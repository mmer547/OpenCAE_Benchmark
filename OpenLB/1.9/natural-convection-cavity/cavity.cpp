// SPDX-License-Identifier: GPL-2.0-or-later
// Uses the OpenLB 1.9r0 thermal/squareCavity2d/laminar coupling and wall models.
#include "case.h"
#include <fstream>
#include <iomanip>

int main(int argc, char** argv) {
  initialize(&argc, &argv);
  using namespace olb::parameters;
  using T = MyCase::value_t;
  using NS = MyCase::descriptor_t_of<NavierStokes>;
  constexpr T nu = 0.001, pr = 0.71, alpha = nu/pr, gravity = 9.81;
  constexpr T beta = 1e4*nu*alpha/gravity;
  MyCase::ParametersD p;
  p.set<RESOLUTION>(32);
  p.set<DOMAIN_EXTENT>({1., 1.});
  p.set<ORIGIN>({0., 0.});
  p.set<PHYS_CHAR_LENGTH>(1.);
  p.set<PHYS_CHAR_VELOCITY>(std::sqrt(gravity*beta));
  p.set<LATTICE_RELAXATION_TIME>(0.65);
  p.set<PHYS_KINEMATIC_VISCOSITY>(nu);
  p.set<PHYS_CHAR_DENSITY>(1.);
  p.set<PHYS_HEAT_CAPACITY>(1.);
  p.set<PHYS_THERMAL_CONDUCTIVITY>(alpha);
  p.set<PHYS_THERMAL_EXPANSION>(beta);
  p.set<GRAVITATIONAL_ACC>(gravity);
  p.set<T_HOT>(300.5);
  p.set<T_COLD>(299.5);
  p.fromCLI(argc, argv);
  const int n = p.get<RESOLUTION>();
  const T dx = 1./n;
  p.set<PHYS_DELTA_X>(dx);
  // Vertical Dirichlet walls on nodes at x=0,1. Horizontal bounce-back walls
  // halfway between nodes at y=0,1. The n samples along a side use midpoint weights.
  IndicatorCuboid2D<T> box({1., 1.+dx}, {0., -dx/2});
  Mesh<T,2> mesh(box, dx, 1);
  mesh.setOverlap(p.get<OVERLAP>());
  MyCase c(p, mesh);
  auto& geom = c.getGeometry();
  geom.rename(0,4);
  geom.rename(4,1,{1,1});
  IndicatorCuboid2D<T> left({dx,1.},{-dx/2,0.});
  IndicatorCuboid2D<T> right({dx,1.},{1.-dx/2,0.});
  geom.rename(4,2,1,left);
  geom.rename(4,3,1,right);
  geom.checkForErrors();
  geom.print();
  prepareLattice(c);
  setInitialValues(c);
  auto& flow = c.getLattice(NavierStokes{});
  auto& heat = c.getLattice(Temperature{});
  const auto& conv = flow.getUnitConverter();
  SuperLatticePhysVelocity2D<T,NS> velocity(flow, conv);
  AnalyticalFfromSuperF2D<T> sample(velocity, true);
  auto& g = geom.getBlockGeometry(0);
  auto& h = heat.getBlock(0);
  std::ofstream history("history.csv");
  history << "step,time_s,Nu_hot,Nu_cold,u_max,v_max\n" << std::setprecision(12);
  std::array<T,4> previous{};
  int stable = 0;
  const int interval = conv.getLatticeTime(5.);
  for (int step=1; step<=conv.getLatticeTime(500.); ++step) {
    flow.collideAndStream();
    heat.collideAndStream();
    c.getOperator("Boussinesq").apply();
    if (step % interval != 0) continue;
    flow.setProcessingContext(ProcessingContext::Evaluation);
    heat.setProcessingContext(ProcessingContext::Evaluation);
    std::array<T,4> values{};
    for (int side=0; side<2; ++side) {
      const int x = side == 0 ? 0 : n, inward = side == 0 ? 1 : -1;
      for (int y=1; y<=n; ++y) {
        if (g.getMaterial(x,y) != side+2) throw std::runtime_error("Wall node mismatch");
        auto temp = [&](int offset) {
          return conv.getPhysTemperature(h.get(x+inward*offset,y).computeRho());
        };
        values[side] += inward*(3*temp(0)-4*temp(1)+temp(2))/(2*dx*n);
      }
    }
    for (int i=0; i<=2*n; ++i) {
      const T pos = i/(2.*n);
      T vertical[2]{0.5,pos}, horizontal[2]{pos,0.5}, u[2]{}, v[2]{};
      if (!sample(u,vertical) || !sample(v,horizontal)) throw std::runtime_error("Invalid velocity sample");
      values[2] = std::max(values[2],u[0]/alpha);
      values[3] = std::max(values[3],v[1]/alpha);
    }
    history << step << ',' << conv.getPhysTime(step);
    bool unchanged = true;
    for (int i=0; i<4; ++i) {
      if (!std::isfinite(values[i])) throw std::runtime_error("Nonfinite observable");
      history << ',' << values[i];
      unchanged &= std::abs(values[i]-previous[i]) < 1.e-6*std::abs(values[i]);
    }
    history << std::endl;
    previous = values;
    stable = unchanged ? stable+1 : 0;
    if (stable >= 4 && conv.getPhysTime(step) >= 80.) return history.good() ? 0 : 1;
    flow.setProcessingContext(ProcessingContext::Simulation);
    heat.setProcessingContext(ProcessingContext::Simulation);
  }
  throw std::runtime_error("Steady criterion not reached by 500 s");
}
