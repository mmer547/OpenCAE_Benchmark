// SPDX-License-Identifier: GPL-2.0-or-later
// Uses the unmodified OpenLB 1.9r0 laminar/cylinder2d case header.
#include "case.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>

int main(int argc, char** argv) {
  initialize(&argc, &argv);
  using namespace olb::parameters;
  using T = MyCase::value_t;
  using D = MyCase::descriptor_t;
  MyCase::ParametersD p;
  p.set<RESOLUTION>(10);
  p.set<REYNOLDS>(20.);
  p.set<DOMAIN_EXTENT>({2.2, 0.41});
  p.set<PHYS_CHAR_VELOCITY>(0.2);
  p.set<LATTICE_RELAXATION_TIME>(0.56);
  p.set<PHYS_CHAR_DENSITY>(1.);
  p.set<RADIUS_CYLINDER>(0.05);
  p.set<CENTER_CYLINDER>({0.2, 0.2});
  p.set<MAX_PHYS_T>(32.);
  p.set<RAMP_UP_END_FRACTION>(0.125);
  p.set<RAMP_UP_UPDATE>(0.01);
  p.set<PHYS_STAT_ITER_T>(0.5);
  p.fromCLI(argc, argv);
  auto mesh = createMesh(p);
  MyCase c(p, mesh);
  prepareGeometry(c);
  prepareLattice(c);
  setInitialValues(c);
  auto& lat = c.getLattice(NavierStokes{});
  auto& geom = c.getGeometry();
  const auto& conv = lat.getUnitConverter();
  const T dx = conv.getPhysDeltaX();
  if (const char* selected = std::getenv("OPENLB_CYLINDER_WALL")) {
    const std::string method(selected);
    if (method != "halfway" && method != "bouzidi")
      throw std::invalid_argument("Unknown cylinder wall method");
    auto& block = lat.getBlock(0);
    auto& geometry = geom.getBlockGeometry(0);
    std::ofstream links("links.csv");
    links << "i,j,cx,cy,q_circle,q_used\n" << std::setprecision(15);
    int count = 0;
    for (int i = 0; i < block.getNx(); ++i) {
      for (int j = 0; j < block.getNy(); ++j) {
        if (geometry.getMaterial(i,j) != 1) continue;
        auto cell = block.get(i,j);
        for (int k = 1; k < D::q; ++k) {
          const int cx = descriptors::c<D>(k,0), cy = descriptors::c<D>(k,1);
          if (geometry.getMaterial(i+cx,j+cy) != 5) continue;
          const T q = cell.getFieldComponent<descriptors::BOUZIDI_DISTANCE>(k);
          if (!(q >= 0 && q <= 1)) throw std::runtime_error("Missing cylinder cut link");
          const T used = method == "halfway" ? T{0.5} : q;
          cell.setFieldComponent<descriptors::BOUZIDI_DISTANCE>(k, used);
          links << i << ',' << j << ',' << cx << ',' << cy << ',' << q << ',' << used << '\n';
          ++count;
        }
      }
    }
    if (!count || !links.good()) throw std::runtime_error("Invalid cylinder link output");
  }
  SuperLatticePhysPressure2D<T,D> pressure(lat, conv);
  AnalyticalFfromSuperF2D<T> sample(pressure, true);
  SuperLatticePhysDrag2D<T,D> drag(lat, geom, 5, conv);
  SuperGeometryFaces2D<T> projected(geom, 5, dx);
  int input[3]{};
  T area[5]{};
  projected(area, input);
  std::ofstream history("history.csv");
  history << "step,time_s,Cd,Cl,dp_Pa\n" << std::setprecision(12);
  const auto end = conv.getLatticeTime(p.get<MAX_PHYS_T>());
  const auto interval = std::max<std::size_t>(1, conv.getLatticeTime(p.get<PHYS_STAT_ITER_T>()));
  for (int step=1; step<=end; ++step) {
    setTemporalValues(c, step-1);
    lat.collideAndStream();
    if (step % interval != 0 && step != end) continue;
    lat.setProcessingContext(ProcessingContext::Evaluation);
    T coeff[2]{};
    drag(coeff, input);
    // Convert the voxel silhouette used by PhysDrag to the physical D times depth dx.
    for (int d=0; d<2; ++d) coeff[d] *= area[d] / (0.1*dx);
    auto fluidPressure = [&](T x) {
      T point[2]{x, 0.2+dx/2}, value{};
      if (!sample(&value, point)) throw std::runtime_error("Pressure sample outside lattice");
      return value;
    };
    // Extrapolate from two fluid columns on each side to the cylinder surface.
    const T front = 2*fluidPressure(0.15-dx) - fluidPressure(0.15-2*dx);
    const T rear  = 2*fluidPressure(0.25+dx) - fluidPressure(0.25+2*dx);
    history << step << ',' << conv.getPhysTime(step) << ',' << coeff[0] << ','
            << coeff[1] << ',' << front-rear << std::endl;
    lat.setProcessingContext(ProcessingContext::Simulation);
  }
  if (std::getenv("OPENLB_WRITE_WAKE")) {
    lat.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysVelocity2D<T,D> velocity(lat,conv);
    auto& block = lat.getBlock(0);
    std::ofstream field("wake.csv");
    field << "x,y,material,u,v\n" << std::setprecision(12);
    for (int j=0; j<block.getNy(); ++j) {
      for (int i=0; i<block.getNx(); ++i) {
        const int material = geom.getBlockGeometry(0).getMaterial(i,j);
        int pos[3]{0,i,j};
        T u[2]{};
        if (material==1 || material==3 || material==4) velocity(u,pos);
        field << i*dx << ',' << (j-0.5)*dx << ',' << material << ',' << u[0] << ',' << u[1] << '\n';
      }
    }
    field.close();
    if (!field) throw std::runtime_error("Failed to write wake field");
  }
  return history.good() ? 0 : 1;
}
