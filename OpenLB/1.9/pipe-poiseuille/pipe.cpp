// SPDX-License-Identifier: GPL-2.0-or-later
#include "case.h"   // OpenLB 1.9r0 examples/laminar/poiseuille3d
#include <cstdlib>
#include <fstream>
#include <iomanip>

int main(int argc, char** argv) {
  initialize(&argc,&argv);
  using namespace olb::parameters;
  using D = MyCase::descriptor_t_of<NavierStokes>;
  const char* name = std::getenv("PIPE_WALL");
  const std::string method = name ? name : "bouzidi";
  if (method != "bouzidi" && method != "halfway") return 2;
  MyCase::ParametersD p;
  setDefaultParameters(p);
  p.set<FLOW_TYPE>(FlowType::FORCED);
  p.set<BOUNDARY_TYPE>(BoundaryType::BOUZIDI);
  p.set<LENGTH>(0.25);
  p.fromCLI(argc,argv);
  auto mesh = createMesh(p);
  MyCase sim(p,mesh);
  prepareGeometry(sim);
  prepareLattice(sim);
  auto& lat = sim.getLattice(NavierStokes{});
  auto& geom = sim.getGeometry();
  auto& b = lat.getBlock(0);
  auto& g = geom.getBlockGeometry(0);
  const auto& c = lat.getUnitConverter();
  const double dx = c.getPhysDeltaX(), radius = p.get<DIAMETER>()/2;
  const double maxU = p.get<PHYS_CHAR_VELOCITY>();
  const int plane = b.getNx()/2;
  std::ofstream links("links.csv");
  links << "y_m,z_m,cx,cy,cz,q_original,q_used\n" << std::setprecision(15);
  for (int i=0; i<b.getNx(); ++i) for (int j=0; j<b.getNy(); ++j) for (int k=0; k<b.getNz(); ++k) {
    if (g.getMaterial(i,j,k)!=1) continue;
    auto cell = b.get(i,j,k);
    for (int pop=1; pop<D::q; ++pop) {
      const double q = cell.getFieldComponent<descriptors::BOUZIDI_DISTANCE>(pop);
      if (q<0) continue;
      const double used = method=="halfway" ? .5 : q;
      cell.setFieldComponent<descriptors::BOUZIDI_DISTANCE>(pop,used);
      if (i==plane) {
        const auto x=g.getPhysR({i,j,k});
        links<<x[1]<<','<<x[2]<<','<<descriptors::c<D>(pop,0)<<','<<descriptors::c<D>(pop,1)
             <<','<<descriptors::c<D>(pop,2)<<','<<q<<','<<used<<'\n';
      }
    }
  }
  // Same axial acceleration as the forced official case, starting from rest.
  const double accel = 4*c.getPhysViscosity()*maxU/(radius*radius);
  Vector<double,3> latticeForce{accel*c.getPhysDeltaT()*c.getPhysDeltaT()/dx,0.,0.};
  fields::set<descriptors::FORCE>(lat,geom.getMaterialIndicator({1,2}),latticeForce);
  AnalyticalConst3D<double,double> one(1.), zero(0.,0.,0.);
  lat.iniEquilibrium(geom.getMaterialIndicator(1),one,zero);
  lat.setParameter<descriptors::OMEGA>(c.getLatticeRelaxationFrequency());
  lat.initialize();
  SuperLatticePhysVelocity3D<double,D> vel(lat,c);
  struct Node {int j,k; double y,z,last;};
  std::vector<Node> nodes;
  for(int j=0;j<b.getNy();++j) for(int k=0;k<b.getNz();++k) {
    if(g.getMaterial(plane,j,k)!=1) continue;
    auto x=g.getPhysR({plane,j,k});
    nodes.push_back({j,k,x[1],x[2],0.});
  }
  std::ofstream history("history.csv");
  history<<"step,t_s,Q_m3s,relative_L2_change\n"<<std::setprecision(15);
  int stable=0;
  const int interval=c.getLatticeTime(.25), end=c.getLatticeTime(20.);
  for(int t=1;t<=end;++t) {
    lat.collideAndStream();
    if(t%interval!=0) continue;
    lat.setProcessingContext(ProcessingContext::Evaluation);
    double sum=0.,diff=0.,norm=0.;
    for(auto& node:nodes) {
      int pos[]{0,plane,node.j,node.k}; double u[3]{};
      vel(u,pos); sum+=u[0]; diff+=std::pow(u[0]-node.last,2); norm+=u[0]*u[0]; node.last=u[0];
    }
    const double change=std::sqrt(diff/norm);
    if(!std::isfinite(change)) throw std::runtime_error("Invalid velocity");
    history<<t<<','<<c.getPhysTime(t)<<','<<sum*dx*dx<<','<<change<<std::endl;
    stable=change<1e-9 ? stable+1 : 0;
    if(stable==3) {
      std::ofstream out("section.csv");
      out<<"y_m,z_m,ux_ms,ux_first_ms,ux_last_ms\n"<<std::setprecision(15);
      for(const auto& node:nodes) {
        int a[]{0,0,node.j,node.k}, bpos[]{0,b.getNx()-1,node.j,node.k};
        double first[3]{}, last[3]{}; vel(first,a); vel(last,bpos);
        out<<node.y<<','<<node.z<<','<<node.last<<','<<first[0]<<','<<last[0]<<'\n';
      }
      return out.good() && history.good() && links.good() ? 0 : 1;
    }
    lat.setProcessingContext(ProcessingContext::Simulation);
  }
  throw std::runtime_error("Steady criterion not met in 20 s");
}
