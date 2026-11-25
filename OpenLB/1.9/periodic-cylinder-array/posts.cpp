// SPDX-License-Identifier: GPL-2.0-or-later
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <numbers>

using namespace olb;
using namespace olb::names;
using D = descriptors::D2Q9<descriptors::FORCE>;
using Posts = Case<NavierStokes,Lattice<double,D>>;

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  if(argc!=4) return 2;
  const int n=std::stoi(argv[1]);
  const double fraction=std::stod(argv[2]);
  const std::string method=argv[3];
  if(n<8 || !(fraction>0 && fraction<std::numbers::pi/4) ||
     (method!="halfway" && method!="bouzidi")) return 2;
  const double center=(n-1)/2., radius=n*std::sqrt(fraction/std::numbers::pi);
  IndicatorCuboid2D<double> box({double(n-1),double(n-1)},{0.,0.});
  Mesh<double,2> mesh(box,1.,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  Posts::ParametersD params;
  Posts problem(params,mesh);
  auto& geometry=problem.getGeometry();
  auto& lat=problem.getLattice(NavierStokes{});
  constexpr double nu=.1;
  const double accel=.5*nu*nu/(double(n)*n*n);
  lat.setUnitConverter(1.,1.,double(n),nu/n,nu,1.);
  IndicatorCircle2D<double> cylinder({center,center},radius);
  geometry.rename(0,1);
  geometry.rename(1,2,cylinder);
  dynamics::set<ForcedBGKdynamics>(lat,geometry.getMaterialIndicator(1));
  setBouzidiBoundary(lat,geometry,2,cylinder);
  auto& block=lat.getBlock(0);
  auto& g=geometry.getBlockGeometry(0);
  if(block.getNx()!=n || block.getNy()!=n) throw std::runtime_error("Incorrect period");
  std::vector<std::pair<int,int>> fluid;
  std::ofstream links("links.csv");
  links<<"i,j,cx,cy,q_geometry,q_used\n"<<std::setprecision(15);
  for(int i=0;i<n;++i) for(int j=0;j<n;++j) {
    if(g.getMaterial(i,j)!=1) continue;
    fluid.emplace_back(i,j);
    auto cell=block.get(i,j);
    for(int p=1;p<D::q;++p) {
      auto q=cell.getFieldComponent<descriptors::BOUZIDI_DISTANCE>(p);
      if(q<0) continue;
      const double used=method=="halfway" ? .5 : q;
      cell.setFieldComponent<descriptors::BOUZIDI_DISTANCE>(p,used);
      links<<i<<','<<j<<','<<descriptors::c<D>(p,0)<<','<<descriptors::c<D>(p,1)<<','<<q<<','<<used<<'\n';
    }
  }
  fields::set<descriptors::FORCE>(lat,geometry.getMaterialIndicator(1),Vector<double,2>{accel,0.});
  AnalyticalConst2D<double,double> rho(1.), zero(0.,0.);
  lat.iniEquilibrium(geometry.getMaterialIndicator(1),rho,zero);
  lat.setParameter<descriptors::OMEGA>(1./(.5+3*nu));
  lat.initialize();
  std::ofstream history("history.csv");
  history<<"step,U_superficial,acceleration,K_over_L2,relative_change\n"<<std::setprecision(15);
  const int interval=std::max(1,n*n/8);
  int stable=0;
  double previous=0.;
  for(int step=1;step<=12*n*n;++step) {
    lat.collideAndStream();
    if(step%interval) continue;
    lat.setProcessingContext(ProcessingContext::Evaluation);
    double sum=0.;
    for(auto [i,j]:fluid) {
      double u[2]{}; block.get(i,j).computeU(u); sum+=u[0];
    }
    const double mean=sum/(n*n);
    const double permeability=nu*mean/(accel*n*n);
    const double change=std::abs(1.-previous/permeability);
    if(!(permeability>0.) || !std::isfinite(change)) throw std::runtime_error("Invalid permeability");
    previous=permeability;
    history<<step<<','<<mean<<','<<accel<<','<<permeability<<','<<change<<std::endl;
    stable=change<1e-8 ? stable+1 : 0;
    if(stable>=3) {
      std::ofstream out("field.csv");
      out<<"i,j,ux,uy\n"<<std::setprecision(15);
      for(auto [i,j]:fluid) {
        double u[2]{}; block.get(i,j).computeU(u);
        out<<i<<','<<j<<','<<u[0]<<','<<u[1]<<'\n';
      }
      return out.good() && history.good() && links.good() ? 0 : 1;
    }
    lat.setProcessingContext(ProcessingContext::Simulation);
  }
  throw std::runtime_error("Steady permeability criterion not met");
}
