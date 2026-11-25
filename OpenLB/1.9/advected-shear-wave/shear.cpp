// SPDX-License-Identifier: GPL-2.0-or-later
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <numbers>
using namespace olb;
using namespace olb::names;
using D=descriptors::D2Q9<>;
using Wave=Case<NavierStokes,Lattice<double,D>>;

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  if (argc!=3) return 2;
  int n=std::stoi(argv[1]);
  double boost=std::stod(argv[2]);
  IndicatorCuboid2D<double> box({double(n-1),3.},{0.,0.});
  Mesh<double,2> mesh(box,1.,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  Wave::ParametersD params;
  Wave problem(params,mesh);
  auto& g=problem.getGeometry();
  auto& lat=problem.getLattice(NavierStokes{});
  lat.setUnitConverter(1.,1.,double(n),.1,.1,1.);
  g.rename(0,1);
  dynamics::set<BGKdynamics>(lat,g.getMaterialIndicator(1));
  constexpr double nu=.1, amp=.001;
  lat.setParameter<descriptors::OMEGA>(1./(.5+3*nu));
  auto& b=lat.getBlock(0);
  if (b.getNx()!=n) throw std::runtime_error("Unexpected wave period");
  const double k=2*std::numbers::pi/n;
  AnalyticalConst2D<double,double> rho(1.);
  AnalyticalFfromCallableF<2,double,double> vel([&](Vector<double,2> x)->Vector<double,2> {
    return {boost,amp*std::sin(k*x[0])};
  });
  auto fluid=g.getMaterialIndicator(1);
  lat.iniEquilibrium(fluid,rho,vel);
  lat.initialize();
  std::ofstream out("wave.csv");
  out<<"step,amplitude,continuum\n"<<std::setprecision(15);
  const int stop=std::ceil(1.5/(nu*k*k));
  for (int t=0;t<=stop;++t) {
    if (t%5==0 || t==stop) {
      double cs=0.,sn=0.;
      for (int i=0;i<n;++i) {
        double u[2]{};
        b.get(i,1).computeU(u);
        cs+=u[1]*std::cos(k*i);
        sn+=u[1]*std::sin(k*i);
      }
      out<<t<<','<<2*std::hypot(cs,sn)/n<<','<<amp*std::exp(-nu*k*k*t)<<'\n';
    }
    if (t!=stop) lat.collideAndStream();
  }
}
