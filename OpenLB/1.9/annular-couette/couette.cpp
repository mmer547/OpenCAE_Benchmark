// SPDX-License-Identifier: GPL-2.0-or-later
#include <olb.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace olb;
using namespace olb::names;
using T = double;
using D = descriptors::D2Q9<>;
using Annulus = Case<NavierStokes, Lattice<T,D>>;

int main(int argc, char** argv) {
  initialize(&argc,&argv);
  Annulus::ParametersD params;
  params.set<parameters::RESOLUTION>(16);
  params.fromCLI(argc,argv);
  const int n = params.get<parameters::RESOLUTION>();
  const T dx = 1./n, half = 2.+2.5*dx;
  IndicatorCuboid2D<T> box({2*half,2*half},{-half,-half});
  Mesh<T,2> mesh(box,dx,1);
  mesh.setOverlap(3);
  Annulus problem(params,mesh);
  auto& geom = problem.getGeometry();
  auto& lat = problem.getLattice(NavierStokes{});
  IndicatorCircle2D<T> inner({0.,0.},1.), outer({0.,0.},2.);
  geom.rename(0,3);
  geom.rename(3,1,outer);
  geom.rename(1,2,inner);
  geom.clean();
  geom.checkForErrors();
  geom.print();
  lat.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,D>>(n,.8,1.,1.,1.,1.);
  const auto& conv = lat.getUnitConverter();
  conv.print();
  dynamics::set<ConstRhoBGKdynamics>(lat,geom.getMaterialIndicator(1));
  setBouzidiBoundary<T,D,BouzidiVelocityPostProcessor>(lat,geom,2,inner);
  setBouzidiBoundary(lat,geom,3,outer);
  if (const char* selected = std::getenv("OPENLB_WALL_SCHEME")) {
    const std::string scheme(selected);
    if (scheme != "halfway" && scheme != "bouzidi")
      throw std::invalid_argument("Unknown wall scheme");
    auto& block = lat.getBlock(0);
    auto& geometry = geom.getBlockGeometry(0);
    std::ofstream links("links.csv");
    links << "x,y,cx,cy,q_circle,q_used\n" << std::setprecision(14);
    for (int i = 0; i < block.getNx(); ++i) {
      for (int j = 0; j < block.getNy(); ++j) {
        if (geometry.getMaterial(i,j) != 1) continue;
        auto cell = block.get(i,j);
        const auto p = geometry.getPhysR({i,j});
        for (int k = 1; k < D::q; ++k) {
          const T q = cell.getFieldComponent<descriptors::BOUZIDI_DISTANCE>(k);
          if (q < 0) continue;
          const T used = scheme == "halfway" ? T{0.5} : q;
          cell.setFieldComponent<descriptors::BOUZIDI_DISTANCE>(k,used);
          links << p[0] << ',' << p[1] << ',' << descriptors::c<D>(k,0) << ','
                << descriptors::c<D>(k,1) << ',' << q << ',' << used << '\n';
        }
      }
    }
    if (!links.good()) throw std::runtime_error("Could not save wall links");
  }
  AnalyticalFfromCallableF<2,T,T> rotation([&](Vector<T,2> p) -> Vector<T,2> {
    return {-p[1]/conv.getConversionFactorVelocity(),p[0]/conv.getConversionFactorVelocity()};
  });
  setBouzidiVelocity(lat,geom.getMaterialIndicator(2),geom.getMaterialIndicator(1),rotation);
  lat.setParameter<descriptors::OMEGA>(conv.getLatticeRelaxationFrequency());
  lat.initialize();
  SuperLatticePhysVelocity2D<T,D> velocity(lat,conv);
  SuperLatticePhysPressure2D<T,D> pressure(lat,conv);
  auto& g = geom.getBlockGeometry(0);
  auto& b = lat.getBlock(0);
  struct Node { int i,j; T x,y; };
  std::vector<Node> fluid;
  for (int j=0;j<b.getNy();++j) {
    for (int i=0;i<b.getNx();++i) {
      if (g.getMaterial(i,j)==1) {
        auto p=g.getPhysR({i,j});
        fluid.push_back({i,j,p[0],p[1]});
      }
    }
  }
  std::vector<Vector<T,2>> previous(fluid.size());
  std::ofstream history("history.csv");
  history << "step,time_s,relative_velocity_change,mean_lattice_density\n" << std::setprecision(12);
  const auto interval=conv.getLatticeTime(.25);
  int stable=0;
  for (std::size_t step=1;step<=conv.getLatticeTime(12.);++step) {
    lat.collideAndStream();
    if (step%interval) continue;
    lat.setProcessingContext(ProcessingContext::Evaluation);
    T change=0., norm=0.;
    for (std::size_t k=0;k<fluid.size();++k) {
      int pos[3]{0,fluid[k].i,fluid[k].j};
      Vector<T,2> u;
      velocity(u.data(),pos);
      const auto diff=u-previous[k];
      change+=diff*diff;
      norm+=u*u;
      previous[k]=u;
    }
    const T relative=std::sqrt(change/norm);
    if (!std::isfinite(relative)) throw std::runtime_error("Invalid velocity field");
    history << step << ',' << conv.getPhysTime(step) << ',' << relative << ','
            << lat.getStatistics().getAverageRho() << std::endl;
    stable=relative<1.e-8 ? stable+1 : 0;
    if (stable>=3 && conv.getPhysTime(step)>=2.) {
      std::ofstream out("field.csv");
      out << "x,y,u,v,p_Pa\n" << std::setprecision(12);
      for (std::size_t k=0;k<fluid.size();++k) {
        const auto& node=fluid[k];
        int pos[3]{0,node.i,node.j};
        T p=0.;
        pressure(&p,pos);
        out << node.x << ',' << node.y << ',' << previous[k][0] << ','
            << previous[k][1] << ',' << p << '\n';
      }
      out.close();
      return out.good() && history.good() ? 0 : 1;
    }
    lat.setProcessingContext(ProcessingContext::Simulation);
  }
  throw std::runtime_error("Steady criterion not reached by 12 s");
}
