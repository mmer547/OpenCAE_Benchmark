// SPDX-License-Identifier: GPL-2.0-or-later
#include <olb.h>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace olb;
using namespace olb::names;
using T = float;
using D = descriptors::D3Q19<>;
using Tunnel = Case<NavierStokes,Lattice<T,D>>;
struct Sample { std::vector<std::pair<LatticeR<3>,double>> nodes; };

int main(int argc,char** argv) {
  initialize(&argc,&argv);
  if(argc!=6) return 2;
  const T dx=std::stof(argv[1]), ul=.065f, uPhys=20.f;
  const bool halfway=std::string(argv[2])=="halfway";
  if(!halfway && std::string(argv[2])!="bouzidi") return 2;
  const int last=std::stoi(argv[5]);
  auto count=[&](T length){return int(std::round(length/dx));};
  const int nx=count(7),ny=count(7),nz=count(14);
  const T dt=dx*ul/uPhys;
  Vector<T,3> origin{-dx*(nx/3),-dx*(ny/2),-dx*(nz/2)};
  IndicatorCuboid3D<T> box({dx*(nx-1),dx*(ny-1),dx*(nz-1)},origin);
  Mesh<T,3> mesh(box,dx,1);
  mesh.setOverlap(2);
  mesh.getCuboidDecomposition().setPeriodicity({true,true,true});
  Tunnel::ParametersD params;
  Tunnel sim(params,mesh);
  auto& geo=sim.getGeometry();
  auto& lat=sim.getLattice(NavierStokes{});
  lat.setUnitConverter(dx,dt,T(.625),uPhys,T(17.7875e-6/1.2318),T(1.2318));
  STLreader<T> blades(argv[3],dx,1,RayMode::FastRayX);
  geo.rename(0,1);
  geo.rename(1,2,blades);
  dynamics::set<SmagorinskyBGKdynamics>(lat,geo.getMaterialIndicator(1));
  setBouzidiBoundary(lat,geo,2,blades);
  auto& block=lat.getBlock(0);
  auto& g=geo.getBlockGeometry(0);
  const auto o=g.getPhysR({0,0,0});
  size_t links=0, offMid=0, solids=0;
  for(int i=0;i<block.getNx();++i) for(int j=0;j<block.getNy();++j) for(int k=0;k<block.getNz();++k) {
    if(g.getMaterial(i,j,k)==2) {++solids;continue;}
    auto xyz=g.getPhysR({i,j,k});
    if(std::abs(xyz[0])>.6 || std::abs(xyz[1])>.85 || std::abs(xyz[2])>5.3) continue;
    auto cell=block.get(i,j,k);
    for(int pop=1;pop<D::q;++pop) {
      T q=cell.getFieldComponent<descriptors::BOUZIDI_DISTANCE>(pop);
      if(q<0) continue;
      ++links;
      if(std::abs(q-.5f)>.001f) ++offMid;
      if(halfway) cell.setFieldComponent<descriptors::BOUZIDI_DISTANCE>(pop,.5f);
    }
  }
  if(links<100 || offMid<links/4 || solids<100) throw std::runtime_error("Invalid surface links");
  std::cout<<"domain "<<block.getNx()<<' '<<block.getNy()<<' '<<block.getNz()<<" dx "<<dx<<" dt "<<dt<<" solid "<<solids<<" links "<<links<<" offmid "<<offMid<<std::endl;
  lat.setParameter<descriptors::OMEGA>(lat.getUnitConverter().getLatticeRelaxationFrequency());
  lat.setParameter<collision::LES::SMAGORINSKY>(T(.18));
  AnalyticalConst3D<T,T> rho(1), flow(ul,0,0);
  lat.iniEquilibrium(geo.getMaterialIndicator(1),rho,flow);
  lat.initialize();
  // One-sided trilinear pressure interpolation along the streamwise axis.
  // Use the first entirely fluid stencil; in particular never average through the blade.
  auto sample=[&](double x,double y,double z,int side) {
    Sample s;
    for(int offset=0;offset<=5;++offset) {
      std::array<double,3> a={(x-o[0])/dx+side*(.5+offset*.5),(y-o[1])/dx,(z-o[2])/dx};
      s.nodes.clear();
      bool ok=true;
      for(int di=0;di<2;++di) for(int dj=0;dj<2;++dj) for(int dk=0;dk<2;++dk){
        LatticeR<3> p{int(std::floor(a[0]))+di,int(std::floor(a[1]))+dj,int(std::floor(a[2]))+dk};
        double weight=1;
        for(int dim=0;dim<3;++dim) weight*=1-std::abs(a[dim]-p[dim]);
        if(weight<1e-8) continue;
        if(g.getMaterial(p)!=1) ok=false;
        s.nodes.emplace_back(p,weight);
      }
      if(ok) return s;
    }
    throw std::runtime_error("No exterior pressure stencil");
  };
  std::ifstream file(argv[4]); std::string line; std::getline(file,line);
  std::vector<Sample> probes;
  while(std::getline(file,line)){
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line); int side; double xc,x,y,z; std::string exp;
    if(!(row>>side>>xc>>exp>>x>>y>>z)) return 3;
    probes.push_back(sample(x,y,z,side));
  }
  if(probes.size()!=124) throw std::runtime_error("Unexpected tap count");
  // Incoming distribution is prescribed by a thin equilibrium reservoir at the
  // periodic streamwise seam. The reference pressure is measured in the inlet stream.
  const auto ref=sample(-1.8,1.8,2.36,-1);
  auto density=[&](const Sample& s){double v=0;for(auto& [p,w]:s.nodes) v+=w*block.get(p).computeRho();return v;};
  std::ofstream out("pressure.csv"); out<<std::setprecision(8)<<"step,time,rho_ref";
  for(size_t i=0;i<probes.size();++i) out<<",p"<<i;
  out<<'\n';
  const T u[3]{ul,0,0};
  for(int t=1;t<=last;++t){
    lat.collideAndStream();
    for(int i=0;i<2;++i) for(int j=0;j<block.getNy();++j) for(int k=0;k<block.getNz();++k)
      block.get(i,j,k).iniEquilibrium(T(1),u);
    if(t%20) continue;
    auto refRho=density(ref);
    if(!std::isfinite(refRho)||std::abs(refRho-1)>.08) throw std::runtime_error("Unstable reference pressure");
    out<<t<<','<<t*dt<<','<<refRho;
    for(auto& probe:probes){
      const double cp=(density(probe)-refRho)*2/(3*ul*ul);
      if(!std::isfinite(cp)||std::abs(cp)>30) throw std::runtime_error("Unstable probe pressure");
      out<<','<<cp;
    }
    out<<std::endl;
  }
}
