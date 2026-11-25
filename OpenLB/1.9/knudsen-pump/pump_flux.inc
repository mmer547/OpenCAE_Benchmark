// At the last simulated time, measure mass crossing fluid-fluid links at the
// mid-pipe plane. At a node, incoming populations encode the just-streamed flux.
{
  using D=MyCase::descriptor_t_of<NavierStokes>;
  auto& lat=myCase.getLattice(NavierStokes{});
  auto& b=lat.getBlock(0);
  auto& g=myCase.getGeometry().getBlockGeometry(0);
  auto& p=myCase.getParameters();
  auto& c=lat.getUnitConverter();
  auto L=p.get<parameters::DOMAIN_EXTENT>()[0];
  double forwardMinusReverse=0,absoluteVelocityFlux=0;
  int pairCount=0;
  for(int i=0;i<b.getNx()-1;++i){
    if(std::abs(g.getPhysR({i,0,0})[0]-L/2)>.49*c.getPhysDeltaX()) continue;
    for(int j=0;j<b.getNy();++j) for(int k=0;k<b.getNz();++k){
      if(g.getMaterial(i,j,k)!=1) continue;
      auto left=b.get(i,j,k);
      double u[3]; left.computeU(u);
      absoluteVelocityFlux+=left.computeRho()*std::abs(u[0]);
      for(int q=1;q<D::q;++q){
        if(descriptors::c<D>(q,0)!=1) continue;
        const int y=j+descriptors::c<D>(q,1),z=k+descriptors::c<D>(q,2);
        if(g.getMaterial(i+1,y,z)!=1) continue;
        auto right=b.get(i+1,y,z);
        forwardMinusReverse+=right[q]-left[descriptors::opposite<D>(q)];
        ++pairCount;
      }
    }
  }
  if(pairCount==0 || !(absoluteVelocityFlux>0)) throw std::runtime_error("Missing pump cross-section");
  const double unit=c.getPhysDensity()*std::pow(c.getPhysDeltaX(),3)/c.getPhysDeltaT();
  std::ofstream out("link_flux.csv");
  out<<"pairs,net_kg_s,abs_macro_kg_s,percent\n"<<std::setprecision(15)
     <<pairCount<<','<<forwardMinusReverse*unit<<','<<absoluteVelocityFlux*unit<<','<<100*forwardMinusReverse/absoluteVelocityFlux<<'\n';
}
