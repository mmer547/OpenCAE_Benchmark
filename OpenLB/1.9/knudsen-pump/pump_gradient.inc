// Impose the known tangential gradient of the prescribed wall temperature in
// lattice units. Away from the shoulders, use the actual radial pipe normal.
{
  auto& g=geometry.getBlockGeometry(0);
  auto& block=lattice.getBlock(0);
  const T mid=(extend[1]-physDeltaX)/2, length=extend[0];
  for(int i=0;i<block.getNx();++i) for(int j=0;j<block.getNy();++j) for(int k=0;k<block.getNz();++k){
    if(g.getMaterial(i,j,k)!=1) continue;
    auto cell=block.get(i,j,k);
    auto normal=cell.template getField<descriptors::NORMAL>();
    if(util::norm<3>(normal)<T(.5)) continue;
    auto x=g.getPhysR({i,j,k});
    if(x[0]>physDeltaX && x[0]<length-physDeltaX){
      T r=std::hypot(x[1]-mid,x[2]-mid);
      if(r>0){
        normal=Vector<T,3>{0,(mid-x[1])/r,(mid-x[2])/r};
        cell.template setField<descriptors::NORMAL>(normal);
      }
    }
    T span=std::clamp(x[0]+physDeltaX/2,T(0),length)-std::clamp(x[0]-physDeltaX/2,T(0),length);
    Vector<T,3> grad{(tempRight-tempLeft)*span/length,0,0};
    cell.template setField<descriptors::TEMPGRADIENT>(grad-(grad*normal)*normal);
  }
}
