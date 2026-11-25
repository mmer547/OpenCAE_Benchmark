// Diagnostics at the same evaluated state as the release pressure measurement.
{
  SuperLatticeVelocity2D<T,DESCRIPTOR> velocity(sLattice1);
  auto& block=sLattice1.getBlock(0);
  T maximum=0.,mass2=0.;
  std::ofstream section;
  if (iT==60000) {
    section.open("diameter.csv");
    section<<"x_lu,fluid2,speed_lu\n"<<std::setprecision(14);
  }
  for (int i=0;i<block.getNx();++i) {
    for (int j=0;j<block.getNy();++j) {
      int pos[3]{0,i,j}; T u[2]{},fraction{};
      velocity(u,pos); c2(&fraction,pos);
      const T speed=std::hypot(u[0],u[1]);
      maximum=std::max(maximum,speed);
      mass2+=fraction;
      if (iT==60000 && j==block.getNy()/2) section<<i<<','<<fraction<<','<<speed<<'\n';
    }
  }
  std::ofstream out("laplace.csv",iT==0?std::ios::out:std::ios::app);
  if (iT==0) out<<"step,dp,dp_ref,max_speed,mass2\n";
  out<<std::setprecision(14)<<iT<<','<<pressureIn-pressureOut<<','<<alpha*(kappa1+kappa2)/(6*radius)
     <<','<<maximum<<','<<mass2<<'\n';
}
