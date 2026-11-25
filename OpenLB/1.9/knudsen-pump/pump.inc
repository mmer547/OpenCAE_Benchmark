{
  using D=MyCase::descriptor_t_of<NavierStokes>;
  auto& lat=myCase.getLattice(NavierStokes{});
  auto& g=myCase.getGeometry().getBlockGeometry(0);
  auto& p=myCase.getParameters();
  const auto& c=lat.getUnitConverter();
  auto extent=p.get<parameters::DOMAIN_EXTENT>();
  const double L=extent[0],H=extent[1],dx=c.getPhysDeltaX(),yc=(H-dx)/2;
  const double dpRef=6*std::pow(c.getPhysDensity()*c.getPhysViscosity(),2)*p.get<parameters::R_SPECIFIC>()
    *(p.get<parameters::TEMPERATURE_RIGHT>()-p.get<parameters::TEMPERATURE_LEFT>())
    /(std::pow(H/2,2)*p.get<parameters::AVERAGE_PRESSURE>()*(1+8*p.get<parameters::SLIPCOEFF>()/H));
  SuperLatticePhysVelocity3D<double,D> velocity(lat,c);
  SuperLatticePhysPressure3D<double,D> pressure(lat,c);
  auto& b=lat.getBlock(0);
  auto sample = [&](bool final) {
    lat.setProcessingContext(ProcessingContext::Evaluation);
    std::ofstream section,line;
    if (final) {
      section.open("pump_section.csv"); line.open("pump_axis.csv");
      section<<"r_over_R,u_ms,rho_kgm3\n"<<std::setprecision(14);
      line<<"x_over_L,p_Pa,p_ref_Pa\n"<<std::setprecision(14);
    }
    double pl=0.,pr=0.,flux=0.,absflux=0.;int nl=0,nr=0;
    for (int i=0;i<b.getNx();++i) for (int j=0;j<b.getNy();++j) for (int k=0;k<b.getNz();++k) {
      if(g.getMaterial(i,j,k)!=1) continue;
      auto x=g.getPhysR({i,j,k});
      int pos[4]{0,i,j,k}; double press{},u[3]{};
      pressure(&press,pos);
      if(x[0]<-.1*L && x[0]>-.4*L) {pl+=press;++nl;}
      if(x[0]>1.1*L && x[0]<1.4*L) {pr+=press;++nr;}
      if(final && std::hypot(x[1]-yc,x[2]-yc)<.72*dx)
        line<<x[0]/L<<','<<press<<','<<p.get<parameters::AVERAGE_PRESSURE>()+dpRef*(std::clamp(x[0]/L,0.,1.)-.5)<<'\n';
      if(std::abs(x[0]-.5*L)<.49*dx) {
        velocity(u,pos);
        double rho=c.getPhysDensity()*b.get(i,j,k).computeRho();
        flux+=rho*u[0]*dx*dx; absflux+=rho*std::abs(u[0])*dx*dx;
        if(final) section<<2*std::hypot(x[1]-yc,x[2]-yc)/H<<','<<u[0]<<','<<rho<<'\n';
      }
    }
    if (!nl || !nr || !absflux) throw std::runtime_error("Missing pump sampling plane");
    return std::array<double,3>{pr/nr-pl/nl,flux,absflux};
  };
  // The example's stopping metric is maximum speed. Check the pressure and
  // signed flux directly over subsequent windows using its same time step.
  auto values=sample(false);
  std::ofstream history("pump_history.csv");
  history<<"extra_s,dp_Pa,net_kgs,abs_kgs\n"<<std::setprecision(14);
  history<<0<<','<<values[0]<<','<<values[1]<<','<<values[2]<<'\n';
  const int window=std::max(1,int(c.getLatticeTime(4.e-7)));
  int stable=0;
  for (int w=1;w<=20 && stable<3;++w) {
    lat.setProcessingContext(ProcessingContext::Simulation);
    for (int t=0;t<window;++t) {
      lat.collideAndStream();
      lat.stripeOffDensityOffset(lat.getStatistics().getAverageRho()-1.);
    }
    const auto next=sample(false);
    stable=(std::abs(next[0]-values[0])<.001*dpRef &&
            std::abs(next[1]-values[1])<.001*next[2]) ? stable+1 : 0;
    values=next;
    history<<c.getPhysTime(w*window)<<','<<values[0]<<','<<values[1]<<','<<values[2]<<'\n';
  }
  values=sample(true);
  std::ofstream out("pump.csv");
  out<<"dp_Pa,ref_Pa,net_kgs,abs_kgs\n"<<std::setprecision(14)<<values[0]<<','<<dpRef<<','<<values[1]<<','<<values[2]<<'\n';
}
