// Final mid-channel section. Re-evaluate the example's slip-law reference
// at exactly the same bulk-node coordinates as the computed velocities.
{
  using D = MyCase::descriptor_t;
  auto& lat = myCase.getLattice(NavierStokes{});
  auto& g = myCase.getGeometry().getBlockGeometry(0);
  auto& p = myCase.getParameters();
  const auto extent = p.get<parameters::DOMAIN_EXTENT>();
  const auto& c = lat.getUnitConverter();
  const double H=extent[1], L=extent[0], dp=p.get<parameters::PRESSURE_DIFFERENCE>();
  const double mu=p.get<parameters::PHYS_CHAR_DENSITY>()*p.get<parameters::PHYS_CHAR_VISCOSITY>();
  const double slip=p.get<parameters::SLIPCOEFF>(), creep=p.get<parameters::CREEPCOEFF>();
  const double dT=p.get<parameters::TEMPERATURE_RIGHT>()-p.get<parameters::TEMPERATURE_LEFT>();
  SuperLatticePhysVelocity2D<double,D> velocity(lat,c);
  lat.setProcessingContext(ProcessingContext::Evaluation);
  std::ofstream out("profile.csv");
  out << "y_over_H,u_ms,ref_ms,H_m\n" << std::setprecision(14);
  for (int i=0;i<lat.getBlock(0).getNx();++i) {
    for (int j=0;j<lat.getBlock(0).getNy();++j) {
      auto x=g.getPhysR({i,j});
      if (g.getMaterial(i,j)!=1 || std::abs(x[0]-L/2)>.49*c.getPhysDeltaX()) continue;
      int pos[3]{0,i,j}; double u[2]{};
      velocity(u,pos);
      double ref=dp/(2*mu*L)*(x[1]*(H-x[1])+slip*H)+creep*dT/L;
      out << x[1]/H << ',' << u[0] << ',' << ref << ',' << H << '\n';
    }
  }
}
