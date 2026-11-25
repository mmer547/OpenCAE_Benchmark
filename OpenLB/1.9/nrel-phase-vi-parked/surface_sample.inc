  // Evaluate pressure at the wall from two exterior stencils one grid
  // spacing apart along x. Both stencils must stay on the fluid side.
  auto wallSample = [&](double x, double y, double z, int side) {
    auto near = sample(x,y,z,side);
    double nearX = 0;
    for (const auto& [p,w] : near.nodes) nearX += w*(o[0]+dx*p[0]);
    double offset = side*(nearX-x)/dx;
    for (const auto& [p,w] : near.nodes) {
      auto farther = p;
      farther[0] += side;
      if (g.getMaterial(farther)!=1)
        throw std::runtime_error("Exterior extrapolation crossed solid");
    }
    Sample result;
    for (const auto& [p,w] : near.nodes) {
      auto farther = p;
      farther[0] += side;
      result.nodes.emplace_back(p,(1+offset)*w);
      result.nodes.emplace_back(farther,-offset*w);
    }
    return result;
  };
