// Tensor product of the one-dimensional populations with moments (1, u, cs^2+u^2).
struct ProductEquilibrium {
  using parameters = meta::list<>;
  static std::string getName() { return "ProductEquilibrium"; }
  template<class DESCRIPTOR, class MOMENTA>
  struct type {
    template<class V> static V factor(int direction,V speed) {
      return direction == 0 ? V(2./3.)-speed*speed :
        V(.5)*(V(1./3.)+speed*speed+V(direction)*speed);
    }
    template<class CELL,class RHO,class U,class FEQ,class V=typename CELL::value_t>
    CellStatistic<V> compute(CELL&,RHO& rho,U& u,FEQ& feq) {
      static_assert(DESCRIPTOR::d==2 && DESCRIPTOR::q==9);
      for(int i=0;i<9;++i)
        feq[i]=rho*factor(descriptors::c<DESCRIPTOR>(i,0),u[0])*
                   factor(descriptors::c<DESCRIPTOR>(i,1),u[1])-descriptors::t<V,DESCRIPTOR>(i);
      return {rho,u[0]*u[0]+u[1]*u[1]};
    }
    template<class CELL,class PARAMETERS,class FEQ,class V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell,PARAMETERS&,FEQ& feq) {
      using M=typename MOMENTA::template type<DESCRIPTOR>;
      V rho,u[2]; M().computeRhoU(cell,rho,u);
      return compute(cell,rho,u,feq);
    }
  };
};
template<class T,class DESCRIPTOR,class MOMENTA=momenta::BulkTuple>
using ProductBGK=dynamics::Tuple<T,DESCRIPTOR,MOMENTA,ProductEquilibrium,collision::BGK>;
