#include <TMB.hpp>
#include <TBM.hpp>

template <class Type>
Type objective_function<Type>::operator()() {
  Type dll = 0.0;
  DATA_VECTOR(init); // initial conds
  DATA_SCALAR(tmax);
  DATA_SCALAR(dt);
  PARAMETER_VECTOR(pars);
  ODE<Type, TB<Type> > ode(init, pars, asDouble(tmax), asDouble(dt));
  matrix<double> out = ode.out();
  REPORT(ode.track);
  REPORT(out);
  return dll;
}
