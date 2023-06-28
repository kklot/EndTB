#include <TMB.hpp>
#include <TBM.hpp>

template <class Type>
vector<Type> Double2Type(vector<double> x){
  vector<Type> y(x.size());
  double * px = x.data();
  for(int i=0;i<x.size();i++) y[i] = Type(*(px + i));
  return y;
}

// https://mc-stan.org/docs/reference-manual/lower-bound-transform.html
// https://discourse.mc-stan.org/t/half-normal-half-cauchy-and-half-t/17314
template <class Type> Type half_normal(Type y, Type mu, Type sd) { return(dnorm(exp(y), mu, sd, true) + y); }
// https://mc-stan.org/docs/stan-users-guide/changes-of-variables.html
template <class T> T log_normal(T y, T mu, T sd) { return(dnorm(y, mu, sd, true) - y); }
// https://github.com/kklot/naomi/blob/94d34246144e4dfcb86161258faf213a7db03268/src/tmb.cpp#L258
template <class T> T prop_beta(T x, T a, T b) { return(log(x) - log(1 - x) + dbeta(x, a, b, true)); }
// beta mean selection
template <class T> T bBeta(T a, T mu) { return((a - mu * a) / mu); }
// Gamma transform
template <class T> T norm2gamma(T x, T a, T b) { return( qgamma(pnorm(x, T(0), T(1)), a, b) ); }
// Lognormal likelihood
template <class Type> Type log_normal_lpdf(Type x, Type meanlog, Type sdlog) {
  return( - pow(log(x + DBL_EPSILON) - meanlog, 2) / (2 * sdlog * sdlog) - log(x * sdlog * sqrt(2*M_PI)) );
}


template <class Type>
Type objective_function<Type>::operator()() {
  Type dll = 0.0;
  // Meta data
  DATA_VECTOR(init); // initial conds
  DATA_SCALAR(tmax);
  DATA_SCALAR(pop1970);
  DATA_SCALAR(year_zero); // maybe estimate this
  DATA_SCALAR(dt);
  double dbdt = asDouble(dt), len_dt = 1/dbdt;
  DATA_IVECTOR(nullid); // indices of parameters in null case
  
  // Target calibration
  DATA_VECTOR(notification_year); // year of the notification in original format eg... 1999
  DATA_VECTOR(notification_meanlog); // notification rate / 100.000
  DATA_VECTOR(notification_sdlog); // notification rate / 100.000
  DATA_VECTOR(mortality_year); // coresponding year - 1990....
  DATA_VECTOR(mortality_meanlog); // mortality rate / 100K from CID
  DATA_VECTOR(mortality_sdlog); // mortality rate / 100K from CID
  DATA_VECTOR(Treat_year); // mortality rate / 100K from CID
  DATA_VECTOR(Treat_qnorm); // mortality rate / 100K from CID
  DATA_VECTOR(Treat_sd); // mortality rate / 100K from CID
  
  DATA_VECTOR(fallback_parameters);
  DATA_IVECTOR(fit_parameters_id);

  PARAMETER_VECTOR(parameters);
  
  vector<Type> pars = fallback_parameters;
  for (int i = 0; i < parameters.size(); i++)
    pars[fit_parameters_id[i]-1] = parameters[i];

  vector<Type> pars_null = pars; // to flexibly change the index
  for (int i = 0; i < nullid.size(); i++) pars_null[nullid[i]] = Type(0); // set zero to health system parameters
  

  ODE<Type, TB<Type> > mod0(init, pars_null, asDouble(tmax), dbdt); // run steady state 
  matrix<double> out0 = mod0.out(); // get the equilibrium
  vector<double> eqVec = out0(Eigen::seqN(1, init.size()), Eigen::last);
  vector<Type> eqVecT = Double2Type<Type>(eqVec); 
  
  ODE<Type, TB<Type>> mod(eqVecT, pars, asDouble(2030 - year_zero), dbdt); // rerun with eq as init, first time point represents 1970: maybe not?

  // extract expected data
  matrix<double> out = mod.out(); // get the equilibrium

  Type lhd = 0;
  vector<Type> // note that the output has an extra row for time
    ept = Double2Type<Type>(out.row(17+1)),  // index of notification model dependent 
    emr = Double2Type<Type>(out.row(18+1)),
    pop = Double2Type<Type>(out.row(22+1))
    ;  // index of mortality
    // eTc = Double2Type<Type>(out.row(11+1)), // index of mortality
    // eTf = Double2Type<Type>(out.row(12+1)); // index of mortality

  for (int i = 0; i < notification_meanlog.size(); i++) {
    int idx = asDouble((notification_year[i] - year_zero) * len_dt);
    Eigen::ArithmeticSequence ii = Eigen::seqN(idx, len_dt);
    vector<Type> annual_v = ept(ii), annual_p = pop(ii);
    Type notification_rate = (annual_v.sum() / (annual_p.sum() / len_dt)) * 1e5;
    dll -= log_normal_lpdf(notification_rate, notification_meanlog[i], notification_sdlog[i]);
  }

  for (int i = 0; i < mortality_meanlog.size(); i++) {
    int idx = asDouble((mortality_year[i] - year_zero) * len_dt);
    Eigen::ArithmeticSequence ii = Eigen::seqN(idx, len_dt);
    vector<Type> annual_v = emr(ii), annual_p = pop(ii);
    Type mortality_rate = (annual_v.sum() / (annual_p.sum() / len_dt)) * 1e5;
    dll -= log_normal_lpdf(mortality_rate, mortality_meanlog[i], mortality_sdlog[i]);
  }

  REPORT(dll);
  REPORT(lhd);  
  REPORT(out0);
  REPORT(out);
  return dll;
}
