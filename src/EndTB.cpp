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
  DATA_VECTOR(nullid); // indices of parameters in null case
  
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
  
  Type prior = 0;

  // transform for sampling
    PARAMETER(beta_s);  // 4.4
    PARAMETER(beta_r); // 12 
    PARAMETER(kappa); // U(0.1 – 10) 
    PARAMETER(b); 
    PARAMETER(mu); 
    PARAMETER(mu_tb); // 0.2
    PARAMETER(theta_s);
    PARAMETER(theta_r);
    PARAMETER(rho);
    PARAMETER(sigma); // 7.9 (95% CrI 3.7-11.8) 
    PARAMETER(delta); // (95% CrI 2.7-11.5)
    PARAMETER(gamma); // 12 (95% CrI 8.5-15) 
    PARAMETER(phi);
    PARAMETER(varepsilon);
    PARAMETER(omega);
    PARAMETER(tau_0);
    PARAMETER(tau_1);
    PARAMETER(chi_s);
    PARAMETER(chi_r);
    PARAMETER(varrho); // 0.05 (95% CrI 0.005-0.09)
    PARAMETER(r_0);
    PARAMETER(r_1);
    PARAMETER(r_2);
    PARAMETER(r_3);
    PARAMETER(varsigma);
    PARAMETER(c_s0);
    PARAMETER(c_r0);
    PARAMETER(c_r1);
    PARAMETER(m_n);
    PARAMETER(m_r);
    PARAMETER(xi);
    vector<Type> pars(32);
    pars << pop1970, beta_s, beta_r, kappa, b, mu, mu_tb, theta_s, theta_r, rho, sigma, delta, gamma, phi, varepsilon, omega, tau_0, tau_1, chi_s, chi_r, varrho, r_0, r_1, r_2, r_3, varsigma, c_s0, c_r0, c_r1, m_n, m_r, xi;

  vector<Type> pars_null = pars; // to flexibly change the index
  for (int i = 0; i < nullid.size(); i++) {
    int idx = asDouble(nullid[i]); // is there a way to cast the index directly
    pars_null[idx] = Type(0); // set zero to health system parameters
  }
  
  pars_null[19] = tmax - year_zero - 1; // move year zero backward for equi. phase

  ODE<Type, TB<Type> > mod0(init, pars_null, asDouble(tmax), dbdt); // run steady state 
  matrix<double> out0 = mod0.out(); // get the equilibrium
  vector<double> eqVec = out0(Eigen::seqN(1, init.size()), Eigen::last);
  vector<Type> eqVecT = Double2Type<Type>(eqVec); 
  eqVecT = eqVecT * (pop1970 / eqVecT.sum()); // adjust to target pop
  
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

  for (int i = 0; i < Treat_qnorm.size(); i++) {
    int idx = asDouble((Treat_year[i] - year_zero) * len_dt);
    vector<Type> 
      annual_success = eTc(Eigen::seqN(idx, len_dt)),
      annual_failed = eTf(Eigen::seqN(idx, len_dt));
    Type 
      prop_f = annual_failed.sum() / (annual_failed.sum() + annual_success.sum()), 
      fail_prop = qnorm(prop_f);
    dll -= dnorm(Treat_qnorm[i], fail_prop, Treat_sd[i], true);
  }
  REPORT(dll);
  REPORT(prior);
  REPORT(lhd);  
  REPORT(out0);
  REPORT(out);
  return dll;
}
