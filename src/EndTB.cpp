#include <TMB.hpp>
#include <TBM.hpp>

template <class Type>
vector<Type> Double2Type(vector<double> x){
  vector<Type> y(x.size());
  double * px = x.data();
  for(int i=0;i<x.size();i++) y[i] = Type(*(px + i));
  return y;
}

template <class Type>
Type objective_function<Type>::operator()() {
  Type dll = 0.0;
  DATA_VECTOR(init); // initial conds
  DATA_SCALAR(tmax);
  DATA_SCALAR(pop1970);
  DATA_SCALAR(year_zero); // maybe estimate this
  DATA_SCALAR(dt);
  DATA_VECTOR(nullid); // indices of parameters in null case
  DATA_VECTOR(noti); // notification rate / 100.000
  DATA_VECTOR(year); // year of the notification in original format eg... 1999
  DATA_VECTOR(mortality); // mortality rate / 100K from CID
  DATA_VECTOR(mortality_year); // coresponding year - 1990....

  // transform for sampling
    PARAMETER(log_beta_s);  // 4.4
    Type beta_s(exp(log_beta_s));
    dll -= dnorm(log_beta_s, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_beta_s;
    
    PARAMETER(log_beta_r); // 12 
    Type beta_r(exp(log_beta_r));
    dll -= dnorm(log_beta_r, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_beta_r;
    
    PARAMETER(log_kappa); // U(0.1 – 10) 
    Type kappa(exp(log_kappa));
    dll -= dnorm(log_kappa, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_kappa;

    PARAMETER(log_b); 
    Type b(exp(log_b));

    PARAMETER(log_mu); 
    Type mu(exp(log_mu));

    PARAMETER(logit_mu_tb); // 0.2
    Type mu_tb(invlogit(logit_mu_tb));
    dll -= log(mu_tb) +  log(1 - mu_tb);
    dll -= dbeta(mu_tb, Type(2.5), Type(10), true);

    PARAMETER(log_theta_s);
    Type theta_s(exp(log_theta_s));

    PARAMETER(log_theta_r);
    Type theta_r(exp(log_theta_r));

    PARAMETER(log_rho);
    Type rho(exp(log_rho));

    PARAMETER(log_sigma); // 7.9 (95% CrI 3.7-11.8) 
    Type sigma(exp(log_sigma));
    dll -= dnorm(log_sigma, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_sigma;

    PARAMETER(log_delta); // (95% CrI 2.7-11.5)
    Type delta(exp(log_delta));
    dll -= dnorm(log_delta, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_delta;
    
    PARAMETER(log_gamma); // 12 (95% CrI 8.5-15) 
    Type gamma(exp(log_gamma));
    dll -= dnorm(log_gamma, Type(log(5) - .5 * .75 * .75), Type(.75), true) + log_gamma;

    PARAMETER(log_phi);
    Type phi(exp(log_phi));

    PARAMETER(log_varepsilon);
    Type varepsilon(exp(log_varepsilon));

    PARAMETER(log_omega);
    Type omega(exp(log_omega));

    PARAMETER(log_tau_0);
    Type tau_0(exp(log_tau_0));

    PARAMETER(log_tau_1);
    Type tau_1(exp(log_tau_1));

    PARAMETER(log_chi_s);
    Type chi_s(exp(log_chi_s));

    PARAMETER(log_chi_r);
    Type chi_r(exp(log_chi_r));

    PARAMETER(logit_varrho); // 0.05 (95% CrI 0.005-0.09)
    Type varrho(invlogit(logit_varrho));
    dll -= log(varrho) +  log(1 - varrho);
    dll -= dbeta(varrho, Type(0.5), Type(10), true);

    PARAMETER(log_r_0);
    Type r_0(exp(log_r_0));

    PARAMETER(log_r_1);
    Type r_1(exp(log_r_1));

    PARAMETER(log_r_2);
    Type r_2(exp(log_r_2));

    PARAMETER(log_r_3);
    Type r_3(exp(log_r_3));

    PARAMETER(log_varsigma);
    Type varsigma(exp(log_varsigma));

    PARAMETER(log_c_s0);
    Type c_s0(exp(log_c_s0));

    PARAMETER(log_c_r0);
    Type c_r0(exp(log_c_r0));

    PARAMETER(log_c_r1);
    Type c_r1(exp(log_c_r1));

    PARAMETER(log_m_n);
    Type m_n(exp(log_m_n));

    PARAMETER(log_m_r);
    Type m_r(exp(log_m_r));

    PARAMETER(log_xi);
    Type xi(exp(log_xi));

    vector<Type> pars(32);
    pars << pop1970, beta_s, beta_r, kappa, b, mu, mu_tb, theta_s, theta_r, rho, sigma, delta, gamma, phi, varepsilon, omega, tau_0, tau_1, chi_s, chi_r, varrho, r_0, r_1, r_2, r_3, varsigma, c_s0, c_r0, c_r1, m_n, m_r, xi;

    PARAMETER(log_sdlog); // half normal - shrink towards zero
    Type sdlog(exp(log_sdlog));
    dll -= dnorm(sdlog, Type(0), Type(2.5), true) + log_sdlog;

    PARAMETER(log_sdlog_m); // half normal - shrink towards zero
    Type sdlog_m(exp(log_sdlog_m));
    dll -= dnorm(sdlog_m, Type(0), Type(2.5), true) + log_sdlog_m;

  vector<Type> pars_null = pars; // to flexibly change the index
  for (int i = 0; i < nullid.size(); i++) {
    int idx = asDouble(nullid[i]); // is there a way to cast the index directly
    pars_null[idx] = Type(0); // set zero to health system parameters
  }
  ODE<Type, TB<Type> > mod0(init, pars_null, asDouble(tmax), asDouble(dt)); // run steady state 
  matrix<double> out0 = mod0.out(); // get the equilibrium
  vector<double> eqVec = out0(Eigen::seqN(1, init.size()), Eigen::last);
  vector<Type> eqVecT = Double2Type<Type>(eqVec); 
  eqVecT = eqVecT * (pop1970 / eqVecT.sum()); // adjust to target pop
  ODE<Type, TB<Type>> mod(eqVecT, pars, asDouble(2030 - year_zero), 0.01); // rerun with eq as init, first time point represents 1970: maybe not?

  // extract expected data
  matrix<double> out = mod.out(); // get the equilibrium
  vector<Type> ept = Double2Type<Type>(out.row(14+17)); // index of notification
  vector<Type> emr = Double2Type<Type>(out.row(14+18)); // index of notification
  // Expected rate per 100.000
  ept = (ept / pop1970) * 1e5; // constant pop
  emr = (emr / pop1970) * 1e5; // constant pop
  
  // calculate the likelihood
  matrix<Type> ill(noti.size(), 3);
  matrix<Type> mll(mortality.size(), 3);

  for (int i = 0; i < noti.size(); i++) {
    int idx = asDouble((year[i] - year_zero) * 100); // 1/0.01
    ill(i, 0) = log(noti[i]);
    ill(i, 1) = log(ept[i]);
    ill(i, 2) = dnorm(log(noti[i]), log(ept[idx]), sdlog, true);
    dll      -= dnorm(log(noti[i]), log(ept[idx]), sdlog, true);
  }

  // calculate the likelihood of mortality
  for (int i = 0; i < mortality.size(); i++) {
    int idx = asDouble((mortality_year[i] - year_zero) * 100); // 1/0.01
    mll(i, 0) = log(mortality[i]);
    mll(i, 1) = log(emr[i]);
    mll(i, 2) = dnorm(log(mortality[i]), log(emr[idx]), sdlog_m, true);
    dll      -= dnorm(log(mortality[i]), log(emr[idx]), sdlog_m, true);
  }
  REPORT(out0);
  REPORT(out);
  return dll;
}
