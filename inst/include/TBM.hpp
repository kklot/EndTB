#include "ode.hpp"

template<class T>
class TB {
  int j = 14;
  double 
  /*0  */ N,
  /*1  */ beta_s,
  /*2  */ beta_r,
  /*3  */ kappa,
  /*4  */ b,
  /*5  */ mu,
  /*6  */ mu_tb,
  /*7  */ theta_s,
  /*8  */ theta_r,
  /*9  */ rho,
  /*10 */ sigma,
  /*11 */ delta,
  /*12 */ gamma,
  /*13 */ phi,
  /*14 */ varepsilon,
  /*15 */ omega,
  /*16 */ tau_0,
  /*17 */ tau_1,
  /*18 */ chi_s,
  /*19 */ chi_r,
  /*20 */ varrho,
  /*21 */ r_0,
  /*22 */ r_1,
  /*23 */ r_2,
  /*24 */ r_3,
  /*25 */ varsigma,
  /*26 */ c_s0,
  /*27 */ c_r0,
  /*28 */ c_r1,
  /*29 */ m_n,
  /*30 */ m_r,
  /*31 */ xi;
public:
  TB(vector<T> pars) :
    N          (asDouble(pars[0])),
    beta_s     (asDouble(pars[1])),
    beta_r     (asDouble(pars[2])),
    kappa      (asDouble(pars[3])),
    b          (asDouble(pars[4])),
    mu         (asDouble(pars[5])),
    mu_tb      (asDouble(pars[6])),
    theta_s    (asDouble(pars[7])),
    theta_r    (asDouble(pars[8])),
    rho        (asDouble(pars[9])),
    sigma      (asDouble(pars[10])),
    delta      (asDouble(pars[11])),
    gamma      (asDouble(pars[12])),
    phi        (asDouble(pars[13])),
    varepsilon (asDouble(pars[14])),
    omega      (asDouble(pars[15])),
    tau_0      (asDouble(pars[16])),
    tau_1      (asDouble(pars[17])),
    chi_s      (asDouble(pars[18])),
    chi_r      (asDouble(pars[19])),
    varrho     (asDouble(pars[20])),
    r_0        (asDouble(pars[21])),
    r_1        (asDouble(pars[22])),
    r_2        (asDouble(pars[23])),
    r_3        (asDouble(pars[24])),
    varsigma   (asDouble(pars[25])),
    c_s0       (asDouble(pars[26])),
    c_r0       (asDouble(pars[27])),
    c_r1       (asDouble(pars[28])),
    m_n        (asDouble(pars[29])),
    m_r        (asDouble(pars[30])),
    xi         (asDouble(pars[31]))
  {}
  void operator()(const state_type &x, state_type &dxdt, const double /* t */) {
    // FOI
    double 
    Nt = std::accumulate(x.begin(), x.end(), 0),
    lambda_s = beta_s / Nt * (x[2] + x[3] + x[7] + x[8] + kappa * (x[4] + x[5] + x[9] + x[10])),
    lambda_r = beta_r / Nt * (x[j+2] + x[j+3] + x[j+7] + x[j+8] + kappa * (x[j+4] + x[j+5] + x[j+9] + x[j+10]));
    // DEs
    dxdt[0]  = b * N - (lambda_r + lambda_s) * x[0] - mu * x[0];
    // Sensitive
    dxdt[1]  = (1 - theta_s) * lambda_s * x[0] - rho * x[1] - mu * x[1];
    dxdt[2]  = theta_s * lambda_s * x[0] + rho * x[1] - sigma * x[2] - mu * x[2] - mu_tb * x[2];
    dxdt[3]  = sigma * x[2] - delta * x[3] - mu * x[3] - mu_tb * x[3];
    dxdt[4]  = delta * x[3] + gamma * x[5] - phi * x[4] - mu * x[4] - mu_tb * x[4];
    dxdt[5]  = (1 - varepsilon * omega) * phi * x[4] - gamma * x[5] - mu * x[5] - mu_tb * x[5];
    dxdt[6]  = varepsilon * omega  * phi * x[4] - (tau_0 + chi_s + varrho) * x[6] - mu * x[6];
    dxdt[7]  = x[12] * r_0 + x[13] * r_1 + x[14] * r_2 - sigma * x[7] - mu * x[7] - mu_tb * x[7];
    dxdt[8]  = sigma * x[7] - delta * x[8] - mu * x[8] - mu_tb * x[8];
    dxdt[9]  = delta * x[8] + gamma * x[10] - phi * x[9] - mu * x[9] - mu_tb * x[9];
    dxdt[10] = (1 - varepsilon * omega) * phi * x[9] - gamma * x[10] - mu * x[10] - mu_tb * x[10];
    dxdt[11] = varepsilon * omega  * phi * x[9] - (tau_0 + chi_s + varrho) * x[11] - mu * x[11];
    dxdt[12] = varsigma * (x[13] + x[14]) - r_0 * x[12] - mu * x[12];
    dxdt[13] =          tau_0 *      c_s0   * (x[6] + x[11]) - varsigma * x[13] - r_1 * x[13] - mu * x[13];
    dxdt[14] = (chi_s + tau_0 * (1 - c_s0)) * (x[6] + x[11]) - varsigma * x[14] - r_2 * x[14] - mu * x[14];
    // MDR
    dxdt[j+1]  = (1 - theta_r) * lambda_r * x[0] - rho * x[j+1] - mu * x[j+1];
    dxdt[j+2]  = theta_r * lambda_r * x[0] + rho * x[j+1] - sigma * x[j+2] - mu * x[j+2] - mu_tb * x[j+2];
    dxdt[j+3]  = sigma * x[j+2] - delta * x[j+3] - mu * x[j+3] -  mu_tb * x[j+3];
    dxdt[j+4]  = delta * x[j+3] + gamma * x[j+5] - phi * x[j+4] - mu * x[j+4] - mu_tb * x[j+4];
    dxdt[j+5]  = (1 - varepsilon * omega) * phi * x[j+4]
               + tau_0 * (1-c_r0) * (1 - xi) * x[j+6] 
               + tau_1 * (1-c_r1) * x[j+15] + chi_r * x[j+15]
               - gamma * x[j+5] - mu * x[j+5] - mu_tb * x[j+5];
    dxdt[j+6]  = varepsilon * omega  * phi * (1 - m_n) * x[j+4] + varrho * x[6] - tau_0 * x[j+6] - mu * x[j+6];
    dxdt[j+7]  = x[j+12] * r_0 + x[j+13] * r_1 + x[j+14] * r_3 - sigma * x[j+7] - mu * x[j+7] - mu_tb * x[j+7];
    dxdt[j+8]  = sigma * x[j+7] - delta * x[j+8] - mu * x[j+8] - mu_tb * x[j+8];
    dxdt[j+9]  = delta * x[j+8] + gamma * x[j+10] - phi * x[j+9] - mu * x[j+9] - mu_tb * x[j+9];
    dxdt[j+10] = (1 - varepsilon * omega) * phi * x[j+9] 
               + tau_0 * (1 - c_r0) * (1 - xi) * x[j+11] 
               + tau_1 * (1 - c_r1) * x[j+16] + chi_r * x[j+16]
               - gamma * x[j+10] - mu * x[j+10] - mu_tb * x[j+10];
    dxdt[j+11] = varepsilon * omega  * phi * (1 - m_r) * x[j+9] + varrho * x[11] - tau_0 * x[j+11] - mu * x[j+11];
    dxdt[j+12] = varsigma * x[j+13] - r_0 * x[j+12] - mu * x[j+12];
    dxdt[j+13] = tau_1 * c_r1 * (x[j+15] + x[j+16]) - varsigma * x[j+13] - r_1 * x[j+13] - mu * x[j+13];
    dxdt[j+14] = tau_0 * c_r0 * (x[j+6] + x[j+11]) - r_3 * x[j+14] - mu * x[j+14];
    dxdt[j+15] = phi * varepsilon * omega * m_n * x[j+4]
               + tau_0 * (1 - c_r0) * xi * x[j+6]
               - (tau_1 + chi_r) * x[j+15] - mu * x[j+15];
    dxdt[j+16] = phi * varepsilon * omega * m_r * x[j+9]
               + tau_0 * (1 - c_r0) * xi * x[j+11]
               - (tau_1 + chi_r) * x[j+16] - mu * x[j+16];
  }
};
