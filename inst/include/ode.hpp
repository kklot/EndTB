#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

template <class T, class M> 
struct ODE {
  int n_state, n_time;
  state_type x;
  runge_kutta4<state_type> stepper; // TODO: let this changable
  M model;
  vector<double> track; // TODO: change from CppAD to std::vector or Eigen's
  ODE(){};
  ODE(vector<T> init, vector<T> pars, double tmax, double dt) :
    n_state(init.size()), // N
    n_time(tmax / dt), // 
    x(n_state), // |state1|state2|...|stateN|
    model(pars),
    track((n_state + 1) * (n_time + 1)) // +1 column to store time steps
  {
    track.setZero();
    // copy starting time
    memcpy(&x[0], &init[0], n_state * sizeof(double));
    *(&track(0)) = 0.0;
    memcpy(&track(0) + 1, &x[0], n_state * sizeof(double));
    // running from time dt
    int mr = 1;
    for (double t = dt; t < tmax; t += dt) 
    {
      stepper.do_step(model, x, t, dt);
      *(&track(0) + mr * (n_state + 1)) = t;
      memcpy(&track(0) + mr * (n_state + 1) + 1, &x[0], n_state * sizeof(double));
      ++mr;
    }
  };
  matrix<double> out () {
    // |time_0|time_1|time_2|...|time_N|
    // |state1|state1|state1|...|state1|
    // ...
    // |stateN|stateN|stateN|...|stateN|
    matrix<double> ans = track.matrix();
    ans.resize(n_state + 1, n_time + 1);
    return ans;
  };
};
