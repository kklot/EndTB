#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

template <class T, class M> 
struct ODE {
  int n_state, n_time;
  state_type x;
  runge_kutta4<state_type> stepper;
  M model;
  vector<double> track;
  ODE(){};
  ODE(vector<T> init, vector<T> pars, double tmax, double dt) :
    n_state(init.size()), 
    n_time(tmax / dt), 
    x(n_state), 
    model(pars),
    track((n_state + 1) * (n_time + 1))
  {
    // copy starting time
    memcpy(&x[0], &init[0], n_state * sizeof(double));
    *(&track(0)) = 0.0;
    memcpy(&track(0) + 1, &x[0], n_state * sizeof(double));
    // running from time
    int mr = 1;
    for (double t = dt; t < tmax; t += dt) {
      stepper.do_step(model, x, t, dt);
      *(&track(0) + mr * (n_state + 1)) = t;
      memcpy(&track(0) + mr * (n_state + 1) + 1, &x[0], n_state * sizeof(double));
      ++mr;
    }
  };
  matrix<double> out () {
    matrix<double> ans = track.matrix();
    ans.resize(n_state + 1, n_time + 1);
    return ans;
  };
};
