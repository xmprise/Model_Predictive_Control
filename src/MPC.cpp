#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
      fg[0] = 0;

      // The part of the cost based on the reference state.
      for (int i = 0; i < N; i++) {
          fg[0] += W_CTE * CppAD::pow(vars[CTE_START + i], 2);
          fg[0] += W_EPSI * CppAD::pow(vars[EPSI_START + i], 2);
          fg[0] += W_V * CppAD::pow(vars[V_START + i] - MAX_V, 2);
      }

      // Minimize the use of actuators.
      for (int i = 0; i < N - 1; i++) {
          fg[0] += W_DELTA * CppAD::pow(vars[DELTA_START + i], 2);
          fg[0] += W_A * CppAD::pow(vars[A_START + i], 2);
      }

      // Minimize the value gap between sequential actuations.
      for (int i = 0; i < N - 2; i++) {
          fg[0] += W_DIFF_DELTA * CppAD::pow(vars[DELTA_START + i + 1] - vars[DELTA_START + i], 2);
          fg[0] += W_DIFF_A * CppAD::pow(vars[A_START + i + 1] - vars[A_START + i], 2);
      }

      // Initial constraints
      // We initialize the model to the initial state.
      // Recall fg[0] is reserved for the cost value, so the other indices are bumped up by 1.
      fg[1 + X_START] = vars[X_START];
      fg[1 + Y_START] = vars[Y_START];
      fg[1 + PSI_START] = vars[PSI_START];
      fg[1 + V_START] = vars[V_START];
      fg[1 + CTE_START] = vars[CTE_START];
      fg[1 + EPSI_START] = vars[EPSI_START];

      // Constraints based on the vehicle mode
      for(int i = 0; i < N - 1; i++) {
          // The state at time t
          AD<double> x0 = vars[X_START + i];
          AD<double> y0 = vars[Y_START + i];
          AD<double> psi0 = vars[PSI_START + i];
          AD<double> v0 = vars[V_START + i];
          AD<double> cte0 = vars[CTE_START + i];
          AD<double> epsi0 = vars[EPSI_START + i];

          // The state at time t+1 .
          AD<double> x1 = vars[X_START + i + 1];
          AD<double> y1 = vars[Y_START + i + 1];
          AD<double> psi1 = vars[PSI_START + i + 1];
          AD<double> v1 = vars[V_START + i + 1];
          AD<double> cte1 = vars[CTE_START + i + 1];
          AD<double> epsi1 = vars[EPSI_START + i + 1];

          // Only consider the actuation at time t.
          AD<double> delta0 = vars[DELTA_START + i];
          AD<double> a0 = vars[A_START + i];

          AD<double> f0 = coeffs[0] + coeffs[1] * x0;
          AD<double> psides0 = CppAD::atan(coeffs[1]);

          fg[2 + X_START + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
          fg[2 + Y_START + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
          fg[2 + PSI_START + i] = psi1 - (psi0 + v0 * (-delta0 / Lf) * dt);
          fg[2 + V_START + i] = v1 - (v0 + a0 * dt);
          fg[2 + CTE_START + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
          fg[2 + EPSI_START + i] = epsi1 - ((psi0 - psides0) + v0 * (-delta0 / Lf) * dt);

      }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 0;
  // TODO: Set the number of constraints
  size_t n_constraints = 0;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  return {};
}
