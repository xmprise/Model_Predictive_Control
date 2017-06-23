#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

class FG_eval {
public:
    // Coefficients of the fitted polynomial.
    Eigen::VectorXd coeffs;
    FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& vars) {
        // `fg` is a vector containing the cost and constraints.
        // `vars` is a vector containing the variable values (state & actuators).

        fg[0] = 0;

        // The part of the cost based on the reference state.
        for (int i = 0; i < N; i++) {
            fg[0] += 2000 * CppAD::pow(vars[CTE_START + i] - 0, 2);
            fg[0] += 2000 * CppAD::pow(vars[EPSI_START + i] - 0, 2);
            fg[0] += CppAD::pow(vars[V_START + i] - 90, 2);
        }

        // Minimize the use of actuators.
        for (int i = 0; i < N -1; i++) {
            fg[0] += 5 * CppAD::pow(vars[DELTA_START + i], 2);
            fg[0] += 5 * CppAD::pow(vars[A_START +i], 2);
        }

        // Minimize the value gap between sequential actuations.
        for (int i = 0; i < N - 2; i++) {
            fg[0] += 200 * CppAD::pow(vars[DELTA_START + i + 1] - vars[DELTA_START + i], 2);
            fg[0] += 10 * CppAD::pow(vars[A_START + i + 1] - vars[A_START + i], 2);
        }

        // setup constraints

        fg[1 + X_START] = vars[X_START];
        fg[1 + Y_START] = vars[Y_START];
        fg[1 + PSI_START] = vars[PSI_START];
        fg[1 + V_START] = vars[V_START];
        fg[1 + CTE_START] = vars[CTE_START];
        fg[1 + EPSI_START] = vars[EPSI_START];

        for (int i = 0; i < N - 1; i++) {

            // time t + 1
            AD<double> x1 = vars[X_START + i + 1];
            AD<double> y1 = vars[Y_START + i + 1];
            AD<double> psi1 = vars[PSI_START + i + 1];
            AD<double> v1 = vars[V_START + i + 1];
            AD<double> cte1 = vars[CTE_START  + i + 1];
            AD<double> epsi1 = vars[EPSI_START + i + 1];

            // time t
            AD<double> x0 = vars[X_START + i];
            AD<double> y0 = vars[Y_START + i];
            AD<double> psi0 = vars[PSI_START + i];
            AD<double> v0 = vars[V_START + i];
            AD<double> cte0 = vars[CTE_START + i];
            AD<double> epsi0 = vars[EPSI_START + i];

            AD<double> delta0 = vars[DELTA_START + i];
            AD<double> a0 = vars[A_START + i];

            AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
            AD<double> psides0 = CppAD::atan(3*coeffs[3]*x0*x0 + 2*coeffs[2]*x0 + coeffs[1]);

            fg[2 + X_START + i] = x1 - (x0 + v0*CppAD::cos(psi0)*dt);
            fg[2 + Y_START + i] = y1 - (y0 + v0*CppAD::sin(psi0)*dt);
            fg[2 + PSI_START + i] = psi1 - (psi0 - v0*delta0 / Lf*dt);
            fg[2 + V_START + i] = v1 - (v0 + a0*dt);
            fg[2 + CTE_START + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
            fg[2 + EPSI_START + i] = epsi1 - ((psi0 - psides0) - v0*delta0 / Lf*dt);

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

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    // Set the number of model variables (includes both states and inputs).
    // For example: If the state is a 4 element vector, the actuators is a 2
    // element vector and there are 10 timesteps. The number of variables is:
    //
    // 4 * 10 + 2 * 9
    size_t n_vars = N*6 + (N-1)*2;
    // Set the number of constraints
    size_t n_constraints = N*6;

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++) {
        vars[i] = 0;
    }

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    // Set lower and upper limits for variables.

    // set all non-actuators limits
    for (int i = 0; i < DELTA_START; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // set limits for delta
    for (int i = DELTA_START; i < A_START; i++) {
        vars_lowerbound[i] = -0.436332 * Lf; // -degrad(25) * Lf
        vars_upperbound[i] = 0.436332 * Lf; // degrad(25) * lf
    }

    // set limits for acceleration/deceleration
    for (int i = A_START; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }

    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (int i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }

    constraints_lowerbound[X_START] = x;
    constraints_lowerbound[Y_START] = y;
    constraints_lowerbound[PSI_START] = psi;
    constraints_lowerbound[V_START] = v;
    constraints_lowerbound[CTE_START] = cte;
    constraints_lowerbound[EPSI_START] = epsi;

    constraints_upperbound[X_START] = x;
    constraints_upperbound[Y_START] = y;
    constraints_upperbound[PSI_START] = psi;
    constraints_upperbound[V_START] = v;
    constraints_upperbound[CTE_START] = cte;
    constraints_upperbound[EPSI_START] = epsi;

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

    // Return the first actuator values. The variables can be accessed with
    // `solution.x[i]`.
    //
    // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
    // creates a 2 element double vector.

    vector<double> result;

    result.push_back(solution.x[DELTA_START]);
    result.push_back(solution.x[A_START]);

    for (int i=0; i < N-1; i++) {
        result.push_back(solution.x[X_START + i + 1]);
        result.push_back(solution.x[Y_START + i + 1]);
    }

    return result;
}