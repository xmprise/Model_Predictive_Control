#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

const size_t N = 10; // Number of timestamp in the horizon
const double dt = 0.1; // How much time elapses between actuations

const double Lf = 2.67; // Length from front of vehicle to Center-of-Gravity

const size_t X_START = 0;
const size_t Y_START = X_START + N;
const size_t PSI_START = Y_START + N;
const size_t V_START = PSI_START + N;
const size_t CTE_START = V_START + N;
const size_t EPSI_START = CTE_START + N;
const size_t DELTA_START = EPSI_START + N;
const size_t A_START = DELTA_START + N - 1;
const double MAX_V = 90;

// Weights for cost function
const double W_CTE = 2000;
const double W_EPSI = 2000;
const double W_V = 1.0;
const double W_DELTA = 5.0;
const double W_A = 5.0;
const double W_DIFF_DELTA = 200.0;
const double W_DIFF_A = 10.0;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
