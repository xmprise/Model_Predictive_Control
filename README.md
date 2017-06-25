# Model Predictive Control

This project is the final task of SDC Term2, which implements the MPC using the location information provided by the simulator and completes the simulation track.

![MPC](https://github.com/xmprise/Model_Predictive_Control/blob/master/pic/img_1.png)

- https://youtu.be/suPJSf6iMJk
- MPC predicted trajectory in green, and reference to the vehicle's coordinate in yellow.
- Tested in Ubuntu 16.04, Intel i7 6700 3.4Ghz

### Model

A simple kinematic model was used. Future improvements can be made using Dynamics models. The state consists of four variables: the x-coordinate, the y-coordinate, the heading direction, and the velocity. Actuator consists of steering angle and acceleration (negative or positive).

**According to wikipedia**
Model predictive controllers rely on dynamic models of the process, most often linear empirical models obtained by system identification. The main advantage of MPC is the fact that it allows the current timeslot to be optimized, while keeping future timeslots in account. This is achieved by optimizing a finite time-horizon, but only implementing the current timeslot. MPC has the ability to anticipate future events and can take control actions accordingly. PID and LQR controllers do not have this predictive ability.

![state and actuations](https://github.com/xmprise/Model_Predictive_Control/blob/master/pic/img_2.png)
### State Variables

#### `px` 
  - The current location in the x-axis in global map coordinate
#### `py`
  - The current location in the y-axis in global map coordinate
#### `psi`
  - The current heading of the vehicle
#### `v`
  - The current velocity of the vehicle


These values are given by the simulator. In addition, provide 'waypoints' that can be used to fit polynomials.
It estimates the curve of the road ahead. The third order polynomial is known to be able to estimate most road curves well. The polynomial function is based on the local coordinates of the vehicle.

### Actuations 

#### `delta`
  - Steering angle represents the angle of turning the vehicle.
  
#### `a`
  - This is the 'throttle' and 'brake' value that indicates the acceleration or deceleration of the vehicle.

### Kinematic Model
Based on physics, there is a simplified version of how the world works. state variable is updated by the elapsed time dt and updated based on the current state and actuators 'delta' and 'a'.

![KinematicModel](https://github.com/xmprise/Model_Predictive_Control/blob/master/pic/img3.png)

### Timestep Length and Elapsed Duration (N & dt)
Generally, it is advantageous to make dt as small as possible and T (N * dt) as large as possible.
To find the appropriate value for N and dt you must set the target speed and keep it constant. The target value for T can then be obtained by fixing dt to 0.1 and increasing the N until the predicted value does not go as far as the horizon and covers most of the roads seen in front of the car. In a few tests, I think T is about 1 second as a reasonable target value. N doubled, dt reduced by half, and vice versa. N = 30, dt = 0.1 made the vibration worse, but N = 20 and dt 0.1 seemed to work better. I then gradually adjusted the parameters and finally got N = 10, dt = 0.1.

### MPC
I refer to the Model Predictive Control quiz to implement MPC. We have added a weight to the cost function to improve the controller. This allows you to individually adjust each part of the cost function to achieve a smooth result. The final weight values are shown below.
```
const double W_CTE = 2000;
const double W_EPSI = 2000;
const double W_V = 1.0;
const double W_DELTA = 5.0;
const double W_A = 5.0;
const double W_DIFF_DELTA = 200.0;
const double W_DIFF_A = 10.0;
```
The cost function is shown below.
```
// The part of the cost based on the reference state.
for (int i = 0; i < N; i++) {
    fg[0] += W_CTE * CppAD::pow(vars[CTE_START + i] - 0, 2);
    fg[0] += W_EPSI * CppAD::pow(vars[EPSI_START + i] - 0, 2);
    fg[0] += CppAD::pow(vars[V_START + i] - MAX_V, 2);
}

// Minimize the use of actuators.
for (int i = 0; i < N -1; i++) {
    fg[0] += W_DELTA * CppAD::pow(vars[DELTA_START + i], 2);
    fg[0] += W_A * CppAD::pow(vars[A_START +i], 2);
}

// Minimize the value gap between sequential actuations.
for (int i = 0; i < N - 2; i++) {
    fg[0] += W_DIFF_DELTA * CppAD::pow(vars[DELTA_START + i + 1] - vars[DELTA_START + i], 2);
    fg[0] += W_DIFF_A * CppAD::pow(vars[A_START + i + 1] - vars[A_START + i], 2);
}
```
The constraints are based on the kinematic model. The formula is shown below.
![KinematicModel](https://github.com/xmprise/Model_Predictive_Control/blob/master/pic/img_5.png)

###Polynomial Fitting and MPC Preprocessing
To estimate the road curve, the waypoints are given in any global coordinate system and must be converted to the vehicle coordinate system.
```
// convert waypoints to vehicle coordinate
for(int  i=0; i < ptsx.size(); ++i){
    // move to origin and rotate around origin to align with axis
    double shift_x = ptsx[i] - px;
    double shift_y = ptsy[i] - py;

    ptsx[i] = shift_x * cos(0-psi) - shift_y * sin(0-psi);
    ptsy[i] = shift_x * sin(0-psi) + shift_y * cos(0-psi);
}
```
I used the third order polynomial as an estimate of the current road curve. It is known that most roads hit well. Using smaller order polynomials can be underfitting, and higher orders are risk of overfitting.

###Model Predictive Control with Latency
I used a kinematic model to predict the position of the car when receiving the command to process the simulated waiting time of 100 ms. To do this, the wait time was converted to seconds and connected to the dt variable. I took this 'delay' into account for all calculations.

```
double x_latency = 0 + v * cos(0) * latency_in_seconds; // x0 + v * cos(psi0) * dt
double y_latency = 0 + v * sin(0) * latency_in_seconds; // y0 + v * sin(psi0) * dt
double psi_latency = 0 - v * steer_value / Lf * latency_in_seconds; // psi0 - v * steer_value / Lf * dt

v += throttle_value * latency_in_seconds; // predicted v after dt
state << x_latency, y_latency, psi_latency, v, cte, epsi;
```
---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.
* Fortran Compiler
  * Mac: `brew install gcc` (might not be required)
  * Linux: `sudo apt-get install gfortran`. Additionall you have also have to install gcc and g++, `sudo apt-get install gcc g++`. Look in [this Dockerfile](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/Dockerfile) for more info.
* [Ipopt](https://projects.coin-or.org/Ipopt)
  * Mac: `brew install ipopt`
  * Linux
    * You will need a version of Ipopt 3.12.1 or higher. The version available through `apt-get` is 3.11.x. If you can get that version to work great but if not there's a script `install_ipopt.sh` that will install Ipopt. You just need to download the source from the Ipopt [releases page](https://www.coin-or.org/download/source/Ipopt/) or the [Github releases](https://github.com/coin-or/Ipopt/releases) page.
    * Then call `install_ipopt.sh` with the source directory as the first argument, ex: `bash install_ipopt.sh Ipopt-3.12.1`. 
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [CppAD](https://www.coin-or.org/CppAD/)
  * Mac: `brew install cppad`
  * Linux `sudo apt-get install cppad` or equivalent.
  * Windows: TODO. If you can use the Linux subsystem and follow the Linux instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).
