//
// Created by Andrei Shikov on 24/11/2017.
//

#include "Trajectory.h"
#include "Eigen-3.3/Core"
#include <random>

int N_SAMPLES = 10;

Trajectory::Trajectory(
        vector<double> &x_coeffs,
        vector<double> &y_coeffs,
        vector<double> &target,
        double &t,
        vector <vector<double>> &predictions) {

    double x = x_coeffs[0];
    double y = y_coeffs[0];
    auto poly = polyfit({x, target[0]}, {y, target[1]}, 3);
    auto steps = t * 5; // 5 steps every second
    auto delta_x = (x_coeffs[0] - target[0]) / steps;
    auto delta_y = (y_coeffs[0] - target[1]) / steps;
    for (auto i = 1; i < steps; i++) {
        this->trajectory_x.push_back(x_coeffs[0] + delta_x * i);
        this->trajectory_y.push_back(y_coeffs[0] + delta_y * i);
    }
}
