//
// Created by Andrei Shikov on 27/11/2017.
//


#ifndef PATH_PLANNING_HELPERS_H
#define PATH_PLANNING_HELPERS_H

#include "Eigen-3.3/Eigen/Core"

Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);

double polyeval(Eigen::VectorXd coeffs, double x);

#endif //PATH_PLANNING_HELPERS_H
