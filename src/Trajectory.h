//
// Created by Andrei Shikov on 24/11/2017.
//

#ifndef PATH_PLANNING_TRAJECTORY_H
#define PATH_PLANNING_TRAJECTORY_H

using namespace std;

class Trajectory {

  public:
    Trajectory(vector<double> &x, vector<double> &y, vector<double> &target, double &t, vector<vector<double>> &predictions);

    vector<double> trajectory_x;
    vector<double> trajectory_y;
};


#endif //PATH_PLANNING_TRAJECTORY_H
