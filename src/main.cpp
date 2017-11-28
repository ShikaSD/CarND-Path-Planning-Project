#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <ctime>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "MPC.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint = (closestWaypoint + 1) % maps_x.size();
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

int car_ahead(nlohmann::json &sensor_fusion, double lane, double car_x, double car_y, double car_yaw, double look_ahead, bool look_back) {
  double current_dist = 1e+10;
  int car_index = -1;

  for (int i = 0; i < sensor_fusion.size(); i++) {
    auto car = sensor_fusion[i];
    double x = car[1];
    double y = car[2];
    double vx = car[3];
    double vy = car[4];
    double s = car[5];
    double d = car[6];

    double limit = look_back ? pi() * .75 : pi() * .5;

    if (d > 4 * lane && d < 4 * (lane + 1)) {
      // We are on the same lane
      double dist = distance(x, y, car_x, car_y);
      current_dist = min(dist, current_dist);
      double angle = atan2(y - car_y, x - car_x) - car_yaw;
      if ((dist < look_ahead) && (fabs(atan2(sin(angle), cos(angle))) < limit)) {
        if (fabs(current_dist - dist) < 0.1) {
          car_index = i;
        }
      }
    }
  }

  return car_index;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    double s;
    double d_x;
    double d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  int target_lane = 1;
  double speed_limit = 49.9;
  double inc_update = (speed_limit / 2.24) * 0.02 / 2;
  double ref_speed = 0;
  double min_speed = inc_update;
  MPC mpc;

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &target_lane, &speed_limit, &inc_update, &ref_speed, &min_speed, &mpc](
      uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
      uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto q1 = hasData(data);

      if (q1 != "") {
        auto j = json::parse(q1);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x     = j[1]["x"];
          double car_y     = j[1]["y"];
          double car_s     = j[1]["s"];
          double car_d     = j[1]["d"];
          double car_yaw   = deg2rad(j[1]["yaw"]);
          double car_speed = j[1]["speed"];
          double car_vx    = car_speed * cos(car_yaw);
          double car_vy    = car_speed * sin(car_yaw);

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          double buf_size = previous_path_x.size() / 2;

          int current_lane = ((int) car_d) / 4 ;

          int empty_limit = 2;
          vector<vector<double>> wpts;
          double ref_x = car_x, ref_y = car_y, ref_yaw = car_yaw;

          if (buf_size < empty_limit) {
            double prev_x = ref_x - cos(car_yaw);
            double prev_y = ref_y - sin(car_yaw);
            wpts.push_back({ prev_x, prev_y });
            wpts.push_back({ ref_x, ref_y });
          } else {
            ref_x = previous_path_x[buf_size - 1];
            ref_y = previous_path_y[buf_size - 1];

            double prev_x = (previous_path_x[buf_size - 2]);
            double prev_y = (previous_path_y[buf_size - 2]);

            ref_yaw = atan2(ref_y - prev_y, ref_x - prev_x);
            wpts.push_back({ prev_x, prev_y });
            wpts.push_back({ ref_x, ref_y });
          }

          double look_ahead = 60;

          int car_ahead_index = car_ahead(sensor_fusion, current_lane, car_x, car_y, car_yaw, look_ahead / 2, false);
          bool is_car_ahead = car_ahead_index > -1;

          if (is_car_ahead) {
            ref_speed -= inc_update;

            auto car_ahead_ = sensor_fusion[car_ahead_index];
            double car_ahead_vx = car_ahead_[3];
            double car_ahead_vy = car_ahead_[4];
            double car_ahead_speed = sqrt(car_ahead_vx * car_ahead_vx + car_ahead_vy * car_ahead_vy);
            double car_dist = distance(car_x, car_y, car_ahead_[1], car_ahead_[2]);
            // Speed update
            min_speed = max(inc_update, car_ahead_speed * 2.24 - (look_ahead / 2 - car_dist));

            // Choose lane if not switching
            if (target_lane == current_lane) {
              // TODO: Implement for current lane a well
              vector<vector<double>> lanes = { { (double) current_lane + 1, -1 }, { (double) current_lane - 1, -1 } };

              int next_lane = current_lane;
              for (auto &pair: lanes) {
                int lane = (int) pair[0];
                if (lane < 0 || lane > 2) {
                  pair[1] = 1e10;
                  continue;
                }

                // TODO: Use both
                auto lane_car_index = car_ahead(sensor_fusion, lane, car_x, car_y, car_yaw, look_ahead, false);
                auto lane_back_index = car_ahead(sensor_fusion, lane, car_x, car_y, car_yaw, look_ahead, true);

//                if (lane_car_index != lane_back_index) lane_car_index = lane_back_index;

                if (lane_car_index != -1) {
                  auto lane_ahead_ = sensor_fusion[lane_car_index];
                  double lane_ahead_vx = lane_ahead_[3];
                  double lane_ahead_vy = lane_ahead_[4];
                  double lane_ahead_speed = sqrt(lane_ahead_vx * lane_ahead_vx + lane_ahead_vy * lane_ahead_vy);
                  double lane_ahead_dist = distance(car_x, car_y, lane_ahead_[1], lane_ahead_[2]);

                  if (lane_ahead_dist > car_dist && car_ahead_speed < lane_ahead_speed) {
                    target_lane = lane;
                    pair[1] = 1000 / (lane_ahead_dist + lane_ahead_speed);
                  } else {
                    pair[1] = 1e10;
                  }
                }
              }
              sort(lanes.begin(), lanes.end(), [&](vector<double> &a, vector<double> &b) { return a[1] < b[1]; });
              if (lanes[0][1] < 1e10) {
                target_lane = (int) lanes[0][0];
              }
            }

          } else if (ref_speed < speed_limit) {
            ref_speed += inc_update;
          }

          bool lane_change = target_lane != current_lane;

          cout <<"car index: " <<car_ahead_index <<" min_speed: " << min_speed <<endl;

          ref_speed = max(min_speed, min(ref_speed, speed_limit));

          double m_per_update = 0.02 * ref_speed / 2.24; // freq * speed_limit / mph to mps

          double dist = lane_change ? 40 : 30;
          for (int i = 0; i < 2; ++i) {
            wpts.push_back(getXY(car_s + (i + 1) * dist, (2 + 4 * target_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y));
          }

          Eigen::VectorXd x_coords(wpts.size());
          Eigen::VectorXd y_coords(wpts.size());
          // Transform to car coords
          for (int i = 0; i < wpts.size(); ++i) {
            auto wpt = wpts[i];
            double shift_x = wpt[0] - ref_x;
            double shift_y = wpt[1] - ref_y;
            x_coords[i] = (shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw));
            y_coords[i] = (shift_y * cos(-ref_yaw) + shift_x * sin(-ref_yaw));
          }

          auto coeffs = polyfit(x_coords, y_coords, 3);

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (int i = 0; i < buf_size; ++i) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = look_ahead;
          double target_y = polyeval(coeffs, target_x);
          double target_dist = distance(0, 0, target_x, target_y);
          double steps = target_dist / m_per_update;
          double step = target_x / steps;

          for (int i = 1; i <= 50 - buf_size - 1; ++i) {
            double x = step * i;
            double y = polyeval(coeffs, x);

            double x_val = (x * cos(ref_yaw) - y * sin(ref_yaw)) + ref_x;
            double y_val = (y * cos(ref_yaw) + x * sin(ref_yaw)) + ref_y;

            next_x_vals.push_back(x_val);
            next_y_vals.push_back(y_val);
          }

//          Eigen::VectorXd x0(6);
//          double cte = polyeval(coeffs[0], 0);
//          x0 <<0, 0, 0, car_speed, polyeval(coeffs[0], 0), -atan(coeffs[1]);
//          mpc.Solve(x0, coeffs, speed_limit / 2.24);
//
//          for (int i = 0; i < mpc.ptsx.size(); ++i) {
//            double x = mpc.ptsx[i];
//            double y = mpc.ptsy[i];
//
//            double x_val = (x * cos(ref_yaw) - y * sin(ref_yaw)) + ref_x;
//            double y_val = (y * cos(ref_yaw) + x * sin(ref_yaw)) + ref_y;
//
//            next_x_vals.push_back(x_val);
//            next_y_vals.push_back(y_val);
//          }

          json msgJson;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\"," + msgJson.dump() + "]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
