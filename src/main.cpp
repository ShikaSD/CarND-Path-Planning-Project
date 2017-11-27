#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <iostream>
#include <thread>
#include <vector>
#include "json.hpp"
#include "spline.h"
#include "Trajectory.h"
#include "helpers.h"

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

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(
    double s,
    double d,
    const vector<double> &maps_s,
    const vector<double> &maps_x,
    const vector<double> &maps_y,
    const vector<double> &maps_dx,
    const vector<double> &maps_dy)
{
	int wp = -1;

	while(s > maps_s[wp+1] && (wp < (int)(maps_s.size()-1) ))
	{
		wp++;
	}

  auto prev_wp = wp - 1;
  if (prev_wp < 0) {
    prev_wp += maps_x.size();
  }

  auto wp2 = (wp+1)%maps_x.size();
  vector<double> x_coords(10), y_coords(10);
  for (int i = 0; i < 10; i++) {
      auto index = (prev_wp + i) % maps_x.size();
      x_coords[i] = maps_x[index] + d * maps_dx[index];
      y_coords[i] = maps_y[index] + d * maps_dy[index];
  }

  auto coeffs = polyfit(Eigen::VectorXd(x_coords.data()), Eigen::VectorXd(y_coords.data()), 3);

//  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
//  double seg_s = (s-maps_s[prev_wp]);

  double seg_s = maps_s[wp2] - maps_s[wp];
  if (seg_s < 0) seg_s = 6945.6 - maps_s[wp];
  s -= maps_s[wp];
  double delta_s = s / seg_s;
  double delta_x = (maps_x[wp2] - maps_x[wp]) * delta_s;

  double x = maps_x[wp] + d * maps_dx[wp] + delta_x;
  double y = polyeval(coeffs, x);


//	double heading = atan2((maps_y[wp2]-maps_y[wp]),(maps_x[wp2]-maps_x[wp]));
//	// the x,y,s along the segment
//	double seg_s = (s-maps_s[wp]);
//
//	double seg_x = maps_x[wp]+seg_s*cos(heading);
//	double seg_y = maps_y[wp]+seg_s*sin(heading);
//
//	double x = seg_x + d*maps_dx[wp2];
//	double y = seg_y + d*maps_dy[wp2];

	return {x,y};

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

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](
      uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
      uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

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

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

//          int next;
//          for (next = 0; next < previous_path_x.size(); next++) {
//              double x = previous_path_x[next];
//              double y = previous_path_y[next];
//              if (car_yaw - fabs(atan2(x - (car_x + 0.2 * car_vx), y - (car_y + 0.2 * car_vy))) < pi() / 4) {
//                  break;
//              }
//          }
//
//          next_x_vals.insert(next_x_vals.end(), previous_path_x.begin() + next, previous_path_x.end());
//          next_y_vals.insert(next_y_vals.end(), previous_path_y.begin() + next, previous_path_y.end());


          // int next = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
          // cout << "Next waypoint: " << next << endl;
          for (int i = 1; i <= 100/* - previous_path_x.size() + next*/; ++i) {
              double s = (previous_path_x.size() == 0) ? car_s : end_path_s;

              auto xy = getXY(
                      car_s + i * 0.5,
                      6,
                      map_waypoints_s,
                      map_waypoints_x,
                      map_waypoints_y,
                      map_waypoints_dx,
                      map_waypoints_dy);

              next_x_vals.push_back(xy[0]);
              next_y_vals.push_back(xy[1]);
          }


          cout <<"s = " <<car_s <<endl;

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
