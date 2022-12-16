

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include <iostream>
#include <vector>
#include <string>


using namespace std;

class Potential {
 public:
  Vec3 Force(Vec3 r, double t) {
    Vec3 f(0, 0, 0);
    return f;
  }
};


int main()
{

  Config params("input.txt");
  
  double D = params.get_parameter<double>("D");
  double gamma = params.get_parameter<double>("gamma");

  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  double system_size_x =
		params.get_parameter<double>("system_size_x");
  double system_size_y =
		params.get_parameter<double>("system_size_y");
  double system_size_z =
		params.get_parameter<double>("system_size_z");

  unsigned int N =
		params.get_parameter<double>("N");

  double dt = params.get_parameter<double>("dt");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");

  double equilibration_time = 
         params.get_parameter<double>("equilibration_time");
  double time_between_samples = 
         params.get_parameter<double>("time_between_samples");
  unsigned int number_of_samples = 
         params.get_parameter<unsigned int>("number_of_samples");

  std::vector<Vec3>  positions = initialize_position(N, 3,
      system_size_x,system_size_y, system_size_z);

  Potential pot;  

  
  SystemEDBD<Potential> system(seed, system_size_x, system_size_y, system_size_z, dt, verlet_list_radius, pot, D, gamma);


  system.SetPositions(positions);

  system.Integrate(equilibration_time);


  for (unsigned int i = 0; i < number_of_samples; ++i) {
    string name = "data/positions_" + to_string(i) + ".dat";
    system.SavePositions(name);
    system.Integrate(time_between_samples);
  } 
 
  string name = "data/positions.dat";
  system.SavePositions(name);



  return 0;
}
