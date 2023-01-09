

#include "config_file.h"
#include "vec3.h"
#include "systemEDBD.h"
#include "initialize_positions.h"

#include "pair_correlation.h"
#include "density.h"
#include "read_positions.h"

#include <iostream>
#include <vector>
#include <string>


using namespace std;

class Potential {
 public:
  Vec3 Force(Vec3 r, double t) {
    Vec3 f(0, 0, -r.z);
    return f;
  }

  bool is_nonzero;
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


  unsigned int number_of_bins =
        params.get_parameter<unsigned int>("number_of_bins");
  double bin_size =
        params.get_parameter<double>("bin_size");
  double bulk_density = N / (system_size_x * system_size_y * system_size_z);


   double zmax = params.get_parameter<double>("zmax");
   double zmin = params.get_parameter<double>("zmin");

   Density density(zmin, zmax, number_of_bins, 'z',
                   system_size_x * system_size_y);

  PairCorrelation pair_corr(number_of_bins, bin_size, bulk_density, system_size_x, system_size_y, system_size_z);

  //std::vector<Vec3>  positions = initialize_position(N, 1.1,
  //    0, system_size_x,0, system_size_y, zmin, zmax);

  std::vector<Vec3> positions = read_positions("positions.dat");

  Potential pot;  
  pot.is_nonzero = true;

  
  SystemEDBD<Potential> system(seed, system_size_x, system_size_y, system_size_z, dt, verlet_list_radius, pot, D, gamma);


  system.SetPositions(positions);

  while (system.GetTime() < equilibration_time) {
    system.Integrate(time_between_samples);
    cout << equilibration_time << "\t"<<system.GetTime() << endl;
  }


  string name; 
  for (unsigned int i = 0; i < number_of_samples; ++i) {
    //string name = "data/positions_" + to_string(i) + ".dat";
    //system.SavePositions(name);
    //pair_corr.sample( system.GetPositions() );
    density.Sample( system.GetPositions() );
    system.Integrate(time_between_samples);
    cout << equilibration_time + number_of_samples * time_between_samples
         << "\t"<<system.GetTime() << endl;
  } 
  
  name = "data/positions.dat";
  system.SavePositions(name);

  //pair_corr.write("data/gr.dat");
  density.Write("rhoz.dat");

  cout << "\n";
  cout << pair_corr.GetNumberOfSamplesIgnored() << endl;
  cout << pair_corr.GetNumberOfSamples() << endl;

  //cout << system.GetNumberOfVerletListUpdates() << endl;
  //cout << system.GetVerletListRadius() << endl;
  //cout << (system.verlet_list_radius_ - 1.0) / 2.0 << endl;
  //cout << system.verlet_list_radius_ << endl;
  return 0;
}
