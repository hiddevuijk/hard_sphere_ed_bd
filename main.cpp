

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
    Vec3 f(0, 0, -1.0 * U * (r.z - z0) );
    return f;
  }
  double U;
  double z0;
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
  //double bin_size =
  //      params.get_parameter<double>("bin_size");
  //double bulk_density = N / (system_size_x * system_size_y * system_size_z);



  bool pbc_x = true; 
  bool pbc_y = true;
  bool pbc_z = false;
   

   Density density(0, system_size_z, number_of_bins, 'z',
                   system_size_x * system_size_y);

  //PairCorrelation pair_corr(number_of_bins, bin_size, bulk_density, system_size_x, system_size_y, system_size_z);

  //std::vector<Vec3>  positions = initialize_position(N, 1.1,
  //    system_size_x, system_size_y, system_size_z);

  std::vector<Vec3> positions = read_positions("positions_init.dat");

  Potential pot;  
  pot.U = 2.0;
  pot.z0 = system_size_z / 2.0;
  pot.is_nonzero = true;
  
  SystemEDBD<Potential> system(seed, system_size_x, system_size_y, system_size_z,pbc_x, pbc_y,pbc_z, dt, verlet_list_radius, pot, D, gamma);


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
    system.Integrate(time_between_samples);
    density.Sample( system.GetPositions() );
    cout << equilibration_time + number_of_samples * time_between_samples
         << "\t"<<system.GetTime() << endl;
  } 
  
  name = "data/positions.dat";
  system.SavePositions(name);

  //pair_corr.write("data/gr.dat");
  density.Write("rhoz.dat");

  //cout << "\n";
  //cout << pair_corr.GetNumberOfSamplesIgnored() << endl;
  //cout << pair_corr.GetNumberOfSamples() << endl;

  //cout << system.GetNumberOfVerletListUpdates() << endl;
  //cout << system.GetVerletListRadius() << endl;
  //cout << (system.verlet_list_radius_ - 1.0) / 2.0 << endl;
  //cout << system.verlet_list_radius_ << endl;
  return 0;
}
