
#include "systemEDBD.h"

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
  unsigned long int seed = 1234567890;
  double system_size_x = 10.0;
  double system_size_y = 10.0;
  double system_size_z = 10.0;
  double dt = 1e-4;
  double verlet_list_radius = 3.0;

  Potential pot;  

  SystemEDBD<Potential> system(seed, system_size_x, system_size_y, system_size_z, dt, verlet_list_radius, pot);


  return 0;
}
