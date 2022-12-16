
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
  double system_size_x = -1;
  double system_size_y = -1;
  double system_size_z = -1;
  double dt = 1e-4;
  double verlet_list_radius = 3.0;

  Potential pot;  

  
  SystemEDBD<Potential> system(seed, system_size_x, system_size_y, system_size_z, dt, verlet_list_radius, pot);


  Vec3 r1(0, 0, 0);
  Vec3 r2(3, 0, 0);

  Vec3 v1, v2;
 
  vector<Vec3> pos;
  pos.push_back(r1);
  pos.push_back(r2);



  system.SetPositions(pos);

  system.SetV(0,  0.0, 0.0, 0.0);
  system.SetV(1, -1.0, 0.0, 0.0);
  
  system.MoveBallistically(1.0);

  r1 = system.GetPosition(0);
  r2 = system.GetPosition(1);
  v1 = system.GetVelocity(0);
  v2 = system.GetVelocity(1);

  cout << "at collision: \n";
  cout << "r: \n";

  cout << "\t" <<  r1 << endl;
  cout << "\t" <<  r2 << endl;

  cout << "v: \n";
  cout <<  "\t" << v1 << endl;
  cout <<  "\t" << v2 << endl;

  system.MakeCollision(0,1);

  r1 = system.GetPosition(0);
  r2 = system.GetPosition(1);
  v1 = system.GetVelocity(0);
  v2 = system.GetVelocity(1);

  cout << "after collision: \n";
  cout << "r: \n";

  cout << "\t" <<  r1 << endl;
  cout << "\t" <<  r2 << endl;

  cout << "v: \n";
  cout <<  "\t" << v1 << endl;
  cout <<  "\t" << v2 << endl;


  return 0;
}
