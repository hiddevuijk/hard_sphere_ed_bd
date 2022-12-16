#ifndef GUARD_INITIALIZE_POSITIONS_H
#define GUARD_INITIALIZE_POSITIONS_H


#include "vec3.h"
#include <vector>


std::vector<Vec3> initialize_position(
      unsigned int N, double d, double Lx, double Ly, double Lz)
{
  std::vector<Vec3> positions(N);
  std::vector<Vec3> lattice;
  int nx = 0;
  int ny = 0;
  Vec3 temp; 
  while (nx * d < Lx ) {
    while (ny * d < Ly ) {
      temp.x = d * nx; 
      temp.y = d * ny; 
      lattice.push_back(temp);
      ny += 1;    
    }
    nx += 1;
  }

  for (unsigned int i = 0; i < N; ++i) {
    positions[i] = lattice[i];
  } 

  return positions;
}



#endif
