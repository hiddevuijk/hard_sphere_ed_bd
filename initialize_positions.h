#ifndef GUARD_INITIALIZE_POSITIONS_H
#define GUARD_INITIALIZE_POSITIONS_H


#include "vec3.h"
#include <vector>


std::vector<Vec3> initialize_position(
      unsigned int N, double d,
      double Lx, double Ly, double Lz)
{
  std::vector<Vec3> positions(N);
  std::vector<Vec3> lattice;
  Vec3 temp; 

  double zmid = Lz / 2.0;
  for (unsigned int zi = 0; zi < floor(Lz/d) / 2; ++zi) {
  for (unsigned int yi = 0; yi < floor(Ly/d); ++yi) {
  for (unsigned int xi = 0; xi < floor(Lx/d); ++xi) {

    temp.x = xi * d; 
    temp.y = yi * d; 
    temp.z = zmid + zi * d; 
    lattice.push_back(temp);
    if (zi > 0) {
      temp.z = zmid - zi * d; 
      lattice.push_back(temp);
    }
  }}}

  for (unsigned int i = 0; i < N; ++i) {
    positions[i] = lattice[i];
  } 

  return positions;
}



#endif
