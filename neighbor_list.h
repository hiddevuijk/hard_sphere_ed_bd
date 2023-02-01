#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H


#include <vector>
#include <list>
#include "vec3.h"

// returns the distance between the two vectors a and b
// taking into account the periodic boundary conditions
// If double_bond == true, then each neighbor pair occurs
//   twice in the neighbor list. If false, each bond 
//   occurs twice.
double dist2(Vec3 a, Vec3 b, double Lx, double Ly, double Lz,
              bool pbc_x, bool pbc_y, bool pbc_z)
{
  a-= b;
  //if (pbc_x) a.x -= Lx * std::round( a.x / Lx);
  //if (pbc_y) a.y -= Ly * std::round( a.y / Ly);
  //if (pbc_z) a.z -= Lz * std::round( a.z / Lz);
  if (pbc_x) {
    a.x = fabs(a.x);
    a.x -= static_cast<int>(a.x / Lx +  0.5) * Lx;
  }
  if (pbc_y) {
    a.y = fabs(a.y);
    a.y -= static_cast<int>(a.y / Ly +  0.5) * Ly;
  }
  if (pbc_z) {
    a.z = fabs(a.z);
    a.z -= static_cast<int>(a.z / Lz +  0.5) * Lz;
  }

  return a.LengthSquared();
}

// Contains the three indices that specify the cell
// that a particle is in
struct TripleIndex {
  TripleIndex() : xi(0), yi(0), zi(0) {};

  TripleIndex(unsigned int xi, unsigned int yi, unsigned int zi)
    : xi(xi), yi(yi), zi(zi) {};

  unsigned int xi, yi, zi;
};

// a vector of vectors of vectors of lists
typedef std::vector<std::vector<std::vector<std::list<unsigned int> > > >
        Vec3List;
// a vectors of vectors of lists
typedef std::vector<std::vector<std::list<unsigned int> > > Vec2List;

//  vector of lists
typedef std::vector<std::list<unsigned int> > VecList;

// returns a Verlet neighbor list
// where each neighbor pair is saved once.
// If i and j are neighbors with i < j, then
//  the i'th element of the neighbor list is a list
//  that contains j, but  the j'th element of the neighbor list
//  does not contain i.
// The inputs are
//   Lx, Ly, Lz: the x,y,z dimensions of the box,
//    with (x=0, y=0,z=0) in a corner of the box.
//   r_verlet: radius of the neighbor list
//   positions: a vector with the positions of the particles 
//
VecList get_neighbor_list(
    double Lx, double Ly, double Lz,
    bool pbc_x, bool pbc_y, bool pbc_z,
    bool double_bond,
    double r_verlet,
    const std::vector<Vec3>& positions)
{

  unsigned int n_particles = positions.size();

  // for the cell dimension, the smalles size is taken
  // that is larger than r_verlet, and fits and integer
  // number of times in the x direction
  // 
  // number of cells along the x directon
  unsigned int n_cells_x = std::floor(Lx / r_verlet);
  // cell size in the x direction
  double delta_x = Lx / n_cells_x;

  // same for the y and z direction
  unsigned int n_cells_y = std::floor(Ly / r_verlet);
  unsigned int n_cells_z = std::floor(Lz / r_verlet);
  double delta_y = Ly / n_cells_y;
  double delta_z = Lz / n_cells_z;

  // list of all neighbor pairs, 
  VecList neighbor_list(n_particles);

  // a three dimensional array containing lists with the 
  // particle indices of the particles in that cell
  //Vec3List cell_list(n_cells_z,
  //          Vec2List(n_cells_y, VecList(n_cells_x) ) );
  Vec3List cell_list(n_cells_x,
            Vec2List(n_cells_y, VecList(n_cells_z) ) );

  // the i'th element contains the cell index
  std::vector<TripleIndex> particles_cell_index(n_particles);
  // generate cell list
  Vec3 p; // position of particle pi
  TripleIndex ci; // cell index of particle  pi
  for (unsigned int pi = 0; pi < n_particles; ++pi) {
    p = positions[pi]; 

    if (pbc_x) p.x -= Lx * floor( p.x / Lx);
    if (pbc_y) p.y -= Ly * floor( p.y / Ly);
    if (pbc_z) p.z -= Lz * floor( p.z / Lz);

    //ci.xi = floor(p.x / delta_x);
    //ci.yi = floor(p.y / delta_y);
    //ci.zi = floor(p.z / delta_z);

    ci.xi = static_cast<unsigned int>(p.x / delta_x);
    ci.yi = static_cast<unsigned int>(p.y / delta_y);
    ci.zi = static_cast<unsigned int>(p.z / delta_z);

    // add particle pi to cell ci
    cell_list[ci.xi][ci.yi][ci.zi].push_back(pi);
    // save cell index ci of particle pi
    particles_cell_index[pi] = ci;

  } // end loop over particles

  // generate neighbor_list from cell list
 
  // index of particle in neighboring cell
  // of particle pi
  unsigned int pj; 

  // loop over particles 
  TripleIndex cj;
  for (unsigned int pi = 0; pi < n_particles; ++pi) {

    // cell index of pi
    ci = particles_cell_index[pi];
    //loop over cells that are neighbors of ci
    for (int dcxi = -1; dcxi <= 1; ++dcxi){
      if  (ci.xi == 0 and dcxi == -1) {
        if (pbc_x) cj.xi = n_cells_x -1;
        else continue;
      } else if(ci.xi == n_cells_x -1 and dcxi == 1) {
        if (pbc_x) cj.xi = 0;
        else continue;
      } else {
        cj.xi = ci.xi + dcxi;
      }

    for (int dcyi = -1; dcyi <= 1; ++dcyi){
      if  (ci.yi == 0 and dcyi == -1) {
        if (pbc_y) cj.yi = n_cells_y -1;
        else continue;
      } else if(ci.yi == n_cells_y -1 and dcyi == 1) {
        if (pbc_y) cj.yi = 0;
        else continue;
      } else {
        cj.yi = ci.yi + dcyi;
      }

    for (int dczi = -1; dczi <= 1; ++dczi){

      if  (ci.zi == 0 and dczi == -1) {
        if (pbc_z) cj.zi = n_cells_z -1;
        else continue;
      } else if(ci.zi == n_cells_z -1 and dczi == 1) {
        if (pbc_z) cj.zi = 0;
        else continue;
      } else {
        cj.zi = ci.zi + dczi;
      }

      // loop over list of particles in neighboring cell ci
      for (std::list<unsigned int>::iterator it = 
           cell_list[cj.xi][cj.yi][cj.zi].begin();
           it != cell_list[cj.xi][cj.yi][cj.zi].end();
           ++it) {
        
          pj = *it; 
          if ((!double_bond and pi >= pj) or (pi == pj)) {

            continue;
          
          }else if (dist2(positions[pj],positions[pi],Lx,Ly,Lz,pbc_x,pbc_y,pbc_z)
                  < r_verlet*r_verlet) {

            neighbor_list[pi].push_back(pj);
          }
      } // end loop over particles in cell cj
          

    }}} // end loop over neighboring cells

  } // end loop over particles (pi)

  return neighbor_list;

} // end function get_neighbor_list


VecList get_verlet_list(
  double Lx, double Ly, double Lz,
  bool pbc_x, bool pbc_y, bool pbc_z,
  bool double_bond,
  double r_verlet,
  const std::vector<Vec3>& positions)
{

  unsigned int n_particles = positions.size();

  VecList verlet_list(n_particles);

  double dist_sq;
  for (unsigned int i = 0;     i < n_particles; ++i) {
  for (unsigned int j = i + 1; j < n_particles; ++j) {

    dist_sq = dist2(positions[i],positions[j],Lx,Ly,Lz,pbc_x,pbc_y,pbc_z);
    if (dist_sq < r_verlet * r_verlet) {

      verlet_list[i].push_back(j);

      if (double_bond) {
        verlet_list[j].push_back(i);
      }

    }
  }} 

  return verlet_list;
}



#endif  // NEIGHBOR_LIST_H
