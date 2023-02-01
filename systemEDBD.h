#ifndef GUARD_SYSTEMEDBD_H
#define GUARD_SYSTEMEDBD_H

/*
 TO DO:

  + Test MakeCollision
  + Test GetNextCollision

  + set max_diff_ in constructor 

*/



#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>

#include "vec3.h"
#include "neighbor_list.h"

namespace systemEDBD_helper {

// calculate square distance between r1 and r2,
// apply periodic boundaires if L > 0
double distance_squared(Vec3 r1, Vec3 r2,
                      double Lx, double Ly, double Lz,
                      bool pbc_x, bool pbc_y, bool pbc_z)
{
  r1 -= r2;

  if (pbc_x) {
    r1.x = fabs(r1.x);
    r1.x -= static_cast<int>(r1.x / Lx + 0.5) * Lx;
  }
  if (pbc_y) {
    r1.y = fabs(r1.y);
    r1.y -= static_cast<int>(r1.y / Ly + 0.5) * Ly;
  }
  if (pbc_z) {
    r1.z = fabs(r1.z);
    r1.z -= static_cast<int>(r1.z / Lz + 0.5) * Lz;
  }

	return r1.LengthSquared();

}

double distance_squared(Vec3 r1, Vec3 r2,
                        double Lx, double Ly, double Lz)
{
     
	r1 -= r2;
  if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);
	return r1.LengthSquared();
}

double distance(Vec3 r1, Vec3 r2, double Lx, double Ly, double Lz)
{
	r1 -= r2;
  if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);
	return r1.Length();
}

};


template <class Potential>
class SystemEDBD {
 public:
  SystemEDBD(unsigned long int seed,
             double system_size_x,
             double system_size_y,
             double system_size_z,
             bool pbc_x, bool pbc_y, bool pbc_z,
             double dt,
             double verlet_list_radius,
             Potential potential,
             double D, double gamma);

  //void SetR(unsigned int i, double x, double y, double z)
  //  {positions_[i].x=x;positions_[i].y=y; positions_[i].z=z;}

  //void SetV(unsigned int i, double vx, double vy, double vz)
  //{velocities_[i].x=vx;velocities_[i].y=vy;velocities_[i].z=vz;}

  Vec3 GetPosition(unsigned int i) const
    { return positions_[i]; }

  std::vector<Vec3> GetPositions() const 
    { return positions_; }

  Vec3 GetVelocity(unsigned int i) const
    { return velocities_[i]; }

  unsigned int GetNColl() const { return Ncoll; }
  double GetVerletListRadius() const { return verlet_list_radius_;}

  void SetPositions(const std::vector<Vec3>& positions);

  void MakeTimeStep(double dt);
  void RetakeTimeStep(double dt);
  void Integrate(double delta_t);

  void SavePositions(std::string name) const;

  Potential& GetPotential() { return potential_;}

  long unsigned int GetNumberOfNeighborListUpdates() const
    { return number_of_neighbor_list_updates_;}

  void ResetTime() { time_ = 0; }
  void SetTime(double new_time) { time_ = new_time;}

  double GetTime() const { return time_; }
  long unsigned int GetNumberOfParticles() const
    { return number_of_particles_; }

	// initialize the particles on a square lattice
	void RandomInit(unsigned int number_of_particles);

  bool CheckOverlaps() const;

 private:


  const boost::normal_distribution<double> normal_distribution_;
  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<double> > random_normal_distribution_;


  void UpdateNeighborList();

  void UpdateVelocities(double dt);

  void GetNextCollision(unsigned int& p1,
        unsigned int& p2, double& dt_collision ) const;

  // returns collision time of p1 and p2
  // returns a negative number if they don't collide
  double PairTime(unsigned int p1, unsigned int p2) const;

  // Move all particles v_i * dt forward
  // This is the only method that moves the system
  // forward in time. 
  // Other methods can only do so by calling MoveBallistically.
  void MoveBallistically(double dt);

  // p1 and p2 are colliding, so reset velocities
  // accoring to a elastic collision
  void MakeCollision(unsigned int p1, unsigned int p2);


  unsigned int number_of_particles_;

  // system_size in the x, y, and z directions
  double system_size_x_;
  double system_size_y_;
  double system_size_z_;

  bool pbc_x_;
  bool pbc_y_;
  bool pbc_z_;

  // time step size
  double dt_;

  // Radius of the Verlet list
  double verlet_list_radius_;

  // particle positions
  std::vector<Vec3> positions_;
  std::vector<Vec3> positions_at_last_update_;
  // particle velocities
  std::vector<Vec3> velocities_;


  std::vector<std::list<unsigned int> > neighbor_list_;

  double max_diff_;

  Potential potential_;

  // keep track of the number of Verlet list updates
  unsigned long int number_of_neighbor_list_updates_;

  // current time
  double time_;

  double D_;
  double gamma_;

  unsigned int Ncoll;
};


template <class Potential>
SystemEDBD<Potential>::SystemEDBD(
  unsigned long int seed,
  double system_size_x,
  double system_size_y,
  double system_size_z,
  bool pbc_x, bool pbc_y, bool pbc_z,
  double dt,
  double verlet_list_radius,
  Potential potential,
  double D, double gamma)
  : normal_distribution_(0.0, 1.0),
    random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_,
                                normal_distribution_),
    system_size_x_(system_size_x),
    system_size_y_(system_size_y),
    system_size_z_(system_size_z),
    pbc_x_(pbc_x), pbc_y_(pbc_y), pbc_z_(pbc_z),
    dt_(dt),
    verlet_list_radius_(verlet_list_radius),
    potential_(potential),
    number_of_neighbor_list_updates_(0),
    time_(0.0),
    D_(D), gamma_(gamma)
    ,Ncoll(0)
{
  max_diff_ = (verlet_list_radius_ - 1.0) / 2.0;
}

template <class Potential>
void SystemEDBD<Potential>::SetPositions(
    const std::vector<Vec3>& positions)
{

  number_of_particles_ = positions.size();
  positions_ = positions;

  positions_at_last_update_ =
    std::vector<Vec3>(number_of_particles_);

  neighbor_list_ = std::vector<std::list<unsigned int> >(number_of_particles_);

  velocities_ = std::vector<Vec3>(number_of_particles_);

  UpdateNeighborList();
}

template <class Potential>
void SystemEDBD<Potential>::Integrate(double delta_time)
{
  while (delta_time > dt_) {
    MakeTimeStep(dt_);
    delta_time -= dt_;
  }
  MakeTimeStep(delta_time);
}


template <class Potential>
void SystemEDBD<Potential>::MakeTimeStep(double dt)
{
  if (dt <= 0) return;

  UpdateVelocities(dt);

  unsigned int p1, p2;
  double dt_collision;

  GetNextCollision(p1, p2, dt_collision);
  
  // if collision happens within dt
  while (dt_collision < dt) {
    // move particles untill next collision
    MoveBallistically(dt_collision);
    // collide particles
    MakeCollision(p1, p2);

    dt -= dt_collision;
    GetNextCollision(p1, p2, dt_collision);
  }

  MoveBallistically(dt);

}

template <class Potential>
void SystemEDBD<Potential>::RetakeTimeStep(double dt)
{
  // Check if the max of dt * v_i is not larger than
  // the max allowed distance between updates of the 
  // Verlet list.
  double maxdr = 0;
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    double dr = dt * velocities_[i].Length();
    if (dr > maxdr) maxdr = dr;
  }

  // increase the Verlet list radius,
  // such that max_diff_ = 1.5 x maxdr
  // because otherwise an infinite loop occurs,
  if (maxdr > max_diff_) {
    verlet_list_radius_ = 1.0 + 3 * maxdr;
    max_diff_ = (verlet_list_radius_ - 1.0) / 2.0;
    std::cout << "\t" << time_ << "\t" << verlet_list_radius_ << "\t" <<  max_diff_ << std::endl;
  }

  // update the neighbor list
  UpdateNeighborList();

  if (dt <= 0) return;

  unsigned int p1, p2;
  double dt_collision;

  GetNextCollision(p1, p2, dt_collision);
  
  // if collision happens within dt
  while (dt_collision < dt) {
    // move particles untill next collision
    MoveBallistically(dt_collision);
    // collide particles
    MakeCollision(p1, p2);

    dt -= dt_collision;
    GetNextCollision(p1, p2, dt_collision);
  }

  MoveBallistically(dt);

}
template <class Potential>
void SystemEDBD<Potential>::UpdateVelocities(double dt)
{

  double sqrt_2_dt = sqrt(2 * D_ / dt);

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    velocities_[i] *= 0;

    velocities_[i] =
         potential_.Force(positions_[i], time_) / gamma_;

    velocities_[i].x +=
        sqrt_2_dt * random_normal_distribution_();
    velocities_[i].y +=
        sqrt_2_dt * random_normal_distribution_();
    velocities_[i].z +=
       sqrt_2_dt * random_normal_distribution_();
  }
  
}

template <class Potential>
void SystemEDBD<Potential>::MoveBallistically(double dt)
{

  bool update_neighbor_list = false;
  Vec3 new_position_i;
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    new_position_i = positions_[i] + velocities_[i] * dt;


    // check if neighbor list needs to be updated
    double dist = systemEDBD_helper::distance_squared(
        new_position_i, positions_at_last_update_[i],
        system_size_x_, system_size_y_, system_size_z_,
        pbc_x_, pbc_y_, pbc_z_);

    if (dist > max_diff_ * max_diff_) {
      update_neighbor_list = true;
      break;
    }

  }
  // If one of the particles moved to far
  // the Verlet list needs to be updated,
  // and the time step needs to be retaken
  // (i.e. recalculate the next collision etc.).
  if (update_neighbor_list) {
    RetakeTimeStep(dt);
  } else {
    // move particles 
    for (unsigned int i = 0; i < number_of_particles_; ++i) {
      positions_[i] += velocities_[i] * dt;
    }

    time_ += dt;
  }
}


template <class Potential>
double SystemEDBD<Potential>::PairTime(unsigned int p1, unsigned int p2) const
{

  Vec3 dr = positions_[p1] - positions_[p2];
  // periodic boundary conditions
  if (system_size_x_ > 0) dr.x -=
        system_size_x_ * round(dr.x/system_size_x_);
	if (system_size_y_ > 0) dr.y -=
        system_size_y_ * round(dr.y/system_size_y_);
	if (system_size_z_ > 0) dr.z -=
        system_size_z_ * round(dr.z/system_size_z_);

  Vec3 dv = velocities_[p1] - velocities_[p2];

  double dr_dv = dr.DotProduct(dv);

  // particles don't collide
  if (dr_dv >= 0) return -1;

  double determinant = dr.LengthSquared() - 1;
  determinant *= -dv.LengthSquared();
  determinant += dr_dv * dr_dv;

  // particles don't collide
  if (determinant <= 0) return -1;

  double dtc=  - (dr_dv + sqrt(determinant)) / dv.LengthSquared();
  return dtc;
}


template <class Potential>
void SystemEDBD<Potential>::GetNextCollision(unsigned int& p1,
        unsigned int& p2, double& dt_collision ) const
{
  // no collision found yet
  // set to dt_, because any collision time longer 
  // than that is not important
  dt_collision = dt_; 

  double dt_pi_pj; // collision time of pi and pj
  for (unsigned int pi = 0; pi < number_of_particles_; ++pi) {
    for (std::list<unsigned int>::const_iterator it=neighbor_list_[pi].begin(); it != neighbor_list_[pi].end(); ++it) {

      dt_pi_pj = PairTime(pi, *it); 
      if (dt_pi_pj < dt_collision and dt_pi_pj > 0) {
        p1 = pi;
        p2 = *it;
        dt_collision = dt_pi_pj;
      }
    }
  }
  //for (unsigned int pi = 0; pi < number_of_particles_; ++pi) {
  //  // pi_n is neighbor number n of particle pi
  //  for (unsigned int pi_n = 0;
  //      pi_n< number_of_neighbors_[pi]; ++pi_n){

  //    pj = verlet_list_[pi][pi_n];
  //    dt_pi_pj = PairTime(pi, pj); 
  //    if (dt_pi_pj < dt_collision and dt_pi_pj > 0) {
  //      p1 = pi;
  //      p2 = pj;
  //      dt_collision = dt_pi_pj;
  //    }
  //  }
  //}
}

template <class Potential>
void SystemEDBD<Potential>::MakeCollision(unsigned int p1, unsigned int p2)
{
  Ncoll += 1; 
  Vec3 dr = positions_[p1] - positions_[p2];
  // periodic boundary conditions
  if (pbc_x_) {
    dr.x -= system_size_x_ * round(dr.x/system_size_x_);
  }
	if (pbc_y_) {
    dr.y -= system_size_y_ * round(dr.y/system_size_y_);
  }
	if (pbc_z_) {
    dr.z -= system_size_z_ * round(dr.z/system_size_z_);
  }

  Vec3 n_perp = dr / dr.Length();

  Vec3 dv = velocities_[p1] - velocities_[p2];

  double a = dv.DotProduct(n_perp);

  velocities_[p1] -= n_perp * a;
  velocities_[p2] += n_perp * a;
}


template <class Potential>
void SystemEDBD<Potential>::SavePositions(std::string name) const
{
  std::ofstream out;
  out.open(name);
  out <<std::setprecision(16);
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
	out << positions_[i].x << '\t'
        << positions_[i].y << '\t'		
        << positions_[i].z << '\n';
  }

	out.close();
}


template <class Potential>
void SystemEDBD<Potential>::UpdateNeighborList()
{
  number_of_neighbor_list_updates_ += 1;


  neighbor_list_ = get_neighbor_list(
          system_size_x_, system_size_y_, system_size_z_,
          pbc_x_, pbc_y_, pbc_z_,true, verlet_list_radius_,
          positions_);

  positions_at_last_update_ = positions_;
}

template <class Potential>
bool SystemEDBD<Potential>::CheckOverlaps() const
{

  bool overlap = false;
  // loop over all particle pairs
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
  for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
    Vec3 dr = positions_[i] - positions_[j];
    // if distance < diameter, there is an overlap
    if (dr.Length() < 1.0)  {
      overlap = true;
      break;
    }
  }}

  return overlap;
}

#endif

