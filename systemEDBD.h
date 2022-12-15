#ifndef GUARD_SYSTEMEDBD_H
#define GUARD_SYSTEMEDBD_H

/*
 TO DO:
  + set max_diff_ in constructor 
  + who should update the Verlet list?????
*/



#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>

#include "vec3.h"

namespace systemEDBD_helper {

// calculate square distance between r1 and r2,
// apply periodic boundaires if L > 0
double distance_squared(Vec3 r1, Vec3 r2, double Lx, double Ly, double Lz)
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
             double dt,
             double verlet_list_radius,
             Potential potential);

  void SetPositions(const std::vector<Vec3>& positions);

  void MakeTimeStep(double dt);
  void Integrate(double delta_t);

  void SavePositions(std::string name) const;

  std::vector<Vec3> GetPositions() const { return positions_;}

  long unsigned int GetNumberOfVerletListUpdates() const
    { return number_of_verlet_list_updates_;}

  void ResetTime() { time_ = 0; }
  void SetTime(double new_time) { time_ = new_time;}

  double GetTime() const { return time_; }
  long unsigned int GetNumberOfParticles() const
    { return number_of_particles_; }

	// initialize the particles on a square lattice
	void RandomInit(unsigned int number_of_particles);

 private:


  const boost::normal_distribution<double> normal_distribution_;
  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<double> > random_normal_distribution_;


  void UpdateVerletList();

  void UpdateVelocities(double dt);

  void GetNextCollision(unsigned int& p1,
        unsigned int& p2, double& dt_collision ) const;
  void MoveBallistically(double dt);
  void MakeCollision(unsigned int p1, unsigned int p2);
  unsigned int number_of_particles_;

  // system_size in the x, y, and z directions
  double system_size_x_;
  double system_size_y_;
  double system_size_z_;

  // time step size
  double dt_;

  // Radius of the Verlet list
  double verlet_list_radius_;

  // particle positions
  std::vector<Vec3> positions_;
  std::vector<Vec3> positions_at_last_update_;
  // particle velocities
  std::vector<Vec3> velocities_;

  // Verlet list
  std::vector<std::vector<unsigned int> > verlet_list_;
  // number of neighbors in the Verlet list
  std::vector<unsigned int> number_of_neighbors_;


  double max_diff_;

  Potential potential_;

  // keep track of the number of Verlet list updates
  unsigned long int number_of_verlet_list_updates_;

  // current time
  double time_;

  double D_;
  double gamma_;

};


template <class Potential>
SystemEDBD<Potential>::SystemEDBD(
  unsigned long int seed,
  double system_size_x,
  double system_size_y,
  double system_size_z,
  double dt,
  double verlet_list_radius,
  Potential potential)
  : normal_distribution_(0.0, 1.0),
    random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_,
                                normal_distribution_),
    system_size_x_(system_size_x),
    system_size_y_(system_size_y),
    system_size_z_(system_size_z),
    dt_(dt),
    verlet_list_radius_(verlet_list_radius),
    potential_(potential),
    number_of_verlet_list_updates_(0),
    time_(0.0)
{

  max_diff_ = 1.0;
}

template <class Potential>
void SystemEDBD<Potential>::SetPositions(
    const std::vector<Vec3>& positions)
{
  number_of_particles_ = positions.size();
  positions_ = positions;
  positions_at_last_update_ =
    std::vector<Vec3>(number_of_particles_);

  verlet_list_ =
  std::vector<std::vector<unsigned int> >(number_of_particles_,
    std::vector<unsigned int>(number_of_particles_));

  number_of_neighbors_ =
    std::vector<unsigned int>(number_of_particles_);

  velocities_ = std::vector<Vec3>(number_of_particles_);
  UpdateVerletList();
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

  UpdateVelocities(dt);

  unsigned int p1, p2;
  double dt_collision;

  GetNextCollision(p1, p2, dt_collision);

  while (dt_collision < dt) {
    MoveBallistically(dt_collision);
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
  bool update_verlet_list = false;

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    velocities_[i] = 0;
    velocities_[i] =
        dt * potential_.Force(positions_[i], time_) / gamma_;

    velocities_[i].x +=
        sqrt_2_dt * random_normal_distribution_();
    velocities_[i].y +=
        sqrt_2_dt * random_normal_distribution_();
    velocities_[i].z +=
        sqrt_2_dt * random_normal_distribution_();

    double dist = systemEDBD_helper::distance_squared(positions_[i],
                  positions_at_last_update_[i], system_size_x_,
                  system_size_y_, system_size_z_);
    if (dist > max_diff_ * max_diff_) update_verlet_list = true;
  }
  
  if (update_verlet_list) UpdateVerletList();
}

template <class Potential>
void SystemEDBD<Potential>::MoveBallistically(double dt)
{
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_[i] += velocities_[i] * dt;
  }

  time_ += dt;
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
void SystemEDBD<Potential>::UpdateVerletList()
{
  number_of_verlet_list_updates_ += 1;
  //std::cout << "update Verlet List\t";
  //std::cout << number_of_verlet_list_updates_ << std::endl;

  std::fill(number_of_neighbors_.begin(),
		    number_of_neighbors_.end(), 0);

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_at_last_update_[i] = positions_[i];
    for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
      if (systemEDBD_helper::distance_squared(positions_[i],
				                    positions_[j], system_size_x_,
                                    system_size_y_, system_size_z_)
			< verlet_list_radius_ * verlet_list_radius_ ) {
        verlet_list_[i][ number_of_neighbors_[i] ] = j;
        ++number_of_neighbors_[i];
        //verlet_list_[j][ number_of_neighbors_[j] ] = i;
        //++number_of_neighbors_[j];
      }
    }
  }

}


#endif

