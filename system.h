/*
 To do:
	- fix initialization
*/


#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <boost/random.hpp>

#include "vec3.h"
#include <string>

double distance(Vec3 r1, Vec3 r2, double L)
{
	r1 -= r2;
	r1.x -= round(L/r1.x);
	r1.y -= round(L/r1.y);
	return r1.Length();
}

class System {
 public:
	System(unsigned long int seed,
		   unsigned int number_of_particles,
		   double system_size_xy,
		   double max_mc_step_size,
		   double verlet_list_radius,
		   double max_diff);
  void SavePositions(std::string name) const;		    
  // attempt an MC move
  void MCMove();
  // attempt an MC move without using the Verlet list
  void MCMoveNoVerlet();		

  // make number_of_paricles_ MC move attempts
  void MCMoveNoVerletFull() {
    for (unsigned int i = 0; i < number_of_particles_; ++i) {
      MCMoveNoVerlet();
	}
  }

 private:
	// uniform distribution [-1,1]
	const boost::uniform_real<double> uniform_distribution_;
    const boost::random::uniform_int_distribution<int>
                                       random_int_distribution_;

	boost::mt19937 random_number_generator_;
	boost::variate_generator<boost::mt19937&,
			boost::uniform_real<double> > random_uniform_distribution_11_;
	boost::variate_generator<boost::mt19937&,
			boost::uniform_real<double> > random_uniform_distribution_01_;

	// initialize the particles on a square lattice
	void RandomInit();

	void UpdateVerletList();



	// private variable

	unsigned int number_of_particles_;
	// system size in the x and y direction
	double system_size_xy_;
	
	// max step size of an MC move
	double max_mc_step_size_;

	// Radius of for the Verlet list
	double verlet_list_radius_;

	// particle positions
	std::vector<Vec3> positions_;
	// particle positions at when the Verlet list was last updated
	std::vector<Vec3> positions_at_last_update_;

	// Verlet list
	std::vector<std::vector<int> > verlet_list_;

	// when distance between position_[i] and position_at_last_update_[i]
	// is larger than max_diff_, the Verlet list needs to be updated
	double max_diff_;	


	// keep track of the performance of the MC algorithm
	unsigned long int number_of_attempted_moves_;
	unsigned long int number_of_accepted_moves_;
	unsigned long int number_of_verlet_list_updates_;
	
};

System::System(
	unsigned long int seed,
	unsigned int number_of_particles,
	double system_size_xy,
	double max_mc_step_size,
	double verlet_list_radius,
	double max_diff)
  : uniform_distribution_(-1,1),
    random_int_distribution_(0, number_of_particles - 1),
	random_number_generator_(seed),
	random_uniform_distribution_11_(random_number_generator_,
							     uniform_distribution_),
	random_uniform_distribution_01_(random_number_generator_,
							     uniform_distribution_),
	number_of_particles_(number_of_particles),
	system_size_xy_(system_size_xy),
	max_mc_step_size_(max_mc_step_size),
	verlet_list_radius_(verlet_list_radius),
    positions_(number_of_particles_),
    positions_at_last_update_(number_of_particles_),
	verlet_list_(number_of_particles_,
				std::vector<int>(number_of_particles)),
	max_diff_(max_diff)
{
  RandomInit();
  UpdateVerletList();
}


void System::RandomInit()
{

  std::vector<Vec3> lattice_positions;

  //if (number_of_particles_ < sqrt(system_size_xy_) ) {
  if (true) {

	int n = ceil( pow(number_of_particles_, 1./2) );
	double dx = system_size_xy_ / n;
	double dy = dx;

	int jx = 0;
	int jy = 0;
	Vec3 r;
	for (int i = 0; i < 1+n * n; ++i) {
      r.x = jx * dx;
	  r.y = jy * dy;
	  r.z = 0;
	  lattice_positions.push_back(r);
	  jx++;
	  if (jx == n) {
        jx = 0;	  
		jy++;
	  }
	}
  } else {
  }

 // shuggffle lattice_positions

  for (unsigned int i = 0; i < number_of_particles_; ++i){
    positions_[i] = lattice_positions[i];
  }
}


void System::UpdateVerletList()
{
	
}



void System::MCMoveNoVerlet()
{
  number_of_attempted_moves_++;

  // pick a random particle to move
  unsigned int i = random_int_distribution_(random_number_generator_);

  // random displacement
  Vec3 new_position(random_uniform_distribution_11_(),
				    random_uniform_distribution_11_(),
				    random_uniform_distribution_11_());
  new_position *= max_mc_step_size_;
  new_position += positions_[i];

  // check for overlap
  bool overlap = false;
  for (unsigned int j = 0; j < number_of_particles_; ++j) {
	if (j == i) continue;
    if (distance(new_position, positions_[j], system_size_xy_) < 1) {
      overlap = true;
      break;
	}
  }
  if (!overlap) {
    double delta_U = positions_[i].z * positions_[i].z;   
	delta_U -= new_position.z * new_position.z;
	if( random_uniform_distribution_01_() < std::exp(- delta_U ) ) {
		positions_[i] = new_position;

		// accepted move
		number_of_accepted_moves_++;
	}	

  }
}


void System::SavePositions(std::string name) const
{
  std::ofstream out;
  out.open(name);
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
	out << positions_[i].x << '\t'
        << positions_[i].y << '\t'		
        << positions_[i].z << '\n';
  }

	out.close();
}

#endif
