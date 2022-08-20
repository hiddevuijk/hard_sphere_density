/*
 To do:
	- fix initialization
    - add Verlet algorithm
    -- check GetPositions
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

namespace system_helper {

double distance_squared(Vec3 r1, Vec3 r2, double L)
{
	r1 -= r2;
	r1.x -= L * round(r1.x/L);
	r1.y -= L * round(r1.y/L);
	return r1.LengthSquared();
}

double distance(Vec3 r1, Vec3 r2, double L)
{
	r1 -= r2;
	r1.x -= L * round(r1.x/L);
	r1.y -= L * round(r1.y/L);
	return r1.Length();
}

};

class System {
 public:
	System(unsigned long int seed,
		   unsigned int number_of_particles,
		   double system_size_xy,
		   double max_mc_step_size,
		   double verlet_list_radius,
		   double A);
  void SavePositions(std::string name) const;		    
  // attempt an MC move
  void MCMove();

  // make n * number_of_particles_ MC move attempts
  void MCMoveFull(int n) {
    for (unsigned int i = 0; i < n * number_of_particles_; ++i) {
      MCMove();
	}
  }

  // attempt an MC move without using the Verlet list
  void MCMoveNoVerlet();		

  // make n * number_of_paricles_ MC move attempts
  void MCMoveNoVerletFull(int n) {
    for (unsigned int i = 0; i < n * number_of_particles_; ++i) {
      MCMoveNoVerlet();
	}
  }

  void SetPotential(double newA) { A_ = newA; }
  double getPotential() const { return A_; }

  std::vector<Vec3> GetPositions() const { return positions_; }

  long unsigned int GetNumberOfAttemptedMoves() const
		{ return number_of_attempted_moves_; }
  long unsigned int GetNumberOfAcceptedMoves() const
		{ return number_of_accepted_moves_; }
  long unsigned int GetNumberOfVerletListUpdates() const
		{ return number_of_verlet_list_updates_; }

 private:
	// uniform distribution [-1,1]
	const boost::uniform_real<double> uniform_distribution_11_;
	const boost::uniform_real<double> uniform_distribution_01_;
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
	std::vector<std::vector<unsigned int> > verlet_list_;
	// number of neighbors in the Verlet list
	std::vector<unsigned int> number_of_neighbors_;

	// when distance between position_[i] and position_at_last_update_[i]
	// is larger than max_diff_, the Verlet list needs to be updated
	double max_diff_;	

	// externale potential: U(z) = A_ * z * z
	double A_;

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
	double A)
  : uniform_distribution_11_(-1,1),
    uniform_distribution_01_(0,1),
    random_int_distribution_(0, number_of_particles - 1),
	random_number_generator_(seed),
	random_uniform_distribution_11_(random_number_generator_,
							     uniform_distribution_11_),
	random_uniform_distribution_01_(random_number_generator_,
							     uniform_distribution_01_),
	number_of_particles_(number_of_particles),
	system_size_xy_(system_size_xy),
	max_mc_step_size_(max_mc_step_size),
	verlet_list_radius_(verlet_list_radius),
    positions_(number_of_particles_),
    positions_at_last_update_(number_of_particles_),
	verlet_list_(number_of_particles_,
				std::vector<unsigned int>(number_of_particles)),
	number_of_neighbors_(number_of_particles_),
	A_(A),
	number_of_attempted_moves_(0),
	number_of_accepted_moves_(0),
	number_of_verlet_list_updates_(0)
	
{
  max_diff_ = (verlet_list_radius - 1) / 2;
  // maximum MC step size sqrt(3) * max Mc step in one dim
  max_diff_ -= sqrt(3.0) * max_mc_step_size_,

  RandomInit();
  UpdateVerletList();
}

void System::RandomInit()
{

  std::vector<Vec3> lattice_positions;

  double dx = 1.25;
  double dy = dx;
  double dz = 1.0;
  int n_per_xy = floor(system_size_xy_ / dx);

  Vec3 temp;
  int iz = 0;
  while (lattice_positions.size() < number_of_particles_) {
    int ix = 0;
    int iy = 0;
    while (iy < n_per_xy) {
      if (iz == 0) {
        temp.x = ix * dx;
        temp.y = iy * dy;
        temp.z = 0.0;
        lattice_positions.push_back(temp);
      } else {
        temp.x = ix * dx;
        temp.y = iy * dy;
        temp.z = iz * dz;
		if ( (iz % 2) == 1) {
		  temp.x += dx * sqrt(3)/2;
		  temp.y += dy * sqrt(3)/2;
        }
        lattice_positions.push_back(temp);

        temp.z = -iz * dz;
        lattice_positions.push_back(temp);
      }

      ix += 1;
      if (ix == n_per_xy) {
        ix = 0;
        iy += 1;
      }
    }
    iz += 1;
  }
  // shuffle positions

  for (unsigned int i = 0; i < lattice_positions.size(); ++i) {
    unsigned int j = random_int_distribution_(random_number_generator_);
	temp = lattice_positions[i];
    lattice_positions[i] = lattice_positions[j];
    lattice_positions[j] = temp;
  }

  for (unsigned int i = 0; i < number_of_particles_; ++i){
    positions_[i] = lattice_positions[i];
  }
}


void System::MCMove()
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
  // nj is the njth neighbor in the Verlet list,
  // its index is j = verlet_list_[i][nj]
  for (unsigned int nj = 0; nj < number_of_neighbors_[i]; ++nj) {
	unsigned int j = verlet_list_[i][nj];
    if (system_helper::distance_squared(new_position, positions_[j], system_size_xy_) < 1) {
      overlap = true;
      break;
	}
  }
  if (!overlap) {
    double delta_U = A_ * new_position.z * new_position.z;
    delta_U -= A_ * positions_[i].z * positions_[i].z;   
	if( random_uniform_distribution_01_() < std::exp(- delta_U ) ) {
		positions_[i] = new_position;

		// accepted move
		number_of_accepted_moves_++;

		// if the total distance it moved since the last update
		// is larger that max_diff, update the Verlet list
		double dist = system_helper::distance_squared(positions_[i], positions_at_last_update_[i], system_size_xy_);
		if (dist > max_diff_ * max_diff_) UpdateVerletList();
	}	

  }
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
    if (system_helper::distance_squared(new_position, positions_[j], system_size_xy_) < 1) {
      overlap = true;
      break;
	}
  }
  if (!overlap) {
    double delta_U = A_ * new_position.z * new_position.z;
    delta_U -= A_ * positions_[i].z * positions_[i].z;   
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
  out <<std::setprecision(16);
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
	out << positions_[i].x << '\t'
        << positions_[i].y << '\t'		
        << positions_[i].z << '\n';
  }

	out.close();
}

void System::UpdateVerletList()
{
  number_of_verlet_list_updates_ += 1;

  std::fill(number_of_neighbors_.begin(),
		    number_of_neighbors_.end(), 0);

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_at_last_update_[i] = positions_[i];
    for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
      if (system_helper::distance_squared(positions_[i],
				  positions_[j], system_size_xy_)
			< verlet_list_radius_ * verlet_list_radius_ ) {
        verlet_list_[i][ number_of_neighbors_[i] ] = j;
        verlet_list_[j][ number_of_neighbors_[j] ] = i;
        ++number_of_neighbors_[i];
        ++number_of_neighbors_[j];
      }
    }
  }

}






#endif
