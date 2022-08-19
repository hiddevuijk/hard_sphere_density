#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/random.hpp>

#include "vec3.h"


template <class Potential>
class System {
 public:
  System(long unsigned int seed, unsigned int number_of_particles, double system_size_xy,
         double max_step_size, double verlet_radius, double max_displacement);

 private:
  // uniform distribution [-1,1]
  const boost::uniform_real<double> uniform_distribution;

  long unsigned int seed;
  boost::mt19937 random_number_generator;
  boost::variate_generator<boost::mt19937&,
      boost::uniform_real<double> > random_uniform_distribution;
 
  // random initialization on a lattice
  void InitRandom();

  // try to move a particle
  void MCMove();
  // try to move a particle using a Verlet list
  void MCMoveVerlet();

  // update Verlet List
  void UpdateVerletList();

  // private data members
  unsigned int number_of_particles_;
  double system_size_xy_;   // systemsize in the x and y direction
  double max_step_size_;   // max move distance

  double verlet_radius_; // radius of the verlet list
  // update verlet list if max displacement is 
  // larger than "max_displacement"
  double max_displacement_;

  std::vector<Vec3> positions_;
  std::vector<Vec3> positions_at_last_verlet_update_;
  std::vector<std::vector<unsigned int> > verlet_list_;


  // keep track of performance
  unsigned long int number_of_attempted_moves_;
  unsigned long int number_of_accepted_moves_;
  unsigned long int number_of_verlet_updates_;
   
};




#endif
