
#include "system.h"

template <class Potential>
System<Potential>::System(long unsigned int seed, unsigned int number_of_particles,
    double system_size_xy, double max_step_size, double verlet_radius,
    double max_displacement)
  : uniform_distribution(-1,1),
    random_number_generator(seed),
    random_uniform_distribution(random_number_generator, uniform_distribution),
    number_of_particles_(number_of_particles),
    system_size_xy_(system_size_xy),
    max_step_size_(max_step_size),
    verlet_radius_(verlet_radius),
    max_displacement_(max_displacement),
    positions_(number_of_particles),
    positions_at_last_verlet_update_(number_of_particles),
    verlet_list_(number_of_particles, std::vector<unsigned int>(number_of_particles)),
    number_of_attempted_moves_(0),
    number_of_accepted_moves_(0),
    number_of_verlet_updates_(0)
{

  // random initialization of positions
  InitRandom();

  // initialize verlet list
  UpdateVerletList();

}


template <class Potential>
void System<Potential>::InitRandom()
{
}

template <class Potential>
void System<Potential>::UpdateVerletList()
{
}

