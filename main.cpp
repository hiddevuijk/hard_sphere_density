#include "config_file.h"
#include "vec3.h"
#include "system.h"

#include <iostream>

using namespace std;


class Potential { 
 public:
 Potential(double a) : a(a) {}
 Potential() : a(0) {}
 Vec3 Force(const Vec3& r1, const Vec3& r2) {
	 return r1 - r2;
}
 double a;
};

int main()
{
  Config params("input.txt");
  unsigned long int seed = params.get_parameter<unsigned long int>("seed");
  unsigned int number_of_particles = params.get_parameter<unsigned int>("number_of_particles");
  double system_size_xy = params.get_parameter<double>("system_size_xy");
  double max_mc_step_size = params.get_parameter<double>("max_mc_step_size");
  double verlet_list_radius = params.get_parameter<double>("verlet_list_radius");
  double max_diff = 0.1;
  Potential potential;

  System<Potential> system(seed, number_of_particles, system_size_xy,
					max_mc_step_size, verlet_list_radius,
					max_diff, potential);

  system.SavePositions("positions.dat"); 

  return 0;
}
