#include "config_file.h"
#include "vec3.h"
#include "system.h"
#include "density.h"

#include <iostream>

using namespace std;

bool CheckOverlap( vector<Vec3> r, double L)
{
  bool overlap = false;
  int N = r.size();
  for(int i = 0; i < N; i++) {
	for(int j = i+1; j < N; j++) {
      double d = system_helper::distance(r[i], r[j], L);
      if ( d < 1 ) {
        overlap = true;
        break;
	  }
	}
    if (overlap) break;
  }

	return overlap;
}

int main()
{
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  unsigned int number_of_particles =
		params.get_parameter<unsigned int>("number_of_particles");

  double system_size_xy =
		params.get_parameter<double>("system_size_xy");

  double max_mc_step_size =
		params.get_parameter<double>("max_mc_step_size");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");


  double A = params.get_parameter<double>("A");
  double Ae = params.get_parameter<double>("Ae");

  long unsigned int MC_moves_per_sample = params.get_parameter<long unsigned int>("MC_moves_per_sample");

  long unsigned int initial_MC_moves = params.get_parameter<long unsigned int>("initial_MC_moves");

  double zlim = params.get_parameter<double>("zlim");
  unsigned int number_of_bins =
		params.get_parameter<double>("number_of_bins");

  unsigned int number_of_samples = params.get_parameter<unsigned int>("number_of_samples");

  System system(seed, number_of_particles, system_size_xy,
					max_mc_step_size, verlet_list_radius, A);
  system.SetPotentialExp(Ae);

  double area = system_size_xy * system_size_xy;
  Density rho_z(-zlim, zlim, number_of_bins, 'z', area);

  string positions_name = "positions0.dat";
  system.SavePositions(positions_name); 


  //system.MCMoveNoVerletFull(initial_MC_moves);
  system.MCMoveFull(initial_MC_moves);

  cout << "Init done\n" << flush;
  cout << (1.0 * system.GetNumberOfAcceptedMoves() ) /
			system.GetNumberOfAttemptedMoves() << endl << flush;



  for (long unsigned int i = 0; i < number_of_samples; i++) {
    //system.MCMoveNoVerletFull(MC_moves_per_sample);
  	system.MCMoveFull(MC_moves_per_sample);
    rho_z.Sample(system.GetPositions());


    cout << number_of_samples << "\t" << i << '\n' << flush;
  }

  string density_name = "rhoz.dat";
  rho_z.Save(density_name);

  cout << "Monte-Carlo move acceptance rate:\n" << flush;
  cout << (1.0 * system.GetNumberOfAcceptedMoves() ) /
			system.GetNumberOfAttemptedMoves() << endl << flush;

  if(CheckOverlap(system.GetPositions(), system_size_xy ) ) cout << "FUCK" << endl << flush;

  positions_name = "positions.dat";
  system.SavePositions(positions_name); 

  cout << "Accepted moves per particle per Verlet list update:\n" << flush;
  cout << system.GetNumberOfAcceptedMoves() * 1.0
       / (system.GetNumberOfVerletListUpdates() * number_of_particles)
	   << endl << flush;

  cout << "Attempted moves per particle per Verlet list update:\n" << flush;
  cout << system.GetNumberOfAttemptedMoves() * 1.0
       / (system.GetNumberOfVerletListUpdates() * number_of_particles)
	   << endl << flush;

  return 0;
}
