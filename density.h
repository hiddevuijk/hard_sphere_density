#ifndef GUARD_DENSITY_H
#define GUARD_DENSITY_H

#include "vec3.h"

#include <vector>
#include <fstream>

class Density {
 public:
  Density(double xmin, double xmax, unsigned int number_of_bins, char xyz);
  void Sample(const std::vector<Vec3>& positions);
 private:

  // range of the density
  double xmin_, xmax_;
  // bin size
  double dx_;

  // total number of particles sampled
  long unsigned int number_of_samples_;

  // if a position is outside xmin_ -- xmax_
  // the samples is excluded
  long unsigned int number_of_excluded_samples_;
};

#endif
