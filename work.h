#ifndef WORK_H
#define WORK_H

#include <string>
#include <vector>

////////////////
// PROTOTYPES //
////////////////

double* activeWork(std::string filename, int const& tau);
  // Returns active work of particles within data file for periods of time
  // indicated by number of frames.

#endif
