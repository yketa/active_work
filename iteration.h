#ifndef ITERATION_H
#define ITERATION_H

#include "particle.h"
#include "maths.h"

////////////////
// PROTOTYPES //
////////////////

void iterate_ABP_WCA(System *system, double const& timeStep);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential..

#endif
