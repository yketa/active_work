#ifndef ITERATION_HPP
#define ITERATION_HPP

#include "particle.hpp"
#include "maths.hpp"

void iterate_ABP_WCA(System *system);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

#endif
