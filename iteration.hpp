#ifndef ITERATION_HPP
#define ITERATION_HPP

#include "particle.hpp"

void iterate_ABP_WCA(System* system);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential, using custom dimensionless parameters
  // relations.
  // [Euler's scheme]

void iterate_ABP_WCA(System0* system);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.
  // [Euler's scheme]

template<class SystemClass> void ABP_WCA(
  SystemClass* system, std::vector<Particle>& newParticles) {
  // Compute interactions with WCA potentials between all particles of the
  // system.

  pairs_ABP<SystemClass>(system,
    [&system, &newParticles](int index1, int index2)
      { system->WCA_force(index1, index2, newParticles); });
}

#endif
