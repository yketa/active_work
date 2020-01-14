#ifndef ITERATION_HPP
#define ITERATION_HPP

#include "particle.hpp"

void iterate_ABP_WCA(System* system, int Niter);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential, using custom dimensionless parameters
  // relations.

void iterate_ABP_WCA(System0* system, int Niter);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
void aligningTorque(System* system);
  // Compute aligning torques between all particles of the system.
#endif

template<class SystemClass> void ABP_WCA(SystemClass* system) {
  // Compute interactions with WCA potentials between all particles of the
  // system.

  pairs_ABP<SystemClass>(system,
    [&system](int index1, int index2) { system->WCA_force(index1, index2); });
}

template<class SystemClass> void ABP_WCA(
  SystemClass* system, std::vector<Particle>& newParticles) {
  // Compute interactions with WCA potentials between all particles of the
  // system.

  pairs_ABP<SystemClass>(system,
    [&system, &newParticles](int index1, int index2)
      { system->WCA_force(index1, index2, newParticles); });
}

#endif
