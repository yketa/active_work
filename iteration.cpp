#include <cmath>
#include <math.h>
#include <vector>

#include "iteration.h"
#include "particle.h"

void iterate_ABP_WCA(System *system) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

  std::vector<Particle> newParticles(system->getNumberParticles());

  // COMPUTATION
  for (int i=0; i < system->getNumberParticles(); i++) {

    // POSITIONS
    for (int dim=0; dim < 2; dim++) {
      newParticles[i].position()[dim] =
        (system->getParticle(i))->position()[dim]; // initialise new positions with previous ones
      newParticles[i].position()[dim] +=
        sqrt(system->getTimeStep()*2/3/system->getPersistenceLength())
          *(system->getRandomGenerator())->gauss_cutoff(); // add noise
      newParticles[i].position()[dim] +=
        system->getTimeStep()*cos(
            (system->getParticle(i))->orientation()[0] - dim*M_PI/2); // add self-propulsion
    }

    // ORIENTATIONS
    newParticles[i].orientation()[0] =
      (system->getParticle(i))->orientation()[0]; // initialise new orientation with previous one
    newParticles[i].orientation()[0] +=
      sqrt(system->getTimeStep()*2/system->getPersistenceLength())
        *(system->getRandomGenerator())->gauss_cutoff(); // add noise
  }

  double force[2];
  for (int i=0; i < system->getNumberParticles(); i++) {
    for (int j=i+1; j < system->getNumberParticles(); j++) {

      // INTERACTIONS
      system->WCA_potential(i, j, &force[0]); // force acting on i from j

      for (int dim=0; dim < 2; dim++){
        newParticles[i].position()[dim] +=
          system->getTimeStep()*force[dim]/3/system->getPersistenceLength(); // add force to i
        newParticles[j].position()[dim] -=
          system->getTimeStep()*force[dim]/3/system->getPersistenceLength(); // substract force to j
      }
    }
  }

  // SAVE AND COPY
  system->saveNewState(&newParticles[0]);
}
