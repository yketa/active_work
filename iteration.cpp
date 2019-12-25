#include <cmath>
#include <math.h>
#include <vector>

#include "iteration.hpp"
#include "particle.hpp"

#ifdef DEBUG
#include <iostream>
#endif

void iterate_ABP_WCA(System* system) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential, using custom dimensionless parameters
  // relations.

  Parameters* parameters = system->getParameters();

  std::vector<Particle> newParticles(parameters->getNumberParticles());

  // COMPUTATION
  for (int i=0; i < parameters->getNumberParticles(); i++) {

    // POSITIONS
    for (int dim=0; dim < 2; dim++) {
      newParticles[i].position()[dim] =
        (system->getParticle(i))->position()[dim]; // initialise new positions with previous ones
      newParticles[i].position()[dim] +=
        sqrt(parameters->getTimeStep()*2/3/parameters->getPersistenceLength())
          *(system->getRandomGenerator())->gauss_cutoff(); // add noise
      newParticles[i].position()[dim] +=
        #if CONTROLLED_DYNAMICS
        (1.0 - 2*parameters->getBiasingParameter()
          /3/parameters->getPersistenceLength())*
        #endif
        parameters->getTimeStep()*cos(
            (system->getParticle(i))->orientation()[0] - dim*M_PI/2); // add self-propulsion

      (system->getParticle(i))->force()[dim] = 0; // initialise force
    }

    // ORIENTATIONS
    newParticles[i].orientation()[0] =
      (system->getParticle(i))->orientation()[0]; // initialise new orientation with previous one
    newParticles[i].orientation()[0] +=
      sqrt(parameters->getTimeStep()*2/parameters->getPersistenceLength())
        *(system->getRandomGenerator())->gauss_cutoff(); // add noise
  }

  // ALIGNING TORQUE
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  double torque;
  for (int i=0; i < parameters->getNumberParticles(); i++) {
    #if 1
    for (int j=i + 1; j < parameters->getNumberParticles(); j++) {
      torque = 2.0*system->getTorqueParameter()/parameters->getNumberParticles()
        *sin((system->getParticle(i))->orientation()[0]
          - (system->getParticle(j))->orientation()[0])
        *parameters->getTimeStep();
      newParticles[i].orientation()[0] += torque;
      newParticles[j].orientation()[0] -= torque;
    }
    #else
    for (int j=0; j < parameters->getNumberParticles(); j++) {
      torque = 2.0*system->getTorqueParameter()/parameters->getNumberParticles()
        *sin((system->getParticle(i))->orientation()[0]
          - (system->getParticle(j))->orientation()[0])
        *parameters->getTimeStep();
      newParticles[i].orientation()[0] += torque;
    }
    #endif
  }
  #endif

  // FORCES
  #if USE_CELL_LIST // with cell list
  cellList_ABP_WCA<System>(system, newParticles);
  #else // brute force
  bruteForce_ABP_WCA<System>(system, newParticles);
  #endif

  // SAVE AND COPY
  system->saveNewState(newParticles);
}


void iterate_ABP_WCA(System0* system) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

  Parameters* parameters = system->getParameters();

  std::vector<Particle> newParticles(system->getNumberParticles());

  // COMPUTATION
  for (int i=0; i < parameters->getNumberParticles(); i++) {

    // DIAMETERS
    newParticles[i].diameter()[0] = (system->getParticle(i))->diameter()[0];

    // POSITIONS
    for (int dim=0; dim < 2; dim++) {
      newParticles[i].position()[dim] =
        (system->getParticle(i))->position()[dim]; // initialise new positions with previous ones
      newParticles[i].position()[dim] +=
        sqrt(parameters->getTimeStep()*2*parameters->getTransDiffusivity())
          *(system->getRandomGenerator())->gauss_cutoff(); // add noise
      newParticles[i].position()[dim] +=
        parameters->getTimeStep()*parameters->getPropulsionVelocity()*cos(
            (system->getParticle(i))->orientation()[0] - dim*M_PI/2); // add self-propulsion

      (system->getParticle(i))->force()[dim] = 0; // initialise force
    }

    // ORIENTATIONS
    newParticles[i].orientation()[0] =
      (system->getParticle(i))->orientation()[0]; // initialise new orientation with previous one
    newParticles[i].orientation()[0] +=
      sqrt(parameters->getTimeStep()*2*parameters->getRotDiffusivity())
        *(system->getRandomGenerator())->gauss_cutoff(); // add noise
  }

  // FORCES
  #if USE_CELL_LIST // with cell list
  cellList_ABP_WCA<System0>(system, newParticles);
  #else
  bruteForce_ABP_WCA<System0>(system, newParticles);
  #endif

  // SAVE AND COPY
  system->saveNewState(newParticles);
}
