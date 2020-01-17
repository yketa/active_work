#include <cmath>
#include <math.h>
#include <vector>

#include "iteration.hpp"
#include "particle.hpp"

#ifdef DEBUG
#include <iostream>
#endif


///////////////////////////////
// ACTIVE BROWNIAN PARTICLES //
///////////////////////////////

void iterate_ABP_WCA(System* system, int Niter) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential, using custom dimensionless parameters
  // relations.

  Parameters* parameters = system->getParameters();

  std::vector<Particle> newParticles(parameters->getNumberParticles());

  double selfPropulsion; // self-propulsion force
  double noise; // noise realisation

  #if HEUN // HEUN'S SCHEME

  double selfPropulsionCorrection; // correction to the self-propulsion force

  std::vector<double> positions (2*parameters->getNumberParticles(), 0.0); // positions backup
  std::vector<double> forces (2*parameters->getNumberParticles(), 0.0); // forces backup
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  std::vector<double> orientations (parameters->getNumberParticles(), 0.0); // orientations backup
  std::vector<double> torques (parameters->getNumberParticles(), 0.0); // torques backup
  #endif
  #endif

  for (int iter=0; iter < Niter; iter++) {

    // COMPUTATION
    for (int i=0; i < parameters->getNumberParticles(); i++) {

      // POSITIONS
      for (int dim=0; dim < 2; dim++) {
        // initialise velocity
        (system->getParticle(i))->velocity()[dim] = 0.0;
        // initialise new positions with previous ones
        newParticles[i].position()[dim] =
          (system->getParticle(i))->position()[dim];
        // add self-propulsion
        selfPropulsion =
          #if CONTROLLED_DYNAMICS
          (1.0 - 2.0*system->getBiasingParameter()
            /3.0/parameters->getPersistenceLength())*
          #endif
          cos((system->getParticle(i))->orientation()[0] - dim*M_PI/2);
        (system->getParticle(i))->velocity()[dim] += selfPropulsion;
        newParticles[i].position()[dim] +=
          parameters->getTimeStep()*selfPropulsion;
        // add noise
        noise = (system->getRandomGenerator())->gauss_cutoff();
        (system->getParticle(i))->velocity()[dim] +=
          sqrt(2.0/3.0/parameters->getPersistenceLength())
          *noise;
        newParticles[i].position()[dim] +=
          sqrt(parameters->getTimeStep()
            *2.0/3.0/parameters->getPersistenceLength())
          *noise;
        // initialise force
        (system->getParticle(i))->force()[dim] = 0.0;
      }

      // ORIENTATIONS
      // initialise new orientation with previous one
      newParticles[i].orientation()[0] =
        (system->getParticle(i))->orientation()[0];
      // add noise
      newParticles[i].orientation()[0] +=
        sqrt(parameters->getTimeStep()*2.0/parameters->getPersistenceLength())
          *(system->getRandomGenerator())->gauss_cutoff();
      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      // initialise torque
      (system->getParticle(i))->torque()[0] = 0.0;
      #endif
    }

    // FORCES AND ALIGNING TORQUES
    ABP_WCA<System>(system); // compute forces
    #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
    aligningTorque<System>(system,
      [&system](int index) {
        return (system->getParticle(index))->orientation(); },
      [&system](int index) {
        return (system->getParticle(index))->torque(); }); // compute torques
    #endif

    for (int i=0; i < parameters->getNumberParticles(); i++) {
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->velocity()[dim] +=
          (system->getParticle(i))->force()[dim]
          /3.0/parameters->getPersistenceLength(); // add force
        newParticles[i].position()[dim] +=
          (system->getParticle(i))->force()[dim]
          *parameters->getTimeStep()/3.0/parameters->getPersistenceLength(); // add force displacement
      }
      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      newParticles[i].orientation()[0] +=
        (system->getParticle(i))->torque()[0]*parameters->getTimeStep(); // add torque rotation
      #endif
    }

    // HEUN'S SCHEME
    #if HEUN
    for (int i=0; i < parameters->getNumberParticles(); i++) {

      for (int dim=0; dim < 2; dim++) {
        // POSITIONS
        positions[2*i + dim] = (system->getParticle(i))->position()[dim]; // save initial position
        (system->getParticle(i))->position()[dim] =
          newParticles[i].position()[dim]; // integrate position as if using Euler's scheme
        // FORCES
        forces[2*i + dim] = (system->getParticle(i))->force()[dim]; // save computed force at initial position
        (system->getParticle(i))->force()[dim] = 0.0; // re-initialise force
      }

      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      // ORIENTATIONS
      orientations[i] = (sytem->getParticle(i))->orientation()[0]; // save initial orientation
      (sytem->getParticle(i))->orientation()[0] =
        newParticles[i].orientation()[0]; // integrate position as if using Euler's scheme
      // TORQUES
      torques[i] = (system->getParticle(i))->torque()[dim]; // save computed force at initial position
      (sytem->getParticle(i))->torque()[0] = 0.0; // re-initialise torque
      #endif
    }

    // FORCES AND ALIGNING TORQUES
    ABP_WCA<System>(system); // re-compute forces
    #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
    aligningTorque<System>(system,
      [&system](int index) {
        return (system->getParticle(index))->orientation(); },
      [&system](int index) {
        return (system->getParticle(index))->torque(); }); // re-compute torques
    #endif

    for (int i=0; i < parameters->getNumberParticles(); i++) {

      // CORRECTION TO INTERPARTICLE FORCE
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->velocity()[dim] +=
          ((system->getParticle(i))->force()[dim] - forces[2*i + dim])
          /3.0/parameters->getPersistenceLength()/2; // velocity
        newParticles[i].position()[dim] +=
          ((system->getParticle(i))->force()[dim] - forces[2*i + dim])
          *parameters->getTimeStep()/3.0/parameters->getPersistenceLength()/2; // position
        (system->getParticle(i))->force()[dim] =
          ((system->getParticle(i))->force()[dim] + forces[2*i + dim])/2; // force
      }

      // CORRECTION TO SELF-PROPULSION FORCE
      for (int dim=0; dim < 2; dim++) {
        selfPropulsionCorrection =
          #if CONTROLLED_DYNAMICS
          (1.0 - 2.0*system->getBiasingParameter()
            /3.0/parameters->getPersistenceLength())*
          #endif
          #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
          (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          - cos(orientations[2*i + dim] - dim*M_PI/2))
          #else
          (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          - cos((system->getParticle(i))->orientation()[0] - dim*M_PI/2))
          #endif
          /2;
        (system->getParticle(i))->velocity()[dim] +=
          selfPropulsionCorrection; // velocity
        newParticles[i].position()[dim] +=
          parameters->getTimeStep()*selfPropulsionCorrection; // position
      }

      // CORRECTION TO TORQUE
      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      newParticles[i].orientation()[0] +=
        (system->getParticle(i))->torque()[0] - torques[i])
        *parameters->getTimeStep()/2; // orientation
      (system->getParticle(i))->torque()[0] =
        ((system->getParticle(i))->torque()[0] + torques[i])/2; // torque
      #endif

      // RESET INITIAL POSITIONS AND ORIENTATION
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->position()[dim] = positions[2*i + dim]; // position
      }
      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      system->getParticle(i))->orientation()[0] = orientations[i]; // orientation
      #endif
    }
    #endif

    // SAVE AND COPY
    system->saveNewState(newParticles);
  }
}

void iterate_ABP_WCA(System0* system, int Niter) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

  Parameters* parameters = system->getParameters();

  std::vector<Particle> newParticles(parameters->getNumberParticles());

  double selfPropulsion; // self-propulsion force
  double noise; // noise realisation

  #if HEUN // HEUN'S SCHEME

  double selfPropulsionCorrection; // correction to the self-propulsion force

  std::vector<double> positions (2*parameters->getNumberParticles(), 0.0); // positions backup
  std::vector<double> forces (2*parameters->getNumberParticles(), 0.0); // forces backup
  #endif

  for (int iter=0; iter < Niter; iter++) {

    // COMPUTATION
    for (int i=0; i < parameters->getNumberParticles(); i++) {

      // POSITIONS
      for (int dim=0; dim < 2; dim++) {
        // initialise velocity
        (system->getParticle(i))->velocity()[dim] = 0.0;
        // initialise new positions with previous ones
        newParticles[i].position()[dim] =
          (system->getParticle(i))->position()[dim];
        // add self-propulsion
        selfPropulsion =
          parameters->getPropulsionVelocity()*
          cos((system->getParticle(i))->orientation()[0] - dim*M_PI/2);
        (system->getParticle(i))->velocity()[dim] += selfPropulsion;
        newParticles[i].position()[dim] +=
          parameters->getTimeStep()*selfPropulsion;
        // add noise
        noise = (system->getRandomGenerator())->gauss_cutoff();
        (system->getParticle(i))->velocity()[dim] +=
          sqrt(2.0*parameters->getTransDiffusivity())
          *noise;
        newParticles[i].position()[dim] +=
          sqrt(parameters->getTimeStep()
            *2.0*parameters->getTransDiffusivity())
          *noise;
        // initialise force
        (system->getParticle(i))->force()[dim] = 0.0;
      }

      // ORIENTATIONS
      // initialise new orientation with previous one
      newParticles[i].orientation()[0] =
        (system->getParticle(i))->orientation()[0];
      // add noise
      newParticles[i].orientation()[0] +=
        sqrt(parameters->getTimeStep()*2.0*parameters->getRotDiffusivity())
          *(system->getRandomGenerator())->gauss_cutoff();
    }

    // FORCES
    ABP_WCA<System0>(system); // compute forces

    for (int i=0; i < parameters->getNumberParticles(); i++) {
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->velocity()[dim] +=
          (system->getParticle(i))->force()[dim]
          *parameters->getPotentialParameter(); // add force
        newParticles[i].position()[dim] +=
          (system->getParticle(i))->force()[dim]
          *parameters->getTimeStep()*parameters->getPotentialParameter(); // add force displacement
      }
    }

    // HEUN'S SCHEME
    #if HEUN
    for (int i=0; i < parameters->getNumberParticles(); i++) {

      for (int dim=0; dim < 2; dim++) {
        // POSITIONS
        positions[2*i + dim] = (system->getParticle(i))->position()[dim]; // save initial position
        (system->getParticle(i))->position()[dim] =
          newParticles[i].position()[dim]; // integrate position as if using Euler's scheme
        // FORCES
        forces[2*i + dim] = (system->getParticle(i))->force()[dim]; // save computed force at initial position
        (system->getParticle(i))->force()[dim] = 0.0; // re-initialise force
      }
    }

    // FORCES
    ABP_WCA<System0>(system); // re-compute forces

    for (int i=0; i < parameters->getNumberParticles(); i++) {

      // CORRECTION TO INTERPARTICLE FORCE
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->velocity()[dim] +=
          ((system->getParticle(i))->force()[dim] - forces[2*i + dim])
          *parameters->getPotentialParameter(); // velocity
        newParticles[i].position()[dim] +=
          ((system->getParticle(i))->force()[dim] - forces[2*i + dim])
          *parameters->getTimeStep()*parameters->getPotentialParameter(); // position
        (system->getParticle(i))->force()[dim] =
          ((system->getParticle(i))->force()[dim] + forces[2*i + dim])/2; // force
      }

      // CORRECTION TO SELF-PROPULSION FORCE
      for (int dim=0; dim < 2; dim++) {
        selfPropulsionCorrection =
          parameters->getPropulsionVelocity()*
          (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          - cos((system->getParticle(i))->orientation()[0] - dim*M_PI/2))
          /2;
        (system->getParticle(i))->velocity()[dim] +=
          selfPropulsionCorrection; // velocity
        newParticles[i].position()[dim] +=
          parameters->getTimeStep()*selfPropulsionCorrection; // position
      }

      // RESET INITIAL POSITIONS
      for (int dim=0; dim < 2; dim++) {
        (system->getParticle(i))->position()[dim] = positions[2*i + dim]; // position
      }
    }
    #endif

    // SAVE AND COPY
    system->saveNewState(newParticles);
  }
}


/////////////////////////////////
// INTERACTING BROWNIAN ROTORS //
/////////////////////////////////

void iterate_rotors(Rotors* rotors, int Niter) {
  // Updates system to next step according to the dynamics of interacting
  // Brownian rotors.

  #if HEUN // HEUN'S SCHEME
  std::vector<double> orientations(rotors->getNumberParticles(), 0.0); // orientations backup
  std::vector<double> torques(rotors->getNumberParticles(), 0.0); // torques backup
  #endif

  for (int iter=0; iter < Niter; iter++) {

    // COMPUTATION
    for (int i=0; i < rotors->getNumberParticles(); i++) {
      // reset torques
      rotors->getTorque(i)[0] = 0.0;
      // add noise
      rotors->getOrientation(i)[0] +=
        sqrt(2.0*rotors->getRotDiffusivity()*rotors->getTimeStep())
        *(rotors->getRandomGenerator())->gauss_cutoff();
    }
    // compute aligning torques
    aligningTorque<Rotors>(rotors,
      [&rotors](int index) {
        return rotors->getOrientation(index); },
      [&rotors](int index) {
        return rotors->getTorque(index); }); // compute torques
    // add torque
    for (int i=0; i < rotors->getNumberParticles(); i++) {
      rotors->getOrientation(i)[0] +=
        rotors->getTorque(i)[0]*rotors->getTimeStep();
    }

    // HEUN'S SCHEME
    #if HEUN
    for (int i=0; i < rotors->getNumberParticles(); i++) {
      // save initial torques
      torques[i] = rotors->getTorque(i)[0];
      // reset torques
      rotors->getTorque(i)[0] = 0;
    }
    // re-compute aligning torques
    aligningTorque<Rotors>(rotors,
      [&rotors](int index) {
        return rotors->getOrientation(index); },
      [&rotors](int index) {
        return rotors->getTorque(index); }); // compute torques
    for (int i=0; i < rotors->getNumberParticles(); i++) {
      // correction to orientations
      rotors->getOrientation(i)[0] +=
        (rotors->getTorque(i)[0] - torques[i])*rotors->getTimeStep()
        /2;
      // correction to torques
      rotors->getTorque(i)[0] =
        (rotors->getTorque(i)[0] + torques[i])
        /2;
    }
    #endif

    // SAVE
    rotors->saveState();
  }
}
