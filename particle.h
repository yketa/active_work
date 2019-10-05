#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <fstream>

#include "maths.h"

/////////////
// CLASSES //
/////////////


// PARTICLE

class Particle {

  public:

    // CONSTRUCTORS

    Particle();
    Particle(double x, double y, double ang);

    // METHODS

    double* position(); // returns pointer to position
    double* orientation(); // returns pointer to orientation

    void copy(Particle *particle);
      // Copy content of other particle.

  private:

    // ATTRIBUTES

    double r[2]; // position (2D)
    double theta; // orientation

};


// SYSTEM

class System {
  /*  Contains all the details to simulate a system of active particles.
   *  (see Phys. Rev. E, 99, 022605 (2019) for details)
   *
   *  Parameters are stored in a binary file with the following structure:

   *  [HEADER (see System::saveInitialState)]
   *  | (int) N | (double) lp | (double) L | (int) seed |
   *
   *  [BODY (see System::saveNewState)] (all double)
   *  ||                              FRAME                               || ...
   *  ||      ||              PARTICLE 1               | ... | PARTICLE N || ...
   *  || STEP || WRAPPED R | UNWRAPPED R | ORIENTATION | ... |     ...    || ...
   *  ||  dt  ||  x  |  y  |   x  |  y   |    theta    | ... |     ...    || ...
   */

  public:

    // CONSTRUCTORS

    System(int N, double lp, double L, int seed, std::string filename);

    // METHODS

    int getNumberParticles() const; // returns number of particles
    double getPersistenceLength() const; // returns dimensionless persistence length
    double getSystemSize() const; // returns system size
    int getRandomSeed() const; // returns random seed
    std::string getOutputFile() const; // returns output file name

    rnd* getRandomGenerator(); // returns pointer to random generator
    Particle* getParticle(int index); // returns pointer to given particle
    std::ofstream* getOutputFileStream(); // returns pointer to output file stream
      // WARNING: Content of file is erased.

    double getDistance(int index1, int index2);
      // Returns distance between two particles in a given system.
    void WCA_potential(int index1, int index2, double *force);
      // Writes WCA force acting on particles[index1] by particles[index2] onto
      // `force`.

    void saveInitialState();
      // Saves initial state of particles to output file.
    void saveNewState(Particle *newParticles, double const& dt);
      // Saves new state of particles to output file then copy it.

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles in the system
    double const persistenceLength; // dimensionless persistence length
    double const systemSize; // system size
    int const randomSeed; // random seed
    std::string const outputFile; // output file name

    rnd randomGenerator; // random number generator
    std::vector<Particle> particles; // vector of particles
    std::ofstream outputFileStream; // output file stream

    std::vector<std::vector<double> > increments; // vector of increments for unwrapped coordinates

};

#endif
