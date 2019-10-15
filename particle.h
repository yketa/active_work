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
   *  (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)
   *
   *  Parameters are stored in a binary file with the following structure:
   *
   *  [HEADER (see System::System)]
   *  | (int) N | (double) lp | (double) phi | (double) L | (int) seed | (double) dt | (int) framesWork |
   *
   *  [BODY (see System::saveInitialState & System::saveNewState)] (all double)
   *  ||                FRAME 1                 || ... || FRAME framesWork ||                    || ...
   *  ||      PARTICLE 1     | ... | PARTICLE N || ... ||        ...       ||                    || ...
   *  ||   R   | ORIENTATION | ... |     ...    || ... ||        ...       || SUMMED ACTIVE WORK || ...
   *  || X | Y |    theta    | ... |     ...    || ... ||        ...       ||          W         || ...
   */

  public:

    // CONSTRUCTORS

    System(
      int N, double lp, double phi, int seed, double dt, std::string filename,
      int nWork, bool dump = true);
    System(
      int N, double lp, double phi, int seed, double dt, std::string filename,
      bool dump = true);

    // DESTRUCTORS

    ~System();

    // METHODS

    int getNumberParticles() const; // returns number of particles
    double getPersistenceLength() const; // returns dimensionless persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    int getRandomSeed() const; // returns random seed
    double getTimeStep() const; // returns time step
    std::string getOutputFile() const; // returns output file name

    rnd* getRandomGenerator(); // returns pointer to random generator
    Particle* getParticle(int index); // returns pointer to given particle
    std::ofstream* getOutputFileStream(); // returns pointer to output file stream

    double getDistance(int index1, int index2);
      // Returns distance between two particles in a given system.
    void WCA_potential(int index1, int index2, double *force);
      // Writes WCA force acting on particles[index1] by particles[index2] onto
      // `force`.

    void saveInitialState();
      // Saves initial state of particles to output file.
    void saveNewState(Particle *newParticles);
      // Saves new state of particles to output file then copy it.

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles in the system
    double const persistenceLength; // dimensionless persistence length
    double const packingFraction; // packing fraction
    double const systemSize; // system size
    int const randomSeed; // random seed
    double const timeStep; // time step
    std::string const outputFile; // output file name

    int const framesWork; // number of frames on which to sum the active work before dumping
      // taken roughly equal to lp/dt
    bool const dumpParticles; // dump positions and orientations to output file

    rnd randomGenerator; // random number generator
    std::vector<Particle> particles; // vector of particles
    std::ofstream outputFileStream; // output file stream
      // WARNING: Content of file is erased.

    int dumpFrame; // index of last frame dumped
    double workSum; // sum of the active works since the last dump

};

#endif
