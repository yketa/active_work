#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <fstream>

#include "maths.h"

/////////////
// CLASSES //
/////////////

class Particle;
class System;
class CellList;


/*  PARTICLE
 *  --------
 *  Store positions and orientation of a given particle.
 */

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


/*  CELL LIST
 *  --------
 *  Speed up computation by storing closest neighbours.
 */

class CellList {

  public:

    // CONSTRUCTORS

    CellList();

    // DESTRUCTORS

    ~CellList();

    // METHODS

    void initialise(System *system, double const& rcut);
      // Initialise cell list.

    void update(System *system);
      // Put particles in the cell list.

    int index(Particle *particle);
      // Index of the box corresponding to a given particle.

    std::vector<int> getNeighbours(Particle *particle);
      // Returns vector of indexes of neighbouring particles.

  private:

    // ATTRIBUTES

    double cutOff; // cut-off radius of the interactions

    int numberBoxes; // number of boxes in each dimension
    double sizeBox; // size of each box
    int dmin; // trick to avoid putting too much neighbours when rcut is large

    std::vector<std::vector<int>> cellList; // cells with indexes of particles

};


/*  SYSTEM
 *  --------
 *  Store physical and integration parameter.
 *  Access to distance and potentials.
 *  Save system state to output file.
 */

class System {
  /*  Contains all the details to simulate a system of active particles.
   *  (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)
   *
   *  Parameters are stored in a binary file with the following structure:
   *
   *  [HEADER (see System::System)]
   *  | (int) N | (double) lp | (double) phi | (double) L | (int) seed | (double) dt | (int) framesWork | (bool) dump | (int) period |
   *
   *  [INITIAL FRAME]
   *  ||                FRAME 0                 ||
   *  ||      PARTICLE 1     | ... | PARTICLE N ||
   *  ||   R   | ORIENTATION | ... |     ...    ||
   *  || X | Y |    theta    | ... |     ...    ||
   *
   *  [BODY (see System::saveInitialState & System::saveNewState)] (all double)
   *  ||             FRAME period*1             || ... || FRAME period*framesWork ||                    || ...
   *  ||      PARTICLE 1     | ... | PARTICLE N || ... ||        ...              ||                    || ...
   *  ||   R   | ORIENTATION | ... |     ...    || ... ||        ...              || SUMMED ACTIVE WORK || ...
   *  || X | Y |    theta    | ... |     ...    || ... ||        ...              ||          W         || ...
   */

  public:

    // CONSTRUCTORS

    System(
      int N, double lp, double phi, int seed, double dt, std::string filename,
      int nWork = 0, bool dump = true, int period = 1);

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

    std::vector<int> getNeighbours(int const& index);
      // Returns vector of indexes of neighbouring particles.
    double diffPeriodic(double const& x1, double const& x2);
      // Returns distance between two pointson a line taking into account periodic
      // boundary condition.
    double getDistance(int const& index1, int const& index2);
      // Returns distance between two particles in a given system.
    void WCA_potential(int const& index1, int const& index2,
      double *force);
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
    int const dumpPeriod; // period of dumping of positions and orientations in number of frames

    rnd randomGenerator; // random number generator
    std::vector<Particle> particles; // vector of particles
    std::ofstream outputFileStream; // output file stream
      // WARNING: Content of file is erased.
    CellList cellList; // cell list

    int dumpFrame; // index of last frame dumped
    double workSum; // sum of the active works since the last dump

};

#endif
