#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <string>
#include <fstream>

#include "maths.hpp"
#include "param.hpp"
#include "write.hpp"

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

    double* force(); // returns pointer to force

  private:

    // ATTRIBUTES

    double r[2]; // position (2D)
    double theta; // orientation

    double f[2]; // force exerted on particle (2D)

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

    int getNumberBoxes(); // return number of boxes in each dimension
    std::vector<int> getCell(int const &index); // return vector of indexes in cell

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
   *  ||             FRAME period*1             || ... || FRAME period*framesWork ||                                                                                 || ...
   *  ||      PARTICLE 1     | ... | PARTICLE N || ... ||        ...              ||                                                                                 || ...
   *  ||   R   | ORIENTATION | ... |     ...    || ... ||        ...              || ACTIVE WORK | ACTIVE WORK (FORCE) | ACTIVE WORK (ORIENTATION) | ORDER PARAMETER || ...
   *  || X | Y |    theta    | ... |     ...    || ... ||        ...              ||      W      |          Wp         |             Wo            |        nu       || ...
   */

  public:

    // CONSTRUCTORS

    System();
    System(
      Parameters* parameters, int seed = 0, std::string filename = "",
      int nWork = 0, bool dump = true, int period = 1);
    System(
      System* system, int seed = 0, std::string filename = "",
      int nWork = 0, bool dump = true, int period = 1);

    // DESTRUCTORS

    ~System();

    // METHODS

    Parameters* getParameters(); // returns pointer to class of parameters

    int getNumberParticles() const; // returns number of particles
    double getPersistenceLength() const; // returns dimensionless persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    double getTimeStep() const; // returns time step

    int getRandomSeed() const; // returns random seed
    rnd* getRandomGenerator(); // returns pointer to random generator

    Particle* getParticle(int index); // returns pointer to given particle
    std::vector<Particle> getParticles(); // returns vector of particles

    CellList* getCellList(); // returns pointer to CellList object

    std::string getOutputFile() const; // returns output file name

    double getWork(); // returns last computed normalised rate of active work
    double getWorkForce(); // returns last computed force part of the normalised rate of active work
    double getWorkOrientation(); // returns last computed orientation part of the normalised rate of active work
    double getOrder(); // returns last computed order parameter
    // NOTE: All these quantities are computed every framesWork*dumpPeriod iterations.

    double diffPeriodic(double const& x1, double const& x2);
      // Returns distance between two pointson a line taking into account periodic
      // boundary condition.
    double getDistance(int const& index1, int const& index2);
      // Returns distance between two particles in a given system.

    void WCA_potential(int const& index1, int const& index2,
      std::vector<Particle>& newParticles);
      // Compute WCA forces between particles[index1] and particles[index2],
      // add to particles[index1].force() and particles[index2].force(), and
      // increments positions in particles[index1].position() and
      // particles[index2].position().

    void copyParticles(std::vector<Particle>& newParticles);
      // Replace vector of particles by newParticles.
    void copyParticles(System* system);
      // Replace vector of particles by the one from system.

    void saveInitialState();
      // Saves initial state of particles to output file.
    void saveNewState(std::vector<Particle>& newParticles);
      // Saves new state of particles to output file then copy it.

  private:

    // ATTRIBUTES

    Parameters param; // class of simulation parameters

    int const randomSeed; // random seed
    rnd randomGenerator; // random number generator

    std::vector<Particle> particles; // vector of particles

    CellList cellList; // cell list

    Output output; // output class

    int const framesWork; // number of frames on which to sum the active work before dumping
      // taken roughly equal to lp/dt
    bool const dumpParticles; // dump positions and orientations to output file
    int const dumpPeriod; // period of dumping of positions and orientations in number of frames

    int dumpFrame; // index of last frame dumped
    double workSum[2]; // sum of the active works since the last dump
    double workForceSum[2]; // sum of the force part of the active works since the last dump
    double workOrientationSum[2]; // sum of the orientation part of the active works since the last dump
    double orderSum[2]; // sum of order parameter norm since the last dump

};

#endif
