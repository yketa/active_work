#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <string>
#include <fstream>

#include "maths.hpp"
#include "write.hpp"

/////////////
// CLASSES //
/////////////

class Particle;
class CellList;
class Parameters;
class System;


/*  PARTICLE
 *  --------
 *  Store diameter, positions and orientation of a given particle.
 */

class Particle {

  public:

    // CONSTRUCTORS

    Particle();
    Particle(double x, double y, double ang);
    Particle(double x, double y, double ang, double d);

    // METHODS

    double* position(); // returns pointer to position
    double* orientation(); // returns pointer to orientation

    double* diameter(); // returns pointer to diameter

    double* force(); // returns pointer to force

  private:

    // ATTRIBUTES

    double r[2]; // position (2D)
    double theta; // orientation

    double sigma; // diameter

    double f[2]; // force exerted on particle (2D)

};


/*  CELL LIST
 *  ---------
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

    template<class SystemClass> void initialise(
      SystemClass* system, double const& rcut) {
      // Initialise cell list.

      // parameters of cell list
      cutOff = rcut;
      numberBoxes = std::max((int) (system->getSystemSize()/rcut), 1);
      sizeBox = system->getSystemSize()/numberBoxes;

      // set size of cell list
      for (int i=0; i < pow(numberBoxes, 2); i++) {
        cellList.push_back(std::vector<int>());
      }

      // set number of neighbours to explore
      if ( numberBoxes == 1 ) { dmin = 1; }
      else if ( numberBoxes == 2 ) { dmin = 0; }
      else { dmin = -1; }

      // put particles in the boxes
      update(system);
    }

    template<class SystemClass> void update(SystemClass* system) {
      // Put particles in the cell list.

      // flush old lists
      for (int i=0; i < (int) cellList.size(); i++) {
        cellList[i].clear();
      }

      // create new lists
      for (int i=0; i < system->getNumberParticles(); i++) {
        cellList[index(system->getParticle(i))].push_back(i); // particles are in increasing order of indexes
      }
    }

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


/*  PARAMETERS
 *  ----------
 *  Store parameters relative to a system.
 */

class Parameters {

  public:

    // CONSTRUCTORS

    Parameters();
    Parameters( // using custom dimensionless parameters relations
      int N, double lp, double phi, double dt);
    Parameters( // defining all parameters independently except system size
      int N, double epsilon, double v0, double D, double Dr, double phi,
      double dt);
    Parameters( // defining all parameters independently
      int N, double epsilon, double v0, double D, double Dr, double phi,
      double L, double dt);
    Parameters( // copy other class, changing system size
      Parameters* parameters);

    // METHODS

    int getNumberParticles() const; // returns number of particles in the system
    double getPotentialParameter() const; // returns coefficient parameter of potential
    double getPropulsionVelocity() const; // returns self-propulsion velocity
    double getTransDiffusivity() const; // returns translational diffusivity
    double getRotDiffusivity() const; // returns rotational diffusivity
    double getPersistenceLength() const; // returns persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    double getTimeStep() const; // returns time step

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles in the system
    double const potentialParameter; // coefficient parameter of potential
    double const propulsionVelocity; // self-propulsion velocity
    double const transDiffusivity; // translational diffusivity
    double const rotDiffusivity; // rotational diffusivity
    double const persistenceLength; // persistence length
    double const packingFraction; // packing fraction
    double const systemSize; // system size
    double const timeStep; // time step

};


/*  SYSTEM
 *  ------
 *  Store physical and integration parameter.
 *  Access to distance and potentials.
 *  Save system state to output file.
 *  Using custom dimensionless parameters relations.
 */

class System {
  /*  Contains all the details to simulate a system of active Brownian
   *  particles, with dimensionless parameters taken from Nemoto et al., PRE 99
   *  022605 (2019).
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

    void setBiasingParameter(double sValue); // set new biasing parameter
    double getBiasingParameter(); // returns biasing parameter

    int getDump(); // returns number of frames dumped since last reset
    void resetDump();
      // Reset time-extensive quantities over trajectory.
    void copyDump(System* system);
      // Copy dumps from other system.
      // WARNING: This also copies the index of last frame dumped. Consistency
      //          has to be checked.

    double getWork(); // returns last computed normalised rate of active work
    double getWorkForce(); // returns last computed force part of the normalised rate of active work
    double getWorkOrientation(); // returns last computed orientation part of the normalised rate of active work
    double getOrder(); // returns last computed averaged integrated order parameter
    // NOTE: All these quantities are computed every framesWork*dumpPeriod iterations.

    double getTotalWork(); // returns computed active work since last reset
    double getTotalWorkForce(); // returns computed force part of the active work since last rest
    double getTotalWorkOrientation(); // returns computed orientation part of the active work since last reset
    double getTotalOrder(); // returns computed integrated order parameter since last reset
    // NOTE: All these quantities are updated every framesWork*dumpPeriod iterations.
    //       All these quantities are extensive in time since last reset.

    #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
    void setTorqueParameter(double g); // set new torque parameter
    double getTorqueParameter(); // returns torque parameter

    double getTorqueIntegral1(); // returns last computed normalised first integral in the expression of the modified active work
    double getTorqueIntegral2(); // returns last computed normalised second integral in the expression of the modified active work
    // NOTE: All these quantities are computed every framesWork*dumpPeriod iterations.

    double getTotalTorqueIntegral1(); // returns computed normalised first integral in the expression of the modified active work since last reset
    double getTotalTorqueIntegral2(); // returns computed normalised second integral in the expression of the modified active work since last reset
    // NOTE: All these quantities are computed every framesWork*dumpPeriod iterations.
    //       All these quantities are extensive in time steps since last reset.
    #endif

    double diffPeriodic(double const& x1, double const& x2);
      // Returns distance between two pointson a line taking into account periodic
      // boundary condition.
    double getDistance(int const& index1, int const& index2);
      // Returns distance between two particles in a given system.

    void WCA_force(int const& index1, int const& index2,
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

    double biasingParameter; // biasing parameter [cloning algorithm]

    int dumpFrame; // number of frames dumped since last reset
    // Quantities
    // (0): sum of quantity since last dump
    // (1): normalised quantity over last dump period
    // (2): time-extensive quantity over trajectory since last reset
    double workSum[3]; // active work
    double workForceSum[3]; //force part of the active work
    double workOrientationSum[3]; // orientation part of the active work
    double orderSum[3]; // integrated order parameter norm (in units of the time step)

    #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
    double torqueParameter; // torque parameter [cloning algorithm]

    // Quantities
    // (0): integral since last dump (in units of the time step)
    // (1): normalised integral over last dump period
    // (2): integral over trajectory since last reset (in units of the time step)
    double torqueIntegral1[3]; // first integral in the expression of the modified active work
    double torqueIntegral2[3]; // second integral in the expression of the modified active work
    #endif

};


/*  SYSTEM0
 *  -------
 *  Store physical and integration parameter.
 *  Access to distance and potentials.
 *  Save system state to output file.
 *  Using all free parameters in the ABP model.
 */

class System0 {
  /*  Contains all the details to simulate a system of active Brownian
   *  particles.
   *  (see https://yketa.github.io/DAMTP_2019_Wiki/#Active%20Brownian%20particles)
   *
   *  Parameters are stored in a binary file with the following structure:
   *
   *  [HEADER (see System::System)]
   *  | (int) N | (double) epsilon | (double) v0 | (double) D | (double) Dr | (double) lp | (double) phi | (double) L | (int) seed | (double) dt | (int) framesWork | (bool) dump | (int) period |
   *  || PARTICLE 1 | ... | PARTICLE N ||
   *  ||  diameter  | ... |  diameter  ||
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

    System0();
    System0(
      Parameters* parameters, int seed = 0, std::string filename = "",
      int nWork = 0, bool dump = true, int period = 1);
    System0(
      Parameters* parameters, std::vector<double>& diameters, int seed = 0,
      std::string filename = "", int nWork = 0, bool dump = true,
      int period = 1);
    System0(
      System0* system, int seed = 0, std::string filename = "",
      int nWork = 0, bool dump = true, int period = 1);
    System0(
      System0* system, std::vector<double>& diameters, int seed = 0,
      std::string filename = "", int nWork = 0, bool dump = true,
      int period = 1);

    // DESTRUCTORS

    ~System0();

    // METHODS

    Parameters* getParameters(); // returns pointer to class of parameters

    int getNumberParticles() const; // returns number of particles
    double getPotentialParameter() const; // returns coefficient parameter of potential
    double getPropulsionVelocity() const; // returns self-propulsion velocity
    double getTransDiffusivity() const; // returns translational diffusivity
    double getRotDiffusivity() const; // returns rotational diffusivity
    double getPersistenceLength() const; // returns persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    double getTimeStep() const; // returns time step

    int getRandomSeed() const; // returns random seed
    rnd* getRandomGenerator(); // returns pointer to random generator

    Particle* getParticle(int index); // returns pointer to given particle
    std::vector<Particle> getParticles(); // returns vector of particles

    CellList* getCellList(); // returns pointer to CellList object

    std::string getOutputFile() const; // returns output file name

    int getDump(); // returns number of frames dumped since last reset
    void resetDump();
      // Reset time-extensive quantities over trajectory.
    void copyDump(System0* system);
      // Copy dumps from other system.
      // WARNING: This also copies the index of last frame dumped. Consistency
      //          has to be checked.

    double getWork(); // returns last computed normalised rate of active work
    double getWorkForce(); // returns last computed force part of the normalised rate of active work
    double getWorkOrientation(); // returns last computed orientation part of the normalised rate of active work
    double getOrder(); // returns last computed averaged integrated order parameter
    // NOTE: All these quantities are computed every framesWork*dumpPeriod iterations.

    double getTotalWork(); // returns computed active work since last reset
    double getTotalWorkForce(); // returns computed force part of the active work since last rest
    double getTotalWorkOrientation(); // returns computed orientation part of the active work since last reset
    double getTotalOrder(); // returns computed integrated order parameter since last reset
    // NOTE: All these quantities are updated every framesWork*dumpPeriod iterations.
    //       All these quantities are extensive in time since last reset.

    double diffPeriodic(double const& x1, double const& x2);
      // Returns distance between two pointson a line taking into account periodic
      // boundary condition.
    double getDistance(int const& index1, int const& index2);
      // Returns distance between two particles in a given system.

    void WCA_force(int const& index1, int const& index2,
      std::vector<Particle>& newParticles);
      // Compute WCA forces between particles[index1] and particles[index2],
      // add to particles[index1].force() and particles[index2].force(), and
      // increments positions in particles[index1].position() and
      // particles[index2].position().

    void copyParticles(std::vector<Particle>& newParticles);
      // Replace vector of particles by newParticles.
    void copyParticles(System0* system);
      // Replace vector of particles by the one from system.
    void copyParticles(System0* system, std::vector<double>& diameters);
      // Replace vector of particles by the one from system.
      // Replace particles' diameters by the ones in diameters.

    void setDiameters(std::vector<double>& diameters);
      // Set all diameters and re-define system size.

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

    int dumpFrame; // number of frames dumped since last reset
    // Quantities
    // (0): sum of quantity since last dump
    // (1): normalised quantity over last dump period
    // (2): time-extensive quantity over trajectory since last reset
    double workSum[3]; // active work
    double workForceSum[3]; //force part of the active work
    double workOrientationSum[3]; // orientation part of the active work
    double orderSum[3]; // integrated order parameter norm (in units of the time step)

};


////////////////
// PROTOTYPES //
////////////////

double getGlobalPhase(std::vector<Particle>& particles);
  // Returns global phase.

std::vector<double> getOrderParameter(std::vector<Particle>& particles);
  // Returns order parameter.

double getOrderParameterNorm(std::vector<Particle>& particles);
  // Returns order parameter norm.

void _WCA_force(
  System* system, int const& index1, int const& index2, double* force);
  // Writes to `force' the force deriving from the WCA potential between
  // particles `index1' and `index2'.

void _WCA_force(
  System0* system, int const& index1, int const& index2, double* force);
  // Writes to `force' the force deriving from the WCA potential between
  // particles `index1' and `index2'.

template<class SystemClass> double WCA_potential(
  SystemClass* system, double const& index1, double const& index2) {
  // Returns WCA potential between two particles in a given system.

  double potential = 0;

  double dist = system->getDistance(index1, index2); // dimensionless distance between particles
  double sigma =
    ((system->getParticle(index1))->diameter()[0]
    + (system->getParticle(index2))->diameter()[0])/2; // equivalent diameter

  if (dist/sigma < pow(2, 1./6.)) { // distance lower than cut-off

    // compute potential
    potential = (system->getParameters())->getPotentialParameter()
      *(4*(1/pow(dist/sigma, 12) - 1/pow(dist/sigma, 6)) - 1);
  }

  return potential;
}

template<class SystemClass> double _diffPeriodic(
  SystemClass* system, double const& x1, double const& x2) {
  // Returns algebraic distance from `x1' to `x2' on a line taking into account
  // periodic boundary condition of the system.

    return algDistPeriod(x1, x2, system->getSystemSize());
}

template<class SystemClass> double _getDistance(
  SystemClass* system, int const& index1, int const& index2) {
  // Returns distance between two particles in a given system.

  return dist2DPeriod(
    (system->getParticle(index1))->position(),
    (system->getParticle(index2))->position(),
    system->getSystemSize());
}

#endif
