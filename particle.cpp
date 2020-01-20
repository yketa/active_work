#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include "dat.hpp"
#include "particle.hpp"
#include "maths.hpp"

/////////////
// CLASSES //
/////////////

/************
 * PARTICLE *
 ************/

// CONSTRUCTORS

Particle::Particle() : r {0, 0}, theta (0), v {0, 0}, sigma (1), f {0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , gamma (0)
  #endif
  {}
Particle::Particle(double x, double y, double ang) :
  r {x, y}, theta (ang), v {0, 0}, sigma (1), f {0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , gamma (0)
  #endif
  {}
Particle::Particle(double x, double y, double ang, double d) :
  r {x, y}, theta (ang), v {0, 0}, sigma (d), f {0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , gamma (0)
  #endif
  {}

// METHODS

double* Particle::position() { return &r[0]; } // returns pointer to position
double* Particle::orientation() { return &theta; } // returns pointer to orientation
double* Particle::velocity() { return &v[0]; } // returns pointer to velocity

double* Particle::diameter() { return &sigma; } // returns pointer to diameter

double* Particle::force() { return &f[0]; }; // returns pointer to force

#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
double* Particle::torque() { return &gamma; } // returns pointer to aligning torque
#endif


/*************
 * CELL LIST *
 *************/

// CONSTRUCTORS

CellList::CellList() {}

// DESTRUCTORS

CellList::~CellList() {}

// METHODS

int CellList::getNumberBoxes() { return numberBoxes; } // return number of boxes in each dimension
std::vector<int>* CellList::getCell(int const &index) {
  return &cellList[index]; } // return pointer to vector of indexes in cell

int CellList::index(Particle *particle) {
  // Index of the box corresponding to a given particle.

  int x = (int) ((particle->position())[0]/sizeBox);
  int y = (int) ((particle->position())[1]/sizeBox);
  return (x == numberBoxes ? 0 : x) + numberBoxes*(y == numberBoxes ? 0 : y);
}

std::vector<int> CellList::getNeighbours(Particle *particle) {
  // Returns vector of indexes of neighbouring particles.

  std::vector<int> neighbours; // vector of neighbouring particles

  int indexParticle = index(particle);
  int x = indexParticle%numberBoxes;
  int y = indexParticle/numberBoxes;

  int neighbourIndex;
  for (int dx=dmin; dx < 2; dx++) {
    for (int dy=dmin; dy < 2; dy++) {
      neighbourIndex =
        (numberBoxes + (x + dx))%numberBoxes
          + numberBoxes*((numberBoxes + (y + dy))%numberBoxes); // index of neighbouring cell
      neighbours.insert(
        std::end(neighbours),
        std::begin(cellList[neighbourIndex]),
        std::end(cellList[neighbourIndex])); // add particle indexes of neighbouring cell
    }
  }

  return neighbours;
}


/**************
 * PARAMETERS *
 **************/

// CONSTRUCTORS

Parameters::Parameters() :
  numberParticles(0), potentialParameter(0), propulsionVelocity(0),
    transDiffusivity(0), rotDiffusivity(0), persistenceLength(0),
    packingFraction(0), systemSize(0), timeStep(0) {}

Parameters::Parameters(int N, double lp, double phi, double dt) :
  numberParticles(N), potentialParameter(1.0), propulsionVelocity(1.0),
    transDiffusivity(1.0/(3.0*lp)), rotDiffusivity(1.0/lp),
    persistenceLength(lp), packingFraction(phi),
    systemSize(sqrt(M_PI*N/phi)/2.0), timeStep(dt) {}

Parameters::Parameters(
  int N, double epsilon, double v0, double D, double Dr, double phi,
  double dt) :
  numberParticles(N), potentialParameter(epsilon), propulsionVelocity(v0),
    transDiffusivity(D), rotDiffusivity(Dr), persistenceLength(v0/Dr),
    packingFraction(phi), systemSize(sqrt(M_PI*N/phi)/2.0), timeStep(dt) {}

Parameters::Parameters(
  int N, double epsilon, double v0, double D, double Dr, double phi, double L,
  double dt) :
  numberParticles(N), potentialParameter(epsilon), propulsionVelocity(v0),
    transDiffusivity(D), rotDiffusivity(Dr), persistenceLength(v0/Dr),
    packingFraction(phi), systemSize(L), timeStep(dt) {}

Parameters::Parameters(Parameters* parameters) :
  numberParticles(parameters->getNumberParticles()),
  potentialParameter(parameters->getPotentialParameter()),
  propulsionVelocity(parameters->getPropulsionVelocity()),
  transDiffusivity(parameters->getTransDiffusivity()),
  rotDiffusivity(parameters->getRotDiffusivity()),
  persistenceLength(parameters->getPersistenceLength()),
  packingFraction(parameters->getPackingFraction()),
  systemSize(parameters->getSystemSize()),
  timeStep(parameters->getTimeStep()) {}

// METHODS

int Parameters::getNumberParticles() { return numberParticles; }
double Parameters::getPotentialParameter() { return potentialParameter; }
double Parameters::getPropulsionVelocity() { return propulsionVelocity; }
double Parameters::getTransDiffusivity() { return transDiffusivity; }
double Parameters::getRotDiffusivity() { return rotDiffusivity; }
double Parameters::getPersistenceLength() { return persistenceLength; }
double Parameters::getPackingFraction() {return packingFraction; }
double Parameters::getSystemSize() { return systemSize; }
double Parameters::getTimeStep() { return timeStep; }


/**********
 * SYSTEM *
 **********/

// CONSTRUCTORS

System::System() :
  param(new Parameters()),
  randomSeed(0), randomGenerator(),
  particles(0),
  cellList(),
  output(""), velocitiesDumps(),
  framesWork(0), dumpParticles(0), dumpPeriod(0),
  biasingParameter(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , torqueParameter(0),
  torqueIntegral1 {0, 0, 0}, torqueIntegral2 {0, 0, 0}
  #endif
  {}

System::System(
  Parameters* parameters, int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(parameters),
  randomSeed(seed), randomGenerator(),
  particles(parameters->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(parameters->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    parameters->getPersistenceLength()/(parameters->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  biasingParameter(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , torqueParameter(0),
  torqueIntegral1 {0, 0, 0}, torqueIntegral2 {0, 0, 0}
  #endif
  {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // put particles on a grid with random orientation
  int gridSize = ceil(sqrt(getNumberParticles())); // size of the grid on which to put the particles
  double gridSpacing = getSystemSize()/gridSize;
  for (int i=0; i < getNumberParticles(); i++) { // loop over particles
    // position on the grid
    particles[i].position()[0] = (i%gridSize)*gridSpacing;
    particles[i].position()[1] = (i/gridSize)*gridSpacing;
    // random orientation
    particles[i].orientation()[0] = 2*M_PI*randomGenerator.random01();
  }

  // initialise cell list
  cellList.initialise<System>(this, pow(2., 1./6.));
}

System::System(
  System* system, int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(system->getParameters()),
  randomSeed(seed), randomGenerator(),
  particles(system->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(system->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    system->getPersistenceLength()/(system->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  biasingParameter(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , torqueParameter(0),
  torqueIntegral1 {0, 0, 0}, torqueIntegral2 {0, 0, 0}
  #endif
  {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // initialise cell list
  cellList.initialise<System>(this, pow(2., 1./6.));
  // copy particles and update cell list
  copyParticles(system);
  // copy dumps
  copyDump(system);
}

System::System(
  std::string inputFilename, int inputFrame, double dt,
  int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(new Parameters()),
  randomSeed(seed), randomGenerator(),
  particles(0),
  cellList(),
  output(filename), velocitiesDumps(0),
  framesWork(nWork), dumpParticles(dump), dumpPeriod(period),
  biasingParameter(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0}
  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  , torqueParameter(0),
  torqueIntegral1 {0, 0, 0}, torqueIntegral2 {0, 0, 0}
  #endif
  {

  // load data
  Dat inputDat(inputFilename, false); // data object
  param = Parameters( // set parameters
    inputDat.getNumberParticles(),
    inputDat.getPersistenceLength(),
    inputDat.getPackingFraction(),
    dt > 0 ? dt : inputDat.getTimeStep());
  particles.resize(getNumberParticles()); // resize vector of particles
  velocitiesDumps.resize(2*getNumberParticles()); // resize vector of locations of velocity dumps
  for (int i=0; i < getNumberParticles(); i++) {
    // set positions
    for (int dim=0; dim < 2; dim++) {
      particles[i].position()[dim] = inputDat.getPosition(inputFrame, i, dim);
    }
    // set orientation
    particles[i].orientation()[0] = inputDat.getOrientation(inputFrame, i);
  }

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // initialise cell list
  cellList.initialise<System>(this, pow(2., 1./6.));
}

// DESTRUCTORS

System::~System() {}

// METHODS

Parameters* System::getParameters() { return &param; }

int System::getNumberParticles() {
  return param.getNumberParticles(); }
double System::getPersistenceLength() {
  return param.getPersistenceLength(); }
double System::getPackingFraction() {
  return param.getPackingFraction(); }
double System::getSystemSize() {
  return param.getSystemSize(); }
double System::getTimeStep() {
  return param.getTimeStep(); }

void System::setTimeStep(double dt) {
  param = Parameters(
    getNumberParticles(), getPersistenceLength(), getPackingFraction(), dt);
}

int System::getRandomSeed() const { return randomSeed; }
rnd* System::getRandomGenerator() { return &randomGenerator; }

Particle* System::getParticle(int index) { return &(particles[index]); }
std::vector<Particle> System::getParticles() { return particles; }

CellList* System::getCellList() { return &cellList; }

std::string System::getOutputFile() const { return output.getOutputFile(); }

void System::setBiasingParameter(double sValue) { biasingParameter = sValue; }
double System::getBiasingParameter() { return biasingParameter; }

int System::getDump() { return dumpFrame; }

void System::resetDump() {
  // Reset time-extensive quantities over trajectory.

  dumpFrame = 0;

  workSum[0] = 0;
  workForceSum[0] = 0;
  workOrientationSum[0] = 0;
  orderSum[0] = 0;

  workSum[2] = 0;
  workForceSum[2] = 0;
  workOrientationSum[2] = 0;
  orderSum[2] = 0;

  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  torqueIntegral1[0] = 0;
  torqueIntegral2[0] = 0;

  torqueIntegral1[2] = 0;
  torqueIntegral2[2] = 0;
  #endif
}

void System::copyDump(System* system) {
  // Copy dumps from other system.
  // WARNING: This also copies the index of last frame dumped. Consistency
  //          has to be checked.

  dumpFrame = system->getDump();

  workSum[2] = system->getTotalWork();
  workForceSum[2] = system->getTotalWorkForce();
  workOrientationSum[2] = system->getTotalWorkOrientation();
  orderSum[2] = system->getTotalOrder();

  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  torqueIntegral1[2] = system->getTotalTorqueIntegral1();
  torqueIntegral2[2] = system->getTotalTorqueIntegral2();
  #endif
}

double System::getWork() { return workSum[1]; }
double System::getWorkForce() { return workForceSum[1]; }
double System::getWorkOrientation() { return workOrientationSum[1]; }
double System::getOrder() { return orderSum[1]; }

double System::getTotalWork() { return workSum[2]; }
double System::getTotalWorkForce() { return workForceSum[2]; }
double System::getTotalWorkOrientation() { return workOrientationSum[2]; }
double System::getTotalOrder() { return orderSum[2]; }

#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
void System::setTorqueParameter(double g) { torqueParameter = g; }
double System::getTorqueParameter() { return torqueParameter; }

double System::getTorqueIntegral1() { return torqueIntegral1[1]; }
double System::getTorqueIntegral2() { return torqueIntegral2[1]; }

double System::getTotalTorqueIntegral1() { return torqueIntegral1[2]; }
double System::getTotalTorqueIntegral2() { return torqueIntegral2[2]; }
#endif

double System::diffPeriodic(double const& x1, double const& x2) {
  // Returns algebraic distance from `x1' to `x2' on a line taking into account
  // periodic boundary condition of the system.

  return _diffPeriodic<System>(this, x1, x2);
}

double System::getDistance(int const& index1, int const& index2) {
  // Returns distance between two particles in a given system.

  return _getDistance<System>(this, index1, index2);
}

void System::WCA_force(int const& index1, int const& index2) {
  // Compute WCA forces between particles[index1] and particles[index2],
  // and add to particles[index1].force() and particles[index2].force().

  if ( index1 != index2 ) { // only consider different particles

    double force[2];
    _WCA_force(this, index1, index2, &force[0]);

    for (int dim=0; dim < 2; dim++) {
      if ( force[dim] != 0 ) {

        // update force arrays
        particles[index1].force()[dim] += force[dim];
        particles[index2].force()[dim] -= force[dim];
      }
    }
  }
}

void System::copyParticles(std::vector<Particle>& newParticles) {
  // Replace vector particles by newParticles.

  particles = newParticles;

  // UPDATING CELL LIST
  cellList.update<System>(this);
}

void System::copyParticles(System* system) {
  // Replace vector of particles by the one from system.

  particles = system->getParticles();

  // UPDATING CELL LIST
  cellList.update<System>(this);
}

void System::saveInitialState() {
  // Saves initial state of particles to output file.

  // output
  if ( dumpParticles ) {

    for (int i=0; i < getNumberParticles(); i++) { // output all particles
      // POSITIONS
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write<double>(particles[i].position()[dim]);
      }
      // ORIENTATIONS
      output.write<double>(particles[i].orientation()[0]); // output orientation
      // VELOCITIES
      velocitiesDumps[i] = output.tellp(); // location to dump velocities at next time step
      for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
        output.write<double>(0.0); // zero by default for initial frame
      }
    }
  }

  // reset dump
  resetDump();
}

void System::saveNewState(std::vector<Particle>& newParticles) {
  // Saves new state of particles to output file then copy it.

  // DUMP FRAME
  dumpFrame++;

  ////////////
  // SAVING //
  ////////////

  for (int i=0; i < getNumberParticles(); i++) { // output all particles

    // ACTIVE WORK and ORDER PARAMETER (computation)
    for (int dim=0; dim < 2; dim++) {
      // active work
      workSum[0] +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *(newParticles[i].position()[dim] - particles[i].position()[dim]) // NOTE: at this stage, newParticles[i].position() are not rewrapped, so this difference is the actual displacement
        /2;
      // force part of the active work
      workForceSum[0] +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *getTimeStep()*particles[i].force()[dim]/3/getPersistenceLength()
        /2;
      // orientation part of the active work
      workOrientationSum[0] +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *getTimeStep()*cos(particles[i].orientation()[0] - dim*M_PI/2)
        /2;
    }

    // WRAPPED COORDINATES
    for (int dim=0; dim < 2; dim++) {
      // keep particles in the box
      newParticles[i].position()[dim] =
        _wrapCoordinate<System>(this, newParticles[i].position()[dim]);
      // output wrapped position in each dimension
      if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
        output.write<double>(newParticles[i].position()[dim]);
      }
    }

    if ( dumpParticles && (dumpFrame - 1) % dumpPeriod == 0 ) {

      // VELOCITIES
      for (int dim=0; dim < 2; dim++) {
        output.write<double>(
          particles[i].velocity()[dim],
          velocitiesDumps[i] + dim*sizeof(double));
      }
    }

    if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {

      // ORIENTATION
      output.write<double>(newParticles[i].orientation()[0]);

      // VELOCITIES
      velocitiesDumps[i] = output.tellp(); // location to dump velocities at next time step
      for (int dim=0; dim < 2; dim++) {
        output.write<double>(0.0); // zero by default until rewrite at next time step
      }
    }


  }

  // ORDER PARAMETER
  orderSum[0] +=
    (getOrderParameterNorm(particles) + getOrderParameterNorm(newParticles))
      /2;

  #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
  // ORDER PARAMETER
  std::vector<double> orderOld = getOrderParameter(particles);
  double orderNormSqOld = pow(orderOld[0], 2) + pow(orderOld[1], 2);
  std::vector<double> orderNew = getOrderParameter(newParticles);
  double orderNormSqNew = pow(orderNew[0], 2) + pow(orderNew[1], 2);
  // GLOBAL PHASE
  double globalPhaseOld = getGlobalPhase(particles);
  double globalPhaseNew = getGlobalPhase(newParticles);
  // FIRST INTEGRAL
  torqueIntegral1[0] += (orderNormSqOld + orderNormSqNew)/2;
  // SECOND INTEGRAL
  for (int i=0; i < getNumberParticles(); i++) {
    #if 1 // method with global phase
    torqueIntegral2[0] += orderNormSqOld
      *pow(sin(particles[i].orientation()[0] - globalPhaseOld), 2)/2.0;
    torqueIntegral2[0] += orderNormSqNew
      *pow(sin(newParticles[i].orientation()[0] - globalPhaseNew), 2)/2.0;
    #else // method with double sum
    double sumSinOld (0.0), sumSinNew (0.0);
    for (int j=0; j < getNumberParticles(); j++) {
      sumSinOld += sin(
        particles[i].orientation()[0] - particles[j].orientation()[0]);
      sumSinNew += sin(
        newParticles[i].orientation()[0] - newParticles[j].orientation()[0]);
    }
    torqueIntegral2[0] += (pow(sumSinOld, 2) + pow(sumSinNew, 2))
      /2.0/pow(getNumberParticles(), 2);
    #endif
  }
  #endif

  // ACTIVE WORK and ORDER PARAMETER (output)
  if ( dumpFrame % (framesWork*dumpPeriod) == 0 ) {
    // compute normalised rates since last dump
    workSum[1] = workSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    workForceSum[1] = workForceSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    workOrientationSum[1] = workOrientationSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    orderSum[1] = orderSum[0]/(
      framesWork*dumpPeriod);
    // output normalised rates
    output.write<double>(workSum[1]);
    output.write<double>(workForceSum[1]);
    output.write<double>(workOrientationSum[1]);
    output.write<double>(orderSum[1]);
    // update time extensive quantities over trajectory since last reset
    workSum[2] += workSum[0]/getNumberParticles();
    workForceSum[2] += workForceSum[0]/getNumberParticles();
    workOrientationSum[2] += workOrientationSum[0]/getNumberParticles();
    orderSum[2] += orderSum[0];
    // reset sums
    workSum[0] = 0;
    workForceSum[0] = 0;
    workOrientationSum[0] = 0;
    orderSum[0] = 0;
    #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
    // compute normalised integrals since last dump
    torqueIntegral1[1] = torqueIntegral1[0]/(
      framesWork*dumpPeriod);
    torqueIntegral2[1] = torqueIntegral2[0]/(
      getNumberParticles()*framesWork*dumpPeriod);
    // update time (in units of time step) extensive integrals since last reset
    torqueIntegral1[2] += torqueIntegral1[0];
    torqueIntegral2[2] += torqueIntegral2[0]/
      getNumberParticles();
    // reset sums
    torqueIntegral1[0] = 0;
    torqueIntegral2[0] = 0;
    #endif
  }

  /////////////
  // COPYING //
  /////////////

  copyParticles(newParticles);
}


/***********
 * SYSTEM0 *
 ***********/

// CONSTRUCTORS

System0::System0() :
  param(new Parameters()),
  randomSeed(0), randomGenerator(),
  particles(0),
  cellList(),
  output(""), velocitiesDumps(),
  framesWork(0), dumpParticles(0), dumpPeriod(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {}

System0::System0(
  Parameters* parameters, int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(parameters),
  randomSeed(seed), randomGenerator(),
  particles(parameters->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(parameters->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    1/(parameters->getRotDiffusivity()*parameters->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // set diameters and system size
  std::vector<double> diameters (getNumberParticles(), 1.0);
  setDiameters(diameters);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPotentialParameter());
  output.write<double>(getPropulsionVelocity());
  output.write<double>(getTransDiffusivity());
  output.write<double>(getRotDiffusivity());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // write particles' diameters
  for (int i=0; i < getNumberParticles(); i++) {
    output.write<double>(getParticle(i)->diameter()[0]);
  }

  // put particles on a grid with random orientation
  int gridSize = ceil(sqrt(getNumberParticles())); // size of the grid on which to put the particles
  double gridSpacing = getSystemSize()/gridSize;
  for (int i=0; i < getNumberParticles(); i++) { // loop over particles
    // position on the grid
    particles[i].position()[0] = (i%gridSize)*gridSpacing;
    particles[i].position()[1] = (i/gridSize)*gridSpacing;
    // random orientation
    particles[i].orientation()[0] = 2*M_PI*randomGenerator.random01();
  }

  // initialise cell list
  double maxDiameter (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    maxDiameter = std::max(maxDiameter, getParticle(i)->diameter()[0]);
  }
  cellList.initialise<System0>(this, pow(2.*maxDiameter, 1./6.));
}

System0::System0(
  Parameters* parameters, std::vector<double>& diameters, int seed,
  std::string filename, int nWork, bool dump, int period) :
  param(parameters),
  randomSeed(seed), randomGenerator(),
  particles(parameters->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(parameters->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    1/(parameters->getRotDiffusivity()*parameters->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // set diameters and system size
  setDiameters(diameters);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPotentialParameter());
  output.write<double>(getPropulsionVelocity());
  output.write<double>(getTransDiffusivity());
  output.write<double>(getRotDiffusivity());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // write particles' diameters
  for (int i=0; i < parameters->getNumberParticles(); i++) {
    output.write<double>(getParticle(i)->diameter()[0]);
  }

  // put particles on a grid with random orientation
  int gridSize = ceil(sqrt(getNumberParticles())); // size of the grid on which to put the particles
  double gridSpacing = getSystemSize()/gridSize;
  for (int i=0; i < getNumberParticles(); i++) { // loop over particles
    // position on the grid
    particles[i].position()[0] = (i%gridSize)*gridSpacing;
    particles[i].position()[1] = (i/gridSize)*gridSpacing;
    // random orientation
    particles[i].orientation()[0] = 2*M_PI*randomGenerator.random01();
  }

  // initialise cell list
  double maxDiameter (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    maxDiameter = std::max(maxDiameter, getParticle(i)->diameter()[0]);
  }
  cellList.initialise<System0>(this, pow(2.0*maxDiameter, 1./6.));
}

System0::System0(
  System0* system, int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(system->getParameters()),
  randomSeed(seed), randomGenerator(),
  particles(system->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(system->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    1/(system->getRotDiffusivity()*system->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // get and re-set diameters
  std::vector<double> diameters (getNumberParticles(), 0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    diameters[i] = (system->getParticle(i))->diameter()[0];
  }
  setDiameters(diameters);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPotentialParameter());
  output.write<double>(getPropulsionVelocity());
  output.write<double>(getTransDiffusivity());
  output.write<double>(getRotDiffusivity());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // write particles' diameters
  for (int i=0; i < getNumberParticles(); i++) {
    output.write<double>(getParticle(i)->diameter()[0]);
  }

  // initialise cell list
  double maxDiameter (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    maxDiameter = std::max(maxDiameter, diameters[i]);
  }
  cellList.initialise<System0>(this, pow(2.0*maxDiameter, 1./6.));
  // copy particles and update cell list
  copyParticles(system);
  // copy dumps
  copyDump(system);
}

System0::System0(
  System0* system, std::vector<double>& diameters, int seed,
  std::string filename, int nWork, bool dump, int period) :
  param(system->getParameters()),
  randomSeed(seed), randomGenerator(),
  particles(system->getNumberParticles()),
  cellList(),
  output(filename), velocitiesDumps(system->getNumberParticles()),
  framesWork(nWork > 0 ? nWork : (int)
    1/(system->getRotDiffusivity()*system->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // get and re-set diameters
  setDiameters(diameters);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPotentialParameter());
  output.write<double>(getPropulsionVelocity());
  output.write<double>(getTransDiffusivity());
  output.write<double>(getRotDiffusivity());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // write particles' diameters
  for (int i=0; i < getNumberParticles(); i++) {
    output.write<double>(getParticle(i)->diameter()[0]);
  }

  // initialise cell list
  double maxDiameter (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    maxDiameter = std::max(maxDiameter, diameters[i]);
  }
  cellList.initialise<System0>(this, pow(2.0*maxDiameter, 1./6.));
  // copy particles and update cell list
  copyParticles(system, diameters);
    // NOTE: this sets once again all particle diameters, something more efficient could be done
  // copy dumps
  copyDump(system);
}

System0::System0(
  std::string inputFilename, int inputFrame, double dt,
  int seed, std::string filename,
  int nWork, bool dump, int period) :
  param(new Parameters()),
  randomSeed(seed), randomGenerator(),
  particles(0),
  cellList(),
  output(filename), velocitiesDumps(0),
  framesWork(nWork), dumpParticles(dump), dumpPeriod(period),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // load data
  Dat0 inputDat(inputFilename, false); // data object
  param = Parameters( // set parameters
    inputDat.getNumberParticles(),
    inputDat.getPotentialParameter(),
    inputDat.getPropulsionVelocity(),
    inputDat.getTransDiffusivity(),
    inputDat.getRotDiffusivity(),
    inputDat.getPackingFraction(),
    inputDat.getSystemSize(),
    dt > 0 ? dt : inputDat.getTimeStep());
  particles.resize(getNumberParticles()); // resize vector of particles
  velocitiesDumps.resize(2*getNumberParticles()); // resize vector of locations of velocity dumps
  std::vector<double> diameters = inputDat.getDiameters(); // diameters of particles
  for (int i=0; i < getNumberParticles(); i++) {
    // set positions
    for (int dim=0; dim < 2; dim++) {
      particles[i].position()[dim] = inputDat.getPosition(inputFrame, i, dim);
    }
    // set orientation
    particles[i].orientation()[0] = inputDat.getOrientation(inputFrame, i);
    // set diameter
    particles[i].diameter()[0] = diameters[i];
  }

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // write header with system parameters to output file
  output.write<int>(getNumberParticles());
  output.write<double>(getPotentialParameter());
  output.write<double>(getPropulsionVelocity());
  output.write<double>(getTransDiffusivity());
  output.write<double>(getRotDiffusivity());
  output.write<double>(getPersistenceLength());
  output.write<double>(getPackingFraction());
  output.write<double>(getSystemSize());
  output.write<int>(randomSeed);
  output.write<double>(getTimeStep());
  output.write<int>(framesWork);
  output.write<bool>(dumpParticles);
  output.write<int>(dumpPeriod);

  // write particles' diameters
  for (int i=0; i < getNumberParticles(); i++) {
    output.write<double>(getParticle(i)->diameter()[0]);
  }

  // initialise cell list
  double maxDiameter (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    maxDiameter = std::max(maxDiameter, getParticle(i)->diameter()[0]);
  }
  cellList.initialise<System0>(this, pow(2.0*maxDiameter, 1./6.));
}

// DESTRUCTORS

System0::~System0() {}

// METHODS

Parameters* System0::getParameters() { return &param; }

int System0::getNumberParticles() {
  return param.getNumberParticles(); }
double System0::getPotentialParameter() {
  return param.getPotentialParameter(); }
double System0::getPropulsionVelocity() {
  return param.getPropulsionVelocity(); }
double System0::getTransDiffusivity() {
  return param.getTransDiffusivity(); }
double System0::getRotDiffusivity() {
  return param.getRotDiffusivity(); }
double System0::getPersistenceLength() {
  return param.getPersistenceLength(); }
double System0::getPackingFraction() {
  return param.getPackingFraction(); }
double System0::getSystemSize() {
  return param.getSystemSize(); }
double System0::getTimeStep() {
  return param.getTimeStep(); }

void System0::setTimeStep(double dt) {
  param = Parameters(
    getNumberParticles(), getPotentialParameter(), getPropulsionVelocity(),
    getTransDiffusivity(), getRotDiffusivity(), getPackingFraction(),
    getSystemSize(), dt);
}

int System0::getRandomSeed() const { return randomSeed; }
rnd* System0::getRandomGenerator() { return &randomGenerator; }

Particle* System0::getParticle(int index) { return &(particles[index]); }
std::vector<Particle> System0::getParticles() { return particles; }

CellList* System0::getCellList() { return &cellList; }

std::string System0::getOutputFile() const { return output.getOutputFile(); }

int System0::getDump() { return dumpFrame; }

void System0::resetDump() {
  // Reset time-extensive quantities over trajectory.

  dumpFrame = 0;

  workSum[0] = 0;
  workForceSum[0] = 0;
  workOrientationSum[0] = 0;
  orderSum[0] = 0;

  workSum[2] = 0;
  workForceSum[2] = 0;
  workOrientationSum[2] = 0;
  orderSum[2] = 0;
}

void System0::copyDump(System0* system) {
  // Copy dumps from other system.
  // WARNING: This also copies the index of last frame dumped. Consistency
  //          has to be checked.

  dumpFrame = system->getDump();

  workSum[2] = system->getTotalWork();
  workForceSum[2] = system->getTotalWorkForce();
  workOrientationSum[2] = system->getTotalWorkOrientation();
  orderSum[2] = system->getTotalOrder();
}

double System0::getWork() { return workSum[1]; }
double System0::getWorkForce() { return workForceSum[1]; }
double System0::getWorkOrientation() { return workOrientationSum[1]; }
double System0::getOrder() { return orderSum[1]; }

double System0::getTotalWork() { return workSum[2]; }
double System0::getTotalWorkForce() { return workForceSum[2]; }
double System0::getTotalWorkOrientation() { return workOrientationSum[2]; }
double System0::getTotalOrder() { return orderSum[2]; }

double System0::diffPeriodic(double const& x1, double const& x2) {
  // Returns algebraic distance from `x1' to `x2' on a line taking into account
  // periodic boundary condition of the system.

  return _diffPeriodic<System0>(this, x1, x2);
}

double System0::getDistance(int const& index1, int const& index2) {
  // Returns distance between two particles in a given system.

  return _getDistance<System0>(this, index1, index2);
}

void System0::WCA_force(int const& index1, int const& index2) {
  // Compute WCA forces between particles[index1] and particles[index2],
  // and add to particles[index1].force() and particles[index2].force().

  if ( index1 != index2 ) { // only consider different particles

    double force[2];
    _WCA_force(this, index1, index2, &force[0]);

    for (int dim=0; dim < 2; dim++) {
      if ( force[dim] != 0 ) {

        // update force arrays
        particles[index1].force()[dim] += force[dim];
        particles[index2].force()[dim] -= force[dim];
      }
    }
  }
}

void System0::copyParticles(std::vector<Particle>& newParticles) {
  // Replace vector particles by newParticles.

  particles = newParticles;

  // UPDATING CELL LIST
  cellList.update<System0>(this);
}

void System0::copyParticles(System0* system) {
  // Replace vector of particles by the one from system.

  // COPYING PARTICLES
  particles = system->getParticles();

  // UPDATING CELL LIST
  cellList.update<System0>(this);
}

void System0::copyParticles(System0* system, std::vector<double>& diameters) {
  // Replace vector of particles by the one from system.
  // Replace particles' diameters by the ones in diameters.

  // COPYING PARTICLES
  particles = system->getParticles();

  // SETTING DIAMETERS
  if ( (int) diameters.size() != getNumberParticles() ) { // check size of diameters array
    std::cout << "Size of diameters array does not match number of particles."; std::cout << std::endl;
    exit(0);
  }
  for (int i=0; i < getNumberParticles(); i++) {
    particles[i].diameter()[0] = diameters[i];
  }

  // UPDATING CELL LIST
  cellList.update<System0>(this);
}

void System0::setDiameters(std::vector<double>& diameters) {
  // Set all diameters and re-define system size.

  if ( (int) diameters.size() != getNumberParticles() ) { // check size of diameters array
    std::cout << "Size of diameters array does not match number of particles.";
    std::cout << std::endl;
    exit(0);
  }
  double totalArea (0.0);
  for (int i=0; i < getNumberParticles(); i++) {
    getParticle(i)->diameter()[0] = diameters[i]; // set diameter
    totalArea += M_PI*pow(diameters[i], 2)/4.0; // add corresponding area
  }
  param = Parameters(
    getNumberParticles(), getPotentialParameter(), getPropulsionVelocity(),
    getTransDiffusivity(), getRotDiffusivity(), getPackingFraction(),
    sqrt(getPackingFraction()/totalArea), getTimeStep());
}

void System0::saveInitialState() {
  // Saves initial state of particles to output file.

  // output
  if ( dumpParticles ) {

    for (int i=0; i < getNumberParticles(); i++) { // output all particles
      // POSITIONS
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write<double>(particles[i].position()[dim]);
      }
      // ORIENTATIONS
      output.write<double>(particles[i].orientation()[0]); // output orientation
      // VELOCITIES
      velocitiesDumps[i] = output.tellp(); // location to dump velocities at next time step
      for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
        output.write<double>(0.0); // zero by default for initial frame
      }
    }
  }

  // reset dump
  resetDump();
}

void System0::saveNewState(std::vector<Particle>& newParticles) {
  // Saves new state of particles to output file then copy it.

  // DUMP FRAME
  dumpFrame++;

  ////////////
  // SAVING //
  ////////////

  for (int i=0; i < getNumberParticles(); i++) { // output all particles

    // ACTIVE WORK and ORDER PARAMETER (computation)
    for (int dim=0; dim < 2; dim++) {
      // active work
      workSum[0] +=
        getPropulsionVelocity()*
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *(newParticles[i].position()[dim] - particles[i].position()[dim]) // NOTE: at this stage, newParticles[i].position() are not rewrapped, so this difference is the actual displacement
        /2;
      // force part of the active work
      workForceSum[0] +=
        getPropulsionVelocity()*
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *getTimeStep()*getPotentialParameter()*particles[i].force()[dim]
        /2;
      // orientation part of the active work
      workOrientationSum[0] +=
        pow(getPropulsionVelocity(), 2)*
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *getTimeStep()*cos(particles[i].orientation()[0] - dim*M_PI/2)
        /2;
    }

    // WRAPPED COORDINATES
    for (int dim=0; dim < 2; dim++) {
      // keep particles in the box
      newParticles[i].position()[dim] =
        _wrapCoordinate<System0>(this, newParticles[i].position()[dim]);
      // output wrapped position in each dimension
      if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
        output.write<double>(newParticles[i].position()[dim]);
      }
    }

    if ( dumpParticles && (dumpFrame - 1) % dumpPeriod == 0 ) {

      // VELOCITIES
      for (int dim=0; dim < 2; dim++) {
        output.write<double>(
          particles[i].velocity()[dim],
          velocitiesDumps[i] + dim*sizeof(double));
      }
    }

    if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {

      // ORIENTATION
      output.write<double>(newParticles[i].orientation()[0]);

      // VELOCITIES
      velocitiesDumps[i] = output.tellp(); // location to dump velocities at next time step
      for (int dim=0; dim < 2; dim++) {
        output.write<double>(0.0); // zero by default until rewrite at next time step
      }
    }
  }

  // ORDER PARAMETER
  orderSum[0] +=
    (getOrderParameterNorm(particles) + getOrderParameterNorm(newParticles))
      /2;

  // ACTIVE WORK and ORDER PARAMETER (output)
  if ( dumpFrame % (framesWork*dumpPeriod) == 0 ) {
    // compute normalised rates since last dump
    workSum[1] = workSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    workForceSum[1] = workForceSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    workOrientationSum[1] = workOrientationSum[0]/(
      getNumberParticles()*getTimeStep()*framesWork*dumpPeriod);
    orderSum[1] = orderSum[0]/(
      framesWork*dumpPeriod);
    // output normalised rates
    output.write<double>(workSum[1]);
    output.write<double>(workForceSum[1]);
    output.write<double>(workOrientationSum[1]);
    output.write<double>(orderSum[1]);
    // update time extensive quantities over trajectory since last reset
    workSum[2] += workSum[0]/getNumberParticles();
    workForceSum[2] += workForceSum[0]/getNumberParticles();
    workOrientationSum[2] += workOrientationSum[0]/getNumberParticles();
    orderSum[2] += orderSum[0];
    // reset sums
    workSum[0] = 0;
    workForceSum[0] = 0;
    workOrientationSum[0] = 0;
    orderSum[0] = 0;
  }

  /////////////
  // COPYING //
  /////////////

  copyParticles(newParticles);
}


/**********
 * ROTORS *
 **********/

// CONSTRUCTORS

Rotors::Rotors(
  int N, double Dr, double g, double dt, int seed,
  std::string filename, int period) :
  numberParticles(N), rotDiffusivity(Dr), torqueParameter(g),
  timeStep(dt), dumpPeriod(period),
  randomSeed(seed), randomGenerator(),
  orientations(N), torques(N, 0.0),
  output(filename), dumpFrame(-1) {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // give random orientations to rotors
  for (int i=0; i < numberParticles; i++) { // loop over rotors
    orientations[i] = 2*M_PI*randomGenerator.random01();
  }

  // write header with system parameters to output file
  output.write<int>(numberParticles);
  output.write<double>(rotDiffusivity);
  output.write<double>(torqueParameter);
  output.write<double>(timeStep);
  output.write<int>(dumpPeriod);
  output.write<int>(seed);
}

// DESTRUCTORS

Rotors::~Rotors() {}

// METHODS

int Rotors::getNumberParticles() const { return numberParticles; }
double Rotors::getRotDiffusivity() const { return rotDiffusivity; }
double Rotors::getTorqueParameter() const { return torqueParameter; }

double Rotors::getTimeStep() const { return timeStep; }

double* Rotors::getOrientation(int index) { return &(orientations[index]); }
double* Rotors::getTorque(int index) { return &(torques[index]); }

rnd* Rotors::getRandomGenerator() { return &randomGenerator; }

void Rotors::saveState() {
  // Save current state.

  dumpFrame++;

  if ( dumpFrame % dumpPeriod == 0 ) {
    for (int i=0; i < numberParticles; i++) {
      output.write<double>(orientations[i]);
    }
  }
}


///////////////
// FUNCTIONS //
///////////////

double getGlobalPhase(std::vector<Particle>& particles) {
  // Returns global phase.

  std::vector<double> order = getOrderParameter(particles);

  return getAngle(order[0]/sqrt(pow(order[0], 2) + pow(order[1], 2)), order[1]);
}

std::vector<double> getOrderParameter(std::vector<Particle>& particles) {
  // Returns order parameter.

  std::vector<double> order(2, 0);

  for (int i=0; i < (int) particles.size(); i++) { // loop over particles
    for (int dim=0; dim < 2; dim++) { // loop over dimensions
      order[dim] += cos(particles[i].orientation()[0] - dim*M_PI/2)
        /particles.size();
    }
  }

  return order;
}

double getOrderParameterNorm(std::vector<Particle>& particles) {
  // Returns order parameter norm.

  std::vector<double> order = getOrderParameter(particles);

  return sqrt(pow(order[0], 2) + pow(order[1], 2));
}

void _WCA_force(
  System* system, int const& index1, int const& index2, double* force) {
  // Writes to `force' the force deriving from the WCA potential between
  // particles `index1' and `index2'.

  force[0] = 0.0;
  force[1] = 0.0;

  double dist = system->getDistance(index1, index2); // dimensionless distance between particles

  if (dist < pow(2., 1./6.)) { // distance lower than cut-off

    // compute force
    double coeff = 48.0/pow(dist, 14.0) - 24.0/pow(dist, 8.0);
    for (int dim=0; dim < 2; dim++) {
      force[dim] = system->diffPeriodic(
          (system->getParticle(index2))->position()[dim],
          (system->getParticle(index1))->position()[dim])
        *coeff;
    }
  }
}

void _WCA_force(
  System0* system, int const& index1, int const& index2, double* force) {
  // Writes to `force' the force deriving from the WCA potential between
  // particles `index1' and `index2'.

  force[0] = 0.0;
  force[1] = 0.0;

  double dist = system->getDistance(index1, index2); // distance between particles
  double sigma =
    ((system->getParticle(index1))->diameter()[0]
    + (system->getParticle(index2))->diameter()[0])/2.0; // equivalent diameter

  if (dist/sigma < pow(2., 1./6.)) { // distance lower than cut-off

    // compute force
    double coeff =
      (48.0/pow(dist/sigma, 14.0) - 24.0/pow(dist/sigma, 8.0))/pow(sigma, 2.0);
    for (int dim=0; dim < 2; dim++) {
      force[dim] = system->diffPeriodic(
          (system->getParticle(index2))->position()[dim],
          (system->getParticle(index1))->position()[dim])
        *coeff;
    }
  }
}
