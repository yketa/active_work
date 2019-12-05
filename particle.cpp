#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

#include "particle.hpp"
#include "maths.hpp"

/////////////
// CLASSES //
/////////////

/************
 * PARTICLE *
 ************/

// CONSTRUCTORS

Particle::Particle() : r {0, 0}, theta (0), f {0, 0} {}
Particle::Particle(double x, double y, double ang) :
  r {x, y}, theta (ang), f {0, 0} {}

// METHODS

double* Particle::position() { return &r[0]; } // returns pointer to position
double* Particle::orientation() { return &theta; } // returns pointer to orientation

double* Particle::force() {return &f[0]; }; // returns pointer to force


/*************
 * CELL LIST *
 *************/

// CONSTRUCTORS

CellList::CellList() {}

// DESTRUCTORS

CellList::~CellList() {}

// METHODS

int CellList::getNumberBoxes() { return numberBoxes; } // return number of boxes in each dimension
std::vector<int> CellList::getCell(int const &index) { return cellList[index]; } // return vector of indexes in cell

void CellList::initialise(System *system, double const& rcut) {
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

void CellList::update(System *system) {
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

int CellList::index(Particle *particle) {
  // Index of the box corresponding to a given particle.

  return (int) ((particle->position())[0]/sizeBox)
    + numberBoxes*((int) ((particle->position())[1]/sizeBox));
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


/**********
 * SYSTEM *
 **********/

// CONSTRUCTORS

System::System() :
  param(new Parameters()),
  randomSeed(0), randomGenerator(),
  particles(0),
  cellList(),
  output(""),
  framesWork(0), dumpParticles(0), dumpPeriod(0),
  biasingParameter(0),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {}

System::System(
  Parameters* parameters, int seed, std::string filename,
  int nWork, bool dump, int period, double sValue) :
  param(parameters),
  randomSeed(seed), randomGenerator(),
  particles(parameters->getNumberParticles()),
  cellList(),
  output(filename),
  framesWork(nWork > 0 ? nWork : (int)
    parameters->getPersistenceLength()/(parameters->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  biasingParameter(sValue),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

    // set seed of random generator
    randomGenerator.setSeed(randomSeed);

    // write header with system parameters to output file
    output.write(getNumberParticles());
    output.write(getPersistenceLength());
    output.write(getPackingFraction());
    output.write(getSystemSize());
    output.write(randomSeed);
    output.write(getTimeStep());
    output.write(framesWork);
    output.write(dumpParticles);
    output.write(dumpPeriod);

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
    cellList.initialise(this, pow(2., 1./6.));
}

System::System(
  System* system, int seed, std::string filename,
  int nWork, bool dump, int period, double sValue) :
  param(system->getParameters()),
  randomSeed(seed), randomGenerator(),
  particles(system->getNumberParticles()),
  cellList(),
  output(filename),
  framesWork(nWork > 0 ? nWork : (int)
    system->getPersistenceLength()/(system->getTimeStep()*period)),
    dumpParticles(dump), dumpPeriod(period),
  biasingParameter(sValue),
  dumpFrame(-1),
  workSum {0, 0, 0}, workForceSum {0, 0, 0}, workOrientationSum {0, 0, 0},
    orderSum {0, 0, 0} {

  // set seed of random generator
  randomGenerator.setSeed(randomSeed);

  // write header with system parameters to output file
  output.write(getNumberParticles());
  output.write(getPersistenceLength());
  output.write(getPackingFraction());
  output.write(getSystemSize());
  output.write(randomSeed);
  output.write(getTimeStep());
  output.write(framesWork);
  output.write(dumpParticles);
  output.write(dumpPeriod);

  // initialise cell list
  cellList.initialise(this, pow(2., 1./6.));
  // copy particles and update cell list
  copyParticles(system);
  // copy dumps
  copyDump(system);
}

// DESTRUCTORS

System::~System() {}

// METHODS

Parameters* System::getParameters() { return &param; }

int System::getNumberParticles() const {
  return param.getNumberParticles(); }
double System::getPersistenceLength() const {
  return param.getPersistenceLength(); }
double System::getPackingFraction() const {
  return param.getPackingFraction(); }
double System::getSystemSize() const {
  return param.getSystemSize(); }
double System::getTimeStep() const {
  return param.getTimeStep(); }

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
}

double System::getWork() { return workSum[1]; }
double System::getWorkForce() { return workForceSum[1]; }
double System::getWorkOrientation() { return workOrientationSum[1]; }
double System::getOrder() { return orderSum[1]; }

double System::getTotalWork() { return workSum[2]; }
double System::getTotalWorkForce() { return workForceSum[2]; }
double System::getTotalWorkOrientation() { return workOrientationSum[2]; }
double System::getTotalOrder() { return orderSum[2]; }

double System::diffPeriodic(double const& x1, double const& x2) {
  // Returns algebraic distance from `x1' to `x2' on a line taking into account
  // periodic boundary condition of the system.

  double diff = x2 - x1;

  if ( fabs(diff) > getSystemSize()/2 ) {
    double diff1 = fabs(x1) + fabs(getSystemSize() - x2);
    double diff2 = fabs(getSystemSize() - x1) + fabs(x2);
    if ( diff1 < diff2 ) { diff = diff1; }
    else { diff = diff2; }
    diff *= (x2 > x1 ? -1 : 1);
  }

  return diff;
}

double System::getDistance(int const& index1, int const& index2) {
  // Returns distance between two particles in a given system.

  return sqrt(
      pow(
        diffPeriodic( // separation in x position
          particles[index1].position()[0],
          particles[index2].position()[0]),
        2)
      + pow( // separation in y position
        diffPeriodic(
          particles[index1].position()[1],
          particles[index2].position()[1]),
        2));
}

void System::WCA_potential(int const& index1, int const& index2,
  std::vector<Particle>& newParticles) {
  // Compute WCA forces between particles[index1] and particles[index2],
  // add to particles[index1].force() and particles[index2].force(), and
  // increments positions in particles[index1].position() and
  // particles[index2].position().

  if ( index1 != index2 ) { // only consider different particles

    double dist = getDistance(index1, index2); // dimensionless distance between particles

    if (dist < pow(2, 1./6.)) { // distance lower than cut-off

      double coeff = 48/pow(dist, 14) - 24/pow(dist, 8);
      double force;

      for (int dim=0; dim < 2; dim++) {

        // compute force
        force = diffPeriodic(
            particles[index2].position()[dim],
            particles[index1].position()[dim])
          *coeff;

        // update force arrays
        particles[index1].force()[dim] += force;
        particles[index2].force()[dim] -= force;

        // increment positions
        newParticles[index1].position()[dim] +=
          getTimeStep()*force/3/getPersistenceLength();
        newParticles[index2].position()[dim] -=
          getTimeStep()*force/3/getPersistenceLength();
      }
    }
  }
}

void System::copyParticles(std::vector<Particle>& newParticles) {
  // Replace vector particles by newParticles.

  particles = newParticles;

  // UPDATING CELL LIST
  cellList.update(this);
}

void System::copyParticles(System* system) {
  // Replace vector of particles by the one from system.

  particles = system->getParticles();

  // UPDATING CELL LIST
  cellList.update(this);
}

void System::saveInitialState() {
  // Saves initial state of particles to output file.

  // output
  if ( dumpParticles ) {

    for (int i=0; i < getNumberParticles(); i++) { // output all particles
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write(particles[i].position()[dim]);
      }
      output.write(particles[i].orientation()[0]); // output orientation
    }
  }

  // reset dump
  resetDump();
}

void System::saveNewState(std::vector<Particle>& newParticles) {
  // Saves new state of particles to output file then copy it.

  // DUMP FRAME
  dumpFrame++;

  // SAVING
  double orderParameterIn[2] {0, 0}, orderParameterFin[2] {0, 0};
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
      // order parameter
      orderParameterIn[dim] +=
        cos(particles[i].orientation()[0] - dim*M_PI/2);
      orderParameterFin[dim] +=
        cos(newParticles[i].orientation()[0] - dim*M_PI/2);
    }

    // WRAPPED COORDINATES
    for (int dim=0; dim < 2; dim++) {
      // keep particles in the box
      newParticles[i].position()[dim] =
        std::remainder(newParticles[i].position()[dim], getSystemSize());
      if (newParticles[i].position()[dim] < 0) {
        newParticles[i].position()[dim] += getSystemSize();
      }
      // output wrapped position in each dimension
      if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
        output.write(newParticles[i].position()[dim]);
      }
    }

    // ORIENTATION
    if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
      output.write(newParticles[i].orientation()[0]);
    }
  }
  orderSum[0] +=
    (sqrt(pow(orderParameterIn[0], 2) + pow(orderParameterIn[1], 2))
      + sqrt(pow(orderParameterFin[0], 2) + pow(orderParameterFin[1], 2)))
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
      getNumberParticles()*framesWork*dumpPeriod);
    // output normalised rates
    output.write(workSum[1]);
    output.write(workForceSum[1]);
    output.write(workOrientationSum[1]);
    output.write(orderSum[1]);
    // update time extensive quantities over trajectory since last reset
    workSum[2] += workSum[0]/getNumberParticles();
    workForceSum[2] += workForceSum[0]/getNumberParticles();
    workOrientationSum[2] += workOrientationSum[0]/getNumberParticles();
    orderSum[2] += orderSum[0]/getNumberParticles();
    // reset sums
    workSum[0] = 0;
    workForceSum[0] = 0;
    workOrientationSum[0] = 0;
    orderSum[0] = 0;
  }

  // COPYING
  copyParticles(newParticles);
}
