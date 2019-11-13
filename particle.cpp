#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

#include "particle.h"
#include "maths.h"

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
  for (int i=0; i < cellList.size(); i++) {
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

System::System(
  int N, double lp, double phi, int seed, double dt, std::string filename,
  int nWork, bool dump, int period) :
  numberParticles(N), persistenceLength(lp), packingFraction(phi),
    systemSize(sqrt(M_PI*N/phi)/2), randomSeed(seed), timeStep(dt),
    outputFile(filename),
  framesWork(nWork > 0 ? nWork : (int) lp/(dt*period)), dumpParticles(dump),
    dumpPeriod(period),
  randomGenerator(), particles(N),
    outputFileStream(filename.c_str(), std::ios::out | std::ios::binary),
    cellList(),
  dumpFrame(-1),
    workSum(0), workForceSum(0), workOrientationSum(0), orderSum(0) {

    // set seed of random generator
    randomGenerator.setSeed(seed);

    // writing header with system parameters to output file
    outputFileStream.write((char*) &numberParticles, sizeof(int));
    outputFileStream.write((char*) &persistenceLength, sizeof(double));
    outputFileStream.write((char*) &packingFraction, sizeof(double));
    outputFileStream.write((char*) &systemSize, sizeof(double));
    outputFileStream.write((char*) &randomSeed, sizeof(int));
    outputFileStream.write((char*) &timeStep, sizeof(double));
    outputFileStream.write((char*) &framesWork, sizeof(int));
    outputFileStream.write((char*) &dumpParticles, sizeof(bool));
    outputFileStream.write((char*) &dumpPeriod, sizeof(int));

    // putting particles on a grid with random orientation
    int gridSize = ceil(sqrt(numberParticles)); // size of the grid on which to put the particles
    double gridSpacing = systemSize/gridSize;
    for (int i=0; i < numberParticles; i++) { // loop over particles
      // position on the grid
      particles[i].position()[0] = (i%gridSize)*gridSpacing;
      particles[i].position()[1] = (i/gridSize)*gridSpacing;
      // random orientation
      particles[i].orientation()[0] = 2*M_PI*randomGenerator.random01();
    }

    // initialise cell list
    cellList.initialise(this, pow(2., 1./6.));
}

// DESTRUCTORS

System::~System() { outputFileStream.close(); }

// METHODS

int System::getNumberParticles() const { return numberParticles; }
double System::getPersistenceLength() const { return persistenceLength; }
double System::getPackingFraction() const { return packingFraction; }
double System::getSystemSize() const { return systemSize; }
int System::getRandomSeed() const { return randomSeed; }
double System::getTimeStep() const { return timeStep; }
std::string System::getOutputFile() const { return outputFile; }

rnd* System::getRandomGenerator() { return &randomGenerator; }
Particle* System::getParticle(int index) { return &(particles[index]); }
std::ofstream* System::getOutputFileStream() { return &outputFileStream; }
CellList* System::getCellList() { return &cellList; } // returns pointer to CellList object

double System::diffPeriodic(double const& x1, double const& x2) {
  // Returns algebraic distance from `x1' to `x2' on a line taking into account
  // periodic boundary condition of the system.

  double diff = x2 - x1;

  if ( fabs(diff) > systemSize/2 ) {
    double diff1 = fabs(x1) + fabs(systemSize - x2);
    double diff2 = fabs(systemSize - x1) + fabs(x2);
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

    double dist = this->getDistance(index1, index2); // dimensionless distance between particles

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
          timeStep*force/3/persistenceLength;
        newParticles[index2].position()[dim] -=
          timeStep*force/3/persistenceLength;
      }
    }
  }
}

void System::copyParticles(std::vector<Particle>& newParticles) {
  // Replace vector particles by newParticles.

  particles = newParticles;
}

void System::saveInitialState() {
  // Saves initial state of particles to output file.

  // output
  if ( dumpParticles ) {

    for (int i=0; i < numberParticles; i++) { // output all particles
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        outputFileStream.write((char*) &(particles[i].position()[dim]),
          sizeof(double));
      }
      outputFileStream.write((char*) particles[i].orientation(),
        sizeof(double)); // output orientation
    }
  }

  // last frame dumped
  dumpFrame = 0;
}

void System::saveNewState(std::vector<Particle>& newParticles) {
  // Saves new state of particles to output file then copy it.

  // DUMP FRAME
  dumpFrame++;

  // SAVING
  double orderParameterIn[2] {0, 0}, orderParameterFin[2] {0, 0};
  for (int i=0; i < numberParticles; i++) { // output all particles

    // ACTIVE WORK and ORDER PARAMETER (computation)
    for (int dim=0; dim < 2; dim++) {
      // active work
      workSum +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *(newParticles[i].position()[dim] - particles[i].position()[dim]) // NOTE: at this stage, newParticles[i].position() are not rewrapped, so this difference is the actual displacement
        /2;
      // force part of the active work
      workForceSum +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *timeStep*particles[i].force()[dim]/3/persistenceLength
        /2;
      // orientation part of the active work
      workOrientationSum +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *timeStep*cos(particles[i].orientation()[0] - dim*M_PI/2)
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
        std::remainder(newParticles[i].position()[dim], systemSize);
      if (newParticles[i].position()[dim] < 0) {
        newParticles[i].position()[dim] += systemSize;
      }
      // output wrapped position in each dimension
      if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
        outputFileStream.write((char*) &(newParticles[i].position()[dim]),
          sizeof(double));
      }
    }

    // ORIENTATION
    if ( dumpParticles && dumpFrame % dumpPeriod == 0 ) {
      outputFileStream.write((char*) newParticles[i].orientation(),
        sizeof(double));
    }
  }
  orderSum +=
    (sqrt(pow(orderParameterIn[0], 2) + pow(orderParameterIn[1], 2))
      + sqrt(pow(orderParameterFin[0], 2) + pow(orderParameterFin[1], 2)))
    /2;

  // ACTIVE WORK and ORDER PARAMETER (output)
  if ( dumpFrame % (framesWork*dumpPeriod) == 0 ) {
    workSum /= numberParticles*timeStep*framesWork*dumpPeriod;
    workForceSum /= numberParticles*timeStep*framesWork*dumpPeriod;
    workOrientationSum /= numberParticles*timeStep*framesWork*dumpPeriod;
    orderSum /= numberParticles*framesWork*dumpPeriod;
    outputFileStream.write((char*) &workSum, sizeof(double));
    outputFileStream.write((char*) &workForceSum, sizeof(double));
    outputFileStream.write((char*) &workOrientationSum, sizeof(double));
    outputFileStream.write((char*) &orderSum, sizeof(double));
    workSum = 0;
    workForceSum = 0;
    workOrientationSum = 0;
    orderSum = 0;
  }

  // COPYING
  copyParticles(newParticles);
  // UPDATING CELL LIST
  cellList.update(this);
}
