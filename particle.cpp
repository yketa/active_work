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


// PARTICLE

// CONSTRUCTORS

Particle::Particle() : r {0, 0}, theta (0) {}
Particle::Particle(double x, double y, double ang) :
  r {x, y}, theta (ang) {}

// METHODS

double* Particle::position() { return &r[0]; } // returns pointer to position
double* Particle::orientation() { return &theta; } // returns pointer to orientation

void Particle::copy(Particle *particle) {
  // Copy content of other particle.

  for (int dim=0; dim < 2; dim++) {
    r[dim] = particle->position()[dim];
  }
  theta = particle->orientation()[0];
}

// SYSTEM

// CONSTRUCTORS

System::System(
  int N, double lp, double phi, int seed, double dt, std::string filename,
  int nWork, bool dump, int period) :
  numberParticles(N), persistenceLength(lp), packingFraction(phi),
    systemSize(sqrt(M_PI*N/phi)/2), randomSeed(seed), timeStep(dt),
    outputFile(filename),
  framesWork(nWork > 0 ? nWork : (int) lp/(dt*period)), dumpParticles(dump),
    dumpPeriod(period),
  randomGenerator(), particles(N), outputFileStream(
    filename.c_str(), std::ios::out | std::ios::binary),
  dumpFrame(-1), workSum(0) {

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

double System::getDistance(int index1, int index2) {
  // Returns distance between two particles in a given system.

  double diff[2];
  double coord[2];
  for (int dim=0; dim < 2; dim++) {
    coord[0] = particles[index1].position()[dim];
    coord[1] = particles[index2].position()[dim];
    diff[dim] = std::min(
      std::min(
        fabs(coord[0]) + fabs(systemSize - coord[1]),
        fabs(systemSize - coord[0]) + fabs(coord[1])),
      fabs(coord[0] - coord[1])); // separation in each direction accounting for periodic boundaries
  }
  return sqrt(pow(diff[0], 2) + pow(diff[1], 2));
}

void System::WCA_potential(int index1, int index2, double *force) {
  // Writes WCA force acting on particles[index1] by particles[index2] onto
  // `force`.

  double dist = this->getDistance(index1, index2); // dimensionless distance between particles

  if (dist < pow(2, 1./6.)) { // distance lower than cut-off
    double coeff = 48/pow(dist, 14) - 24/pow(dist, 8);
    for (int dim=0; dim < 2; dim++) {
      force[dim] =
        remainder(
          particles[index1].position()[dim] - particles[index2].position()[dim],
          systemSize/2)
        *coeff;
    }
  }
  else { // distance greater than cut-off
    force[0] = 0;
    force[1] = 0;
  }
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

void System::saveNewState(Particle *newParticles) {
  // Saves new state of particles to output file then copy it.

  // DUMP FRAME
  dumpFrame++;

  // SAVING
  for (int i=0; i < numberParticles; i++) { // output all particles

    // ACTIVE WORK (computation)
    for (int dim=0; dim < 2; dim++) {
      workSum +=
        (cos(newParticles[i].orientation()[0] - dim*M_PI/2)
          + cos(particles[i].orientation()[0] - dim*M_PI/2))
        *(newParticles[i].position()[dim] - particles[i].position()[dim]) // NOTE: at this stage, newPartices[i].position() are not rewrapped, so this difference is the actual displacement
        /2;
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

  // ACTIVE WORK (output)
  if ( dumpFrame % (framesWork*dumpPeriod) == 0 ) {
    workSum /= numberParticles*timeStep*framesWork*dumpPeriod;
    outputFileStream.write((char*) &workSum, sizeof(double));
    workSum = 0;
  }

  // COPYING
  for (int i=0; i < numberParticles; i++) {
    particles[i].copy(&(newParticles[i]));
  }
}
