#include <cmath>
#include <vector>
#include <string>

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

System::System(int N, double lp, double L, int seed, std::string filename) :
  numberParticles(N), persistenceLength(lp), systemSize(L), randomSeed(seed),
    outputFile(filename),
  particles(N), randomGenerator(), outputFileStream(
    filename.c_str(), std::ios::out | std::ios::binary),
  increments(N, std::vector<double>(2, 0)) {

    randomGenerator.setSeed(seed); // set random seed

    // writing header with system parameters to output file
    outputFileStream.write((char*) &N, sizeof(int));
    outputFileStream.write((char*) &lp, sizeof(double));
    outputFileStream.write((char*) &L, sizeof(double));
    outputFileStream.write((char*) &seed, sizeof(int));
}

// METHODS

int System::getNumberParticles() const { return numberParticles; }
double System::getPersistenceLength() const { return persistenceLength; }
double System::getSystemSize() const { return systemSize; }
int System::getRandomSeed() const { return randomSeed; }
std::string System::getOutputFile() const { return outputFile; }

rnd* System::getRandomGenerator() { return &randomGenerator; }
Particle* System::getParticle(int index) { return &(particles[index]); }
std::ofstream* System::getOutputFileStream() { return &outputFileStream; }

double System::getDistance(int index1, int index2) {
  // Returns distance between two particles in a given system.

  double diff[2];
  for (int i=0; i < 2; i++) {
    diff[i] = fmod(
      fabs(particles[index1].position()[i] - particles[index2].position()[i]),
      systemSize/2
    ); // separation in each direction accounting for periodic boundaries
  }
  return sqrt(pow(diff[0], 2) + pow(diff[1], 2));
}

void System::WCA_potential(int index1, int index2, double *force) {
  // Writes WCA force acting on particles[index1] by particles[index2] onto
  // `force`.

  double dist = this->getDistance(index1, index2); // dimensionless distance between particles

  if (dist < pow(2, 1./6.)) { // distance lower than cut-off radius
    double coeff = 48/pow(dist, 14) - 24/pow(dist, 8);
    for (int i=0; i < 2; i++) {
      force[i] =
        (particles[index1].position()[i] - particles[index2].position()[i])
          *coeff;
    }
  }
}

void System::saveInitialState() {
  // Saves initial state of particles to output file.

  double dt = 0;
  outputFileStream.write((char*) &dt, sizeof(double));

  for (int i=0; i < numberParticles; i++) { // output all particles
    for (int j=0; j < 2; j++) { // output twice from wrapped and unwrapped coordinates
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        outputFileStream.write((char*) &(particles[i].position()[dim]),
          sizeof(double));
      }
    }
    outputFileStream.write((char*) particles[i].orientation(),
      sizeof(double)); // output orientation
  }
}

void System::saveNewState(Particle *newParticles, double const& timeStep) {
  // Saves new state of particles to output file then copy it.

  // SAVING

  outputFileStream.write((char*) &timeStep, sizeof(double));

  double newPosition;

  for (int i=0; i < numberParticles; i++) { // output all particles

    // WRAPPED COORDINATES
    for (int dim=0; dim < 2; dim++) {
      // output wrapped position in each dimension
      outputFileStream.write((char*) &(newParticles[i].position()[dim]),
        sizeof(double));
    }

    // UNWRAPPED COORDINATES
    for (int dim=0; dim < 2; dim++) {
      // updating increments
      if (
        (particles[i].position()[dim] - systemSize/2)*
          (newParticles[i].position()[dim] - systemSize/2) < 0 // position switches "sign"
        && fabs(particles[i].position()[dim] - systemSize/2) > systemSize/4 // particle not in the centre of the box
      ) {
        increments[i][dim] +=
          ((particles[i].position()[dim] > systemSize/2)
            - (particles[i].position()[dim] < systemSize/2))* // sign of position
          systemSize;
      }
      // output unwrapped coordinates in each dimension
      newPosition = newParticles[i].position()[dim] + increments[i][dim];
      outputFileStream.write((char*) &newPosition,
        sizeof(double));
    }

    // ORIENATION
    outputFileStream.write((char*) newParticles[i].orientation(),
      sizeof(double));
  }

  // COPYING

  for (int i=0; i < numberParticles; i++) {
    particles[i].copy(&(newParticles[i]));
  }
}
