#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "dat.hpp"

/*******
 * DAT *
 *******/

// CONSTRUCTORS

Dat::Dat(std::string filename) :
  numberParticles(), persistenceLength(), packingFraction(), systemSize(),
    randomSeed(), timeStep(), framesWork(), dumpParticles(), dumpPeriod(),
  input(filename) {

  // HEADER INFORMATION
  input.read<const int>(&numberParticles);
  input.read<const double>(&persistenceLength);
  input.read<const double>(&packingFraction);
  input.read<const double>(&systemSize);
  input.read<const int>(&randomSeed);
  input.read<const double>(&timeStep);
  input.read<const int>(&framesWork);
  input.read<const bool>(&dumpParticles);
  input.read<const int>(&dumpPeriod);

  // FILE PARTS LENGTHS
  headerLength = input.tellg();
  particleLength = 5*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;
  workLength = 4*sizeof(double);

  // ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER PARAMETER SUMS AND FRAMES
  numberWork = (input.getFileSize() - headerLength - frameLength)/(
    framesWork*frameLength*dumpParticles + workLength);
  frames = !dumpParticles ? 0 :
    (input.getFileSize() - headerLength - numberWork*workLength)/frameLength;

  // FILE CORRUPTION CHECK
  if ( input.getFileSize() !=
    headerLength + frames*frameLength + numberWork*workLength ) {
    std::cout << "Invalid file size." << std::endl;
    exit(0);
  }

  // ACTIVE WORK AND ORDER PARAMETER
  double work;
  for (int i=0; i < numberWork; i++) {
    input.read<double>(&work,
      headerLength                     // header
      + frameLength                    // frame with index 0
      + (1 + i)*framesWork*frameLength // all following packs of framesWork frames
      + i*workLength);                 // previous values of the active work
    activeWork.push_back(work);
    input.read<double>(&work);
    activeWorkForce.push_back(work);
    input.read<double>(&work);
    activeWorkOri.push_back(work);
    input.read<double>(&work);
    orderParameter.push_back(work);
  }
}

// DESTRUCTORS

Dat::~Dat() {}

// METHODS

int Dat::getNumberParticles() const { return numberParticles; }
double Dat::getPersistenceLength() const { return persistenceLength; }
double Dat::getPackingFraction() const { return packingFraction; }
double Dat::getSystemSize() const { return systemSize; }
int Dat::getRandomSeed() const { return randomSeed; }
double Dat::getTimeStep() const { return timeStep; }
int Dat::getFramesWork() const { return framesWork; }

long int Dat::getNumberWork() const { return numberWork; }
long int Dat::getFrames() const { return frames; }

std::vector<double> Dat::getActiveWork() { return activeWork; }
std::vector<double> Dat::getActiveWorkForce() { return activeWorkForce; }
std::vector<double> Dat::getActiveWorkOri() { return activeWorkOri; }
std::vector<double> Dat::getOrderParameter() { return orderParameter; }

double Dat::getPosition(
  int const& frame, int const& particle, int const& dimension) {
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + dimension*sizeof(double));                     // dimension
}

double Dat::getOrientation(int const& frame, int const& particle){
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double));                             // positions
}

double Dat::getVelocity(
  int const& frame, int const& particle, int const& dimension) {
  // Returns velocity of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 3*sizeof(double)                               // positions and orientation
    + dimension*sizeof(double));                     // dimension
}


/********
 * DAT0 *
 ********/

// CONSTRUCTORS

Dat0::Dat0(std::string filename) :
  numberParticles(), potentialParameter(), propulsionVelocity(),
    transDiffusivity(), rotDiffusivity(), persistenceLength(),
    packingFraction(), systemSize(), randomSeed(), timeStep(), framesWork(),
    dumpParticles(), dumpPeriod(),
  input(filename) {

  // HEADER INFORMATION
  input.read<const int>(&numberParticles);
  input.read<const double>(&potentialParameter);
  input.read<const double>(&propulsionVelocity);
  input.read<const double>(&transDiffusivity);
  input.read<const double>(&rotDiffusivity);
  input.read<const double>(&persistenceLength);
  input.read<const double>(&packingFraction);
  input.read<const double>(&systemSize);
  input.read<const int>(&randomSeed);
  input.read<const double>(&timeStep);
  input.read<const int>(&framesWork);
  input.read<const bool>(&dumpParticles);
  input.read<const int>(&dumpPeriod);

  // DIAMETERS
  double diameter;
  for (int i=0; i < numberParticles; i++) {
    input.read<double>(&diameter);
    diameters.push_back(diameter);
  }

  // FILE PARTS LENGTHS
  headerLength = input.tellg();
  particleLength = 5*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;
  workLength = 4*sizeof(double);

  // ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER PARAMETER SUMS AND FRAMES
  numberWork = (input.getFileSize() - headerLength - frameLength)/(
    framesWork*frameLength*dumpParticles + workLength);
  frames = !dumpParticles ? 0 :
    (input.getFileSize() - headerLength - numberWork*workLength)/frameLength;

  // FILE CORRUPTION CHECK
  if ( input.getFileSize() !=
    headerLength + frames*frameLength + numberWork*workLength ) {
    std::cout << "Invalid file size." << std::endl;
    exit(0);
  }

  // ACTIVE WORK AND ORDER PARAMETER
  double work;
  for (int i=0; i < numberWork; i++) {
    input.read<double>(&work,
      headerLength                     // header
      + frameLength                    // frame with index 0
      + (1 + i)*framesWork*frameLength // all following packs of framesWork frames
      + i*workLength);                 // previous values of the active work
    activeWork.push_back(work);
    input.read<double>(&work);
    activeWorkForce.push_back(work);
    input.read<double>(&work);
    activeWorkOri.push_back(work);
    input.read<double>(&work);
    orderParameter.push_back(work);
  }
}

// DESTRUCTORS

Dat0::~Dat0() {}

// METHODS

int Dat0::getNumberParticles() const { return numberParticles; }
double Dat0::getPotentialParameter() const { return potentialParameter; }
double Dat0::getPropulsionVelocity() const { return propulsionVelocity; }
double Dat0::getTransDiffusivity() const { return transDiffusivity; }
double Dat0::getRotDiffusivity() const { return rotDiffusivity; }
double Dat0::getPersistenceLength() const { return persistenceLength; }
double Dat0::getPackingFraction() const { return packingFraction; }
double Dat0::getSystemSize() const { return systemSize; }
int Dat0::getRandomSeed() const { return randomSeed; }
double Dat0::getTimeStep() const { return timeStep; }
int Dat0::getFramesWork() const { return framesWork; }

long int Dat0::getNumberWork() const { return numberWork; }
long int Dat0::getFrames() const { return frames; }

std::vector<double> Dat0::getDiameters() { return diameters; }

std::vector<double> Dat0::getActiveWork() { return activeWork; }
std::vector<double> Dat0::getActiveWorkForce() { return activeWorkForce; }
std::vector<double> Dat0::getActiveWorkOri() { return activeWorkOri; }
std::vector<double> Dat0::getOrderParameter() { return orderParameter; }

double Dat0::getPosition(
  int const& frame, int const& particle, int const& dimension) {
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + dimension*sizeof(double));                     // dimension
}

double Dat0::getOrientation(int const& frame, int const& particle){
  // Returns position of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double));                             // positions
}

double Dat0::getVelocity(
  int const& frame, int const& particle, int const& dimension) {
  // Returns velocity of a given particle at a given frame.

  return input.read<double>(
    headerLength                                     // header
    + frame*frameLength                              // other frames
    + particle*particleLength                        // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 3*sizeof(double)                               // positions and orientation
    + dimension*sizeof(double));                     // dimension
}
