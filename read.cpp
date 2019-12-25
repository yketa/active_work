#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "read.hpp"

/*******
 * DAT *
 *******/

// CONSTRUCTORS

Dat::Dat(std::string filename) :
  numberParticles(), persistenceLength(), packingFraction(), systemSize(),
    randomSeed(), timeStep(), framesWork(), dumpParticles(), dumpPeriod(),
  inputFileStream(filename.c_str(), std::ios::in | std::ios::binary) {

  // HEADER INFORMATION
  inputFileStream.read((char*) &numberParticles, sizeof(int));
  inputFileStream.read((char*) &persistenceLength, sizeof(double));
  inputFileStream.read((char*) &packingFraction, sizeof(double));
  inputFileStream.read((char*) &systemSize, sizeof(double));
  inputFileStream.read((char*) &randomSeed, sizeof(int));
  inputFileStream.read((char*) &timeStep, sizeof(double));
  inputFileStream.read((char*) &framesWork, sizeof(int));
  inputFileStream.read((char*) &dumpParticles, sizeof(bool));
  inputFileStream.read((char*) &dumpPeriod, sizeof(int));

  // FILE PARTS LENGTHS
  headerLength = inputFileStream.tellg();
  particleLength = 3*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;
  workLength = 4*sizeof(double);

  // ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER PARAMETER SUMS AND FRAMES
  inputFileStream.seekg(0, std::ios_base::end);
  fileSize = inputFileStream.tellg();
  numberWork = (fileSize - headerLength - frameLength)/(
    framesWork*frameLength*dumpParticles + workLength);
  frames = !dumpParticles ? 0 :
    (fileSize - headerLength - numberWork*workLength)/frameLength;

  // FILE CORRUPTION CHECK
  if ( fileSize !=
    headerLength + frames*frameLength + numberWork*workLength ) {
    std::cout << "Invalid file size." << std::endl;
    exit(0);
  }

  // ACTIVE WORK AND ORDER PARAMETER
  double work;
  for (int i=0; i < numberWork; i++) {
    inputFileStream.seekg(
      headerLength // header
      + frameLength // frame with index 0
      + (1 + i)*framesWork*frameLength // all following packs of framesWork frames
      + i*workLength); // previous values of the active work
    inputFileStream.read((char*) &work, sizeof(double));
    activeWork.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    activeWorkForce.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    activeWorkOri.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    orderParameter.push_back(work);
  }
}

// DESTRUCTORS

Dat::~Dat() { inputFileStream.close(); }

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

  inputFileStream.seekg(
    headerLength // header
    + frame*frameLength // other frames
    + particle*particleLength // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + dimension*sizeof(double)); // dimension
  double position;
  inputFileStream.read((char*) &position, sizeof(double));

  return position;
}

double Dat::getOrientation(int const& frame, int const& particle){
  // Returns position of a given particle at a given frame.

  inputFileStream.seekg(
    headerLength // header
    + frame*frameLength // other frames
    + particle*particleLength // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double)); // positions
  double orientation;
  inputFileStream.read((char*) &orientation, sizeof(double));

  return orientation;
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
  inputFileStream(filename.c_str(), std::ios::in | std::ios::binary) {

  // HEADER INFORMATION
  inputFileStream.read((char*) &numberParticles, sizeof(int));
  inputFileStream.read((char*) &potentialParameter, sizeof(double));
  inputFileStream.read((char*) &propulsionVelocity, sizeof(double));
  inputFileStream.read((char*) &transDiffusivity, sizeof(double));
  inputFileStream.read((char*) &rotDiffusivity, sizeof(double));
  inputFileStream.read((char*) &persistenceLength, sizeof(double));
  inputFileStream.read((char*) &packingFraction, sizeof(double));
  inputFileStream.read((char*) &systemSize, sizeof(double));
  inputFileStream.read((char*) &randomSeed, sizeof(int));
  inputFileStream.read((char*) &timeStep, sizeof(double));
  inputFileStream.read((char*) &framesWork, sizeof(int));
  inputFileStream.read((char*) &dumpParticles, sizeof(bool));
  inputFileStream.read((char*) &dumpPeriod, sizeof(int));

  // DIAMETERS
  double diameter;
  for (int i=0; i < numberParticles; i++) {
    inputFileStream.read((char*) &diameter, sizeof(double));
    diameters.push_back(diameter);
  }

  // FILE PARTS LENGTHS
  headerLength = inputFileStream.tellg();
  particleLength = 3*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;
  workLength = 4*sizeof(double);

  // ESTIMATION OF NUMBER OF COMPUTED WORK AND ORDER PARAMETER SUMS AND FRAMES
  inputFileStream.seekg(0, std::ios_base::end);
  fileSize = inputFileStream.tellg();
  numberWork = (fileSize - headerLength - frameLength)/(
    framesWork*frameLength*dumpParticles + workLength);
  frames = !dumpParticles ? 0 :
    (fileSize - headerLength - numberWork*workLength)/frameLength;

  // FILE CORRUPTION CHECK
  if ( fileSize !=
    headerLength + frames*frameLength + numberWork*workLength ) {
    std::cout << "Invalid file size." << std::endl;
    exit(0);
  }

  // ACTIVE WORK AND ORDER PARAMETER
  double work;
  for (int i=0; i < numberWork; i++) {
    inputFileStream.seekg(
      headerLength // header
      + frameLength // frame with index 0
      + (1 + i)*framesWork*frameLength // all following packs of framesWork frames
      + i*workLength); // previous values of the active work
    inputFileStream.read((char*) &work, sizeof(double));
    activeWork.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    activeWorkForce.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    activeWorkOri.push_back(work);
    inputFileStream.read((char*) &work, sizeof(double));
    orderParameter.push_back(work);
  }
}

// DESTRUCTORS

Dat0::~Dat0() { inputFileStream.close(); }

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

  inputFileStream.seekg(
    headerLength // header
    + frame*frameLength // other frames
    + particle*particleLength // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + dimension*sizeof(double)); // dimension
  double position;
  inputFileStream.read((char*) &position, sizeof(double));

  return position;
}

double Dat0::getOrientation(int const& frame, int const& particle){
  // Returns position of a given particle at a given frame.

  inputFileStream.seekg(
    headerLength // header
    + frame*frameLength // other frames
    + particle*particleLength // other particles
    + (std::max(frame - 1, 0)/framesWork)*workLength // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double)); // positions
  double orientation;
  inputFileStream.read((char*) &orientation, sizeof(double));

  return orientation;
}
