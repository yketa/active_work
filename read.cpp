#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "read.h"

/////////////
// CLASSES //
/////////////


// DAT

// CONSTRUCTORS

Dat::Dat(std::string filename) :
  inputFileStream(filename.c_str(), std::ios::in | std::ios::binary),
  numberParticles(), persistenceLength(), packingFraction(), systemSize(),
    randomSeed(), timeStep(), framesWork(), dumpParticles() {

  // HEADER INFORMATION
  inputFileStream.read((char*) &numberParticles, sizeof(int));
  inputFileStream.read((char*) &persistenceLength, sizeof(double));
  inputFileStream.read((char*) &packingFraction, sizeof(double));
  inputFileStream.read((char*) &systemSize, sizeof(double));
  inputFileStream.read((char*) &randomSeed, sizeof(int));
  inputFileStream.read((char*) &timeStep, sizeof(double));
  inputFileStream.read((char*) &framesWork, sizeof(int));
  inputFileStream.read((char*) &dumpParticles, sizeof(bool));

  // FILE PARTS LENGTHS
  headerLength = inputFileStream.tellg();
  particleLength = 3*sizeof(double)*dumpParticles;
  frameLength = numberParticles*particleLength;

  // ESTIMATION OF NUMBER OF COMPUTED WORK SUMS AND FRAMES
  inputFileStream.seekg(0, std::ios_base::end);
  fileSize = inputFileStream.tellg();
  numberWork = (fileSize - headerLength - frameLength)/(
    framesWork*frameLength + sizeof(double));
  frames = !dumpParticles ? 0 :
    (fileSize - headerLength - numberWork*sizeof(double))/frameLength;

  // FILE CORRUPTION CHECK
  if ( fileSize !=
    headerLength + frames*frameLength + numberWork*sizeof(double) ) {
    std::cout << "Invalid file size." << std::endl;
    exit(0);
  }

  // ACTIVE WORK
  double work;
  for (int i=0; i < numberWork; i++) {
    inputFileStream.seekg(
      headerLength // header
      + frameLength // frame with index 0
      + (1 + i)*framesWork*frameLength // all following packs of framesWork frames
      + i*sizeof(double)); // previous values of the active work
    inputFileStream.read((char*) &work, sizeof(double));
    activeWork.push_back(work);
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

double Dat::getPosition(
  int const& frame, int const& particle, int const& dimension) {
  // Returns position of a given particle at a given frame.

  inputFileStream.seekg(
    headerLength // header
    + frame*frameLength // other frames
    + particle*particleLength // other particles
    + (std::max(frame - 1, 0)/framesWork)*sizeof(double) // active work sums (taking into account the frame with index 0)
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
    + (std::max(frame - 1, 0)/framesWork)*sizeof(double) // active work sums (taking into account the frame with index 0)
    + 2*sizeof(double)); // positions
  double orientation;
  inputFileStream.read((char*) &orientation, sizeof(double));

  return orientation;
}
