#include <string>

#include "read.h"

/////////////
// CLASSES //
/////////////


// DAT

// CONSTRUCTORS

Dat::Dat(std::string filename) :
  inputFileStream(filename.c_str(), std::ios::in | std::ios::binary),
  numberParticles(), persistenceLength(), packingFraction(), systemSize(),
    randomSeed() {

  // reading header
  inputFileStream.read((char*) &numberParticles, sizeof(int));
  inputFileStream.read((char*) &persistenceLength, sizeof(double));
  inputFileStream.read((char*) &packingFraction, sizeof(double));
  inputFileStream.read((char*) &systemSize, sizeof(double));
  inputFileStream.read((char*) &randomSeed, sizeof(int));

  // header length
  headerLength = inputFileStream.tellg();

  // number of frames
  inputFileStream.seekg(0, inputFileStream.end);
  frames = ((int) inputFileStream.tellg() - headerLength)/
    ((1 + 6*numberParticles)*sizeof(double));
}

// DESTRUCTORS

Dat::~Dat() { inputFileStream.close(); }

// METHODS

int Dat::getNumberParticles() const { return numberParticles; }
double Dat::getPersistenceLength() const { return persistenceLength; }
double Dat::getPackingFraction() const { return packingFraction; }
double Dat::getSystemSize() const { return systemSize; }
int Dat::getRandomSeed() const { return randomSeed; }

int Dat::getFrames() const { return frames; }

void Dat::getActiveWork(double* work, int const& frame, int const& particle) {
  // Returns work of particle between frames `frame' - 1 and `frame'.

  inputFileStream.seekg(headerLength + sizeof(double)*(
    frame*(1 + 6*numberParticles) // frames
    + 1 // time step
    + particle*6 // other particles
  ));
  inputFileStream.read((char*) work, sizeof(double));
}
double Dat::getActiveWork(int const& frame, int const& particle) {
  // Returns work of particle between frames `frame' - 1 and `frame'.

  inputFileStream.seekg(headerLength + sizeof(double)*(
    frame*(1 + 6*numberParticles) // frames
    + 1 // time step
    + particle*6 // other particles
  ));
  double work;
  inputFileStream.read((char*) &work, sizeof(double));

  return work;
}

double Dat::getPosition(
  int const& frame, int const& particle, int const& dimension,
  bool const& wrapped) {
  // Returns position of a given particle at a given frame.

  inputFileStream.seekg(headerLength + sizeof(double)*(
    frame*(1 + 6*numberParticles) // frames
    + 1 // time step
    + particle*6 // other particles
    + 1 // active work
    + (1 - wrapped)*2 // wrapped coordinates
    + dimension // dimension
  ));
  double position;
  inputFileStream.read((char*) &position, sizeof(double));

  return position;
}

double Dat::getOrientation(int const& frame, int const&particle){
  // Returns position of a given particle at a given frame.

  inputFileStream.seekg(headerLength + sizeof(double)*(
    frame*(1 + 6*numberParticles) // frames
    + 1 // time step
    + particle*6 // other particles
    + 1 // active work
    + 4 // position coordinates
  ));
  double orientation;
  inputFileStream.read((char*) &orientation, sizeof(double));

  return orientation;
}

double Dat::getTimeStep(int const& frame) {
  // Returns time step at a given frame.

  inputFileStream.seekg(headerLength + sizeof(double)*(
    frame*(1 + 6*numberParticles) // frames
  ));
  double timeStep;
  inputFileStream.read((char*) &timeStep, sizeof(double));

  return timeStep;
}
