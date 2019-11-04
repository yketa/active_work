#ifndef READ_H
#define READ_H

#include <fstream>
#include <string>
#include <vector>

// DAT

class Dat {
  /*  Read files as defined by the System class (see particle.h).
   */

  public:

    // CONSTRUCTORS

    Dat(std::string filename);

    // DESTRUCTORS

    ~Dat();

    // METHODS

    int getNumberParticles() const; // returns number of particles
    double getPersistenceLength() const; // returns persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    int getRandomSeed() const; // returns random seed
    double getTimeStep() const; // returns time step
    int getFramesWork() const; // returns number of frames on which to sum the active work before dumping

    long int getNumberWork() const; // returns number of computed work sums
    long int getFrames() const; // returns number of frames

    std::vector<double> getActiveWork();
      // Returns vector of computed active work sums.
    std::vector<double> getActiveWorkForce();
      // Returns vector of computed active work (force) sums.
    std::vector<double> getActiveWorkOri();
      // Returns vector of computed active work (orientation) sums.

    double getPosition(
      int const& frame, int const& particle, int const& dimension);
      // Returns position of a given particle at a given frame.
    double getOrientation(int const& frame, int const& particle);
      // Returns position of a given particle at a given frame.

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles
    double const persistenceLength; // persistence length
    double const packingFraction; // packing fraction
    double const systemSize; // size of the system
    int const randomSeed; // random seed
    double const timeStep; // time step
    int const framesWork; // number of frames on which to sum the active work before dumping
    bool const dumpParticles; // positions and orientations dumped in file
    int const dumpPeriod; // period of dumping of positions and orientations in number of frames

    std::ifstream inputFileStream; // input file stream
    long int fileSize; // size of file

    long int headerLength; // length of header in input file
    long int particleLength; // length the data of a single particle takes in a frame
    long int frameLength; // length the data of a single frame takes in a file
    long int workLength; // length the data of a single work dump take in a file

    long int numberWork; // number of computed work sums
    long int frames; // number of frames

    std::vector<double> activeWork; // computed active work sums
    std::vector<double> activeWorkForce; // computed active work (force) sums
    std::vector<double> activeWorkOri; // computed active work (orientation) sums

};

#endif
