#ifndef READ_H
#define READ_H

#include <fstream>
#include <string>

/////////////
// CLASSES //
/////////////


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

    int getFrames() const; // returns total number of frames

    void getActiveWork(double* work, int const& frame, int const& particle);
      // Returns work of particle between frames `frame' - 1 and `frame'.
    double getActiveWork(int const& frame, int const& particle);
      // Returns work of particle between frames `frame' - 1 and `frame'.
    double getPosition(
      int const& frame, int const& particle, int const& dimension,
      bool const& wrapped = true);
      // Returns position of a given particle at a given frame.
    double getOrientation(int const& frame, int const& particle);
      // Returns position of a given particle at a given frame.
    double getTimeStep(int const& frame);
      // Returns time step at a given frame.

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles
    double const persistenceLength; // persistence length
    double const packingFraction; // packing fraction
    double const systemSize; // size of the system
    int const randomSeed; // random seed

    std::ifstream inputFileStream; // input file stream
    int headerLength; // length of header in input file
    int frames; // number of frames

};

#endif
