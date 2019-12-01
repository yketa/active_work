#include <math.h>

#include "param.hpp"

/////////////
// CLASSES //
/////////////

/**************
 * PARAMETERS *
 **************/

// CONSTRUCTORS

Parameters::Parameters() :
  numberParticles(0), persistenceLength(0), packingFraction(0), systemSize(0),
    timeStep(0) {}

Parameters::Parameters(int N, double lp, double phi, double dt) :
  numberParticles(N), persistenceLength(lp), packingFraction(phi),
    systemSize(sqrt(M_PI*N/phi)/2), timeStep(dt) {}

Parameters::Parameters(Parameters* parameters) :
  numberParticles(parameters->getNumberParticles()),
  persistenceLength(parameters->getPersistenceLength()),
  packingFraction(parameters->getPackingFraction()),
  systemSize(parameters->getSystemSize()),
  timeStep(parameters->getTimeStep()) {}

// METHODS

int Parameters::getNumberParticles() const { return numberParticles; }
double Parameters::getPersistenceLength() const {return persistenceLength; }
double Parameters::getPackingFraction() const { return packingFraction; }
double Parameters::getSystemSize() const { return systemSize; }
double Parameters::getTimeStep() const { return timeStep; }