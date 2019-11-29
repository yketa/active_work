#include <math.h>

#include "param.h"

/////////////
// CLASSES //
/////////////

/**************
 * PARAMETERS *
 **************/

// CONSTRUCTORS

Parameters::Parameters(int N, double lp, double phi, double dt) :
  numberParticles(N), persistenceLength(lp), packingFraction(phi),
    systemSize(sqrt(M_PI*N/phi)/2), timeStep(dt) {}

// METHODS

int Parameters::getNumberParticles() const { return numberParticles; }
double Parameters::getPersistenceLength() const {return persistenceLength; }
double Parameters::getPackingFraction() const { return packingFraction; }
double Parameters::getSystemSize() const { return systemSize; }
double Parameters::getTimeStep() const { return timeStep; }
