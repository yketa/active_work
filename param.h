#ifndef PARAM_H
#define PARAM_H

/////////////
// CLASSES //
/////////////

class Parameters;

/*  PARAMETERS
 *  ----------
 *  Store parameters relative to a system.
 */

class Parameters {

  public:

    // CONSTRUCTORS

    Parameters(
      int N, double lp, double phi, double dt);

    // METHODS

    int getNumberParticles() const; // returns number of particles
    double getPersistenceLength() const; // returns dimensionless persistence length
    double getPackingFraction() const; // returns packing fraction
    double getSystemSize() const; // returns system size
    double getTimeStep() const; // returns time step

  private:

    // ATTRIBUTES

    int const numberParticles; // number of particles in the system
    double const persistenceLength; // dimensionless persistence length
    double const packingFraction; // packing fraction
    double const systemSize; // system size
    double const timeStep; // time step

};

#endif
