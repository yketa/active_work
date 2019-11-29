#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "param.h"
#include "particle.h"
#include "iteration.h"
#include "env.h"

int main() {

  // VARIABLE DEFINITION

  // system
  int N = getEnvInt("N", 1); // number of particles in the system
  double lp = getEnvDouble("LP", 1); // dimensionless persistence length
  double phi = getEnvDouble("PHI", 0.02); // packing fraction
  int seed = getEnvInt("SEED", 1); // random seed
  std::string filename = getEnvString("FILE", "out.dat"); // output file name

  // simulation
  double dt = getEnvDouble("DT", 1); // time step
  int Niter = getEnvInt("NITER", 1000000); // number of iterations

  // active work computation
  int nWork = getEnvInt("NWORK", 0); // number of frames on which to compute active work

  // particle output
  bool dump = getEnvBool("DUMP", 1); // dump positions and orientations to output file
  int period = getEnvInt("PERIOD", 1); // period of dumping of positions and orientations in number of frames

  // PARAMETERS

  Parameters parameters(N, lp, phi, dt); // class of simulation parameters

  // SYSTEM

  System system(&parameters, seed, filename, nWork, dump, period);
  system.saveInitialState(); // save first frame

  // ITERATION

  for (int i=0; i < Niter; i++) { iterate_ABP_WCA(&system); }

}
