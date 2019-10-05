#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "particle.h"
#include "iteration.h"
#include "env.h"

int main() {

  // VARIABLE DEFINITION

  // system
  double N = getEnvDouble("N", 1);
  double lp = getEnvDouble("LP", 1);
  double L = getEnvDouble("L", 2);
  int seed = getEnvInt("SEED", 1);
  std::string filename = getEnvString("FILE", "out.dat");

  // simulation
  double dt = getEnvDouble("DT", 1);
  int Niter = getEnvInt("NITER", 1000000);

  // SYSTEM

  System system(N, lp, L, seed, filename);
  system.saveInitialState();

  // ITERATION

  for (int i=0; i < Niter; i++) {
    iterate_ABP_WCA(&system, dt);
  }

}
