#include <algorithm>

#include "env.hpp"
#include "iteration.hpp"
#include "particle.hpp"

int main() {

  // VARIABLE DEFINITION

  // system
  int N = getEnvInt("N", 1); // number of rotors in the system
  double Dr = getEnvDouble("DR", 1.0/2.0); // rotational diffusivity
  double g = getEnvDouble("G", -Dr/2); // aligning torque parameter
  int seed = getEnvInt("SEED", 1); // random seed
  std::string filename = getEnvString("FILE", "out.datR"); // output file name

  // simulation
  double dt = getEnvDouble("DT", 1e-3); // time step
  int Niter = getEnvInt("NITER", 1000000); // number of iterations

  // rotors output
  int period = getEnvInt("PERIOD", 1); // period of dumping of orientations in number of frames

  // SYSTEM

  Rotors rotors(N, Dr, g, dt, seed, filename, period); // define system
  rotors.saveState(); // save first frame

  // ITERATION

  iterate_rotors(&rotors, Niter);

}
