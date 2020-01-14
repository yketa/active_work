#include <algorithm>

#include "env.hpp"
#include "fire.hpp"
#include "iteration.hpp"
#include "maths.hpp"
#include "particle.hpp"

int main() {

  // VARIABLE DEFINITION

  // system
  int N = getEnvInt("N", 1); // number of particles in the system
  double Dr = getEnvDouble("DR", 1.0/2.0); // rotational diffusivity
  double epsilon = getEnvDouble("EPSILON", Dr/3.0); // coefficient parameter of potential
  double v0 = getEnvDouble("V0", 1); // self-propulsion velocity
  double D = getEnvDouble("D", epsilon); // translational diffusivity
  double phi = getEnvDouble("PHI", 0.02); // packing fraction
  int seed = getEnvInt("SEED", 1); // random seed
  std::string filename = getEnvString("FILE", "out.dat0"); // output file name

  // diameters
  double I = getEnvDouble("I", 0); // polydispersity index
  std::vector<double> diameters (N); // array of diameters
  for (int i=0; i < N; i++) {
    diameters[i] = 1 - I + 2*I*i/(N - 1);
  }
  // randomisation of diameters order
  rnd randomGenerator;
  randomGenerator.setSeed(seed);
  std::random_shuffle(diameters.begin(), diameters.end(),
    [&randomGenerator](int max) { return randomGenerator.randomInt(max); });

  // simulation
  double dt = getEnvDouble("DT", 1e-3); // time step
  int Niter = getEnvInt("NITER", 1000000); // number of iterations

  // active work computation
  int nWork = getEnvInt("NWORK", 0); // number of frames on which to compute active work

  // particle output
  bool dump = getEnvBool("DUMP", 1); // dump positions and orientations to output file
  int period = getEnvInt("PERIOD", 1); // period of dumping of positions and orientations in number of frames

  // PARAMETERS

  Parameters parameters(N, epsilon, v0, D, Dr, phi, dt); // class of simulation parameters

  // SYSTEM

  System0 system(&parameters, diameters, seed, filename, nWork, dump, period); // define system
  FIRE_WCA(&system, // FIRE minimisation algorithm
    getEnvDouble("EMIN", 1),
    getEnvInt("ITERMAX", (int) 100.0/dt),
    getEnvDouble("DTMIN", dt*1e-3),
    getEnvDouble("DT0", dt*1e-1),
    getEnvDouble("DTMAX", dt));
  system.saveInitialState(); // save first frame

  // ITERATION

  iterate_ABP_WCA(&system, Niter);

}
