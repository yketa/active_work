/*  Cloning algorithm from RLJ.
 *
 */

#include <random>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cloningserial.hpp"
#include "env.hpp"
#include "param.hpp"
#include "particle.hpp"
#include "write.hpp"

using namespace std;

int main() {

	// cloning parameters
	double tmax = getEnvDouble("TMAX", 1); // dimensionless time to simulate
	int nc = getEnvInt("NC", 10); // number of clones
	double sValue = getEnvDouble("SVALUE", 0); // biasing parameter
	int seed = getEnvInt("SEED", 0); // master random seed
	int nRuns = getEnvInt("NRUNS", 1); // number of different runs
	int cloneMethod = 2; // should keep this set to 2 (!!)
	int initSim = getEnvInt("INITSIM", 1); // number of initial elementary number of iterations to "randomise" the systems

	// openMP parameters
	#ifdef _OPENMP
	int threads = getEnvInt("THREADS", -1); // number of threads
	printf("# compiled with openMP\n");
  if ( threads > 0 ) {
    printf("# setting threads %d\n",threads);
    omp_set_num_threads(threads);
  }
  printf("# running on %d threads\n",omp_get_max_threads());
	#endif

	// physical parameters
	int N = getEnvInt("N", 100); // number of particles in the system
  double lp = getEnvDouble("LP", 5); // dimensionless persistence length
  double phi = getEnvDouble("PHI", 0.65); // packing fraction

	// simulation parameters
	int tau = getEnvInt("TAU", 100); // elementary number of steps
	double dt = getEnvDouble("DT", 0.001); // time step

	// output to file
	std::string filename = getEnvString("FILE", ""); // output file name
	Output output(filename); // output class
	output.write(tmax);
	output.write(nc);
	output.write(sValue);
	output.write(seed);
	output.write(nRuns);
	output.write(cloneMethod);
	output.write(initSim);
	output.write(N);
	output.write(lp);
	output.write(phi);
	output.write(tau);
	output.write(dt);

	// parameters class
	Parameters parameters(N, lp, phi, dt);
	// dummy system
	System dummy(&parameters);

  printf("## CloningSerial Code: tmax %.3e numClones %d runs %d s %.3e tau %d Delta.t %.3e\n",tmax,nc,nRuns,sValue,tau,dt);
	#if CONTROLLED_DYNAMICS
	std::cout << "## Using controlled dynamics." << std::endl;
	#endif

  // cloning object (default cloning method is [eq])
  CloningSerial clones(tau);
  clones.setCloneMethod(cloneMethod);

  // set up the clones etc, using dummySystem to get system sizes, hop rates, etc
  clones.Init(nc,sValue,&dummy,seed);
  std::cout << "## master seed " << seed << std::endl;

  for (int run = 0; run<nRuns;run++) {

    // go! (this includes generating "random" [different] initial conditions for the clones)

    clones.doCloning(tmax,initSim);

    // output. Note OP here is the mean escape rate of the modified dynamics, is that the same as K(s) ?
    std::cout << "#psi_OP_t_s " << clones.outputPsi << " "
                         << clones.outputOP[0] << " "
                         << clones.outputOP[1] << " "
												 << clones.outputOP[2] << " "
												 << clones.outputOP[3] << " "
                         << clones.outputWalltime << " "
                         << sValue << std::endl;

		// output to file
		output.write(clones.outputPsi);
		output.write(clones.outputOP[0]);
		output.write(clones.outputOP[1]);
		output.write(clones.outputOP[2]);
		output.write(clones.outputOP[3]);
		output.write(clones.outputWalltime);
  }

}
