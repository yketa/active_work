#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cloningserial.hpp"
#include "env.hpp"
#include "particle.hpp"
#include "readwrite.hpp"

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
	Write output(filename); // output class
	output.write<double>(tmax);
	output.write<int>(nc);
	output.write<double>(sValue);
	output.write<int>(seed);
	output.write<int>(nRuns);
	output.write<int>(cloneMethod);
	output.write<int>(initSim);
	output.write<int>(N);
	output.write<double>(lp);
	output.write<double>(phi);
	output.write<int>(tau);
	output.write<double>(dt);
	output.close();

	// save trajectories
	std::string clonesDirectory = getEnvString("CLONES_DIRECTORY", ""); // if different than "" then clones trajectories are saved to this directory

	// parameters class
	Parameters parameters(N, lp, phi, dt);
	// dummy system
	System dummy(&parameters);

	printf("## CloningSerial Code: tmax %.3e numClones %d runs %d s %.3e tau %d Delta.t %.3e\n",tmax,nc,nRuns,sValue,tau,dt);
	#if CONTROLLED_DYNAMICS == 1
	std::cout << "## Modified translational EOM." << std::endl;
	#endif
	#if CONTROLLED_DYNAMICS == 2
	std::cout << "## Modified translational and rotational EOM with 1st method." << std::endl;
	#endif
	#if CONTROLLED_DYNAMICS == 3
	std::cout << "## Modified translational and rotational EOM with 2nd method." << std::endl;
	#endif
	// see cloning.py for more info

	// cloning object (default cloning method is [eq])
	CloningSerial<System> clones(nc, tau, cloneMethod, clonesDirectory);

	// set up the clones etc, using dummySystem to get system sizes, hop rates, etc
	clones.init(&dummy, seed);
	std::cout << "## master seed " << seed << std::endl;

  double sFactor = N*tau*dt;

	for (int run = 0; run<nRuns;run++) {

		#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
		#ifdef TORQUE_DUMP
		Write torqueDump(getEnvString("TORQUE_DUMP_FILE", "torque.dump"));
		torqueDump.close();
		#endif
		#endif

		// go! (this includes generating "random" [different] initial conditions for the clones)
		clones.doCloning(tmax,sValue,initSim,

			// ITERATION FUNCTION
			[](System* system, int Niter) { iterate_ABP_WCA(system, Niter); },

			// GET WEIGHT FUNCTION
			[&sValue, &sFactor](System* system) {

			double sWeight;

			#if CONTROLLED_DYNAMICS
			sWeight = sValue*(                       // sw = s(
				1.0 - system->getBiasingParameter()/   // 1 - s/
				(3.0*system->getPersistenceLength())   // (3*lp)
				+ system->getWorkForce());             // + w_f))
			#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
			sWeight += system->getTorqueParameter()* // sw += g
				(1.0/system->getNumberParticles()      // (1/N
				- system->getTorqueIntegral1()         // - I_1
				- system->getTorqueParameter()*        // - g
				system->getPersistenceLength()*        // lp
				system->getTorqueIntegral2());         // I_2)
			#endif
			#else
			sWeight = sValue*system->getWork();
			#endif

			return sFactor*sWeight;
			},

			// CONTROL FUNCTION
			[&nc
			#if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
			#ifdef TORQUE_DUMP
			, &torqueDump
			#endif
			#endif
			](std::vector<System*>& systems, int pullOffset, int pushOffset) {

				// TORQUE PARAMETER VALUE

				// ORDER PARAMETER METHOD
				#if CONTROLLED_DYNAMICS == 2

				// order parameter squared
				double nusq = 0;
				#ifdef _OPENMP
				#pragma omp parallel for reduction (+:nusq)
				#endif
				for (int i=0; i<nc; i++) {
					nusq += systems[pushOffset + i]->getTotalTorqueIntegral1()
						/systems[pushOffset + i]->getDump();
				}
				nusq /= nc;

				// define and set g
				double g = (1.0/(systems[0]->getNumberParticles()*nusq) - 1.0)
					/systems[0]->getPersistenceLength();
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for (int i=0; i<nc; i++) {
					systems[pullOffset + i]->setTorqueParameter(g);
				}

				// output
				#ifdef TORQUE_DUMP
				torqueDump.open();
				torqueDump.write<double>(systems[pullOffset]->getTorqueParameter());
				torqueDump.close();
				#endif

				#endif

				// POLYNOMIAL METHOD
				#if CONTROLLED_DYNAMICS == 3

				// polynomial coefficients
				double torqueIntegral1 (0.0), torqueIntegral2 (0.0); // torque integrals
				double workForce (0.0); // force part of the normalised rate of active work
				#ifdef _OPENMP
				#pragma omp parallel for reduction (+:torqueIntegral1,torqueIntegral2,workForce)
				#endif
				for (int i=0; i<nc; i++) {
					torqueIntegral1 += systems[pushOffset + i]->getTotalTorqueIntegral1()
						/systems[pullOffset + i]->getDump();
					torqueIntegral2 += systems[pushOffset + i]->getTotalTorqueIntegral2()
						/systems[pullOffset + i]->getDump();
					workForce += systems[pushOffset + i]->getTotalWorkForce()
						/(systems[pushOffset + i]->getTimeStep()
							*systems[pushOffset + i]->getDump());
				}
				torqueIntegral1 /= nc;
				torqueIntegral2 /= nc;
				workForce /= nc;

				// define and set g
				double g = -(torqueIntegral1 - 1.0/systems[0]->getNumberParticles())
					/(2*systems[0]->getPersistenceLength()*torqueIntegral2);
				#ifdef _OPENMP
				#pragma omp parallel for
				#endif
				for (int i=0; i<nc; i++) {
					systems[pullOffset + i]->setTorqueParameter(g);
				}

				// output
				#ifdef TORQUE_DUMP
				torqueDump.open();
				torqueDump.write<double>(systems[pullOffset]->getTorqueParameter());
				torqueDump.close();
				#endif

				#endif
			},

			// RUN INDEX
			run
		);

		clones.outputOP.assign(4, 0.0);
		for (int i=0; i < nc; i++) {
			clones.outputOP[0] += (clones.finalSystem(i))->getTotalWork()
				/((clones.finalSystem(i))->getTimeStep()
					*(clones.finalSystem(i))->getDump()); // normalised rate of active work
			clones.outputOP[1] += (clones.finalSystem(i))->getTotalWorkForce()
				/((clones.finalSystem(i))->getTimeStep()
					*(clones.finalSystem(i))->getDump()); // force part of the normalised rate of active work
			clones.outputOP[2] += (clones.finalSystem(i))->getTotalWorkOrientation()
				/((clones.finalSystem(i))->getTimeStep()
					*(clones.finalSystem(i))->getDump()); // orientation part of the normalised rate of active work
			clones.outputOP[3] += (clones.finalSystem(i))->getTotalOrder()
				/(clones.finalSystem(i))->getDump(); // order parameter
		}

		for (unsigned int j=0;j<4;j++) { clones.outputOP[j] /= nc; }

		std::cout << std::endl;
		std::cout << "##s "    << sValue << std::endl
		          << "#SCGF "  << clones.outputPsi/N << std::endl
		          << "#w "     << clones.outputOP[0] << std::endl
		          << "#wf "    << clones.outputOP[1] << std::endl
							<< "#wo "		 << clones.outputOP[2] << std::endl
						  << "#nu "    << clones.outputOP[3] << std::endl << std::endl
		          << "##time " << clones.outputWalltime << std::endl;

		// output to file
		output.open();
		output.write<double>(clones.outputPsi);
		output.write<double>(clones.outputOP[0]);
		output.write<double>(clones.outputOP[1]);
		output.write<double>(clones.outputOP[2]);
		output.write<double>(clones.outputOP[3]);
		output.write<double>(clones.outputWalltime);
		output.close();
	}

}
