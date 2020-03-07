// Adapted cloning algorithm from RLJ.

#ifndef CLONINGSERIAL_HPP
#define CLONINGSERIAL_HPP

#include <random>
#include <vector>
#include <chrono>
#include <experimental/filesystem>

// note this is called CloningSerial even if it also supports openMP
// all loops over clones should be parallelised via openMP if supported
#ifdef _OPENMP
#include <omp.h>
#endif

#include "dat.hpp"
#include "particle.hpp"
#include "iteration.hpp"

/////////////
// CLASSES //
/////////////


/*  CLONINGSERIAL
 *  -------------
 *  Cloning algorithm.
 */

template<class SystemClass> class CloningSerial {
  /*  Contains all the parameters relevant to the cloning algorithm and the
   *  clones themselves.
   *  (see https://yketa.github.io/DAMTP_2019_Wiki/#Cloning%20algorithm)
   */

  public:

    // CONSTRUCTORS

    // lightweight constructor
    CloningSerial(
      int nClones, int nWork, int method = 2, std::string cDir = "") :
      nc (nClones), cloneMethod (method), tau (nWork), clonesDirectory(cDir)
      { systems.resize(2*nc); }

    // DESTRUCTORS

    // simple destructor, we delete the systems but the vector deals with itself
    ~CloningSerial() { deleteClones(); }

    // METHODS

    void init(SystemClass* dummy, int masterSeed); // initialise list of clones

    template<typename F, typename G, typename H> void doCloning(
      double tmax, double sValue, int initSim,
      F iterate, G getSWeight, H control,
      int runIndex = 0); // this runs the cloning for total time tmax
    void writeTrajFiles(int iter, Write& clonesLog, int runIndex);

    void selectClones(int newClones[], double key[], int pullOffset);

    int binsearch(double key[], double val, int keylength); // this is a binary search used in clone selection

    void deleteClones() { // delete clones
      if ( nc>0 ) { for (int i=0; i<2*nc; i++) delete systems[i]; }
    }

    SystemClass* finalSystem(int index) { return systems[finalOffset + index]; }

    // ATTRIBUTES

    std::vector<SystemClass*> systems; // these are the clones, the vector has size (2.nc)
    int const nc; // how many clones
    int const cloneMethod; // this determines which clone selection method to use (default 2?)

    int const tau; // number of simulation steps between cloning steps

    std::string const clonesDirectory; // if != "": directory where trajectories are saved

    double outputPsi; // this is the estimate of Psi
    std::vector<double> outputOP; // averages over the trajectories of the different clones of (0) active work (1) force part of the active work (2) orientation part of the active work (3) order parameter
    double outputWalltime; // time taken
    int finalOffset; // used to access the final population at the end of the run

    std::mt19937 cloneTwister; // random numbers

};

template<class SystemClass> void CloningSerial<SystemClass>::init(
  SystemClass* dummy, int masterSeed) {
  // initialise systems array with 2.nc copies of the "dummy" input
  // and gives a new seed to the random number generator on each copy
  // (also makes sure to clean up any old systems that are already in the array)
  // ... but note this does not initialise the state of the actual systems

  // std::cout << "#init CloningSerial" << std::endl;

  // delete any existing clones
  deleteClones();

  cloneTwister.seed(masterSeed);
  // std::uniform_int_distribution<int> sdist(0, 10000000);  // 7 digit random number
  std::uniform_int_distribution<int> sdist(0, 0x7fffffff);  // 31-bit integer, avoids any ambiguity with signed/unsigned

  int processSeeds[2*nc];  // BIG array of seeds(!), generate them here in order to ensure reproducibility
  for (int i=0;i<2*nc;i++) processSeeds[i] = sdist(cloneTwister);

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0;i<2*nc;i++) {
    systems[i] = new SystemClass(dummy, processSeeds[i], tau,
      clonesDirectory == "" ? "" : std::experimental::filesystem::path(
        std::experimental::filesystem::path(clonesDirectory) /
        [](int index)
          { return std::string(6 - std::to_string(index).length(), '0')
            + std::to_string(index) + std::string(".tmp.dat"); }
          (i)
      ).u8string()); // create new system from copy of dummy, with random seed from processSeeds, computing active work and order parameter for every tau iterations
  }

  // just for extra safety we can opt to throw our random number generators out of sync here
  // (the risk is that in a very big population, several clones might have the same seed,
  //  the effect here is to increase the effective range of seeds by safetyFactor
  #if 1
  for (int i=0;i<nc;i++) {
    int dum = 0;
    int safetyFactor = 1000;
    //int k = i % nc;
    if ( i==0 )
      std::cout << "#seed safetyFactor " << safetyFactor << std::endl;
    for (int r=0; r < i % safetyFactor; r++ ) {
      dum += (systems[i]->getRandomGenerator())->random01();
      dum += (systems[i+nc]->getRandomGenerator())->random01();
    }
  }
  #endif

  for (int i=0;i<2*nc;i++) systems[i]->saveInitialState(); // important if trajectories are actually saved

}

template<class SystemClass> template<typename F, typename G, typename H>
  void CloningSerial<SystemClass>::doCloning(
    double tmax, double sValue, int initSim,
    F iterate, G getSWeight, H control,
    int runIndex) {
  // this runs the cloning for total time tmax
  // input functions:
  // void iterate(SystemClass* system, int Niter): iterate system for Niter iterations
  // double getSWeight(SystemClass* system): returns product of biasing parameter and trajectory weight over last cloning step
  // void control(std::vector<SystemClass*> systems, int pullOffset, int pushOffset): modifications to controlled dynamics at the end of each cloning step

  // !! this is the main cloning algorithm

  // std::cout << "#cloning: cloneMethod is " << cloneMethod << std::endl;

  // clones log
  Write clonesLog(clonesDirectory == "" ? "" :
    std::experimental::filesystem::path
      (std::experimental::filesystem::path(clonesDirectory) / "clones.log")
    .u8string());

  // slightly fancy c++ clock
  std::chrono::system_clock::time_point startTime =
    std::chrono::system_clock::now();

  // random initial condition for the nc systems that form the current population
  int arrswitch = 0;
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0;i<nc;i++) {
    systems[i]->setBiasingParameter(0); // setting 0 biasing parameter (unmodified dynamics to sample initial configurations)
    iterate(systems[i], tau*initSim); // simulate an elementary number of steps
    systems[i]->resetDump(); // reset dumps: important between different runs and to only count the relevant quantities within the cloning framework
  }

  // biasing parameter
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0; i<2*nc; i++) systems[i]->setBiasingParameter(sValue); // setting desired biasing parameter

  double lnX = 0.0;  // this is used in our final estimate of psi

  std::vector<double> sWeight;
  sWeight.resize(nc);

  // this is the main loop
  int iter;
  for (iter = 0; iter < tmax / (tau*systems[0]->getTimeStep()); iter++) // for each iteration of the algorithm
  {
    std::vector<double> upsilon(nc);  // these are the cloning factors

    // pushOffset refers to the current population, pullOffset will be the new population
    // (the systems array is of size 2.nc, only one half is "active" at each time)
    int pushOffset =  arrswitch   *nc;
    int pullOffset = (1-arrswitch)*nc;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < nc; i++) // for each clone in the current population
    {
        iterate(systems[pushOffset+i], tau); // run dynamics

        sWeight[i] = getSWeight(systems[pushOffset+i]);
        upsilon[i] = exp(-sWeight[i]);
    }

    #ifdef DEBUG
      // diagnostic printing (for debug)
      std::cout << "#logUps";
      for (int i=0;i<nc;i++) std::cout << " " << log( upsilon[i] );
      std::cout << std::endl;
      // std::cout << "#s w_mod";
      // for (int i=0; i<nc; i++) std::cout << " " << sWeight[i];
      // std::cout << std::endl;
    #endif

    // construct key based on the upsilon params
    std::vector<double> key(nc+1);
    key[0] = 0.0;
    for (int i = 0; i < nc; i++) {
      key[i + 1] = key[i] + upsilon[i];
    }
    double totups = key[nc]; //Extract total of upsilons

    //Calculate cloning factor and store as log to avoid very large values
    double X = double(totups) / double(nc);
    lnX = lnX + log(X);

    // decide who to copy (from which old clone does each new clone "pull" its state)
    int newClones[nc];
    selectClones(newClones,key.data(),pullOffset);

    #ifdef DEBUG
      // diagnostic printing (for debug)
      std::cout << "#pull";
      for (int i=0;i<nc;i++) std::cout << " " << newClones[i] ;
      std::cout << std::endl;
    #endif
    for (int i=0;i<nc;i++) clonesLog.write<int>(newClones[i]);

    // this is the actual cloning step
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<nc; i++) {
      systems[pullOffset + i]->copyState(systems[ pushOffset + newClones[i] ]); // clone particles
      systems[pullOffset + i]->copyDump(systems[ pushOffset + newClones[i] ]); // clone dumps
    }

    // CONTROLLED DYNAMICS
    control(systems, pullOffset, pushOffset);

    arrswitch = 1 - arrswitch; //Set the other set of systems to be used in next time loop
  }

  finalOffset = arrswitch * nc;

  outputPsi = double(lnX) / (iter*tau*systems[0]->getTimeStep());

  // C++ clocks again
  std::chrono::system_clock::time_point endTime =
    std::chrono::system_clock::now();

  // this horrible syntax just computes elapsed time in seconds, at microsecond resolution
  outputWalltime = 1e-6*std::chrono::duration_cast<std::chrono::microseconds>
    (endTime - startTime).count();

  // WRITE TRAJECTORY FILES
  if ( clonesDirectory != "" ) writeTrajFiles(iter, clonesLog, runIndex);
}

template<class SystemClass> void CloningSerial<SystemClass>::writeTrajFiles(
  int iter, Write& clonesLog, int runIndex) {}

template<> void CloningSerial<Rotors>::writeTrajFiles(
  int iter, Write& clonesLog, int runIndex) {}

template<> void CloningSerial<System>::writeTrajFiles(
  int iter, Write& clonesLog, int runIndex) {
  // Write trajectory files for cloning loop.

  int period = 1;

  // BUILD TRAJECTORIES

  clonesLog.flush();
  Read log(clonesLog.getOutputFile());

  std::vector<std::vector<int>> trajectories (nc);
  std::vector<int> parents (nc);
  for (int i=0; i < nc; i++) { parents[i] = i; }

  for (int t=0; t < iter; t++) {
    for (int i=0; i < nc; i++) {
      int newParent = log.read<int>(
        -((t + 1)*nc - parents[i])*sizeof(int),
        std::ios_base::end);
      trajectories[i].insert(trajectories[i].begin(), newParent);
      parents[i] = newParent;
    }
  }

  // WRITE .dat FILES

  std::vector<Dat*> dat;
  std::vector<std::vector<double>> activeWork;
  std::vector<std::vector<double>> activeWorkForce;
  std::vector<std::vector<double>> activeWorkOri;
  std::vector<std::vector<double>> orderParameter;
  std::vector<std::vector<double>> torqueIntegral1;
  std::vector<std::vector<double>> torqueIntegral2;
  for (int i=0; i < 2*nc; i++) {
    systems[i]->flushOutputFile();
    dat.push_back(new Dat(systems[i]->getOutputFile(), true));
    activeWork.push_back(dat[i]->getActiveWork());
    activeWorkForce.push_back(dat[i]->getActiveWorkForce());
    activeWorkOri.push_back(dat[i]->getActiveWorkOri());
    orderParameter.push_back(dat[i]->getOrderParameter());
    torqueIntegral1.push_back(dat[i]->getTorqueIntegral1());
    torqueIntegral2.push_back(dat[i]->getTorqueIntegral2());
  }

  for (int i=0; i < nc; i++) {

    Write output(std::experimental::filesystem::path(
      std::experimental::filesystem::path(clonesDirectory) /
      [](int index, int rIndex)
        { return
          std::string(4 - std::to_string(rIndex).length(), '0')
            + std::to_string(rIndex) + std::string(".")
          + std::string(6 - std::to_string(index).length(), '0')
            + std::to_string(index) + std::string(".dat"); }
        (i, runIndex)
      ).u8string());

    // header
    output.write<int>(systems[i]->getNumberParticles());
    output.write<double>(systems[i]->getPersistenceLength());
    output.write<double>(systems[i]->getPackingFraction());
    output.write<double>(systems[i]->getSystemSize());
    output.write<double>(systems[i]->getTorqueParameter());
    output.write<int>(-1);
    output.write<double>(systems[i]->getTimeStep());
    output.write<int>(1);
    output.write<bool>(true);
    output.write<int>(tau);

    // initial frame
    int initFrame = (dat[i]->getFrames() - 1) - period*((int) (iter + 1)/2);
    int initClone = trajectories[i][0];
    for (int p=0; p < systems[i]->getNumberParticles(); p++) { // output all particles
      // POSITIONS
      for (int dim=0; dim < 2; dim++) { // output position in each dimension
        output.write<double>(dat[initClone]->getPosition(initFrame, p, dim));
      }
      // ORIENTATIONS
      output.write<double>(dat[initClone]->getOrientation(initFrame, p)); // output orientation
      // VELOCITIES
      for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
        output.write<double>(dat[initClone]->getVelocity(initFrame, p, dim)); // output velocity
      }
    }

    // other frames
    for (int frame=1; frame <= iter; frame++) {

      int currClone = trajectories[i][frame - 1] + nc*(1 - (frame%2));
      int currFrame = (dat[currClone]->getFrames() - 1)
        - ((int) (iter - frame)/2)*period;

      for (int p=0; p < systems[i]->getNumberParticles(); p++) { // output all particles
        // POSITIONS
        for (int dim=0; dim < 2; dim++) { // output position in each dimension
          output.write<double>(
            dat[currClone]->getPosition(currFrame, p, dim));
        }
        // ORIENTATIONS
        output.write<double>(
          dat[currClone]->getOrientation(currFrame, p)); // output orientation
        // VELOCITIES
        for (int dim=0; dim < 2; dim++) { // output velocity in each dimension
          output.write<double>(
            dat[currClone]->getVelocity(currFrame, p, dim)); // output velocity
        }
      }

      // active work, polarisation, and torque integrals
      output.write<double>(activeWork[currClone][currFrame/period - 1]);
      output.write<double>(activeWorkForce[currClone][currFrame/period - 1]);
      output.write<double>(activeWorkOri[currClone][currFrame/period - 1]);
      output.write<double>(orderParameter[currClone][currFrame/period - 1]);
      output.write<double>(torqueIntegral1[currClone][currFrame/period - 1]);
      output.write<double>(torqueIntegral2[currClone][currFrame/period - 1]);
    }

    // close file
    output.close();
  }
}

template<class SystemClass> void CloningSerial<SystemClass>::selectClones(
  int newClones[], double key[], int pullOffset) {

  std::uniform_real_distribution<double> randc(0,1.0);  // random double in [0,1]

  // this offset determines who is pulling (it is either 0 or nc)
  double totups = key[nc];

  if ( cloneMethod == 1 ) {
    // this is the iid method
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<nc; i++) {
      // RLJ: we use the system rng of the destination system to decide where to "pull" from
      double rr = (systems[pullOffset+i]->getRandomGenerator())->random01()
        *totups;
      //cout << rr << " " << totups << " " << offset+i << endl;
      newClones[i] = binsearch(key, rr, nc+1);
    }
  }
  else if ( cloneMethod == 2 ) {
    // this is the eq method (should be default)
    double alpha = (systems[pullOffset]->getRandomGenerator())->random01()
      *totups/nc;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<nc; i++) {
      double rr = alpha + (i*totups/(double)nc);
      newClones[i] = binsearch(key, rr, nc+1);
    }
  }
  // if not method 1 or 2 then just copy everyone equally (should never happen)
  else { for(int i=0; i<nc; i++) newClones[i] = i; }
}

template<class SystemClass> int CloningSerial<SystemClass>::binsearch(
  double *key, double val, int keylength) {
  // this is a binary search used in clone selection

  int l = 0; /*Left hand limit*/
  int r = keylength - 1; /*Right hand limit*/
  int m; /*Midpoint*/
  int k; /*Element containing value*/
  bool elementfound = false; /*True when the element containing the value has been found*/
  while (elementfound == false)
  {
      m = int(floor(float(l + r) / 2)); /*Calculate midpoint*/
      if (val<key[m]) /*If value lower than midpoint, shift the right hand limit*/
      {
          r = m;
      }
      else /*Otherwise shift the left hand limit*/
      {
          l = m;
      }
      if (l == r - 1) /*Value lies in element between limits*/
      {
          elementfound = true;
      }
  }
  k = r - 1; /*Element index is index of upper limit when disregarding the index of the 0 value at the beginning of the key*/
  return k;
}

#endif
