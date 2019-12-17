/*  Cloning algorithm from RLJ.
 *
 */

#ifndef CLONINGSERIAL_HPP
#define CLONINGSERIAL_HPP

#include <random>
#include <vector>
#include <chrono>

// System class
#include "particle.hpp"
// iterate_ABP_WCA function
#include "iteration.hpp"

// note this is called CloningSerial even if it also supports openMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class CloningSerial {
public:
  vector<System*> systems;  // these are the clones, the vector has size (2.nc)
  int nc;                 // how many clones
  int cloneMethod;        // this determines which clone selection method to use (default 2?)

  int tau; // number of time steps to simulate between each cloning step

  double outputPsi;          // this is the estimate of Psi
  vector<double> outputOP;   // averages over the trajectories of the different clones of (0) active work (1) force part of the active work (2) orientation part of the active work (3) order parameter
  double outputWalltime;     // time taken
  int finalOffset;           // used to access the final population at the end of the run

  // for output: this function is a pointer to system i in the final population (for this thread)
  System* finalSys(int i) { return systems[ finalOffset + i ]; }

  void setCloneMethod(int s) {
    cloneMethod = s;
    cout << "#cloningSerial: setting cloneMethod " << s << endl;
  }

  // random numbers
  mt19937 cloneTwister;

  // lightweight constructor
  CloningSerial(int nWork) :  nc (0) , cloneMethod (2), tau(nWork) {;} // the elementary number of step that we want to simulate (tau) has to be given here
  // simple destructor, we delete the systems but the vector deals with itself
  ~CloningSerial() { if (nc>0) { for (int i=0;i<2*nc;i++) delete systems[i]; } }

  // initialise systems array with 2.nc copies of the "dummy" input
  //   and gives a new seed to the random number generator on each copy
  // (also makes sure to clean up any old systems that are already in the array)
  // ... but note this does not initialise the state of the actual systems
  // tau corresponds
  void Init(int _nc, System* dummy, int masterSeed);

  // this runs the cloning for total time tmax
  void doCloning(double tmax, double sValue, int initSim = 1);

  // this part of the algorithm is a bit tricky so put it in a separate function
  void selectClones(int newClones[], double key[], int pullOffset);

  // this function is used in clone selection, could put it in an external class...
  int binsearch(double key[], double val, int keylength);
};

// NOTE: in the following all loops over clones should be parallelised via openMP

void CloningSerial::selectClones(int newClones[], double key[], int pullOffset) {
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
          double rr = (systems[pullOffset+i]->getRandomGenerator())->random01() * totups;
          //cout << rr << " " << totups << " " << offset+i << endl;
          newClones[i] = binsearch(key, rr, nc+1);
        }
      }
      else if ( cloneMethod == 2 ) {
        // this is the eq method (should be default)
        double alpha = (systems[pullOffset]->getRandomGenerator())->random01() * totups / nc;

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

// initialise list of clones
void CloningSerial::Init(int _nc, System* dummy, int masterSeed) {

  //cout << "#init CloningSerial" << endl;

  // delete any existing clones
  if (nc>0) { for (int i=0;i<2*nc;i++) delete systems[i]; }

  nc = _nc;
  systems.resize(2*nc); // = new SystemClass* [2*nc];

  cloneTwister.seed(masterSeed);
  // std::uniform_int_distribution<int> sdist(0, 10000000);  // 7 digit random number
  std::uniform_int_distribution<int> sdist(0, 0x7fffffff);  // 31-bit integer, avoids any ambiguity with signed/unsigned

  int processSeeds[2*nc];  // BIG array of seeds(!), generate them here in order to ensure reproducibility
  for (int i=0;i<2*nc;i++) processSeeds[i] = sdist(cloneTwister);

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for (int i=0;i<2*nc;i++) {
    systems[i] = new System(dummy, processSeeds[i], "", tau, false, 1, 0); // create new system from copy of dummyRun, with random seed from processSeeds, computing active work and order parameter for every tau iterations, not dumping to output file, and with 0 biasing parameter
    systems[i]->saveInitialState(); // save first frame (important for frame counting)
  }

  // just for extra safety we can opt to throw our random number generators out of sync here
  // (the risk is that in a very big population, several clones might have the same seed,
  //  the effect here is to increase the effective range of seeds by safetyFactor
  if (1) {
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
  }


}

// void CloningSerial::doCloning (horrible c++ template syntax)
void CloningSerial::doCloning(double tmax, double sValue, int initSim) {

  // !! this is the main cloning algorithm

  //cout << "#cloning: cloneMethod is " << cloneMethod << endl;

  // slightly fancy c++ clock
  std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

    // random initial condition for the nc systems that form the current population
    int arrswitch = 0;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0;i<nc;i++) {
      systems[i]->setBiasingParameter(0); // setting 0 biasing parameter (unmodified dynamics to sample initial configurations)
      for (int j=0; j < tau*initSim; j++) { iterate_ABP_WCA(systems[i]); } // simulate an elementary number of steps
      systems[i]->resetDump(); // reset dumps: important between different runs and to only count the relevant quantities within the cloning framework
    }

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<2*nc; i++) systems[i]->setBiasingParameter(sValue); // setting desired biasing parameter

    double lnX = 0.0;  // this is used in our final estimate of psi

    double actualTau = tau*systems[0]->getTimeStep();
    double sFactor = systems[0]->getNumberParticles()*actualTau;

    std::vector<double> weight;
    weight.resize(nc);

    // this is the main loop
    double iter;
    for (iter = 0; iter < tmax / (tau*systems[0]->getTimeStep()); iter += 1.0 ) // for each iteration of the algorithm
    {
      vector<double> upsilon(nc);  // these are the cloning factors

      // pushOffset refers to the current population, pullOffset will be the new population
      // (the systems array is of size 2.nc, only one half is "active" at each time)
      int pushOffset =  arrswitch   *nc;
      int pullOffset = (1-arrswitch)*nc;

      #ifdef DEBUG
      #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
      std::cout << "##torqueParameter " << systems[pushOffset]->getTorqueParameter() << std::endl;
      #endif
      #endif

      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int i = 0; i < nc; i++) //For each lattice in the current population
      {
          for (int k=0; k < tau; k++) { iterate_ABP_WCA(systems[pushOffset+i]); } // run dynamics
          #if CONTROLLED_DYNAMICS
          weight[i] = 1 - systems[pushOffset+i]->getBiasingParameter()/           // w = 1 - s/
            (3*systems[pushOffset+i]->getPersistenceLength())                     // (3*lp)
            + systems[pushOffset+i]->getWorkForce();                              // + w_f
          #if CONTROLLED_DYNAMICS == 2 || CONTROLLED_DYNAMICS == 3
          if (systems[pushOffset+i]->getTorqueParameter() != 0) {
            weight[i] += (systems[pushOffset+i]->getTorqueParameter()/            // w += (g
              systems[pushOffset+i]->getBiasingParameter())*                      // /s)
              (systems[pushOffset+i]->getPersistenceLength()*                     // (lp
              systems[pushOffset+i]->getTorqueIntegral0()                         // I_0
              - systems[pushOffset+i]->getTorqueIntegral1()                       // - I_1
              - systems[pushOffset+i]->getTorqueParameter()*                      // - g
              systems[pushOffset+i]->getPersistenceLength()*                      // lp
              systems[pushOffset+i]->getTorqueIntegral2());                       // I_2)
          }
          #endif
          upsilon[i] = exp(                                                       // cloning factor = exp(
            -systems[pushOffset+i]->getBiasingParameter() * sFactor               // - s*N*tau
            * weight[i]);                                                         // w)
          #else
          weight[i] = systems[pushOffset+i]->getWork();
          upsilon[i] = exp(                                                       // cloning factor = exp(
            -systems[pushOffset+i]->getBiasingParameter() * sFactor               // - s*N*tau
            * weight[i]);                                                         // w)
          #endif
      }

      #ifdef DEBUG
        // diagnostic printing (for debug)
        cout << "#logUps";
        for (int i=0;i<nc;i++) cout << " " << log( upsilon[i] );
        cout << endl;
        std::cout << "#w_mod";
        for (int i=0; i<nc; i++) std::cout << " " << weight[i];
        std::cout << std::endl;
      #endif

      // construct key based on the upsilon params
      vector<double> key(nc+1);
      key[0] = 0.0;
      for (int i = 0; i < nc; i++) {
        key[i + 1] = key[i] + upsilon[i];
      }
      double totups = key[nc]; //Extract total of upsilons

      //Calculate cloning factor and store as log to avoid very large values
      double X = double(totups) / double(nc);
      lnX = lnX + log(X);

      // RLJ: this uses the function pointer from above
      //selectclones(c, key, nc, totups, infos, arrswitch);

      // decide who to copy (from which old clone does each new clone "pull" its state)
      int newClones[nc];
      selectClones(newClones,key.data(),pullOffset);

      #ifdef DEBUG
        // diagnostic printing (for debug)
        cout << "#pull";
        for (int i=0;i<nc;i++) cout << " " << newClones[i] ;
        cout << endl;
      #endif

      // this is the actual cloning step
      // (at this point we use the = operator for the data member of SystemClass)
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i<nc; i++) {
        systems[pullOffset + i]->copyParticles(systems[ pushOffset + newClones[i] ]); // clone particles
        systems[pullOffset + i]->copyDump(systems[ pushOffset + newClones[i] ]); // clone dumps
      }

      // TORQUE PARAMETER VALUE

      // ORDER PARAMETER METHOD
      #if CONTROLLED_DYNAMICS == 2

      // order parameter squared
      double nusq = 0;
      #ifdef DEBUG
      double torqueIntegral0 = 0;
      #endif
      #ifdef _OPENMP
      #ifdef DEBUG
      #pragma omp parallel for reduction (+:nusq,torqueIntegral0)
      #else
      #pragma omp parallel for reduction (+:nusq)
      #endif
      #endif
      for (int i=0; i<nc; i++) {
        nusq += systems[pullOffset + i]->getTotalTorqueIntegral1()
          /systems[pullOffset + i]->getDump();
        #ifdef DEBUG
        torqueIntegral0 += systems[pullOffset + i]->getTotalTorqueIntegral0()
          /(systems[pullOffset + i]->getTimeStep()*systems[pullOffset + i]->getDump());
        #endif
      }
      #ifdef DEBUG
      torqueIntegral0 /= nc;
      std::cout << "##order parameter^2: " << nusq << " <I_0> " << torqueIntegral0 << std::endl;
      #endif

      // define and set g
      double g = 0;
      if (sValue != 0) { g = (1/(systems[0]->getNumberParticles()*nusq) - 1)
       /systems[0]->getPersistenceLength(); }
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i<nc; i++) {
        systems[pullOffset + i]->setTorqueParameter(g);
      }

      #endif

      // POLYNOMIAL METHOD
      #if CONTROLLED_DYNAMICS == 3

      // polynomial coefficients
      double torqueIntegral0 (0.0), torqueIntegral1 (0.0), torqueIntegral2 (0.0); // torque integrals
      double workForce (0.0); // force part of the normalised rate of active work
      #ifdef _OPENMP
      #pragma omp parallel for reduction (+:torqueIntegral0,torqueIntegral1,torqueIntegral2,workForce)
      #endif
      for (int i=0; i<nc; i++) {
        torqueIntegral0 += systems[pullOffset + i]->getTotalTorqueIntegral0()
          /(systems[pullOffset + i]->getTimeStep()*systems[pullOffset + i]->getDump());
        torqueIntegral1 += systems[pullOffset + i]->getTotalTorqueIntegral1()
          /systems[pullOffset + i]->getDump();
        torqueIntegral2 += systems[pullOffset + i]->getTotalTorqueIntegral2()
          /systems[pullOffset + i]->getDump();
        workForce += systems[pullOffset + i]->getTotalWorkForce()
          /(systems[pullOffset + i]->getTimeStep()*systems[pullOffset + i]->getDump());
      }
      torqueIntegral0 /= nc;
      torqueIntegral1 /= nc;
      torqueIntegral2 /= nc;
      workForce /= nc;
      double psi = double(lnX) / ((iter + 1.0)*actualTau) / systems[0]->getNumberParticles();
      #ifdef DEBUG
      std::cout << "##<I_0> = " << torqueIntegral0 << " <I_1> = " << torqueIntegral1 << " <I_2> = " << torqueIntegral2 << std::endl;
      std::cout << "##<w_f> = " << workForce << " SCGF = " << psi << std::endl;
      #endif
      double polyA = systems[0]->getPersistenceLength()*torqueIntegral2;
      double polyB = torqueIntegral1 - systems[0]->getPersistenceLength()*torqueIntegral0;
      double polyC = sValue*(sValue/systems[0]->getPersistenceLength()/3.0 - workForce - 1.0) - psi;
      double Delta = pow(polyB, 2) - 4.0*polyA*polyC;

      #ifdef DEBUG
      std::cout << "##torque polynimial: " << polyA << " g^2 + " << polyB << " g + " << polyC << std::endl;
      std::cout << "##order parameter^2: " << polyB + 1.0/systems[0]->getNumberParticles() << std::endl;
      if (Delta > 0) { std::cout << "#roots: " << (-polyB-sqrt(Delta))/(2.0*polyA) << " " << (-polyB+sqrt(Delta))/(2.0*polyA) << std::endl; }
      else if (Delta == 0) { std::cout << "#root: " << -polyB/(2.0*polyA) << std::endl; }
      else { std::cout << "#no real root" << std::endl; }
      #endif

      // define and set g
      double g = 0;
      #if 1 // highest root
      if (sValue != 0 && Delta >= 0) g = (-polyB + sqrt(Delta))/(2*polyA);
      #endif
      #if 0 // lowest root
      if (sValue != 0 && Delta >= 0) g = (-polyB - sqrt(Delta))/(2*polyA);
      #endif
      #if 0 // negative root closer to 0
      if (sValue != 0 && Delta >= 0) {
        g = -polyB;
        if (polyB > 0) {
          if (-polyB + sqrt(Delta) > 0) g = -polyB - sqrt(Delta);
          else g = -polyB + sqrt(Delta);
        }
        else {
          if (-polyB - sqrt(Delta) < 0) g = -polyB - sqrt(Delta);
        }
        g /= 2.0*polyA;
        g = (-polyB - sqrt(Delta))/2.0/polyA;
      }
      #endif
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int i=0; i<nc; i++) {
        systems[pullOffset + i]->setTorqueParameter(g);
      }

      #endif

      // we maintain outputOP as our estimate of the (time-integrated) OP over the current population
      // this is a bit inefficient but it should be a small overhead

      arrswitch = 1 - arrswitch; //Set the other set of systems to be used in next time loop

      #ifdef DEBUG
      std::cout << std::endl;
      #endif
    }

  finalOffset = arrswitch * nc;
  outputOP.assign(4, 0.0);
  for (int i=0; i < nc; i++) {
    outputOP[0] += systems[ finalOffset + i]->getTotalWork()
      /(systems[ finalOffset + i]->getTimeStep()*systems[ finalOffset + i]->getDump()); // normalised rate of active work
    outputOP[1] += systems[ finalOffset + i]->getTotalWorkForce()
      /(systems[ finalOffset + i]->getTimeStep()*systems[ finalOffset + i]->getDump()); // force part of the normalised rate of active work
    outputOP[2] += systems[ finalOffset + i]->getTotalWorkOrientation()
      /(systems[ finalOffset + i]->getTimeStep()*systems[ finalOffset + i]->getDump()); // orientation part of the normalised rate of active work
    outputOP[3] += systems[ finalOffset + i]->getTotalOrder()
      /systems[ finalOffset + i]->getDump(); // order parameter
  }

  for (unsigned int j=0;j<4;j++) { outputOP[j] /= nc; }

  outputPsi = double(lnX) / (iter*tau*systems[0]->getTimeStep());

  // C++ clocks again
  std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();

  // this horrible syntax just computes elapsed time in seconds, at microsecond resolution
  outputWalltime = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
}

// this is a binary search used in clone selection
int CloningSerial::binsearch(double *key, double val, int keylength) {
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
