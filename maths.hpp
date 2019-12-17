#ifndef MATHS_HPP
#define MATHS_HPP

#include <random>
#include <cmath>

/////////////
// CLASSES //
/////////////

class rnd {
  // This simple rnd class is a wrapper for the built-in c++ random number
  // generator.
  // (from RLJ, modified for Gaussian with cut-off)

  private:
    // Nuts and bolts.. should not need to touch this.

    // ATTRIBUTES
    std::default_random_engine generator;
    int max;
    std::uniform_int_distribution<int> *intmax;
    std::uniform_real_distribution<double> *real01;
    std::normal_distribution<double> *normal;
    double g_cutoff; // cut-off for Gaussian white noise

  public:

    // CONTRUCTOR
    rnd(double g_co = 3) {
      max = 0x7fffffff;
      intmax = new std::uniform_int_distribution<int>(0, max);
      real01 = new std::uniform_real_distribution<double>(0.0, 1.0);
      normal = new std::normal_distribution<double>( 0.0 , 1.0 );
      g_cutoff = g_co;
    }
    // DESTRUCTOR
    ~rnd() { delete intmax; delete real01; delete normal; }

    // METHODS
    // set the random seed
    void   setSeed(int seed) { generator.seed(seed); }
    // member functions for generating random double in [0,1] and random integer in [0,max-1]
    double random01() { return (*real01)(generator); }
    int    randomInt(int max) { return (*intmax)(generator) % max; }
    double gauss() { return (*normal)(generator); }
    double gauss_cutoff() {
      double g = this->gauss();
      while (fabs(g) > g_cutoff) {
        g = this->gauss();
      }
      return g;
    }

};


////////////////
// PROTOTYPES //
////////////////

double getAngle(double cosinus, double signSinus);
  // Returns angle in radians from its cosinus and sign of its sinus.

double algDistPeriod(double const& x1, double const& x2, double const& length);
  // Returns algebraic distance from `x1' to `x2' on a line of length `length'
  // taking into account periodic boundary condition.

double dist2DPeriod(double* pos0, double* pos1, double const& length);
  // Returns distance between points on a plane, with positions `pos0' and
  // `pos1' taking into account period boundary condition in a square system
  // of size `length'.

#endif
