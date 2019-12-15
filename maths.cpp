#include "maths.hpp"

#include <math.h>

double getAngle(double cosinus, double signSinus) {
  // Returns angle in radians from its cosinus and sign of its sinus.

  double angle = acos(cosinus);
  angle *= signSinus > 0 ? 1 : -1;
  return angle;
}
