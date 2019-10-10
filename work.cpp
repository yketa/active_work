#include <string>
#include <cmath>
#include <algorithm>

#include "work.h"
#include "read.h"

double* activeWork(std::string filename, int const& tau) {
  // Returns active work of particles within data file for periods of time
  // indicated by number of frames `tau'.

  Dat dat(filename); // data class to read from file

  double* workSums = (double*) malloc((dat.getFrames() - tau)*sizeof(double));
  for (int i=0; i < dat.getFrames() - tau; i++) { workSums[i] = 0; }

  double work;
  for (int frame=0; frame < dat.getFrames() - 1; frame++) { // loop over all frames
    for (int i=0; i < dat.getNumberParticles(); i++) { // loop over al particles

      // add to all relevant work sums
      dat.getActiveWork(&work, frame + 1, i);
      for (int t=std::max(0, frame + 1 - tau);
        t < std::min(frame, dat.getFrames() - tau);
        t++) {
          workSums[t] += work;
      }
    }
  }

  return workSums;
}
