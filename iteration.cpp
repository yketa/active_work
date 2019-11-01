#include <cmath>
#include <math.h>
#include <vector>

#include "iteration.h"
#include "particle.h"

void iterate_ABP_WCA(System *system) {
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

  std::vector<Particle> newParticles(system->getNumberParticles());

  // COMPUTATION
  for (int i=0; i < system->getNumberParticles(); i++) {

    // POSITIONS
    for (int dim=0; dim < 2; dim++) {
      newParticles[i].position()[dim] =
        (system->getParticle(i))->position()[dim]; // initialise new positions with previous ones
      newParticles[i].position()[dim] +=
        sqrt(system->getTimeStep()*2/3/system->getPersistenceLength())
          *(system->getRandomGenerator())->gauss_cutoff(); // add noise
      newParticles[i].position()[dim] +=
        system->getTimeStep()*cos(
            (system->getParticle(i))->orientation()[0] - dim*M_PI/2); // add self-propulsion

      (system->getParticle(i))->force()[dim] = 0; // initialise force
    }

    // ORIENTATIONS
    newParticles[i].orientation()[0] =
      (system->getParticle(i))->orientation()[0]; // initialise new orientation with previous one
    newParticles[i].orientation()[0] +=
      sqrt(system->getTimeStep()*2/system->getPersistenceLength())
        *(system->getRandomGenerator())->gauss_cutoff(); // add noise
  }

  // FORCES
  #if USE_CELL_LIST // with cell list

  int index1, index2; // index of the couple of particles
  int i, j; // indexes of the cells
  int k, l; // indexes of the particles in the cell
  std::vector<int> cell1, cell2;
  int numberBoxes = (system->getCellList())->getNumberBoxes();
  for (i=0; i < pow(numberBoxes, 2); i++) { // loop over cells

    cell1 = (system->getCellList())->getCell(i); // indexes of particles in the first cell
    for (k=0; k < cell1.size(); k++) { // loop over particles in the first cell
      index1 = cell1[k];

      // interactions with particles in the same cell
      for (l=k+1; l < cell1.size(); l++) { // loop over particles in the first cell
        index2 = cell1[l];
        system->WCA_potential(index1, index2, newParticles);
      }

      if ( numberBoxes == 1 ) { continue; } // only one cell

      // interactions with particles in other cells
      if ( numberBoxes == 2 ) { // 2 x 2 cells

        for (j=0; j < 4; j++) {
          if ( i == j ) { continue; } // same cell
          cell2 = (system->getCellList())->getCell(j); // indexes of particles in the second cell

          for (l=0; l < cell2.size(); l++) { // loop over particles in the second cell
            index2 = cell2[l];
            if ( index1 < index2 ) { // only count once each couple
              system->WCA_potential(index1, index2, newParticles);
            }
          }
        }
      }
      else { // 3 x 3 cells or more

        int x = i%numberBoxes;
        int y = i/numberBoxes;
        for (int dx=0; dx <= 1; dx++) {
          for (int dy=-1; dy < 2*dx; dy++) { // these two loops correspond to (dx, dy) = {0, -1}, {1, -1}, {1, 0}, {1, 1}, so that half of the neighbouring cells are explored
            j = (numberBoxes + (x + dx))%numberBoxes
              + numberBoxes*((numberBoxes + (y + dy))%numberBoxes); // index of neighbouring cell
            cell2 = (system->getCellList())->getCell(j); // indexes of particles in the second cell

            for (l=0; l < cell2.size(); l++) { // loop over particles in the second cell
              index2 = cell2[l];
              system->WCA_potential(index1, index2, newParticles);
            }
          }
        }
      }
    }
  }

  #else // without cell list

  for (int index1=0; index1 < system->getNumberParticles(); index1++) {
    for (int index2=index1+1; index2 < system->getNumberParticles(); index2++) {
      system->WCA_potential(index1, index2, newParticles);
    }
  }

  #endif

  // SAVE AND COPY
  system->saveNewState(newParticles);
}
