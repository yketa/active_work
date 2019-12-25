#ifndef ITERATION_HPP
#define ITERATION_HPP

#include "particle.hpp"

void iterate_ABP_WCA(System* system);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential, using custom dimensionless parameters
  // relations.

void iterate_ABP_WCA(System0* system);
  // Updates system to next step according to the dynamics of active Brownian
  // particles with WCA potential.

template<class SystemClass> void cellList_ABP_WCA(
  SystemClass* system, std::vector<Particle>& newParticles) {
  // Compute interactions with WCA potentials between all particles of the
  // system, using its cell list.

  int index1, index2; // index of the couple of particles
  int i, j; // indexes of the cells
  int k, l; // indexes of the particles in the cell
  std::vector<int> cell1, cell2;
  int numberBoxes = (system->getCellList())->getNumberBoxes();
  for (i=0; i < pow(numberBoxes, 2); i++) { // loop over cells

    cell1 = (system->getCellList())->getCell(i); // indexes of particles in the first cell
    for (k=0; k < (int) cell1.size(); k++) { // loop over particles in the first cell
      index1 = cell1[k];

      // interactions with particles in the same cell
      for (l=k+1; l < (int) cell1.size(); l++) { // loop over particles in the first cell
        index2 = cell1[l];
        system->WCA_potential(index1, index2, newParticles);
      }

      if ( numberBoxes == 1 ) { continue; } // only one cell

      // interactions with particles in other cells
      if ( numberBoxes == 2 ) { // 2 x 2 cells

        for (j=0; j < 4; j++) {
          if ( i == j ) { continue; } // same cell
          cell2 = (system->getCellList())->getCell(j); // indexes of particles in the second cell

          for (l=0; l < (int) cell2.size(); l++) { // loop over particles in the second cell
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

            for (l=0; l < (int) cell2.size(); l++) { // loop over particles in the second cell
              index2 = cell2[l];
              system->WCA_potential(index1, index2, newParticles);
            }
          }
        }
      }
    }
  }
}

template<class SystemClass> void bruteForce_ABP_WCA(
  SystemClass* system, std::vector<Particle>& newParticles) {
  // Compute interactions with WCA potentials between all particles of the
  // system, using a double loop.

  for (int index1=0; index1 < system->getNumberParticles(); index1++) {
    for (int index2=index1+1; index2 < system->getNumberParticles(); index2++) {
      system->WCA_potential(index1, index2, newParticles);
    }
  }
}

#endif
