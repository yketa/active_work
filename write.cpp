#include <string>

#include "write.h"

/////////////
// CLASSES //
/////////////

/**************
 * OUTPUT *
 **************/

// CONSTRUCTORS

Output::Output(std::string filename) :
  outputFile(filename),
  outputStream(filename.c_str(), std::ios::out | std::ios::binary) {}

// DESTRUCTORS

Output::~Output() { outputStream.close(); }

// METHODS

void Output::write(bool out) {
  if ( outputStream ) {
    outputStream.write((char*) &out, sizeof(bool));
  }
}
void Output::write(int out) {
  if ( outputStream ) {
    outputStream.write((char*) &out, sizeof(int));
  }
}
void Output::write(double out) {
  if ( outputStream ) {
    outputStream.write((char*) &out, sizeof(double));
  }
}

std::string Output::getOutputFile() const { return outputFile; } // returns output file name
