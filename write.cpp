#include <string>

#include "write.hpp"

/////////////
// CLASSES //
/////////////

/**********
 * OUTPUT *
 **********/

// CONSTRUCTORS

Output::Output(std::string filename) :
  outputFile(filename),
  outputStream(filename.c_str(), std::ios::out | std::ios::binary) {}

// DESTRUCTORS

Output::~Output() { outputStream.close(); }

// METHODS

std::string Output::getOutputFile() const { return outputFile; } // returns output file name

long int Output::tellp() { return outputStream.tellp(); } // returns position in output stream
