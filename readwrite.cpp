#include <string>

#include "readwrite.hpp"

/********
 * READ *
 ********/

// CONSTRUCTORS

Read::Read(std::string filename) :
  inputFile(filename),
  inputStream(filename.c_str(), std::ios::in | std::ios::binary),
  fileSize (0) {

  // get file size
  if ( inputStream ) {
    inputStream.seekg(0, std::ios_base::end);
    fileSize = tellg();
    inputStream.seekg(0, std::ios_base::beg);
  }
}

// DESTRUCTORS

Read::~Read() { inputStream.close(); }

// METHODS

std::string Read::getInputFile() const { return inputFile; } // returns input file name

long int Read::getFileSize() { return fileSize; } // returns size of file
long int Read::tellg() { return inputStream.tellg(); } // returns position in input stream


/*********
 * WRITE *
 *********/

// CONSTRUCTORS

Write::Write(std::string filename) :
  outputFile(filename),
  outputStream(filename.c_str(), std::ios::out | std::ios::binary) {}

// DESTRUCTORS

Write::~Write() { close(); }

// METHODS

std::string Write::getOutputFile() const { return outputFile; } // returns output file name

long int Write::tellp() { return outputStream.tellp(); } // returns position in output stream

void Write::flush() { outputStream.flush(); } // flush output file stream
void Write::close() { outputStream.close(); } // close output file stream
