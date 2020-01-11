#ifndef WRITE_HPP
#define WRITE_HPP

#include <fstream>
#include <string>

/////////////
// CLASSES //
/////////////

class Output;


/*  OUTPUT
 *  ------
 *  Write to binary output file.
 */

class Output {

  public:

    // CONSTRUCTORS

    Output(std::string filename);

    // DESTRUCTORS

    ~Output();

    // METHODS

    template<class OutClass> void write(OutClass out) {
      // Write to output file.

      if ( outputStream ) {
        outputStream.write((char*) &out, sizeof(OutClass));
      }
    }

    template<class OutClass> void write(OutClass out,
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Write to output file at a given position given by offset relative to
      // way.

      if ( outputStream ) {
        outputStream.seekp(offset, way); // set position in output stream
        write<OutClass>(out); // write
        outputStream.seekp(0, std::ios_base::end); // set position to end of file
      }
    }

    std::string getOutputFile() const; // returns output file name

    long int tellp(); // returns position in output stream

  private:

    // ATTRIBUTES

    std::string const outputFile; // output file name
    std::ofstream outputStream; // output binary file stream
      // WARNING: Content of file is erased.

};

#endif
