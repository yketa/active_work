#ifndef READWRITE_HPP
#define READWRITE_HPP

#include <fstream>
#include <string>

/////////////
// CLASSES //
/////////////

class Read;
class Write;


/*  READ
 *  ----
 *  Read from binary input file.
 */

class Read {

  public:

    // CONSTRUCTORS

    Read(std::string filename);

    // DESTRUCTORS

    ~Read();

    // METHODS

    // READ WITHOUT OFFSET
    template<class InClass> void read(InClass* in) {
      // Read from input file.

      if ( inputStream ) {
        inputStream.read((char*) in, sizeof(InClass));
      }
    }
    template<class InClass> InClass read() {
      // Read from input file.

      InClass in (0);
      read<InClass>(&in);
      return in;
    }

    // READ WITH OFFSET
    template<class InClass> void read(InClass* in,
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Read from itput file at a given position given by offset relative to
      // way.

      if ( inputStream ) {
        inputStream.seekg(offset, way); // set position in input stream
        read<InClass>(in); // read
      }
    }
    template<class InClass> InClass read(
      long int offset, std::ios_base::seekdir way = std::ios_base::beg) {
      // Read from itput file at a given position given by offset relative to
      // way.

      InClass in (0);
      read<InClass>(&in, offset, way);
      return in;
    }

    std::string getInputFile() const; // returns input file name

    long int getFileSize(); // returns size of file
    long int tellg(); // returns position in input stream

  private:

    // ATTRIBUTES

    std::string const inputFile; // input file name
    std::ifstream inputStream; // input binary file stream
    long int fileSize; // size of file

};


/*  WRITE
 *  -----
 *  Write to binary output file.
 */

class Write {

  public:

    // CONSTRUCTORS

    Write(std::string filename);

    // DESTRUCTORS

    ~Write();

    // METHODS

    // READ WITH OFFSET
    template<class OutClass> void write(OutClass out) {
      // Write to output file.

      if ( outputStream ) {
        outputStream.write((char*) &out, sizeof(OutClass));
      }
    }

    // READ WITHOUT OFFSET
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
