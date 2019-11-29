#ifndef WRITE_H
#define WRITE_H

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

    void write(bool out);
    void write(int out);
    void write(double out);

    std::string getOutputFile() const; // returns output file name

  private:

    // ATTRIBUTES

    std::string const outputFile; // output file name
    std::ofstream outputStream; // output binary file stream
      // WARNING: Content of file is erased.

};

#endif
