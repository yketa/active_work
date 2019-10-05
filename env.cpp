#include <cstdlib>
#include <string>

#include "env.h"

double getEnvDouble(std::string const& envName, double const& defaultValue){
  // Returns environment variable as double.
  // (from https://stackoverflow.com/questions/5866134/how-to-read-linux-environment-variables-in-c)

  const char* val = std::getenv(envName.c_str());
  if ( val == 0 ){
    return defaultValue;
  }
  else {
    return std::atof(val);
  }
}

int getEnvInt(std::string const& envName, int const& defaultValue) {
  // Returns environment variable as integer.
  // (from https://stackoverflow.com/questions/5866134/how-to-read-linux-environment-variables-in-c)

  const char* val = std::getenv(envName.c_str());
  if ( val == 0 ){
    return defaultValue;
  }
  else {
//    return std::atoi(val);
    return (int) std::atof(val);
  }
}

std::string getEnvString(std::string const& envName,
  std::string const& defaultValue) {
  // Returns environment variable as string.
  // (from https://stackoverflow.com/questions/5866134/how-to-read-linux-environment-variables-in-c)

  const char* val = std::getenv(envName.c_str());
  if ( val == 0 ){
    return defaultValue;
  }
  else {
    return std::string(val);
  }
}
