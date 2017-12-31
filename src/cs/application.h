// Copyright 2009, Andreas Biegert

#ifndef CS_APPLICATION_H_
#define CS_APPLICATION_H_

#include <stdio.h>
#include <stdlib.h>

#include <string>

#include "globals.h"
#include "getopt_pp.h"

using namespace GetOpt;

namespace cs {

// Basic (abstract) application class.
// Defines the high level behavior of an application. A new application is
// written by deriving a class from Application and writing an implementation
// of the Run() and maybe some other (like ParseOptions() etc.) methods.
class Application {
 public:
  // Register the application instance.
  Application();

  // Clean up the application settings, flush the diagnostic stream.
  virtual ~Application();

  // Main function (entry point) for the application.
  int main(int argc,                 // argc in a regular main
           char* argv[],             // argv in a regular main
           FILE* fout,               // output stream
           const std::string& name   // application name
           );

 protected:
  // Version number for usage output.
  static const char* kVersionNumber;
  // Copyright string for usage output.
  static const char* kCopyright;

  // Runs the application and return exit code. To be implemented by derived
  // classes.
  virtual int Run() = 0;
  // Parses command line options.
  virtual void ParseOptions(GetOpt_pp& /* options */) {};
  // Prints options summary to stream.
  virtual void PrintOptions() const {};
  // Prints usage banner to stream.
  virtual void PrintUsage() const {};
  // Prints short application description.
  virtual void PrintBanner() const {};
  // Prints copyright notification.
  void PrintHelp() const;

  static Application* instance_;       // current app. instance
  std::string         app_name_;       // application name
  std::string         log_level_;      // log reporting level
  std::string         log_file_;       // name of logfile
  FILE*               log_fp_;         // file pointer to logfile
  FILE*               out_;            // file pointer to output stream
};

}  // namespace cs

#endif  // CS_APPLICATION_H_
