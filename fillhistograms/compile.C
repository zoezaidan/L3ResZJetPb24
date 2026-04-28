#include "TSystem.h"
#include <iostream>

void compile(bool clean = true) {

  if (clean) {
    std::cout << "Cleaning pre-existing compilation artifacts..." << std::endl;

    // Remove histograms.C artifacts
    gSystem->Unlink("histograms_C.so");
    gSystem->Unlink("histograms_C.d");
    gSystem->Unlink("histograms_C_ACLiC_dict_rdict.pcm");

    // Remove eventhistograms.C artifacts
    gSystem->Unlink("eventhistograms_C.so");
    gSystem->Unlink("eventhistograms_C.d");
    gSystem->Unlink("eventhistograms_C_ACLiC_dict_rdict.pcm");

    std::cout << "Cleanup complete." << std::endl;
  }

  std::cout << "Compiling histograms.C..." << std::endl;
  gROOT->ProcessLine(".L histograms.C+");

  std::cout << "Compiling eventhistograms.C..." << std::endl;
  gROOT->ProcessLine(".L eventhistograms.C+");

  std::cout << "Compilation complete!" << std::endl;
  gROOT->ProcessLine(".q");
}
