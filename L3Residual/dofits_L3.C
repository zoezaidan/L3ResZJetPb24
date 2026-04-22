// Purpose: orchestration wrapper for the split L3 residual workflow.
//
// The active implementation now lives in:
//   1. L3Residual/L3Res.C
//   2. L3Residual/createL2L3ResTextFile.C
//
// This wrapper is optional: users can call the split macros directly, but it
// keeps one macro entry point for the full fit-plus-export chain.

#include "TString.h"

#include <iostream>
#include <string>

#include "L3Res.C"
#include "createL2L3ResTextFile.C"

using namespace std;

void dofits_L3(TString inFileL3Derived = "L3_derived.root",
               TString sampleTypesCSV = "photonjet",
               TString inputPtRangesCSV = "",
               string outfilename = "L3Res_photonjet",
               string runLabel = "2024ppRef",
               string lumiLabel = "pp 480.4 pb^{-1}",
               int refAlphaBin = 5,
               double fitAlphaMin = 0.0,
               double fitAlphaMax = 0.4,
               string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
               string outBaseDir = "L3Residual") {
  string normalizedBaseDir = outBaseDir;
  if (!normalizedBaseDir.empty() && normalizedBaseDir[0] != '/' && normalizedBaseDir[0] != '.') {
    normalizedBaseDir = "./" + normalizedBaseDir;
  }

  L3Res(inFileL3Derived,
        sampleTypesCSV,
        inputPtRangesCSV,
        outfilename,
        runLabel,
        lumiLabel,
        refAlphaBin,
        fitAlphaMin,
        fitAlphaMax,
        normalizedBaseDir);

  if (l2ResidualFile.empty()) {
    cout << "ERROR: The split export step requires an L2Residual text file." << endl;
    return;
  }

  TString fitRootPath;
  if (normalizedBaseDir.empty() || normalizedBaseDir == "." || normalizedBaseDir == "./") {
    fitRootPath = Form("%s_fit.root", outfilename.c_str());
  } else {
    fitRootPath = Form("%s/%s/%s_fit.root", normalizedBaseDir.c_str(), outfilename.c_str(), outfilename.c_str());
  }

  // Require an explicit local path; do not attempt silent fallbacks.
  cout << "INFO: Opening fit ROOT file: " << fitRootPath.Data() << endl;
  createL2L3ResTextFile(fitRootPath, l2ResidualFile);
}