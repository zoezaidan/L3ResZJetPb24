// Purpose: orchestration wrapper for the split L3 residual workflow.
//
// The active implementation now lives in:
//   1. L3Residual/L3Res.C
//   2. L3Residual/createL2L3ResTextFile.C
//
// This wrapper is optional: users can call the split macros directly, but it
// keeps one macro entry point for the full fit-plus-export chain.

#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>

using namespace std;

static TString escapeRootString(const TString& value) {
  TString escaped = value;
  escaped.ReplaceAll("\\", "\\\\");
  escaped.ReplaceAll("\"", "\\\"");
  return escaped;
}

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

  gROOT->LoadMacro("L3Residual/L3Res.C");
  gROOT->LoadMacro("L3Residual/createL2L3ResTextFile.C");

  const TString fitCommand = Form("L3Res(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%.17g,%.17g,\"%s\")",
                                  escapeRootString(inFileL3Derived).Data(),
                                  escapeRootString(sampleTypesCSV).Data(),
                                  escapeRootString(inputPtRangesCSV).Data(),
                                  escapeRootString(outfilename.c_str()).Data(),
                                  escapeRootString(runLabel.c_str()).Data(),
                                  escapeRootString(lumiLabel.c_str()).Data(),
                                  refAlphaBin,
                                  fitAlphaMin,
                                  fitAlphaMax,
                                  escapeRootString(outBaseDir.c_str()).Data());
  gROOT->ProcessLine(fitCommand);

  if (l2ResidualFile.empty()) {
    cout << "ERROR: The split export step requires an L2Residual text file." << endl;
    return;
  }

  const TString exportCommand = Form("createL2L3ResTextFile(\"%s\",\"%s\")",
                                     escapeRootString(Form("%s/%s/%s_fit.root", outBaseDir.c_str(), outfilename.c_str(), outfilename.c_str())).Data(),
                                     escapeRootString(l2ResidualFile.c_str()).Data());
  gROOT->ProcessLine(exportCommand);
}