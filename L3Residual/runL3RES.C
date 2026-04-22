// Centralized wrapper for the full split L3 residual workflow.
// Similar in spirit to L2Residual/runL2RES.C.
//
// This wrapper runs:
// 1) Optional raw photon+jet response plots for input Data/MC files
// 2) Derivation of Data/MC response ratio vs pTref
// 3) Final L3 fit workflow + L3/L2L3 text export
//
// Output is written to:
//   L3Residual/L3_derived_2026_04_21_photonjet_full/
//   (including plots in /plots subdirectory)

#include "deriveL3_from_photonjet.C"
#include "plotresponse_L3.C"
#include "dofits_L3.C"

#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <string>

using namespace std;

void runL3RES(
  TString mcInputFile = "/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L3ResPhotonJet/2026_04_13_QCDPhoton_full.root",
  TString dataInputFile = "/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L3ResPhotonJet/2026_04_13_2024ppRefHP_full.root",
  string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
  bool makeInputPlots = true,
  int refAlphaBin = 5,
  double fitAlphaMin = 0.0,
  double fitAlphaMax = 0.4,
  TString sampleTypesCSV = "photonjet",
  TString inputPtRangesCSV = "60-300",
  string runLabel = "2024ppRef",
  string lumiLabel = "pp 480.4 pb^{-1}") {

  const string outTag = "L3_derived_2026_04_21_photonjet_full";
  const string outBaseDir = "L3Residual";
  const string outDir = outBaseDir + "/" + outTag;
  const string plotsDir = outDir + "/plots";
  const TString derivedFile = Form("%s/%s.root", outDir.c_str(), outTag.c_str());

  gSystem->mkdir(outDir.c_str(), kTRUE);
  gSystem->mkdir(plotsDir.c_str(), kTRUE);

  cout << "============================================" << endl;
  cout << "Running centralized L3 workflow" << endl;
  cout << "MC input:   " << mcInputFile << endl;
  cout << "Data input: " << dataInputFile << endl;
  cout << "Output dir: " << outDir << endl;
  cout << "Derived:    " << derivedFile << endl;
  cout << "============================================" << endl;

  if (makeInputPlots) {
    cout << "[1/3] Plotting input photon+jet response distributions" << endl;
    plotresponse_L3(dataInputFile, "Data_2026_04_21", false, "", runLabel.c_str(), lumiLabel.c_str(), plotsDir.c_str());
    plotresponse_L3(mcInputFile, "MC_2026_04_21", true, "", runLabel.c_str(), lumiLabel.c_str(), plotsDir.c_str());
    plotresponse_L3(dataInputFile, "DataVsMC_2026_04_21", false, mcInputFile, runLabel.c_str(), lumiLabel.c_str(), plotsDir.c_str());
  }

  cout << "[2/3] Deriving L3 response inputs" << endl;
  deriveL3_from_photonjet(mcInputFile, dataInputFile, derivedFile, refAlphaBin, false, true);

  cout << "[3/3] Running L3 fits and text export" << endl;
  dofits_L3(derivedFile,
            sampleTypesCSV,
            inputPtRangesCSV,
            outTag,
            runLabel,
            lumiLabel,
            refAlphaBin,
            fitAlphaMin,
            fitAlphaMax,
            l2ResidualFile,
            outBaseDir);

  cout << "============================================" << endl;
  cout << "L3 workflow finished." << endl;
  cout << "All outputs are in: " << outDir << endl;
  cout << "============================================" << endl;
}
