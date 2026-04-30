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

#include "L3Res.C"
#include "createL2L3ResTextFile.C"
#include "deriveL3.C"
#include "plotresponse_L3.C"

#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include <string>

using namespace std;

void runL3RES(

    TString mcInputFile =
        "/eos/home-z/zzaidanc/L3ResZJetpp24/Outputs/jtpt40_Z40/2026_04_29_ZJet_MC_ak4_jtpt40_Z40.root",
    TString dataInputFile =
        "/eos/home-z/zzaidanc/L3ResZJetpp24/Outputs/jtpt40_Z40/2026_04_29_ZJet_Data_ak4_jtpt40_Z40.root",
    string l2ResidualFile =
        "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
    bool makeInputPlots = true, int refAlphaBin = 5, double fitAlphaMin = 0.0,
    double fitAlphaMax = 0.4, TString sampleTypesCSV = "photonjet",
    TString inputPtRangesCSV = "60-300", string runLabel = "2024ppRef",
    string lumiLabel = "pp 480.4 pb^{-1}") {

  // Make paths absolute to avoid ROOT path resolution issues
  TString currentDir = gSystem->pwd();
  const string currentDirStr = string(currentDir.Data());

  const string outTag = "L3_derived_2026_04_30_photonjet_final";
  const string outBaseDir = currentDirStr + "/L3Residual";
  const string outDir = outBaseDir + "/" + outTag;
  const string plotsDir = outDir + "/plots";

  gSystem->mkdir(outDir.c_str(), kTRUE);
  gSystem->mkdir(plotsDir.c_str(), kTRUE);

  const TString derivedFile =
      Form("%s/%s.root", outDir.c_str(), outTag.c_str());

  cout << "============================================" << endl;
  cout << "Running centralized L3 workflow" << endl;
  cout << "MC input:   " << mcInputFile << endl;
  cout << "Data input: " << dataInputFile << endl;
  cout << "Output dir: " << outDir << endl;
  cout << "Derived:    " << derivedFile << endl;
  cout << "============================================" << endl;

  if (makeInputPlots) {
    cout << "[1/3] Plotting input photon+jet response distributions" << endl;
    plotresponse_L3(dataInputFile, "Data", false, "", runLabel.c_str(),
                    lumiLabel.c_str(), plotsDir.c_str());
    plotresponse_L3(mcInputFile, "MC", true, "", runLabel.c_str(),
                    lumiLabel.c_str(), plotsDir.c_str());
    plotresponse_L3(dataInputFile, "DataVsMC", false, mcInputFile,
                    runLabel.c_str(), lumiLabel.c_str(), plotsDir.c_str());
  }

  cout << "[2/3] Deriving L3 response inputs" << endl;
  deriveL3(mcInputFile, dataInputFile, derivedFile, refAlphaBin,
                          false, true);

  cout << "[3/3] Running L3 fits and text export" << endl;
  L3Res(derivedFile, sampleTypesCSV, inputPtRangesCSV, outTag, runLabel,
        lumiLabel, refAlphaBin, fitAlphaMin, fitAlphaMax, outBaseDir);

  if (l2ResidualFile.empty()) {
    cout << "ERROR: The split export step requires an L2Residual text file."
         << endl;
    return;
  }

  string normalizedBaseDir = outBaseDir;
  if (!normalizedBaseDir.empty() && normalizedBaseDir[0] != '/' &&
      normalizedBaseDir[0] != '.') {
    normalizedBaseDir = "./" + normalizedBaseDir;
  }

  TString fitRootPath;
  if (normalizedBaseDir.empty() || normalizedBaseDir == "." ||
      normalizedBaseDir == "./") {
    fitRootPath = Form("%s_fit.root", outTag.c_str());
  } else {
    fitRootPath = Form("%s/%s/%s_fit.root", normalizedBaseDir.c_str(),
                       outTag.c_str(), outTag.c_str());
  }

  // Make the path absolute to avoid ROOT path resolution issues
  if (!fitRootPath.BeginsWith("/")) {
    fitRootPath = currentDir + "/" + fitRootPath;
  }

  cout << "INFO: Opening fit ROOT file: " << fitRootPath.Data() << endl;
  createL2L3ResTextFile(fitRootPath, l2ResidualFile);

  cout << "============================================" << endl;
  cout << "L3 workflow finished." << endl;
  cout << "All outputs are in: " << outDir << endl;
  cout << "============================================" << endl;
}
