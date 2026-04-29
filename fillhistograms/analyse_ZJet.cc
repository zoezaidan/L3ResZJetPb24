#include <iostream>
using std::cout;
using std::endl;

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <typeinfo>

#include "configurations.h"
#include "histograms.h"
#include "settings.h"

#include "JetCorrector.h"
#include "eventhistograms.h"
#include "helpers.h"
#include "input_config.h"
#include "chain_builder.h"

R__LOAD_LIBRARY(histograms_C.so)
R__LOAD_LIBRARY(eventhistograms_C.so)

float rho = 0.;

bool debug = false;
bool applyjetvetomap = true;

namespace {

std::string formatPercent(double numerator, double denominator) {
  if (!(denominator > 0.0)) {
    return "n/a";
  }

  std::ostringstream stream;
  stream << std::fixed << std::setprecision(2)
         << 100.0 * numerator / denominator << "%";
  return stream.str();
}

} // namespace

// Z candidate structure
struct ZCandidate {
  float pt, eta, phi, mass, rapidity;
  int lepton1_idx, lepton2_idx;
  float lepton1_pt, lepton1_eta, lepton1_phi;
  float lepton2_pt, lepton2_eta, lepton2_phi;
  bool isValid;

  ZCandidate() : pt(-1), eta(0), phi(0), mass(0), rapidity(0),
                 lepton1_idx(-1), lepton2_idx(-1),
                 lepton1_pt(0), lepton1_eta(0), lepton1_phi(0),
                 lepton2_pt(0), lepton2_eta(0), lepton2_phi(0),
                 isValid(false) {}
};

// Reconstruct Z from dimuons
ZCandidate ReconstructZFromMuons(Int_t nReco,
                                 std::vector<float>* recoPt,
                                 std::vector<float>* recoEta,
                                 std::vector<float>* recoPhi,
                                 std::vector<int>* recoCharge,
                                 std::vector<bool>* recoIDTight,
                                 std::vector<float>* recoPFChIso,
                                 std::vector<float>* recoPFPhoIso,
                                 std::vector<float>* recoPFNeuIso,
                                 std::vector<float>* recoPFPUIso) {
  ZCandidate z;
  if (nReco < 2) return z;
  float bestMassDiff = 999;
  float targetMass = 91.1880;  // Z mass Average in GeV (PDG 2024)

  // Loop over all opposite-sign muon pairs
  for (int i = 0; i < nReco; i++) {
    // Muon quality cuts
    if ((*recoPt)[i] < 20.0) continue;  // pT threshold
    if (fabs((*recoEta)[i]) > 2.4) continue;  // Acceptance
    // Tight muon ID check
    if (!(*recoIDTight)[i]) continue;
    // Muon 'i' Isolation Check ---
    float iso_i = (*recoPFChIso)[i] + std::max(0.0f, (*recoPFNeuIso)[i] + (*recoPFPhoIso)[i] - 0.5f * (*recoPFPUIso)[i]);
    if ((iso_i / (*recoPt)[i]) > 0.15) continue;
    for (int j = i+1; j < nReco; j++) {
      if ((*recoPt)[j] < 20.0) continue;
      if (fabs((*recoEta)[j]) > 2.4) continue;
      if (!(*recoIDTight)[j]) continue;
      // --- Muon 'j' Isolation Check ---
      float iso_j = (*recoPFChIso)[j] + std::max(0.0f, (*recoPFNeuIso)[j] + (*recoPFPhoIso)[j] - 0.5f * (*recoPFPUIso)[j]);
      if ((iso_j / (*recoPt)[j]) > 0.15) continue;
      // Opposite sign requirement
      if ((*recoCharge)[i] * (*recoCharge)[j] >= 0) continue;
      // Reconstruct Z 4-vector
      TLorentzVector mu1, mu2, Z;
      const float muon_mass = 0.1056583745;  // muon mass in GeV
      mu1.SetPtEtaPhiM((*recoPt)[i], (*recoEta)[i], (*recoPhi)[i], muon_mass);
      mu2.SetPtEtaPhiM((*recoPt)[j], (*recoEta)[j], (*recoPhi)[j], muon_mass);
      Z = mu1 + mu2;
      // Mass window: 60 < M < 120 GeV (wide for inclusive Z selection)
      if (Z.M() < 60.0 || Z.M() > 120.0) continue;
      // Select pair closest to Z mass
      float massDiff = fabs(Z.M() - targetMass);
      if (massDiff < bestMassDiff) {
        bestMassDiff = massDiff;
        z.pt = Z.Pt();
        z.eta = Z.Eta();
        z.phi = Z.Phi();
        z.mass = Z.M();
        z.rapidity = Z.Rapidity();
        z.lepton1_idx = i;
        z.lepton2_idx = j;
        z.lepton1_pt = (*recoPt)[i];
        z.lepton1_eta = (*recoEta)[i];
        z.lepton1_phi = (*recoPhi)[i];
        z.lepton2_pt = (*recoPt)[j];
        z.lepton2_eta = (*recoEta)[j];
        z.lepton2_phi = (*recoPhi)[j];
        z.isValid = true;
      }
    }
  }
  return z;
}

// Reconstruct Z from dielectrons
ZCandidate ReconstructZFromElectrons(Int_t nEle,
                                     std::vector<float>* elePt,
                                     std::vector<float>* eleEta,
                                     std::vector<float>* elePhi,
                                     std::vector<int>* eleCharge,
                                     std::vector<int>* eleCutIdWP80, 
                                     std::vector<float>* elePFRelIsoWithEA) {
  ZCandidate z;
  if (nEle < 2) return z;
  float bestMassDiff = 999;
  float targetMass = 91.1880;  // Z mass in GeV

  // Loop over all opposite-sign electron pairs
  for (int i = 0; i < nEle; i++) {
    // Electron quality cuts
    if ((*elePt)[i] < 20.0) continue;  // pT threshold
    if (fabs((*eleEta)[i]) > 2.5) continue;  // Acceptance (slightly wider than muons)

    // Tight electron ID check (WP80 = tight)
    //if (eleCutIdWP80 && (*eleCutIdWP80)[i] != 1) continue;

    // --- Electron 'i' Isolation Check ---
    if (elePFRelIsoWithEA && (*elePFRelIsoWithEA)[i] > 0.15) continue;
    for (int j = i+1; j < nEle; j++) {
      if ((*elePt)[j] < 20.0) continue;
      if (fabs((*eleEta)[j]) > 2.5) continue;
      //if (eleCutIdWP80 && (*eleCutIdWP80)[j] != 1) continue;
      if (elePFRelIsoWithEA && (*elePFRelIsoWithEA)[j] > 0.15) continue;
      // Opposite sign requirement
      if ((*eleCharge)[i] * (*eleCharge)[j] >= 0) continue;
      // Reconstruct Z 4-vector
      TLorentzVector ele1, ele2, Z;
      const float electron_mass = 0.000510998950;  // electron mass in GeV
      ele1.SetPtEtaPhiM((*elePt)[i], (*eleEta)[i], (*elePhi)[i], electron_mass);
      ele2.SetPtEtaPhiM((*elePt)[j], (*eleEta)[j], (*elePhi)[j], electron_mass);
      Z = ele1 + ele2;
      // Mass window: 60 < M < 120 GeV
      if (Z.M() < 60.0 || Z.M() > 120.0) continue;
      // Select pair closest to Z mass
      float massDiff = fabs(Z.M() - targetMass);
      if (massDiff < bestMassDiff) {
        bestMassDiff = massDiff;
        z.pt = Z.Pt();
        z.eta = Z.Eta();
        z.phi = Z.Phi();
        z.mass = Z.M();
        z.rapidity = Z.Rapidity();
        z.lepton1_idx = i;
        z.lepton2_idx = j;
        z.lepton1_pt = (*elePt)[i];
        z.lepton1_eta = (*eleEta)[i];
        z.lepton1_phi = (*elePhi)[i];
        z.lepton2_pt = (*elePt)[j];
        z.lepton2_eta = (*eleEta)[j];
        z.lepton2_phi = (*elePhi)[j];
        z.isValid = true;
      }
    }
  }
  return z;
}

ZCandidate ReconstructGenZFromMuons(Int_t nGen,
                                    std::vector<int>* genPID,
                                    std::vector<int>* genStatus,
                                    std::vector<float>* genPt,
                                    std::vector<float>* genEta,
                                    std::vector<float>* genPhi,
                                    std::vector<int>* genMotherID) {
  ZCandidate genZ;
  
  // 1. First check if the pointers exist
  if (!genPID || !genStatus || !genPt) return genZ;
  
  // 2. FORCE nGen to be the actual vector size, ignoring the garbage tree variable
  int safe_nGen = genPID->size();
  
  // 3. Now safely check if we have enough particles
  if (safe_nGen < 2) return genZ;

  float bestMassDiff = 999;
  float targetMass = 91.1880;
  const float muon_mass = 0.1056583745;

  // Find two opposite-sign final-state muons
  for (int i = 0; i < safe_nGen; i++) {
    // Must be muon (PID = ±13)
    if (abs((*genPID)[i]) != 13) continue;

    // Must be final state (status = 1)
    if ((*genStatus)[i] != 1) continue;

    // Optional: check mother is Z or intermediate muon
    if (genMotherID && i < (int)genMotherID->size()) {
      int momPID = abs((*genMotherID)[i]);
      // Accept if mother is Z(23), muon(13), or initial state
      if (momPID != 23 && momPID != 13 && momPID != 0 && momPID != 2212) continue;
    }

    if ((*genPt)[i] < 10.0) continue;  // Minimum pT

    for (int j = i+1; j < nGen; j++) {
      if (abs((*genPID)[j]) != 13) continue;
      if ((*genStatus)[j] != 1) continue;
      if ((*genPt)[j] < 10.0) continue;

      // Opposite sign check
      if ((*genPID)[i] * (*genPID)[j] >= 0) continue;

      // Reconstruct Z
      TLorentzVector mu1, mu2, Z;
      mu1.SetPtEtaPhiM((*genPt)[i], (*genEta)[i], (*genPhi)[i], muon_mass);
      mu2.SetPtEtaPhiM((*genPt)[j], (*genEta)[j], (*genPhi)[j], muon_mass);
      Z = mu1 + mu2;

      // Mass window check
      if (Z.M() < 60.0 || Z.M() > 120.0) continue;

      float massDiff = fabs(Z.M() - targetMass);
      if (massDiff < bestMassDiff) {
        bestMassDiff = massDiff;
        genZ.pt = Z.Pt();
        genZ.eta = Z.Eta();
        genZ.phi = Z.Phi();
        genZ.mass = Z.M();
        genZ.rapidity = Z.Rapidity();
        genZ.lepton1_idx = i;
        genZ.lepton2_idx = j;
        genZ.lepton1_pt = (*genPt)[i];
        genZ.lepton1_eta = (*genEta)[i];
        genZ.lepton1_phi = (*genPhi)[i];
        genZ.lepton2_pt = (*genPt)[j];
        genZ.lepton2_eta = (*genEta)[j];
        genZ.lepton2_phi = (*genPhi)[j];
        genZ.isValid = true;
      }
    }
  }

  return genZ;
}

// Reconstruct generator-level Z from gen-level electrons (ggHiNtuplizer tree)
ZCandidate ReconstructGenZFromElectrons(Int_t nMC,
                                        std::vector<int>* mcPID,
                                        std::vector<int>* mcStatus,
                                        std::vector<float>* mcPt,
                                        std::vector<float>* mcEta,
                                        std::vector<float>* mcPhi,
                                        std::vector<int>* mcMomPID) {
  ZCandidate genZ;
  if (nMC < 2 || !mcPID || !mcStatus || !mcPt) return genZ;

  float bestMassDiff = 999;
  float targetMass = 91.1880;
  
  // Identify if we are looking for muons (13) or electrons (11)
  // Based on the first lepton found in the event
  int targetPID = 11; 
  for(int k=0; k<nMC; k++) {
      if (abs((*mcPID)[k]) == 13) { targetPID = 13; break; }
      if (abs((*mcPID)[k]) == 11) { targetPID = 11; break; }
  }

  for (int i = 0; i < nMC; i++) {
    if (abs((*mcPID)[i]) != targetPID) continue;
    if ((*mcStatus)[i] != 1) continue;
    if ((*mcPt)[i] < 10.0) continue;

    for (int j = i+1; j < nMC; j++) {
      if (abs((*mcPID)[j]) != targetPID) continue; // FIXED: used to be hardcoded 11
      if ((*mcStatus)[j] != 1) continue;
      if ((*mcPt)[j] < 10.0) continue;
      if ((*mcPID)[i] * (*mcPID)[j] >= 0) continue;

      TLorentzVector l1, l2, Z;
      float m = (targetPID == 13) ? 0.105658 : 0.000511;
      l1.SetPtEtaPhiM((*mcPt)[i], (*mcEta)[i], (*mcPhi)[i], m);
      l2.SetPtEtaPhiM((*mcPt)[j], (*mcEta)[j], (*mcPhi)[j], m);
      Z = l1 + l2;

      if (Z.M() < 60.0 || Z.M() > 120.0) continue;

      float massDiff = fabs(Z.M() - targetMass);
      if (massDiff < bestMassDiff) {
        bestMassDiff = massDiff;
        genZ.pt = Z.Pt(); genZ.eta = Z.Eta(); genZ.phi = Z.Phi();
        genZ.mass = Z.M(); genZ.rapidity = Z.Rapidity();
        genZ.lepton1_pt = (*mcPt)[i]; genZ.lepton2_pt = (*mcPt)[j];
        genZ.isValid = true;
      }
    }
  }
  return genZ;
}

// Helper function to extract jet radius from tree path
// Examples: "ak4PFJetAnalyzer/t" -> 0.4, "ak2PFJetAnalyzer/t" -> 0.2
float ExtractJetRadius(const std::string& jetPath) {
    // Look for "ak" followed optionally by other chars, then a digit
    // Handles: ak4, akCs4, akPu4, etc.
    size_t pos = jetPath.find("ak");
    if (pos == std::string::npos) {
    log(LOG_WARNING, "Could not extract jet radius from path " + jetPath +
               "; using default radius = 0.4");
        return 0.4;
    }

    // Move past "ak"
    pos += 2;

    // Skip any additional characters (like "Cs", "Pu", etc.) until we find a digit
    while (pos < jetPath.size() && !isdigit(jetPath[pos])) {
        pos++;
    }

    if (pos >= jetPath.size() || !isdigit(jetPath[pos])) {
      log(LOG_WARNING, "Invalid jet radius in path " + jetPath +
                 "; using default radius = 0.4");
        return 0.4;
    }

    int radiusTimes10 = jetPath[pos] - '0';  // Convert char to int
    float radius = radiusTimes10 / 10.0;

    log(LOG_INFO, "Detected jet radius R = " + std::to_string(radius) +
              " from path " + jetPath);
    return radius;
}

// Z+Jet analysis for L3 residual corrections
// jetTree: jet tree path, e.g. "ak4PFJetAnalyzer/t"
// analysisType: ZJET_MUMU, ZJET_EE, or ZJET (combined)
void analyse_ZJet(string input = "ZJETHP",
                  string outputfiletag = "AK4_zjet",
                  bool isMC = false, bool checkjetid = false,
                  string inputType = "era", int maxFiles = -1,
                  int maxEvents = -1, string outputDir = "",
                  int batchIndex = -1, int totalBatches = 1,
                  string jetPath = "ak4PFJetAnalyzer/t",
                  AnalysisType analysisType = AnalysisType::ZJET_MUMU,
                  float jtptlimitforalpha = 15) {

  bool usecalotrig = false;
  bool checkvalidjet = false; // this is for checking valid jet range after applying l2. now for
                              // tightly limited range. TODO: do something smarter

  if (debug && g_verbosity < LOG_TRACE) {
    g_verbosity = LOG_TRACE;
  }
  
  if (!isZJetAnalysisType(analysisType)) {
      log(LOG_ERROR, "Invalid AnalysisType for analyse_ZJet: " +
             std::string(analysisTypeName(analysisType)));
      log(LOG_ERROR,
       "Valid options: AnalysisType::ZJET_MUMU, AnalysisType::ZJET_EE, AnalysisType::ZJET");
    return;
  }

  const bool useMuonFlavor = usesMuonZJetFlavor(analysisType);
  const bool useElectronFlavor = usesElectronZJetFlavor(analysisType);

  log(LOG_INFO, "Running Z+jet analysis with AnalysisType: " +
                    std::string(analysisTypeName(analysisType)));

  // Extract jet radius from tree path
  float jetRadius = ExtractJetRadius(jetPath);
  // Derived parameters based on jet radius
  float leptonJetDeltaR = jetRadius;  // dR cleaning between jet and leptons

  // Build input configuration
  InputConfig config;
  config.maxFiles = maxFiles;
  config.maxEvents = maxEvents;
  config.outputTag = outputfiletag;
  config.batchIndex = batchIndex;
  config.totalBatches = totalBatches;
  config.skipFiles = 0;

  // Set default output directory
  if (outputDir.empty()) {
    config.outputDir = "../output";
  } else {
    config.outputDir = outputDir;
  }

  // Determine input type and path
  if (inputType == "era") {
    auto it = filenames.find(input);
    if (it != filenames.end()) {
      config.type = InputType::FILE;
      config.path = it->second;
    } else {
      log(LOG_ERROR, "Unknown era: " + input);
      return;
    }
  } else if (inputType == "directory") {
    config.type = InputType::DIRECTORY;
    config.path = input;
  } else if (inputType == "filelist") {
    config.type = InputType::FILELIST;
    config.path = input;
  } else {
    config.type = InputType::FILE;
    config.path = input;
  }

  // Calculate skip for batch mode
  if (batchIndex >= 0 && totalBatches > 0 && maxFiles > 0) {
    config.skipFiles = batchIndex * maxFiles;
  }

  // Generate output filename
  string outputfilename;
  string inputName = input;
  // For directory/filelist, use last component of path as name
  // First trim any trailing '/' characters to avoid embedding the full path
  while (!inputName.empty() && inputName.back() == '/') inputName.pop_back();
  if (inputType != "era") {
    size_t lastSlash = inputName.find_last_of("/");
    if (lastSlash != string::npos && lastSlash < inputName.size() - 1) {
      inputName = inputName.substr(lastSlash + 1);
    }
  }
  // For filelist, strip the .txt extension if present
  if (inputType == "filelist") {
    size_t extPos = inputName.rfind(".txt");
    if (extPos != string::npos && extPos == inputName.size() - 4) {
      inputName = inputName.substr(0, extPos);
    }
  }

  // Update output tag to include jet radius
  string radiusTag = Form("_ak%.0f", 10*jetRadius);
  string fullOutputTag = outputfiletag + radiusTag;  // e.g., "zjet_R0.4"
  if (batchIndex >= 0) {
    // For batch mode, use just the outputfiletag (cleaner naming)
    outputfilename = Form("%s/%s_batch%d_of_%d.root",
                         config.outputDir.c_str(),
                         fullOutputTag.c_str(), batchIndex, totalBatches);
  } else {
    outputfilename = Form("%s/%s_%s.root", config.outputDir.c_str(),
                         inputName.c_str(), fullOutputTag.c_str());
  }
  if (debug)
    outputfilename = "test.root";

  log(LOG_INFO, "Output file: " + outputfilename);

  // Define tree paths
  std::string evtPath = "hiEvtAnalyzer/HiTree";
  std::string triggerPath = "hltanalysis/HltTree";
  std::string skimPath = "skimanalysis/HltTree";
  std::string electronPath = "ggHiNtuplizer/EventTree";
  std::string muonPath = "muonAnalyzer/MuonTree";

  log(LOG_INFO, "Using jet tree: " + jetPath);
  log(LOG_INFO, "Building input chains...");

  bool needElectronTree = (useElectronFlavor || isMC);
  bool needMuonTree = useMuonFlavor;

  TreeChains *chains = BuildChainsFromConfig(config, jetPath, needElectronTree, needMuonTree);
  if (!chains || chains->nEntries == 0) {
    log(LOG_ERROR, "No entries found in input");
    return;
  }

  auto evtTree = chains->evtChain;
  auto electronTree = chains->photonChain;  // For electrons
  auto muonTree = chains->muonChain;      // For muons
  auto jetTree = chains->jetChain;
  auto triggerTree = chains->triggerChain;
  auto skimTree = chains->skimChain;

  // Cuts and weights from event tree
  Int_t hiBin = -1;
  Float_t weight = 1, vz = 0, evtwt = 1;

  // Muon variables (vectors)
  Int_t nReco = 0;
  std::vector<float> *recoPt = 0;
  std::vector<float> *recoEta = 0;
  std::vector<float> *recoPhi = 0;
  std::vector<int> *recoCharge = 0;
  std::vector<bool> *recoIDTight = 0;  // Muon ID bits

  std::vector<float> *recoPFChIso = 0;
  std::vector<float> *recoPFPhoIso = 0;
  std::vector<float> *recoPFNeuIso = 0;
  std::vector<float> *recoPFPUIso = 0;

  // Electron variables (vectors)
  Int_t nEle = 0;
  std::vector<float> *elePt = 0;
  std::vector<float> *eleEta = 0;
  std::vector<float> *elePhi = 0;
  std::vector<int> *eleCharge = 0;
  std::vector<int> *eleCutIdWP80 = 0;  // Tight Electron ID
  std::vector<float> *elePFRelIsoWithEA = 0;

  // MC truth matching branches (only for MC)
  // Electron gen-info (from ggHiNtuplizer/EventTree)
  std::vector<int> *mcPID = 0;
  std::vector<int> *mcStatus = 0;
  std::vector<float> *mcPt = 0;
  std::vector<float> *mcEta = 0;
  std::vector<float> *mcPhi = 0;
  std::vector<float> *mcMass = 0;
  std::vector<int> *mcMomPID = 0;
  // Muon gen-info (from muonAnalyzer/MuonTree)
  Int_t nGen = 0;
  std::vector<int> *genPID = 0;
  std::vector<int> *genStatus = 0;
  std::vector<float> *genPt = 0;
  std::vector<float> *genEta = 0;
  std::vector<float> *genPhi = 0;
  std::vector<int> *genMotherID = 0;

  // Event tree setup
  // Now enable only the branches we need
  evtTree->SetBranchStatus("*", 0);
  evtTree->SetBranchStatus("hiBin", 1);
  evtTree->SetBranchStatus("vz", 1);
  if (isMC) {
    evtTree->SetBranchStatus("weight", 1);
  }
  // Set branch addresses BEFORE SetBranchStatus (like analyse.cc does)
  evtTree->SetBranchAddress("hiBin", &hiBin);
  evtTree->SetBranchAddress("vz", &vz);
  if (isMC) {
    evtTree->SetBranchAddress("weight", &weight);
  }

  //// EVENT FILTERS
  // auto skimTree = (TTree*)inFile->Get(skimPath.c_str());
  // if (!isMC) skimTree->SetBranchStatus("*",1);

  // Int_t pprimaryVertexFilter = 1;
  // if (!isMC) skimTree->SetBranchAddress("pprimaryVertexFilter",
  // &pprimaryVertexFilter);

  // Trigger branches
  Int_t trigger = 0;
  Int_t HLT_SingleMu = 1;
  Int_t HLT_SingleEle = 1;

  // auto triggerTree = (TTree*)inFile->Get(triggerPath.c_str());

  if (!isMC) {
    triggerTree->SetBranchStatus("*", 0);
    if (useMuonFlavor) {
      triggerTree->SetBranchStatus("HLT_PPRefL2SingleMu7_v1", 1);
      triggerTree->SetBranchAddress("HLT_PPRefL2SingleMu7_v1", &HLT_SingleMu); // SetBranchAddress does not take wildcards, so we need to specify the exact trigger name.
    }
    if (useElectronFlavor) {
      // TODO: Update trigger name based on actual trigger in forest
      log(LOG_WARNING,
          "Using placeholder dielectron trigger name - update this after checking the forest branches");
      // triggerTree->SetBranchStatus("HLT_SingleEle*", 1);
      // triggerTree->SetBranchAddress("HLT_SingleEle*", &HLT_SingleEle);
    }
  }

  // JETS
  // Disable all branches first, then enable only what we need
  jetTree->SetMakeClass(1); // Debug why this is necessary
  jetTree->SetBranchStatus("*", 0);

  Int_t evt;
  // Reconstruted jet information
  Int_t nref;
  Float_t jtpt[MAXJETS];
  Float_t jtpt_uncorr[MAXJETS];
  Float_t jteta[MAXJETS];
  Float_t jtphi[MAXJETS];

  Float_t jtnhf[MAXJETS];
  Float_t jtchf[MAXJETS];
  Float_t jtnef[MAXJETS];
  Float_t jtcef[MAXJETS];
  Float_t jtmuf[MAXJETS];
  Int_t jtchm[MAXJETS]; // charged multiplicity

  Float_t jtPfNHM[MAXJETS]; // Neutral Hadron Multiplicity
  Float_t jtPfCEM[MAXJETS]; // Charged EM Multiplicity
  Float_t jtPfNEM[MAXJETS]; // Neutral EM Multiplicity
  Float_t jtPfMUM[MAXJETS]; // Muon Multiplicity

  //Int_t jtn[MAXJETS];

  jetTree->SetBranchStatus("nref", 1);
  jetTree->SetBranchStatus("evt", 1);
  jetTree->SetBranchStatus("rawpt", 1);
  jetTree->SetBranchStatus("jteta", 1);
  jetTree->SetBranchStatus("jtphi", 1);

  jetTree->SetBranchAddress("evt", &evt);
  jetTree->SetBranchAddress("nref", &nref);
  jetTree->SetBranchAddress("rawpt", jtpt); // we want uncorrected rawpt, jtpt might have some JEC already applied
  jetTree->SetBranchAddress("jteta", jteta);
  jetTree->SetBranchAddress("jtphi", jtphi);

  jetTree->SetBranchStatus("jtPfNHF", 1);
  jetTree->SetBranchStatus("jtPfCHF", 1);
  jetTree->SetBranchStatus("jtPfNEF", 1);
  jetTree->SetBranchStatus("jtPfCEF", 1);
  jetTree->SetBranchStatus("jtPfMUF", 1);
  jetTree->SetBranchStatus("jtPfCHM", 1);

  jetTree->SetBranchAddress("jtPfNHF", jtnhf);
  jetTree->SetBranchAddress("jtPfCHF", jtchf);
  jetTree->SetBranchAddress("jtPfNEF", jtnef);
  jetTree->SetBranchAddress("jtPfCEF", jtcef);
  jetTree->SetBranchAddress("jtPfMUF", jtmuf);
  jetTree->SetBranchAddress("jtPfCHM", jtchm);

  jetTree->SetBranchStatus("jtPfNHM", 1);
  jetTree->SetBranchStatus("jtPfCEM", 1);
  jetTree->SetBranchStatus("jtPfNEM", 1);
  jetTree->SetBranchStatus("jtPfMUM", 1);

  jetTree->SetBranchAddress("jtPfNHM", jtPfNHM);
  jetTree->SetBranchAddress("jtPfCEM", jtPfCEM);
  jetTree->SetBranchAddress("jtPfNEM", jtPfNEM);
  jetTree->SetBranchAddress("jtPfMUM", jtPfMUM);

  // Gen level jet information
  Float_t jtpt_gen[MAXJETS];
  Float_t jteta_gen[MAXJETS];
  Float_t jtphi_gen[MAXJETS];
  Float_t refdrjt[MAXJETS];

  if (isMC) {
    jetTree->SetBranchStatus("refpt", 1);
    jetTree->SetBranchStatus("refeta", 1);
    jetTree->SetBranchStatus("refphi", 1);
    jetTree->SetBranchStatus("refdrjt", 1);

    jetTree->SetBranchAddress("refpt", jtpt_gen);
    jetTree->SetBranchAddress("refeta", jteta_gen);
    jetTree->SetBranchAddress("refphi", jtphi_gen);
    jetTree->SetBranchAddress("refdrjt", refdrjt);
  }

  // Set muon branch addresses (from muonAnalyzer/MuonTree)
  if (useMuonFlavor) {
    muonTree->SetBranchStatus("*", 0);
    muonTree->SetBranchStatus("nReco", 1);
    muonTree->SetBranchStatus("recoPt", 1);
    muonTree->SetBranchStatus("recoEta", 1);
    muonTree->SetBranchStatus("recoPhi", 1);
    muonTree->SetBranchStatus("recoCharge", 1);
    muonTree->SetBranchStatus("recoIDTight", 1);

    muonTree->SetBranchAddress("nReco", &nReco);
    muonTree->SetBranchAddress("recoPt", &recoPt);
    muonTree->SetBranchAddress("recoEta", &recoEta);
    muonTree->SetBranchAddress("recoPhi", &recoPhi);
    muonTree->SetBranchAddress("recoCharge", &recoCharge);
    muonTree->SetBranchAddress("recoIDTight", &recoIDTight);

    muonTree->SetBranchStatus("recoPFChIso", 1);
    muonTree->SetBranchStatus("recoPFPhoIso", 1);
    muonTree->SetBranchStatus("recoPFNeuIso", 1);
    muonTree->SetBranchStatus("recoPFPUIso", 1);

    muonTree->SetBranchAddress("recoPFChIso", &recoPFChIso);
    muonTree->SetBranchAddress("recoPFPhoIso", &recoPFPhoIso);
    muonTree->SetBranchAddress("recoPFNeuIso", &recoPFNeuIso);
    muonTree->SetBranchAddress("recoPFPUIso", &recoPFPUIso);
    // Gen-level muon branches (MC only)
    if (isMC) {
      muonTree->SetBranchStatus("nGen", 1);
      muonTree->SetBranchStatus("genPID", 1);
      muonTree->SetBranchStatus("genStatus", 1);
      muonTree->SetBranchStatus("genPt", 1);
      muonTree->SetBranchStatus("genEta", 1);
      muonTree->SetBranchStatus("genPhi", 1);
      muonTree->SetBranchStatus("genMotherID", 1);

      muonTree->SetBranchAddress("nGen", &nGen);
      muonTree->SetBranchAddress("genPID", &genPID);
      muonTree->SetBranchAddress("genStatus", &genStatus);
      muonTree->SetBranchAddress("genPt", &genPt);
      muonTree->SetBranchAddress("genEta", &genEta);
      muonTree->SetBranchAddress("genPhi", &genPhi);
      muonTree->SetBranchAddress("genMotherID", &genMotherID);
    }
  }

  // Set electron branch addresses (from ggHiNtuplizer/EventTree)
  if (useElectronFlavor) {
    electronTree->SetBranchStatus("*", 0);
    electronTree->SetBranchStatus("nEle", 1);
    electronTree->SetBranchStatus("elePt", 1);
    electronTree->SetBranchStatus("eleEta", 1);
    electronTree->SetBranchStatus("elePhi", 1);
    electronTree->SetBranchStatus("eleCharge", 1);
    electronTree->SetBranchStatus("eleCutIdWP80", 1);

    electronTree->SetBranchAddress("nEle", &nEle);
    electronTree->SetBranchAddress("elePt", &elePt);
    electronTree->SetBranchAddress("eleEta", &eleEta);
    electronTree->SetBranchAddress("elePhi", &elePhi);
    electronTree->SetBranchAddress("eleCharge", &eleCharge);
    electronTree->SetBranchAddress("eleCutIdWP80", &eleCutIdWP80);

    electronTree->SetBranchStatus("elePFRelIsoWithEA", 1);
    electronTree->SetBranchAddress("elePFRelIsoWithEA", &elePFRelIsoWithEA);
  }

  if (isMC) {
    // MC truth is in electronTree (ggHiNtuplizer/EventTree) - for electrons
    if (electronTree && (useElectronFlavor || isMC)) {
      if (!useElectronFlavor) { electronTree->SetBranchStatus("*", 0); }
      electronTree->SetBranchStatus("mcPID", 1);
      electronTree->SetBranchStatus("mcStatus", 1);
      electronTree->SetBranchStatus("mcMomPID", 1);
      electronTree->SetBranchStatus("mcPt", 1);
      electronTree->SetBranchStatus("mcEta", 1);
      electronTree->SetBranchStatus("mcPhi", 1);
      // Note: mcMass might not exist in all forests
      // electronTree->SetBranchStatus("mcMass", 1);

      electronTree->SetBranchAddress("mcPID", &mcPID);
      electronTree->SetBranchAddress("mcStatus", &mcStatus);
      electronTree->SetBranchAddress("mcMomPID", &mcMomPID);
      electronTree->SetBranchAddress("mcPt", &mcPt);
      electronTree->SetBranchAddress("mcEta", &mcEta);
      electronTree->SetBranchAddress("mcPhi", &mcPhi);
      // electronTree->SetBranchAddress("mcMass", &mcMass);
    }
  }

  TTreeCache::SetLearnEntries(10);
  evtTree->AddBranchToCache("*", true);
  triggerTree->AddBranchToCache("*", true);
  if (electronTree)
    electronTree->AddBranchToCache("*", true);
  if (muonTree)
    muonTree->AddBranchToCache("*", true);
  jetTree->AddBranchToCache("*", true);

  TFile *outfile = new TFile(outputfilename.c_str(), "RECREATE");

  // Local map for histogram storage (not global to avoid ROOT cleanup issues)
  map<string, vector<histograms *>> _histos;

  // Create folders for centrality bins. Mostly a placeholder in case PbPb
  // MC/data is checked.
  for (int j = 0; j < nhibins; ++j) {
    if (hibins[j] < hibins[j + 1]) {
      string name = Form("hibin_%.1f_%.1f", hibins[j], hibins[j + 1]);
      outfile->mkdir(name.c_str());

      TDirectory *dir = outfile->GetDirectory(name.c_str());
      assert(dir);
      dir->cd();

      for (int i = 0; i < netabins; ++i) {
        if (etaedges[i] < etaedges[i + 1]) {
          string name2 = Form("eta_%.1f_%.1f", etaedges[i], etaedges[i + 1]);
          dir->mkdir(name2.c_str());

          TDirectory *dir2 = dir->GetDirectory(name2.c_str());
          assert(dir2);
          dir2->cd();

          // Use ZJET analysis type
          histograms *h = new histograms(dir2, etaedges[i], etaedges[i + 1],
                                         hibins[j], hibins[j + 1], isMC, AnalysisType::ZJET);
          _histos[name2.c_str()].push_back(h);
        }
      }
    }
  }

  outfile->mkdir("event");
  TDirectory *dir = outfile->GetDirectory("event");
  assert(dir);
  dir->cd();

  eventhistograms *eh = new eventhistograms(dir, isMC);

#if REDOJES == 1
  log(LOG_INFO, "Applying MC JEC from file " + jecfile);
  vector<string> JECFiles;
  JECFiles.push_back(jecfile.c_str());
  if (!isMC) {
    log(LOG_INFO, "Applying L2 Residual from file " + l2file);
    JECFiles.push_back(l2file.c_str());
  }
  JetCorrector JEC(JECFiles);
#endif

  // JER to be added for L3 residual analysis

  // Jet veto map
  TFile *mapfile = nullptr;
  TH2D *vetomap = nullptr;
  if (applyjetvetomap) {
    mapfile = TFile::Open(vetomapFile, "READ");
    if (!mapfile || mapfile->IsZombie()) {
      log(LOG_WARNING,
          "Veto map unavailable, disabling veto map selection");
      if (mapfile) {
        mapfile->Close();
        delete mapfile;
        mapfile = nullptr;
      }
      applyjetvetomap = false;
    } else {
      vetomap = dynamic_cast<TH2D *>(mapfile->Get("jetvetomap_all"));
      if (!vetomap) {
        log(LOG_WARNING,
            "jetvetomap_all not found in veto map file, disabling veto map selection");
        mapfile->Close();
        delete mapfile;
        mapfile = nullptr;
        applyjetvetomap = false;
      }
    }
  }

  log(LOG_INFO, "Number of entries: " + std::to_string(chains->nEntries));
  Long64_t nentries = chains->nEntries;
  if (config.maxEvents > 0 && config.maxEvents < nentries) {
    nentries = config.maxEvents;
  }

  log(LOG_INFO, "Processing " + std::to_string(nentries) + " events");

  // Counters for monitoring
  Long64_t nEventsWithMuonZ = 0;
  Long64_t nEventsWithElectronZ = 0;
  Long64_t nEventsPassedSelection = 0;
  Long64_t nEventsWithGenZ = 0;
  Long64_t nEventsPassedWithGenZ = 0;

  Long64_t nEvents_hasJets = 0;
  Long64_t nEvents_ZptCut = 0;
  Long64_t nEvents_hasAwayJet = 0;
  Long64_t nEvents_passVeto = 0;

  Long64_t i_processed = 0;
  for (Long64_t i = 0; i < nentries; ++i) {
    if (debug) {
      log(LOG_DEBUG, "Event " + std::to_string(i) + ": begin");
    }
    evtTree->GetEntry(i);
    triggerTree->GetEntry(i);
    // Get entries from the trees required for the active AnalysisType.
    if (muonTree && useMuonFlavor) { muonTree->GetEntry(i); }
    if (electronTree && (useElectronFlavor || isMC)) { electronTree->GetEntry(i); }
    if (debug) {
      log(LOG_TRACE, "Event " + std::to_string(i) + ": lepton trees loaded");
    }
    
    log_progress_every(i + 1, nentries, 10000, LOG_INFO);
    if (i == 40000000) break;
    i_processed++;

    // Trigger logic
    if (!isMC) {
      if (analysisType == AnalysisType::ZJET_MUMU) {
        trigger = HLT_SingleMu;
      } else if (analysisType == AnalysisType::ZJET_EE) {
        trigger = HLT_SingleEle;
      } else if (analysisType == AnalysisType::ZJET) {
        trigger = (HLT_SingleMu || HLT_SingleEle);
      }
      if (!trigger) continue;
    } else {
      trigger = true;  // MC: no trigger requirement
    }

    evtwt = 1;
    if (isMC) {
      evtwt *= weight; // Just use generator weight, no pthat reweighting
    }

    // cout << weight << " " << evtwt << endl;

    // BASIC EVENT FILTERS
    //  if (!isMC) {
    //    skimTree->GetEntry(i);
    //    if (pprimaryVertexFilter != 1) continue;
    //  }
    jetTree->GetEntry(i);
    if (debug) {
      log(LOG_TRACE, "Event " + std::to_string(i) + ": jet tree loaded with nref=" +
                         std::to_string(nref) + ", nReco=" +
                         std::to_string(nReco) + ", nEle=" +
                         std::to_string(nEle));
    }

    // Reconstruct Z candidate(s)
    ZCandidate zMuon, zElectron, selectedZ;
    bool hasZMuon = false, hasZElectron = false;

    if (useMuonFlavor) {
      zMuon = ReconstructZFromMuons(nReco, recoPt, recoEta, recoPhi, recoCharge, recoIDTight, 
                              recoPFChIso, recoPFPhoIso, recoPFNeuIso, recoPFPUIso);
      hasZMuon = zMuon.isValid;
      if (hasZMuon) nEventsWithMuonZ++;
    }

    if (useElectronFlavor) {
      zElectron = ReconstructZFromElectrons(nEle, elePt, eleEta, elePhi, eleCharge, eleCutIdWP80, elePFRelIsoWithEA);
      hasZElectron = zElectron.isValid;
      if (hasZElectron) nEventsWithElectronZ++;
    }

    // Select which Z to use
    if (analysisType == AnalysisType::ZJET_MUMU) {
      if (!hasZMuon) continue;
      selectedZ = zMuon;
    } else if (analysisType == AnalysisType::ZJET_EE) {
      if (!hasZElectron) continue;
      selectedZ = zElectron;
    } else if (analysisType == AnalysisType::ZJET) {
      if (!hasZMuon && !hasZElectron) continue;
      if (hasZMuon && !hasZElectron) selectedZ = zMuon;
      else if (!hasZMuon && hasZElectron) selectedZ = zElectron;
      else {
        // Both valid - pick closest to Z mass
        float muonDiff = fabs(zMuon.mass - 91.1880);
        float eleDiff = fabs(zElectron.mass - 91.1880);
        selectedZ = (muonDiff < eleDiff) ? zMuon : zElectron;
      }
    }
    if (debug) {
      log(LOG_DEBUG, "Event " + std::to_string(i) + ": reco Z selected pT=" +
                         std::to_string(selectedZ.pt) + ", mass=" +
                         std::to_string(selectedZ.mass));
    }

    // ========================================
    // GEN-LEVEL Z MATCHING (MC only)
    // ========================================

    bool hasGenZ = false;
    ZCandidate genZ;
    float deltaR_recoGen = 999;
    float deltaPt_recoGen = 999;
    float deltaMass_recoGen = 999;

    if (isMC) {
      if (debug) {
        log(LOG_TRACE, "Event " + std::to_string(i) + ": entering gen-Z matching");
      }
      // For BOTH muon and electron channels, the gen info is in the mcPID branches (ggHiNtuplizer)
      //! TO BE UPDATED with the muonAnalyzer gen branches. 
      if (electronTree) {
        int nMC = mcPID ? mcPID->size() : 0;
        
        if (analysisType == AnalysisType::ZJET_MUMU) {
            // Re-using the electron function logic but looking for muons (PID 13)
            // We'll create a temp function or just modify yours:
            genZ = ReconstructGenZFromElectrons(nMC, mcPID, mcStatus, mcPt, mcEta, mcPhi, mcMomPID);
            
            /* Note: If your ReconstructGenZFromElectrons function is strictly hardcoded 
               to PID 11, you should change that function's loop to:
              int targetPID =
                 (analysisType == AnalysisType::ZJET_MUMU) ? 13 : 11;
               if (abs((*mcPID)[i]) != targetPID) continue; 
            */
        } else {
            genZ = ReconstructGenZFromElectrons(nMC, mcPID, mcStatus, mcPt, mcEta, mcPhi, mcMomPID);
        }
      }

      // Calculate matching metrics if gen Z was found
      if (genZ.isValid) {
        hasGenZ = true;

        // Calculate ΔR between reco and gen Z
        float dEta = selectedZ.eta - genZ.eta;
        float dPhi = selectedZ.phi - genZ.phi;
        if (dPhi > TMath::Pi()) dPhi -= 2 * TMath::Pi();
        if (dPhi < -TMath::Pi()) dPhi += 2 * TMath::Pi();
        deltaR_recoGen = sqrt(dEta * dEta + dPhi * dPhi);

        deltaPt_recoGen = selectedZ.pt - genZ.pt;
        deltaMass_recoGen = selectedZ.mass - genZ.mass;

        // Count gen-matched events
        nEventsWithGenZ++;
      }
    }
    if (debug) {
      log(LOG_DEBUG, "Event " + std::to_string(i) + ": gen matching done hasGenZ=" +
                         std::to_string(hasGenZ) + ", deltaR=" +
                         std::to_string(deltaR_recoGen));
    }

    // Optional: Apply gen-matching requirement for MC validation
    // Uncomment these lines if you want to require gen-matching:
    if (isMC && !hasGenZ) continue;
    if (isMC && deltaR_recoGen > 0.3) continue;  // Require dR < 0.3

    // Debug output for first few events
    if (debug && isMC && i < 5) {  // Print first 5 MC events
      cout << "\n=== Event " << i << " Gen-Matching Debug ===" << endl;
      cout << "AnalysisType: " << analysisTypeName(analysisType) << endl;
      cout << "Reco Z: pT=" << selectedZ.pt << " GeV, eta=" << selectedZ.eta 
           << ", phi=" << selectedZ.phi << ", mass=" << selectedZ.mass << " GeV" << endl;

      if (hasGenZ) {
        cout << "Gen Z:  pT=" << genZ.pt << " GeV, eta=" << genZ.eta 
             << ", phi=" << genZ.phi << ", mass=" << genZ.mass << " GeV" << endl;
        cout << "Matching: ΔR=" << deltaR_recoGen << ", ΔpT=" << deltaPt_recoGen 
             << " GeV, ΔM=" << deltaMass_recoGen << " GeV" << endl;
      } else {
        cout << "Gen Z: NOT FOUND" << endl;

        // Debug: show what gen leptons are available
        if (analysisType == AnalysisType::ZJET_MUMU && nGen > 0) {
          cout << "Gen-level particles in muonTree (nGen=" << nGen << "):" << endl;
          for (int ig = 0; ig < nGen && ig < (int)genPID->size(); ig++) {
            cout << "  [" << ig << "] PID=" << (*genPID)[ig] 
                 << " status=" << (*genStatus)[ig]
                 << " pT=" << (*genPt)[ig] << " GeV";
            if (genMotherID && ig < (int)genMotherID->size()) {
              cout << " momPID=" << (*genMotherID)[ig];
            }
            cout << endl;
          }
        } else if (analysisType == AnalysisType::ZJET_EE && mcPID && mcPID->size() > 0) {
          cout << "Gen-level particles in ggHiNtuplizer (nMC=" << mcPID->size() << "):" << endl;
          for (size_t ig = 0; ig < mcPID->size() && ig < 10; ig++) {
            cout << "  [" << ig << "] PID=" << (*mcPID)[ig]
                 << " status=" << (*mcStatus)[ig]
                 << " pT=" << (*mcPt)[ig] << " GeV";
            if (mcMomPID && ig < mcMomPID->size()) {
              cout << " momPID=" << (*mcMomPID)[ig];
            }
            cout << endl;
          }
        }
      }
    }


    // Require at least 1 jet
    if (nref < 1) continue;
    if (jtpt[0] < jtptmin) continue;

    nEvents_hasJets++;
    if (debug) {
      log(LOG_TRACE, "Event " + std::to_string(i) + ": passed leading-jet preselection");
    }

    eh->event_vz->Fill(vz, evtwt);

    // Z+jet variables
    double Z_pt = selectedZ.pt;
    double Z_eta = selectedZ.eta;
    double Z_phi = selectedZ.phi;
    double Z_mass = selectedZ.mass;
    double Z_rapidity = selectedZ.rapidity;
    double jet_pt, jet_eta, jet_phi, ptavgtp, alpha, balance;
    double asymmtp;
    //double djrespasymm;

    // Z pt selection
    if (Z_pt < 60) continue;

    nEvents_ZptCut++;

    if (debug && nEvents_ZptCut <= 5) {
      cout << "\n=== Z+Jet Event " << nEvents_ZptCut << " (event index " << i << ") ===" << endl;
      cout << "Z: pT=" << Z_pt << " GeV, eta=" << Z_eta << ", phi=" << Z_phi << ", mass=" << Z_mass << " GeV" << endl;
      cout << "nref=" << nref << " jets (before selection)" << endl;
    }

    if (checkvalidjet) { // This is based on validity of JEC. - obsolete?
      for (int j = 0; j < nref; ++j) {
        if (abs(jteta[j]) > 2.964)
          jtpt[j] = 0; // Always invalid jets
        else if (abs(jteta[j]) > 2.5 && jtpt[j] > 120)
          jtpt[j] = 0; // ?
        else if (abs(jteta[j]) > 1.93 && jtpt[j] > 170)
          jtpt[j] = 0; // ?
        if (jtpt[j] < 80)
          jtpt[j] = 0;
        // if 80-120 abseta < 2.964
        // 120-170 < 2.5
        // 170-1000 < 1.93
      }
    }

    // JET ID - this is 2023 AK4CHS jet selection 12/2024
    // fill passjteta for all jets in the event?
    std::vector<bool> passjetid(nref, false);
    for (int j = 0; j < nref; ++j) {
      if (checkjetid) {
        float eta = fabs(jteta[j]);
        
        // 1. Calculate Multiplicities
        // Total Multiplicity
        int m = jtchm[j] + jtPfNHM[j] + jtPfCEM[j] + jtPfNEM[j] + jtPfMUM[j];
        // Charged Multiplicity
        int cm = jtchm[j] + jtPfCEM[j] + jtPfMUM[j];
        // Neutral Multiplicity
        int nm = jtPfNHM[j] + jtPfNEM[j];

        // 2. Apply pp-specific ID cuts based on eta range
        if (eta <= 2.6) {
          if ((jtnhf[j] < 0.9) && (jtnef[j] < 0.9) && (m > 1) && (jtmuf[j] < 0.8) && 
              (jtchf[j] > 0.01) && (cm > 0) && (jtcef[j] < 0.8)) {
            passjetid[j] = true;
          }
        } 
        else if (eta > 2.6 && eta <= 2.7) {
          if ((jtnhf[j] < 0.9) && (jtnef[j] < 0.99) && (jtmuf[j] < 0.8) && 
              (cm > 0) && (jtcef[j] < 0.8)) {
            passjetid[j] = true;
          }
        } 
        else if (eta > 2.7 && eta <= 3.0) {
          if ((jtnef[j] < 0.99) && (nm > 1)) {
            passjetid[j] = true;
          }
        } 
        else if (eta > 3.0 && eta <= 5.0) {
          if ((jtnhf[j] > 0.2) && (jtnef[j] < 0.9) && (nm > 10)) {
            passjetid[j] = true;
          }
        }
      } else {
        passjetid[j] = true; // If checkjetid is false, everyone passes
        //       if (jtpt[j] > jtptmin)    cout << passjetid[j] << endl;
        //  if (passjetid[j] < 2 and nref > 2  and jtpt[j] > 70  and jtpt[1] >
        //  40) cout << "Pass jetid: " << passjetid[j] << " pt: " << jtpt[j] <<
        //  " " << jteta[j] << " " << j << " " << i <<  endl;
      }
    }

    // Apply JEC
    if (debug) {
      log(LOG_TRACE, "Event " + std::to_string(i) + ": entering JEC loop");
    }
    for (int j = 0; j < nref; ++j) {
      jtpt_uncorr[j] = jtpt[j];

#if REDOJES == 1
      // cout << "Applying JES" << endl;
  JEC.SetJetPT(jtpt[j]);
  JEC.SetJetEta(jteta[j]);
  JEC.SetJetPhi(jtphi[j]);
  jtpt[j] = JEC.GetCorrectedPT();
#endif
      // JER not applied for photon+jet L3 residual analysis
    }
    if (debug) {
      log(LOG_TRACE, "Event " + std::to_string(i) + ": finished JEC loop");
    }

    // ========================================
    // Z+JET SELECTION
    // ========================================

    // Find leading and subleading away-side jets
    // identify all jets back-to-back with Z boson (dphi > 2.7)// No just apply a small dR requirement
    std::vector<int> awayJetIndices;
    awayJetIndices.reserve(nref);
    const float probeJetPtMin = jtptmin; // Minimum pT for the leading away-side jet (probe jet)
    const float alphaJetPtFloor = std::min(probeJetPtMin, jtptlimitforalpha);

    for (int j = 0; j < nref; j++) {
      if (checkjetid && passjetid[j] == 0) continue; // Apply jet ID
      if (jtpt[j] < alphaJetPtFloor) continue;

      // Calculate delta-phi with Z
      float dphi = abs(jtphi[j] - Z_phi);
      if (dphi > TMath::Pi())
        dphi = 2 * TMath::Pi() - dphi;

      // Back-to-back requirement
      if (dphi < 2.7) continue;

      // Check overlap with BOTH Z decay leptons
      float deta1 = jteta[j] - selectedZ.lepton1_eta;
      float dphi1 = abs(jtphi[j] - selectedZ.lepton1_phi);
      if (dphi1 > TMath::Pi()) dphi1 = 2 * TMath::Pi() - dphi1;
      float dR1 = sqrt(deta1 * deta1 + dphi1 * dphi1);

      float deta2 = jteta[j] - selectedZ.lepton2_eta;
      float dphi2 = abs(jtphi[j] - selectedZ.lepton2_phi);
      if (dphi2 > TMath::Pi()) dphi2 = 2 * TMath::Pi() - dphi2;
      float dR2 = sqrt(deta2 * deta2 + dphi2 * dphi2);

      if (dR1 < leptonJetDeltaR || dR2 < leptonJetDeltaR) continue;  // Reject if overlaps with either lepton

      // This jet is away-side
      awayJetIndices.push_back(j);
    }

    const int nAwayJets = awayJetIndices.size();
    if (debug) {
      log(LOG_DEBUG, "Event " + std::to_string(i) + ": away-side jet candidates=" +
                         std::to_string(nAwayJets));
    }

    if (nAwayJets < 1) continue; // Need at least one away-side jet

    // The leading away-side jet is still required to satisfy the nominal
    // probe-jet threshold.
    if (jtpt[awayJetIndices[0]] < probeJetPtMin) continue;

    nEvents_hasAwayJet++;

    // Debug: Print jet selection details
    if (debug && nEvents_hasAwayJet <= 5) {
      cout << "  Found " << nAwayJets << " away-side jets" << endl;
      for (int ii = 0; ii < nAwayJets && ii < 3; ii++) {
        int jidx = awayJetIndices[ii];
        cout << "    Jet " << ii << ": pT=" << jtpt[jidx] << " GeV, eta=" << jteta[jidx] 
             << ", phi=" << jtphi[jidx];

        float dphi_debug = abs(jtphi[jidx] - Z_phi);
        if (dphi_debug > TMath::Pi()) dphi_debug = 2 * TMath::Pi() - dphi_debug;
        cout << ", Δφ(Z,jet)=" << dphi_debug << endl;
      }
      cout << "  Leading jet: pT=" << jtpt[awayJetIndices[0]] << " GeV" << endl;
      cout << "  Alpha=" << alpha << ", Balance=" << balance << endl;
    }


    // Sort away-side jets by pT (descending)
    for (int ii = 0; ii < nAwayJets - 1; ii++) {
      for (int jj = ii + 1; jj < nAwayJets; jj++) {
        if (jtpt[awayJetIndices[jj]] > jtpt[awayJetIndices[ii]]) {
          int temp = awayJetIndices[ii];
          awayJetIndices[ii] = awayJetIndices[jj];
          awayJetIndices[jj] = temp;
        }
      }
    }

    // Leading away-side jet (probe jet)
    int awayJetIdx = awayJetIndices[0];

    // Calculate alpha (using 2nd away-side jet if available)
    alpha = 0;
    if (nAwayJets >= 2) {
      int secondAwayJetIdx = awayJetIndices[1];
      if (jtpt[secondAwayJetIdx] >= jtptlimitforalpha) {
        float ptavg_temp = Z_pt;
        alpha = jtpt[secondAwayJetIdx] / ptavg_temp;
      }
    } else {
      alpha = 0; // Only one away-side jet
    }

    // Calculate Z+jet variables
    jet_pt = jtpt[awayJetIdx];
    jet_eta = jteta[awayJetIdx];
    jet_phi = jtphi[awayJetIdx];

    float dphi_zjet = abs(jet_phi - Z_phi);
    if (dphi_zjet > TMath::Pi())
      dphi_zjet = 2 * TMath::Pi() - dphi_zjet;

    ptavgtp = Z_pt;
    balance = jet_pt / Z_pt;  // Response
    asymmtp = balance;            // For compatibility with histogram filling

    // Apply jet veto map to Z and leading/subleading away-side jets
    if (applyjetvetomap && vetomap) {
      bool passVetoMap = true;

      // Check leading away-side jet
      int jet_bin = vetomap->FindBin(jet_eta, jet_phi);
      if (vetomap->GetBinContent(jet_bin) > 0) passVetoMap = false;

      // Check subleading away-side jet if available
      if (passVetoMap && nAwayJets >= 2) {
        int secondAwayJetIdx = awayJetIndices[1];
        int subjet_bin = vetomap->FindBin(jteta[secondAwayJetIdx], jtphi[secondAwayJetIdx]);
        if (vetomap->GetBinContent(subjet_bin) > 0) passVetoMap = false;
      }

      if (!passVetoMap) continue;
    }

    nEvents_passVeto++;

    nEventsPassedSelection++;
    if (isMC && hasGenZ) {
      nEventsPassedWithGenZ++;
    }
    if (debug) {
      log(LOG_DEBUG, "Event " + std::to_string(i) + ": passed final selection, filling histograms");
    }
    // cout << "TP:" << tagpt << " " << probept << " " << alpha << endl;
    // ========================================
    // FILL Z+JET HISTOGRAMS
    // ========================================

    for (auto &histrange : _histos) {
      for (auto &h : histrange.second) {

        // Check if jet eta and centrality are in range for this histogram set
        if (jet_eta >= h->etamin && jet_eta < h->etamax &&
            hiBin >= h->hibinmin && hiBin < h->hibinmax) {

          // Z properties (using new Z+Jet specific histograms)
          h->z_pt->Fill(Z_pt, evtwt);
          h->z_eta->Fill(Z_eta, evtwt);
          h->z_phi->Fill(Z_phi, evtwt);
          h->z_mass->Fill(Z_mass, evtwt);
          // (Rapidity removed as we didn't add a specific histogram for it)

          // Gen-level Z histograms (MC only)
          if (isMC && hasGenZ) {
            if (h->genz_pt) h->genz_pt->Fill(genZ.pt, evtwt);
            if (h->genz_eta) h->genz_eta->Fill(genZ.eta, evtwt);
            if (h->genz_phi) h->genz_phi->Fill(genZ.phi, evtwt);

            // Z response: reco_pT / gen_pT
            if (h->zresponse && genZ.pt > 0) {
              h->zresponse->Fill(genZ.pt, Z_pt / genZ.pt, evtwt);
            }

            // Z pT resolution: (reco - gen) / gen
            if (h->z_ptres && genZ.pt > 0) {
              h->z_ptres->Fill((Z_pt - genZ.pt) / genZ.pt, evtwt);
            }
          }

          // Away-side jet properties
          h->awayside_jet_pt->Fill(jet_pt, evtwt);
          h->awayside_jet_eta->Fill(jet_eta, evtwt);
          h->awayside_jet_phi->Fill(jet_phi, evtwt);
          h->awayside_jet_uncorr_pt->Fill(jtpt_uncorr[awayJetIdx], evtwt);

          // Z+Jet system
          h->zjet_dphi->Fill(dphi_zjet, evtwt);
          h->zjet_balance->Fill(balance, evtwt);
          h->zjet_ptavg->Fill(ptavgtp, evtwt);
          h->zjet_alpha->Fill(alpha, evtwt);

          // Trigger histograms
          if (trigger) {
            if (h->HLT_Z) h->HLT_Z->Fill(1, evtwt);
            if (h->HLT_Z_ptav) h->HLT_Z_ptav->Fill(ptavgtp, evtwt);
          }

          // Alpha-dependent balance profiles  (analogous to dijet asymmetry)
          if (alpha < 0.1) {
            h->zjet_balance_a01->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a01->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.2) {
            h->zjet_balance_a02->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a02->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.3) {
            h->zjet_balance_a03->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a03->Fill(ptavgtp, jet_eta, balance, evtwt);

            // Jet composition for alpha < 0.3 (like dijets)
            h->jet_nef->Fill(jet_pt, jtnef[awayJetIdx], evtwt);
            h->jet_cef->Fill(jet_pt, jtcef[awayJetIdx], evtwt);
            h->jet_nhf->Fill(jet_pt, jtnhf[awayJetIdx], evtwt);
            h->jet_chf->Fill(jet_pt, jtchf[awayJetIdx], evtwt);
            h->jet_muf->Fill(jet_pt, jtmuf[awayJetIdx], evtwt);
          }
          if (alpha < 0.4) {
            h->zjet_balance_a04->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a04->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.5) {
            h->zjet_balance_a05->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a05->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.6) {
            h->zjet_balance_a06->Fill(ptavgtp, balance, evtwt);
            h->zjet_balance2D_a06->Fill(ptavgtp, jet_eta, balance, evtwt);
          }

          // 3D balance profiles (KEY for L3 residual derivation)
          // Fill with CUMULATIVE alpha cuts (matching dijet analyse.cc pattern)
          // Alpha bins are read from histograms::alphavalues array
          // Each event with alpha < threshold is filled into the bin corresponding to threshold
          // Only fill in the wide eta bin since these have internal eta binning
          if ((h->etamin - h->etamax) < -10) {
            // Loop over alpha thresholds from histograms::alphavalues (skip first bin which is 0)
            for (unsigned int ia = 1; ia <= histograms::nalphavalues; ++ia) {
              double alphaThreshold = histograms::alphavalues[ia];
              double alphaFillValue = alphaThreshold - 0.0001;

              if (alpha < alphaThreshold) {
                // Fill balance profiles (weighted)
                h->zjet_balance3D->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->zjet_balance3Dwide->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->zjet_balance3Dnarrow->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->zjet_balance3Dabseta->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->zjet_balance3Dabsetawide->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->zjet_balance3Dabsetanarrow->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);

                // Fill counts histograms (UNWEIGHTED - just count entries)
                if (h->zjet_balance3D_counts) h->zjet_balance3D_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->zjet_balance3Dwide_counts) h->zjet_balance3Dwide_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->zjet_balance3Dnarrow_counts) h->zjet_balance3Dnarrow_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->zjet_balance3Dabseta_counts) h->zjet_balance3Dabseta_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->zjet_balance3Dabsetawide_counts) h->zjet_balance3Dabsetawide_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->zjet_balance3Dabsetanarrow_counts) h->zjet_balance3Dabsetanarrow_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
              }
            }
          }
        }
      }
    }

    // Additional jet histograms for all jets in the event
    for (int j = 0; j < nref; ++j) {
      for (auto &histrange : _histos) {
        for (auto &h : histrange.second) {
          if (jteta[j] >= h->etamin && jteta[j] < h->etamax &&
              hiBin >= h->hibinmin && hiBin < h->hibinmax) {

            if (checkjetid && passjetid[j] == 0) continue;

            h->jetetaphi->Fill(jteta[j], jtphi[j], weight);

            if (j == 0 && passjetid[0]) {
              // trigger check; leading jet pt
              if (trigger)
                if (h->HLT_Z) h->HLT_Z->Fill(jtpt[awayJetIndices[0]], evtwt);
            }

            // These are actually obsolete after all the selections
            /* if (j == 0 and nref > 1 and dphitp > 2.7) { // Fill dijet
              system based on leading jet pT
              h->dijetasymmetry->Fill(abs(djetasymm),evtwt);
              h->dijetasymmetry_now->Fill(abs(djetasymm));
              h->dijetdeltaphi->Fill(dphi,evtwt);
              h->dijetdeltaeta->Fill(ddeta,evtwt);
              } */

            h->jet_pt->Fill(jtpt[j], evtwt);
            h->jet_pt_now->Fill(jtpt[j], 1);
            h->jet_uncorr_pt->Fill(jtpt_uncorr[j], evtwt);
            h->jet_pt_genweight->Fill(jtpt[j], weight);
            h->jet_eta->Fill(jteta[j], evtwt);
            h->jet_phi->Fill(jtphi[j], evtwt);

            if (isMC) {
              h->genjet_pt->Fill(jtpt_gen[j], evtwt);
              h->genjet_eta->Fill(jteta_gen[j], evtwt);
              h->genjet_phi->Fill(jtphi_gen[j], evtwt);

              h->jetresponse->Fill(jtpt_gen[j], jtpt[j] / jtpt_gen[j], evtwt);
              h->ptres->Fill((jtpt[j] - jtpt_gen[j]) / jtpt_gen[j], evtwt);
              h->ptgenvsptreco->Fill(jtpt_gen[j], jtpt[j], evtwt);
              h->ptrecovsweight->Fill(jtpt[j], weight);
              h->ptgenvsweight->Fill(jtpt_gen[j], weight);

              if (h->responses3D) {
                h->responses3D->Fill(jtpt_gen[j], jteta_gen[j],
                                     jtpt[j] / jtpt_gen[j], evtwt);
              }
              if (h->phiresponse) {
                h->phiresponse->Fill(jtpt_gen[j], jteta_gen[j],
                                     jtphi[j] - jtphi_gen[j], evtwt);
              }
              if (h->etaresponse) {
                h->etaresponse->Fill(jtpt_gen[j], jteta_gen[j],
                                     jteta[j] - jteta_gen[j], evtwt);
              }
            }
          }
        }
      }
    }
  } // end of event loop

  log(LOG_INFO, "========================================");
  log(LOG_INFO, "Z+Jet Analysis Summary");
  log(LOG_INFO, "========================================");
  log(LOG_INFO, "Total events processed: " +
                    std::to_string((i_processed < nentries) ? i_processed : nentries));

  log(LOG_INFO, "--- Z Reconstruction ---");
  if (useMuonFlavor) {
    log(LOG_INFO, "Events with muon Z: " + std::to_string(nEventsWithMuonZ) +
                      " (" + formatPercent(nEventsWithMuonZ, nentries) + ")");
  }
  if (useElectronFlavor) {
    log(LOG_INFO, "Events with electron Z: " +
                      std::to_string(nEventsWithElectronZ) + " (" +
                      formatPercent(nEventsWithElectronZ, nentries) + ")");
  }
  
  if (isMC) {
    log(LOG_INFO, "Events with gen-level Z: " +
                      std::to_string(nEventsWithGenZ) + " (" +
                      formatPercent(nEventsWithGenZ, nentries) + ")");
  }
  
  log(LOG_INFO, "--- Selection Flow ---");
  Long64_t nWithZ = (analysisType == AnalysisType::ZJET_MUMU) ? nEventsWithMuonZ : 
                    (analysisType == AnalysisType::ZJET_EE) ? nEventsWithElectronZ :
                    (nEventsWithMuonZ + nEventsWithElectronZ);
  
  log(LOG_INFO, "Events with Z:              " + std::to_string(nWithZ) +
                    " (" + (nWithZ > 0 ? std::string("100%") : std::string("n/a")) + ")");
  log(LOG_INFO, "  + has jets:               " + std::to_string(nEvents_hasJets) +
                    " (" + formatPercent(nEvents_hasJets, nWithZ) + ")");
  log(LOG_INFO, "  + Z_pT > 60 GeV:          " + std::to_string(nEvents_ZptCut) +
                    " (" + formatPercent(nEvents_ZptCut, nWithZ) + ")");
  log(LOG_INFO, "  + has away-side jet:      " +
                    std::to_string(nEvents_hasAwayJet) + " (" +
                    formatPercent(nEvents_hasAwayJet, nWithZ) + ")");
  log(LOG_INFO, "  + passes veto map:        " +
                    std::to_string(nEventsPassedSelection) + " (" +
                    formatPercent(nEventsPassedSelection, nWithZ) + ")");
  
  if (isMC) {
    log(LOG_INFO, "--- Gen-Matching (MC only) ---");
    log(LOG_INFO, "Events passed with gen-Z: " +
                      std::to_string(nEventsPassedWithGenZ));
    if (nEventsPassedSelection > 0) {
      float genMatchEff = 100.0 * nEventsPassedWithGenZ / nEventsPassedSelection;
      std::ostringstream genMatchLine;
      genMatchLine << std::fixed << std::setprecision(2)
                   << "Gen-matching efficiency: " << genMatchEff << "%";
      log(LOG_INFO, genMatchLine.str());
    }
  }
  
  log(LOG_INFO, "========================================");
  log(LOG_INFO, "FINAL: " + std::to_string(nEventsPassedSelection) +
                    " events passed all selections");
  log(LOG_INFO, "========================================");

  // Write output histograms
  for (auto &histrange : _histos) {
    for (auto &h : histrange.second) {
      h->Write();
    }
  }
  eh->Write();

  // Write and close output file
  outfile->Write();
  log(LOG_INFO, "Wrote " + outputfilename);

  if (mapfile) {
    mapfile->Close();
    delete mapfile;
    mapfile = nullptr;
  }

  // Close file properly - TFile destructor handles all owned histogram cleanup
  outfile->Close();
  delete outfile;
}
