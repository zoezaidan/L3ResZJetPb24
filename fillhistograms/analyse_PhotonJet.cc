#include <iostream>
using std::cout;
using std::endl;

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TTreeCache.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <typeinfo>
#include <vector>

#include "configurations.h"
#include "histograms.h"
#include "settings.h"

#include "eventhistograms.h"
#include "helpers.h"
#include "input_config.h"
#include "chain_builder.h"
#include "JetCorrector.h"

R__LOAD_LIBRARY(histograms_C.so)
R__LOAD_LIBRARY(eventhistograms_C.so)

std::mt19937 _mersennetwister;
std::uint32_t _seed = 4;
//_seed = 4;

bool debug = false;
bool applyjetvetomap = true;

// Helper function to load pthat weights from file
std::map<float, double> LoadPthatWeights(const std::string& weightsFile) {
    std::map<float, double> weights;
    std::ifstream fin(weightsFile);
    if (!fin) {
    log(LOG_ERROR, "Could not open pthat weights file: " + weightsFile);
        return weights;
    }
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        float bin; double w;
        if (iss >> bin >> w) {
            weights[bin] = w;
        }
    }
    return weights;
}

// Helper function to get the pthat bin for a given value
float GetPthatBin(float pthat, const std::vector<float>& bins) {
    float result = bins.front();
    for (size_t i = 0; i < bins.size(); ++i) {
        if (pthat >= bins[i]) result = bins[i];
        else break;
    }
    return result;
}

// Photon+Jet analysis for L3 residual corrections
// jetTree: jet tree path, e.g. "ak4PFJetAnalyzer/t" or "ak4PFJetAnalyzerSDZcut1/t"
void analyse_PhotonJet(string input = "PHOTONHP",
                       string outputfiletag = "AK4_photonjet",
                       bool isMC = false, bool checkjetid = false,
                       string inputType = "era", int maxFiles = -1,
                       int maxEvents = -1, string outputDir = "",
                       int batchIndex = -1, int totalBatches = 1,
                       string jetPath = "ak4PFJetAnalyzer/t") {

  if (debug && g_verbosity < LOG_TRACE) {
    g_verbosity = LOG_TRACE;
  }
  log(LOG_INFO, "Initialising Photon+Jet analysis");

  bool usecalotrig = false;
  bool checkvalidjet =
      false; // this is for checking valid jet range after applying l2. now for
             // tightly limited range. TODO: do something smarter

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
    config.outputDir =
        "/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L3ResPhotonJet";
  } else {
    config.outputDir = outputDir;
  }
  log(LOG_DEBUG, "Configured analysis inputs and output directory");

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

  if (batchIndex >= 0) {
    // For batch mode, use just the outputfiletag (cleaner naming)
    outputfilename = Form("%s/%s_batch%d_of_%d.root",
                         config.outputDir.c_str(),
                         outputfiletag.c_str(), batchIndex, totalBatches);
  } else {
    outputfilename = Form("%s/%s_%s.root", config.outputDir.c_str(),
                         inputName.c_str(), outputfiletag.c_str());
  }
  if (debug)
    outputfilename = "test.root";

  log(LOG_INFO, "Output file: " + outputfilename);

  // Define tree paths
  std::string evtPath = "hiEvtAnalyzer/HiTree";
  std::string triggerPath = "hltanalysis/HltTree";
  std::string skimPath = "skimanalysis/HltTree";
  std::string photonPath = "ggHiNtuplizer/EventTree";

  log(LOG_INFO, "Using jet tree: " + jetPath);
  log(LOG_INFO, "Building input chains");
  TreeChains *chains = BuildChainsFromConfig(config, jetPath, true);
  if (!chains || chains->nEntries == 0) {
    log(LOG_ERROR, "No entries found in input");
    return;
  }
  log(LOG_INFO,
      "Input chain building complete with " +
          std::to_string(chains->nEntries) + " entries");

  auto evtTree = chains->evtChain;
  auto photonTree = chains->photonChain;
  auto jetTree = chains->jetChain;
  auto triggerTree = chains->triggerChain;
  auto skimTree = chains->skimChain;

  // Cuts and weights from event tree
  Int_t hiBin = -1;
  Float_t weight = 1, vz = 0, pthat = 0, evtwt = 1;

  // Photon variables (vectors)
  Int_t nPho;
  std::vector<float> *phoEt = 0;
  std::vector<float> *phoEta = 0;
  std::vector<float> *phoPhi = 0;
  std::vector<float> *phoE = 0;
  std::vector<float> *phoSCEta = 0;
  std::vector<float> *phoSCPhi = 0;
  std::vector<float> *phoHoverE = 0;
  std::vector<float> *phoSigmaIEtaIEta = 0;
  std::vector<float> *phoR9 = 0;
  std::vector<float> *pfcIso3subUEec = 0;  // PF charged hadron isolation
  std::vector<float> *pfnIso3subUEec = 0;  // PF neutral hadron isolation
  std::vector<float> *pfpIso3subUEec = 0;  // PF photon isolation

  // MC truth matching branches (only used for MC)
  std::vector<int> *mcPID = 0;
  std::vector<int> *mcMomPID = 0;
  std::vector<float> *mcCalIsoDR04 = 0;
  std::vector<float> *mcPt = 0;
  std::vector<float> *mcEta = 0;
  std::vector<float> *mcPhi = 0;
  std::vector<int> *pho_genMatchedIndex = 0;

  

  // Now enable only the branches we need
  evtTree->SetBranchStatus("*", 0);
  evtTree->SetBranchStatus("hiBin", 1);
  evtTree->SetBranchStatus("vz", 1);
  if (isMC) {
    evtTree->SetBranchStatus("weight", 1);
    evtTree->SetBranchStatus("pthat", 1);
  }
  // Set branch addresses BEFORE SetBranchStatus (like analyse.cc does)
  evtTree->SetBranchAddress("hiBin", &hiBin);
  evtTree->SetBranchAddress("vz", &vz);
  if (isMC) {
    evtTree->SetBranchAddress("weight", &weight);
    evtTree->SetBranchAddress("pthat", &pthat);
  }

  //// EVENT FILTERS
  // auto skimTree = (TTree*)inFile->Get(skimPath.c_str());
  // if (!isMC) skimTree->SetBranchStatus("*",1);

  // Int_t pprimaryVertexFilter = 1;
  // if (!isMC) skimTree->SetBranchAddress("pprimaryVertexFilter",
  // &pprimaryVertexFilter);

  Int_t trigger = 0;

  // Photon trigger
  Int_t HLT_Photon30 = 1;

  // auto triggerTree = (TTree*)inFile->Get(triggerPath.c_str());

  if (!isMC) {
  //   cout << "Use Photon trigger: HLT_PPRefGEDPhoton30_v6" << endl;
    triggerTree->SetBranchStatus("*", 0);
    triggerTree->SetBranchStatus("HLT_PPRefGEDPhoton30_v6", 1);
    triggerTree->SetBranchAddress("HLT_PPRefGEDPhoton30_v6", &HLT_Photon30);
  }
  log(LOG_DEBUG, "Configured event, trigger, jet, and photon branches");

  // Simple cache setup for remote xrootd reads.
  TTreeCache::SetLearnEntries(10);
  evtTree->AddBranchToCache("*", true);
  triggerTree->AddBranchToCache("*", true);
  photonTree->AddBranchToCache("*", true);
  jetTree->AddBranchToCache("*", true);

  // JETS
  // Disable all branches first, then enable only what we need
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

  Int_t jtn[MAXJETS];

  jetTree->SetBranchAddress("evt", &evt);
  jetTree->SetBranchAddress("nref", &nref);
  jetTree->SetBranchAddress("rawpt", &jtpt); // we want uncorrected rawpt, jtpt might have some JEC already applied
  jetTree->SetBranchAddress("jteta", &jteta);
  jetTree->SetBranchAddress("jtphi", &jtphi);
  
  jetTree->SetBranchStatus("evt", 1);
  jetTree->SetBranchStatus("nref", 1);
  jetTree->SetBranchStatus("rawpt", 1);
  jetTree->SetBranchStatus("jteta", 1);
  jetTree->SetBranchStatus("jtphi", 1);

  jetTree->SetBranchAddress("jtPfNHF", &jtnhf);
  jetTree->SetBranchAddress("jtPfCHF", &jtchf);
  jetTree->SetBranchAddress("jtPfNEF", &jtnef);
  jetTree->SetBranchAddress("jtPfCEF", &jtcef);
  jetTree->SetBranchAddress("jtPfMUF", &jtmuf);
  jetTree->SetBranchAddress("jtPfCHM", &jtchm);
  jetTree->SetBranchStatus("jtPfNHF", 1);
  jetTree->SetBranchStatus("jtPfCHF", 1);
  jetTree->SetBranchStatus("jtPfNEF", 1);
  jetTree->SetBranchStatus("jtPfCEF", 1);
  jetTree->SetBranchStatus("jtPfMUF", 1);
  jetTree->SetBranchStatus("jtPfCHM", 1);

  // Gen level jet information
  Float_t jtpt_gen[MAXJETS];
  Float_t jteta_gen[MAXJETS];
  Float_t jtphi_gen[MAXJETS];
  Float_t refdrjt[MAXJETS];

  if (isMC) {
    jetTree->SetBranchAddress("refpt", &jtpt_gen);
    jetTree->SetBranchAddress("refeta", &jteta_gen);
    jetTree->SetBranchAddress("refphi", &jtphi_gen);
    jetTree->SetBranchAddress("refdrjt", &refdrjt);
    jetTree->SetBranchStatus("refpt", 1);
    jetTree->SetBranchStatus("refeta", 1);
    jetTree->SetBranchStatus("refphi", 1);
    jetTree->SetBranchStatus("refdrjt", 1);
  }

  // Set photon branch addresses (vectors)
  photonTree->SetBranchStatus("*", 0);
  photonTree->SetBranchStatus("nPho", 1);
  photonTree->SetBranchStatus("phoEt", 1);
  photonTree->SetBranchStatus("phoEta", 1);
  photonTree->SetBranchStatus("phoPhi", 1);
  photonTree->SetBranchStatus("phoE", 1);
  photonTree->SetBranchStatus("phoSCEta", 1);
  photonTree->SetBranchStatus("phoSCPhi", 1);
  photonTree->SetBranchStatus("phoHoverE", 1);
  photonTree->SetBranchStatus("phoSigmaIEtaIEta_2012", 1);
  photonTree->SetBranchStatus("phoR9_2012", 1);
  photonTree->SetBranchStatus("pfcIso3subUEec", 1);
  photonTree->SetBranchStatus("pfnIso3subUEec", 1);
  photonTree->SetBranchStatus("pfpIso3subUEec", 1);
  if (isMC) {
    photonTree->SetBranchStatus("mcPID", 1);
    photonTree->SetBranchStatus("mcMomPID", 1);
    photonTree->SetBranchStatus("mcCalIsoDR04", 1);
    photonTree->SetBranchStatus("mcPt", 1);
    photonTree->SetBranchStatus("mcEta", 1);
    photonTree->SetBranchStatus("mcPhi", 1);
    photonTree->SetBranchStatus("pho_genMatchedIndex", 1);
  }
  photonTree->SetBranchAddress("nPho", &nPho);
  photonTree->SetBranchAddress("phoEt", &phoEt);
  photonTree->SetBranchAddress("phoEta", &phoEta);
  photonTree->SetBranchAddress("phoPhi", &phoPhi);
  photonTree->SetBranchAddress("phoE", &phoE);
  photonTree->SetBranchAddress("phoSCEta", &phoSCEta);
  photonTree->SetBranchAddress("phoSCPhi", &phoSCPhi);
  photonTree->SetBranchAddress("phoHoverE", &phoHoverE);
  photonTree->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta);
  photonTree->SetBranchAddress("phoR9_2012", &phoR9);
  photonTree->SetBranchAddress("pfcIso3subUEec", &pfcIso3subUEec);
  photonTree->SetBranchAddress("pfnIso3subUEec", &pfnIso3subUEec);
  photonTree->SetBranchAddress("pfpIso3subUEec", &pfpIso3subUEec);
  if (isMC) {
    photonTree->SetBranchAddress("mcPID", &mcPID);
    photonTree->SetBranchAddress("mcMomPID", &mcMomPID);
    photonTree->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04);
    photonTree->SetBranchAddress("mcPt", &mcPt);
    photonTree->SetBranchAddress("mcEta", &mcEta);
    photonTree->SetBranchAddress("mcPhi", &mcPhi);
    photonTree->SetBranchAddress("pho_genMatchedIndex", &pho_genMatchedIndex);
  }

  //TODO: Add electron veto for photons

  auto pthatWeights = LoadPthatWeights("jecfiles/2024_PP_30_170_pthat_weights.txt");
  std::vector<float> pthatBins;
  for (const auto& kv : pthatWeights) pthatBins.push_back(kv.first);
  std::sort(pthatBins.begin(), pthatBins.end());
  if (isMC && pthatWeights.empty()) {
    log(LOG_WARNING, "No pthat weights were loaded; MC event weights may be zero");
  } else if (isMC) {
    log(LOG_DEBUG,
      "Loaded " + std::to_string(pthatWeights.size()) + " pthat bins");
  }

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

          // Use PHOTONJET analysis type - only creates photon+jet specific histograms
          histograms *h = new histograms(dir2, etaedges[i], etaedges[i + 1],
                                         hibins[j], hibins[j + 1], isMC, AnalysisType::PHOTONJET);
          _histos[name2.c_str()].push_back(h);
        }
      }
    }
  }
  log(LOG_INFO, "Created centrality and eta histogram directories");

  outfile->mkdir("event");
  TDirectory *dir = outfile->GetDirectory("event");
  assert(dir);
  dir->cd();

  eventhistograms *eh = new eventhistograms(dir, isMC);

#if REDOJES == 1
  log(LOG_INFO, "Applying MC JEC from file " + jecfile);
  vector<string> JECFiles;
  // This is MCTruth (L2Relative)  
  JECFiles.push_back(jecfile.c_str());
  // L2 residual for data only
  if (!isMC) {
    log(LOG_INFO, "Applying L2 Residual from file " + l2file);
    JECFiles.push_back(l2file.c_str());
  }
   
  JetCorrector JEC(JECFiles);
  log(LOG_INFO, "Initialised jet corrector inputs");
#endif

  // JER not needed for photon+jet L3 residual analysis

  // Jet veto map
  auto mapfile = new TFile("jecfiles/Summer24Prompt24_RunBCDEFGHI.root","READ"); 
  auto vetomap = (TH2D*)mapfile->Get("jetvetomap_all");
  
  // Null check for veto map
  if (applyjetvetomap && (!mapfile || mapfile->IsZombie() || !vetomap)) {
    log(LOG_WARNING, "Veto map unavailable, disabling veto map selection");
    applyjetvetomap = false;
  } else if (applyjetvetomap) {
    log(LOG_DEBUG, "Loaded jet veto map successfully");
  }

  log(LOG_INFO, "Number of entries: " + std::to_string(chains->nEntries));
  Long64_t nentries = chains->nEntries;
  if (config.maxEvents > 0 && config.maxEvents < nentries) {
    nentries = config.maxEvents;
  }
  if (debug && nentries > 1000)
    nentries = 1000;

  log(LOG_INFO, "Processing " + std::to_string(nentries) + " events");
  for (Long64_t i = 0; i < nentries; ++i) {
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + "/" +
                         std::to_string(nentries) + ": reading event tree");

    evtTree->GetEntry(i);
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + "/" +
                         std::to_string(nentries) + ": reading trigger tree");

    triggerTree->GetEntry(i);
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + "/" +
                         std::to_string(nentries) + ": reading photon tree");

    photonTree->GetEntry(i);
    log_progress_every(i + 1, nentries, 1000, LOG_INFO);

    // Photon trigger logic
     if (!isMC) {
       trigger = HLT_Photon30;
     }
     if (isMC) trigger = true; // MC: no trigger requirement

     if (!trigger) continue;
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": passed trigger selection");


    auto get_weight = [pthatWeights, pthatBins](float pthat) {
        float bin = GetPthatBin(pthat, pthatBins);
        auto it = pthatWeights.find(bin);
        if (it != pthatWeights.end()) return static_cast<float>(it->second);
        return 0.f;
    };

    evtwt = 1;
    if (isMC) {
      evtwt *= weight*get_weight(pthat);
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) +
                         ": computed event weight = " + std::to_string(evtwt));


    // cout << weight << " " << evtwt << endl;

    // BASIC EVENT FILTERS
    //  if (!isMC) {
    //    skimTree->GetEntry(i);
    //    if (pprimaryVertexFilter != 1) continue;
    //  }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": reading jet tree");

    jetTree->GetEntry(i);

    if (debug) {
      log(LOG_DEBUG, "Processing event " + std::to_string(i) +
                         ", evt number: " + std::to_string(evt));
    }

    // Photon+jet: need at least 1 photon and 1 jet
    if (nPho < 1)
      continue;
    if (nref < 1)
      continue;
    if (jtpt[0] < jtptmin) {
      continue;
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) +
                         ": passed basic photon and jet multiplicity selection");


    eh->event_vz->Fill(vz, evtwt);
    if (isMC)
      eh->event_pthatwsgenweight->Fill(pthat, weight);

    // This is photon+jet analysis
    double photon_pt, photon_eta, photon_phi, jet_pt, jet_eta, jet_phi, ptavgtp,
        alpha, balance;
    double asymmtp;
    double djrespasymm;

    if (checkvalidjet) { // This is based on validity of JEC. - obsolete?
      for (int j = 0; j < nref; ++j) {
        if (abs(jteta[j]) > 2.964)
          jtpt[j] = 0; // Always invalid jets

        else if (abs(jteta[j]) > 2.5 and jtpt[j] > 120)
          jtpt[j] = 0; // ?
        else if (abs(jteta[j]) > 1.93 and jtpt[j] > 170)
          jtpt[j] = 0; // ?

        if (jtpt[j] < 80)
          jtpt[j] = 0;
        // if 80-120 abseta < 2.964
        // 120-170 < 2.5
        // 170-1000 < 1.93
      }
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": evaluated valid jet range");


    // JET ID - this is 2023 AK4CHS jet selection 12/2024
    // fill passjteta for all jets in the event?
    bool passjetid[nref];
    for (int j = 0; j < nref; ++j) {
      passjetid[j] = true;
      if (checkjetid) {
        if (abs(jteta[j]) <= 2.6) {
          if (jtnhf[j] >= 0.99)
            passjetid[j] = false;
          if (jtnef[j] >= 0.9)
            passjetid[j] = false;
          if (jtchf[j] <= 0.01)
            passjetid[j] = false;
          if (jtcef[j] >= 0.8)
            passjetid[j] = false;
          if (jtmuf[j] >= 0.8)
            passjetid[j] = false;
          if (jtchm[j] <= 0)
            passjetid[j] = false;
        } else if (abs(jteta[j]) <= 2.7) {
          if (jtnhf[j] >= 0.9)
            passjetid[j] = false;
          if (jtnef[j] >= 0.99)
            passjetid[j] = false;
          if (jtmuf[j] >= 0.8)
            passjetid[j] = false;
          if (jtcef[j] >= 0.8)
            passjetid[j] = false;

        } else if (abs(jteta[j]) <= 3.0) {
          if (jtnhf[j] >= 0.99)
            passjetid[j] = false;
          if (jtnef[j] >= 0.99)
            passjetid[j] = false;
        } else if (abs(jteta[j]) <= 5.0) {
          if (jtnef[j] >= 0.4)
            passjetid[j] = false;
        }
        //	   if (jtpt[j] > jtptmin)	   cout << passjetid[j] << endl;
        //  if (passjetid[j] < 2 and nref > 2  and jtpt[j] > 70  and jtpt[1] >
        //  40) cout << "Pass jetid: " << passjetid[j] << " pt: " << jtpt[j] <<
        //  " " << jteta[j] << " " << j << " " << i <<  endl;
      }
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": finished jet ID selection");


    // Apply JEC
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
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": applied jet energy corrections");


    // ========================================
    // PHOTON+JET SELECTION
    // ========================================

    // 1. Find leading photon passing selection
    int leadPhotonIdx = -1;
    float leadPhotonPt = 0;
    int leadPhotonGenIdx = -1;

    for (int ipho = 0; ipho < nPho; ipho++) {
      int currentGenIdx = -1;
      // Kinematic cuts
      if ((*phoEt)[ipho] < 60.0)
        continue; // Trigger threshold
      if (abs((*phoEta)[ipho]) > 1.3)
        continue; // Barrel only
      
      if ((*phoSigmaIEtaIEta)[ipho] < 0.002)
        continue;

      // MC gen-matching selection: require matched generator photon with
      // reasonable mother PID and low calorimeter isolation. Only for MC.
      if (isMC) {
        if (!pho_genMatchedIndex) continue;
        if (ipho >= (int)pho_genMatchedIndex->size()) continue;
        int genIdx = (*pho_genMatchedIndex)[ipho];
        if (genIdx == -1) continue;
        if (!mcPID || genIdx >= (int)mcPID->size()) continue;
        if ((*mcPID)[genIdx] != 22) continue;
        if (!mcMomPID || genIdx >= (int)mcMomPID->size()) continue;
        int mom = (*mcMomPID)[genIdx];
        int absMom = std::abs(mom);
        if (!(absMom <= 22 || mom == -999)) continue;
        if (!mcCalIsoDR04 || genIdx >= (int)mcCalIsoDR04->size()) continue;
        if (!((*mcCalIsoDR04)[genIdx] < 3.0)) continue;
        currentGenIdx = genIdx;
      }

      // Find highest pT photon
      if ((*phoEt)[ipho] > leadPhotonPt) {
        leadPhotonPt = (*phoEt)[ipho];
        leadPhotonIdx = ipho;
        if (isMC) leadPhotonGenIdx = currentGenIdx;
      }
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": finished leading photon search");


    // Stop early if no photon passes the kinematic preselection
    if (leadPhotonIdx < 0)
      continue;

    // Photon ID cuts
    if ((*phoHoverE)[leadPhotonIdx] > 0.129991)
      continue;
    if ((*phoSigmaIEtaIEta)[leadPhotonIdx] > 0.0114521)
      continue;

    if ((*pfcIso3subUEec)[leadPhotonIdx] > 1.88518)
      continue;
    if ((*pfnIso3subUEec)[leadPhotonIdx] > 2.0)
      continue;
    if ((*pfpIso3subUEec)[leadPhotonIdx] > 2.0)
      continue;
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": passed photon ID cuts");


    // 2. Find leading and subleading away-side jets
    // First: identify all jets back-to-back with photon (dphi > 2.7)// No just apply a small dR requirement
    int awayJetIndices[MAXJETS];
    int nAwayJets = 0;

    for (int j = 0; j < nref; j++) {
      // Apply jet ID
      if (checkjetid && passjetid[j] == 0)
        continue;

      // Jet kinematic cuts
      if (jtpt[j] < 40.0)
        continue; // Minimum jet pT

      // Calculate delta-phi with photon
      float dphi = abs(jtphi[j] - (*phoPhi)[leadPhotonIdx]);
      if (dphi > TMath::Pi())
        dphi = 2 * TMath::Pi() - dphi;

      // Back-to-back requirement
      if (dphi < 2.7488935)
        continue; // 2*pi/3 = 2.0943951

      // Calculate delta-R (reject jets close to photon)
      float deta = jteta[j] - (*phoEta)[leadPhotonIdx];
      float deltaR = sqrt(deta * deta + dphi * dphi);
      if (deltaR < 0.4)
        continue; // Isolation cone

      // This jet is away-side
      awayJetIndices[nAwayJets] = j;
      nAwayJets++;
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": identified " +
                         std::to_string(nAwayJets) + " away-side jets");


    // Need at least one away-side jet
    if (nAwayJets < 1)
      continue;

    // Sort away-side jets by pT (descending)
    for (int i = 0; i < nAwayJets - 1; i++) {
      for (int j = i + 1; j < nAwayJets; j++) {
        if (jtpt[awayJetIndices[j]] > jtpt[awayJetIndices[i]]) {
          int temp = awayJetIndices[i];
          awayJetIndices[i] = awayJetIndices[j];
          awayJetIndices[j] = temp;
        }
      }
    }

    // Leading away-side jet (probe jet)
    int awayJetIdx = awayJetIndices[0];

    // 3. Calculate alpha (using 2nd away-side jet if available)
    alpha = 0;
    if (nAwayJets >= 2) {
      int secondAwayJetIdx = awayJetIndices[1];
      float ptavg_temp =(*phoEt)[leadPhotonIdx];
      alpha = jtpt[secondAwayJetIdx] / ptavg_temp;
    } else {
      alpha = 0; // Only one away-side jet
    }

    // 4. Calculate photon+jet variables
    photon_pt = (*phoEt)[leadPhotonIdx];
    photon_eta = (*phoEta)[leadPhotonIdx];
    photon_phi = (*phoPhi)[leadPhotonIdx];

    float genPhotonPt = -1;
    float genPhotonEta = 0;
    float genPhotonPhi = 0;
    bool hasGenPhoton = false;
    if (isMC && leadPhotonGenIdx >= 0 && mcPt && mcEta && mcPhi) {
      if (leadPhotonGenIdx < (int)mcPt->size() && leadPhotonGenIdx < (int)mcEta->size() && leadPhotonGenIdx < (int)mcPhi->size()) {
        genPhotonPt = (*mcPt)[leadPhotonGenIdx];
        genPhotonEta = (*mcEta)[leadPhotonGenIdx];
        genPhotonPhi = (*mcPhi)[leadPhotonGenIdx];
        hasGenPhoton = true;
      }
    }

    jet_pt = jtpt[awayJetIdx];
    jet_eta = jteta[awayJetIdx];
    jet_phi = jtphi[awayJetIdx];

    float dphi_photonjet = abs(jet_phi - photon_phi);
    if (dphi_photonjet > TMath::Pi())
      dphi_photonjet = 2 * TMath::Pi() - dphi_photonjet;

    ptavgtp = photon_pt;
    balance = jet_pt / photon_pt; // Response
    asymmtp = balance;            // For compatibility with histogram filling
          log(LOG_TRACE, "Event " + std::to_string(i + 1) +
                 ": computed photon-jet observables");


    // Apply jet veto map to photon and leading/subleading away-side jets
    if (applyjetvetomap && vetomap) {
      bool passVetoMap = true;
      
      // Check photon position - Not required for photon since the vetomap is mostly due to pixel failures
      // int pho_bin = vetomap->FindBin(photon_eta, photon_phi);
      // if (vetomap->GetBinContent(pho_bin) > 0) passVetoMap = false;
      
      // Check leading away-side jet (probe)
      if (passVetoMap) {
        int jet_bin = vetomap->FindBin(jet_eta, jet_phi);
        if (vetomap->GetBinContent(jet_bin) > 0) passVetoMap = false;
      }
      
      // Check subleading away-side jet if available
      if (passVetoMap && nAwayJets >= 2) {
        int secondAwayJetIdx = awayJetIndices[1];
        int subjet_bin = vetomap->FindBin(jteta[secondAwayJetIdx], jtphi[secondAwayJetIdx]);
        if (vetomap->GetBinContent(subjet_bin) > 0) passVetoMap = false;
      }
      
      if (!passVetoMap) continue;
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": passed veto map selection");


    // cout << "TP:" << tagpt << " " << probept << " " << alpha << endl;
    // ========================================
    // FILL PHOTON+JET HISTOGRAMS
    // ========================================

    for (auto &histrange : _histos) {
      for (auto &h : histrange.second) {

        // Check if jet eta and centrality are in range for this histogram set
        if (jet_eta >= h->etamin && jet_eta < h->etamax &&
            hiBin >= h->hibinmin && hiBin < h->hibinmax) {

          // Photon properties
          h->photon_pt->Fill(photon_pt, evtwt);
          h->photon_eta->Fill(photon_eta, evtwt);
          h->photon_phi->Fill(photon_phi, evtwt);
          h->photon_HoverE->Fill((*phoHoverE)[leadPhotonIdx], evtwt);
          h->photon_sigmaIetaIeta->Fill((*phoSigmaIEtaIEta)[leadPhotonIdx],
                                        evtwt);

          if (isMC && hasGenPhoton) {
            if (h->genphoton_pt) h->genphoton_pt->Fill(genPhotonPt, evtwt);
            if (h->genphoton_eta) h->genphoton_eta->Fill(genPhotonEta, evtwt);
            if (h->genphoton_phi) h->genphoton_phi->Fill(genPhotonPhi, evtwt);
            if (h->photonresponse && genPhotonPt > 0) h->photonresponse->Fill(genPhotonPt, photon_pt / genPhotonPt, evtwt);
            if (h->photon_ptres && genPhotonPt > 0) h->photon_ptres->Fill((photon_pt - genPhotonPt) / genPhotonPt, evtwt);
          }

          // Away-side jet properties
          h->awayside_jet_pt->Fill(jet_pt, evtwt);
          h->awayside_jet_eta->Fill(jet_eta, evtwt);
          h->awayside_jet_phi->Fill(jet_phi, evtwt);
          h->awayside_jet_uncorr_pt->Fill(jtpt_uncorr[awayJetIdx], evtwt);

          // Photon+Jet system
          h->photonjet_dphi->Fill(dphi_photonjet, evtwt);
          h->photonjet_balance->Fill(balance, evtwt);
          h->photonjet_ptavg->Fill(ptavgtp, evtwt);
          h->photonjet_alpha->Fill(alpha, evtwt);

          // Trigger histograms
          if (HLT_Photon30) {
            h->HLTPhoton30->Fill(1, evtwt);
            h->HLTPhoton30_ptav->Fill(ptavgtp, evtwt);
          }

          // Alpha-dependent balance profiles (analogous to dijet asymmetry)
          if (alpha < 0.1) {
            h->photonjet_balance_a01->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a01->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.2) {
            h->photonjet_balance_a02->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a02->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.3) {
            h->photonjet_balance_a03->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a03->Fill(ptavgtp, jet_eta, balance, evtwt);

            // Jet composition for alpha < 0.3 (like dijets)
            h->jet_nef->Fill(jet_pt, jtnef[awayJetIdx], evtwt);
            h->jet_cef->Fill(jet_pt, jtcef[awayJetIdx], evtwt);
            h->jet_nhf->Fill(jet_pt, jtnhf[awayJetIdx], evtwt);
            h->jet_chf->Fill(jet_pt, jtchf[awayJetIdx], evtwt);
            h->jet_muf->Fill(jet_pt, jtmuf[awayJetIdx], evtwt);
          }
          if (alpha < 0.4) {
            h->photonjet_balance_a04->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a04->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.5) {
            h->photonjet_balance_a05->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a05->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.6) {
            h->photonjet_balance_a06->Fill(ptavgtp, balance, evtwt);
            h->photonjet_balance2D_a06->Fill(ptavgtp, jet_eta, balance, evtwt);
          }

          // 3D balance profiles (KEY HISTOGRAMS for L3 residual derivation)
          // Fill with CUMULATIVE alpha cuts (matching dijet analyse.cc pattern)
          // Alpha bins are read from histograms::alphavalues array
          // Each event with alpha < threshold is filled into the bin corresponding to threshold
          // Only fill in the wide eta bin since these have internal eta binning
          if ((h->etamin - h->etamax) < -10) {
            // Loop over alpha thresholds from histograms::alphavalues (skip first bin which is 0)
            for (unsigned int ia = 1; ia <= histograms::nalphavalues; ++ia) {
              double alphaThreshold = histograms::alphavalues[ia];
              double alphaFillValue = alphaThreshold - 0.0001;  // Fill just below threshold to land in correct bin
              
              if (alpha < alphaThreshold) {
                // Fill balance profiles (weighted)
                h->photonjet_balance3D->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->photonjet_balance3Dwide->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->photonjet_balance3Dnarrow->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->photonjet_balance3Dabseta->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->photonjet_balance3Dabsetawide->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->photonjet_balance3Dabsetanarrow->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3D_jetpt) h->photonjet_balance3D_jetpt->Fill(jet_pt, jet_eta, alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3Dwide_jetpt) h->photonjet_balance3Dwide_jetpt->Fill(jet_pt, jet_eta, alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3Dnarrow_jetpt) h->photonjet_balance3Dnarrow_jetpt->Fill(jet_pt, jet_eta, alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3Dabseta_jetpt) h->photonjet_balance3Dabseta_jetpt->Fill(jet_pt, abs(jet_eta), alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3Dabsetawide_jetpt) h->photonjet_balance3Dabsetawide_jetpt->Fill(jet_pt, abs(jet_eta), alphaFillValue, balance, evtwt);
                if (h->photonjet_balance3Dabsetanarrow_jetpt) h->photonjet_balance3Dabsetanarrow_jetpt->Fill(jet_pt, abs(jet_eta), alphaFillValue, balance, evtwt);
                
                // Fill counts histograms (UNWEIGHTED - just count entries)
                if (h->photonjet_balance3D_counts) h->photonjet_balance3D_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->photonjet_balance3Dwide_counts) h->photonjet_balance3Dwide_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->photonjet_balance3Dnarrow_counts) h->photonjet_balance3Dnarrow_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->photonjet_balance3Dabseta_counts) h->photonjet_balance3Dabseta_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->photonjet_balance3Dabsetawide_counts) h->photonjet_balance3Dabsetawide_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->photonjet_balance3Dabsetanarrow_counts) h->photonjet_balance3Dabsetanarrow_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                
                // Fill balance distribution (photon_pT, alpha, balance_value)
                // Fill with weight for each cumulative alpha cut
                if (h->photonjet_balance_dist) h->photonjet_balance_dist->Fill(ptavgtp, alphaFillValue, balance, evtwt);
              }
            }
          }
        }
      }
    }
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": filled photon-jet histograms");


    // Additional jet histograms for all jets in the event
    for (int j = 0; j < nref; ++j) {

      for (auto &histrange : _histos) {
        for (auto &h : histrange.second) {

          if (jteta[j] >= h->etamin and jteta[j] < h->etamax and
              hiBin >= h->hibinmin and hiBin < h->hibinmax) {

            if (checkjetid and passjetid[j] == 0)
              continue;

            h->jetetaphi->Fill(jteta[j], jtphi[j], weight);

            if (j == 0 and passjetid[0]) {
              // Photon trigger check; leading jet pt
              if (HLT_Photon30)
                h->HLTPhoton30->Fill(jtpt[awayJetIndices[0]], evtwt);
            }

            // These are actually obsolete after all the selections
            /*	       if (j == 0 and nref > 1 and dphitp > 2.7) { // Fill dijet
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
          log(LOG_TRACE, "Event " + std::to_string(i + 1) + ": filled inclusive jet histograms");

  } // end of event loop
  log(LOG_INFO, "Finished processing all events");

  // Write output histograms
  log(LOG_INFO, "Writing output histograms");

  for (auto &histrange : _histos) {
    for (auto &h : histrange.second) {
      h->Write();
    }
  }
  eh->Write();
  log(LOG_INFO, "Persisted all histograms to output file");

  // Write and close output file
  outfile->Write();
  log(LOG_INFO, "Wrote " + outputfilename);
  
  // Close file properly - TFile destructor handles all owned histogram cleanup
  outfile->Close();
  delete outfile;
  log(LOG_INFO, "Closed output file and released resources");
}
