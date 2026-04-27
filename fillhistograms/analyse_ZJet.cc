#include <iostream>
using std::cout;
using std::endl;

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iterator>
#include <typeinfo>

#include "configurations.h"
#include "histograms.h"
#include "settings.h"

#include "eventhistograms.h"
#include "helpers.h"
#include "input_config.h"
#include "chain_builder.h"

R__LOAD_LIBRARY(histograms_C.so)
R__LOAD_LIBRARY(eventhistograms_C.so)

#include "JetMETCorrections/Modules/interface/JetResolution.h"
JME::JetResolution *_jer(0);
JME::JetResolutionScaleFactor *_jer_sf(0);
float rho = 0.;

std::mt19937 _mersennetwister;
std::uint32_t _seed = 4;
//_seed = 4;

#if REDOJES == 1
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#endif

bool debug = false;
bool applyjetvetomap = true;

// Helper function to load pthat weights from file
std::map<float, double> LoadPthatWeights(const std::string& weightsFile) {
    std::map<float, double> weights;
    std::ifstream fin(weightsFile);
    if (!fin) {
        cout << "Could not open pthat weights file: " + weightsFile << endl;
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

//I will modify this code that it's originally for Photons+Jets to Muon+Jets
// .
// .
// .
// .
// .
// MC Matt produced:
// /store/user/mnguyen//JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/MuonjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root
// /store/user/mnguyen//JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/MuonjetMadraph/260206_124018/0001/merged_HiForestMiniAOD.root
// /store/user/mnguyen//JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/MuonjetMadraph/260206_124018/0002/merged_HiForestMiniAOD.root
//and the otherones
// /store/user/mnguyen//JEC/PPRefSingleMuon0/SingleMuon0/merged_HiForestMiniAOD.root
// /store/user/mnguyen//JEC/PPRefSingleMuon1/SingleMuon1/merged_HiForestMiniAOD.root
// /store/user/mnguyen//JEC/PPRefSingleMuon2/SingleMuon2/merged_HiForestMiniAOD.root
// /store/user/mnguyen//JEC/PPRefSingleMuon3/SingleMuon3/merged_HiForestMiniAOD.root


// Muon+Jet analysis for L3 residual corrections
// jetTree: jet tree path, e.g. "ak4PFJetAnalyzer/t" or "ak4PFJetAnalyzerSDMuoncut1/t"


//here I need to change a few things (ask)
void analyse_ZJet(string input = "ZSM0",
                       string outputfiletag = "AK4_Zjet",
                       bool isMC = false, bool checkjetid = false,
                       string inputType = "era", int maxFiles = -1,
                       int maxEvents = -1, string outputDir = "",
                       int batchIndex = -1, int totalBatches = 1,
                       string jetPath = "ak4PFJetAnalyzer/t") {

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
        "/eos/home-z/zzaidanc/JetMinPOG/L3ResZJet/";
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
      cerr << "ERROR: Unknown era: " << input << endl;
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
    outputfilename = Form("%s/%s_%s_batch%d_of_%d.root",
                         config.outputDir.c_str(),
                         inputName.c_str(),
                         outputfiletag.c_str(), batchIndex, totalBatches);
  } else {
    outputfilename = Form("%s/%s_%s.root", config.outputDir.c_str(),
                         inputName.c_str(), outputfiletag.c_str());
  }
  if (debug)
    outputfilename = "test.root";

  cout << "Output file: " << outputfilename << endl;

  // Define tree paths
  std::string evtPath = "hiEvtAnalyzer/HiTree";
  std::string triggerPath = "hltanalysis/HltTree";
  std::string skimPath = "skimanalysis/HltTree";
  std::string MuonPath = "muonAnalyzer/MuonTree";

  cout << "Using jet tree: " << jetPath << endl;
  cout << "Building input chains..." << endl;
  TreeChains *chains = BuildChainsFromConfig(config, jetPath, false, true);
  if (!chains || chains->nEntries == 0) {
    cerr << "ERROR: No entries found in input!" << endl;
    return;
  }

  auto evtTree = chains->evtChain;
  auto MuonTree = chains->muonChain;
  auto jetTree = chains->jetChain;
  auto triggerTree = chains->triggerChain;
  auto skimTree = chains->skimChain;

  // Read-ahead cache for faster EOS/xrootd I/O
  const int cacheSize = 50*1024*1024; // 50 MB per chain
  evtTree->SetCacheSize(cacheSize);
  jetTree->SetCacheSize(cacheSize);
  triggerTree->SetCacheSize(cacheSize);
  MuonTree->SetCacheSize(cacheSize);

  // Cuts and weights from event tree
  Int_t hiBin = -1;
  Float_t weight = 1, vz = 0, pthat = 0, evtwt = 1;

  // Muon variables (from muonAnalyzer/MuonTree)
  Int_t nReco = 0;
  std::vector<float> *recoPt = 0;
  std::vector<float> *recoEta = 0;
  std::vector<float> *recoPhi = 0;
  std::vector<int>   *recoCharge = 0;
  std::vector<int>   *recoIDTight = 0;
  std::vector<float> *recoPFChIso = 0;
  std::vector<float> *recoPFPhoIso = 0;
  std::vector<float> *recoPFNeuIso = 0;
  std::vector<float> *recoPFPUIso = 0;

  // MC gen-level muon branches (only used for MC)
  std::vector<int>   *genPID = 0;
  std::vector<int>   *genStatus = 0;
  std::vector<float> *genPt = 0;
  std::vector<float> *genEta = 0;
  std::vector<float> *genPhi = 0;
  std::vector<int>   *genMotherID = 0;

  

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

  // Muon trigger
  Int_t HLT_PPRefL2SingleMu12 = 1;

  //auto triggerTree = (TTree*)inFile->Get(triggerPath.c_str());

  if (!isMC) {
    cout << "Use Muon trigger: HLT_PPRefL2SingleMu12_v6" << endl;
    triggerTree->SetBranchStatus("*", 0);
    triggerTree->SetBranchStatus("HLT_PPRefL2SingleMu12_v6", 1);
    triggerTree->SetBranchAddress("HLT_PPRefL2SingleMu12_v6", &HLT_PPRefL2SingleMu12);
  }

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

  // Set Muon branch addresses (muonAnalyzer/MuonTree)
  MuonTree->SetBranchStatus("*", 0);
  MuonTree->SetBranchStatus("nReco", 1);
  MuonTree->SetBranchStatus("recoPt", 1);
  MuonTree->SetBranchStatus("recoEta", 1);
  MuonTree->SetBranchStatus("recoPhi", 1);
  MuonTree->SetBranchStatus("recoCharge", 1);
  MuonTree->SetBranchStatus("recoIDTight", 1);
  MuonTree->SetBranchStatus("recoPFChIso", 1);
  MuonTree->SetBranchStatus("recoPFPhoIso", 1);
  MuonTree->SetBranchStatus("recoPFNeuIso", 1);
  MuonTree->SetBranchStatus("recoPFPUIso", 1);
  if (isMC) {
    MuonTree->SetBranchStatus("genPID", 1);
    MuonTree->SetBranchStatus("genStatus", 1);
    MuonTree->SetBranchStatus("genPt", 1);
    MuonTree->SetBranchStatus("genEta", 1);
    MuonTree->SetBranchStatus("genPhi", 1);
    MuonTree->SetBranchStatus("genMotherID", 1);
  }
  MuonTree->SetBranchAddress("nReco", &nReco);
  MuonTree->SetBranchAddress("recoPt", &recoPt);
  MuonTree->SetBranchAddress("recoEta", &recoEta);
  MuonTree->SetBranchAddress("recoPhi", &recoPhi);
  MuonTree->SetBranchAddress("recoCharge", &recoCharge);
  MuonTree->SetBranchAddress("recoIDTight", &recoIDTight);
  MuonTree->SetBranchAddress("recoPFChIso", &recoPFChIso);
  MuonTree->SetBranchAddress("recoPFPhoIso", &recoPFPhoIso);
  MuonTree->SetBranchAddress("recoPFNeuIso", &recoPFNeuIso);
  MuonTree->SetBranchAddress("recoPFPUIso", &recoPFPUIso);
  if (isMC) {
    MuonTree->SetBranchAddress("genPID", &genPID);
    MuonTree->SetBranchAddress("genStatus", &genStatus);
    MuonTree->SetBranchAddress("genPt", &genPt);
    MuonTree->SetBranchAddress("genEta", &genEta);
    MuonTree->SetBranchAddress("genPhi", &genPhi);
    MuonTree->SetBranchAddress("genMotherID", &genMotherID);
  }

  //TODO: Add electron veto for Zs

  auto pthatWeights = LoadPthatWeights("jecfiles/2024_PP_private_test_weights.txt");
  std::vector<float> pthatBins;
  for (const auto& kv : pthatWeights) pthatBins.push_back(kv.first);
  std::sort(pthatBins.begin(), pthatBins.end());

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
  cout << "Applying MC JEC from file " << jecfile.c_str() << endl;
  FactorizedJetCorrector *corr;
  vector<JetCorrectorParameters> vpar;
  // This is MCTruth (L2Relative)
  vpar.push_back(JetCorrectorParameters(jecfile.c_str()));
  // L2 residual for data only
  if (!isMC) {
    cout << "Applying L2 Residual from file " << l2file.c_str() << endl;
    vpar.push_back(JetCorrectorParameters(l2file.c_str()));
  }
  corr = new FactorizedJetCorrector(vpar);
#endif

  // JER not needed for Z+jet L3 residual analysis

  // Jet veto map
  auto mapfile = new TFile("jecfiles/Summer24Prompt24_RunBCDEFGHI.root","READ"); 
  auto vetomap = (TH2D*)mapfile->Get("jetvetomap_all");
  
  // Null check for veto map
  if (applyjetvetomap && (!mapfile || mapfile->IsZombie() || !vetomap)) {
    cerr << "WARNING: Veto map unavailable, disabling veto map selection" << endl;
    applyjetvetomap = false;
  }

  cout << "Number of entries :" << chains->nEntries << endl;
  Long64_t nentries = chains->nEntries;
  if (config.maxEvents > 0 && config.maxEvents < nentries) {
    nentries = config.maxEvents;
  }
  if (debug && nentries > 1000)
    nentries = 1000;

  // Event-level batch splitting: split events within a single file across jobs
  Long64_t firstEvent = 0;
  if (batchIndex >= 0 && totalBatches > 1 && maxFiles <= 0) {
    Long64_t eventsPerBatch = nentries / totalBatches;
    firstEvent = batchIndex * eventsPerBatch;
    Long64_t lastEvent = (batchIndex == totalBatches - 1) ? nentries : firstEvent + eventsPerBatch;
    nentries = lastEvent - firstEvent;
    cout << "Batch " << batchIndex << "/" << totalBatches
         << ": processing events " << firstEvent << " to " << firstEvent + nentries << endl;
  }

  cout << "Processing " << nentries << " events" << endl;
  for (Long64_t i = firstEvent; i < firstEvent + nentries; ++i) {
    if ((i - firstEvent) % 1000000 == 0)
      cout << "Event " << (i - firstEvent) << " / " << nentries << endl;
    evtTree->GetEntry(i);
    triggerTree->GetEntry(i);
    MuonTree->GetEntry(i);

    // Z trigger logic
     if (!isMC) {
       trigger = HLT_PPRefL2SingleMu12;
     }
     if (isMC) trigger = true; // MC: no trigger requirement

     if (!trigger) continue;

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

    // cout << weight << " " << evtwt << endl;

    // BASIC EVENT FILTERS
    //  if (!isMC) {
    //    skimTree->GetEntry(i);
    //    if (pprimaryVertexFilter != 1) continue;
    //  }
    jetTree->GetEntry(i);

    if (debug) {
      cout << "Processing event " << i << ", evt number: " << evt << endl;
    }

    // Z+jet: need at least 1 muon and 1 jet
    if (nReco < 1)
      continue;
    if (nref < 1)
      continue;
    if (jtpt[0] < jtptmin) {
      continue;
    }

    eh->event_vz->Fill(vz, evtwt);
    if (isMC)
      eh->event_pthatwsgenweight->Fill(pthat, weight);

    // This is Z+jet analysis
    double Muon_pt, Muon_eta, Muon_phi, jet_pt, jet_eta, jet_phi, ptavgtp,
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

    // Apply JEC
    for (int j = 0; j < nref; ++j) {
      jtpt_uncorr[j] = jtpt[j];

#if REDOJES == 1
      // cout << "Applying JES" << endl;
      corr->setJetPt(jtpt[j]);
      // corr->setJetE(jteu[jetidx]);
      corr->setJetEta(jteta[j]);
      // 	 corr->setJetPhi(jthpi[j]);

      vector<float> v = corr->getSubCorrections();
      float jes = v.back();

      //	 cout << "New jes correction jet pt: " << jtpt[j] << " " <<
      //jteta[j] << " "  << jes << endl;
      jtpt[j] *= jes;
#endif

      // JER not applied for photon+jet L3 residual analysis
    }

    // ========================================
    // PHOTON+JET SELECTION
    // ========================================

    // 1. Find leading photon passing selection
    int leadMuonIdx = -1;
    float leadMuonPt = 0;
    int leadMuonGenIdx = -1;

    for (int iMuon = 0; iMuon < nReco; iMuon++) {
      int currentGenIdx = -1;
      // Kinematic cuts
      if ((*recoPt)[iMuon] < 20.0)
        continue; // Minimum muon pT
      if (fabs((*recoEta)[iMuon]) > 2.4)
        continue; // Muon acceptance

      // Muon ID and isolation
      if (!(*recoIDTight)[iMuon])
        continue;
      float iso = (*recoPFChIso)[iMuon] +
                  std::max(0.f, (*recoPFNeuIso)[iMuon] + (*recoPFPhoIso)[iMuon]
                                - 0.5f * (*recoPFPUIso)[iMuon]);
      if (iso / (*recoPt)[iMuon] > 0.15)
        continue;

      // MC gen-matching: find closest gen muon (PID=13) within dR < 0.3
      if (isMC) {
        if (!genPt) continue;
        float minDR = 0.3;
        for (int ig = 0; ig < (int)genPt->size(); ig++) {
          if (!genPID || !genStatus) continue;
          if (abs((*genPID)[ig]) != 13) continue;
          if ((*genStatus)[ig] != 1) continue;
          if ((*genPt)[ig] < 10.0) continue;
          float deta_gm = (*recoEta)[iMuon] - (*genEta)[ig];
          float dphi_gm = (*recoPhi)[iMuon] - (*genPhi)[ig];
          if (dphi_gm > TMath::Pi()) dphi_gm -= 2 * TMath::Pi();
          if (dphi_gm < -TMath::Pi()) dphi_gm += 2 * TMath::Pi();
          float dR_gm = sqrt(deta_gm * deta_gm + dphi_gm * dphi_gm);
          if (dR_gm < minDR) {
            minDR = dR_gm;
            currentGenIdx = ig;
          }
        }
        if (currentGenIdx < 0) continue;
      }

      // Find highest pT muon
      if ((*recoPt)[iMuon] > leadMuonPt) {
        leadMuonPt = (*recoPt)[iMuon];
        leadMuonIdx = iMuon;
        if (isMC) leadMuonGenIdx = currentGenIdx;
      }
    }

    // Stop early if no Z passes the kinematic preselection
    if (leadMuonIdx < 0)
      continue;

    // Muon ID and isolation already applied in selection loop above

    // 2. Find leading and subleading away-side jets
    // First: identify all jets back-to-back with photon (dphi > 2.7)// No just apply a small dR requirement
    int awayJetIndices[MAXJETS];
    int nAwayJets = 0;

    for (int j = 0; j < nref; j++) {
      // Apply jet ID
      if (checkjetid && passjetid[j] == 0)
        continue;

      // Jet kinematic cuts
      if (jtpt[j] < 20.0)
        continue; // Minimum jet pT

      // Calculate delta-phi with Muon
      float dphi = abs(jtphi[j] - (*recoPhi)[leadMuonIdx]);
      if (dphi > TMath::Pi())
        dphi = 2 * TMath::Pi() - dphi;

      // Back-to-back requirement
      if (dphi < 2.7488935)
        continue; // 2*pi/3 = 2.0943951

      // Calculate delta-R (reject jets close to Muon)
      float deta = jteta[j] - (*recoEta)[leadMuonIdx];
      float deltaR = sqrt(deta * deta + dphi * dphi);
      if (deltaR < 0.4)
        continue; // Isolation cone

      // This jet is away-side
      awayJetIndices[nAwayJets] = j;
      nAwayJets++;
    }

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
      float ptavg_temp = (*recoPt)[leadMuonIdx];
      alpha = jtpt[secondAwayJetIdx] / ptavg_temp;
    } else {
      alpha = 0; // Only one away-side jet
    }

    // 4. Calculate Muon+jet variables
    Muon_pt = (*recoPt)[leadMuonIdx];
    Muon_eta = (*recoEta)[leadMuonIdx];
    Muon_phi = (*recoPhi)[leadMuonIdx];

    float genMuonPt = -1;
    float genMuonEta = 0;
    float genMuonPhi = 0;
    bool hasGenMuon = false;
    if (isMC && leadMuonGenIdx >= 0 && genPt && genEta && genPhi) {
      if (leadMuonGenIdx < (int)genPt->size()) {
        genMuonPt = (*genPt)[leadMuonGenIdx];
        genMuonEta = (*genEta)[leadMuonGenIdx];
        genMuonPhi = (*genPhi)[leadMuonGenIdx];
        hasGenMuon = true;
      }
    }

    jet_pt = jtpt[awayJetIdx];
    jet_eta = jteta[awayJetIdx];
    jet_phi = jtphi[awayJetIdx];

    float dphi_Muonjet = abs(jet_phi - Muon_phi);
    if (dphi_Muonjet > TMath::Pi())
      dphi_Muonjet = 2 * TMath::Pi() - dphi_Muonjet;

    ptavgtp = Muon_pt;
    balance = jet_pt / Muon_pt; // Response
    asymmtp = balance;            // For compatibility with histogram filling

    // Apply jet veto map to photon and leading/subleading away-side jets
    if (applyjetvetomap && vetomap) {
      bool passVetoMap = true;
      
      // Check Muon position - Not required for Muon since the vetomap is mostly due to pixel failures
      // int Muon_bin = vetomap->FindBin(Muon_eta, Muon_phi);
      // if (vetomap->GetBinContent(Muon_bin) > 0) passVetoMap = false;
      
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

    // cout << "TP:" << tagpt << " " << probept << " " << alpha << endl;
    // ========================================
    // FILL PHOTON+JET HISTOGRAMS
    // ========================================

    for (auto &histrange : _histos) {
      for (auto &h : histrange.second) {

        // Check if jet eta and centrality are in range for this histogram set
        if (jet_eta >= h->etamin && jet_eta < h->etamax &&
            hiBin >= h->hibinmin && hiBin < h->hibinmax) {

          // Muon properties
          h->Muon_pt->Fill(Muon_pt, evtwt);
          h->Muon_eta->Fill(Muon_eta, evtwt);
          h->Muon_phi->Fill(Muon_phi, evtwt);
          float leadRelIso = ((*recoPFChIso)[leadMuonIdx] +
                              std::max(0.f, (*recoPFNeuIso)[leadMuonIdx] +
                                           (*recoPFPhoIso)[leadMuonIdx] -
                                           0.5f * (*recoPFPUIso)[leadMuonIdx])) /
                             (*recoPt)[leadMuonIdx];
          h->Muon_HoverE->Fill(leadRelIso, evtwt);         // reused for rel. isolation
          h->Muon_sigmaIetaIeta->Fill(leadRelIso, evtwt);  // reused for rel. isolation

          if (isMC && hasGenMuon) {
            if (h->genMuon_pt) h->genMuon_pt->Fill(genMuonPt, evtwt);
            if (h->genMuon_eta) h->genMuon_eta->Fill(genMuonEta, evtwt);
            if (h->genMuon_phi) h->genMuon_phi->Fill(genMuonPhi, evtwt);
            if (h->Muonresponse && genMuonPt > 0) h->Muonresponse->Fill(genMuonPt, Muon_pt / genMuonPt, evtwt);
            if (h->Muon_ptres && genMuonPt > 0) h->Muon_ptres->Fill((Muon_pt - genMuonPt) / genMuonPt, evtwt);
          }

          // Away-side jet properties
          h->awayside_jet_pt->Fill(jet_pt, evtwt);
          h->awayside_jet_eta->Fill(jet_eta, evtwt);
          h->awayside_jet_phi->Fill(jet_phi, evtwt);
          h->awayside_jet_uncorr_pt->Fill(jtpt_uncorr[awayJetIdx], evtwt);

          // Muon+Jet system
          h->Muonjet_dphi->Fill(dphi_Muonjet, evtwt);
          h->Muonjet_balance->Fill(balance, evtwt);
          h->Muonjet_ptavg->Fill(ptavgtp, evtwt);
          h->Muonjet_alpha->Fill(alpha, evtwt);

          // Trigger histograms
          if (HLT_PPRefL2SingleMu12) {
            h->HLT_PPRefL2SingleMu12->Fill(1, evtwt);
            h->HLT_PPRefL2SingleMu12_ptav->Fill(ptavgtp, evtwt);
          }

          // Alpha-dependent balance profiles (analogous to dijet asymmetry)
          if (alpha < 0.1) {
            h->Muonjet_balance_a01->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a01->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.2) {
            h->Muonjet_balance_a02->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a02->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.3) {
            h->Muonjet_balance_a03->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a03->Fill(ptavgtp, jet_eta, balance, evtwt);
            // Jet composition for alpha < 0.3 (like dijets)
            h->jet_nef->Fill(jet_pt, jtnef[awayJetIdx], evtwt);
            h->jet_cef->Fill(jet_pt, jtcef[awayJetIdx], evtwt);
            h->jet_nhf->Fill(jet_pt, jtnhf[awayJetIdx], evtwt);
            h->jet_chf->Fill(jet_pt, jtchf[awayJetIdx], evtwt);
            h->jet_muf->Fill(jet_pt, jtmuf[awayJetIdx], evtwt);
          }
          if (alpha < 0.4) {
            h->Muonjet_balance_a04->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a04->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.5) {
            h->Muonjet_balance_a05->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a05->Fill(ptavgtp, jet_eta, balance, evtwt);
          }
          if (alpha < 0.6) {
            h->Muonjet_balance_a06->Fill(ptavgtp, balance, evtwt);
            h->Muonjet_balance2D_a06->Fill(ptavgtp, jet_eta, balance, evtwt);
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
                h->Muonjet_balance3D->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->Muonjet_balance3Dwide->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->Muonjet_balance3Dnarrow->Fill(ptavgtp, jet_eta, alphaFillValue, balance, evtwt);
                h->Muonjet_balance3Dabseta->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->Muonjet_balance3Dabsetawide->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                h->Muonjet_balance3Dabsetanarrow->Fill(ptavgtp, abs(jet_eta), alphaFillValue, balance, evtwt);
                
                // Fill counts histograms (UNWEIGHTED - just count entries)
                if (h->Muonjet_balance3D_counts) h->Muonjet_balance3D_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->Muonjet_balance3Dwide_counts) h->Muonjet_balance3Dwide_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->Muonjet_balance3Dnarrow_counts) h->Muonjet_balance3Dnarrow_counts->Fill(ptavgtp, jet_eta, alphaFillValue);
                if (h->Muonjet_balance3Dabseta_counts) h->Muonjet_balance3Dabseta_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->Muonjet_balance3Dabsetawide_counts) h->Muonjet_balance3Dabsetawide_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                if (h->Muonjet_balance3Dabsetanarrow_counts) h->Muonjet_balance3Dabsetanarrow_counts->Fill(ptavgtp, abs(jet_eta), alphaFillValue);
                
                // Fill balance distribution (Muon_pT, alpha, balance_value)
                // Fill with weight for each cumulative alpha cut
                if (h->Muonjet_balance_dist) h->Muonjet_balance_dist->Fill(ptavgtp, alphaFillValue, balance, evtwt);
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

          if (jteta[j] >= h->etamin and jteta[j] < h->etamax and
              hiBin >= h->hibinmin and hiBin < h->hibinmax) {

            if (checkjetid and passjetid[j] == 0)
              continue;

            h->jetetaphi->Fill(jteta[j], jtphi[j], weight);

            if (j == 0 and passjetid[0]) {
              // Muon trigger check; leading jet pt
              if (HLT_PPRefL2SingleMu12)
                h->HLT_PPRefL2SingleMu12->Fill(jtpt[awayJetIndices[0]], evtwt);
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
  } // end of event loop

  // Write output histograms

  for (auto &histrange : _histos) {
    for (auto &h : histrange.second) {
      h->Write();
    }
  }
  eh->Write();

  // Write and close output file
  outfile->Write();
  cout << "Wrote " << outputfilename.c_str() << endl;
  
  // Close file properly - TFile destructor handles all owned histogram cleanup
  outfile->Close();
  delete outfile;
}
