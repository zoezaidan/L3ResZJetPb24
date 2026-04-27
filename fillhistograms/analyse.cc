#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include <iterator>
#include "TMath.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include "TFile.h"
#include <typeinfo>
#include "TTree.h"

#include "histograms.h"
#include "settings.h"
#include "configurations.h"

#include "eventhistograms.h"
#include "helpers.h"
#include "input_config.h"
#include "chain_builder.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
JME::JetResolution *_jer(0);
JME::JetResolutionScaleFactor *_jer_sf(0);
float rho = 0.;

std::mt19937 _mersennetwister;
std::uint32_t _seed = 4;

#if REDOJES == 1
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#endif

map<string, vector<histograms*> > _histos;

bool debug = false;
bool applyjetvetomap = false;


//void analyse(string era = "RERECOMC", string outputfiletag = "AK4_nojetid", bool isMC = true, bool checkjetid = false, bool iszb = false, bool dol2res = false, bool dojer = false, bool fillforJER = false) {
void analyse(string input = "RERECOHP", string outputfiletag = "AK4_nojetid", bool isMC = false, bool checkjetid = false, bool iszb = false, bool dol2res = false, bool dojer = false, bool fillforJER = false,float jtptlimitforalpha = 15, string inputType = "era", int maxFiles = -1, int maxEvents = -1, string outputDir = "", int batchIndex = -1, int totalBatches = 1, string jetPath = "ak4PFJetAnalyzer/t") {

  bool usecalotrig = false;
  bool checkvalidjet = false; // this is for checking valid jet range after applying l2. now for tightly limited range. TODO: do something smarter

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
    config.outputDir = "/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/L2ResDiJet";
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
  if (inputType != "era") {
    size_t lastSlash = input.find_last_of("/");
    if (lastSlash != string::npos && lastSlash < input.size() - 1) {
      inputName = input.substr(lastSlash + 1);
    }
  }

  if (batchIndex >= 0) {
    outputfilename = Form("%s/%s_%s_batch%d_of_%d.root", config.outputDir.c_str(), inputName.c_str(), outputfiletag.c_str(), batchIndex, totalBatches);
  } else {
    outputfilename = Form("%s/%s_%s.root", config.outputDir.c_str(), inputName.c_str(), outputfiletag.c_str());
  }
  //string outputfilename = Form("/eos/user/l/lamartik/HIJEC_rereco_results_HI2023MCTruth_chs/%s_%s.root",era.c_str(),outputfiletag.c_str());
  if (debug) outputfilename = "test.root";
  
  TRandom3 r;
  // Define and activate branches
  std::string evtPath = "hiEvtAnalyzer/HiTree";
  std::string triggerPath = "hltanalysis/HltTree";
  std::string skimPath = "skimanalysis/HltTree";
//   std::string egmPath = "ggHiNtuplizer/EventTree"; 

  cout << "Using jet tree: " << jetPath << endl;

  cout << "Building input chains..." << endl;
  TreeChains *chains = BuildChainsFromConfig(config, jetPath, false);
  if (!chains || chains->nEntries == 0) {
    cerr << "ERROR: No entries found in input!" << endl;
    return;
  }

  auto evtTree = chains->evtChain;
  auto triggerTree = chains->triggerChain;
  auto skimTree = chains->skimChain;
  auto jetTree = chains->jetChain;

  const bool hasTriggerTree = (triggerTree && triggerTree->GetNtrees() > 0);
  if (!isMC && !hasTriggerTree) {
    cerr << "ERROR: Data input is missing hltanalysis/HltTree" << endl;
    return;
  }
  
  // Cuts and weights from event tree
  Int_t       hiBin = -1;
  Float_t     weight = 1, vz = 0, pthat = 0, evtwt = 1;
  
  evtTree->SetBranchAddress("hiBin", &hiBin);
  evtTree->SetBranchAddress("vz", &vz);
  if (isMC) evtTree->SetBranchAddress("weight", &weight);
  if (isMC) evtTree->SetBranchAddress("pthat", &pthat);

  evtTree->SetBranchStatus("*",0);
  evtTree->SetBranchStatus("hiBin",1);
  evtTree->SetBranchStatus("vz",1);
  if (isMC) evtTree->SetBranchStatus("weight",1);
  if (isMC) evtTree->SetBranchStatus("pthat",1);

  //// EVENT FILTERS
  if (!isMC) skimTree->SetBranchStatus("*",1);

  Int_t pprimaryVertexFilter = 1;
  if (!isMC) skimTree->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  
  Int_t trigger = 0;
  bool useTriggerSelection = false;

  // Triggger paths in the files
  Int_t HLT_ZB = 0, HLT_40 = 0, HLT_60 = 0, HLT_80 = 0, HLT_100 = 0, HLT_120 = 0;
 
  //  if (isMC) triggerTree->SetBranchAddress("HLT_PPRefZeroBias_v1",&HLT_ZB);
  
  if (!usecalotrig and !isMC) {
    bool hasPFTriggers =
      triggerTree->GetBranch("HLT_AK4PFJet40_v1") &&
      triggerTree->GetBranch("HLT_AK4PFJet60_v1") &&
      triggerTree->GetBranch("HLT_AK4PFJet100_v1") &&
      triggerTree->GetBranch("HLT_AK4PFJet120_v1");

    if (hasPFTriggers) {
      cout << "Use PF triggers" << endl;
      useTriggerSelection = true;
  
      triggerTree->SetBranchAddress("HLT_AK4PFJet40_v1",&HLT_40);
      triggerTree->SetBranchAddress("HLT_AK4PFJet60_v1",&HLT_60);
      if (triggerTree->GetBranch("HLT_AK4PFJet80_v1")) triggerTree->SetBranchAddress("HLT_AK4PFJet80_v1",&HLT_80);
      triggerTree->SetBranchAddress("HLT_AK4PFJet100_v1",&HLT_100);
      triggerTree->SetBranchAddress("HLT_AK4PFJet120_v1",&HLT_120);

      triggerTree->SetBranchStatus("*",0);
      if (triggerTree->GetBranch("HLT_PPRefZeroBias_v1")) triggerTree->SetBranchStatus("HLT_PPRefZeroBias_v1",1);
      triggerTree->SetBranchStatus("HLT_AK4PFJet40_v1",1);
      triggerTree->SetBranchStatus("HLT_AK4PFJet60_v1",1);
      if (triggerTree->GetBranch("HLT_AK4PFJet80_v1")) triggerTree->SetBranchStatus("HLT_AK4PFJet80_v1",1);
      triggerTree->SetBranchStatus("HLT_AK4PFJet100_v1",1);
      triggerTree->SetBranchStatus("HLT_AK4PFJet120_v1",1);
    } else {
      cout << "WARNING: Expected PF trigger branches not found; disabling trigger selection for this input." << endl;
    }
    
  }
  if (usecalotrig and !isMC) {
    bool hasCaloTriggers =
      triggerTree->GetBranch("HLT_AK4CaloJet40_v1") &&
      triggerTree->GetBranch("HLT_AK4CaloJet60_v1") &&
      triggerTree->GetBranch("HLT_AK4CaloJet100_v1") &&
      triggerTree->GetBranch("HLT_AK4CaloJet120_v1");

    if (hasCaloTriggers) {
      cout << "Use Calo triggers" << endl;
      useTriggerSelection = true;
    
    triggerTree->SetBranchAddress("HLT_AK4CaloJet40_v1",&HLT_40);
    triggerTree->SetBranchAddress("HLT_AK4CaloJet60_v1",&HLT_60); 
    //triggerTree->SetBranchAddress("HLT_AK4CaloJet80_v1",&HLT_80); 
    triggerTree->SetBranchAddress("HLT_AK4CaloJet100_v1",&HLT_100);
    triggerTree->SetBranchAddress("HLT_AK4CaloJet120_v1",&HLT_120);
    
    triggerTree->SetBranchStatus("*",0);
    
    if (triggerTree->GetBranch("HLT_PPRefZeroBias_v1")) triggerTree->SetBranchStatus("HLT_PPRefZeroBias_v1",1);
    triggerTree->SetBranchStatus("HLT_AK4CaloJet40_v1",1);
    triggerTree->SetBranchStatus("HLT_AK4CaloJet60_v1",1);
    //triggerTree->SetBranchStatus("HLT_AK4CaloJet80_v1",1);
    triggerTree->SetBranchStatus("HLT_AK4CaloJet100_v1",1);
    triggerTree->SetBranchStatus("HLT_AK4CaloJet120_v1",1);
    } else {
      cout << "WARNING: Expected Calo trigger branches not found; disabling trigger selection for this input." << endl;
    }
    
  }

  // JETS 
  jetTree->SetBranchStatus("*",1);    

  Int_t     evt;
  
  // Reconstruted jet information
  Int_t     nref;
  Float_t   jtpt[MAXJETS];
  Float_t   jtpt_uncorr[MAXJETS];
  Float_t   jteta[MAXJETS];
  Float_t   jtphi[MAXJETS];

  Float_t   jtnhf[MAXJETS];
  Float_t   jtchf[MAXJETS];
  Float_t   jtnef[MAXJETS];
  Float_t   jtcef[MAXJETS];
  Float_t   jtmuf[MAXJETS];

  Int_t   jtchm[MAXJETS]; // charged multiplicity

  Int_t   jtn[MAXJETS];
 
  jetTree->SetBranchAddress("evt", &evt);
  jetTree->SetBranchAddress("nref", &nref);
  jetTree->SetBranchAddress("jtpt", &jtpt);
  jetTree->SetBranchAddress("jteta", &jteta);
  jetTree->SetBranchAddress("jtphi", &jtphi);

  jetTree->SetBranchAddress("jtPfNHF", &jtnhf);
  jetTree->SetBranchAddress("jtPfCHF", &jtchf);
  jetTree->SetBranchAddress("jtPfNEF", &jtnef);
  jetTree->SetBranchAddress("jtPfCEF", &jtcef);
  jetTree->SetBranchAddress("jtPfMUF", &jtmuf);

  jetTree->SetBranchAddress("jtPfCHM", &jtchm);
  
  // Gen level jet information
  Float_t   jtpt_gen[MAXJETS];
  Float_t   jteta_gen[MAXJETS];
  Float_t   jtphi_gen[MAXJETS];
  Float_t   refdrjt[MAXJETS];

  if (isMC) {
    jetTree->SetBranchAddress("refpt", &jtpt_gen);
    jetTree->SetBranchAddress("refeta", &jteta_gen);
    jetTree->SetBranchAddress("refphi", &jtphi_gen);
    jetTree->SetBranchAddress("refdrjt", &refdrjt);
  }

  TFile *outfile = new TFile(outputfilename.c_str(),"RECREATE");

  // Create folders for centrality bins. Mostly a placeholder in case PbPb MC/data is checked.
  for (int j = 0; j < nhibins; ++j) { 
  
      if (hibins[j] < hibins[j+1]) {
	string name = Form("hibin_%.1f_%.1f",hibins[j],hibins[j+1]);
	outfile->mkdir(name.c_str());
	
	TDirectory *dir = outfile->GetDirectory(name.c_str()); assert(dir);
	dir->cd();
    
	for (int i = 0; i < netabins ; ++i) { 
	  if (etaedges[i] < etaedges[i+1]) {
	    string name2 = Form("eta_%.1f_%.1f",etaedges[i],etaedges[i+1]);
	    dir->mkdir(name2.c_str());
	
	    TDirectory *dir2 = dir->GetDirectory(name2.c_str()); assert(dir2);
	    dir2->cd();
	    
	    histograms *h = new histograms(dir2, etaedges[i], etaedges[i+1], hibins[j],hibins[j+1], isMC);
	    _histos[name2.c_str()].push_back(h);
	  }
	}
      }
  }

  outfile->mkdir("event");
  TDirectory *dir = outfile->GetDirectory("event"); assert(dir);
  dir->cd();

  eventhistograms *eh = new eventhistograms(dir, isMC);
   

#if REDOJES == 1
  cout << "Applying MC JEC from file " << jecfile.c_str() << endl;
  FactorizedJetCorrector* corr;
  vector<JetCorrectorParameters> vpar;
  // This is MCTruth
  vpar.push_back(JetCorrectorParameters(jecfile.c_str()));
  // L2 residual
  if (!isMC and dol2res) vpar.push_back(JetCorrectorParameters(l2file.c_str()));
  if (dol2res) cout << "Applying L2 residual" << endl;
  corr = new FactorizedJetCorrector(vpar);
#endif

  if (dojer) {
    cout << "Applying JER SF" << endl;

    _jer = new JME::JetResolution(resolutionFile);
    _jer_sf =  new JME::JetResolutionScaleFactor(scaleFactorFile);
    
  }

  // Jet veto map
//   auto mapfile = new TFile("jecfiles/Summer23BPixPrompt23_RunD_v1.root","READ");
//   auto vetomap = (TH2D*)mapfile->Get("jetvetomap_all");

  
   cout << "Number of entries :" <<  chains->nEntries  << endl; 
   Long64_t nentries = chains->nEntries;
   if (config.maxEvents > 0 && config.maxEvents < nentries) {
     nentries = config.maxEvents;
   }
   if (debug && nentries > 1000) nentries = 1000;

   cout << "Processing " << nentries << " events" << endl;
   for (Long64_t i = 0; i < nentries; ++i) {
    evtTree->GetEntry(i);
    if (!isMC) triggerTree->GetEntry(i);

     //trigger = HLT_ZB or HLT_40 or HLT_60;

     trigger = true;
     if (!isMC && useTriggerSelection) {
       if (iszb) trigger = HLT_ZB;
       else trigger = HLT_60;
     }

     if (!trigger) continue;
     
     evtwt = 1;
     if (isMC) {
       evtwt *=  weight;
     }

     // cout << weight << " " << evtwt << endl;
          
     // BASIC EVENT FILTERS
     if (!isMC) {
       skimTree->GetEntry(i);
       if (pprimaryVertexFilter != 1) continue;
     }
     jetTree->GetEntry(i);

     // Filter out events without jets with pt > 10 GeV
     if (nref < 2) continue;
     if (jtpt[0] < jtptmin) {
       continue;
     }

     eh->event_vz->Fill(vz, evtwt);
     if (isMC)  eh->event_pthatwsgenweight->Fill(pthat,weight);

     // This is dijet with tag and probe
     double tagpt, probept, tageta, probeeta, pt3, ptavgtp, alpha, asymmtp;
     double tagpt_gen, probept_gen, tageta_gen, probeeta_gen, pt3_gen, ptavgtp_gen, alpha_gen, asymmtp_gen;
     double djrespasymm; 

     if (checkvalidjet) {   // This is based on validity of JEC. - obsolete?
       for (int j = 0; j < nref; ++j ) {
	 if (abs(jteta[j]) > 2.964) jtpt[j] = 0;   // Always invalid jets
	 
	 else if (abs(jteta[j]) > 2.5 and jtpt[j] > 120) jtpt[j] = 0;  // ?
	 else if (abs(jteta[j]) > 1.93 and jtpt[j] > 170) jtpt[j] = 0;  // ? 
	 
	 if (jtpt[j] < 80) jtpt[j] = 0;
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
	        if (jtnhf[j] >= 0.99) passjetid[j] = false;
		if (jtnef[j] >= 0.9) passjetid[j] = false;
		if (jtchf[j] <= 0.01) passjetid[j] = false;
		if (jtcef[j] >= 0.8) passjetid[j] = false;
		if (jtmuf[j] >= 0.8) passjetid[j] = false;
		if (jtchm[j] <= 0) passjetid[j] = false;
	  }
	  else if (abs(jteta[j]) <= 2.7) {
             if (jtnhf[j] >= 0.9) passjetid[j] = false;
	     if (jtnef[j] >= 0.99) passjetid[j] = false;
	     if (jtmuf[j] >= 0.8) passjetid[j] = false;
	     if (jtcef[j] >= 0.8) passjetid[j] = false;

          }  
          else if (abs(jteta[j]) <= 3.0) {
               if (jtnhf[j] >= 0.99) passjetid[j] = false;
	       if (jtnef[j] >= 0.99) passjetid[j] = false;
          }
          else if (abs(jteta[j]) <= 5.0)  {
	    if (jtnef[j] >= 0.4) passjetid[j] = false;

          }
	   //	   if (jtpt[j] > jtptmin)	   cout << passjetid[j] << endl;
	   //  if (passjetid[j] < 2 and nref > 2  and jtpt[j] > 70  and jtpt[1] > 40) cout << "Pass jetid: " << passjetid[j] << " pt: " << jtpt[j] << " " << jteta[j] << " " << j << " " << i <<  endl;
	 }
       }

     
     // Apply JEC
     for (int j = 0; j < nref; ++j ) {
 	 jtpt_uncorr[j] = jtpt[j];

#if REDOJES == 1
	 // cout << "Applying JES" << endl;
	 corr->setJetPt(jtpt[j]);
	 // corr->setJetE(jteu[jetidx]);
	 corr->setJetEta(jteta[j]);
	 // 	 corr->setJetPhi(jthpi[j]);

	 vector<float> v = corr->getSubCorrections();
	 float jes = v.back();

	 //	 cout << "New jes correction jet pt: " << jtpt[j] << " " << jteta[j] << " "  << jes << endl;
	 jtpt[j] *= jes;
#endif

	 if (isMC and dojer) {
	   double jet_resolution = _jer->getResolution({{JME::Binning::JetPt, jtpt[j]}, {JME::Binning::JetEta, jteta[j]}, {JME::Binning::Rho, rho}});
	   double jer_sf = _jer_sf->getScaleFactor({{JME::Binning::JetEta, jteta[j]}, {JME::Binning::JetPt, jtpt[j]}}, Variation::NOMINAL);

	   float jersfcorr = 1;
	   // TODO: not hardcoded radius
	   // Also check resolution
	   if (refdrjt[j] != -999 and refdrjt[j] < 0.2 and abs((jtpt[j]-jtpt_gen[j])) < 3*jet_resolution*jtpt[j] ) {
	     jersfcorr += (jer_sf-1)*(jtpt[j]-jtpt_gen[j])/jtpt[j];
	     //	     cout << "We use scaling method " <<  abs((jtpt[j]-jtpt_gen[j])) << endl;
	   }
	   else {
	     //	     cout << "We use stochastic method " <<  abs((jtpt[j]-jtpt_gen[j])) << endl;
	     _mersennetwister = std::mt19937(_seed);
	     //	     double sigma = std::sqrt(std::max(jer_sf*jer_sf - 1,0)); // technically should be max(sf*sf-1,0)
	     std::normal_distribution<> d(0, jet_resolution);
	     if (jer_sf*jer_sf > 1) jersfcorr += d(_mersennetwister)*std::sqrt(jer_sf*jer_sf-1);
	     else jersfcorr += d(_mersennetwister)*std::sqrt(0);

	   }
	   
	   //	   cout << jet_resolution << " SF " << jer_sf << " eta: " << jteta[j] << " corr:  " << jersfcorr << " reco: " << jtpt[j] << " gen " << jtpt_gen[j] << " matching " << refdrjt[j] <<  endl;
	   jtpt[j] *= jersfcorr;
	 }

     }
     int ind1 = -1, ind2 = -1, ind3 = -1;
     //     int ind[3] = {0, 1, -1};
     int ind[10] = {0, 1, -1, -1, -1, -1, -1, -1, -1, -1};
     
     if (nref > 1) {
   
       // Find highest pt
       float highestpt = 0;
       for (int scan = 0; scan < nref; scan ++) {
	 if (passjetid[scan] == 1)  {
	   if (jtpt[scan] > highestpt) {
	     highestpt = jtpt[scan];
	     ind[0] = scan;
	   }
	 }
       }

       // Second highest
       float sechighestpt = 0;
       for (int scan = 0; scan < nref; scan ++) {
	 if (passjetid[scan] == 1)  {
	   if (jtpt[scan] < highestpt and jtpt[scan] > sechighestpt) {
	     sechighestpt = jtpt[scan];
	     ind[1] = scan;
	   }
	 }
       }

       // Third
       float thirhighestpt = 0;
       for (int scan = 0; scan < nref; scan ++) {
	 if (passjetid[scan] == 1)  {
	   if (jtpt[scan] < sechighestpt and jtpt[scan] > thirhighestpt) {
	     thirhighestpt = jtpt[scan];
	     ind[2] = scan;
	   }
	 }
       }

      
       if (ind[1] == -1) continue; // Ditch events with only 1 good jet

       /*       cout << "Good corr:" << jtpt[ind[0]] << " " << jtpt[ind[1]] << " " << jtpt[ind[2]] << " " << jtpt[ind[3]] << " " << jtpt[ind[4]] << " nref " << nref << endl;
       cout << "Good corr:" << highestpt << " " << sechighestpt << " " << thirhighestpt << " nref " << nref << endl;
       cout << "Good uncorr:" << jtpt_uncorr[ind[0]] << " " << jtpt_uncorr[ind[1]] << " " << jtpt_uncorr[ind[2]] << " " << jtpt_uncorr[ind[3]] << " " << jtpt_uncorr[ind[4]] << " nref " << nref << endl;
       */

       if (jtpt_uncorr[ind[0]] < jtptmin or jtpt_uncorr[ind[1]] < jtptmin) continue; // Didn't find dijets with good pt

       int probeind = 0;
       int tagind = -1;

       if (applyjetvetomap) {
	//  if (vetomap->GetBinContent(vetomap->FindBin(jteta[0],jtphi[0])) > 0 or vetomap->GetBinContent(vetomap->FindBin(jteta[1],jtphi[1])) > 0) continue;
       }

       
       if (fillforJER) {
          const auto rand = r.Rndm();
          if (rand < 0.5) { tagind = ind[1]; probeind = ind[0]; }
          else  { tagind = ind[0];  probeind = ind[1];}
          
       }
       else {
        if (abs(jteta[ind[0]]) > 1.3 and abs(jteta[ind[1]]) <= 1.3) { tagind = ind[1]; probeind = ind[0]; }
        else if (abs(jteta[ind[1]]) > 1.3 and abs(jteta[ind[0]]) <= 1.3)  { tagind = ind[0];  probeind = ind[1]; }

        else if (abs(jteta[ind[0]]) <= 1.3 and abs(jteta[ind[1]]) <= 1.3)  {
          const auto rand = r.Rndm();
          if (rand < 0.5) { tagind = ind[1]; probeind = ind[0]; }
          else  { tagind = ind[0];  probeind = ind[1];}
        }
       }

       if (tagind < 0 || probeind < 0) continue;

       tagpt = jtpt[tagind];
       probept = jtpt[probeind];

       tageta = jteta[tagind];
       probeeta = jteta[probeind];

       if (isMC) {
         tagpt_gen = jtpt_gen[tagind];
         probept_gen = jtpt_gen[probeind];

         tageta_gen = jteta_gen[tagind];
         probeeta_gen = jteta_gen[probeind];
       }

       float dphitp = DPhi(jtphi[tagind],jtphi[probeind]);

       ptavgtp = 0.5*(tagpt  + probept);
       asymmtp = probept - tagpt;

       if (isMC) {
         ptavgtp_gen = 0.5*(tagpt_gen  + probept_gen);
         asymmtp_gen = probept_gen - tagpt_gen;
       }
       
       
       if (ind[2] != -1) alpha = jtpt[ind[2]]/ptavgtp;
       else alpha = 0; // In case only two jets

      if (ind[2] != -1 && jtpt[ind[2]] < jtptlimitforalpha) continue;
     
       for (auto &histrange : _histos) { 
	   for (auto &h : histrange.second) {
	     
	     if (tageta >= h->etamin and tageta < h->etamax and probeeta >= h->etamin and probeeta < h->etamax and hiBin >= h->hibinmin and hiBin < h->hibinmax and dphitp > 2.7  and nref >= 2 and tagind > -1) {
	       
	       // This is the full eta range
	       if ((h->etamin - h->etamax) < -10) {
		 h->alphas->Fill(alpha,evtwt);
		 
		 h->probe_pt->Fill(probept,evtwt);
		 h->probe_eta->Fill(probeeta,evtwt);
		 h->tag_pt->Fill(tagpt,evtwt);
		 h->tag_eta->Fill(tageta,evtwt);
		 
		 if (HLT_ZB) h->HLTZB_ptav->Fill(ptavgtp, evtwt);
		 if (HLT_40) h->HLT40_ptav->Fill(ptavgtp, evtwt);
		 if (HLT_60) h->HLT60_ptav->Fill(ptavgtp, evtwt);
		 if (HLT_100) h->HLT100_ptav->Fill(ptavgtp, evtwt);
		 if (HLT_120) h->HLT120_ptav->Fill(ptavgtp, evtwt);
		   
		 h->asymmdist3D->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		 h->absasymmdist3D->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);

		 
		 int probebin = h->asymmdist3D_a10->FindBin(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp);
		 int tagbin = h->asymmdist3D_a10->FindBin(ptavgtp, abs(tageta), asymmtp/2./ptavgtp);
		 //		 cout << " probebin " << probebin << " tagbin " << tagbin << " " << (probebin==tagbin) << endl;
		 //              cout << abs(probeeta) << " " << abs(tageta) << endl;
		 bool tagincorrecteta = (probebin == tagbin);    // For JER we want tag and probe in the same eta bin
		 
		 if (alpha < 0.1)   {
		   if (tagincorrecteta) h->asymmdist3D_a10->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a10->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a10) h->absasymmdist3D_gen_a10->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);
		   
		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.1-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.1-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.1-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.1-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.1-0.0001, asymmtp/2./ptavgtp, evtwt);
 
		   h->dijetasymmetry2D_a01->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		 }

		 if (alpha < 0.15) {
		   if (tagincorrecteta) h->asymmdist3D_a15->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a15->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a15) h->absasymmdist3D_gen_a15->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.15-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.15-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.15-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.15-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.15-0.0001, asymmtp/2./ptavgtp, evtwt);
}
		 if (alpha < 0.2)  {
		   if (tagincorrecteta) h->asymmdist3D_a20->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a20->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a20) h->absasymmdist3D_gen_a20->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.2-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.2-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.2-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.2-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.2-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry2D_a02->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		 }

		 if (alpha < 0.25) {
		   if (tagincorrecteta) h->asymmdist3D_a25->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a25->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a25) h->absasymmdist3D_gen_a25->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.25-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.25-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.25-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.25-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.25-0.0001, asymmtp/2./ptavgtp, evtwt);

		 }
		 if (alpha < 0.3)  {
		   if (tagincorrecteta) h->asymmdist3D_a30->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a30->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a30) h->absasymmdist3D_gen_a30->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetbalance_a03->Fill(asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry2D_a03->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry_a03->Fill(ptavgtp, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.3-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.3-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.3-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.3-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.3-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->jet_nef->Fill(probept,jtnef[probeind],evtwt);
		   h->jet_cef->Fill(probept,jtcef[probeind],evtwt);
		   h->jet_nhf->Fill(probept,jtnhf[probeind],evtwt);
		   h->jet_chf->Fill(probept,jtchf[probeind],evtwt);
		   h->jet_muf->Fill(probept,jtmuf[probeind],evtwt);       

		 }
		 //else if (alpha < 0.35)   h->dijetasymmetry2D_a035->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		 if (alpha < 0.35) {
		   if (tagincorrecteta) h->asymmdist3D_a35->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a35->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a35) h->absasymmdist3D_gen_a35->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.35-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.35-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.35-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.35-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.35-0.0001, asymmtp/2./ptavgtp, evtwt);

		 }
		 if (alpha < 0.4)   {
		   if (tagincorrecteta) h->asymmdist3D_a40->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a40->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a40) h->absasymmdist3D_gen_a40->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry2D_a04->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.4-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.4-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.4-0.0001, asymmtp/2./ptavgtp, evtwt); 
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.4-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.4-0.0001, asymmtp/2./ptavgtp, evtwt);
		 }

		 if (alpha < 0.45) {
		   if (tagincorrecteta) h->asymmdist3D_a45->Fill(ptavgtp, abs(probeeta), asymmtp/2./ptavgtp, evtwt);
		   if (tagincorrecteta) h->absasymmdist3D_a45->Fill(ptavgtp, abs(probeeta), abs(asymmtp/2./ptavgtp), evtwt);
       if (isMC && tagincorrecteta && h->absasymmdist3D_gen_a45) h->absasymmdist3D_gen_a45->Fill(ptavgtp_gen, abs(probeeta_gen), abs(asymmtp_gen/2./ptavgtp_gen), evtwt);

		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.45-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.45-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.45-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.45-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.45-0.0001, asymmtp/2./ptavgtp, evtwt);
		 }
		 
		 if (alpha < 0.5) {
		   h->dijetasymmetry2D_a05->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.5-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.5-0.0001, asymmtp/2./ptavgtp, evtwt);

		   h->dijetasymmetry3Dnarrow->Fill(ptavgtp, probeeta, 0.5-0.0001, asymmtp/2./ptavgtp, evtwt); 
		   h->dijetasymmetry3Dabsetanarrow->Fill(ptavgtp, abs(probeeta), 0.5-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.5-0.0001, asymmtp/2./ptavgtp, evtwt);

		 }
		 if (alpha < 0.6) {
		   //		   h->dijetasymmetry2D_a06->Fill(ptavgtp, probeeta, asymmtp/2./ptavgtp, evtwt);
		   //		   h->dijetasymmetry3D->Fill(ptavgtp, probeeta, 0.6-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabseta->Fill(ptavgtp, abs(probeeta), 0.6-0.0001, asymmtp/2./ptavgtp, evtwt);
		   h->dijetasymmetry3Dabsetawide->Fill(ptavgtp, abs(probeeta), 0.6-0.0001, asymmtp/2./ptavgtp, evtwt);
		 }
	       }
		 // Second: alpha < 1
		 h->dijetbalance_a1->Fill(asymmtp/2./ptavgtp, evtwt);
		 h->dijetasymmetry_a1->Fill(ptavgtp, asymmtp/2./ptavgtp, evtwt);
	     }
	   }
	 } 
     }

	 for (int j = 0; j < nref; ++j ) {

	 for (auto &histrange : _histos) {
	   for (auto &h : histrange.second) {

	     if (jteta[j] >= h->etamin and jteta[j] < h->etamax and hiBin >= h->hibinmin and hiBin < h->hibinmax) {

	       if (checkjetid and passjetid[j] == 0) continue;
	       
	       h->jetetaphi->Fill(jteta[j],jtphi[j],weight);

	       
	       if (j == 0 and passjetid[0]) {

		 // Trigger checks; leading jet pt
		 if (HLT_ZB) h->HLTZB->Fill(jtpt[ind[0]],evtwt);

		 if (HLT_40) h->HLT40->Fill(jtpt[ind[0]],evtwt);
		 if (HLT_60) h->HLT60->Fill(jtpt[ind[0]],evtwt); //  if (jtpt[0] > 100.) cout << jtpt[0] << "  " << HLT_100 << endl; }
		 if (HLT_80) h->HLT80->Fill(jtpt[ind[0]],evtwt);
		 if (HLT_100) h->HLT100->Fill(jtpt[ind[0]],evtwt);
		 if (HLT_120) { h->HLT120->Fill(jtpt[ind[0]],evtwt);}

	       }

	       h->jet_pt->Fill(jtpt[j],evtwt);
	       h->jet_pt_now->Fill(jtpt[j],1);
	       h->jet_uncorr_pt->Fill(jtpt_uncorr[j],evtwt);
	       h->jet_pt_genweight->Fill(jtpt[j],weight);
	       h->jet_eta->Fill(jteta[j],evtwt);
	       h->jet_phi->Fill(jtphi[j],evtwt);

	       if (isMC) {
	     
		 h->genjet_pt->Fill(jtpt_gen[j],evtwt);
		 h->genjet_eta->Fill(jteta_gen[j],evtwt);
		 h->genjet_phi->Fill(jtphi_gen[j],evtwt);
		 
		 h->jetresponse->Fill(jtpt_gen[j],jtpt[j]/jtpt_gen[j],evtwt);
	     
		 h->ptres->Fill((jtpt[j]-jtpt_gen[j])/jtpt_gen[j],evtwt);

		 h->ptgenvsptreco->Fill(jtpt_gen[j],jtpt[j],evtwt);
		 h->ptrecovsweight->Fill(jtpt[j],weight);
		 h->ptgenvsweight->Fill(jtpt_gen[j],weight);

		 h->responses3D->Fill(jtpt_gen[j], jteta_gen[j], jtpt[j]/jtpt_gen[j], evtwt);
		 h->phiresponse->Fill(jtpt_gen[j], jteta_gen[j], jtphi[j]-jtphi_gen[j], evtwt);
		 h->etaresponse->Fill(jtpt_gen[j], jteta_gen[j], jteta[j]-jteta_gen[j], evtwt);
		 
	       }
	     }
	   }
	 }
     }    
  }

  // Write output histograms
  
  for (auto &histrange : _histos) {
    for (auto &h : histrange.second) {
      h->Write();
    }  
  }
  eh->Write();
 
  cout << "Wrote " << outputfilename.c_str() << endl;

 }
