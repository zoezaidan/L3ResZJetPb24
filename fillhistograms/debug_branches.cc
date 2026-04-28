#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

using namespace std;

void debug_branches() {
    const char* filepath = "root://cms-xrd-global.cern.ch//store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root";
    
    cout << "Opening file..." << endl;
    TFile* f = TFile::Open(filepath, "READ");
    if (!f || f->IsZombie()) {
        cout << "Failed to open file." << endl;
        return;
    }

    // ==========================================
    // 1. EVENT TREE
    // ==========================================
    cout << "\n--- Testing hiEvtAnalyzer/HiTree ---" << endl;
    TTree* evtTree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    evtTree->SetBranchStatus("*", 0);

    Int_t hiBin = -999; evtTree->SetBranchStatus("hiBin", 1); evtTree->SetBranchAddress("hiBin", &hiBin);
    Float_t vz = -999; evtTree->SetBranchStatus("vz", 1); evtTree->SetBranchAddress("vz", &vz);
    Float_t weight = -999; evtTree->SetBranchStatus("weight", 1); evtTree->SetBranchAddress("weight", &weight);
    Float_t pthat = -999; evtTree->SetBranchStatus("pthat", 1); evtTree->SetBranchAddress("pthat", &pthat);

    cout << "Reading event tree entry 0..." << endl;
    evtTree->GetEntry(0);
    cout << "  hiBin  = " << hiBin << endl;
    cout << "  vz     = " << vz << endl;
    cout << "  weight = " << weight << endl;
    cout << "  pthat  = " << pthat << endl;

    // ==========================================
    // 2. JET TREE
    // ==========================================
    cout << "\n--- Testing ak4PFJetAnalyzer/t ---" << endl;
    TTree* jetTree = (TTree*)f->Get("ak4PFJetAnalyzer/t");
    jetTree->SetBranchStatus("*", 0);

    Int_t nref = -999; jetTree->SetBranchStatus("nref", 1); jetTree->SetBranchAddress("nref", &nref);
    const int MAXJETS = 1000; 
    Float_t rawpt[MAXJETS]; jetTree->SetBranchStatus("rawpt", 1); jetTree->SetBranchAddress("rawpt", rawpt);
    Float_t jteta[MAXJETS]; jetTree->SetBranchStatus("jteta", 1); jetTree->SetBranchAddress("jteta", jteta);
    Float_t jtPfNHF[MAXJETS]; jetTree->SetBranchStatus("jtPfNHF", 1); jetTree->SetBranchAddress("jtPfNHF", jtPfNHF);
    Int_t jtchm[MAXJETS]; jetTree->SetBranchStatus("jtPfCHM", 1); jetTree->SetBranchAddress("jtPfCHM", jtchm);

    cout << "Reading jet tree entry 0..." << endl;
    jetTree->GetEntry(0);
    cout << "  nref = " << nref << endl;
    if (nref > 0) {
        cout << "  rawpt[0]   = " << rawpt[0] << endl;
        cout << "  jteta[0]   = " << jteta[0] << endl;
        cout << "  jtPfNHF[0] = " << jtPfNHF[0] << endl;
        cout << "  jtchm[0]   = " << jtchm[0] << endl;
    }

    // ==========================================
    // 3. MUON TREE
    // ==========================================
    cout << "\n--- Testing muonAnalyzer/MuonTree ---" << endl;
    TTree* muonTree = (TTree*)f->Get("muonAnalyzer/MuonTree");
    muonTree->SetBranchStatus("*", 0);

    Int_t nReco = -999; muonTree->SetBranchStatus("nReco", 1); muonTree->SetBranchAddress("nReco", &nReco);
    vector<float>* recoPt = nullptr; muonTree->SetBranchStatus("recoPt", 1); muonTree->SetBranchAddress("recoPt", &recoPt);
    vector<bool>* recoIDTight = nullptr; muonTree->SetBranchStatus("recoIDTight", 1); muonTree->SetBranchAddress("recoIDTight", &recoIDTight);
    
    Int_t nGen = -999; muonTree->SetBranchStatus("nGen", 1); muonTree->SetBranchAddress("nGen", &nGen);
    vector<int>* genPID = nullptr; muonTree->SetBranchStatus("genPID", 1); muonTree->SetBranchAddress("genPID", &genPID);

    cout << "Reading muon tree entry 0..." << endl;
    muonTree->GetEntry(0);
    cout << "  nReco = " << nReco << endl;
    if (recoPt) cout << "  recoPt size = " << recoPt->size() << endl;
    if (recoIDTight) cout << "  recoIDTight size = " << recoIDTight->size() << endl;

    cout << "  nGen = " << nGen << endl;
    if (genPID) cout << "  genPID size = " << genPID->size() << endl;

    cout << "\nDEBUG COMPLETE: All variables read and printed successfully." << endl;
}
