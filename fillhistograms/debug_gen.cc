#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

void debug_gen() {
    TFile* f = TFile::Open("root://cms-xrd-global.cern.ch//store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root");
    
    int testEntry = 10000;

    // 1. Check muonAnalyzer
    TTree* t1 = (TTree*)f->Get("muonAnalyzer/MuonTree");
    std::vector<int>* genPID1 = 0;
    if(t1) {
        t1->SetBranchAddress("genPID", &genPID1);
        t1->GetEntry(testEntry);
        std::cout << "muonAnalyzer/MuonTree: genPID size = " << (genPID1 ? genPID1->size() : -1) << std::endl;
    }

    // 2. Check ggHiNtuplizer
    TTree* t2 = (TTree*)f->Get("ggHiNtuplizer/EventTree");
    std::vector<int>* mcPID = 0;
    if(t2) {
        t2->SetBranchAddress("mcPID", &mcPID);
        t2->GetEntry(testEntry);
        std::cout << "ggHiNtuplizer/EventTree: mcPID size = " << (mcPID ? mcPID->size() : -1) << std::endl;
    }

    // 3. Check HiGenParticleAna
    TTree* t3 = (TTree*)f->Get("HiGenParticleAna/hi");
    std::vector<int>* genPID3 = 0;
    if(t3) {
        // Note: variable name might be 'mcPID' or 'pid' in this tree
        t3->SetBranchAddress("mcPID", &genPID3); 
        t3->GetEntry(testEntry);
        std::cout << "HiGenParticleAna/hi: mcPID size = " << (genPID3 ? genPID3->size() : -1) << std::endl;
    }
    
    // 4. Check Jet Tree to see if jets exist at Entry 10000
    TTree* tj = (TTree*)f->Get("ak4PFJetAnalyzer/t");
    Int_t nref;
    tj->SetBranchAddress("nref", &nref);
    tj->GetEntry(testEntry);
    std::cout << "ak4PFJetAnalyzer/t: nref = " << nref << std::endl;
}
