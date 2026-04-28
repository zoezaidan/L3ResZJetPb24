#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

void inspect_particles() {
    TFile* f = TFile::Open("root://cms-xrd-global.cern.ch//store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root");
    TTree* t = (TTree*)f->Get("ggHiNtuplizer/EventTree");
    TTree* tj = (TTree*)f->Get("ak4PFJetAnalyzer/t");

    std::vector<int>* mcPID = 0;
    std::vector<float>* mcPt = 0;
    Int_t nref;
    Float_t rawpt[1000];

    t->SetBranchAddress("mcPID", &mcPID);
    t->SetBranchAddress("mcPt", &mcPt);
    tj->SetBranchAddress("nref", &nref);
    tj->SetBranchAddress("rawpt", rawpt);

    for (int i = 0; i < 10; ++i) {
        t->GetEntry(i);
        tj->GetEntry(i);
        std::cout << "--- Event " << i << " ---" << std::endl;
        std::cout << "  Gen Particles (mcPID): ";
        if (mcPID) {
            for (size_t j = 0; j < mcPID->size(); ++j) 
                std::cout << (*mcPID)[j] << "(" << (*mcPt)[j] << " GeV) ";
        }
        std::cout << "\n  Jets: nref=" << nref << " Leading Pt=" << (nref > 0 ? rawpt[0] : 0) << std::endl;
    }
}
