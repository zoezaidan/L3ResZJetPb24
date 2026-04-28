import ROOT

file_path = "root://cms-xrd-global.cern.ch//store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root"
f = ROOT.TFile.Open(file_path)
t = f.Get("ggHiNtuplizer/EventTree")

print("--- Electron Branches in ggHiNtuplizer ---")
for b in t.GetListOfBranches():
    name = b.GetName()
    # Filter for anything that might be an ID
    if "ele" in name.lower() and ("id" in name.lower() or "bit" in name.lower() or "mva" in name.lower() or "cut" in name.lower()):
        # Print the branch name and its type
        dtype = b.GetClassName() if b.GetClassName() else b.GetTitle()
        print(f"{name:25} : {dtype}")
