import ROOT
import sys

# Once Matthew tells you the site, replace the server name here.
# For example, if it's at LLR, it might be "polgrid4.in2p3.fr"
server = "cms-xrd-global.cern.ch" 

# The exact path Matthew provided in the chat
lfn = "/store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root"

full_path = f"root://{server}/{lfn}"

print(f"Attempting to connect to: {full_path}")
print("This might take a few seconds...\n")

# Attempt to open the file over XRootD
test_file = ROOT.TFile.Open(full_path)

# Check if the file successfully opened and isn't corrupted
if not test_file or test_file.IsZombie():
    print("[ERROR] Failed to open the file.")
    print("Reasons this might happen:")
    print(" 1. The server name is wrong.")
    print(" 2. The file genuinely does not exist at that site.")
    print(" 3. The site is currently down or blocking your access.")
    sys.exit(1)

print("[SUCCESS] File opened successfully!")

# Let's peek inside to prove we can actually read the data
print("\n--- File Contents ---")
test_file.ls()
print("---------------------\n")

# Always close your files
test_file.Close()
print("Test complete.")

import ROOT

file_path = "root://cms-xrd-global.cern.ch//store/user/mnguyen/JEC/DYto2L-2Jets_MLL-50_TuneCP5_5p36TeV_amcatnloFXFX-pythia8/ZjetMadraph/260206_124018/0000/merged_HiForestMiniAOD.root"
f = ROOT.TFile.Open(file_path)

trees = {
    "hiEvtAnalyzer/HiTree": ["hiBin", "vz", "weight", "pthat"],
    "ak4PFJetAnalyzer/t": [
        "nref", "rawpt", "jteta", "jtphi", "jtPfNHF", "jtPfCHF", 
        "jtPfNEF", "jtPfCEF", "jtPfMUF", "jtPfCHM", 
        "refpt", "refeta", "refphi", "refdrjt"
    ],
    "muonAnalyzer/MuonTree": [
        "nReco", "recoPt", "recoEta", "recoPhi", "recoCharge", "recoIDTight", 
        "nGen", "genPID", "genStatus", "genPt", "genEta", "genPhi", "genMotherID"
    ]
}

print("--- FILE INSPECTION vs C++ EXPECTATIONS ---")
for tree_path, branches in trees.items():
    t = f.Get(tree_path)
    print(f"\n======================================")
    print(f"Tree: {tree_path}")
    if not t:
        print("  -> ERROR: Tree not found!")
        continue

    for b in branches:
        br = t.GetBranch(b)
        if not br:
            print(f"  [MISSING] {b:15}")
        else:
            dtype = br.GetClassName() if br.GetClassName() else br.GetTitle()
            print(f"  [FOUND]   {b:15} : {dtype}")

    # Explicitly check event 0 limits to see if we overflow MAXJETS
    if tree_path == "ak4PFJetAnalyzer/t":
        try:
            t.GetEntry(0)
            print(f"\n  [EVENT 0 CHECK] nref = {t.nref}")
        except Exception as e:
            print(f"\n  [EVENT 0 CHECK] CRASHED while reading: {e}")
