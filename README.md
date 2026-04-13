# Residual Analysis Framework

Residual-analysis workflows for L2 residuals, JER, and L3 residuals using ppRef HiForest inputs.

The repository now uses a single top-level README for setup and quick-start commands. Detailed step-by-step documentation lives under `docs/`.

## Clone and environment

Tested inside CMSSW `CMSSW_15_1_0_patch3`.

```bash
cmsrel CMSSW_15_1_0_patch3
cd CMSSW_15_1_0_patch3/src
cmsenv
```

Create a fork of this repository for your own use and development or directly use this repository.

```bash
git clone -b L3ResPhotonJet https://gitlab.cern.ch/bharikri/residualanalysis.git
cd residualanalysis
```

If ROOT is not already available in your shell and you are not using CMSSW, use the LCG fallback noted in the detailed docs.

## Build helper classes

Run this when you change histogram classes, binning, or ACLiC-built helper code:

```bash
cd fillhistograms
root -l -b -q compile.C
cd ..
```


## Quick workflows

### L2 residuals

Histogram production is driven by `fillhistograms/analyse.cc`. The derivation and fitting stages are in `L2Residual/`.

```bash
cd fillhistograms
root -l -b -q 'analyse.cc("RERECOMC","l2_mc",true,true,false,false,false,false,15,"era",-1,10000,"../test_output/L2",-1,1,"ak4PFJetAnalyzer/t")'
root -l -b -q 'analyse.cc("RERECOHP","l2_data",false,true,false,false,false,false,15,"era",-1,10000,"../test_output/L2",-1,1,"ak4PFJetAnalyzer/t")'
cd ..

root -l -b -q 'L2Residual/deriveL2_from3D.C("test_output/L2/PHOTONMC_l2_mc.root","test_output/L2/PHOTONHP_l2_data.root","test_output/L2/L2_derived.root",5,true,false)'
root -l -b -q 'L2Residual/dofits.C("test_output/L2/L2_derived.root","test_output/L2/L2_derived.root",0.15,0.35,"kfactor_test","L2fits",true)'
root -l -b -q 'L2Residual/fit_pt_param.C("test_output/L2/L2_derived.root","test_output/L2/L2_derived.root",60.,700.,"ptparam_test","ptfits",true)'
root -l -b -q 'L2Residual/doTxt.C("L2fits/kfactor_test.root","test_output/L2/L2Residual_test.txt")'
```

### JER and JER scale factors

JER depends on dijet histogram production with the asymmetry and response 3D histograms enabled.

```bash
root -l -b -q 'JER/JERSF_RMS.C("JER/JERSF_sigmas_RMS.root","mc_forjer.root","zb_forjer.root","hp_forjer.root")'
root -l -b -q 'JER/JERSF_fits.C("JER/JERSF_sigmas_fits.root","mc_forjer.root","zb_forjer.root","hp_forjer.root")'
root -l -b -q 'JER/JERSF_fits_vsalpha.C("JER/JERSF_sigmas_fits.root","JER/JERSFs_fromfits.root",false)'
root -l -b -q 'JER/JERSF_printtxt.C("JER/JERSFs_fromfits.root","JER/JERSF_fromfits.txt")'
```

### L3 residuals

Photon+jet histogram production is in `fillhistograms/analyse_PhotonJet.cc`, followed by `L3Residual/deriveL3_from_photonjet.C` and `L3Residual/dofits_L3.C`.

```bash
cd fillhistograms
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_mc.txt","photonjet_mc",true,true,"filelist",-1,-1,"/output/dir")'
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_data.txt","photonjet_data",false,true,"filelist",-1,-1,"/output/dir")'
cd ..

root -l -b -q 'L3Residual/deriveL3_from_photonjet.C("/output/dir/filelist_mc_photonjet_mc.root","/output/dir/filelist_data_photonjet_data.root","L3Residual/L3_derived_photonjet.root",true,5,false,true)'
root -l -b -q 'L3Residual/dofits_L3.C("L3Residual/L3_derived_photonjet.root","60-300","L3Res_photonjet","2024ppRef","pp 480.4 pb^{-1}",true,5,0.0,0.4,"fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt","L3Residual","")'
```

For a combined photon+jet and Z+jet fit, pass a comma-separated input list and matching fit windows:

```bash
root -l -b -q 'L3Residual/dofits_L3.C("L3Residual/L3_derived_photonjet.root,L3Residual/L3_derived_zjet.root","60-400,80-300","L3Res_combined","2024ppRef","pp 480.4 pb^{-1}",true,5,0.0,0.4,"fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt","L3Residual","#gamma+jet,Z+jet")'
```

## Documentation map

- `docs/L2Residual.md`: dijet histogram filling, data/MC propagation, alpha fits, eta fits, and text outputs.
- `docs/JER.md`: MC truth resolution, JER scale-factor derivation, alpha extrapolation, and text exports.
- `docs/L3Residual.md`: photon+jet cuts, histogram contracts, derivation, kFSR handling, and final L3 text output.
- `docs/Systematics.md`: L2 systematic-uncertainty production and text-file exports.
- `docs/Batch.md`: batch inputs, submission helpers, and merge flow.
- `docs/residualanalysis.wiki/residualanalysis.md`: GitLab wiki landing page linking the same material.

## Wiki Sync

The GitLab wiki is a separate git repository. Clone it wherever you want, then point the sync script at that checkout.

Example setup:

```bash
git clone <your-project-url>.wiki.git /path/to/residualanalysis.wiki
python3 docs/sync_wiki.py \
	--wiki-dir /path/to/residualanalysis.wiki \
	--project-url <your-project-url>
cd /path/to/residualanalysis.wiki
git status
git add .
git commit -m "Sync wiki from docs"
git push
```

If the wiki checkout lives in `docs/residualanalysis.wiki`, `python3 docs/sync_wiki.py` is enough.

## Batch processing

The submission helpers stay in `batch/`, but the user-facing instructions are consolidated in `docs/Batch.md`.

## Notes

- `triggerstudy/plottriggereff.C` contains the trigger-efficiency plotting utility.
- Existing example outputs in `L2fits/`, `L3Residual/`, and `test_output/` are not part of the documentation flow.
