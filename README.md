# residualanalysis

Repository for deriving jet energy calibration and resolution products from HiForest ntuples and merged histogram ROOT files.

Active workflows in this repository:

- L2 residual JEC from dijet asymmetry
- L3 residual JEC from photon+jet balance
- MC truth JER fits and JER scale factors
- batch submission for histogram production

The codebase started from a 2023 ppRef dijet L2/JER workflow and now also contains the in-progress 2024 ppRef and OO L2/L3 workflows. The main correction paths are still direct-balance based.

Legacy 2023 input forests used in the cleaned-up repository lived under:

`/eos/cms/store/group/phys_heavyions/lamartik/DIJET_JEC_FORESTS/`

## Repository layout

| Path | Purpose |
| --- | --- |
| `fillhistograms/` | histogram production for dijet, photon+jet, and JER inputs |
| `L2Residual/` | dijet L2 residual derivation, alpha fits, pT fits, and text export |
| `L3Residual/` | photon+jet L3 residual derivation, shared fits, and L2L3 export |
| `JER/` | MC truth resolution and JER scale-factor derivations |
| `batch/` | HTCondor submission and output merging for histogram production |
| `docs/` | current workflow documentation |
| `triggerstudy/` | trigger turn-on studies and related plotting |

## Environment setup

You need a CMS environment with ROOT and the JetMET dependencies available. Any compatible CMSSW release can work for the histogram-filling step; the current docs use `CMSSW_15_1_0_patch3` as an example.

Typical setup:

```bash
cmsrel CMSSW_15_1_0_patch3
cd CMSSW_15_1_0_patch3/src
cmsenv
git clone <your-project-url>
cd residualanalysis
```

Compile helper classes after changing histogram definitions or shared histogram headers:

```bash
cd fillhistograms
root -l -b -q compile.C
cd ..
```

## Before you run

- `fillhistograms/settings.h` controls the JEC and JER payloads applied during histogram filling.
- `fillhistograms/jecfiles/` stores payloads used both during histogram production and in exported text files.
- The cleaned-up 2023 code was still written for a fairly narrow set of datasets. If you adapt the workflows to a new dataset, verify trigger logic, alpha and pT binning, eta selections, and downstream post-processing macros rather than assuming they are fully generic.
- In the dijet filler, the most common hard-coded switches to check are `applyjetvetomap`, `isrun3jersf`, and `usecalotrig`.
- Batch mode only fills histograms. Merge the shard outputs first, then run the downstream L2, JER, or L3 derivation and fit macros on the merged ROOT file.

## Workflow overview

### 1. Histogram production

Main entry points:

- `fillhistograms/analyse.cc` for dijet L2 residual and JER inputs
- `fillhistograms/analyse_PhotonJet.cc` for photon+jet L3 inputs
- `batch/submit_condor.py` and `batch/merge_outputs.sh` for Condor production

The local fillers support the same basic input modes used throughout the repository: `era`, `file`, `directory`, and `filelist`.

### 2. L2 residuals

Minimal chain:

1. Fill dijet histograms with `fillhistograms/analyse.cc` using `dol2res=true`.
2. Derive MC/data response ratios with `L2Residual/deriveL2_from3D.C`.
3. Fit the alpha dependence with `L2Residual/dofits.C`.
4. Optionally fit the pT dependence with `L2Residual/fit_pt_param.C`.
5. Export payloads with one of:
	- `L2Residual/doTxt.C`
	- `L2Residual/L2res_param_txt.C`
	- `L2Residual/L2res_Run3param_txt.C`

The detailed step-by-step interface, argument tables, and example commands are in `docs/L2Residual.md`.

### 3. JER and JER scale factors

For histogram filling in the dijet chain:

- use `dojer=true` for truth-response studies
- use `fillforJER=true` for asymmetry distributions used in scale-factor extraction

MC truth resolution chain:

- `JER/MCJER.C`
- `JER/MCJPR.C`
- `JER/MCRESP.C`
- `JER/doTxtMCJER.C`

JER scale-factor chain:

1. width extraction with either:
	- `JER/JERSF_fits.C`
	- `JER/JERSF_RMS.C`
2. alpha extrapolation with `JER/JERSF_fits_vsalpha.C`
3. text export with `JER/JERSF_printtxt.C`

Apply the intended L2 residual correction before deriving JER scale factors.

Detailed signatures and outputs are documented in `docs/JER.md`.

### 4. L3 residuals

Active split workflow:

`deriveL3_from_photonjet.C -> L3Res.C -> createL2L3ResTextFile.C`

Minimal chain:

1. Fill photon+jet histograms with `fillhistograms/analyse_PhotonJet.cc`.
2. Derive the response-ratio inputs with `L3Residual/deriveL3_from_photonjet.C`.
3. Run the shared pTref fit with `L3Residual/L3Res.C`.
4. Export standalone L3 and combined L2L3 text payloads with `L3Residual/createL2L3ResTextFile.C`.

For the current photon+jet path, `L3Residual/runL3RES.C` can also be used as a wrapper that runs input plotting, derivation, fitting, and text export in one go.

Detailed signatures, sample combinations, and example commands are in `docs/L3Residual.md`.

### 5. Batch workflow

Batch helpers only cover histogram production:

1. submit with `batch/submit_condor.py`
2. merge with `batch/merge_outputs.sh`
3. run the downstream L2, JER, or L3 macros on the merged ROOT file

See `docs/Batch.md` for the full submission interface and options.

## Documentation

The current workflow documentation lives in `docs/`:

- `docs/WikiHome.md`
- `docs/L2Residual.md`
- `docs/JER.md`
- `docs/L3Residual.md`
- `docs/Systematics.md`
- `docs/Batch.md`

## Notes

- Several older 2023 payloads and helper files are still kept for comparison, legacy reprocessing, or cross-checks.
- If you change histogram definitions or branch handling, rerun `fillhistograms/compile.C` before launching new production.
