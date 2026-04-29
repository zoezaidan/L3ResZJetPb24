# L3Residual

This page documents the active photon+jet L3 residual workflow. The split fit/export chain is:

- [fillhistograms/analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc)
- [L3Residual/deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C)
- [L3Residual/L3Res.C](../L3Residual/L3Res.C)
- [L3Residual/createL2L3ResTextFile.C](../L3Residual/createL2L3ResTextFile.C)

[L3Residual/dofits_L3.C](../L3Residual/dofits_L3.C) is now an optional thin wrapper for the split chain above; the active implementation lives in the split macros themselves.

The repository also contains a Z+jet histogram filler and a shared enum-based
derivation helper for mixed photon/Z workflows:

- [fillhistograms/analyse_ZJet.cc](../fillhistograms/analyse_ZJet.cc)
- [L3Residual/deriveL3.C](../L3Residual/deriveL3.C)

## At a glance

| Stage | Macro | Required inputs | Main outputs |
| --- | --- | --- | --- |
| Histogram filling | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc) | photon+jet trees and event filters | `photonjet_balance3D*`, `photonjet_balance3D*_counts`, `photonjet_balance_dist` |
| Derivation | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C) | MC and data `photonjet_balance3D*` histograms | `ratio_vsptref_alphaN`, `ratio_norm_vsptref_alphaN`, `ratio_vsjetpt_alphaN`, `ratio_norm_vsjetpt_alphaN`, `L3Res_vsa*`, `balance3D_mc`, `balance3D_data` |
| Shared pTref fit | [L3Res.C](../L3Residual/L3Res.C) | one or more derived ROOT files, explicit sample list, fit windows | per-input alpha diagnostics, `kFSR`, `corr_vsptref`, `corr_vsjetpt`, combined sequential pTref fit, combined direct JetPt graph, fit ROOT file, validation plots |
| Text export | [createL2L3ResTextFile.C](../L3Residual/createL2L3ResTextFile.C) | fit ROOT file from `L3Res.C`, L2Residual text payload | local L3 text file, exported L3 text file, combined L2L3 text file, shared full pTref export plot |

## Macro signatures

### `analyse_ZJet.cc`

```cpp
void analyse_ZJet(string input = "ZJETHP",
                  string outputfiletag = "AK4_zjet",
                  bool isMC = false,
                  bool checkjetid = false,
                  string inputType = "era",
                  int maxFiles = -1,
                  int maxEvents = -1,
                  string outputDir = "",
                  int batchIndex = -1,
                  int totalBatches = 1,
                  string jetPath = "ak4PFJetAnalyzer/t",
                  AnalysisType analysisType = AnalysisType::ZJET_MUMU,
                  float jtptlimitforalpha = 15)
```

Flavor selection is now encoded in `AnalysisType`:

- `AnalysisType::ZJET_MUMU`: dimuon Z+jet
- `AnalysisType::ZJET_EE`: dielectron Z+jet
- `AnalysisType::ZJET`: combined ee+mumu Z+jet

The macro now uses the shared `log()` helper for normal status, warning, and
summary output.

### `deriveL3_from_photonjet.C`

```cpp
void deriveL3_from_photonjet(TString mcFile,
                             TString dataFile,
                             TString outfilename = "L3Residual_PhotonJet.root",
                             int refAlphaBin = 5,
                             bool useabs = true,
                             bool usewideabs = false)
```

Purpose:

- read the photon+jet 3D balance profiles from MC and data
- preserve the input pT binning directly from the histograms
- build the derived `pTref` and JetPt ratio series for every cumulative alpha bin
- write the normalized alpha series used later for `kFSR` extraction

Important notes:

- The macro requires the direct JetPt balance profiles already stored by the filler. It does not fall back to a photon-pT-derived approximation.
- `refAlphaBin` is the cumulative alpha bin used for the eta maps and for the normalized alpha series.
- Weighted photon+jet profiles must keep the ROOT-stored `TProfile3D` bin errors when collapsing over eta or pT. The macro does not approximate uncertainties with `1/sqrt(entries)`.

### `deriveL3.C`

```cpp
void deriveL3_from_photonjet(TString mcFile,
                             TString dataFile,
                             TString outfilename = "L3Residual.root",
                             int refAlphaBin = 5,
                             bool useabs = true,
                             bool usewideabs = false,
                             AnalysisType analysisType = AnalysisType::ZJET)
```

This shared derivation helper now selects photon+jet vs Z+jet through
`AnalysisType` instead of a string mode. It accepts `AnalysisType::PHOTONJET`
or any of the Z+jet enum values.

Unlike the older strict behavior, direct JetPt balance profiles in this helper
are optional diagnostic inputs. The main pTref derivation continues even when
those JetPt profiles are missing, and the JetPt diagnostic outputs are skipped
in that case.

### `L3Res.C`

```cpp
void L3Res(TString inFileL3Derived = "L3_derived.root",
           TString sampleTypesCSV = "",
           TString inputPtRangesCSV = "",
           string outfilename = "L3Res_photonjet",
           string runLabel = "2024ppRef",
           string lumiLabel = "pp 480.4 pb^{-1}",
           int refAlphaBin = 5,
           double fitAlphaMin = 0.0,
           double fitAlphaMax = 0.4,
           string outBaseDir = "L3Residual")
```

Purpose:

- run the alpha fits in each selected input independently
- build `kFSR` and the `alpha -> 0` corrected `corr_vsptref` histogram for each input
- build the corresponding `alpha -> 0` corrected `corr_vsjetpt` histogram from the direct JetPt ratio series
- combine the selected `pTref` points from all inputs into one shared fit
- combine the selected JetPt points from all inputs into one shared direct JetPt graph
- apply the template-style sequential fit chain in `pTref`
- save the fit ROOT file consumed by `createL2L3ResTextFile.C`

Sample values:

- `photonjet`
- `zjet`

The sample list is mandatory for multi-input fits. The macro does not guess labels or tags from filenames.

### `createL2L3ResTextFile.C`

```cpp
void createL2L3ResTextFile(TString fitRootFile = "L3Residual/L3Res_photonjet/L3Res_photonjet_fit.root",
                           string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt")
```

Purpose:

- reuse the stored full direct `pTref` final fit when available, or refit the shared `pTref` graph with the same full model for backward compatibility
- clone each eta row from the input L2Residual text payload
- append the same global direct-`pTref` L3 correction block to every cloned row
- write the standalone L3 text file and the combined L2L3 text file

## Minimal workflow

### 1. Produce photon+jet input histograms

```bash
cd fillhistograms
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_mc.txt","photonjet_mc",true,true,"filelist",-1,-1,"/output/dir")'
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_data.txt","photonjet_data",false,true,"filelist",-1,-1,"/output/dir")'
cd ..
```

### 2. Derive the L3 input histograms

```bash
root -l -b -q 'L3Residual/deriveL3_from_photonjet.C("/output/dir/filelist_mc_photonjet_mc.root","/output/dir/filelist_data_photonjet_data.root","L3Residual/L3_derived_photonjet.root",5,false,true)'
```

### 3. Fit the shared pTref response

```bash
root -l -b -q 'L3Residual/L3Res.C("L3Residual/L3_derived_photonjet.root","photonjet","60-400","L3Res_photonjet","2024ppRef","pp 480.4 pb^{-1}",5,0.0,0.4,"L3Residual")'
```

### 4. Export the text payloads

```bash
root -l -b -q 'L3Residual/createL2L3ResTextFile.C("L3Residual/L3Res_photonjet/L3Res_photonjet_fit.root","fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt")'
```

### Combined photon+jet and Z+jet fit

```bash
root -l -b -q 'L3Residual/L3Res.C("L3Residual/L3_derived_photonjet.root,L3Residual/L3_derived_zjet.root","photonjet,zjet","60-400,80-300","L3Res_combined","2024ppRef","pp 480.4 pb^{-1}",5,0.0,0.4,"L3Residual")'
root -l -b -q 'L3Residual/createL2L3ResTextFile.C("L3Residual/L3Res_combined/L3Res_combined_fit.root","fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt")'
```

## What the stages do

### Derivation

The derivation stage stores the full alpha series in both `pTref` and JetPt:

- `ratio_vsptref_alphaN`
- `ratio_norm_vsptref_alphaN`
- `ratio_vsjetpt_alphaN`
- `ratio_norm_vsjetpt_alphaN`
- `L3Res_vsa_<ptbin>_<etabin>`
- `L3Res_vsa_norm_<ptbin>_<etabin>`

All `pTref` binning comes directly from the input histogram axis. The eta output histograms use the axis definitions from the loaded balance histogram branch and the repository bin arrays in [fillhistograms/histograms.h](../fillhistograms/histograms.h).

For uncertainty handling, the derivation stage uses the ROOT-provided profile or histogram bin errors wherever those exist. The only manual step is the normalization to the reference alpha bin, because ROOT does not track the covariance between different cumulative-alpha ratio histograms. That normalization therefore uses conservative independent numerator/denominator propagation, and the reference-alpha normalized point is fixed to `1 +/- 0` by identity.

### Shared pTref fit

`L3Res.C` fits the normalized alpha series linearly in `alpha` to extract `kFSR` in each `pTref` bin, forms the corrected `corr_vsptref`, and then applies the template-style sequential fit chain:

1. constant
2. log-linear
3. log-linear + `1/x`
4. quadratic in `log10(0.01*x)`
5. quadratic + `1/x`
6. final exported reference fit: full direct `pTref` L3 model
    `1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x)`

The alpha extrapolation uses `TGraphErrors`, so the stored `ratio_norm_*` bin errors are used as fit weights. The corrected `corr_vsptref` and `corr_vsjetpt` histograms keep the reference-alpha ratio uncertainty on the final response; the fitted `kFSR` intercept uncertainty is stored separately in `kFSR` and is not multiplied back into `corr_*`, because the normalized alpha-fit points share the same reference denominator and ROOT does not model that covariance.

Only the final shared fit is shown in the main pTref summary plot. The per-input alpha fits and `kFSR` summaries are written under:

```text
<outBaseDir>/<outfilename>/inputs/<token>/pdf/
<outBaseDir>/<outfilename>/inputs/<token>/alpha_extrap/
```

### Text export

`createL2L3ResTextFile.C` keeps the L2 text payload as the reference row layout and writes one global direct-`pTref` L3Residual function for the full sample. The same fitted L3 block is repeated for every eta row so the combined file has the form

```text
L2Residual(eta, JetPt) * L3Residual(pTref-derived global fit)
```

The final text payload is anchored to the direct global `pTref` fit stored upstream by `L3Res.C`; no JetPt remap or second JetPt fit is applied in the active path. Each exported row uses the intersection of the chosen fit window and the source L2 payload validity range. The output folder contains:

- `<outfilename>.txt`: local standalone L3 payload
- `L3Residuals_<runLabel>_<tag>_AK4PF.txt`: exported standalone L3 payload
- `L2L3Residuals_<runLabel>_<tag>_AK4PF.txt`: combined L2L3 payload

The export validation plots are written in:

```text
<outBaseDir>/<outfilename>/pdf/
```

The key export figure is:

- `L3Res_<runLabel>_ptref_export_fit.png`: the shared `pTref` response points and the full direct `pTref` export fit that is written to the final text payloads.

## Raw-distribution checks

To inspect the raw input distributions before the fit, run:

```bash
root -l -b -q 'L3Residual/plotresponse_L3.C("/path/to/data_histograms.root","Data")'
root -l -b -q 'L3Residual/plotresponse_L3.C("/path/to/mc_histograms.root","MC",true)'
```

These plots are the fastest way to validate the starting point before running the derived fit chain.