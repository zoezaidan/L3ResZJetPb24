# L3Residual

This page documents the photon+jet L3 residual workflow from histogram filling to the final text payloads. The source macros are:

- [fillhistograms/analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc)
- [L3Residual/deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C)
- [L3Residual/dofits_L3.C](../L3Residual/dofits_L3.C)

## At a glance

| Stage | Macro | Required inputs | Main outputs |
| --- | --- | --- | --- |
| Histogram filling | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc) | photon tree, jet tree, event filters, alpha thresholds | `photonjet_balance3D*`, `photonjet_balance3D*_counts`, `photonjet_balance_dist` |
| Derivation | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C) | MC and data `photonjet_balance3D*` | `ratio_vsphotonpt_alphaN`, `ratio_norm_vsphotonpt_alphaN`, `L3Res_vsa*`, `balance3D_mc`, `balance3D_data` |
| Final fit | [dofits_L3.C](../L3Residual/dofits_L3.C) | one or more derived ROOT files, optional L2 text file | per-input alpha diagnostics, final L3 fit, L3 text file, optional combined L2L3 text file |

## Runnable Macros

### 1. Histogram production with `analyse_PhotonJet.cc`

Signature:

```cpp
void analyse_PhotonJet(string input = "PHOTONHP",
                       string outputfiletag = "AK4_photonjet",
                       bool isMC = false,
                       bool checkjetid = false,
                       string inputType = "era",
                       int maxFiles = -1,
                       int maxEvents = -1,
                       string outputDir = "",
                       int batchIndex = -1,
                       int totalBatches = 1,
                       string jetPath = "ak4PFJetAnalyzer/t")
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `input` | era key, file path, directory path, or filelist path depending on `inputType` |
| `outputfiletag` | tag appended to the produced ROOT filename |
| `isMC` | switches MC-only branches and pThat weighting |
| `checkjetid` | enables the jet ID requirement |
| `inputType` | one of `era`, `file`, `directory`, or `filelist` |
| `maxFiles` | maximum number of input files to read, `-1` means all |
| `maxEvents` | maximum number of events to process, `-1` means all |
| `outputDir` | directory where the histogram ROOT file is written |
| `batchIndex` | current batch slot for Condor splitting |
| `totalBatches` | total number of batch slots |
| `jetPath` | jet tree path such as `ak4PFJetAnalyzer/t` or `ak4PFJetAnalyzerSDZcut1/t` |

Minimal examples:

```bash
cd fillhistograms
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_mc.txt","photonjet_mc",true,true,"filelist",-1,-1,"/output/dir")'
root -l -b -q 'analyse_PhotonJet.cc("/path/to/filelist_data.txt","photonjet_data",false,true,"filelist",-1,-1,"/output/dir")'
cd ..
```

### 2. Derivation with `deriveL3_from_photonjet.C`

Signature:

```cpp
void deriveL3_from_photonjet(TString mcFile,
                             TString dataFile,
                             TString outfilename = "L3Residual_PhotonJet.root",
                             bool dodt = true,
                             int alphabin = 5,
                             bool useabs = true,
                             bool usewideabs = false)
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `mcFile` | photon+jet histogram ROOT file from MC |
| `dataFile` | photon+jet histogram ROOT file from data |
| `outfilename` | output derived ROOT file |
| `dodt` | retained for compatibility, not used to branch the current derivation logic |
| `alphabin` | reference cumulative alpha bin used for the eta maps |
| `useabs` | use the absolute-eta histogram branch |
| `usewideabs` | use the wide absolute-eta histogram branch; this overrides `useabs` |

Example:

```bash
root -l -b -q 'L3Residual/deriveL3_from_photonjet.C("/output/dir/filelist_mc_photonjet_mc.root","/output/dir/filelist_data_photonjet_data.root","L3Residual/L3_derived_photonjet.root",true,5,false,true)'
```

### 3. Final fit with `dofits_L3.C`

Signature:

```cpp
void dofits_L3(TString inFileL3Derived = "L3_derived.root",
               TString inputPtRangesCSV = "",
               string outfilename = "L3Res_photonjet",
               string runLabel = "2024ppRef",
               string lumiLabel = "pp 480.4 pb^{-1}",
               bool saveAlphaExtrap = true,
               int refAlphaBin = 5,
               double fitAlphaMin = 0.0,
               double fitAlphaMax = 0.4,
               string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
               string outBaseDir = "L3Residual",
               TString inputLabelsCSV = "")
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `inFileL3Derived` | one derived ROOT file or a comma-separated list of derived ROOT files |
| `inputPtRangesCSV` | fit window per input, for example `80-300` or `60-400,80-300`; empty means use the full correction-histogram range for each input |
| `outfilename` | output tag used for the fit ROOT file and local text output |
| `runLabel` | run label used in plot names and exported text filenames |
| `lumiLabel` | luminosity label passed to the plot stamp |
| `saveAlphaExtrap` | if `true`, run or reuse per-input alpha extrapolation and save diagnostics under `inputs/<token>/alpha_extrap/` |
| `refAlphaBin` | reference cumulative alpha bin used as the denominator for alpha normalization |
| `fitAlphaMin` | lower alpha boundary for the linear alpha extrapolation |
| `fitAlphaMax` | upper alpha boundary for the linear alpha extrapolation |
| `l2ResidualFile` | optional L2 text file used to reuse eta bins and to build the combined L2L3 text file |
| `outBaseDir` | parent output directory |
| `inputLabelsCSV` | optional comma-separated labels used in the combined plot legend and per-input output folder names |

Single-input example:

```bash
root -l -b -q 'L3Residual/dofits_L3.C("L3Residual/L3_derived_photonjet.root","60-300","L3Res_photonjet","2024ppRef","pp 480.4 pb^{-1}",true,5,0.0,0.4,"fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt","L3Residual","")'
```

Combined-input example:

```bash
root -l -b -q 'L3Residual/dofits_L3.C("L3Residual/L3_derived_photonjet.root,L3Residual/L3_derived_zjet.root","60-400,80-300","L3Res_combined","2024ppRef","pp 480.4 pb^{-1}",true,5,0.0,0.4,"fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt","L3Residual","#gamma+jet,Z+jet")'
```

## Essential Histograms

### Filled in `analyse_PhotonJet.cc`

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| `photonjet_balance3D*` | balance profile vs photon pT, eta, and cumulative alpha | `deriveL3_from_photonjet.C` |
| `photonjet_balance3Dabseta*` | same payload in absolute eta | `deriveL3_from_photonjet.C` |
| `photonjet_balance3D*_counts` | event counts per `(pT, eta, alpha)` bin | `deriveL3_from_photonjet.C`, QA |
| `photonjet_balance_dist` | balance distribution vs photon pT and alpha | optional distribution overlays |

### Produced by `deriveL3_from_photonjet.C`

| Output | Meaning |
| --- | --- |
| `L3Res_vsa_<ptbin>_<etabin>` | MC/data ratio vs alpha in one pT and eta bin |
| `L3Res_vsa_norm_<ptbin>_<etabin>` | the same ratio normalized to the reference alpha bin |
| `ratio_vsphotonpt_alphaN` | MC/data ratio vs photon pT for cumulative alpha bin `N` |
| `ratio_norm_vsphotonpt_alphaN` | `ratio_vsphotonpt_alphaN / ratio_vsphotonpt_alphaRef` |
| `balance3D_mc`, `balance3D_data` | copied 3D profiles used by `dofits_L3.C` if raw alpha extrapolation is rerun |
| `counts3D_mc`, `counts3D_data` | copied 3D count histograms |

## Minimal Chain Detail

### 1. Histogram filling

The response variable propagated through the chain is

$$
\mathrm{balance} = \frac{p_{T,\mathrm{jet}}}{p_{T,\gamma}}.
$$

The second away-side jet defines

$$
\alpha = \frac{p_{T,\mathrm{2nd\ away\ jet}}}{p_{T,\gamma}}.
$$

The filler writes cumulative alpha bins. Each stored alpha bin already represents `alpha < threshold`.

### 2. Derivation

The derivation stage forms MC/data balance ratios by eta, alpha, and photon pT. The reference-alpha normalized outputs are the direct inputs for the alpha extrapolation in the final fit macro.

### 3. Final fit in `dofits_L3.C`

The fit macro runs the same flow for one input or several inputs:

1. prepare each input independently
2. derive or reuse `alpha -> 0` corrections for each input
3. concatenate the selected pT points from all inputs
4. fit one shared L3 correction function

If `saveAlphaExtrap=true`, each input gets its own output directory:

```text
<outBaseDir>/<outfilename>/inputs/<token>/pdf/
<outBaseDir>/<outfilename>/inputs/<token>/alpha_extrap/
<outBaseDir>/<outfilename>/inputs/<token>/textfiles/
```

Input preparation rules:

| Input content | Behavior |
| --- | --- |
| `ratio_vsphotonpt_alphaN` present | recompute the alpha extrapolation from the stored ratio histograms |
| archived `gAlphaNorm_pt*`, `fAlpha_pt*_col`, `kFSR_vsPt`, and `ratio_vspT_alpha0` present | reuse the archived alpha diagnostics and correction histogram |
| only `corr_vspT` or `ratio_vspT_alpha0` present | skip alpha refit and use the stored correction directly |

The final fit function is

$$
C_{L3}(x) = \frac{1}{p_0 + \frac{p_1}{x} + \frac{p_2 \log(x)}{x} + p_3 \frac{(x/p_4)^{p_5} - 1}{(x/p_4)^{p_5} + 1} + p_6 x^{-0.3051} + p_7 x},
$$

with $x = p_T^\gamma$.

The pT validity written into the final text file comes from the effective fit windows:

| Mode | Written pT validity |
| --- | --- |
| single input | the selected window for that input |
| multiple inputs | the union of the selected per-input windows used in the combined fit |

## Text Outputs

### 1. L3 residual text file

The exported header is

```text
{ 1 JetEta 1 JetPt 1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x) Correction L2Relative}
```

The correction name is kept as `L2Relative` to match the payload convention already used in this repository.

Each row is written as

```text
etaMin etaMax 10 ptMin ptMax p0 p1 p2 p3 p4 p5 p6 p7
```

Column meaning:

| Column block | Meaning |
| --- | --- |
| `etaMin etaMax` | eta validity range |
| `10` | payload size, equal to `2` range values plus `8` fit parameters |
| `ptMin ptMax` | pT validity range used by the final fit |
| `p0 ... p7` | parameters of the final L3 fit function |

If an L2 text file is supplied, the L3 writer reuses the eta binning from that L2 file. Otherwise one eta range is written using the input histogram eta support.

Files written:

| Path | Meaning |
| --- | --- |
| `<outBaseDir>/<outfilename>/textfiles/<outfilename>.txt` | local copy of the L3 result |
| `<outBaseDir>/<outfilename>/textfiles/L3Residuals_<run>_<tag>_AK4PF.txt` | JetMET-style export |

For a single default photon+jet input, `<tag>` is `photonjet`. For a combined fit, `<tag>` is `combined`.

### 2. Combined L2L3 residual text file

If the L2 text file can be read, `dofits_L3.C` also writes a combined L2L3 payload. The writer:

1. reads the original L2 formula and eta binning
2. appends the shifted L3 fit expression to that L2 formula
3. writes the original L2 parameters followed by the fitted L3 parameters

The combined header is therefore the original L2 expression multiplied by the shifted L3 expression, with the original correction name preserved.

Each row contains:

| Block | Meaning |
| --- | --- |
| eta range and pT range | copied from the L2 input record |
| original L2 parameters | copied from the L2 input record |
| appended L3 parameters | the `p0 ... p7` values from the final fit |

Written file:

```text
<outBaseDir>/<outfilename>/textfiles/L2L3Residuals_<run>_<tag>_AK4PF.txt
```

## Final Files Worth Tracking

| Product | Producer |
| --- | --- |
| photon+jet histogram ROOT file | `analyse_PhotonJet.cc` |
| derived L3 ROOT file | `deriveL3_from_photonjet.C` |
| per-input alpha diagnostic ROOT files | `dofits_L3.C` |
| final fit ROOT file | `dofits_L3.C` |
| final L3 text file | `dofits_L3.C` |
| optional combined L2L3 text file | `dofits_L3.C` |