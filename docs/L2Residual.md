# L2Residual

This page documents the dijet L2 residual workflow from histogram filling to the final text exports. The main macros are:

- [fillhistograms/analyse.cc](../fillhistograms/analyse.cc)
- [L2Residual/deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C)
- [L2Residual/dofits.C](../L2Residual/dofits.C)
- [L2Residual/fit_pt_param.C](../L2Residual/fit_pt_param.C)
- [L2Residual/doTxt.C](../L2Residual/doTxt.C)
- [L2Residual/L2res_param_txt.C](../L2Residual/L2res_param_txt.C)
- [L2Residual/L2res_Run3param_txt.C](../L2Residual/L2res_Run3param_txt.C)

## At a glance

| Stage | Macro | Required inputs | Main outputs |
| --- | --- | --- | --- |
| Histogram filling | [analyse.cc](../fillhistograms/analyse.cc) | dijet event trees, jet ID settings, alpha thresholds | `dijetasymmetry3D*`, `dijetasymmetry2D_a*`, `asymmdist3D_a*`, `absasymmdist3D_a*` |
| Derivation | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C) | MC and data `dijetasymmetry3D*` | `ratio_pt*_alpha*`, `Respvsa_*`, `Respvsa_norm_*` |
| Alpha fits | [dofits.C](../L2Residual/dofits.C) | derived response-vs-alpha histograms | `ratio`, `corrections_<ptrange>` |
| pT fits | [fit_pt_param.C](../L2Residual/fit_pt_param.C) | `corrections_<ptrange>` | `loglin_etaN`, `const_etaN`, `run3_etaN` |
| Text export | [doTxt.C](../L2Residual/doTxt.C), [L2res_param_txt.C](../L2Residual/L2res_param_txt.C), [L2res_Run3param_txt.C](../L2Residual/L2res_Run3param_txt.C) | alpha-fit or pT-fit ROOT files | JetMET-style L2 residual text files |

## Runnable Macros

### 1. Histogram production with `analyse.cc`

Signature:

```cpp
void analyse(string input = "RERECOHP",
             string outputfiletag = "AK4_nojetid",
             bool isMC = false,
             bool checkjetid = false,
             bool iszb = false,
             bool dol2res = false,
             bool dojer = false,
             bool fillforJER = false,
             float jtptlimitforalpha = 15,
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
| `input` | era key or path depending on `inputType` |
| `outputfiletag` | output ROOT filename tag |
| `isMC` | enable MC-specific branches and weights |
| `checkjetid` | apply the jet ID selection |
| `iszb` | zero-bias specific handling |
| `dol2res` | enable the L2 residual histogram set |
| `dojer` | enable the truth-response JER histograms |
| `fillforJER` | enable the asymmetry-distribution histograms used by JER scale factors |
| `jtptlimitforalpha` | threshold used when constructing the alpha observable |
| `inputType` | one of `era`, `file`, `directory`, or `filelist` |
| `maxFiles` | maximum input files, `-1` means all |
| `maxEvents` | maximum events, `-1` means all |
| `outputDir` | output directory |
| `batchIndex` | batch slot for Condor splitting |
| `totalBatches` | total number of batch slots |
| `jetPath` | jet tree path |

Example commands for the L2 branch:

```bash
cd fillhistograms
root -l -b -q 'analyse.cc("RERECOMC","l2_mc",true,true,false,true,false,false,15,"era",-1,-1,"../test_output/L2",-1,1,"ak4PFJetAnalyzer/t")'
root -l -b -q 'analyse.cc("RERECOHP","l2_data",false,true,false,true,false,false,15,"era",-1,-1,"../test_output/L2",-1,1,"ak4PFJetAnalyzer/t")'
cd ..
```

### 2. Derivation with `deriveL2_from3D.C`

Signature:

```cpp
void deriveL2_from3D(TString inFileName,
                     TString inFileNameDT,
                     TString outfilename = "RERUNALL/L2residuals.root",
                     int alphabin = 5,
                     bool useabs = true,
                     bool usewideabs = false)
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `inFileName` | MC histogram ROOT file |
| `inFileNameDT` | data histogram ROOT file |
| `outfilename` | derived output ROOT file |
| `alphabin` | reference cumulative alpha bin used for the eta maps |
| `useabs` | use absolute eta binning |
| `usewideabs` | use wide absolute eta binning |

Example:

```bash
root -l -b -q 'L2Residual/deriveL2_from3D.C("test_output/L2/PHOTONMC_l2_mc.root","test_output/L2/PHOTONHP_l2_data.root","test_output/L2/L2_derived.root",5,true,false)'
```

### 3. Alpha fits with `dofits.C`

Signature:

```cpp
void dofits(TString inzb,
            TString inHP,
            float fitmin = 0.15,
            float fitmax = 0.35,
            string outfilename = "kfactor",
            string outfolder = "L2fits",
            bool doabseta = true)
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `inzb` | first derived input file |
| `inHP` | second derived input file |
| `fitmin` | lower alpha bound for the linear fit |
| `fitmax` | upper alpha bound for the linear fit |
| `outfilename` | output ROOT file tag |
| `outfolder` | output directory |
| `doabseta` | use absolute eta bins |

Example:

```bash
root -l -b -q 'L2Residual/dofits.C("test_output/L2/L2_derived.root","test_output/L2/L2_derived.root",0.15,0.35,"kfactor_test","L2fits",true)'
```

### 4. pT fits with `fit_pt_param.C`

Signature:

```cpp
void fit_pt_param(TString inzb,
                  TString inHP,
                  float fitmin = 25.,
                  float fitmax = 1000.,
                  string outfilename = "testing_pt_dep",
                  string outfolder = "ptfits",
                  bool doabseta = true)
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `inzb` | first alpha-fit ROOT file |
| `inHP` | second alpha-fit ROOT file |
| `fitmin` | lower pT bound for the fit |
| `fitmax` | upper pT bound for the fit |
| `outfilename` | output ROOT file tag |
| `outfolder` | output directory |
| `doabseta` | use absolute eta bins |

Example:

```bash
root -l -b -q 'L2Residual/fit_pt_param.C("L2fits/kfactor_test.root","L2fits/kfactor_test.root",60.,700.,"ptparam_test","ptfits",true)'
```

### 5. Text writers

Signatures:

```cpp
void doTxt(TString inFileName = "../RERECORESULTS/kfactor.root",
           TString outFileName = "L2residual.txt")

void L2res_param_txt(TString inFileName = "ptfits/testing_pt_dep.root",
                     TString kfactFile = "L2fits/correctionfile.root",
                     TString outFileName = "L2residual-param.txt",
                     bool loglin = true,
                     bool kfsr = false)

void L2res_Run3param_txt(TString inFileName = "ptfits/testing_pt_dep.root",
                         TString kfactFile = "L2fits/correctionfile.root",
                         TString outFileName = "L2residual-param.txt",
                         bool loglin = true,
                         bool kfsr = false)
```

Parameters:

| Macro | Parameters |
| --- | --- |
| `doTxt` | alpha-fit ROOT file and output text filename |
| `L2res_param_txt` | pT-fit ROOT file, k-factor ROOT file, output text filename, `loglin` switch, `kfsr` switch |
| `L2res_Run3param_txt` | same inputs as `L2res_param_txt`, but writes the Run-3-style functional form |

Examples:

```bash
root -l -b -q 'L2Residual/doTxt.C("L2fits/kfactor_test.root","test_output/L2/L2Residual_test.txt")'
root -l -b -q 'L2Residual/L2res_param_txt.C("ptfits/ptparam_test.root","L2fits/kfactor_test.root","test_output/L2/L2Residual_param.txt",true,false)'
root -l -b -q 'L2Residual/L2res_Run3param_txt.C("ptfits/ptparam_test.root","L2fits/kfactor_test.root","test_output/L2/L2Residual_run3.txt",true,false)'
```

## Essential Histograms

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| `dijetasymmetry3D*` | main asymmetry profiles vs pT, eta, and alpha | `deriveL2_from3D.C` |
| `dijetasymmetry2D_a*` | QA profiles at fixed cumulative alpha cuts | plotting and checks |
| `asymmdist3D_a*` | asymmetry distributions for Gaussian JER fits | `JERSF_fits.C` |
| `absasymmdist3D_a*` | absolute asymmetry distributions for RMS JER fits | `JERSF_RMS.C` |

## Minimal Chain Detail

### 1. Histogram filling

The asymmetry observable is

$$
A = \frac{p_{T,\mathrm{probe}} - p_{T,\mathrm{tag}}}{2 p_{T,\mathrm{avg}}}.
$$

The histogram filler stores this cumulatively in alpha. The `dijetasymmetry3D*` family is the required input for the L2 derivation.

### 2. Derivation in `deriveL2_from3D.C`

For each `(pT, eta, alpha)` bin the macro converts asymmetry to response and forms the MC/data ratio:

$$
R = \frac{1 + A}{1 - A}, \qquad \mathrm{L2Res} = \frac{R_{MC}}{R_{Data}}.
$$

Main outputs:

| Output | Meaning |
| --- | --- |
| `ratio_pt<lo>to<hi>_alpha<cut>` | eta-dependent MC/data ratio at one cumulative alpha cut |
| `Respvsa_<ptbin>_<etabin>` | raw response ratio vs alpha |
| `Respvsa_norm_<ptbin>_<etabin>` | response ratio normalized to the reference alpha bin |

### 3. Alpha fits in `dofits.C`

The alpha-fit stage uses the normalized response-vs-alpha histograms and extracts an eta-dependent correction factor.

| Output | Meaning |
| --- | --- |
| `ratio` | final eta-dependent k-factor histogram |
| `corrections_<ptrange>` | alpha-corrected eta histogram for each pT slice |

### 4. pT fits in `fit_pt_param.C`

The pT-fit stage takes one eta bin across all `corrections_<ptrange>` histograms and builds a new histogram vs pT. It then fits that histogram in each eta bin.

Functions written to the output ROOT file:

| Name | Expression |
| --- | --- |
| `loglin_etaN` | `[0] + [1]*log(x)` |
| `const_etaN` | `pol0` |
| `run3_etaN` | `1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))` |

## Text Outputs

### 1. Binned export from `doTxt.C`

Header:

```text
{2 JetEta JetPt 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}
```

### 2. Parametrized export from `L2res_param_txt.C`

Headers:

```text
{1 JetEta 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}
```

and

```text
{1 JetEta 1 JetPt 1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))*[3] Correction L2Relative}
```

The last parameter is the eta-dependent k-factor if `kfsr=true`.

## Final Files Worth Tracking

| Product | Producer |
| --- | --- |
| histogram ROOT files | `analyse.cc` |
| derived L2 ROOT file | `deriveL2_from3D.C` |
| alpha-fit ROOT file | `dofits.C` |
| pT-fit ROOT file | `fit_pt_param.C` |
| L2 residual text files | `doTxt.C`, `L2res_param_txt.C`, `L2res_Run3param_txt.C` |