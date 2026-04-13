# JER

This page documents the truth-resolution and JER scale-factor workflows. The main macros are:

- [JER/MCJER.C](../JER/MCJER.C)
- [JER/MCJPR.C](../JER/MCJPR.C)
- [JER/MCRESP.C](../JER/MCRESP.C)
- [JER/JERSF_RMS.C](../JER/JERSF_RMS.C)
- [JER/JERSF_fits.C](../JER/JERSF_fits.C)
- [JER/JERSF_fits_vsalpha.C](../JER/JERSF_fits_vsalpha.C)
- [JER/JERSF_printtxt.C](../JER/JERSF_printtxt.C)
- [JER/doTxtMCJER.C](../JER/doTxtMCJER.C)

## At a glance

| Branch | Macro sequence | Required inputs | Main outputs |
| --- | --- | --- | --- |
| MC truth resolution | `MCJER.C` or `MCJPR.C` or `MCRESP.C` then `doTxtMCJER.C` | truth-response histograms from `analyse.cc` | fitted resolution functions and text files |
| JER scale factors | `JERSF_RMS.C` or `JERSF_fits.C` then `JERSF_fits_vsalpha.C` then `JERSF_printtxt.C` | asymmetry-distribution histograms from `analyse.cc` | `SF` histogram and JER scale-factor text file |

## Runnable Macros

### 1. MC truth resolution macros

Signatures:

```cpp
void MCJER(TString inFileName,
           int minpt = 28,
           string dirname = "MCJER",
           TString outFileName = "testingMCjer-rereco.root")

void MCJPR(TString inFileName,
           bool doeta = 0,
           string dirname = "MCJPR",
           TString outFileName = "testingMCpfires-rereco.root")

void MCRESP(string inFileName,
            int minpt = 15,
            string dirname = "MCJER")
```

Parameters:

| Macro | Parameters |
| --- | --- |
| `MCJER` | input histogram file, minimum pT, output directory, output ROOT file |
| `MCJPR` | input histogram file, `doeta` switch for eta or phi resolution, output directory, output ROOT file |
| `MCRESP` | input histogram file, minimum pT, output directory |

Examples:

```bash
root -l -b -q 'JER/MCJER.C("mc_forjer.root",28,"JER","JER/MCJER.root")'
root -l -b -q 'JER/MCJPR.C("mc_forjer.root",true,"JER","JER/MCJPR_eta.root")'
root -l -b -q 'JER/MCRESP.C("mc_forjer.root",15,"JER")'
```

### 2. JER scale-factor width extraction

Signatures:

```cpp
void JERSF_RMS(TString outFileName = "JERSF_sigmas_RMS.root",
               string inFileName = "mc.root",
               string inFileNameZB = "zb.root",
               string inFileNameDT = "data.root")

void JERSF_fits(string outfilename = "JERSF_sigmas_fits.root",
                string inFileName = "mc.root",
                string inFileNameZB = "zb.root",
                string inFileNameDT = "data.root")
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `outFileName` or `outfilename` | output ROOT file with sigma histograms or graphs |
| `inFileName` | MC asymmetry histogram ROOT file |
| `inFileNameZB` | zero-bias asymmetry histogram ROOT file |
| `inFileNameDT` | data asymmetry histogram ROOT file |

Examples:

```bash
root -l -b -q 'JER/JERSF_RMS.C("JER/JERSF_sigmas_RMS.root","mc_forjer.root","zb_forjer.root","hp_forjer.root")'
root -l -b -q 'JER/JERSF_fits.C("JER/JERSF_sigmas_fits.root","mc_forjer.root","zb_forjer.root","hp_forjer.root")'
```

### 3. Alpha extrapolation and final SF export

Signatures:

```cpp
void JERSF_fits_vsalpha(string filein = "JERSF_sigmas_fits.root",
                        string outfilename = "JERSFs_fromfits.root",
                        bool RMS = false)

void JERSF_printtxt(string filein = "JERSFs_fromRMS.root",
                    string filename = "JERSF.txt")
```

Parameters:

| Parameter | Meaning |
| --- | --- |
| `filein` | input ROOT file from the Gaussian or RMS sigma stage |
| `outfilename` | output ROOT file containing the `SF` histogram |
| `RMS` | if `true`, interpret the inputs as RMS-based widths |
| `filename` | output JER scale-factor text file |

Examples:

```bash
root -l -b -q 'JER/JERSF_fits_vsalpha.C("JER/JERSF_sigmas_fits.root","JER/JERSFs_fromfits.root",false)'
root -l -b -q 'JER/JERSF_printtxt.C("JER/JERSFs_fromfits.root","JER/JERSF_fromfits.txt")'
```

## Essential Histograms

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| `responses3D` | reco/gen jet response distribution vs pT and eta | `MCJER.C`, `MCRESP.C` |
| `etaresponse` | eta residual distribution | `MCJPR.C` with `doeta=true` |
| `phiresponse` | phi residual distribution | `MCJPR.C` with `doeta=false` |
| `asymmdist3D_a*` | asymmetry distributions per alpha cut | `JERSF_fits.C` |
| `absasymmdist3D_a*` | absolute asymmetry distributions per alpha cut | `JERSF_RMS.C` |

## Minimal Chain Detail

### 1. MC truth resolution

`MCJER.C` fits the reco/gen jet response widths in bins of pT and eta. `MCJPR.C` fits eta or phi residual widths. `MCRESP.C` extracts the mean response only.

Fit forms used later by the text writer:

| Case | Expression |
| --- | --- |
| jet pT resolution | `sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])` |
| eta or phi residual resolution | `sqrt(pow([0],2)+pow([1],2)/x+pow([2]/x,2)+pow([3]/x,3))` |

### 2. JER scale factors from asymmetry widths

The scale-factor branch starts from the asymmetry histograms and ends with one scale factor per eta bin.

| Stage | Main output |
| --- | --- |
| `JERSF_fits.C` | sigma graphs from Gaussian widths |
| `JERSF_RMS.C` | sigma graphs from truncated RMS widths |
| `JERSF_fits_vsalpha.C` | `ratio_<etaBin>` and final `SF` histogram |
| `JERSF_printtxt.C` | JetMET scale-factor text file |

`JERSF_fits_vsalpha.C` does the following:

1. fits data and MC widths vs alpha with `pol1`
2. takes the alpha intercept as the `alpha -> 0` width
3. forms the data/MC ratio vs pT
4. fits that ratio with `pol0`
5. stores the constant in the `SF` histogram

## Text Outputs

### 1. MC resolution text from `doTxtMCJER.C`

Headers:

```text
{1 JetEta 1 JetPt (sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])) Resolution}
```

or, for eta or phi residuals,

```text
{1 JetEta 1 JetPt (sqrt(pow([0],2)+pow([1],2)/x+pow([2]/x,2)+pow([3]/x,3))) Resolution}
```

Each row contains eta range, parameter count, pT validity range, and the fitted parameters.

### 2. JER scale-factor text from `JERSF_printtxt.C`

Header:

```text
{ 2 JetEta JetPt 0 None ScaleFactor }
```

Each row contains eta range, pT range, payload size `3`, and the same scale factor repeated three times.

## Final Files Worth Tracking

| Product | Producer |
| --- | --- |
| MC resolution ROOT files | `MCJER.C`, `MCJPR.C`, `MCRESP.C` |
| MC resolution text file | `doTxtMCJER.C` |
| sigma ROOT file | `JERSF_RMS.C` or `JERSF_fits.C` |
| final SF ROOT file with `SF` | `JERSF_fits_vsalpha.C` |
| JER scale-factor text file | `JERSF_printtxt.C` |