# JER

This page keeps the essential JER chain: which histograms matter, which macro consumes them, and how the final MC resolution and scale-factor text files are structured.

## At a glance

| Branch | Macro sequence | Essential inputs | Essential outputs |
| --- | --- | --- | --- |
| MC truth resolution | [MCJER.C](../JER/MCJER.C#L5-154) or [MCJPR.C](../JER/MCJPR.C#L8-150) -> [doTxtMCJER.C](../JER/doTxtMCJER.C#L12-52) | `responses3D`, `etaresponse`, `phiresponse` | `fit_eta_*` functions and resolution text files |
| JER scale factors | [JERSF_fits.C](../JER/JERSF_fits.C#L1-188) or [JERSF_RMS.C](../JER/JERSF_RMS.C#L1-236) -> [JERSF_fits_vsalpha.C](../JER/JERSF_fits_vsalpha.C#L5-149) -> [JERSF_printtxt.C](../JER/JERSF_printtxt.C#L1-34) | `asymmdist3D_a*` or `absasymmdist3D_a*` | `SF` histogram and scale-factor text file |

## Essential histograms

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| [`responses3D`](../fillhistograms/analyse.cc#L802) | reco/gen jet response distribution vs pT and eta | [MCJER.C](../JER/MCJER.C#L5-154), [MCRESP.C](../JER/MCRESP.C#L1-85) |
| [`etaresponse`](../fillhistograms/analyse.cc#L804) | `eta_reco - eta_gen` distribution | [MCJPR.C](../JER/MCJPR.C#L8-150) |
| [`phiresponse`](../fillhistograms/analyse.cc#L803) | `phi_reco - phi_gen` distribution | [MCJPR.C](../JER/MCJPR.C#L8-150) |
| [`asymmdist3D_a10` ... `a45`](../fillhistograms/analyse.cc#L613) | asymmetry distributions per alpha cut | [JERSF_fits.C](../JER/JERSF_fits.C#L1-188) |
| [`absasymmdist3D_a10` ... `a45`](../fillhistograms/analyse.cc#L614) | absolute asymmetry distributions per alpha cut | [JERSF_RMS.C](../JER/JERSF_RMS.C#L1-236) |

## Minimal chain detail

### 1. MC truth resolution

| Macro | What it extracts | Output object |
| --- | --- | --- |
| [MCJER.C](../JER/MCJER.C#L5-154) | Gaussian width of reco/gen response in each `(pT_gen, eta)` bin | `widths_<etaBin>`, `fit_eta_<etaRange>` |
| [MCJPR.C](../JER/MCJPR.C#L8-150) | Gaussian width of eta or phi residuals | `widths_<etaBin>`, `fit_eta_<etaRange>` |
| [MCRESP.C](../JER/MCRESP.C#L1-85) | mean response only | `widths_<etaBin>` as response means |

Fit forms used for the text export:

| Case | Expression |
| --- | --- |
| jet pT resolution | `sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])` |
| eta or phi resolution | `sqrt(pow([0],2)+pow([1],2)/x+pow([2]/x,2)+pow([3]/x,3))` |

### 2. JER scale factors from asymmetry widths

The scale-factor branch starts from the asymmetry histograms and ends with one scale factor per eta bin.

| Stage | Key output |
| --- | --- |
| [JERSF_fits.C](../JER/JERSF_fits.C#L1-188) | `gsigmasMC*`, `gsigmasDT*` from Gaussian widths |
| [JERSF_RMS.C](../JER/JERSF_RMS.C#L1-236) | `gsigmasMC*`, `gsigmasDT*` from truncated RMS |
| [JERSF_fits_vsalpha.C](../JER/JERSF_fits_vsalpha.C#L5-149) | `ratio_<etaBin>` and final `SF` histogram |
| [JERSF_printtxt.C](../JER/JERSF_printtxt.C#L1-34) | JetMET scale-factor text file |

The final SF extraction in [JERSF_fits_vsalpha.C](../JER/JERSF_fits_vsalpha.C#L96-149) is:

1. fit data and MC widths vs alpha with `pol1`
2. take the alpha intercept as the `alpha -> 0` width
3. build `Data/MC` vs pT
4. fit that ratio with `pol0`
5. store the constant in the `SF` histogram

## Text outputs

### 1. MC resolution text from `doTxtMCJER.C`

Headers:

```text
{1 JetEta 1 JetPt (sqrt([0]*[0]/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])) Resolution}
```

or, for eta/phi residuals,

```text
{1 JetEta 1 JetPt (sqrt(pow([0],2)+pow([1],2)/x+pow([2]/x,2)+pow([3]/x,3))) Resolution}
```

Row layout written by [doTxtMCJER.C](../JER/doTxtMCJER.C#L12-52):

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3 | number of parameters |
| 4-5 | pT min, pT max |
| 6-9 | `p0 p1 p2 p3` |

### 2. JER scale-factor text from `JERSF_printtxt.C`

Header:

```text
{ 2 JetEta JetPt 0 None ScaleFactor }
```

Row layout written by [JERSF_printtxt.C](../JER/JERSF_printtxt.C#L1-34):

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3-4 | pT min, pT max |
| 5 | number of payload entries, fixed to `3` |
| 6-8 | scale factor repeated three times |

The file is mirrored to negative eta and positive eta explicitly.

## Final files worth tracking

| Product | Produced by |
| --- | --- |
| MC resolution ROOT file with `fit_eta_*` | [MCJER.C](../JER/MCJER.C#L5-154) or [MCJPR.C](../JER/MCJPR.C#L8-150) |
| MC resolution text file | [doTxtMCJER.C](../JER/doTxtMCJER.C#L12-52) |
| sigma ROOT file from Gaussian or RMS branch | [JERSF_fits.C](../JER/JERSF_fits.C#L1-188) or [JERSF_RMS.C](../JER/JERSF_RMS.C#L1-236) |
| final SF ROOT file containing `SF` | [JERSF_fits_vsalpha.C](../JER/JERSF_fits_vsalpha.C#L5-149) |
| JER scale-factor text file | [JERSF_printtxt.C](../JER/JERSF_printtxt.C#L1-34) |