# L2Residual

This page keeps only the essential L2 chain: which histograms are filled, which outputs are passed to the next macro, and how the final text files are structured.

## At a glance

| Stage | Macro | Essential input histograms or objects | Essential output passed downstream |
| --- | --- | --- | --- |
| Histogram filling | [analyse.cc](../fillhistograms/analyse.cc#L314-807) | event trees, jets, alpha thresholds, dijet selection | `dijetasymmetry3D*`, `dijetasymmetry2D_a*`, `asymmdist3D_a*`, `absasymmdist3D_a*` |
| Derivation | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L1-259) | `dijetasymmetry3D*` from MC and data | `mc_pt*_alpha*`, `dt_pt*_alpha*`, `ratio_pt*_alpha*`, `Respvsa_*`, `Respvsa_norm_*` |
| Alpha fits | [dofits.C](../L2Residual/dofits.C#L33-301) | `Respvsa_norm_*`, `ratio_pt*_alpha0.3` | `ratio`, `corrections_<ptrange>` |
| pT fits | [fit_pt_param.C](../L2Residual/fit_pt_param.C#L28-169) | `corrections_<ptrange>` | `loglin_etaN`, `const_etaN`, `run3_etaN` |
| Text export | [doTxt.C](../L2Residual/doTxt.C#L31-101), [L2res_param_txt.C](../L2Residual/L2res_param_txt.C#L1-97), [L2res_Run3param_txt.C](../L2Residual/L2res_Run3param_txt.C#L1-47) | `corrections_<ptrange>` or fitted functions plus k-factors | JetMET-format L2 residual text files |

## Essential histograms

### Filled in `analyse.cc`

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| [`dijetasymmetry3D*`](../fillhistograms/analyse.cc#L617) | main asymmetry profiles vs pT, eta, alpha | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L1-259) |
| [`dijetasymmetry2D_a01` ... `a05`](../fillhistograms/analyse.cc#L624) | quick QA profiles at fixed cumulative alpha cuts | plotting and validation |
| [`asymmdist3D_a10` ... `a45`](../fillhistograms/analyse.cc#L613) | asymmetry distributions for JER Gaussian-width workflow | [JERSF_fits.C](../JER/JERSF_fits.C#L1-188) |
| [`absasymmdist3D_a10` ... `a45`](../fillhistograms/analyse.cc#L614) | absolute asymmetry distributions for JER RMS workflow | [JERSF_RMS.C](../JER/JERSF_RMS.C#L1-236) |

Selection block references:

| Item | Code |
| --- | --- |
| Jet ID logic | [analyse.cc](../fillhistograms/analyse.cc#L399-423) |
| Back-to-back requirement | [analyse.cc](../fillhistograms/analyse.cc#L566-585) |
| Cumulative asymmetry fills | [analyse.cc](../fillhistograms/analyse.cc#L612-748) |
| Binning arrays | [histograms.h](../fillhistograms/histograms.h#L298-319) |

## Minimal chain detail

### 1. Histogram filling

`analyse.cc` builds the asymmetry payload

$$
A = \frac{p_{T,\mathrm{probe}} - p_{T,\mathrm{tag}}}{2 p_{T,\mathrm{avg}}}
$$

and stores it cumulatively in alpha. The L2 derivation only needs the `dijetasymmetry3D*` family.

### 2. Derivation in `deriveL2_from3D.C`

For each `(pT, eta, alpha)` bin the macro converts asymmetry to response and forms the MC/data ratio.

$$
R = \frac{1 + A}{1 - A}, \qquad \mathrm{L2Res} = \frac{R_{MC}}{R_{Data}}
$$

Outputs to keep in mind:

| Output | Meaning |
| --- | --- |
| `ratio_pt<lo>to<hi>_alpha<cut>` | eta-dependent MC/data ratio at one alpha cut |
| `Respvsa_<ptbin>_<etabin>` | raw response ratio vs alpha |
| `Respvsa_norm_<ptbin>_<etabin>` | response ratio normalized to the reference alpha bin |

References: [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L183), [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L259)

### 3. Fit order and the actual objects passed downstream

The implemented order is:

| Order | Object | Produced in | Consumed in |
| --- | --- | --- | --- |
| 1 | `ratio_pt<lo>to<hi>_alpha0.3` | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L183) | normalization reference inside [dofits.C](../L2Residual/dofits.C#L49) |
| 2 | `Respvsa_norm_<ptbin>_<etabin>` | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L259) | alpha fit inside [dofits.C](../L2Residual/dofits.C#L107) |
| 3 | `ratio` histogram of alpha-fit intercepts vs eta | [dofits.C](../L2Residual/dofits.C#L196) | multiplied into the reference ratio histograms in [dofits.C](../L2Residual/dofits.C#L280) |
| 4 | `corrections_<ptrange>` | [dofits.C](../L2Residual/dofits.C#L293) | pT fit input in [fit_pt_param.C](../L2Residual/fit_pt_param.C#L131) |

The important clarification is that the pT fits do not use the alpha-fit functions directly. They use the already alpha-corrected eta-bin values stored in `corrections_<ptrange>`.

### 4. Alpha fits in `dofits.C`

This stage uses the normalized response-vs-alpha histograms and extracts the eta-dependent correction factor.

| Fit | Expression | Purpose |
| --- | --- | --- |
| alpha fit | `pol1` | fit `Respvsa_norm_*` vs alpha and store the intercept in `ratio` |
| eta helper | `[0]+[1]*cosh(x)/(1+[2]*cosh(x))` | defined in the macro but not currently applied to `ratio` |

Outputs to keep:

| Output | Meaning |
| --- | --- |
| `ratio` | final eta-dependent k-factor histogram |
| `corrections_<ptrange>` | `ratio_pt*_alpha0.3` multiplied by `ratio`, one eta histogram per pT slice |

Reference: [dofits.C](../L2Residual/dofits.C#L196), [dofits.C](../L2Residual/dofits.C#L280), [dofits.C](../L2Residual/dofits.C#L293)

### 5. pT fits in `fit_pt_param.C`

`fit_pt_param.C` fits `corrections_<ptrange>` as a function of pT separately in each eta bin.

The pT-fit input histogram is assembled by taking the same eta bin across all `corrections_<ptrange>` histograms and filling a new histogram vs pT: [fit_pt_param.C](../L2Residual/fit_pt_param.C#L131)

| Function name in output ROOT | Expression |
| --- | --- |
| `loglin_etaN` | `[0] + [1]*log(x)` |
| `const_etaN` | `pol0` |
| `run3_etaN` | `1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))` |

Reference: [fit_pt_param.C](../L2Residual/fit_pt_param.C#L118-169)

## Text outputs

### 1. Binned export from `doTxt.C`

Header:

```text
{2 JetEta JetPt 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}
```

Each row written by [doTxt.C](../L2Residual/doTxt.C#L36-101) is:

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3-4 | validity pT min, validity pT max |
| 5 | number of function parameters |
| 6-7 | formula pT min, formula pT max |
| 8-10 | `p0 p1 p2` |

In the current implementation the last two parameters are written as `0 1`, so the binned export is effectively a piecewise constant eta correction over each pT block.

### 2. Parametrized export from `L2res_param_txt.C`

Log-linear header:

```text
{1 JetEta 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}
```

Run-3-style header:

```text
{1 JetEta 1 JetPt 1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))*[3] Correction L2Relative}
```

Row layout:

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3 | number of parameters |
| 4-5 | pT min, pT max |
| 6+ | function parameters in header order |

The last parameter is the eta-dependent k-factor from `ratio` when `kfsr` is enabled.

### 3. Which text output should I use?

| Need | Output |
| --- | --- |
| one value per eta and pT block | `doTxt.C` |
| compact log-linear parametrization | `L2res_param_txt.C` |
| Run-3-style parametrization matching current 2024 files | `L2res_Run3param_txt.C` |

## Final files worth tracking

| Product | Produced by |
| --- | --- |
| derived ROOT file | [deriveL2_from3D.C](../L2Residual/deriveL2_from3D.C#L1-259) |
| alpha-fit and eta-correction ROOT file | [dofits.C](../L2Residual/dofits.C#L33-301) |
| pT-fit ROOT file | [fit_pt_param.C](../L2Residual/fit_pt_param.C#L28-169) |
| L2 residual text file | one of the three text exporters above |