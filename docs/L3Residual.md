# L3Residual

This page keeps the essential L3 chain: which photon+jet histograms matter, which objects are passed through the derivation and fit steps, and how the final text files are written.

## At a glance

| Stage | Macro | Essential input histograms or objects | Essential output passed downstream |
| --- | --- | --- | --- |
| Histogram filling | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L70-930) | photons, away-side jets, alpha thresholds | `photonjet_balance3D*`, `photonjet_balance*_a*`, `photonjet_balance_dist` |
| Derivation | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L250-470) | `photonjet_balance3D*` from MC and data | `L3Res_vsa*`, `L3Res_vsa_norm_*`, `ratio_vsphotonpt_alphaN`, `ratio_norm_vsphotonpt_alphaN` |
| Final fits | [dofits_L3.C](../L3Residual/dofits_L3.C#L228-1368) | derived ratio histograms, optional external L2 file | `kFSR_vsPt`, final L3 fit, optional L2L3 fit |
| Text export | [dofits_L3.C](../L3Residual/dofits_L3.C#L1256-1364) | final fitted function, optional L2 records | L3Residual and optional L2L3Residual text files |

## Essential histograms

### Filled in `analyse_PhotonJet.cc`

| Histogram family | Meaning | Used by |
| --- | --- | --- |
| [`photonjet_balance3D*`](../fillhistograms/analyse_PhotonJet.cc#L863) | balance profile vs photon pT, eta, alpha | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L250-470) |
| [`photonjet_balance3Dabseta*`](../fillhistograms/analyse_PhotonJet.cc#L866) | same payload in absolute eta | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L250-470) |
| [`photonjet_balance3D*_counts`](../fillhistograms/analyse_PhotonJet.cc#L874) | event counts for QA and uncertainty handling | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L626-633) |
| [`photonjet_balance_dist`](../fillhistograms/analyse_PhotonJet.cc#L880) | balance distribution vs photon pT and alpha | distribution overlays and checks |

Selection references:

| Item | Code |
| --- | --- |
| photon preselection and ID | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L600-652) |
| MC gen matching | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L613-628) |
| away-side jet selection | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L656-704) |
| veto-map application | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L748-769) |
| cumulative alpha fills | [analyse_PhotonJet.cc](../fillhistograms/analyse_PhotonJet.cc#L817-880) |
| binning arrays | [histograms.h](../fillhistograms/histograms.h#L298-319) |

## Minimal chain detail

### 1. Histogram filling

The propagated response variable is

$$
\mathrm{balance} = \frac{p_{T,\mathrm{jet}}}{p_{T,\gamma}}.
$$

The second away-side jet defines

$$
\alpha = \frac{p_{T,\mathrm{2nd\ away\ jet}}}{p_{T,\gamma}}.
$$

Only the `photonjet_balance3D*` family is required for the main L3 derivation.

### 2. Derivation in `deriveL3_from_photonjet.C`

This stage forms MC/data balance ratios by eta and by photon pT.

| Output | Meaning |
| --- | --- |
| `L3Res_vsa_<ptbin>_<etabin>` | MC/data ratio vs alpha in one pT and eta bin |
| `L3Res_vsa_norm_<ptbin>_<etabin>` | same ratio normalized to a reference alpha bin |
| `ratio_vsphotonpt_alphaN` | MC/data ratio vs photon pT for one cumulative alpha cut |
| `ratio_norm_vsphotonpt_alphaN` | same ratio divided by the reference-alpha ratio |

Reference: [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L250-470)

### 3. Final fits in `dofits_L3.C`

The final L3 fit function is

$$
[0] + [1] \log_{10}(0.01 x)
$$

defined in [dofits_L3.C](../L3Residual/dofits_L3.C#L228-240).

Optional `kFSR` path:

| Output | Meaning | Code |
| --- | --- | --- |
| `kFSR_vsPt` | alpha-intercept of normalized ratio vs alpha | [dofits_L3.C](../L3Residual/dofits_L3.C#L380-382), [dofits_L3.C](../L3Residual/dofits_L3.C#L594-752) |
| `ratio_vspT_alpha0_fromkFSR` | reference-alpha ratio promoted to `alpha -> 0` | [dofits_L3.C](../L3Residual/dofits_L3.C#L860-923) |

Combined-input mode:

- if several derived files are provided, each input should already contain `corr_vspT` or `ratio_vspT_alpha0`
- the macro does not recompute `kFSR` independently per input in combined mode

Reference: [dofits_L3.C](../L3Residual/dofits_L3.C#L303-318)

## Text outputs

### 1. L3 residual text file

Header written by [dofits_L3.C](../L3Residual/dofits_L3.C#L1256-1279):

```text
{ 1 JetEta 1 JetPt [0]+[1]*log10(0.01*x) Correction L3Residual}
```

Row layout:

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3 | payload size, fixed to `4` in this writer |
| 4-5 | pT min, pT max |
| 6-7 | fit parameters `p0 p1` |

If an L2 file is provided, its eta binning is reused for the L3 text output. Otherwise one wide eta bin is written.

Two copies are written:

| Path style | Purpose |
| --- | --- |
| `<tag>/textfiles/<tag>.txt` | local fit output inside the result directory |
| `jecfiles/L3Residuals_<run>_<tag>_AK4PF.txt` | JEC-style export |

### 2. Combined L2L3 residual text file

If `writeL2L3=true`, [dofits_L3.C](../L3Residual/dofits_L3.C#L1307-1364) reads the L2 file, samples the product

$$
\mathrm{L2}(p_T) \times \mathrm{L3}(p_T)
$$

and refits it to the same functional form already used by the L2 file.

Important behavior:

| Item | Meaning |
| --- | --- |
| header | copied from the L2 file, with `L2Relative` renamed to `L2L3Residual` |
| eta bins | copied from the L2 file |
| pT range | copied from each L2 record |
| parameters | refitted combined parameters per eta bin |

## Final files worth tracking

| Product | Produced by |
| --- | --- |
| derived ROOT file | [deriveL3_from_photonjet.C](../L3Residual/deriveL3_from_photonjet.C#L250-470) |
| fit ROOT file and plots | [dofits_L3.C](../L3Residual/dofits_L3.C#L228-1368) |
| L3 residual text file | [dofits_L3.C](../L3Residual/dofits_L3.C#L1256-1279) |
| optional L2L3 residual text file | [dofits_L3.C](../L3Residual/dofits_L3.C#L1307-1364) |