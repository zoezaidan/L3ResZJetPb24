# Systematics

The centralized systematic machinery currently exists only for the L2 chain, under `L2Residual/SYSTEMATICS/`.

## At a glance

| Source | Function | Inputs | Output |
| --- | --- | --- | --- |
| statistical | `stat_uncertainty` | `ratio` histogram from [dofits.C](../L2Residual/dofits.C#L33-301) | text file with eta-bin statistical uncertainty expanded over pT |
| parametrization | `param_uncertainty` | `loglin_etaN` or `run3_etaN` plus `const_etaN` from [fit_pt_param.C](../L2Residual/fit_pt_param.C#L118-169) | text file with model-choice uncertainty |
| JER variation | `symm_uncertainty` | nominal, up, and down fit files plus k-factors | text file with symmetrized JER uncertainty |

## Runnable Macro

The entry point is:

```cpp
void L2res_all_uncertainties()
```

This macro has no runtime parameters. The file paths, output names, eta binning, and pT support are configured inside [L2Residual/SYSTEMATICS/L2res_all_uncertainties.C](../L2Residual/SYSTEMATICS/L2res_all_uncertainties.C).

To run it:

```bash
root -l -b -q 'L2Residual/SYSTEMATICS/L2res_all_uncertainties.C()'
```

Before running it, update the configured input filenames inside the macro to point at:

| Needed input | Produced by |
| --- | --- |
| nominal `ratio` histogram | `L2Residual/dofits.C` |
| nominal `loglin_etaN`, `const_etaN`, or `run3_etaN` functions | `L2Residual/fit_pt_param.C` |
| JER up and down alternatives if used | corresponding variation fit files |

## Text output format

All three sources write a JetMET-style payload with header

```text
{1 JetEta 1 JetPt "" Correction JECSource}
```

Row layout:

| Column | Meaning |
| --- | --- |
| 1-2 | eta min, eta max |
| 3-4 | pT min, pT max |
| 5 | uncertainty value at that eta and pT point |

The script writes negative eta first and then positive eta.

## Source-specific notes

| Source | Minimal explanation | Code |
| --- | --- | --- |
| statistical | takes the error on `ratio` and repeats it over a fixed pT support grid | [L2res_all_uncertainties.C](../L2Residual/SYSTEMATICS/L2res_all_uncertainties.C#L51-94) |
| parametrization | compares `loglin_etaN` or `run3_etaN` against `const_etaN` and writes half the difference times the k-factor | [L2res_all_uncertainties.C](../L2Residual/SYSTEMATICS/L2res_all_uncertainties.C#L96-151) |
| JER variation | evaluates nominal, up, and down fits and symmetrizes the shift | [L2res_all_uncertainties.C](../L2Residual/SYSTEMATICS/L2res_all_uncertainties.C#L153-235) |

There is no equivalent repository-level L3 systematic producer in the current tree.