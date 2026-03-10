# TODO Index

This file tracks every in-code `TODO` currently present in the repository so they can be worked through systematically.

## fillhistograms/analyse.cc

- [analyse.cc TODO at line 47](fillhistograms/analyse.cc#L47): smarter handling for the valid-jet range after applying L2.
- [analyse.cc TODO at line 455](fillhistograms/analyse.cc#L455): remove the hardcoded radius.

## fillhistograms/analyse_PhotonJet.cc

- [analyse_PhotonJet.cc TODO at line 87](fillhistograms/analyse_PhotonJet.cc#L87): smarter handling for the tightly limited valid-jet range.
- [analyse_PhotonJet.cc TODO at line 366](fillhistograms/analyse_PhotonJet.cc#L366): add the electron veto for photons.

## fillhistograms/analyse_JER.cc

- [analyse_JER.cc TODO at line 27](fillhistograms/analyse_JER.cc#L27): clarify local-environment handling.
- [analyse_JER.cc TODO at line 45](fillhistograms/analyse_JER.cc#L45): smarter valid-jet range handling after L2.
- [analyse_JER.cc TODO at line 58](fillhistograms/analyse_JER.cc#L58): add safety checks for file opening.
- [analyse_JER.cc TODO at line 59](fillhistograms/analyse_JER.cc#L59): add safety checks for file opening.
- [analyse_JER.cc TODO at line 265](fillhistograms/analyse_JER.cc#L265): improve the current placeholder logic.
- [analyse_JER.cc TODO at line 329](fillhistograms/analyse_JER.cc#L329): invert the logic to start with `true` and disable only when needed.
- [analyse_JER.cc TODO at line 338](fillhistograms/analyse_JER.cc#L338): implement constituent handling.
- [analyse_JER.cc TODO at line 400](fillhistograms/analyse_JER.cc#L400): review `passjetid` handling.
- [analyse_JER.cc TODO at line 432](fillhistograms/analyse_JER.cc#L432): use the above-defined tag and probe indices.
- [analyse_JER.cc TODO at line 661](fillhistograms/analyse_JER.cc#L661): add tag-and-probe versions for the fill without tag and probe.

## fillhistograms/eventhistograms.h

- [eventhistograms.h TODO at line 19](fillhistograms/eventhistograms.h#L19): add `rho` support.

## L2Residual/deriveL2_from3D.C

- [deriveL2_from3D.C TODO at line 23](L2Residual/deriveL2_from3D.C#L23): decide whether and how to save the full response information from different alpha values.

## L2Residual/dofits.C

- [dofits.C TODO at line 73](L2Residual/dofits.C#L73): copy histogram contents into `TGraph`s and exclude the reference alpha.
- [dofits.C TODO at line 74](L2Residual/dofits.C#L74): remove the reference alpha from the fit.
- [dofits.C TODO at line 83](L2Residual/dofits.C#L83): revisit the bin-range debug snippet.
- [dofits.C TODO at line 146](L2Residual/dofits.C#L146): avoid including the reference point in the graph fit.
- [dofits.C TODO at line 151](L2Residual/dofits.C#L151): replace this with a `TMultiGraph` fitter if needed.
- [dofits.C TODO at line 158](L2Residual/dofits.C#L158): correct the reference label in the y-axis title.
- [dofits.C TODO at line 169](L2Residual/dofits.C#L169): improve the color handling.
- [dofits.C TODO at line 199](L2Residual/dofits.C#L199): debug the current section.
- [dofits.C TODO at line 204](L2Residual/dofits.C#L204): correct the signed-eta bin label.
- [dofits.C TODO at line 205](L2Residual/dofits.C#L205): correct the absolute-eta bin label.

## JER/JERSF_printtxt.C

- [JERSF_printtxt.C TODO at line 2](JER/JERSF_printtxt.C#L2): add uncertainty handling to the output.

## JER/JERSF_polfits.C

- [JERSF_polfits.C TODO at line 4](JER/JERSF_polfits.C#L4): fix the input naming.
- [JERSF_polfits.C TODO at line 5](JER/JERSF_polfits.C#L5): fit using `TGraph`s.
- [JERSF_polfits.C TODO at line 10](JER/JERSF_polfits.C#L10): decide whether to fit data and MC separately or fit the ratio directly.
- [JERSF_polfits.C TODO at line 99](JER/JERSF_polfits.C#L99): fix the latex text drawing on the top pad.
- [JERSF_polfits.C TODO at line 110](JER/JERSF_polfits.C#L110): take ratios and extract SF per eta with possible pT dependence.

## JER/JERSF_fits_vsalpha.C

- [JERSF_fits_vsalpha.C TODO at line 125](JER/JERSF_fits_vsalpha.C#L125): decide whether the pT dependence should be fit with `pol1` instead of `pol0`.

## JER/JERSF_fits.C

- [JERSF_fits.C TODO at line 50](JER/JERSF_fits.C#L50): make the HardProbes versus ZeroBias choice more flexible for 2024.
- [JERSF_fits.C TODO at line 63](JER/JERSF_fits.C#L63): decide whether the Gaussian should always be fit twice.
- [JERSF_fits.C TODO at line 176](JER/JERSF_fits.C#L176): avoid including the reference alpha point in the MC graph.
- [JERSF_fits.C TODO at line 179](JER/JERSF_fits.C#L179): avoid including the reference alpha point in the data graph.

## JER/JERSF.C

- [JERSF.C TODO at line 22](JER/JERSF.C#L22): decide whether to store sigmas in histograms or vectors.

## JER/JERSF_RMS.C

- [JERSF_RMS.C TODO at line 80](JER/JERSF_RMS.C#L80): make the HardProbes versus ZeroBias choice more flexible for 2024.
- [JERSF_RMS.C TODO at line 134](JER/JERSF_RMS.C#L134): replace the current loop logic with a clearer `while` loop.
- [JERSF_RMS.C TODO at line 155](JER/JERSF_RMS.C#L155): improve the line maximum used for the truncation marker.
- [JERSF_RMS.C TODO at line 225](JER/JERSF_RMS.C#L225): avoid including the reference alpha point in the MC graph.
- [JERSF_RMS.C TODO at line 228](JER/JERSF_RMS.C#L228): avoid including the reference alpha point in the data graph.