// This can be used to run the whole workflow for L2 residuals after the
// histograms are filled

#include "deriveL2_from3D.C"
#include "doTxt.C"
#include "dofits.C"
#include "plotresponses.C"

void runL2RES() {

  // void deriveL2_from3D(string inFileName =
  // "HIJEC_results/rerunall/MC_AK4_jetid.root",string inFileNameDT =
  // "HIJEC_results/rerunall/zerobiasall_jetid.root", string outfilename =
  // "RERUNALL/L2residuals_pbpbreco_rerunall_zerobias_jetid.root", bool dodt =
  // true,   int alphabin = 5, bool useabs = true, bool usewideabs = false) {

  TString dir = "/home/laura/Code/jec/HIJEC_rereco_results/";
  ///// Initial responses
  // deriveL2_from3D(dir+"RERECOMC_AK4_PFTRIG_jetid.root",dir+"RERECO_ZB_ALL_PFTRIG_jetid.root","L2residuals_pbpbreco_rereco_zb_jetid.root",
  // true, 5, true, false);
  // deriveL2_from3D(dir+"RERECOMC_AK4_PFTRIG_jetid.root",dir+"RERECOHP_AK4_PFTRIG_jetid.root","L2residuals_pbpbreco_rereco_hp_jetid.root",
  // true, 5, true, false);

  // Plot responses
  plotresponses("L2residuals_pbpbreco_rereco_zb_jetid.root", "ZB");
  plotresponses("L2residuals_pbpbreco_rereco_hp_jetid.root", "HP");

  // Fit
  //  dofits("L2residuals_pbpbreco_rereco_zb_jetid.root",
  //  "L2residuals_pbpbreco_rereco_hp_jetid.root",  0.15,  0.35,
  //  "correctionfile", true);
  // Produce txt
  // doTxt("L2fits/correctionfile.root", "L2residual_corrections.txt");

  ////// L2 closure without JER

  // deriveL2_from3D(dir+"RERECOMC_AK4_PFTRIG_jetid.root",dir+"RERECO_ZB_ALL_PFTRIG_l2closure.root","L2residuals_pbpbreco_rereco_zb_jetid_l2jecclosure.root",
  // true, 5, true, false);
  // deriveL2_from3D(dir+"RERECOMC_AK4_PFTRIG_jetid.root",dir+"RERECOHP_AK4_PFTRIG_jetid_l2jecclosure.root","L2residuals_pbpbreco_rereco_hp_jetid_l2jecclosure.root",
  // true, 5, true, false);

  // Plot responses
  plotresponses("L2residuals_pbpbreco_rereco_zb_jetid_l2jecclosure.root",
                "ZB_l2resclosure");
  plotresponses("L2residuals_pbpbreco_rereco_hp_jetid_l2jecclosure.root",
                "HP_l2resclosure");

  //////// L2 closure with JER SF
}
