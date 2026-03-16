#include "histograms.h"
histograms::histograms(TDirectory *dir, float etamin, float etamax, float hibinmin, float hibinmax, bool ismc) {
 //histograms::histograms(float ptmin, float ptmax, bool ismc) {
  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;
  
  this->etamin = etamin;
  this->etamax = etamax;
  this->isMC = ismc;
  this->hibinmin = hibinmin;
  this->hibinmax = hibinmax;
  
// Jet histograms
  jet_pt = new TH1D("reco jet pT", "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_pt_now = new TH1D("reco jet pT, no evt w", "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_uncorr_pt = new TH1D("reco jet pT uncorr", "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_pt_genweight = new TH1D("reco jet pT, gen w", "reco jet p_{T}, gen weight only; reco jet p_{T};", 100, 100, 1000);
  jet_eta = new TH1D("reco jet eta"," reco jet #eta; reco jet #eta;", 40, -5.2, 5.2);
  jet_phi = new TH1D("reco jet phi"," reco jet #phi; reco jet #phi;", 25, -3.1415926535, 3.1415926535);

  tag_pt = new TH1D("tag jet pT", "tag jet p_{T}; tag jet p_{T};", 100, 15, 1000);
  tag_eta = new TH1D("tag jet eta","tag jet #eta; tag jet #eta;", 40, -5.2, 5.2);
  tag_phi = new TH1D("tag jet phi","tag jet #phi; tag jet #phi;", 25, -3.1415926535, 3.1415926535);

  probe_pt = new TH1D("probe jet pT", "probe jet p_{T}; probe jet p_{T};", 100, 15, 1000);
  probe_eta = new TH1D("probe jet eta"," probe jet #eta; probe jet #eta;", 40, -5.2, 5.2);
  probe_phi = new TH1D("probe jet phi"," probe jet #phi; probe jet #phi;", 25, -3.1415926535, 3.1415926535);

  alphas = new TH1D("alpha","alpha;alpha", 10, 0.0, 0.5);

  // Trigger checks
  HLTZB = new TH1D("HLTZB", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);
  HLT40 = new TH1D("HLT40", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);
  HLT60 = new TH1D("HLT60", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);
  HLT80 = new TH1D("HLT80", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);
  HLT100 = new TH1D("HLT100", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);
  HLT120 = new TH1D("HLT120", "leading jet p_{T}; leading jet p_{T};", 100, 0, 1000);

  HLTZB_ptav = new TH1D("HLTZB_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
  HLT40_ptav = new TH1D("HLT40_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
  HLT60_ptav = new TH1D("HLT60_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
  HLT80_ptav = new TH1D("HLT80_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
  HLT100_ptav = new TH1D("HLT100_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
  HLT120_ptav = new TH1D("HLT120_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);

  // PF energy fraction composition
  jet_nhf = new TProfile("reco jet nhf", "reco jet nhf; reco jet p_{T};", 100, 15, 1000);
  jet_chf = new TProfile("reco jet chf", "reco jet chf; reco jet p_{T};", 100, 15, 1000);
  jet_nef = new TProfile("reco jet nef", "reco jet nef; reco jet p_{T};", 100, 15, 1000); 
  jet_cef = new TProfile("reco jet cef", "reco jet cef; reco jet p_{T};", 100, 15, 1000);
  jet_muf = new TProfile("reco jet muf", "reco jet muf; reco jet p_{T};", 100, 15, 1000);
  jetetaphi = new TH2D("eta-phi distribution",";#eta; #phi",netas,etarange,nphis,phirange);
  
  if (ismc) {
    genjet_pt = new TH1D("gen jet pT", "gen jet p_{T}; gen jet p_{T};", 100, 100, 1000);
    genjet_eta = new TH1D("gen jet eta"," gen jet #eta; gen jet #eta;", 20, -2.5, 2.5);
    genjet_phi = new TH1D("gen jet phi"," gen jet #phi; gen jet #phi;", 20, -2.5, 2.5);
    jetresponse = new TProfile("response","",100,100,1000);

// "Residuals"
// Definition: (reco-true)/true
    ptres = new TH1D("pT res"," pT ; (pT_reco-pT_gen)/pT_gen;", 40, -2, 2);
 
  }

  // Dijets
  dijetasymmetry = new TH1D("dijetasymmetry","  ; asymmetry;", 40, 0, 1);
  dijetasymmetry_now = new TH1D("dijetasymmetry_now","  ; asymmetry;", 40, 0, 1);
  dijetdeltaphi  = new TH1D("dijetdeltaphi"," ; delta phi;", 40, 0, 3.1415926535);
  dijetdeltaeta  = new TH1D("dijetdeltaeta"," ; delta eta;", 40, 0, 5.2);
  //  dijetbalance_a01 = new TH1D("dijetbalance_a01","  ; ;", 20, -2, 2);
  dijetasymmetry_a01 = new TProfile("dijetasymmetry_a01","  ; ;",  nptforjec, &ptforjec[0]);
  // dijetbalance_a02 = new TH1D("dijetbalance_a02","  ; ;", 20, -2, 2);
  dijetasymmetry_a02 = new TProfile("dijetasymmetry_a02","  ; ;",  nptforjec, &ptforjec[0]);
  dijetbalance_a03 = new TH1D("dijetbalance_a03","  ; ;", 20, -2, 2);


  dijetbalance_a03_pt30to40 = new TH2D("dijetbal_pt30to40", "", nwabsetas, wabsetarange, nasym, asymrange);
  dijetbalance_a03_pt40to80 = new TH2D("dijetbal_pt40to80", "", nwabsetas, wabsetarange, nasym, asymrange);

  dijetasymmetry_a03 = new TProfile("dijetasymmetry_a03","  ; ;",  nptforjec, &ptforjec[0]);
  //  dijetbalance_a035 = new TH1D("dijetbalance_a035","  ; ;", 20, -2, 2);
  dijetasymmetry_a035 = new TProfile("dijetasymmetry_a035","  ; ;",  nptforjec, &ptforjec[0]);
  
  //  dijetbalance_a04 = new TH1D("dijetbalance_a04","  ; ;", 20, -2, 2);
  dijetasymmetry_a04 = new TProfile("dijetasymmetry_a04","  ; ;",  nptforjec, &ptforjec[0]);
  //  dijetbalance_a05 = new TH1D("dijetbalance_a05","  ; ;", 20, -2, 2);
  dijetasymmetry_a05 = new TProfile("dijetasymmetry_a05","  ; ;",  nptforjec, &ptforjec[0]);
  //  dijetbalance_a06 = new TH1D("dijetbalance_a06","  ; ;", 20, -2, 2);
  dijetasymmetry_a06 = new TProfile("dijetasymmetry_a06","  ; ;",  nptforjec, &ptforjec[0]);
  
  dijetbalance_a1= new TH1D("dijetbalance_a1","  ; ;", 20, -2, 2);
  dijetasymmetry_a1 = new TProfile("dijetasymmetry_a1","  ; ;",  nptforjec, &ptforjec[0]);

  if ((this->etamin - this->etamax) < -10) {

    dijetasymmetry3D = new TProfile3D("dijetasymmetry3D", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0], nalphavalues, &alphavalues[0]); 
    dijetasymmetry3Dabseta = new TProfile3D("dijetasymmetry3Dabseta", ";;", nptforjec, &ptforjec[0], nwabsetas, &wabsetarange[0], nalphavalues, &alphavalues[0]);
    dijetasymmetry3Dabsetawide = new TProfile3D("dijetasymmetry3Dabsetawide", ";;", nptforjec, &ptforjec[0], ndwabsetas, &dwabsetarange[0], nalphavalues, &alphavalues[0]); 
    dijetasymmetry3Dnarrow = new TProfile3D("dijetasymmetry3Dnarrow", ";;", nptforjec, &ptforjec[0], netas, &etarange[0], nalphavalues, &alphavalues[0]); 
    dijetasymmetry3Dabsetanarrow = new TProfile3D("dijetasymmetry3Dabsetanarrow", ";;", nptforjec, &ptforjec[0], nabsetas, &absetarange[0], nalphavalues, &alphavalues[0]); 
    dijetasymmetry2D_a01 = new TProfile2D("dijetasymmetry2D_a01", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    dijetasymmetry2D_a02 = new TProfile2D("dijetasymmetry2D_a02", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    dijetasymmetry2D_a03 = new TProfile2D("dijetasymmetry2D_a03", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    dijetasymmetry2D_a04 = new TProfile2D("dijetasymmetry2D_a04", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    dijetasymmetry2D_a05 = new TProfile2D("dijetasymmetry2D_a05", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    dijetasymmetry2D_a06 = new TProfile2D("dijetasymmetry2D_a06", ";;", nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
  }
  ptgenvsptreco = new TH2D("ptgenvsptreco","",200,0,1500,200,0,1500);
  ptrecovsweight = new TH2D("ptrecovsweight","",200,0,1500,200,0,0.1);
  ptgenvsweight = new TH2D("ptgenvsweight","",200,0,1500,200,0,0.1);

  // JER STUFF:
  vector<double> x(61);
  for (unsigned int i = 0; i != x.size(); ++i) x[i] = 0.05*i;
  const int nx = x.size()-1;

    vector<double> y(161);
  for (unsigned int i = 0; i != y.size(); ++i) y[i] = -1 + 0.0125*i;
  const int ny = y.size()-1;
  if (ismc) {
    responses3D = new TH3D("responses3D",";;", nptforJER, &ptforJER[0], netaforjer, &etaforjer[0], nx, &x[0]);
    phiresponse = new TH3D("phiresponses3D",";;", nptforJER, &ptforJER[0], netaforjer, &etaforjer[0], ny, &y[0]);
    etaresponse = new TH3D("etaresponses3D",";;", nptforJER, &ptforJER[0], netaforjer, &etaforjer[0], ny, &y[0]);
  }
  absasymmdist3D = new TH3D("absasymmdist3D",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D = new TH3D("asymmdist3D",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);

  // JER needs asymmetries as function of alpha
  asymmdist3D_a10 = new TH3D("asymmdist3D_a10",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a10  = new TH3D("absasymmdist3D_a10",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a15 = new TH3D("asymmdist3D_a15",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a15 = new TH3D("absasymmdist3D_a15",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a20 = new TH3D("asymmdist3D_a20",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a20 = new TH3D("absasymmdist3D_a20",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a25 = new TH3D("asymmdist3D_a25",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a25 = new TH3D("absasymmdist3D_a25",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a30 = new TH3D("asymmdist3D_a30",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a30 = new TH3D("absasymmdist3D_a30",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a35 = new TH3D("asymmdist3D_a35",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a35 = new TH3D("absasymmdist3D_a35",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a40 = new TH3D("asymmdist3D_a40",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a40 = new TH3D("absasymmdist3D_a40",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  asymmdist3D_a45 = new TH3D("asymmdist3D_a45",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], ny, &y[0]);
  absasymmdist3D_a45 = new TH3D("absasymmdist3D_a45",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);

  if (ismc) {
    absasymmdist3D_gen_a10 = new TH3D("absasymmdist3D_gen_a10",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a15 = new TH3D("absasymmdist3D_gen_a15",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a20 = new TH3D("absasymmdist3D_gen_a20",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a25 = new TH3D("absasymmdist3D_gen_a25",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a30 = new TH3D("absasymmdist3D_gen_a30",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a35 = new TH3D("absasymmdist3D_gen_a35",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a40 = new TH3D("absasymmdist3D_gen_a40",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
    absasymmdist3D_gen_a45 = new TH3D("absasymmdist3D_gen_a45",";;", nptforjec, &ptforjec[0], njeretas, &jeretarange[0], nx, &x[0]);
  }
  
}
void histograms::Write() {
  dir->cd();
  dir->Write();
  }
// HistosBasic::~HistosBasic() {
//  Write();/
// };
