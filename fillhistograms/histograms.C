#include "histograms.h"

// Helper function to initialize all pointers to nullptr
void histograms::initializePointers() {
  // Common jet histograms
  jet_pt = nullptr;
  jet_pt_now = nullptr;
  jet_uncorr_pt = nullptr;
  jet_pt_genweight = nullptr;
  jet_eta = nullptr;
  jet_phi = nullptr;
  tag_pt = nullptr;
  tag_eta = nullptr;
  tag_phi = nullptr;
  probe_pt = nullptr;
  probe_eta = nullptr;
  probe_phi = nullptr;

  // Triggers
  HLTZB = nullptr;
  HLT40 = nullptr;
  HLT60 = nullptr;
  HLT80 = nullptr;
  HLT100 = nullptr;
  HLT120 = nullptr;
  HLTZB_ptav = nullptr;
  HLT40_ptav = nullptr;
  HLT60_ptav = nullptr;
  HLT80_ptav = nullptr;
  HLT100_ptav = nullptr;
  HLT120_ptav = nullptr;

  // PF composition
  jet_nhf = nullptr;
  jet_chf = nullptr;
  jet_nef = nullptr;
  jet_cef = nullptr;
  jet_muf = nullptr;
  jetetaphi = nullptr;

  // Gen jet histograms
  genjet_pt = nullptr;
  genjet_eta = nullptr;
  genjet_phi = nullptr;
  genjetdyn_kt = nullptr;
  genjetdyn_deltaR = nullptr;
  genjetdyn_z = nullptr;
  genjetdyn_z_cutdeltaR = nullptr;
  genphoton_pt = nullptr;
  genphoton_eta = nullptr;
  genphoton_phi = nullptr;
  plane_inclusive = nullptr;
  jetresponse = nullptr;
  photonresponse = nullptr;
  ptres = nullptr;
  photon_ptres = nullptr;
  ktres = nullptr;
  deltaRres = nullptr;
  zres = nullptr;

  // Dijet histograms
  dijetasymmetry = nullptr;
  dijetasymmetry_now = nullptr;
  dijetdeltaphi = nullptr;
  dijetdeltaeta = nullptr;
  dijetasymmetry_a01 = nullptr;
  dijetasymmetry_a02 = nullptr;
  dijetasymmetry_a03 = nullptr;
  dijetasymmetry_a035 = nullptr;
  dijetasymmetry_a04 = nullptr;
  dijetasymmetry_a05 = nullptr;
  dijetasymmetry_a06 = nullptr;
  dijetasymmetry_a1 = nullptr;
  dijetbalance_a03 = nullptr;
  dijetbalance_a1 = nullptr;
  dijetasymmetry2D_a01 = nullptr;
  dijetasymmetry2D_a02 = nullptr;
  dijetasymmetry2D_a03 = nullptr;
  dijetasymmetry2D_a04 = nullptr;
  dijetasymmetry2D_a05 = nullptr;
  dijetasymmetry2D_a06 = nullptr;
  dijetasymmetry3D = nullptr;
  dijetasymmetry3Dwide = nullptr;
  dijetasymmetry3Dnarrow = nullptr;
  dijetasymmetry3Dabseta = nullptr;
  dijetasymmetry3Dabsetawide = nullptr;
  dijetasymmetry3Dabsetanarrow = nullptr;

  // Photon+Jet histograms
  photon_pt = nullptr;
  photon_eta = nullptr;
  photon_phi = nullptr;
  photon_HoverE = nullptr;
  photon_sigmaIetaIeta = nullptr;
  photon_SwissCrx = nullptr;
  photon_SeedTime = nullptr;
  awayside_jet_pt = nullptr;
  awayside_jet_eta = nullptr;
  awayside_jet_phi = nullptr;
  awayside_jet_uncorr_pt = nullptr;
  photonjet_dphi = nullptr;
  photonjet_balance = nullptr;
  photonjet_ptavg = nullptr;
  photonjet_alpha = nullptr;
  photonjet_balance_a01 = nullptr;
  photonjet_balance_a02 = nullptr;
  photonjet_balance_a03 = nullptr;
  photonjet_balance_a04 = nullptr;
  photonjet_balance_a05 = nullptr;
  photonjet_balance_a06 = nullptr;
  photonjet_balance2D_a01 = nullptr;
  photonjet_balance2D_a02 = nullptr;
  photonjet_balance2D_a03 = nullptr;
  photonjet_balance2D_a04 = nullptr;
  photonjet_balance2D_a05 = nullptr;
  photonjet_balance2D_a06 = nullptr;
  photonjet_balance3D = nullptr;
  photonjet_balance3Dwide = nullptr;
  photonjet_balance3Dnarrow = nullptr;
  photonjet_balance3Dabseta = nullptr;
  photonjet_balance3Dabsetawide = nullptr;
  photonjet_balance3Dabsetanarrow = nullptr;
  photonjet_balance3D_jetpt = nullptr;
  photonjet_balance3Dwide_jetpt = nullptr;
  photonjet_balance3Dnarrow_jetpt = nullptr;
  photonjet_balance3Dabseta_jetpt = nullptr;
  photonjet_balance3Dabsetawide_jetpt = nullptr;
  photonjet_balance3Dabsetanarrow_jetpt = nullptr;
  photonjet_balance3D_counts = nullptr;
  photonjet_balance3Dwide_counts = nullptr;
  photonjet_balance3Dnarrow_counts = nullptr;
  photonjet_balance3Dabseta_counts = nullptr;
  photonjet_balance3Dabsetawide_counts = nullptr;
  photonjet_balance3Dabsetanarrow_counts = nullptr;
  HLTPhoton30 = nullptr;
  HLTPhoton30_ptav = nullptr;

  // JER histograms
  responses3D = nullptr;
  phiresponse = nullptr;
  etaresponse = nullptr;
  asymmdist3D = nullptr;
  absasymmdist3D = nullptr;
  asymmdist3D_a10 = nullptr;
  absasymmdist3D_a10 = nullptr;
  asymmdist3D_a15 = nullptr;
  absasymmdist3D_a15 = nullptr;
  asymmdist3D_a20 = nullptr;
  absasymmdist3D_a20 = nullptr;
  asymmdist3D_a25 = nullptr;
  absasymmdist3D_a25 = nullptr;
  asymmdist3D_a30 = nullptr;
  absasymmdist3D_a30 = nullptr;
  asymmdist3D_a35 = nullptr;
  absasymmdist3D_a35 = nullptr;
  asymmdist3D_a40 = nullptr;
  absasymmdist3D_a40 = nullptr;
  asymmdist3D_a45 = nullptr;
  absasymmdist3D_a45 = nullptr;

  // Weight histograms
  ptgenvsptreco = nullptr;
  ptrecovsweight = nullptr;
  ptgenvsweight = nullptr;
}

// Original constructor - backward compatible, creates all histograms
histograms::histograms(TDirectory *dir, float etamin, float etamax,
                       float hibinmin, float hibinmax, bool ismc)
    : histograms(dir, etamin, etamax, hibinmin, hibinmax, ismc,
                 AnalysisType::ALL) {}

// Overloaded constructor with analysis type selection
histograms::histograms(TDirectory *dir, float etamin, float etamax,
                       float hibinmin, float hibinmax, bool ismc,
                       AnalysisType type) {
  TDirectory *curdir = gDirectory;
  bool enter = dir->cd();
  assert(enter);
  this->dir = dir;

  this->etamin = etamin;
  this->etamax = etamax;
  this->isMC = ismc;
  this->hibinmin = hibinmin;
  this->hibinmax = hibinmax;
  this->analysisType = type;

  // Initialize all pointers to nullptr first
  initializePointers();

  bool createDijet = (type == AnalysisType::DIJET || type == AnalysisType::ALL);
  bool createPhotonJet =
      (type == AnalysisType::PHOTONJET || type == AnalysisType::ALL);

  // ============================================
  // COMMON HISTOGRAMS (created for all analyses)
  // ============================================

  // Basic jet histograms
  jet_pt =
      new TH1D("reco jet pT", "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_pt_now = new TH1D("reco jet pT, no evt w",
                        "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_uncorr_pt = new TH1D("reco jet pT uncorr",
                           "reco jet p_{T}; reco jet p_{T};", 100, 15, 1000);
  jet_pt_genweight = new TH1D(
      "reco jet pT, gen w", "reco jet p_{T}, gen weight only; reco jet p_{T};",
      100, 100, 1000);
  jet_eta =
      new TH1D("reco jet eta", " reco jet #eta; reco jet #eta;", 40, -5.2, 5.2);
  jet_phi = new TH1D("reco jet phi", " reco jet #phi; reco jet #phi;", 25,
                     -3.1415926535, 3.1415926535);

  // PF energy fraction composition
  jet_nhf = new TProfile("reco jet nhf", "reco jet nhf; reco jet p_{T};", 100,
                         15, 1000);
  jet_chf = new TProfile("reco jet chf", "reco jet chf; reco jet p_{T};", 100,
                         15, 1000);
  jet_nef = new TProfile("reco jet nef", "reco jet nef; reco jet p_{T};", 100,
                         15, 1000);
  jet_cef = new TProfile("reco jet cef", "reco jet cef; reco jet p_{T};", 100,
                         15, 1000);
  jet_muf = new TProfile("reco jet muf", "reco jet muf; reco jet p_{T};", 100,
                         15, 1000);
  jetetaphi = new TH2D("eta-phi distribution", ";#eta; #phi", netas, etarange,
                       nphis, phirange);

  // MC-only common histograms
  if (ismc) {
    genjet_pt =
        new TH1D("gen jet pT", "gen jet p_{T}; gen jet p_{T};", 100, 100, 1000);
    genjet_eta =
        new TH1D("gen jet eta", " gen jet #eta; gen jet #eta;", 20, -2.5, 2.5);
    genjet_phi =
        new TH1D("gen jet phi", " gen jet #phi; gen jet #phi;", 20, -2.5, 2.5);
    jetresponse = new TProfile("response", "", 100, 100, 1000);
    ptres = new TH1D("pT res", " pT ; (pT_reco-pT_gen)/pT_gen;", 40, -2, 2);
    genphoton_pt = new TH1D("gen photon pT",
                            "gen photon p_{T}; gen photon p_{T};", 100, 0, 500);
    genphoton_eta = new TH1D(
        "gen photon eta", "gen photon #eta; gen photon #eta;", 50, -2.5, 2.5);
    genphoton_phi = new TH1D(
        "gen photon phi", "gen photon #phi; gen photon #phi;", nphis, phirange);
    photonresponse = new TProfile("photon response", "", 100, 0, 500);
    photon_ptres = new TH1D("photon pT res",
                            "photon p_{T}; (p_{T,reco}-p_{T,gen})/p_{T,gen};",
                            40, -0.5, 0.5);
  }

  // Weight histograms (common)
  ptgenvsptreco = new TH2D("ptgenvsptreco", "", 200, 0, 1500, 200, 0, 1500);
  ptrecovsweight = new TH2D("ptrecovsweight", "", 200, 0, 1500, 200, 0, 0.1);
  ptgenvsweight = new TH2D("ptgenvsweight", "", 200, 0, 1500, 200, 0, 0.1);

  // ============================================
  // DIJET-SPECIFIC HISTOGRAMS
  // ============================================
  if (createDijet) {
    tag_pt =
        new TH1D("tag jet pT", "tag jet p_{T}; tag jet p_{T};", 100, 15, 1000);
    tag_eta =
        new TH1D("tag jet eta", "tag jet #eta; tag jet #eta;", 40, -5.2, 5.2);
    tag_phi = new TH1D("tag jet phi", "tag jet #phi; tag jet #phi;", 25,
                       -3.1415926535, 3.1415926535);

    probe_pt = new TH1D("probe jet pT", "probe jet p_{T}; probe jet p_{T};",
                        100, 15, 1000);
    probe_eta = new TH1D("probe jet eta", " probe jet #eta; probe jet #eta;",
                         40, -5.2, 5.2);
    probe_phi = new TH1D("probe jet phi", " probe jet #phi; probe jet #phi;",
                         25, -3.1415926535, 3.1415926535);

    alphas = new TH1D("alpha", "alpha;alpha", 10, 0.0, 0.5);

    // Dijet triggers
    HLTZB = new TH1D("HLTZB", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                     1000);
    HLT40 = new TH1D("HLT40", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                     1000);
    HLT60 = new TH1D("HLT60", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                     1000);
    HLT80 = new TH1D("HLT80", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                     1000);
    HLT100 = new TH1D("HLT100", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                      1000);
    HLT120 = new TH1D("HLT120", "leading jet p_{T}; leading jet p_{T};", 100, 0,
                      1000);

    HLTZB_ptav = new TH1D("HLTZB_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
    HLT40_ptav = new TH1D("HLT40_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
    HLT60_ptav = new TH1D("HLT60_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
    HLT80_ptav = new TH1D("HLT80_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
    HLT100_ptav =
        new TH1D("HLT100_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);
    HLT120_ptav =
        new TH1D("HLT120_ptav", "p_{T,avg}; p_{T,avg};", 100, 0, 1000);

    // Dijet asymmetry histograms
    dijetasymmetry = new TH1D("dijetasymmetry", "  ; asymmetry;", 40, 0, 1);
    dijetasymmetry_now =
        new TH1D("dijetasymmetry_now", "  ; asymmetry;", 40, 0, 1);
    dijetdeltaphi =
        new TH1D("dijetdeltaphi", " ; delta phi;", 40, 0, 3.1415926535);
    dijetdeltaeta = new TH1D("dijetdeltaeta", " ; delta eta;", 40, 0, 5.2);

    dijetasymmetry_a01 =
        new TProfile("dijetasymmetry_a01", "  ; ;", nptforjec, &ptforjec[0]);
    dijetasymmetry_a02 =
        new TProfile("dijetasymmetry_a02", "  ; ;", nptforjec, &ptforjec[0]);
    dijetbalance_a03 = new TH1D("dijetbalance_a03", "  ; ;", 20, -2, 2);
    dijetasymmetry_a03 =
        new TProfile("dijetasymmetry_a03", "  ; ;", nptforjec, &ptforjec[0]);
    dijetasymmetry_a035 =
        new TProfile("dijetasymmetry_a035", "  ; ;", nptforjec, &ptforjec[0]);
    dijetasymmetry_a04 =
        new TProfile("dijetasymmetry_a04", "  ; ;", nptforjec, &ptforjec[0]);
    dijetasymmetry_a05 =
        new TProfile("dijetasymmetry_a05", "  ; ;", nptforjec, &ptforjec[0]);
    dijetasymmetry_a06 =
        new TProfile("dijetasymmetry_a06", "  ; ;", nptforjec, &ptforjec[0]);
    dijetbalance_a1 = new TH1D("dijetbalance_a1", "  ; ;", 20, -2, 2);
    dijetasymmetry_a1 =
        new TProfile("dijetasymmetry_a1", "  ; ;", nptforjec, &ptforjec[0]);

    // Dijet 3D profiles (only for wide eta bin)
    if ((this->etamin - this->etamax) < -10) {
      dijetasymmetry3D =
          new TProfile3D("dijetasymmetry3D", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0], nalphavalues, &alphavalues[0]);
      dijetasymmetry3Dabseta = new TProfile3D(
          "dijetasymmetry3Dabseta", ";;", nptforjec, &ptforjec[0], nwabsetas,
          &wabsetarange[0], nalphavalues, &alphavalues[0]);
      dijetasymmetry3Dabsetawide = new TProfile3D(
          "dijetasymmetry3Dabsetawide", ";;", nptforjec, &ptforjec[0],
          ndwabsetas, &dwabsetarange[0], nalphavalues, &alphavalues[0]);
      dijetasymmetry3Dnarrow = new TProfile3D(
          "dijetasymmetry3Dnarrow", ";;", nptforjec, &ptforjec[0], netas,
          &etarange[0], nalphavalues, &alphavalues[0]);
      dijetasymmetry3Dabsetanarrow = new TProfile3D(
          "dijetasymmetry3Dabsetanarrow", ";;", nptforjec, &ptforjec[0],
          nabsetas, &absetarange[0], nalphavalues, &alphavalues[0]);
      dijetasymmetry2D_a01 =
          new TProfile2D("dijetasymmetry2D_a01", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
      dijetasymmetry2D_a02 =
          new TProfile2D("dijetasymmetry2D_a02", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
      dijetasymmetry2D_a03 =
          new TProfile2D("dijetasymmetry2D_a03", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
      dijetasymmetry2D_a04 =
          new TProfile2D("dijetasymmetry2D_a04", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
      dijetasymmetry2D_a05 =
          new TProfile2D("dijetasymmetry2D_a05", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
      dijetasymmetry2D_a06 =
          new TProfile2D("dijetasymmetry2D_a06", ";;", nptforjec, &ptforjec[0],
                         nwetas, &wetarange[0]);
    }

    // JER histograms (dijet-specific)
    vector<double> x(61);
    for (unsigned int i = 0; i != x.size(); ++i)
      x[i] = 0.05 * i;
    const int nx = x.size() - 1;

    vector<double> y(161);
    for (unsigned int i = 0; i != y.size(); ++i)
      y[i] = -1 + 0.0125 * i;
    const int ny = y.size() - 1;

    if (ismc) {
      responses3D = new TH3D("responses3D", ";;", nptforJER, &ptforJER[0],
                             netaforjer, &etaforjer[0], nx, &x[0]);
      phiresponse = new TH3D("phiresponses3D", ";;", nptforJER, &ptforJER[0],
                             netaforjer, &etaforjer[0], ny, &y[0]);
      etaresponse = new TH3D("etaresponses3D", ";;", nptforJER, &ptforJER[0],
                             netaforjer, &etaforjer[0], ny, &y[0]);
    }

    absasymmdist3D = new TH3D("absasymmdist3D", ";;", nptforjec, &ptforjec[0],
                              njeretas, &jeretarange[0], nx, &x[0]);
    asymmdist3D = new TH3D("asymmdist3D", ";;", nptforjec, &ptforjec[0],
                           njeretas, &jeretarange[0], ny, &y[0]);

    // JER needs asymmetries as function of alpha
    asymmdist3D_a10 = new TH3D("asymmdist3D_a10", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a10 =
        new TH3D("absasymmdist3D_a10", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a15 = new TH3D("asymmdist3D_a15", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a15 =
        new TH3D("absasymmdist3D_a15", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a20 = new TH3D("asymmdist3D_a20", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a20 =
        new TH3D("absasymmdist3D_a20", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a25 = new TH3D("asymmdist3D_a25", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a25 =
        new TH3D("absasymmdist3D_a25", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a30 = new TH3D("asymmdist3D_a30", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a30 =
        new TH3D("absasymmdist3D_a30", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a35 = new TH3D("asymmdist3D_a35", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a35 =
        new TH3D("absasymmdist3D_a35", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a40 = new TH3D("asymmdist3D_a40", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a40 =
        new TH3D("absasymmdist3D_a40", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);
    asymmdist3D_a45 = new TH3D("asymmdist3D_a45", ";;", nptforjec, &ptforjec[0],
                               njeretas, &jeretarange[0], ny, &y[0]);
    absasymmdist3D_a45 =
        new TH3D("absasymmdist3D_a45", ";;", nptforjec, &ptforjec[0], njeretas,
                 &jeretarange[0], nx, &x[0]);

    if (ismc) {
      absasymmdist3D_gen_a10 =
          new TH3D("absasymmdist3D_gen_a10", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a15 =
          new TH3D("absasymmdist3D_gen_a15", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a20 =
          new TH3D("absasymmdist3D_gen_a20", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a25 =
          new TH3D("absasymmdist3D_gen_a25", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a30 =
          new TH3D("absasymmdist3D_gen_a30", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a35 =
          new TH3D("absasymmdist3D_gen_a35", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a40 =
          new TH3D("absasymmdist3D_gen_a40", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
      absasymmdist3D_gen_a45 =
          new TH3D("absasymmdist3D_gen_a45", ";;", nptforjec, &ptforjec[0],
                   njeretas, &jeretarange[0], nx, &x[0]);
    }
  }

  // ============================================
  // PHOTON+JET SPECIFIC HISTOGRAMS
  // ============================================
  if (createPhotonJet) {
    photon_pt = new TH1D("photon_pt", "Photon p_{T}", 100, 0, 500);
    photon_eta = new TH1D("photon_eta", "Photon #eta", 50, -2.5, 2.5);
    photon_phi = new TH1D("photon_phi", "Photon #phi", nphis, phirange);
    photon_HoverE = new TH1D("photon_HoverE", "Photon H/E", 50, 0, 0.5);
    photon_sigmaIetaIeta = new TH1D("photon_sigmaIetaIeta",
                                    "Photon #sigma_{i#etai#eta}", 50, 0, 0.05);
    photon_SwissCrx =
        new TH1D("photon_SwissCrx", "Photon Swiss Cross", 50, 0.8, 1.0);
    photon_SeedTime =
        new TH1D("photon_SeedTime", "Photon Seed Time", 50, -10, 10);

    awayside_jet_pt =
        new TH1D("awayside_jet_pt", "Away-side Jet p_{T}", 100, 0, 500);
    awayside_jet_eta =
        new TH1D("awayside_jet_eta", "Away-side Jet #eta", netas, etarange);
    awayside_jet_phi =
        new TH1D("awayside_jet_phi", "Away-side Jet #phi", nphis, phirange);
    awayside_jet_uncorr_pt = new TH1D(
        "awayside_jet_uncorr_pt", "Away-side Jet p_{T} (uncorr)", 100, 0, 500);

    photonjet_dphi =
        new TH1D("photonjet_dphi", "#Delta#phi(#gamma,jet)", 50, 0, 3.14159);
    photonjet_balance =
        new TH1D("photonjet_balance", "p_{T}^{jet}/p_{T}^{#gamma}", 50, 0, 2);
    photonjet_ptavg = new TH1D("photonjet_ptavg", "p_{T,avg}", 100, 0, 500);
    photonjet_alpha =
        new TH1D("photonjet_alpha", "#alpha (3rd jet fraction)", 50, 0, 1);

    // 1D profiles
    photonjet_balance_a01 =
        new TProfile("photonjet_balance_a01", "Balance (#alpha<0.1)", nptforjec,
                     &ptforjec[0]);
    photonjet_balance_a02 =
        new TProfile("photonjet_balance_a02", "Balance (#alpha<0.2)", nptforjec,
                     &ptforjec[0]);
    photonjet_balance_a03 =
        new TProfile("photonjet_balance_a03", "Balance (#alpha<0.3)", nptforjec,
                     &ptforjec[0]);
    photonjet_balance_a04 =
        new TProfile("photonjet_balance_a04", "Balance (#alpha<0.4)", nptforjec,
                     &ptforjec[0]);
    photonjet_balance_a05 =
        new TProfile("photonjet_balance_a05", "Balance (#alpha<0.5)", nptforjec,
                     &ptforjec[0]);
    photonjet_balance_a06 =
        new TProfile("photonjet_balance_a06", "Balance (#alpha<0.6)", nptforjec,
                     &ptforjec[0]);

    // 2D profiles
    photonjet_balance2D_a01 =
        new TProfile2D("photonjet_balance2D_a01", "Balance (#alpha<0.1)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    photonjet_balance2D_a02 =
        new TProfile2D("photonjet_balance2D_a02", "Balance (#alpha<0.2)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    photonjet_balance2D_a03 =
        new TProfile2D("photonjet_balance2D_a03", "Balance (#alpha<0.3)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    photonjet_balance2D_a04 =
        new TProfile2D("photonjet_balance2D_a04", "Balance (#alpha<0.4)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    photonjet_balance2D_a05 =
        new TProfile2D("photonjet_balance2D_a05", "Balance (#alpha<0.5)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);
    photonjet_balance2D_a06 =
        new TProfile2D("photonjet_balance2D_a06", "Balance (#alpha<0.6)",
                       nptforjec, &ptforjec[0], nwetas, &wetarange[0]);

    // Photon trigger
    HLTPhoton30 = new TH1D("HLTPhoton30", "HLT Photon30 events", 2, 0, 2);
    HLTPhoton30_ptav =
        new TH1D("HLTPhoton30_ptav", "p_{T,avg} (HLT Photon30)", 100, 0, 500);

    // 3D profiles (only for wide eta bin)
    if ((this->etamin - this->etamax) < -10) {
      photonjet_balance3D = new TProfile3D(
          "photonjet_balance3D", "Balance vs p_{T,avg}, #eta_{jet}, #alpha",
          nptforjec, &ptforjec[0], nwetas, &wetarange[0], nalphavalues,
          &alphavalues[0]);
      photonjet_balance3Dwide = new TProfile3D(
          "photonjet_balance3Dwide", "Balance (wide #eta bins)", nptforjec,
          &ptforjec[0], nwetas, &wetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dnarrow = new TProfile3D(
          "photonjet_balance3Dnarrow", "Balance (narrow #eta bins)", nptforjec,
          &ptforjec[0], netas, &etarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabseta = new TProfile3D(
          "photonjet_balance3Dabseta",
          "Balance vs p_{T,avg}, |#eta_{jet}|, #alpha", nptforjec, &ptforjec[0],
          nwabsetas, &wabsetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabsetawide = new TProfile3D(
          "photonjet_balance3Dabsetawide", "Balance (wide |#eta| bins)",
          nptforjec, &ptforjec[0], ndwabsetas, &dwabsetarange[0], nalphavalues,
          &alphavalues[0]);
      photonjet_balance3Dabsetanarrow = new TProfile3D(
          "photonjet_balance3Dabsetanarrow", "Balance (narrow |#eta| bins)",
          nptforjec, &ptforjec[0], nabsetas, &absetarange[0], nalphavalues,
          &alphavalues[0]);
      photonjet_balance3D_jetpt = new TProfile3D(
          "photonjet_balance3D_jetpt",
          "Balance vs p_{T}^{jet}, #eta_{jet}, #alpha", nptforjec, &ptforjec[0],
          nwetas, &wetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dwide_jetpt = new TProfile3D(
          "photonjet_balance3Dwide_jetpt",
          "Balance vs p_{T}^{jet} (wide #eta bins)", nptforjec, &ptforjec[0],
          nwetas, &wetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dnarrow_jetpt = new TProfile3D(
          "photonjet_balance3Dnarrow_jetpt",
          "Balance vs p_{T}^{jet} (narrow #eta bins)", nptforjec, &ptforjec[0],
          netas, &etarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabseta_jetpt =
          new TProfile3D("photonjet_balance3Dabseta_jetpt",
                         "Balance vs p_{T}^{jet}, |#eta_{jet}|, #alpha",
                         nptforjec, &ptforjec[0], nwabsetas, &wabsetarange[0],
                         nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabsetawide_jetpt = new TProfile3D(
          "photonjet_balance3Dabsetawide_jetpt",
          "Balance vs p_{T}^{jet} (wide |#eta| bins)", nptforjec, &ptforjec[0],
          ndwabsetas, &dwabsetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabsetanarrow_jetpt =
          new TProfile3D("photonjet_balance3Dabsetanarrow_jetpt",
                         "Balance vs p_{T}^{jet} (narrow |#eta| bins)",
                         nptforjec, &ptforjec[0], nabsetas, &absetarange[0],
                         nalphavalues, &alphavalues[0]);
      photonjet_balance3D_counts = new TH3D(
          "photonjet_balance3D_counts",
          "Entries vs p_{T,avg}, #eta_{jet}, #alpha", nptforjec, &ptforjec[0],
          nwetas, &wetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dwide_counts =
          new TH3D("photonjet_balance3Dwide_counts", "Entries (wide #eta bins)",
                   nptforjec, &ptforjec[0], nwetas, &wetarange[0], nalphavalues,
                   &alphavalues[0]);
      photonjet_balance3Dnarrow_counts =
          new TH3D("photonjet_balance3Dnarrow_counts",
                   "Entries (narrow #eta bins)", nptforjec, &ptforjec[0], netas,
                   &etarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabseta_counts = new TH3D(
          "photonjet_balance3Dabseta_counts",
          "Entries vs p_{T,avg}, |#eta_{jet}|, #alpha", nptforjec, &ptforjec[0],
          nwabsetas, &wabsetarange[0], nalphavalues, &alphavalues[0]);
      photonjet_balance3Dabsetawide_counts = new TH3D(
          "photonjet_balance3Dabsetawide_counts", "Entries (wide |#eta| bins)",
          nptforjec, &ptforjec[0], ndwabsetas, &dwabsetarange[0], nalphavalues,
          &alphavalues[0]);
      photonjet_balance3Dabsetanarrow_counts =
          new TH3D("photonjet_balance3Dabsetanarrow_counts",
                   "Entries (narrow |#eta| bins)", nptforjec, &ptforjec[0],
                   nabsetas, &absetarange[0], nalphavalues, &alphavalues[0]);

      // Balance distribution: (photon_pT, alpha, balance_value)
      // Match the standard binning used elsewhere (histograms.h):
      //   X: ptforjec (variable bins)
      //   Y: alphavalues (variable bins)
      //   Z: balance (50 uniform bins from 0 to 2)
      // NOTE: In this ROOT build, TH3D does not provide a constructor for
      //       (variable x bins, variable y bins, uniform z range). We therefore
      //       pass an explicit z-edge array with uniform spacing (0..2).
      static double balanceEdges[51];
      static bool balanceEdgesInit = false;
      if (!balanceEdgesInit) {
        for (int i = 0; i <= 50; ++i) {
          balanceEdges[i] = 0.0 + (2.0 / 50.0) * i;
        }
        balanceEdgesInit = true;
      }

      photonjet_balance_dist =
          new TH3D("photonjet_balance_dist",
                   "Balance distribution;p_{T}^{#gamma} "
                   "(GeV);#alpha;p_{T}^{jet}/p_{T}^{#gamma}",
                   nptforjec, &ptforjec[0], nalphavalues, &alphavalues[0], 50,
                   balanceEdges);
    }
  }
}

void histograms::Write() {
  dir->cd();
  dir->Write();
}
// HistosBasic::~HistosBasic() {
//  Write();/
// };
