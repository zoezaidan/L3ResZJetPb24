// Derive L3 residual corrections from photon+jet balancing using 3D profiles
// Balance = jet_pT / photon_pT
// L3 Residual = balance_MC / balance_data (corrects data to match MC response)
// 
// This script produces:
// 1. Ratio of pT response (Data/MC) as functions of:
//    - Photon pT (collapsed over eta)
//    - Leading jet pT (derived from balance = jet_pT/photon_pT)
// 2. Individual Data and MC histograms for each quantity
// 3. All calculations for all alpha bins

#include <iostream>
#include <fstream>
#include "../fillhistograms/histograms.h"

void deriveL3_from_photonjet(
    TString mcFile = "/eos/home-z/zzaidanc/JetMinPOG/L3ResZJet/ZJETMC_AK4_zjet_ak4.root",
    TString dataFile = "/eos/home-z/zzaidanc/JetMinPOG/L3ResZJet/ZJETSM_AK4_zjet_ak4.root",
    TString outfilename = "L3Residual_ZJet.root",
    bool dodt = true,
    int alphabin = 5,
    bool useabs = true,
    bool usewideabs = false) {

  // Open MC file
  TFile *inFileMC = TFile::Open(mcFile);
  if (!inFileMC || inFileMC->IsZombie()) {
    cout << "ERROR: Cannot open MC file: " << mcFile << endl;
    return;
  }

  // Open Data file
  TFile *inFileDT = TFile::Open(dataFile);
  if (!inFileDT || inFileDT->IsZombie()) {
    cout << "ERROR: Cannot open Data file: " << dataFile << endl;
    return;
  }

  // These are bins to be processed
  vector<string> etabins = {"eta_-5.2_5.2"};
  int i = 0;

  map<string, TProfile3D*> mc3d, data3d;
  map<string, TH3D*> counts_mc3d, counts_data3d;
  map<string, TH1D*> responses;
  map<int, TH1D*> respETA;

  TFile *outfile = new TFile(outfilename, "RECREATE");

  // Histograms for alpha dependence study
  TH1D* vsalpha_mc = new TH1D("vsalpha_mc", "MC balance vs alpha; alpha; balance", histograms::nalphavalues, &histograms::alphavalues[0]);
  TH1D* vsalpha_data = new TH1D("vsalpha_data", "Data balance vs alpha; alpha; balance", histograms::nalphavalues, &histograms::alphavalues[0]);

  // Histograms for eta dependence - declare pointers
  TH1D* vseta_mc(0);
  TH1D* vseta_data(0);
  TH1D* aerrormc(0);
  TH1D* aerrordt(0);

  // Get the appropriate 3D balance profiles based on eta binning choice
  if (useabs) {
    vseta_mc = new TH1D("vseta_mc", "MC balance; |#eta_{jet}|; balance", histograms::nwabsetas, &histograms::wabsetarange[0]);
    vseta_data = new TH1D("vseta_data", "Data balance; |#eta_{jet}|; balance", histograms::nwabsetas, &histograms::wabsetarange[0]);
    aerrormc = new TH1D("aerrormc", "MC error; |#eta_{jet}|; error", histograms::nwabsetas, &histograms::wabsetarange[0]);
    aerrordt = new TH1D("aerrordt", "Data error; |#eta_{jet}|; error", histograms::nwabsetas, &histograms::wabsetarange[0]);

    mc3d[etabins[i].c_str()] = (TProfile3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabseta");
    data3d[etabins[i].c_str()] = (TProfile3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabseta");
    counts_mc3d[etabins[i].c_str()] = (TH3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabseta_counts");
    counts_data3d[etabins[i].c_str()] = (TH3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabseta_counts");
  }
  else if (usewideabs) {
    vseta_mc = new TH1D("vseta_mc", "MC balance; |#eta_{jet}|; balance", histograms::ndwabsetas, &histograms::dwabsetarange[0]);
    vseta_data = new TH1D("vseta_data", "Data balance; |#eta_{jet}|; balance", histograms::ndwabsetas, &histograms::dwabsetarange[0]);
    aerrormc = new TH1D("aerrormc", "MC error; |#eta_{jet}|; error", histograms::ndwabsetas, &histograms::dwabsetarange[0]);
    aerrordt = new TH1D("aerrordt", "Data error; |#eta_{jet}|; error", histograms::ndwabsetas, &histograms::dwabsetarange[0]);

    mc3d[etabins[i].c_str()] = (TProfile3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabsetawide");
    data3d[etabins[i].c_str()] = (TProfile3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabsetawide");
    counts_mc3d[etabins[i].c_str()] = (TH3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabsetawide_counts");
    counts_data3d[etabins[i].c_str()] = (TH3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3Dabsetawide_counts");
  }
  else {
    vseta_mc = new TH1D("vseta_mc", "MC balance; #eta_{jet}; balance", histograms::nwetas, &histograms::wetarange[0]);
    vseta_data = new TH1D("vseta_data", "Data balance; #eta_{jet}; balance", histograms::nwetas, &histograms::wetarange[0]);
    aerrormc = new TH1D("aerrormc", "MC error; #eta_{jet}; error", histograms::nwetas, &histograms::wetarange[0]);
    aerrordt = new TH1D("aerrordt", "Data error; #eta_{jet}; error", histograms::nwetas, &histograms::wetarange[0]);

    mc3d[etabins[i].c_str()] = (TProfile3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3D");
    data3d[etabins[i].c_str()] = (TProfile3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3D");
    counts_mc3d[etabins[i].c_str()] = (TH3D*)inFileMC->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3D_counts");
    counts_data3d[etabins[i].c_str()] = (TH3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/Muonjet_balance3D_counts");
  }

  // Check that histograms were found
  if (!mc3d[etabins[i].c_str()]) {
    cout << "ERROR: Cannot find MC muon+jet balance histogram (Muonjet_balance3D*)" << endl;
    return;
  }
  if (!data3d[etabins[i].c_str()]) {
    cout << "ERROR: Cannot find Data muon+jet balance histogram (Muonjet_balance3D*)" << endl;
    return;
  }

  cout << etabins[i] << " MC nptbins: " << mc3d[etabins[i].c_str()]->GetXaxis()->GetNbins() << endl;
  cout << etabins[i] << " Data nptbins: " << data3d[etabins[i].c_str()]->GetXaxis()->GetNbins() << endl;
  cout << "Alpha cut bin: " << alphabin << endl;
  // Alpha cut: sum bins 1 to alphabin (cumulative cut alpha < threshold)
  float alphaCutValue = mc3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphabin+1);
  cout << "Alpha cut: alpha < " << alphaCutValue << endl;
  cout << "(Summing alpha bins 1 to " << alphabin << ")" << endl;

  string alphastr = Form("alpha%.2f", alphaCutValue);

  ///////////////// Balance vs eta in bins of pT (for a given alpha cut)
  // NOTE: With cumulative filling in analyse_PhotonJet.cc, each alpha bin 
  // already represents "alpha < threshold", so we read directly from alphabin

  for (int ptbin = 1; ptbin <= mc3d[etabins[i].c_str()]->GetXaxis()->GetNbins(); ++ptbin) {

    cout << "\nGetting L3 corrections as function of eta" << endl;
    cout << "pT bin edges: " << mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin)
         << " to " << mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin+1) << endl;

    string ptstr = Form("%.0fto%.0f",
                        mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin),
                        mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin+1));

    // Process MC: get balance vs eta for this pT bin at the specified alpha cut
    // With cumulative filling, alphabin directly corresponds to "alpha < threshold"
    for (int etabin = 1; etabin <= mc3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
      double val = mc3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphabin);
      double err = mc3d[etabins[i].c_str()]->GetBinError(ptbin, etabin, alphabin);

      vseta_mc->SetBinContent(etabin, val);
      vseta_mc->SetBinError(etabin, (TMath::IsNaN(err) ? 0. : err));
      aerrormc->SetBinContent(etabin, (TMath::IsNaN(err) ? 0. : err));
    }

    // Process Data: get balance vs eta for this pT bin at the specified alpha cut
    for (int etabin = 1; etabin <= data3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
      double val = data3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphabin);
      double err = data3d[etabins[i].c_str()]->GetBinError(ptbin, etabin, alphabin);

      vseta_data->SetBinContent(etabin, val);
      vseta_data->SetBinError(etabin, (TMath::IsNaN(err) ? 0. : err));
      aerrordt->SetBinContent(etabin, (TMath::IsNaN(err) ? 0. : err));
    }

    // Clone and save MC balance
    TH1D* balance_eta_mc = (TH1D*)vseta_mc->Clone(Form("mc_balance_pt%s_%s", ptstr.c_str(), alphastr.c_str()));
    balance_eta_mc->SetLineColor(kBlue);

    // Clone and save Data balance
    TH1D* balance_eta_data = (TH1D*)vseta_data->Clone(Form("dt_balance_pt%s_%s", ptstr.c_str(), alphastr.c_str()));
    balance_eta_data->SetLineColor(kRed);

    // L3 Residual correction = balance_MC / balance_data
    // This corrects data jets to match MC response
    // If data jets are lower (balance_data < balance_mc), correction > 1
    TH1D* l3res_eta = (TH1D*)balance_eta_mc->Clone(Form("l3res_pt%s_%s", ptstr.c_str(), alphastr.c_str()));
    l3res_eta->Divide(balance_eta_data);

    // Propagate errors properly
    for (int bin = 1; bin <= l3res_eta->GetXaxis()->GetNbins(); ++bin) {
      double mc_val = balance_eta_mc->GetBinContent(bin);
      double mc_err = balance_eta_mc->GetBinError(bin);
      double dt_val = balance_eta_data->GetBinContent(bin);
      double dt_err = balance_eta_data->GetBinError(bin);

      if (dt_val > 0 && mc_val > 0) {
        double ratio = mc_val / dt_val;
        double rel_err = sqrt(pow(mc_err/mc_val, 2) + pow(dt_err/dt_val, 2));
        l3res_eta->SetBinError(bin, ratio * rel_err);
      } else {
        l3res_eta->SetBinContent(bin, 1.0);
        l3res_eta->SetBinError(bin, 0.0);
      }
    }

    // Store for later use
    responses[Form("mc_pt%s_%s", ptstr.c_str(), alphastr.c_str())] = balance_eta_mc;
    responses[Form("dt_pt%s_%s", ptstr.c_str(), alphastr.c_str())] = balance_eta_data;
    respETA[ptbin] = (TH1D*)l3res_eta->Clone(Form("l3res_ptbin%d", ptbin));

    // Write histograms
    balance_eta_mc->Write();
    balance_eta_data->Write();
    l3res_eta->Write();
    respETA[ptbin]->Write(Form("ratio_pt%s_%s", ptstr.c_str(), alphastr.c_str()));
  }

  ///////////////// Balance vs alpha in bins of pT and eta (for ISR/FSR studies)

  cout << "\n\nNumber of alpha bins: " << mc3d[etabins[i].c_str()]->GetZaxis()->GetNbins() << endl;

  for (int ptbin = 1; ptbin <= mc3d[etabins[i].c_str()]->GetXaxis()->GetNbins(); ++ptbin) {
    cout << "\nNEW PT BIN " << ptbin << endl;

    for (int etabin = 1; etabin <= data3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
      cout << "NEW ETA BIN " << data3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin)
           << " to " << data3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin+1) << endl;

      // Loop over alpha bins
      for (int abin = 1; abin <= mc3d[etabins[i].c_str()]->GetZaxis()->GetNbins(); ++abin) {
        // MC
        double val_mc = mc3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, abin);
        double err_mc = mc3d[etabins[i].c_str()]->GetBinError(ptbin, etabin, abin);
        vsalpha_mc->SetBinContent(abin, val_mc);
        vsalpha_mc->SetBinError(abin, (TMath::IsNaN(err_mc) ? 0. : err_mc));

        // Data
        double val_dt = data3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, abin);
        double err_dt = data3d[etabins[i].c_str()]->GetBinError(ptbin, etabin, abin);
        vsalpha_data->SetBinContent(abin, val_dt);
        vsalpha_data->SetBinError(abin, (TMath::IsNaN(err_dt) ? 0. : err_dt));
      }

      // Write balance vs alpha for MC and Data
      vsalpha_mc->Write(Form("Balance_vsa_mc_%d_%d", ptbin, etabin));
      vsalpha_data->Write(Form("Balance_vsa_data_%d_%d", ptbin, etabin));

      // Compute ratio (L3Res vs alpha)
      TH1D* vsalpha_ratio = (TH1D*)vsalpha_mc->Clone(Form("L3Res_vsa_%d_%d", ptbin, etabin));
      vsalpha_ratio->Divide(vsalpha_data);
      vsalpha_ratio->Write();

      // Normalize to the chosen alpha bin value
      TH1D* vsalpha_norm = (TH1D*)vsalpha_ratio->Clone(Form("L3Res_vsa_norm_%d_%d", ptbin, etabin));
      double norm = respETA[ptbin]->GetBinContent(etabin);
      if (norm > 0) {
        for (int bin = 1; bin <= vsalpha_norm->GetXaxis()->GetNbins(); ++bin) {
          double val = vsalpha_norm->GetBinContent(bin);
          double err = vsalpha_norm->GetBinError(bin);
          vsalpha_norm->SetBinContent(bin, val / norm);
          vsalpha_norm->SetBinError(bin, err / norm);
        }
      }
      vsalpha_norm->Write();
    }
  }

  ///////////////// Balance vs Photon pT (collapsed over eta) for alpha cuts

  cout << "\n\n===== Balance vs Photon pT for cumulative alpha cuts =====" << endl;

  int nPtBins = mc3d[etabins[i].c_str()]->GetXaxis()->GetNbins();
  int nEtaBins = mc3d[etabins[i].c_str()]->GetYaxis()->GetNbins();
  int nAlphaBins = mc3d[etabins[i].c_str()]->GetZaxis()->GetNbins();

  // Get binning from the 3D profile for photon pT (X-axis)
  TAxis* ptaxis = mc3d[etabins[i].c_str()]->GetXaxis();

  // Extract variable bin edges from the pT axis
  std::vector<double> ptBinEdges(nPtBins + 1);
  for (int bin = 1; bin <= nPtBins + 1; ++bin) {
    ptBinEdges[bin-1] = ptaxis->GetBinLowEdge(bin);
  }

  // Create histograms for photon pT dependence for each cumulative alpha cut
  map<int, TH1D*> balance_vsphotonpt_mc;
  map<int, TH1D*> balance_vsphotonpt_data;
  map<int, TH1D*> ratio_vsphotonpt;

  // Reference alpha bin MC/Data ratio (for normalization to extract kFSR)
  // Will be filled first, then used to normalize other alpha bins
  TH1D* ratio_ref = nullptr;
  
  for (int alphaCutBin = 1; alphaCutBin <= nAlphaBins; ++alphaCutBin) {
    // Cumulative alpha cut: alpha < upper edge of this bin
    float alpha_cut = mc3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphaCutBin+1);

    // Use variable binning from the input histogram
    TH1D* h_mc = new TH1D(Form("balance_vsphotonpt_mc_alpha%d", alphaCutBin),
                          Form("MC Balance vs Photon pT (#alpha < %.2f);Photon p_{T} (GeV);Balance",
                               alpha_cut),
                          nPtBins, ptBinEdges.data());
    h_mc->SetLineColor(kBlue);
    h_mc->SetMarkerColor(kBlue);

    TH1D* h_data = new TH1D(Form("balance_vsphotonpt_data_alpha%d", alphaCutBin),
                            Form("Data Balance vs Photon pT (#alpha < %.2f);Photon p_{T} (GeV);Balance",
                                 alpha_cut),
                            nPtBins, ptBinEdges.data());
    h_data->SetLineColor(kRed);
    h_data->SetMarkerColor(kRed);

    // Collapse over eta bins, reading from alphaCutBin directly (cumulative fill)
    // With cumulative filling, bin N already represents "alpha < threshold_N"
    for (int ptbin = 1; ptbin <= nPtBins; ++ptbin) {
      double sum_mc = 0., sum_mc_entries = 0.;
      double sum_data = 0., sum_data_entries = 0.;

      for (int etabin = 1; etabin <= nEtaBins; ++etabin) {
        // Read directly from alphaCutBin (cumulative fill, no need to sum)
        double val_mc = mc3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphaCutBin);
        double entries_mc = mc3d[etabins[i].c_str()]->GetBinEntries(
            mc3d[etabins[i].c_str()]->GetBin(ptbin, etabin, alphaCutBin));
        double val_data = data3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphaCutBin);
        double entries_data = data3d[etabins[i].c_str()]->GetBinEntries(
            data3d[etabins[i].c_str()]->GetBin(ptbin, etabin, alphaCutBin));

        if (val_mc > 0 && entries_mc > 0 && !TMath::IsNaN(val_mc)) {
          sum_mc += val_mc * entries_mc;
          sum_mc_entries += entries_mc;
        }
        if (val_data > 0 && entries_data > 0 && !TMath::IsNaN(val_data)) {
          sum_data += val_data * entries_data;
          sum_data_entries += entries_data;
        }
      }

      if (sum_mc_entries > 0) {
        double avg_mc = sum_mc / sum_mc_entries;
        double err_mc = avg_mc / sqrt(sum_mc_entries);  // Statistical error estimate
        h_mc->SetBinContent(ptbin, avg_mc);
        h_mc->SetBinError(ptbin, err_mc);
      }
      if (sum_data_entries > 0) {
        double avg_data = sum_data / sum_data_entries;
        double err_data = avg_data / sqrt(sum_data_entries);  // Statistical error estimate
        h_data->SetBinContent(ptbin, avg_data);
        h_data->SetBinError(ptbin, err_data);
      }
    }

    balance_vsphotonpt_mc[alphaCutBin] = h_mc;
    balance_vsphotonpt_data[alphaCutBin] = h_data;

    h_mc->Write();
    h_data->Write();

    // Compute ratio
    TH1D* h_ratio = (TH1D*)h_mc->Clone(Form("ratio_vsphotonpt_alpha%d", alphaCutBin));
    h_ratio->Divide(h_data);
    h_ratio->SetTitle(Form("Balance Ratio (MC/Data) vs Photon pT (#alpha < %.2f)", alpha_cut));
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);

    // Propagate errors
    for (int bin = 1; bin <= h_ratio->GetXaxis()->GetNbins(); ++bin) {
      double mc_val = h_mc->GetBinContent(bin);
      double mc_err = h_mc->GetBinError(bin);
      double dt_val = h_data->GetBinContent(bin);
      double dt_err = h_data->GetBinError(bin);

      if (dt_val > 0 && mc_val > 0) {
        double ratio = mc_val / dt_val;
        double rel_err = sqrt(pow(mc_err/mc_val, 2) + pow(dt_err/dt_val, 2));
        h_ratio->SetBinError(bin, ratio * rel_err);
      } else {
        h_ratio->SetBinContent(bin, 1.0);
        h_ratio->SetBinError(bin, 0.0);
      }
    }

    ratio_vsphotonpt[alphaCutBin] = h_ratio;
    h_ratio->Write();
    
    // Store the reference alpha bin ratio for normalization
    if (alphaCutBin == alphabin) {
      ratio_ref = (TH1D*)h_ratio->Clone("ratio_ref_vsphotonpt");
      ratio_ref->Write();
    }
  }
  
  // Second pass: create normalized ratios (divide by reference alpha bin)
  // This is for kFSR extraction: normalized_ratio = (MC/Data at α) / (MC/Data at ref α)
  if (ratio_ref) {
    cout << "\n===== Creating normalized ratios for kFSR extraction =====" << endl;
    cout << "Reference alpha bin: " << alphabin << " (alpha < " << alphaCutValue << ")" << endl;
    
    for (int alphaCutBin = 1; alphaCutBin <= nAlphaBins; ++alphaCutBin) {
      float alpha_cut = mc3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphaCutBin+1);
      
      TH1D* h_ratio_norm = (TH1D*)ratio_vsphotonpt[alphaCutBin]->Clone(
          Form("ratio_norm_vsphotonpt_alpha%d", alphaCutBin));
      h_ratio_norm->SetTitle(Form("Normalized Ratio vs Photon pT (#alpha < %.2f / ref #alpha < %.2f)", 
                                   alpha_cut, alphaCutValue));
      
      // Divide by reference
      for (int bin = 1; bin <= h_ratio_norm->GetXaxis()->GetNbins(); ++bin) {
        double val = h_ratio_norm->GetBinContent(bin);
        double err = h_ratio_norm->GetBinError(bin);
        double ref_val = ratio_ref->GetBinContent(bin);
        double ref_err = ratio_ref->GetBinError(bin);
        
        if (ref_val > 0 && val > 0) {
          double norm_val = val / ref_val;
          double rel_err = sqrt(pow(err/val, 2) + pow(ref_err/ref_val, 2));
          h_ratio_norm->SetBinContent(bin, norm_val);
          h_ratio_norm->SetBinError(bin, norm_val * rel_err);
        } else {
          h_ratio_norm->SetBinContent(bin, 0);
          h_ratio_norm->SetBinError(bin, 0);
        }
      }
      
      h_ratio_norm->Write();
      
      if (alphaCutBin == alphabin) {
        cout << "  Alpha bin " << alphaCutBin << " (reference): all values should be 1.0" << endl;
      }
    }
  }

  ///////////////// Balance vs Leading Jet pT (derived from balance = jet_pT/photon_pT)
  ///////////////// for cumulative alpha cuts

  cout << "\n===== Balance vs Leading Jet pT for cumulative alpha cuts =====" << endl;

  // For jet pT, we create histograms by reweighting based on jet pT = balance * photon pT
  // We use the same photon pT bins but relabel them as jet pT bins

  map<int, TH1D*> balance_vsjetpt_mc;
  map<int, TH1D*> balance_vsjetpt_data;
  map<int, TH1D*> ratio_vsjetpt;

  for (int alphaCutBin = 1; alphaCutBin <= nAlphaBins; ++alphaCutBin) {
    // Cumulative alpha cut: alpha < upper edge of this bin
    float alpha_cut = mc3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphaCutBin+1);

    // Create histograms with variable binning (same as photon pT bins)
    // The binning represents the derived jet pT values
    TH1D* h_mc = new TH1D(Form("balance_vsjetpt_mc_alpha%d", alphaCutBin),
                          Form("MC Balance vs Derived Jet pT (#alpha < %.2f);Derived Jet p_{T} (GeV);Balance",
                               alpha_cut),
                          nPtBins, ptBinEdges.data());
    h_mc->SetLineColor(kBlue);
    h_mc->SetMarkerColor(kBlue);

    TH1D* h_data = new TH1D(Form("balance_vsjetpt_data_alpha%d", alphaCutBin),
                            Form("Data Balance vs Derived Jet pT (#alpha < %.2f);Derived Jet p_{T} (GeV);Balance",
                                 alpha_cut),
                            nPtBins, ptBinEdges.data());
    h_data->SetLineColor(kRed);
    h_data->SetMarkerColor(kRed);

    // For each photon pT bin, collect jet pT values weighted by balance from all eta bins
    // jet_pT = balance * photon_pT, using cumulative alpha cut (read directly from alphaCutBin)
    for (int ptbin = 1; ptbin <= nPtBins; ++ptbin) {
      double photon_pt_low = ptaxis->GetBinLowEdge(ptbin);
      double photon_pt_high = ptaxis->GetBinLowEdge(ptbin + 1);
      double photon_pt_center = ptaxis->GetBinCenter(ptbin);

      double sum_mc = 0., sum_mc_entries = 0.;
      double sum_data = 0., sum_data_entries = 0.;

      // Collapse over eta bins, reading directly from alphaCutBin (cumulative fill)
      for (int etabin = 1; etabin <= nEtaBins; ++etabin) {
        double val_mc = mc3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphaCutBin);
        double entries_mc = mc3d[etabins[i].c_str()]->GetBinEntries(
            mc3d[etabins[i].c_str()]->GetBin(ptbin, etabin, alphaCutBin));
        double val_data = data3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, alphaCutBin);
        double entries_data = data3d[etabins[i].c_str()]->GetBinEntries(
            data3d[etabins[i].c_str()]->GetBin(ptbin, etabin, alphaCutBin));

        if (val_mc > 0 && entries_mc > 0 && !TMath::IsNaN(val_mc)) {
          sum_mc += val_mc * entries_mc;
          sum_mc_entries += entries_mc;
        }
        if (val_data > 0 && entries_data > 0 && !TMath::IsNaN(val_data)) {
          sum_data += val_data * entries_data;
          sum_data_entries += entries_data;
        }
      }

      // Compute weighted average balance
      if (sum_mc_entries > 0) {
        double avg_balance_mc = sum_mc / sum_mc_entries;
        double avg_jet_pt_mc = avg_balance_mc * photon_pt_center;  // Derived jet pT
        double err_mc = avg_balance_mc / sqrt(sum_mc_entries);

        h_mc->SetBinContent(ptbin, avg_balance_mc);
        h_mc->SetBinError(ptbin, err_mc);
      }

      if (sum_data_entries > 0) {
        double avg_balance_data = sum_data / sum_data_entries;
        double avg_jet_pt_data = avg_balance_data * photon_pt_center;  // Derived jet pT
        double err_data = avg_balance_data / sqrt(sum_data_entries);

        h_data->SetBinContent(ptbin, avg_balance_data);
        h_data->SetBinError(ptbin, err_data);
      }
    }

    balance_vsjetpt_mc[alphaCutBin] = h_mc;
    balance_vsjetpt_data[alphaCutBin] = h_data;

    h_mc->Write();
    h_data->Write();

    // Compute ratio
    TH1D* h_ratio = (TH1D*)h_mc->Clone(Form("ratio_vsjetpt_alpha%d", alphaCutBin));
    h_ratio->Divide(h_data);
    h_ratio->SetTitle(Form("Balance Ratio (MC/Data) vs Derived Jet pT (#alpha < %.2f)", alpha_cut));
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);

    // Propagate errors
    for (int bin = 1; bin <= h_ratio->GetXaxis()->GetNbins(); ++bin) {
      double mc_val = h_mc->GetBinContent(bin);
      double mc_err = h_mc->GetBinError(bin);
      double dt_val = h_data->GetBinContent(bin);
      double dt_err = h_data->GetBinError(bin);

      if (dt_val > 0 && mc_val > 0) {
        double ratio = mc_val / dt_val;
        double rel_err = sqrt(pow(mc_err/mc_val, 2) + pow(dt_err/dt_val, 2));
        h_ratio->SetBinError(bin, ratio * rel_err);
      } else {
        h_ratio->SetBinContent(bin, 1.0);
        h_ratio->SetBinError(bin, 0.0);
      }
    }

    ratio_vsjetpt[alphaCutBin] = h_ratio;
    h_ratio->Write();
  }

  ///////////////// Write L3 residual corrections in JEC text format

  TString txtfilename = outfilename;
  txtfilename.ReplaceAll(".root", ".txt");

  ofstream outtext(txtfilename.Data());
  outtext << "# L3 Residual Corrections from Photon+Jet Balancing" << endl;
  outtext << "# Format: {eta_min eta_max} N pt_min pt_max correction ..." << endl;
  outtext << "# Alpha cut: alpha < " << alphaCutValue << " (summing bins 1 to " << alphabin << ")" << endl;
  outtext << "# MC file: " << mcFile << endl;
  outtext << "# Data file: " << dataFile << endl;

  // Note: nEtaBins and nPtBins already defined above
  for (int etabin = 1; etabin <= nEtaBins; ++etabin) {
    float eta_min = mc3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin);
    float eta_max = mc3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin+1);

    outtext << "{" << eta_min << " " << eta_max << "} " << (2 + nPtBins * 3) << " ";

    // Get pT range for this eta bin
    float pt_min_overall = mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(1);
    float pt_max_overall = mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(nPtBins+1);
    outtext << pt_min_overall << " " << pt_max_overall << " ";

    for (int ptbin = 1; ptbin <= nPtBins; ++ptbin) {
      float pt_min = mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin);
      float pt_max = mc3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin+1);

      // Get L3 residual correction from stored histogram
      float correction = respETA[ptbin]->GetBinContent(etabin);

      // Sanity check
      if (correction < 0.5 || correction > 2.0 || TMath::IsNaN(correction)) {
        cout << "WARNING: eta [" << eta_min << ", " << eta_max
             << "], pT [" << pt_min << ", " << pt_max
             << "], L3Res = " << correction << endl;
        correction = 1.0;
      }

      outtext << pt_min << " " << pt_max << " " << correction << " ";
    }
    outtext << endl;
  }

  outtext.close();
  cout << "\nL3 residual corrections written to: " << txtfilename << endl;

  ///////////////// Create summary 2D maps (using cumulative alpha cut)

  // Extract variable bin edges from the eta axis
  TAxis* etaaxis = mc3d[etabins[i].c_str()]->GetYaxis();
  std::vector<double> etaBinEdges(nEtaBins + 1);
  for (int bin = 1; bin <= nEtaBins + 1; ++bin) {
    etaBinEdges[bin-1] = etaaxis->GetBinLowEdge(bin);
  }

  TH2D* balanceMap_mc = new TH2D("balanceMap_mc", Form("MC Balance (#alpha < %.2f);p_{T}^{#gamma} (GeV);|#eta_{jet}|;Balance", alphaCutValue),
                                  nPtBins, ptBinEdges.data(),
                                  nEtaBins, etaBinEdges.data());

  TH2D* balanceMap_data = new TH2D("balanceMap_data", Form("Data Balance (#alpha < %.2f);p_{T}^{#gamma} (GeV);|#eta_{jet}|;Balance", alphaCutValue),
                                    nPtBins, ptBinEdges.data(),
                                    nEtaBins, etaBinEdges.data());

  TH2D* l3resMap = new TH2D("l3resMap", Form("L3 Residual (MC/Data) (#alpha < %.2f);p_{T}^{#gamma} (GeV);|#eta_{jet}|;L3Res", alphaCutValue),
                            nPtBins, ptBinEdges.data(),
                            nEtaBins, etaBinEdges.data());

  for (int etabin = 1; etabin <= nEtaBins; ++etabin) {
    for (int ptbin = 1; ptbin <= nPtBins; ++ptbin) {
      // Sum over alpha bins 1 to alphabin (cumulative cut)
      double sum_mc = 0., entries_mc = 0.;
      double sum_dt = 0., entries_dt = 0.;
      
      for (int abin = 1; abin <= alphabin; ++abin) {
        double val_mc = mc3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, abin);
        double ent_mc = mc3d[etabins[i].c_str()]->GetBinEntries(
            mc3d[etabins[i].c_str()]->GetBin(ptbin, etabin, abin));
        double val_dt = data3d[etabins[i].c_str()]->GetBinContent(ptbin, etabin, abin);
        double ent_dt = data3d[etabins[i].c_str()]->GetBinEntries(
            data3d[etabins[i].c_str()]->GetBin(ptbin, etabin, abin));
        
        if (val_mc > 0 && ent_mc > 0 && !TMath::IsNaN(val_mc)) {
          sum_mc += val_mc * ent_mc;
          entries_mc += ent_mc;
        }
        if (val_dt > 0 && ent_dt > 0 && !TMath::IsNaN(val_dt)) {
          sum_dt += val_dt * ent_dt;
          entries_dt += ent_dt;
        }
      }
      
      double mc_val = (entries_mc > 0) ? sum_mc / entries_mc : 0.;
      double dt_val = (entries_dt > 0) ? sum_dt / entries_dt : 0.;

      balanceMap_mc->SetBinContent(ptbin, etabin, mc_val);
      balanceMap_data->SetBinContent(ptbin, etabin, dt_val);

      if (dt_val > 0) {
        l3resMap->SetBinContent(ptbin, etabin, mc_val / dt_val);
      } else {
        l3resMap->SetBinContent(ptbin, etabin, 1.0);
      }
    }
  }

  balanceMap_mc->Write();
  balanceMap_data->Write();
  l3resMap->Write();

  // Save original 3D histograms for reference
  mc3d[etabins[i].c_str()]->Write("balance3D_mc");
  data3d[etabins[i].c_str()]->Write("balance3D_data");

  // Save counts histograms if available
  if (counts_mc3d[etabins[i].c_str()])
    counts_mc3d[etabins[i].c_str()]->Write("counts3D_mc");
  if (counts_data3d[etabins[i].c_str()])
    counts_data3d[etabins[i].c_str()]->Write("counts3D_data");

  outfile->Close();
  cout << "Output ROOT file: " << outfilename << endl;

  inFileMC->Close();
  inFileDT->Close();
}
