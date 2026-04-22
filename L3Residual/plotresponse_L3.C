// Plot photon+jet distributions from analyse_PhotonJet.cc output
// Shows raw distributions (Data or MC) with optional normalized MC vs Data comparison
//
// Input: ROOT files produced by analyse_PhotonJet.cc
//
// Usage:
//   root -l -b -q 'plotresponse_L3.C("data.root", "Data")'
//   root -l -b -q 'plotresponse_L3.C("data.root", "Data", false, "mc.root")'  // with MC comparison

#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <vector>

#include "tdrStyle.C"
#include "CMS_lumi.C"

using namespace std;

static void NormalizeToUnityWidth(TH1* h) {
  if (!h) return;
  const double integralWidth = h->Integral("width");
  if (integralWidth <= 0) return; // assume non-zero for normal inputs
  h->Scale(1.0 / integralWidth, "width");
}

static void DrawSelectionText(double ptMin, double ptMax, double alphaCutMax,
                              double meanData, bool hasMC, double meanMC) {
  TLatex latex;
  latex.SetNDC(true);
  latex.SetTextFont(42);
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.18, 0.86, Form("%.0f < p_{T}^{#gamma} < %.0f GeV", ptMin, ptMax));
  latex.DrawLatex(0.18, 0.81, Form("#alpha < %.2f", alphaCutMax));
  latex.DrawLatex(0.18, 0.76, "B_{#gamma+jet} = p_{T}^{jet}/p_{T}^{#gamma}");
  if (hasMC) {
    latex.DrawLatex(0.18, 0.71, Form("#LTB_{#gamma+jet}^{Data}#GT = %.3f", meanData));
    latex.DrawLatex(0.18, 0.66, Form("#LTB_{#gamma+jet}^{MC}#GT = %.3f", meanMC));
  } else {
    latex.DrawLatex(0.18, 0.71, Form("#LTB_{#gamma+jet}#GT = %.3f", meanData));
  }
}

static void StyleAxes(TH1* hist, double yTitleOffset = 1.35) {
  if (!hist) return;
  hist->GetXaxis()->SetTitleSize(0.042);
  hist->GetYaxis()->SetTitleSize(0.042);
  hist->GetXaxis()->SetLabelSize(0.032);
  hist->GetYaxis()->SetLabelSize(0.032);
  hist->GetXaxis()->SetTitleOffset(1.05);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
}

void plotresponse_L3(TString inputFile = "",
                     TString tag = "Data",
                     bool isMC = false,
                     TString mcFile = "",
                     TString runLabel = "2024ppRef",
                     TString lumiLabel = "pp 480.4 pb^{-1}",
                     TString outputDir = "") {
  
  if (inputFile.IsNull()) {
    cout << "ERROR: No input file specified!" << endl;
    cout << "Usage: plotresponse_L3.C(\"data.root\", \"Data\", false)" << endl;
    cout << "       plotresponse_L3.C(\"data.root\", \"Data\", false, \"mc.root\")  // for comparison" << endl;
    return;
  }
  
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  writeExtraText = true;
  extraText = "Preliminary";
  cmsTextSize = 0.78;
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.Data());
  
  TFile* inFile = TFile::Open(inputFile, "READ");
  if (!inFile || inFile->IsZombie()) {
    cout << "ERROR: Cannot open input file: " << inputFile << endl;
    return;
  }
  
  TFile* mcFilePtr = nullptr;
  if (!mcFile.IsNull()) {
    mcFilePtr = TFile::Open(mcFile, "READ");
    if (!mcFilePtr || mcFilePtr->IsZombie()) {
      cout << "ERROR: Cannot open MC file: " << mcFile << endl;
      return;
    }
  }
  
  cout << "============================================" << endl;
  cout << "Photon+Jet Distribution Plotting" << endl;
  cout << "Input file: " << inputFile << endl;
  cout << "Tag: " << tag << endl;
  if (mcFilePtr) cout << "MC comparison file: " << mcFile << endl;
  cout << "============================================" << endl;
  
  // Create output directory
  string outfolder = outputDir.IsNull() ? Form("L3plots_%s", tag.Data()) : Form("%s/%s", outputDir.Data(), tag.Data());
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  
  // Access histograms - navigate to hibin and eta directories
  string histoPath = "hibin_-1.0_0.0/eta_-5.2_5.2";
  TDirectory* histDir = (TDirectory*)inFile->Get(histoPath.c_str());
  if (!histDir) {
    cout << "ERROR: Cannot find histogram directory: " << histoPath << endl;
    cout << "Available keys:" << endl;
    inFile->ls();
    return;
  }
  
  TDirectory* mcHistDir = nullptr;
  if (mcFilePtr) {
    mcHistDir = (TDirectory*)mcFilePtr->Get(histoPath.c_str());
    if (!mcHistDir) {
      cout << "WARNING: Cannot find MC histogram directory: " << histoPath << endl;
      mcFilePtr = nullptr;
    }
  }
  
  // List of histograms to plot with axis labels
  vector<pair<string, pair<string, string>>> histConfigs = {
    {"photon_pt", {"Photon p_{T} (GeV)", "Events"}},
    {"photon_eta", {"Photon #eta", "Events"}},
    {"photon_phi", {"Photon #phi (rad)", "Events"}},
    {"awayside_jet_pt", {"Jet p_{T} (GeV)", "Events"}},
    {"awayside_jet_eta", {"Jet #eta", "Events"}},
    {"awayside_jet_phi", {"Jet #phi (rad)", "Events"}},
    {"photonjet_dphi", {"#Delta#phi(#gamma,jet) (rad)", "Events"}},
    {"photonjet_alpha", {"#alpha", "Events"}},
    {"photonjet_ptavg", {"p_{T,avg} (GeV)", "Events"}}
  };
  
  TCanvas* c = new TCanvas("c", "Distribution", 800, 600);
  c->SetLogy();
  
  int nPlots = 0;
  
  for (const auto& config : histConfigs) {
    string histName = config.first;
    string xLabel = config.second.first;
    string yLabel = config.second.second;
    
    TH1D* h1 = (TH1D*)histDir->Get(histName.c_str());
    if (!h1) {
      cout << "WARNING: Cannot find histogram: " << histName << endl;
      continue;
    }
    
    TH1D* hMC = nullptr;
    if (mcHistDir) {
      hMC = (TH1D*)mcHistDir->Get(histName.c_str());
      if (!hMC) {
        cout << "WARNING: MC file missing histogram: " << histName << endl;
      }
    }
    
    // In comparison mode (Data+MC), only produce the overlay comparison plots.
    // In single-file mode, only produce the raw distributions.
    if (!mcFilePtr) {
      h1->SetLineColor(kBlue + 1);
      h1->SetLineWidth(2);
      h1->SetMarkerStyle(20);
      h1->SetMarkerColor(kBlue + 1);
      h1->SetMarkerSize(0.6);
      h1->GetXaxis()->SetTitle(xLabel.c_str());
      h1->GetYaxis()->SetTitle(yLabel.c_str());
      StyleAxes(h1);
      h1->Draw("E");

      CMS_lumi(c, 0, 0);
      c->Print(Form("%s/%s_%s.png", outfolder.c_str(), histName.c_str(), tag.Data()));

      nPlots++;
      continue;
    }
    
    // If MC file provided, plot normalized comparison
    if (hMC) {
      const int wasLogy = c->GetLogy();
      c->SetLogy(0);

      // Normalize to unit area (width-aware for variable binning)
      TH1D* h1_norm = (TH1D*)h1->Clone(Form("h1_norm_%s", histName.c_str()));
      TH1D* hMC_norm = (TH1D*)hMC->Clone(Form("hMC_norm_%s", histName.c_str()));

      NormalizeToUnityWidth(h1_norm);
      NormalizeToUnityWidth(hMC_norm);
      if (h1_norm->Integral() <= 0 || hMC_norm->Integral() <= 0) {
        delete h1_norm;
        delete hMC_norm;
        continue;
      }
      
      // Plot comparison
      double maxY = TMath::Max(h1_norm->GetMaximum(), hMC_norm->GetMaximum());
      h1_norm->SetMaximum(maxY * 1.3);
      
      h1_norm->SetLineColor(kRed + 1);
      h1_norm->SetLineWidth(2);
      h1_norm->SetMarkerColor(kRed + 1);
      h1_norm->GetXaxis()->SetTitle(xLabel.c_str());
      h1_norm->GetYaxis()->SetTitle("1/N dN/dX");
      StyleAxes(h1_norm, 1.40);
      h1_norm->Draw("E");
      
      hMC_norm->SetLineColor(kBlue + 1);
      hMC_norm->SetLineWidth(2);
      hMC_norm->SetMarkerColor(kBlue + 1);
      hMC_norm->Draw("E SAME");
      
      TLegend* leg = new TLegend(0.65, 0.7, 0.88, 0.85);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.035);
      leg->AddEntry(h1_norm, "Data", "lp");
      leg->AddEntry(hMC_norm, "MC", "lp");
      leg->Draw();
      
      CMS_lumi(c, 0, 0);
      c->Print(Form("%s/%s_comparison.png", outfolder.c_str(), histName.c_str()));

      c->SetLogy(wasLogy);
      
      delete leg;
      delete h1_norm;
      delete hMC_norm;

      nPlots++;
    }
  }
  
  // Plot balance distributions from TH3D
  TH3D* balance_dist = (TH3D*)histDir->Get("photonjet_balance_dist");
  TH3D* balance_dist_mc = nullptr;
  if (mcHistDir) {
    balance_dist_mc = (TH3D*)mcHistDir->Get("photonjet_balance_dist");
  }
  if (balance_dist) {
    cout << "\n=== Plotting balance distributions ===" << endl;
    cout << "Balance distribution entries: " << balance_dist->GetEntries() << endl;
    
    // Get binning from histogram axes (no fixed binning)
    int nPtBins = balance_dist->GetXaxis()->GetNbins();
    int nAlphaBins = balance_dist->GetYaxis()->GetNbins();
    
    cout << "pT bins: " << nPtBins << ", Alpha bins: " << nAlphaBins << endl;
    
    c->SetLogy(0);
    
    int nBalancePlots = 0;
    
    // Loop over all pT bins
    for (int iPt = 1; iPt <= nPtBins; iPt++) {
      double ptMin = balance_dist->GetXaxis()->GetBinLowEdge(iPt);
      double ptMax = balance_dist->GetXaxis()->GetBinUpEdge(iPt);
      
      // Use cumulative alpha cuts (alpha < alphaCutMax), not alpha-bin slices
      for (int iAlpha = 1; iAlpha <= nAlphaBins; iAlpha++) {
        const double alphaCutMax = balance_dist->GetYaxis()->GetBinUpEdge(iAlpha);
        const double alphaEps = 1e-6;
        const int alphaMaxBin = balance_dist->GetYaxis()->FindBin(alphaCutMax - alphaEps);

        // Project to balance axis for this (pT bin, alpha < cut)
        balance_dist->GetXaxis()->SetRange(iPt, iPt);
        balance_dist->GetYaxis()->SetRange(1, alphaMaxBin);
        TH1D* h_bal = (TH1D*)balance_dist->Project3D("z");
        h_bal->SetDirectory(nullptr);
        h_bal->SetName(Form("balance_pt%.0f_%.0f_alphalt%.3f", ptMin, ptMax, alphaCutMax));
        balance_dist->GetXaxis()->SetRange(0, 0);
        balance_dist->GetYaxis()->SetRange(0, 0);
        
        if (h_bal->GetEntries() < 1) {
          delete h_bal;
          continue;
        }
        
        h_bal->SetLineColor(kBlue + 1);
        h_bal->SetLineWidth(2);
        h_bal->SetMarkerStyle(20);
        h_bal->SetMarkerColor(kBlue + 1);
        h_bal->SetMarkerSize(0.6);
        h_bal->GetXaxis()->SetTitle("B_{#gamma+jet}");
        h_bal->GetYaxis()->SetTitle("Events");
        StyleAxes(h_bal, 1.40);
        h_bal->SetTitle("");
        const double meanData = h_bal->GetMean();

        // In single-file mode, write the raw balance slices
        if (!mcFilePtr) {
          h_bal->Draw("E");
          DrawSelectionText(ptMin, ptMax, alphaCutMax, meanData, false, 0.0);

          TLegend* legBal = new TLegend(0.65, 0.78, 0.88, 0.88);
          legBal->SetBorderSize(0);
          legBal->SetFillStyle(0);
          legBal->SetTextSize(0.035);
          legBal->AddEntry(h_bal, "Data", "lp");
          legBal->Draw();

          CMS_lumi(c, 0, 0);
          c->Print(Form("%s/balance_pt%.0f-%.0f_alphalt%.3f_%s.png",
                        outfolder.c_str(), ptMin, ptMax, alphaCutMax, tag.Data()));

          delete legBal;
          delete h_bal;
          nBalancePlots++;
          continue;
        }

        // In comparison mode, only make normalized Data vs MC overlay plots
        if (balance_dist_mc) {
          // Use the same bin indices for MC (same binning assumed)
          balance_dist_mc->GetXaxis()->SetRange(iPt, iPt);
          balance_dist_mc->GetYaxis()->SetRange(1, iAlpha);
          TH1D* h_bal_mc = (TH1D*)balance_dist_mc->Project3D("z");
          h_bal_mc->SetDirectory(nullptr);
          h_bal_mc->SetName(Form("balance_mc_pt%.0f_%.0f_alphalt%.3f", ptMin, ptMax, alphaCutMax));
          balance_dist_mc->GetXaxis()->SetRange(0, 0);
          balance_dist_mc->GetYaxis()->SetRange(0, 0);

          const double meanMC = h_bal_mc->GetMean();

          TH1D* h_bal_norm = (TH1D*)h_bal->Clone(Form("balance_norm_pt%.0f_%.0f_alphalt%.3f", ptMin, ptMax, alphaCutMax));
          TH1D* h_bal_mc_norm = (TH1D*)h_bal_mc->Clone(Form("balance_mc_norm_pt%.0f_%.0f_alphalt%.3f", ptMin, ptMax, alphaCutMax));

          NormalizeToUnityWidth(h_bal_norm);
          NormalizeToUnityWidth(h_bal_mc_norm);
          const double maxYBal = TMath::Max(h_bal_norm->GetMaximum(), h_bal_mc_norm->GetMaximum());
          h_bal_norm->SetMaximum(maxYBal * 1.3);
          h_bal_norm->SetLineColor(kRed + 1);
          h_bal_norm->SetMarkerColor(kRed + 1);
          h_bal_norm->SetLineWidth(2);
          h_bal_norm->GetXaxis()->SetTitle("B_{#gamma+jet}");
          h_bal_norm->GetYaxis()->SetTitle("1/N dN/dB_{#gamma+jet}");
          StyleAxes(h_bal_norm, 1.40);
          h_bal_norm->Draw("E");

          h_bal_mc_norm->SetLineColor(kBlue + 1);
          h_bal_mc_norm->SetMarkerColor(kBlue + 1);
          h_bal_mc_norm->SetLineWidth(2);
          h_bal_mc_norm->Draw("E SAME");

          DrawSelectionText(ptMin, ptMax, alphaCutMax, meanData, true, meanMC);

          TLegend* legBalComp = new TLegend(0.65, 0.74, 0.88, 0.88);
          legBalComp->SetBorderSize(0);
          legBalComp->SetFillStyle(0);
          legBalComp->SetTextSize(0.035);
          legBalComp->AddEntry(h_bal_norm, "Data", "lp");
          legBalComp->AddEntry(h_bal_mc_norm, "MC", "lp");
          legBalComp->Draw();

          CMS_lumi(c, 0, 0);
          c->Print(Form("%s/balance_pt%.0f-%.0f_alphalt%.3f_%s_DataVsMC.png",
                        outfolder.c_str(), ptMin, ptMax, alphaCutMax, tag.Data()));

          delete legBalComp;        

          delete h_bal_norm;
          delete h_bal_mc_norm;
          delete h_bal_mc;
        }

        delete h_bal;
        nBalancePlots++;
      }
    }
    
    cout << "Created " << nBalancePlots << " balance distribution plots" << endl;
    nPlots += nBalancePlots;
    
    c->SetLogy();
  } else {
    cout << "WARNING: photonjet_balance_dist not found" << endl;
  }
  
  delete c;
  
  inFile->Close();
  if (mcFilePtr) mcFilePtr->Close();
  
  cout << "\n============================================" << endl;
  cout << "Plots saved in: " << outfolder << "/" << endl;
  cout << "Total plots: " << nPlots << endl;
  cout << "============================================" << endl;
}
