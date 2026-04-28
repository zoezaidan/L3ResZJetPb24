// Do fits for pt-parametization of L2residuals

#include "../fillhistograms/histograms.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace {
std::string edgeLabel(double x) {
  double xi = std::round(x);
  if (std::fabs(x - xi) < 1e-6)
    return std::to_string(static_cast<int>(xi));
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

std::string ptBinLabel(int ptbin) {
  return edgeLabel(histograms::ptforjec[ptbin - 1]) + "to" +
         edgeLabel(histograms::ptforjec[ptbin]);
}

TH1D *getH1(TFile *f, const std::string &name) {
  if (!f)
    return nullptr;
  return dynamic_cast<TH1D *>(f->Get(name.c_str()));
}
} // namespace

void fit_pt_param(TString inzb = "L2residuals_pbpbreco_rereco_zb_jetid.root",
                  TString inHP = "L2residuals_pbpbreco_rereco_hp_jetid.root",
                  float fitmin = 25., float fitmax = 1000,
                  string outfilename = "testing_pt_dep",
                  string outfolder = "ptfits", bool doabseta = true) {

  // input file is from 3D derivation
  gStyle->SetOptStat(0);

  TFile *inFilezb = new TFile(inzb, "READ");
  TFile *inFile = new TFile(inHP, "READ");

  gSystem->mkdir(outfolder.c_str());
  TFile *outfile = new TFile(
      Form("%s/%s.root", outfolder.c_str(), outfilename.c_str()), "RECREATE");

  std::vector<int> availablePtBins;
  for (int ptbin = 1; ptbin <= histograms::nptforjec; ++ptbin) {
    const std::string lbl = ptBinLabel(ptbin);
    const std::string ratioName = "ratio_pt" + lbl + "_alpha0.3";
    if (getH1(inFilezb, ratioName) || getH1(inFile, ratioName))
      availablePtBins.push_back(ptbin);
  }

  if (availablePtBins.empty()) {
    std::cerr << "No ratio_pt* histograms found in inputs." << std::endl;
    return;
  }

  TH1D *factors = nullptr;
  for (int ptbin : availablePtBins) {
    const std::string lbl = ptBinLabel(ptbin);
    const std::string ratioName = "ratio_pt" + lbl + "_alpha0.3";
    factors = getH1(inFile, ratioName);
    if (!factors)
      factors = getH1(inFilezb, ratioName);
    if (factors)
      break;
  }

  if (!factors) {
    std::cerr << "Could not initialize factors histogram from ratio_pt* inputs."
              << std::endl;
    return;
  }
  factors->Reset();

  int colours[] = {209, 226, 213, 51, 206, 209};

  TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
  c3->SetLogx();

  // Fro mresponses
  map<string, TH1D *> histos;
  std::vector<std::string> ptbins;
  map<int, TH1D *> histosvspt;

  // Pick up histograms for different pT bins from ZB/HP with legacy low/high
  // preference
  for (int ptbin : availablePtBins) {
    const std::string lbl = ptBinLabel(ptbin);
    const std::string ratioName = "ratio_pt" + lbl + "_alpha0.3";
    const double highEdge = histograms::ptforjec[ptbin];
    const bool preferZB = (highEdge <= 80.0);

    TH1D *h = nullptr;
    if (preferZB)
      h = getH1(inFilezb, ratioName);
    if (!h)
      h = getH1(inFile, ratioName);
    if (!h && !preferZB)
      h = getH1(inFilezb, ratioName);
    if (!h)
      continue;

    histos[lbl] = h;
    ptbins.push_back(lbl);
  }

  if (ptbins.empty()) {
    std::cerr << "No usable ratio_pt* histograms available for pt-param fits."
              << std::endl;
    return;
  }

  std::vector<double> ptEdges;
  ptEdges.reserve(ptbins.size() + 1);
  bool first = true;
  for (int ptbin : availablePtBins) {
    const std::string lbl = ptBinLabel(ptbin);
    if (!histos.count(lbl))
      continue;
    if (first) {
      ptEdges.push_back(histograms::ptforjec[ptbin - 1]);
      first = false;
    }
    ptEdges.push_back(histograms::ptforjec[ptbin]);
  }

  auto txt = new TLatex();
  txt->SetTextSize(0.03);
  txt->SetNDC();

  TF1 *f1 = new TF1("f1", "[0] + [1]*log(x)", fitmin, fitmax);
  TF1 *f2 = new TF1("f2", "pol0", fitmin, fitmax);
  TF1 *f3 = new TF1("f3", "1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))", fitmin,
                    fitmax); // Run3 parametrization

  // Histograms vs. pT
  for (int ebin = 1; ebin < 12; ++ebin) { // Keep barrel-focused default range
    auto leg2 = new TLegend(0.12, 0.10, 0.35, 0.36); //  x, y, x, y
    leg2->SetTextSize(0.03);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    histosvspt[ebin] =
        new TH1D(Form("ebin_%d", ebin), "", ptEdges.size() - 1, &ptEdges[0]);

    // etabins for
    for (int ptbin = 1; ptbin <= static_cast<int>(ptbins.size()); ++ptbin) {
      histosvspt[ebin]->SetBinContent(
          ptbin, histos[ptbins[ptbin - 1]]->GetBinContent(ebin));
      histosvspt[ebin]->SetBinError(
          ptbin, histos[ptbins[ptbin - 1]]->GetBinError(ebin));
    }

    histosvspt[ebin]->SetMaximum(1.1);
    histosvspt[ebin]->SetMinimum(0.9);
    histosvspt[ebin]->GetXaxis()->SetTitle("p_{T,avg} (GeV)");
    histosvspt[ebin]->Draw();

    f1->SetParameters(1, 0);
    histosvspt[ebin]->Fit("f1", "R");

    f1->Draw("same");

    f2->SetParameters(1);
    f2->SetLineStyle(kDashed);
    histosvspt[ebin]->Fit("f2", "R");
    f2->Draw("same");

    f3->SetLineStyle(kDashed);
    f3->SetLineColor(kBlue + 1);
    histosvspt[ebin]->Fit("f3", "R");
    f3->Draw("same");

    leg2->AddEntry(
        f1, Form("a + b*log(p_{T}), a = %.3f #pm %.3f, b = %.3f #pm %.3f",
                 f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(1),
                 f1->GetParError(1)));
    leg2->AddEntry(
        f1, Form("#chi2/ndof = %.3f / %d", f1->GetChisquare(), f1->GetNDF()));
    leg2->AddEntry(f2, Form("p0 = %.3f #pm %.3f", f2->GetParameter(0),
                            f2->GetParError(0)));
    leg2->AddEntry(
        f2, Form("#chi2/ndof = %.3f / %d", f2->GetChisquare(), f2->GetNDF()));
    leg2->AddEntry(f3, "1./(a+b*log10(0.01*p_{T})+c/(p_{T}/10.))");
    leg2->AddEntry(
        f3, Form(" a = %.3f #pm %.3f, b = %.3f #pm %.3f, c = %.3f #pm %.3f",
                 f3->GetParameter(0), f3->GetParError(0), f3->GetParameter(1),
                 f3->GetParError(1), f3->GetParameter(2), f3->GetParError(2)));
    leg2->AddEntry(
        f3, Form("#chi2/ndof = %.3f / %d", f3->GetChisquare(), f3->GetNDF()));

    leg2->Draw();

    txt->DrawLatex(0.5, 0.8,
                   Form("%.3f < |#eta| < %.3f ",
                        histos[ptbins[0]]->GetBinLowEdge(ebin),
                        histos[ptbins[0]]->GetBinLowEdge(ebin + 1)));

    f1->Write(Form("loglin_eta%d", ebin));
    f2->Write(Form("const_eta%d", ebin));
    f3->Write(Form("run3_eta%d", ebin));
    c3->Print(Form("%s/ebin_%d.png", outfolder.c_str(), ebin));
  }

  outfile->Close();
}
