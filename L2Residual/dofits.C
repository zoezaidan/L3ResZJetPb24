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

void dofits(TString inzb, TString inHP, float fitmin = 0.15,
            float fitmax = 0.35,
            string outfilename = "kfactor_rerunall_combined_allpts",
            string outfolder = "L2fits", bool doabseta = true) {
  // input file is from 3D derivation
  gStyle->SetOptStat(0);

  TFile *inFilezb = new TFile(inzb, "READ");
  TFile *inFile = new TFile(inHP, "READ");

  gSystem->mkdir(outfolder.c_str());
  TFile *outfile = new TFile(
      Form("%s/%s.root", outfolder.c_str(), outfilename.c_str()), "RECREATE");

  // ROOT::Fit::DataRange range(0.3,0.4);
  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange range;
  range.SetRange(fitmin, fitmax);

  //  mg->Add(g2,"L");

  const int nPtBins = histograms::nptforjec;
  const int etaMax = histograms::nwabsetas;

  std::vector<int> availablePtBins;
  for (int ptbin = 1; ptbin <= nPtBins; ++ptbin) {
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

  // TODO: change this so that histogram contents are copied into TGraphs and
  // the reference alpha is excluded
  // TODO: remove reference alpha from the fit
  int colours[] = {209, 226, 213, 51, 206, 209};
  for (int etabin = 1; etabin <= etaMax; ++etabin) {

    ROOT::Fit::BinData data(opt, range);
    map<int, TH1D *> histos;
    map<int, TGraphErrors *> graphs;
    auto multifit = new TMultiGraph();

    // for (int etabin = 18; etabin < 19; ++etabin) { // TODO bins
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    auto leg = new TLegend(0.57, 0.7, 0.85, 0.85); //  x, y, x, y
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // int etabin = 18;
    //  pick histos for different pT:s
    /*
      for ( int ptbin = 1; ptbin <=2; ++ptbin)  histos[ptbin] =
     (TH1D*)inFilezb->Get(Form("Respvsa_norm_%d_%d",ptbin,etabin)); for ( int
     ptbin = 3; ptbin <=4; ++ptbin)  histos[ptbin] =
     (TH1D*)inFile->Get(Form("Respvsa_norm_%d_%d",ptbin,etabin));


     for ( int ptbin = 1; ptbin <=2; ++ptbin)  ratios[ptbin] =
     (TH1D*)inFilezb->Get(Form("Respvsa_%d_%d",ptbin,etabin)); for ( int ptbin =
     3; ptbin <=4; ++ptbin)  ratios[ptbin] =
     (TH1D*)inFile->Get(Form("Respvsa_%d_%d",ptbin,etabin));
    */

    for (int ptbin : availablePtBins) {
      const double highEdge = histograms::ptforjec[ptbin];
      const bool preferZB =
          (highEdge <= 80.0); // Preserve legacy low-pT from ZB preference
      TH1D *h = nullptr;

      if (preferZB)
        h = getH1(inFilezb, Form("Respvsa_norm_%d_%d", ptbin, etabin));
      if (!h)
        h = getH1(inFile, Form("Respvsa_norm_%d_%d", ptbin, etabin));
      if (!h && !preferZB)
        h = getH1(inFilezb, Form("Respvsa_norm_%d_%d", ptbin, etabin));

      if (h)
        histos[ptbin] = h;
    }

    if (histos.empty()) {
      delete c1;
      delete multifit;
      continue;
    }

    auto first = histos.begin()->second;
    first->SetMaximum(1.15);
    first->SetMinimum(0.85);
    first->Draw("AXIS");

    for (int ptbin : availablePtBins) {
      if (!histos.count(ptbin))
        continue;
      if (ptbin == availablePtBins.front())
        continue;
      //       if (etabin >= 14 and ptbin == 6) break;
      //       if (etabin >= 17 and ptbin == 5) break;
      //  if (etabin > 8 and ptbin == 6) break;
      // if (etabin > 11 and ptbin == 5) break;
      //     if (etabin > 11 and ptbin == 4) break;

      Double_t results[histograms::nalphavaluesgraph];
      Double_t errors[histograms::nalphavaluesgraph];
      Double_t errorsx[histograms::nalphavaluesgraph];

      for (int abin = 0; abin < histograms::nalphavaluesgraph; ++abin) {
        results[abin] = histos[ptbin]->GetBinContent(abin + 1);
        errors[abin] = histos[ptbin]->GetBinError(abin + 1);
        errorsx[abin] = 0;
      }

      graphs[ptbin] = new TGraphErrors(
          histograms::nalphavaluesgraph, histograms::alphavaluesgraph, results,
          errorsx, errors); // TODO: do not include reference?
      //      cout << graphs[ptbin] << endl

      // histograms::alphavalues[]
      ROOT::Fit::FillData(data, histos[ptbin]);
      // ROOT::Fit::FillData(data, graphs[ptbin]);   // TODO: Replace with
      // TMultiGraph fitter
      multifit->Add(graphs[ptbin], "P");

      histos[ptbin]->SetMaximum(1.15);
      histos[ptbin]->SetMinimum(0.85);

      histos[ptbin]->GetXaxis()->SetTitle("#alpha");
      histos[ptbin]->GetYaxis()->SetTitle(
          " < MC / Data > / < MC / Data >_{#alpha < 0.3}"); // TODO: correc
                                                            // reference label

      histos[ptbin]->SetLineColor(colours[ptbin - 1]);
      // histos[ptbin]->Draw("same E1");

      graphs[ptbin]->SetMarkerStyle(kFullCircle);
      graphs[ptbin]->SetMarkerColor(colours[ptbin - 1]);
      graphs[ptbin]->Draw("sameP");

      leg->AddEntry(histos[ptbin], Form("%s", ptBinLabel(ptbin).c_str()));
      // TODO: colours
    }

    leg->Draw("same");

    /* TF1 * f1 = new TF1("f1","pol1",fitmin,fitmax);
    f1->SetParameters(1,0);

     ROOT::Math::WrappedTF1 wf(*f1);

    ROOT::Fit::Fitter fitter;
    fitter.SetFunction(wf);

    fitter.Fit(data);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

    cout << result.Chi2() <<  " " << result.Ndf() << " " << result.Parameter(0)
    << " " << result.ParError(0) << endl; */

    // f1->Draw("same");

    TF1 *f1 = new TF1("f1", "pol1", fitmin, fitmax);
    f1->SetParameters(1, 0);
    multifit->Fit("f1", "R");
    f1->Draw("same");

    factors->SetBinContent(etabin, f1->GetParameter(0));
    factors->SetBinError(etabin, f1->GetParError(0));

    // TODO: debug this
    TLatex *txt = new TLatex();
    txt->SetTextSize(0.03);
    txt->SetNDC();

    if (!doabseta)
      txt->DrawLatex(0.2, 0.8,
                     Form("%.2f < #eta < %.2f",
                          histograms::wetarange[etabin - 1],
                          histograms::wetarange[etabin])); // TODO: correct bin
    else
      txt->DrawLatex(
          0.2, 0.8,
          Form("%.2f < |#eta| < %.2f", histograms::wabsetarange[etabin - 1],
               histograms::wabsetarange[etabin])); // TODO: correct bin

    // txt->DrawLatex( 0.2, 0.85, Form("p0 = %.5f #pm %.5f",result.Parameter(0),
    // result.ParError(0)));
    txt->DrawLatex(
        0.2, 0.85,
        Form("p0 = %.5f #pm %.5f", f1->GetParameter(0), f1->GetParError(0)));

    c1->Print(Form("%s/fits_eta_%d.png", outfolder.c_str(), etabin));

    //    TCanvas *c3 = new TCanvas("c3","c3",800,600);
    //  graphs[3]->Draw("");
  }

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  factors->SetMaximum(1.15);
  factors->SetMinimum(0.95);

  // TF1 * f2 = new TF1("f2","pol1",0.0,3.0);
  TF1 *f2 = new TF1("f2", "[0]+[1]*cosh(x)/(1+[2]*cosh(x))", 0.0, 3.0);
  //  factors->Fit(f2);

  if (doabseta)
    factors->GetXaxis()->SetTitle(" |#eta| ");
  else
    factors->GetXaxis()->SetTitle(" #eta ");

  factors->Draw();
  c2->Print(Form("%s/kfactors.png", outfolder.c_str()));

  factors->Write();
  factors->Write("kfactors");

  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);

  auto leg2 = new TLegend(0.57, 0.7, 0.85, 0.85); //  x, y, x, y
  leg2->SetTextSize(0.03);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);

  map<string, TH1D *> histos;

  // Save non-parametrized correction; response ratio multiplied by kFSR
  // 2023 configuration

  /*  histos["30to80"] = (TH1D*)inFilezb->Get("ratio_pt30to80_alpha0.3");
  histos["80to120"] = (TH1D*)inFile->Get("ratio_pt80to120_alpha0.3");
  histos["120to1000"] = (TH1D*)inFile->Get("ratio_pt120to1000_alpha0.3"); */

  std::vector<std::string> ptbins;
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

  int col = 0;
  for (auto bin : ptbins) {

    for (int hb = 1; hb < histos[bin.c_str()]->GetXaxis()->GetNbins(); ++hb) {
      histos[bin.c_str()]->SetBinError(hb, 0.);
      if (bin == "120to1000" and hb > 11)
        histos[bin.c_str()]->SetBinContent(hb, 0.);
    }

    histos[bin.c_str()]->Multiply(factors);
    histos[bin.c_str()]->SetLineColor(colours[col]);
    histos[bin.c_str()]->SetMarkerStyle(kFullCircle);
    histos[bin.c_str()]->SetMarkerColor(colours[col]);

    if (doabseta)
      histos[bin.c_str()]->GetXaxis()->SetTitle(" |#eta| ");
    else
      histos[bin.c_str()]->GetXaxis()->SetTitle(" #eta ");
    histos[bin.c_str()]->GetYaxis()->SetTitle(" L2 residual correction ");

    histos[bin.c_str()]->SetMaximum(1.3);
    histos[bin.c_str()]->SetMinimum(0.9);

    histos[bin.c_str()]->Draw("same");
    histos[bin.c_str()]->Write(Form("corrections_%s", bin.c_str()));
    leg2->AddEntry(histos[bin.c_str()], Form("%s", bin.c_str()));

    col++;
  }
  leg2->Draw();

  c3->Print(Form("%s/corrections-closure.png", outfolder.c_str()));
  c3->Print(Form("%s/corrections-closure.pdf", outfolder.c_str()));

  outfile->Close();
}
