// Purpose: fit the active L3 residual response workflow in pTref.
//
// Step 1. Read one or more derived ROOT files from deriveL3_from_photonjet.C.
// Step 2. Fit the alpha dependence in every pTref bin to extract kFSR and the
//         alpha->0 corrected response for each input independently.
// Step 3. Combine the selected pTref points from all inputs and run the
//         final full direct pTref fit used by the export step.
// Step 4. Save the fit ROOT file, metadata, and validation plots needed by
//         createL2L3ResTextFile.C.

#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TParameter.h"
#include "TProfile3D.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "CMS_lumi.C"
#include "L3ResCommon.h"
#include "tdrStyle.C"

using namespace std;

static void useCmsPlotStyle(const string &lumiLabel) {
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  writeExtraText = true;
  extraText = "Preliminary";
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.c_str());
}

static void drawCmsLabel(TCanvas &canvas) {
  canvas.cd();
  CMS_lumi(&canvas, 0, 0);
  canvas.Modified();
  canvas.Update();
}

static void styleFrameAxes(TH1D *frame, double yTitleOffset = 1.35) {
  if (!frame)
    return;
  frame->GetXaxis()->SetTitleSize(0.046);
  frame->GetYaxis()->SetTitleSize(0.046);
  frame->GetXaxis()->SetLabelSize(0.034);
  frame->GetYaxis()->SetLabelSize(0.034);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(yTitleOffset);
}

static TString balanceAxisLabel(const PreparedInput &prepared) {
  return Form("B_{%s}", prepared.label.Data());
}

static TString responseAxisLabel(const PreparedInput &prepared) {
  return Form("R_{%s}", prepared.label.Data());
}

static TString alphaRatioAxisLabel(const PreparedInput &prepared) {
  return Form("r_{%s}(#alpha)", prepared.label.Data());
}

struct SequentialFitResult {
  TF1 *f0 = nullptr;
  TF1 *f1 = nullptr;
  TF1 *f2 = nullptr;
  TF1 *f3 = nullptr;
  TF1 *f4 = nullptr;
  TF1 *finalFit = nullptr;
  int finalStatus = -1;
  bool valid = false;
};

static double fitChi2Ndf(const TF1 *fit) {
  if (!fit)
    return -1.0;
  const int ndf = fit->GetNDF();
  if (ndf <= 0)
    return -1.0;
  return fit->GetChisquare() / ndf;
}

static bool fitNeedsReview(const TF1 *fit, int fitStatus,
                           double chi2Limit = 5.0) {
  if (fitStatus != 0)
    return true;
  const double chi2PerNdf = fitChi2Ndf(fit);
  return (chi2PerNdf > chi2Limit);
}

static const char *finalPtRefFitExpr() {
  return "1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/"
         "[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x)";
}

static void initFinalPtRefFit(TF1 *fit, double ySeed = 1.0) {
  if (!fit)
    return;

  const double safeSeed = (ySeed > 0.0 ? ySeed : 1.0);
  fit->SetParameters(0.9722 / safeSeed, 0.7944, 2.14069, 0.10229, 10.52, 1.5550,
                     -0.72222, -2.25e-06);
  fit->SetParLimits(0, 0.5, 2.0);
  fit->SetParLimits(4, 1.0, 100.0);
  fit->SetParLimits(5, 0.1, 5.0);
}

// Purpose: keep the historical sequential pTref fits for diagnostics, but use
// the full direct pTref model as the final shared fit.
static SequentialFitResult runSequentialTemplateFit(TGraphErrors *graph,
                                                    double fitMin,
                                                    double fitMax,
                                                    const TString &fitTag) {
  SequentialFitResult result;
  if (!graph || graph->GetN() < 3 || fitMax <= fitMin)
    return result;

  const double ySeed = graphMeanY(graph);
  const double oneOverXLimit = std::max(0.5, 0.40 * fitMin / 10.0);

  result.f0 =
      new TF1(Form("f_ptref_f0_%s", fitTag.Data()), "[0]", fitMin, fitMax);
  result.f1 = new TF1(Form("f_ptref_f1_%s", fitTag.Data()),
                      "[0]+[1]*log10(0.01*x)", fitMin, fitMax);
  result.f2 = new TF1(Form("f_ptref_f2_%s", fitTag.Data()),
                      "[0]+[1]*log10(0.01*x)+[2]/(x/10.)", fitMin, fitMax);
  result.f3 =
      new TF1(Form("f_ptref_f3_%s", fitTag.Data()),
              "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)", fitMin, fitMax);
  result.f4 =
      new TF1(Form("f_ptref_f4_%s", fitTag.Data()),
              "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)+[3]/(x/10.)",
              fitMin, fitMax);
  result.finalFit = new TF1(Form("f_ptref_final_%s", fitTag.Data()),
                            finalPtRefFitExpr(), fitMin, fitMax);

  result.f0->SetParameter(0, ySeed);
  result.f0->SetParLimits(0, 0.5, 1.5);
  graph->Fit(result.f0, "QRN");

  result.f1->SetParameters(result.f0->GetParameter(0), -0.01);
  result.f1->SetParLimits(0, 0.5, 1.5);
  result.f1->SetParLimits(1, -0.20, 0.20);
  graph->Fit(result.f1, "QRN");

  result.f2->SetParameters(result.f1->GetParameter(0),
                           result.f1->GetParameter(1), 0.02);
  result.f2->SetParLimits(0, 0.5, 1.5);
  result.f2->SetParLimits(1, -0.20, 0.20);
  result.f2->SetParLimits(2, -oneOverXLimit, oneOverXLimit);
  graph->Fit(result.f2, "QRN");

  result.f3->SetParameters(result.f1->GetParameter(0),
                           result.f1->GetParameter(1), 0.005);
  result.f3->SetParLimits(0, 0.5, 1.5);
  result.f3->SetParLimits(1, -0.20, 0.20);
  result.f3->SetParLimits(2, -0.15, 0.15);
  graph->Fit(result.f3, "QRN");

  result.f4->SetParameters(
      result.f3->GetParameter(0), result.f3->GetParameter(1),
      result.f3->GetParameter(2), result.f2->GetParameter(2));
  result.f4->SetParLimits(0, 0.5, 1.5);
  result.f4->SetParLimits(1, -0.20, 0.20);
  result.f4->SetParLimits(2, -0.15, 0.15);
  result.f4->SetParLimits(3, -oneOverXLimit, oneOverXLimit);
  graph->Fit(result.f4, "QRN");

  initFinalPtRefFit(result.finalFit, ySeed);
  result.finalStatus = (int)graph->Fit(result.finalFit, "QRN");
  result.finalStatus = (int)graph->Fit(result.finalFit, "QRN");

  result.valid = std::isfinite(result.finalFit->GetParameter(0)) &&
                 std::isfinite(result.finalFit->GetParameter(7));
  return result;
}

void L3Res(TString inFileL3Derived = "L3Residual/L3Residual_ZJet_jtpt30_Z40.root",
           TString sampleTypesCSV = "Zjet", TString inputPtRangesCSV = "",
           string outfilename = "L3Res_Z40jet30",
           string runLabel = "2024ppRef", string lumiLabel = "pp 480.4 pb^{-1}",
           int refAlphaBin = 5, double fitAlphaMin = 0.0,
           double fitAlphaMax = 0.4, string outBaseDir = "L3Residual") {

  useCmsPlotStyle(lumiLabel);

  const string outfolder = outBaseDir + "/" + outfilename;
  const string pdfFolder = outfolder + "/pdf";
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  gSystem->mkdir(pdfFolder.c_str(), kTRUE);

  vector<TString> inputFiles = splitCsvTokens(inFileL3Derived);
  if (inputFiles.empty())
    inputFiles.push_back(inFileL3Derived);

  vector<SampleType> sampleTypes;
  string sampleError;
  if (!parseSampleTypesCSV(sampleTypesCSV, inputFiles.size(), sampleTypes,
                           sampleError)) {
    cout << "ERROR: " << sampleError << endl;
    return;
  }

  const vector<TString> fitWindowTokens = splitCsvTokens(inputPtRangesCSV);
  const int colors[] = {kBlue + 1,    kRed + 1,    kGreen + 2,
                        kMagenta + 2, kOrange + 7, kCyan + 2};
  const int markers[] = {kFullCircle,  kFullSquare,       kFullTriangleUp,
                         kFullDiamond, kFullTriangleDown, kFullCross};
  const int colorCount = sizeof(colors) / sizeof(colors[0]);
  const int markerCount = sizeof(markers) / sizeof(markers[0]);

  TFile *outputFile =
      new TFile(Form("%s/%s_fit.root", outfolder.c_str(), outfilename.c_str()),
                "RECREATE");

  // Purpose: clone a required histogram and stop the macro with a clear error
  // if the contract is broken.
  auto cloneRequiredHist = [](TFile *inputFile, const TString &name,
                              const TString &context) -> TH1D * {
    TH1D *hist = inputFile ? (TH1D *)inputFile->Get(name) : nullptr;
    if (!hist) {
      cout << "ERROR: Missing required histogram '" << name << "' in "
           << context << endl;
      return nullptr;
    }
    TH1D *clone = (TH1D *)hist->Clone(Form("%s_clone", name.Data()));
    clone->SetDirectory(nullptr);
    return clone;
  };

  // Purpose: load the full alpha series for a required histogram family.
  auto loadSeries = [&](TFile *inputFile, const TString &prefix,
                        const TString &context,
                        bool required = true) -> map<int, TH1D *> {
    map<int, TH1D *> series;
    for (int alphaBin = 1; alphaBin <= 50; ++alphaBin) {
      TH1D *sourceHist =
          inputFile
              ? (TH1D *)inputFile->Get(Form("%s%d", prefix.Data(), alphaBin))
              : nullptr;
      TH1D *hist = nullptr;
      if (sourceHist) {
        hist = (TH1D *)sourceHist->Clone(
            Form("%s%d_clone", prefix.Data(), alphaBin));
        hist->SetDirectory(nullptr);
      }
      if (!hist) {
        if (alphaBin == 1 && required) {
          cout << "ERROR: Missing required histogram '"
               << Form("%s%d", prefix.Data(), alphaBin) << "' in " << context
               << endl;
          series.clear();
        }
        break;
      }
      series[alphaBin] = hist;
    }
    return series;
  };

  auto normalizeSeriesByReference =
      [&](const map<int, TH1D *> &ratioSeries, int referenceAlphaBin,
          const TString &prefix, const TString &context) -> map<int, TH1D *> {
    map<int, TH1D *> normalizedSeries;
    auto refIt = ratioSeries.find(referenceAlphaBin);
    if (refIt == ratioSeries.end() || !refIt->second) {
      cout << "ERROR: Missing required reference alpha histogram for "
              "normalized series '"
           << prefix << "' in " << context << endl;
      return normalizedSeries;
    }

    TH1D *referenceHist = refIt->second;
    for (const auto &entry : ratioSeries) {
      TH1D *hist = (TH1D *)entry.second->Clone(
          Form("%s%d_clone", prefix.Data(), entry.first));
      hist->SetDirectory(nullptr);

      for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        const double value = hist->GetBinContent(bin);
        const double error = hist->GetBinError(bin);
        const double referenceValue = referenceHist->GetBinContent(bin);
        const double referenceError = referenceHist->GetBinError(bin);

        if (!std::isfinite(value) || !std::isfinite(referenceValue) ||
            !(referenceValue > 0.0) || !(value > 0.0)) {
          hist->SetBinContent(bin, 0.0);
          hist->SetBinError(bin, 0.0);
          continue;
        }

        if (entry.first == referenceAlphaBin) {
          hist->SetBinContent(bin, 1.0);
          hist->SetBinError(bin, 0.0);
          continue;
        }

        // ROOT tracks the per-histogram bin errors, but not the covariance
        // between cumulative-alpha ratios from different histograms. When an
        // older derived file does not already store ratio_norm_*, use a
        // conservative independent numerator/denominator propagation here.
        const double normalizedValue = value / referenceValue;
        const double relativeError = sqrt(
            pow(error / value, 2) + pow(referenceError / referenceValue, 2));
        hist->SetBinContent(bin, normalizedValue);
        hist->SetBinError(bin, normalizedValue * relativeError);
      }

      normalizedSeries[entry.first] = hist;
    }

    return normalizedSeries;
  };

  // Purpose: draw the per-alpha ratio overlay used to inspect the stored pTref
  // contract.
  auto saveRatioOverlay = [&](const PreparedInput &prepared,
                              const map<int, TH1D *> &ratioSeries,
                              const vector<double> &alphaEdges,
                              const string &inputPdfFolder) {
    TCanvas canvas(Form("cRatioOverlay_%s", prepared.token.Data()),
                   Form("cRatioOverlay_%s", prepared.token.Data()), 900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hRatioOverlay_%s", prepared.token.Data()),
        Form(";p_{T}^{ref} (GeV);%s", responseAxisLabel(prepared).Data()), 100,
        prepared.refRatioHist->GetXaxis()->GetBinLowEdge(1),
        prepared.refRatioHist->GetXaxis()->GetBinLowEdge(
            prepared.refRatioHist->GetNbinsX() + 1));
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TLegend legend(0.56, 0.55, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.028);

    int colorIndex = 0;
    for (const auto &entry : ratioSeries) {
      TH1D *hist = entry.second;
      hist->SetMarkerStyle(kFullCircle);
      hist->SetMarkerSize(0.8);
      hist->SetMarkerColor(colors[colorIndex % colorCount]);
      hist->SetLineColor(colors[colorIndex % colorCount]);
      hist->Draw("PE1 SAME");

      const int alphaBin = entry.first;
      const double alphaCut =
          (alphaBin >= 1 && alphaBin < (int)alphaEdges.size())
              ? alphaEdges[alphaBin]
              : alphaBin;
      legend.AddEntry(hist, Form("#alpha < %.2f", alphaCut), "PE");
      ++colorIndex;
    }

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    if (prepared.haveEtaRange) {
      label.DrawLatex(
          0.18, 0.81,
          Form("%.3f < #eta < %.3f", prepared.etaMin, prepared.etaMax));
    }
    label.DrawLatex(0.18, 0.76,
                    Form("%s = B_{%s}^{Data} / B_{%s}^{MC}",
                         responseAxisLabel(prepared).Data(),
                         prepared.label.Data(), prepared.label.Data()));
    legend.Draw();
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_ratio_overlay_vsptref.png",
                       inputPdfFolder.c_str(), prepared.token.Data()));

    delete frame;
  };

  // Purpose: draw the reference-alpha raw MC/Data response used in the alpha
  // fit.
  auto saveRefRawOverlay = [&](const PreparedInput &prepared, TH1D *rawMc,
                               TH1D *rawDt, const string &inputPdfFolder) {
    TCanvas canvas(Form("cRawOverlay_%s", prepared.token.Data()),
                   Form("cRawOverlay_%s", prepared.token.Data()), 900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hRawOverlay_%s", prepared.token.Data()),
        Form(";p_{T}^{ref} (GeV);%s", balanceAxisLabel(prepared).Data()), 100,
        rawMc->GetXaxis()->GetBinLowEdge(1),
        rawMc->GetXaxis()->GetBinLowEdge(rawMc->GetNbinsX() + 1));
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.3);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    rawMc->SetMarkerStyle(kFullCircle);
    rawMc->SetMarkerColor(kBlue + 1);
    rawMc->SetLineColor(kBlue + 1);
    rawMc->Draw("PE1 SAME");

    rawDt->SetMarkerStyle(kFullCircle);
    rawDt->SetMarkerColor(kRed + 1);
    rawDt->SetLineColor(kRed + 1);
    rawDt->Draw("PE1 SAME");

    TLegend legend(0.60, 0.74, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);
    legend.AddEntry(rawMc, "MC", "PE");
    legend.AddEntry(rawDt, "Data", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("#alpha_{ref} < %.2f (bin %d)", prepared.refAlphaValue,
                         refAlphaBin));
    label.DrawLatex(0.18, 0.76,
                    Form("B_{%s}^{Data,MC} = p_{T}^{jet} / p_{T}^{ref}",
                         prepared.label.Data()));
    if (prepared.haveEtaRange) {
      label.DrawLatex(
          0.18, 0.71,
          Form("%.3f < #eta < %.3f", prepared.etaMin, prepared.etaMax));
    }
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_refalpha_raw_vsptref.png",
                       inputPdfFolder.c_str(), prepared.token.Data()));

    delete frame;
  };

  // Purpose: draw the alpha fit in one pTref bin and keep the fitted intercept
  // visible.
  auto saveAlphaFitPlot = [&](const PreparedInput &prepared, int ptBin,
                              double ptLow, double ptHigh, TGraph *graph,
                              TF1 *fitLine, double kfsr, double kfsrError,
                              const string &alphaFolder) {
    TCanvas canvas(Form("cAlphaFit_%s_%d", prepared.token.Data(), ptBin),
                   Form("cAlphaFit_%s_%d", prepared.token.Data(), ptBin), 800,
                   600);
    TH1D *frame =
        new TH1D(Form("hAlphaFit_%s_%d", prepared.token.Data(), ptBin),
                 Form(";#alpha;%s", alphaRatioAxisLabel(prepared).Data()), 100,
                 0.0, 0.5);
    frame->SetMinimum(0.97);
    frame->SetMaximum(1.05);
    styleFrameAxes(frame, 1.30);
    frame->Draw();

    TLine unity(0.0, 1.0, 0.5, 1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlue + 1);
    graph->SetLineColor(kBlue + 1);
    graph->Draw("PZ SAME");

    fitLine->SetLineColor(kRed + 1);
    fitLine->SetLineWidth(2);
    fitLine->Draw("SAME");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("p_{T}^{ref}: %.0f-%.0f GeV", ptLow, ptHigh));
    label.DrawLatex(0.18, 0.76,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    label.DrawLatex(0.18, 0.71,
                    Form("%s = %s(#alpha) / %s(#alpha_{ref})",
                         alphaRatioAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data()));
    label.DrawLatex(
        0.18, 0.66,
        Form("f(#alpha) = k_{FSR} + m#alpha,  k_{FSR} = %.4f #pm %.4f", kfsr,
             kfsrError));
    label.DrawLatex(0.18, 0.61,
                    Form("m = %.4f #pm %.4f", fitLine->GetParameter(1),
                         fitLine->GetParError(1)));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_kFSR_alphaFit_pt%.0fto%.0f.png",
                       alphaFolder.c_str(), prepared.token.Data(), ptLow,
                       ptHigh));

    delete frame;
  };

  // Purpose: draw the per-input kFSR summary after the alpha fits finish.
  auto saveKfsrSummary = [&](const PreparedInput &prepared,
                             const string &alphaFolder) {
    TCanvas canvas(Form("cKfsr_%s", prepared.token.Data()),
                   Form("cKfsr_%s", prepared.token.Data()), 900, 700);
    canvas.SetLogx();

    TH1D *frame =
        new TH1D(Form("hKfsr_%s", prepared.token.Data()),
                 Form(";p_{T}^{ref} (GeV);k_{FSR}^{%s}", prepared.label.Data()),
                 100, prepared.kfsrHist->GetXaxis()->GetBinLowEdge(1),
                 prepared.kfsrHist->GetXaxis()->GetBinLowEdge(
                     prepared.kfsrHist->GetNbinsX() + 1));
    frame->SetMinimum(0.9);
    frame->SetMaximum(1.1);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    prepared.kfsrHist->SetMarkerStyle(kFullCircle);
    prepared.kfsrHist->SetMarkerColor(kBlue + 1);
    prepared.kfsrHist->SetLineColor(kBlue + 1);
    prepared.kfsrHist->Draw("PE1 SAME");

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_kFSR_vsptref.png", alphaFolder.c_str(),
                       prepared.token.Data()));

    delete frame;
  };

  // Purpose: draw the reference-alpha ratio and the alpha->0 corrected response
  // on the same axis.
  auto saveCorrSummary = [&](const PreparedInput &prepared,
                             const string &alphaFolder) {
    TCanvas canvas(Form("cCorr_%s", prepared.token.Data()),
                   Form("cCorr_%s", prepared.token.Data()), 900, 700);
    canvas.SetLogx();

    TH1D *frame = new TH1D(
        Form("hCorr_%s", prepared.token.Data()),
        Form(";p_{T}^{ref} (GeV);%s", responseAxisLabel(prepared).Data()), 100,
        prepared.corrHist->GetXaxis()->GetBinLowEdge(1),
        prepared.corrHist->GetXaxis()->GetBinLowEdge(
            prepared.corrHist->GetNbinsX() + 1));
    frame->SetMinimum(0.7);
    frame->SetMaximum(1.5);
    frame->GetXaxis()->SetMoreLogLabels(kTRUE);
    frame->GetXaxis()->SetNoExponent(kTRUE);
    styleFrameAxes(frame);
    frame->Draw();

    TLine unity(frame->GetXaxis()->GetXmin(), 1.0, frame->GetXaxis()->GetXmax(),
                1.0);
    unity.SetLineStyle(kDashed);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    prepared.refRatioHist->SetMarkerStyle(kOpenCircle);
    prepared.refRatioHist->SetMarkerColor(kBlue + 1);
    prepared.refRatioHist->SetLineColor(kBlue + 1);
    prepared.refRatioHist->Draw("PE1 SAME");

    prepared.corrHist->SetMarkerStyle(kFullCircle);
    prepared.corrHist->SetMarkerColor(kRed + 1);
    prepared.corrHist->SetLineColor(kRed + 1);
    prepared.corrHist->Draw("PE1 SAME");

    TLegend legend(0.54, 0.74, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.032);
    legend.AddEntry(prepared.refRatioHist,
                    Form("Ref #alpha < %.2f", prepared.refAlphaValue), "PE");
    legend.AddEntry(prepared.corrHist, "#alpha #rightarrow 0", "PE");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.032);
    label.DrawLatex(0.18, 0.86, prepared.label.Data());
    label.DrawLatex(0.18, 0.81,
                    Form("%s(0) = k_{FSR} #times %s(#alpha_{ref})",
                         responseAxisLabel(prepared).Data(),
                         responseAxisLabel(prepared).Data()));
    label.DrawLatex(0.18, 0.76,
                    Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                         prepared.refAlphaValue, fitAlphaMin, fitAlphaMax));
    drawCmsLabel(canvas);
    canvas.SaveAs(Form("%s/L3Res_%s_corr_alpha0_vsptref.png",
                       alphaFolder.c_str(), prepared.token.Data()));

    delete frame;
  };

  vector<PreparedInput> preparedInputs;
  preparedInputs.reserve(inputFiles.size());

  for (size_t inputIndex = 0; inputIndex < inputFiles.size(); ++inputIndex) {
    PreparedInput prepared;
    prepared.sample = sampleTypes[inputIndex];
    prepared.path = inputFiles[inputIndex];
    prepared.label = sampleLabel(prepared.sample);
    prepared.token = sampleToken(prepared.sample, (int)inputIndex + 1);
    prepared.jecTag = sampleTag(prepared.sample);

    const string inputFolder =
        outfolder + "/inputs/" + string(prepared.token.Data());
    const string inputPdfFolder = inputFolder + "/pdf";
    const string alphaFolder = inputFolder + "/alpha_extrap";
    gSystem->mkdir(inputFolder.c_str(), kTRUE);
    gSystem->mkdir(inputPdfFolder.c_str(), kTRUE);
    gSystem->mkdir(alphaFolder.c_str(), kTRUE);

    TFile *inputFile = TFile::Open(prepared.path, "READ");
    if (!inputFile || inputFile->IsZombie()) {
      cout << "ERROR: Cannot open derived input '" << prepared.path << "'."
           << endl;
      return;
    }

    TProfile3D *balance3D = (TProfile3D *)inputFile->Get("balance3D_mc");
    if (!balance3D) {
      cout << "ERROR: Missing required object 'balance3D_mc' in "
           << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    vector<double> alphaEdges(balance3D->GetZaxis()->GetNbins() + 1, 0.0);
    for (int bin = 1; bin <= balance3D->GetZaxis()->GetNbins() + 1; ++bin) {
      alphaEdges[bin - 1] = balance3D->GetZaxis()->GetBinLowEdge(bin);
    }
    if (refAlphaBin < 1 || refAlphaBin >= (int)alphaEdges.size()) {
      cout << "ERROR: Reference alpha bin " << refAlphaBin
           << " is outside the input alpha range for " << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.refAlphaValue = alphaEdges[refAlphaBin];
    prepared.etaMin = balance3D->GetYaxis()->GetBinLowEdge(1);
    prepared.etaMax = balance3D->GetYaxis()->GetBinLowEdge(
        balance3D->GetYaxis()->GetNbins() + 1);
    prepared.haveEtaRange = (prepared.etaMax > prepared.etaMin);

    map<int, TH1D *> ratioSeries =
        loadSeries(inputFile, "ratio_vsptref_alpha", prepared.path);
    map<int, TH1D *> jetRatioSeries =
        loadSeries(inputFile, "ratio_vsjetpt_alpha", prepared.path);
    if (ratioSeries.empty() || jetRatioSeries.empty()) {
      cout
          << "ERROR: Missing required ratio_vsptref or ratio_vsjetpt series in "
          << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }
    map<int, TH1D *> ratioNormSeries =
        loadSeries(inputFile, "ratio_norm_vsptref_alpha", prepared.path, false);
    if (ratioNormSeries.empty()) {
      cout << "INFO: ratio_norm_vsptref_alpha* missing in " << prepared.path
           << "; building conservative normalized series in memory." << endl;
      ratioNormSeries = normalizeSeriesByReference(
          ratioSeries, refAlphaBin, "ratio_norm_vsptref_alpha", prepared.path);
    }
    map<int, TH1D *> jetRatioNormSeries =
        loadSeries(inputFile, "ratio_norm_vsjetpt_alpha", prepared.path, false);
    if (jetRatioNormSeries.empty()) {
      cout << "INFO: ratio_norm_vsjetpt_alpha* missing in " << prepared.path
           << "; building conservative normalized series in memory." << endl;
      jetRatioNormSeries =
          normalizeSeriesByReference(jetRatioSeries, refAlphaBin,
                                     "ratio_norm_vsjetpt_alpha", prepared.path);
    }
    if (ratioNormSeries.empty() || jetRatioNormSeries.empty()) {
      inputFile->Close();
      delete inputFile;
      return;
    }
    if (!ratioSeries.count(refAlphaBin) || !jetRatioSeries.count(refAlphaBin)) {
      cout << "ERROR: Missing the requested reference alpha bin " << refAlphaBin
           << " in " << prepared.path << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.refRatioHist = cloneRequiredHist(
        inputFile, Form("ratio_vsptref_alpha%d", refAlphaBin), prepared.path);
    prepared.refJetRatioHist = cloneRequiredHist(
        inputFile, Form("ratio_vsjetpt_alpha%d", refAlphaBin), prepared.path);
    TH1D *rawMc = cloneRequiredHist(
        inputFile, Form("balance_vsptref_mc_alpha%d", refAlphaBin),
        prepared.path);
    TH1D *rawDt = cloneRequiredHist(
        inputFile, Form("balance_vsptref_data_alpha%d", refAlphaBin),
        prepared.path);
    if (!prepared.refRatioHist || !prepared.refJetRatioHist || !rawMc ||
        !rawDt) {
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.fitWindow = parseFitWindow(
        inputIndex < fitWindowTokens.size() ? fitWindowTokens[inputIndex] : "",
        prepared.refRatioHist);
    saveRatioOverlay(prepared, ratioSeries, alphaEdges, inputPdfFolder);
    saveRefRawOverlay(prepared, rawMc, rawDt, inputPdfFolder);

    prepared.kfsrHist = (TH1D *)prepared.refRatioHist->Clone(
        Form("kFSR_vsptref_%s", prepared.token.Data()));
    prepared.kfsrHist->Reset();
    prepared.kfsrHist->SetDirectory(nullptr);
    prepared.corrHist = (TH1D *)prepared.refRatioHist->Clone(
        Form("corr_vsptref_%s", prepared.token.Data()));
    prepared.corrHist->Reset();
    prepared.corrHist->SetDirectory(nullptr);
    prepared.kfsrJetHist = (TH1D *)prepared.refJetRatioHist->Clone(
        Form("kFSR_vsjetpt_%s", prepared.token.Data()));
    prepared.kfsrJetHist->Reset();
    prepared.kfsrJetHist->SetDirectory(nullptr);
    prepared.corrJetHist = (TH1D *)prepared.refJetRatioHist->Clone(
        Form("corr_vsjetpt_%s", prepared.token.Data()));
    prepared.corrJetHist->Reset();
    prepared.corrJetHist->SetDirectory(nullptr);

    TDirectory *inputDir =
        outputFile->mkdir(Form("input_%s", prepared.token.Data()));
    if (!inputDir) {
      cout << "ERROR: Could not create output directory for input '"
           << prepared.label << "'." << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    inputDir->cd();
    prepared.refRatioHist->Write("ratio_vsptref_refalpha");
    prepared.refJetRatioHist->Write("ratio_vsjetpt_refalpha");
    rawMc->Write("balance_vsptref_mc_refalpha");
    rawDt->Write("balance_vsptref_data_refalpha");

    auto fillCorrectedResponse =
        [&](TH1D *refHist, const map<int, TH1D *> &normalizedAxisRatioSeries,
            TH1D *kfsrHist, TH1D *corrHist, const TString &graphStem,
            bool saveAlphaPlots) {
          for (int ptBin = 1; ptBin <= refHist->GetNbinsX(); ++ptBin) {
            vector<double> alphaAllX;
            vector<double> alphaAllY;
            vector<double> alphaAllErr;
            vector<double> alphaFitX;
            vector<double> alphaFitY;
            vector<double> alphaFitErr;

            const double refValue = refHist->GetBinContent(ptBin);
            const double refError = refHist->GetBinError(ptBin);
            if (!(refValue > 0.0) || !std::isfinite(refValue)) {
              kfsrHist->SetBinContent(ptBin, 0.0);
              kfsrHist->SetBinError(ptBin, 0.0);
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
              continue;
            }

            for (const auto &entry : normalizedAxisRatioSeries) {
              const int alphaBin = entry.first;
              TH1D *hist = entry.second;
              const double value = hist->GetBinContent(ptBin);
              const double error = hist->GetBinError(ptBin);
              if (!(value > 0.0) || !std::isfinite(value))
                continue;

              const double alphaCenter =
                  0.5 * (alphaEdges[alphaBin - 1] + alphaEdges[alphaBin]);
              alphaAllX.push_back(alphaCenter);
              alphaAllY.push_back(value);
              alphaAllErr.push_back(
                  (error > 0.0 && std::isfinite(error)) ? error : 0.0);

              if (alphaBin == refAlphaBin)
                continue;
              if (alphaCenter < fitAlphaMin || alphaCenter > fitAlphaMax)
                continue;
              if (!(error > 0.0) || !std::isfinite(error))
                continue;
              alphaFitX.push_back(alphaCenter);
              alphaFitY.push_back(value);
              alphaFitErr.push_back(error);
            }

            if (alphaAllX.size() < 2 || alphaFitX.size() < 2) {
              kfsrHist->SetBinContent(ptBin, 0.0);
              kfsrHist->SetBinError(ptBin, 0.0);
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
              continue;
            }

            TGraphErrors *alphaGraph = new TGraphErrors((int)alphaAllX.size());
            for (int point = 0; point < (int)alphaAllX.size(); ++point) {
              alphaGraph->SetPoint(point, alphaAllX[point], alphaAllY[point]);
              alphaGraph->SetPointError(point, 0.0, alphaAllErr[point]);
            }
            alphaGraph->SetName(Form("%s_pt%d", graphStem.Data(), ptBin));

            TGraphErrors fitGraph((int)alphaFitX.size());
            for (int point = 0; point < (int)alphaFitX.size(); ++point) {
              fitGraph.SetPoint(point, alphaFitX[point], alphaFitY[point]);
              fitGraph.SetPointError(point, 0.0, alphaFitErr[point]);
            }
            TF1 alphaLine(Form("fAlpha_%s_pt%d", graphStem.Data(), ptBin),
                          "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
            alphaLine.SetParameters(1.0, 0.0);
            fitGraph.Fit(&alphaLine, "QRN");

            const double kfsr = alphaLine.GetParameter(0);
            const double kfsrError = alphaLine.GetParError(0);
            kfsrHist->SetBinContent(ptBin, kfsr);
            kfsrHist->SetBinError(ptBin, kfsrError);

            if (kfsr > 0.0) {
              corrHist->SetBinContent(ptBin, refValue * kfsr);
              // Keep the reference-ratio uncertainty on the corrected response.
              // The fitted kFSR uncertainty is preserved separately in
              // kfsrHist; multiplying it back in here would double-count the
              // same reference denominator that already appears in the
              // normalized alpha series.
              corrHist->SetBinError(ptBin, refError * kfsr);
            } else {
              corrHist->SetBinContent(ptBin, 0.0);
              corrHist->SetBinError(ptBin, 0.0);
            }

            TF1 *alphaLineOut = (TF1 *)alphaLine.Clone(
                Form("fAlpha_%s_pt%d", graphStem.Data(), ptBin));
            alphaGraph->Write();
            alphaLineOut->Write();

            if (saveAlphaPlots) {
              saveAlphaFitPlot(
                  prepared, ptBin, refHist->GetXaxis()->GetBinLowEdge(ptBin),
                  refHist->GetXaxis()->GetBinLowEdge(ptBin + 1), alphaGraph,
                  alphaLineOut, kfsr, kfsrError, alphaFolder);
            }
          }
        };

    fillCorrectedResponse(prepared.refRatioHist, ratioNormSeries,
                          prepared.kfsrHist, prepared.corrHist,
                          "gAlphaNormPtRef", true);
    fillCorrectedResponse(prepared.refJetRatioHist, jetRatioNormSeries,
                          prepared.kfsrJetHist, prepared.corrJetHist,
                          "gAlphaNormJetPt", false);

    prepared.kfsrHist->Write("kFSR_vsptref");
    prepared.corrHist->Write("corr_vsptref");
    prepared.kfsrJetHist->Write("kFSR_vsjetpt");
    prepared.corrJetHist->Write("corr_vsjetpt");
    saveKfsrSummary(prepared, alphaFolder);
    saveCorrSummary(prepared, alphaFolder);

    prepared.displayGraph = buildWindowedGraph(
        prepared.corrHist, prepared.fitWindow, true,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    prepared.fitGraph = buildWindowedGraph(
        prepared.corrHist, prepared.fitWindow, false,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    if (!prepared.fitGraph || prepared.fitGraph->GetN() < 3) {
      cout << "ERROR: Not enough valid pTref points in the selected fit window "
              "for input '"
           << prepared.label << "'." << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }

    prepared.displayGraph->Write("g_ptref_display");
    prepared.fitGraph->Write("g_ptref_fit");
    FitWindow jetWindow = parseFitWindow("", prepared.refJetRatioHist);
    prepared.displayJetGraph = buildWindowedGraph(
        prepared.corrJetHist, jetWindow, true,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    prepared.fitJetGraph = buildWindowedGraph(
        prepared.corrJetHist, jetWindow, true,
        markers[inputIndex % markerCount], colors[inputIndex % colorCount]);
    if (!prepared.fitJetGraph || prepared.fitJetGraph->GetN() < 3) {
      cout << "ERROR: Not enough valid JetPt points for input '"
           << prepared.label << "'." << endl;
      inputFile->Close();
      delete inputFile;
      return;
    }
    prepared.displayJetGraph->Write("g_jetpt_display");
    prepared.fitJetGraph->Write("g_jetpt_fit");
    outputFile->cd();
    preparedInputs.push_back(prepared);

    for (auto &entry : ratioSeries)
      delete entry.second;
    for (auto &entry : jetRatioSeries)
      delete entry.second;
    for (auto &entry : ratioNormSeries)
      delete entry.second;
    for (auto &entry : jetRatioNormSeries)
      delete entry.second;
    delete rawMc;
    delete rawDt;
    inputFile->Close();
    delete inputFile;
  }

  if (preparedInputs.empty()) {
    cout << "ERROR: No valid derived inputs were prepared for the combined L3 "
            "fit."
         << endl;
    outputFile->Close();
    delete outputFile;
    return;
  }

  double fitMinOut = preparedInputs.front().fitWindow.min;
  double fitMaxOut = preparedInputs.front().fitWindow.max;
  double xMin =
      preparedInputs.front().refRatioHist->GetXaxis()->GetBinLowEdge(1);
  double xMax = preparedInputs.front().refRatioHist->GetXaxis()->GetBinLowEdge(
      preparedInputs.front().refRatioHist->GetNbinsX() + 1);
  double jetXMin =
      preparedInputs.front().refJetRatioHist->GetXaxis()->GetBinLowEdge(1);
  double jetXMax =
      preparedInputs.front().refJetRatioHist->GetXaxis()->GetBinLowEdge(
          preparedInputs.front().refJetRatioHist->GetNbinsX() + 1);
  double etaMinOut = preparedInputs.front().etaMin;
  double etaMaxOut = preparedInputs.front().etaMax;
  bool haveEtaRange = preparedInputs.front().haveEtaRange;

  TGraphErrors *combinedGraph = new TGraphErrors();
  int combinedPoint = 0;
  double yMin = 1e9;
  double yMax = -1e9;
  TGraphErrors *combinedJetGraph = new TGraphErrors();
  int combinedJetPoint = 0;

  for (size_t inputIndex = 0; inputIndex < preparedInputs.size();
       ++inputIndex) {
    fitMinOut = std::min(fitMinOut, preparedInputs[inputIndex].fitWindow.min);
    fitMaxOut = std::max(fitMaxOut, preparedInputs[inputIndex].fitWindow.max);
    xMin = std::min(
        xMin,
        preparedInputs[inputIndex].refRatioHist->GetXaxis()->GetBinLowEdge(1));
    xMax = std::max(
        xMax,
        preparedInputs[inputIndex].refRatioHist->GetXaxis()->GetBinLowEdge(
            preparedInputs[inputIndex].refRatioHist->GetNbinsX() + 1));
    jetXMin = std::min(
        jetXMin,
        preparedInputs[inputIndex].refJetRatioHist->GetXaxis()->GetBinLowEdge(
            1));
    jetXMax = std::max(
        jetXMax,
        preparedInputs[inputIndex].refJetRatioHist->GetXaxis()->GetBinLowEdge(
            preparedInputs[inputIndex].refJetRatioHist->GetNbinsX() + 1));

    if (preparedInputs[inputIndex].haveEtaRange) {
      etaMinOut = std::min(etaMinOut, preparedInputs[inputIndex].etaMin);
      etaMaxOut = std::max(etaMaxOut, preparedInputs[inputIndex].etaMax);
      haveEtaRange = true;
    }

    for (int point = 0; point < preparedInputs[inputIndex].fitGraph->GetN();
         ++point) {
      combinedGraph->SetPoint(
          combinedPoint, preparedInputs[inputIndex].fitGraph->GetX()[point],
          preparedInputs[inputIndex].fitGraph->GetY()[point]);
      combinedGraph->SetPointError(
          combinedPoint, preparedInputs[inputIndex].fitGraph->GetEX()[point],
          preparedInputs[inputIndex].fitGraph->GetEY()[point]);
      yMin = std::min(yMin, preparedInputs[inputIndex].fitGraph->GetY()[point]);
      yMax = std::max(yMax, preparedInputs[inputIndex].fitGraph->GetY()[point]);
      ++combinedPoint;
    }

    for (int point = 0; point < preparedInputs[inputIndex].displayGraph->GetN();
         ++point) {
      yMin = std::min(yMin,
                      preparedInputs[inputIndex].displayGraph->GetY()[point]);
      yMax = std::max(yMax,
                      preparedInputs[inputIndex].displayGraph->GetY()[point]);
    }

    for (int point = 0; point < preparedInputs[inputIndex].fitJetGraph->GetN();
         ++point) {
      combinedJetGraph->SetPoint(
          combinedJetPoint,
          preparedInputs[inputIndex].fitJetGraph->GetX()[point],
          preparedInputs[inputIndex].fitJetGraph->GetY()[point]);
      combinedJetGraph->SetPointError(
          combinedJetPoint,
          preparedInputs[inputIndex].fitJetGraph->GetEX()[point],
          preparedInputs[inputIndex].fitJetGraph->GetEY()[point]);
      ++combinedJetPoint;
    }
  }

  if (combinedGraph->GetN() < 3) {
    cout << "ERROR: The combined fit graph has fewer than three points."
         << endl;
    outputFile->Close();
    delete outputFile;
    delete combinedGraph;
    delete combinedJetGraph;
    return;
  }

  if (combinedJetGraph->GetN() < 3) {
    cout << "ERROR: The combined JetPt graph has fewer than three points."
         << endl;
    outputFile->Close();
    delete outputFile;
    delete combinedGraph;
    delete combinedJetGraph;
    return;
  }

  const double yLowerMargin = 0.08;
  const double yUpperMargin = 0.30;
  yMin = std::max(0.45, yMin - yLowerMargin);
  yMax = std::min(1.8, yMax + yUpperMargin);

  SequentialFitResult fitResult = runSequentialTemplateFit(
      combinedGraph, fitMinOut, fitMaxOut, outfilename.c_str());
  if (!fitResult.valid || !fitResult.finalFit) {
    cout << "ERROR: The sequential pTref fit failed." << endl;
    outputFile->Close();
    delete outputFile;
    delete combinedGraph;
    return;
  }

  TCanvas finalCanvas("cL3ResFinal", "cL3ResFinal", 900, 700);
  finalCanvas.SetLogx();
  TH1D *frame = new TH1D("hL3ResFinalFrame",
                         ";p_{T}^{ref} (GeV);R_{Data/MC} = B^{Data}/B^{MC}",
                         100, xMin, xMax);
  frame->SetMinimum(yMin);
  frame->SetMaximum(yMax);
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  styleFrameAxes(frame);
  frame->Draw();

  TLine unity(xMin, 1.0, xMax, 1.0);
  unity.SetLineStyle(kDashed);
  unity.SetLineColor(kGray + 1);
  unity.Draw("SAME");

  TLegend legend(0.52, 0.68, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextFont(42);
  legend.SetTextSize(0.028);

  for (size_t inputIndex = 0; inputIndex < preparedInputs.size();
       ++inputIndex) {
    preparedInputs[inputIndex].displayGraph->Draw("PE SAME");
    legend.AddEntry(preparedInputs[inputIndex].displayGraph,
                    Form("%s (%.0f-%.0f GeV)",
                         preparedInputs[inputIndex].label.Data(),
                         preparedInputs[inputIndex].fitWindow.min,
                         preparedInputs[inputIndex].fitWindow.max),
                    "PE");
  }

  fitResult.finalFit->SetLineColor(kBlack);
  fitResult.finalFit->SetLineWidth(2);
  fitResult.finalFit->Draw("SAME");

  TLatex label;
  label.SetNDC();
  label.SetTextFont(42);
  label.SetTextSize(0.032);
  label.DrawLatex(0.18, 0.86, Form("Run label: %s", runLabel.c_str()));
  if (haveEtaRange) {
    label.DrawLatex(0.18, 0.81,
                    Form("%.3f < #eta < %.3f", etaMinOut, etaMaxOut));
  }
  label.DrawLatex(0.18, 0.76,
                  Form("#alpha_{ref} < %.2f, fit %.2f #leq #alpha #leq %.2f",
                       preparedInputs.front().refAlphaValue, fitAlphaMin,
                       fitAlphaMax));
  label.DrawLatex(0.18, 0.71,
                  Form("Fit window: %.0f-%.0f GeV", fitMinOut, fitMaxOut));
  label.DrawLatex(0.18, 0.66, "Direct global p_{T}^{ref} fit to R_{Data/MC}");
  label.DrawLatex(0.18, 0.61,
                  "C_{L3}(x)=1./(p_{0}+p_{1}/x+p_{2}log(x)/"
                  "x+p_{3}T(x)+p_{6}x^{-0.3051}+p_{7}x)");
  label.DrawLatex(0.18, 0.56,
                  "T(x)=((x/p_{4})^{p_{5}}-1)/((x/p_{4})^{p_{5}}+1)");
  label.DrawLatex(0.18, 0.51,
                  Form("p_{0}=%.4f, p_{1}=%.4f, p_{2}=%.4f, p_{3}=%.4f",
                       fitResult.finalFit->GetParameter(0),
                       fitResult.finalFit->GetParameter(1),
                       fitResult.finalFit->GetParameter(2),
                       fitResult.finalFit->GetParameter(3)));
  label.DrawLatex(0.18, 0.46,
                  Form("p_{4}=%.4f, p_{5}=%.4f, p_{6}=%.4f, p_{7}=%.6f",
                       fitResult.finalFit->GetParameter(4),
                       fitResult.finalFit->GetParameter(5),
                       fitResult.finalFit->GetParameter(6),
                       fitResult.finalFit->GetParameter(7)));
  label.DrawLatex(0.18, 0.41,
                  Form("#chi^{2}/ndf = %.1f/%d",
                       fitResult.finalFit->GetChisquare(),
                       fitResult.finalFit->GetNDF()));
  legend.Draw();
  drawCmsLabel(finalCanvas);
  finalCanvas.SaveAs(
      Form("%s/L3Res_%s_ptref_final.png", pdfFolder.c_str(), runLabel.c_str()));

  const string outputTag = (preparedInputs.size() > 1
                                ? string("combined")
                                : string(preparedInputs.front().jecTag.Data()));
  outputFile->cd();
  combinedGraph->SetName("g_ptref_combined");
  combinedGraph->Write();
  combinedJetGraph->SetName("g_jetpt_combined");
  combinedJetGraph->Write();
  fitResult.f0->Write("f_ptref_f0");
  fitResult.f1->Write("f_ptref_f1");
  fitResult.f2->Write("f_ptref_f2");
  fitResult.f3->Write("f_ptref_f3");
  fitResult.f4->Write("f_ptref_f4");
  fitResult.finalFit->Write("f_ptref_final");

  TNamed runLabelObj("run_label", runLabel.c_str());
  runLabelObj.Write();
  TNamed outputNameObj("output_name", outfilename.c_str());
  outputNameObj.Write();
  TNamed outputTagObj("output_tag", outputTag.c_str());
  outputTagObj.Write();
  TNamed lumiLabelObj("lumi_label", lumiLabel.c_str());
  lumiLabelObj.Write();
  TParameter<int> refAlphaBinObj("ref_alpha_bin", refAlphaBin);
  refAlphaBinObj.Write();
  TParameter<double> refAlphaValueObj("ref_alpha_value",
                                      preparedInputs.front().refAlphaValue);
  refAlphaValueObj.Write();
  TParameter<double> fitMinObj("pt_fit_min", fitMinOut);
  fitMinObj.Write();
  TParameter<double> fitMaxObj("pt_fit_max", fitMaxOut);
  fitMaxObj.Write();
  TParameter<double> etaMinObj("eta_min", haveEtaRange ? etaMinOut : -5.191);
  etaMinObj.Write();
  TParameter<double> etaMaxObj("eta_max", haveEtaRange ? etaMaxOut : 5.191);
  etaMaxObj.Write();
  TParameter<int> inputCountObj("n_inputs", (int)preparedInputs.size());
  inputCountObj.Write();

  outputFile->Close();

  cout << "============================================" << endl;
  cout << "Wrote fit ROOT file: " << outfolder << "/" << outfilename
       << "_fit.root" << endl;
  cout << "Plots saved in: " << pdfFolder << endl;
  cout << "NOTE: L3Res.C stops at the p_{T}^{ref} fit stage." << endl;
  cout << "      Run createL2L3ResTextFile.C on the fit ROOT file, or use "
          "dofits_L3.C, to produce direct p_{T}^{ref} text payloads."
       << endl;
  cout << "============================================" << endl;
}