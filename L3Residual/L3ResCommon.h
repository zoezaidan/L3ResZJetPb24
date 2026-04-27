#ifndef L3RESIDUAL_L3RESCOMMON_H
#define L3RESIDUAL_L3RESCOMMON_H

#include "TGraphErrors.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

enum class SampleType {
  PhotonJet = 0,
  ZJet = 1,
};

struct FitWindow {
  double min = 0.0;
  double max = 0.0;
};

struct PreparedInput {
  SampleType sample = SampleType::PhotonJet;
  TString path;
  TString label;
  TString token;
  TString jecTag;
  FitWindow fitWindow;
  double refAlphaValue = 0.0;
  double etaMin = 0.0;
  double etaMax = 0.0;
  bool haveEtaRange = false;
  TH1D *refRatioHist = nullptr;
  TH1D *refJetRatioHist = nullptr;
  TH1D *corrHist = nullptr;
  TH1D *corrJetHist = nullptr;
  TH1D *kfsrHist = nullptr;
  TH1D *kfsrJetHist = nullptr;
  TGraphErrors *displayGraph = nullptr;
  TGraphErrors *fitGraph = nullptr;
  TGraphErrors *displayJetGraph = nullptr;
  TGraphErrors *fitJetGraph = nullptr;
};

// Purpose: split a comma-separated macro argument into clean tokens.
inline std::vector<TString> splitCsvTokens(const TString &csv) {
  std::vector<TString> tokens;
  TString cleaned = csv;
  cleaned.ReplaceAll("\n", ",");
  cleaned.ReplaceAll("\t", ",");
  cleaned.ReplaceAll(";", ",");

  TObjArray *parts = cleaned.Tokenize(",");
  if (!parts)
    return tokens;

  for (int index = 0; index < parts->GetEntriesFast(); ++index) {
    TObject *obj = parts->At(index);
    if (!obj)
      continue;
    TString token = obj->GetName();
    token = token.Strip(TString::kBoth);
    if (token.Length() > 0)
      tokens.push_back(token);
  }

  parts->Delete();
  delete parts;
  return tokens;
}

// Purpose: parse a user-supplied sample token into the active enum.
inline bool parseSampleToken(const TString &token, SampleType &sample) {
  TString lower = token;
  lower.ToLower();
  if (lower == "photonjet" || lower == "photon+jet" || lower == "gammajet" ||
      lower == "gamma+jet") {
    sample = SampleType::PhotonJet;
    return true;
  }
  if (lower == "zjet" || lower == "z+jet") {
    sample = SampleType::ZJet;
    return true;
  }
  return false;
}

// Purpose: parse and validate the explicit sample list passed to the fit macro.
inline bool parseSampleTypesCSV(const TString &csv, std::size_t expectedCount,
                                std::vector<SampleType> &samples,
                                std::string &errorMessage) {
  samples.clear();
  const std::vector<TString> tokens = splitCsvTokens(csv);
  if (tokens.size() != expectedCount) {
    errorMessage = "Expected " + std::to_string(expectedCount) +
                   " sample token(s) in sampleTypesCSV.";
    return false;
  }

  for (const TString &token : tokens) {
    SampleType sample;
    if (!parseSampleToken(token, sample)) {
      errorMessage = "Unsupported sample token '" + std::string(token.Data()) +
                     "'. Use photonjet or zjet.";
      samples.clear();
      return false;
    }
    samples.push_back(sample);
  }
  return true;
}

// Purpose: provide the canonical plot label for one sample.
inline TString sampleLabel(SampleType sample) {
  switch (sample) {
  case SampleType::PhotonJet:
    return "#gamma+jet";
  case SampleType::ZJet:
    return "Z+jet";
  }
  return "Input";
}

// Purpose: provide the canonical short tag used in text-file names.
inline TString sampleTag(SampleType sample) {
  switch (sample) {
  case SampleType::PhotonJet:
    return "photonjet";
  case SampleType::ZJet:
    return "zjet";
  }
  return "input";
}

// Purpose: provide the canonical per-sample folder token.
inline TString sampleToken(SampleType sample, int indexOneBased) {
  switch (sample) {
  case SampleType::PhotonJet:
    return Form("%02d_photonplusjet", indexOneBased);
  case SampleType::ZJet:
    return Form("%02d_zplusjet", indexOneBased);
  }
  return Form("%02d_input", indexOneBased);
}

// Purpose: convert a textual fit window like 60-400 into numeric limits.
inline FitWindow parseFitWindow(const TString &token,
                                const TH1D *referenceHist) {
  FitWindow window;
  if (!referenceHist)
    return window;

  window.min = referenceHist->GetXaxis()->GetBinLowEdge(1);
  window.max =
      referenceHist->GetXaxis()->GetBinLowEdge(referenceHist->GetNbinsX() + 1);

  TString cleaned = token;
  cleaned.ReplaceAll(" ", "");
  if (!cleaned.Contains("-"))
    return window;

  const int dashIndex = cleaned.First('-');
  const double userMin = TString(cleaned(0, dashIndex)).Atof();
  const double userMax =
      TString(cleaned(dashIndex + 1, cleaned.Length() - dashIndex - 1)).Atof();

  if (userMin > 0.0)
    window.min = std::max(window.min, userMin);
  if (userMax > 0.0)
    window.max = std::min(window.max, userMax);
  if (window.max <= window.min) {
    window.min = referenceHist->GetXaxis()->GetBinLowEdge(1);
    window.max = referenceHist->GetXaxis()->GetBinLowEdge(
        referenceHist->GetNbinsX() + 1);
  }
  return window;
}

// Purpose: turn a correction histogram into a TGraphErrors used for drawing or
// fitting.
inline TGraphErrors *buildWindowedGraph(const TH1D *hist,
                                        const FitWindow &window,
                                        bool useFullRange, int marker,
                                        int color) {
  if (!hist)
    return nullptr;

  TGraphErrors *graph = new TGraphErrors();
  int point = 0;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double lowEdge = hist->GetXaxis()->GetBinLowEdge(bin);
    const double highEdge = hist->GetXaxis()->GetBinLowEdge(bin + 1);
    const double value = hist->GetBinContent(bin);
    const double error = hist->GetBinError(bin);

    if (!useFullRange && (lowEdge < window.min || highEdge > window.max))
      continue;
    if (!(value > 0.0) || !std::isfinite(value))
      continue;
    if (!(error > 0.0) || !std::isfinite(error))
      continue;
    if (value < 0.5 || value > 1.5)
      continue;

    graph->SetPoint(point, hist->GetBinCenter(bin), value);
    graph->SetPointError(point, 0.5 * hist->GetBinWidth(bin), error);
    ++point;
  }

  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  return graph;
}

// Purpose: return the mean y value of a graph for stable fit seeds.
inline double graphMeanY(const TGraphErrors *graph) {
  if (!graph || graph->GetN() <= 0)
    return 1.0;
  double sum = 0.0;
  for (int point = 0; point < graph->GetN(); ++point)
    sum += graph->GetY()[point];
  return sum / graph->GetN();
}

#endif