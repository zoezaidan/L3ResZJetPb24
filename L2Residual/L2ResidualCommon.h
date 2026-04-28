#pragma once

#include "TAxis.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace l2residual {

struct RunReport {
  std::string logPath;
  std::vector<std::string> warnings;
  std::vector<std::string> errors;
  std::vector<std::string> fitAlerts;
};

inline RunReport *&activeRunReport() {
  static RunReport *report = nullptr;
  return report;
}

inline void recordWarning(const std::string &message) {
  if (auto *report = activeRunReport()) {
    report->warnings.push_back(message);
    std::cout << "WARNING: " << message << std::endl;
    return;
  }
  std::cout << "WARNING: " << message << std::endl;
}

inline void recordError(const std::string &message) {
  if (auto *report = activeRunReport()) {
    report->errors.push_back(message);
    std::cerr << "ERROR: " << message << std::endl;
    return;
  }
  std::cerr << "ERROR: " << message << std::endl;
}

inline void recordFitAlert(const std::string &message) {
  if (auto *report = activeRunReport()) {
    report->fitAlerts.push_back(message);
    std::cout << "FIT ALERT: " << message << std::endl;
    return;
  }
  std::cout << "FIT ALERT: " << message << std::endl;
}

inline std::string runLogFileName() { return "runL2RES.log"; }

inline std::string runSummaryFileName() { return "runL2RES_summary.txt"; }

inline double fitAlertChi2NdfThreshold() { return 5.0; }

inline int minimumFitPoints() { return 2; }

inline double textPtMin() { return 30.0; }

inline double textPtMax() { return 1000.0; }

inline double binnedTextEtaMax() { return 2.9; }

inline double parametrizedTextEtaMax() { return 2.5; }

inline double binnedTextLowPtFloor() { return 0.001; }

inline double binnedTextLowPtThreshold() { return 25.0; }

inline std::string joinPath(const std::string &left, const std::string &right);
inline std::string absolutePath(const std::string &path);
inline std::string parentDirectory(const std::string &path);
inline std::string defaultOutputDirectory(const std::string &inputPath,
                                          const std::string &fallback = ".");

inline std::string edgeLabel(double x) {
  const double xi = std::round(x);
  if (std::fabs(x - xi) < 1e-6)
    return std::to_string(static_cast<int>(xi));
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

inline std::string defaultSubdir() { return "hibin_70.0_100.0/eta_-5.2_5.2"; }

struct RunConfig {
  std::string selectionSubdir;
  std::string runDir;
  std::string mcFile;
  std::string dataFile;
  std::string plotDirName;
  std::string kfsrDirName;
  std::string ptParamDirName;
  std::string logsDirName;
  std::string textDirName;
  int alphaBin;
  bool plotPerAlphaDijetAsymmetry;
  bool useAbs;
  bool useWideAbs;
};

inline RunConfig defaultRunConfig() {
  return {defaultSubdir(),
          absolutePath("test_output/L2Residual_OO_2026_04_21"),
          absolutePath("/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/"
                       "L2ResOO/Template_OO_pythia.root"),
          absolutePath("/eos/cms/store/group/phys_heavyions/bharikri/JetMinPOG/"
                       "L2ResOO/Template_OO_Data.root"),
          "plotresponses",
          "kfsr",
          "ptparam",
          "logs",
          "text",
          5,
          false,
          false,
          true};
}

inline const std::vector<std::string> &centralityBins() {
  static const std::vector<std::string> bins = {
      "hibin_50.0_100.0", "hibin_50.0_80.0", "hibin_70.0_100.0",
      "hibin_70.0_90.0",  "hibin_70.0_80.0", "hibin_80.0_100.0",
      "hibin_80.0_90.0",
  };
  return bins;
}

inline const std::vector<std::string> &etaBins() {
  static const std::vector<std::string> bins = {"eta_-5.2_5.2"};
  return bins;
}

inline std::vector<std::string> selectionSubdirs() {
  std::vector<std::string> subdirs;
  for (const auto &centbin : centralityBins()) {
    for (const auto &etabin : etaBins()) {
      subdirs.push_back(centbin + "/" + etabin);
    }
  }
  return subdirs;
}

struct CorrBin {
  double lowPt;
  double highPt;
  TH1D *hist;
};

struct PtLabelRange {
  std::string label;
  double lowPt;
  double highPt;
};

inline bool hasSelectionSubdir(TFile *file, const std::string &subdir) {
  if (!file || subdir.empty())
    return false;
  return file->GetDirectory(subdir.c_str()) != nullptr;
}

inline bool validateSelectionSubdir(TFile *file, const std::string &subdir,
                                    const std::string &context) {
  if (!file) {
    recordError("Null input file in " + context);
    return false;
  }
  if (subdir.empty()) {
    recordError("Empty selection subdir in " + context);
    return false;
  }
  if (hasSelectionSubdir(file, subdir)) {
    return true;
  }

  std::ostringstream message;
  message << "Missing selection subdir " << subdir << " in " << file->GetName()
          << " for " << context;
  recordError(message.str());
  return false;
}

inline bool validateSelectionSubdir(const TString &fileName,
                                    const std::string &subdir,
                                    const std::string &context) {
  TFile *file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::ostringstream message;
    message << "Could not open " << fileName << " in " << context;
    recordError(message.str());
    if (file)
      file->Close();
    return false;
  }

  const bool ok = validateSelectionSubdir(file, subdir, context);
  file->Close();
  return ok;
}

inline bool validateSelectionSubdir(TFile *first, TFile *second,
                                    const std::string &subdir,
                                    const std::string &context) {
  const bool firstOk = validateSelectionSubdir(first, subdir, context);
  const bool secondOk = validateSelectionSubdir(second, subdir, context);
  return firstOk && secondOk;
}

inline std::string selectionPath(const std::string &subdir,
                                 const std::string &name) {
  return subdir.empty() ? name : subdir + "/" + name;
}

inline std::string alphaSuffix(double alpha) {
  return Form("_alpha%.1f", alpha);
}

inline std::vector<double> axisEdges(const TAxis *axis) {
  std::vector<double> edges;
  if (!axis)
    return edges;

  edges.reserve(axis->GetNbins() + 1);
  for (int bin = 1; bin <= axis->GetNbins(); ++bin) {
    edges.push_back(axis->GetBinLowEdge(bin));
  }
  edges.push_back(axis->GetBinUpEdge(axis->GetNbins()));
  return edges;
}

inline TH1D *makeHistFromAxis(const std::string &name, const TAxis *axis) {
  const auto edges = axisEdges(axis);
  if (!axis || edges.size() < 2)
    return nullptr;
  return new TH1D(name.c_str(), "  ; ;", axis->GetNbins(), edges.data());
}

inline std::string formatEtaBinLabel(const TAxis *axis, int bin,
                                     bool absolute) {
  if (!axis || bin < 1 || bin > axis->GetNbins())
    return absolute ? "|#eta| bin unavailable" : "#eta bin unavailable";
  const char *axisTitle = absolute ? "|#eta|" : "#eta";
  return Form("%.3f < %s < %.3f", axis->GetBinLowEdge(bin), axisTitle,
              axis->GetBinUpEdge(bin));
}

inline bool parsePtLabel(const std::string &label, double &lowPt,
                         double &highPt) {
  const size_t sep = label.find("to");
  if (sep == std::string::npos || sep == 0 || sep + 2 >= label.size())
    return false;

  try {
    lowPt = std::stod(label.substr(0, sep));
    highPt = std::stod(label.substr(sep + 2));
  } catch (...) {
    return false;
  }

  return true;
}

inline std::vector<PtLabelRange>
collectPtLabelRanges(TFile *file, const std::string &subdir,
                     const std::string &prefix,
                     const std::string &suffix = "") {
  std::vector<PtLabelRange> ranges;
  if (!file)
    return ranges;

  TDirectory *dir = subdir.empty() ? file : file->GetDirectory(subdir.c_str());
  if (!dir)
    return ranges;

  std::set<std::string> seen;
  TIter next(dir->GetListOfKeys());
  while (auto *obj = next()) {
    auto *key = dynamic_cast<TKey *>(obj);
    if (!key)
      continue;

    const std::string name = key->GetName();
    if (name.rfind(prefix, 0) != 0)
      continue;

    std::string label;
    if (!suffix.empty()) {
      if (name.size() <= prefix.size() + suffix.size())
        continue;
      if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
        continue;
      label = name.substr(prefix.size(),
                          name.size() - prefix.size() - suffix.size());
    } else {
      label = name.substr(prefix.size());
    }

    if (label.empty() || seen.count(label))
      continue;

    double lowPt = 0.;
    double highPt = 0.;
    if (!parsePtLabel(label, lowPt, highPt))
      continue;

    seen.insert(label);
    ranges.push_back({label, lowPt, highPt});
  }

  std::sort(ranges.begin(), ranges.end(),
            [](const PtLabelRange &left, const PtLabelRange &right) {
              if (left.lowPt != right.lowPt)
                return left.lowPt < right.lowPt;
              return left.highPt < right.highPt;
            });

  return ranges;
}

inline std::string joinPath(const std::string &left, const std::string &right) {
  if (left.empty())
    return right;
  if (right.empty())
    return left;
  return gSystem->ConcatFileName(left.c_str(), right.c_str());
}

inline std::string absolutePath(const std::string &path) {
  if (path.empty() || gSystem->IsAbsoluteFileName(path.c_str()))
    return path;

  const char *logicalPwd = gSystem->Getenv("PWD");
  if (logicalPwd && logicalPwd[0] != '\0') {
    return joinPath(logicalPwd, path);
  }

  return joinPath(gSystem->WorkingDirectory(), path);
}

inline std::string parentDirectory(const std::string &path) {
  if (path.empty())
    return {};

  const char *dirname = gSystem->DirName(path.c_str());
  if (!dirname)
    return {};

  std::string result = dirname;

  if (result == "." || result.empty()) {
    return {};
  }
  return result;
}

inline std::string defaultOutputDirectory(const std::string &inputPath,
                                          const std::string &fallback) {
  const std::string parent = parentDirectory(inputPath);
  if (!parent.empty()) {
    return parent;
  }
  return fallback;
}

inline double alphaFitMin() { return 0.0; }

inline double alphaFitMax() { return 0.4; }

inline double referenceAlpha() { return 0.3; }

inline double ptFitMin() { return 60.0; }

inline double ptFitMax() { return 300.0; }

inline std::string cmsExtraText() { return "Internal"; }

inline std::string cmsCollisionText() { return "OO 2025"; }

inline std::string cmsLumiText() { return "6.1 nb^{-1}"; }

inline bool isReferenceAlpha(double alpha) {
  return std::fabs(alpha - referenceAlpha()) < 1e-6;
}

inline std::string referenceAlphaSelectionText() {
  return Form("Reference: #alpha = %.2f excluded", referenceAlpha());
}

inline std::string referenceAlphaCutText() {
  return Form("#alpha < %.2f", referenceAlpha());
}

inline std::string alphaFitRangeText() {
  return Form("Fit: %.2f #leq #alpha #leq %.2f", alphaFitMin(), alphaFitMax());
}

inline std::string ptFitRangeText() {
  return Form("Fit: %.0f #leq p_{T}^{avg} #leq %.0f GeV", ptFitMin(),
              ptFitMax());
}

inline void writeRunSummary(const RunReport &report,
                            const std::string &summaryPath) {
  std::ofstream out(summaryPath);
  if (!out)
    return;

  out << "runL2RES log: " << report.logPath << '\n';
  out << "errors: " << report.errors.size() << '\n';
  for (const auto &message : report.errors)
    out << "ERROR: " << message << '\n';
  out << "warnings: " << report.warnings.size() << '\n';
  for (const auto &message : report.warnings)
    out << "WARNING: " << message << '\n';
  out << "fit alerts: " << report.fitAlerts.size() << '\n';
  for (const auto &message : report.fitAlerts)
    out << "FIT ALERT: " << message << '\n';
}

inline void printRunSummary(const RunReport &report,
                            const std::string &summaryPath) {
  std::cout << "runL2RES log: " << report.logPath << std::endl;
  std::cout << "runL2RES summary: " << summaryPath << std::endl;

  if (report.errors.empty() && report.warnings.empty() &&
      report.fitAlerts.empty()) {
    std::cout << "No warnings, errors, or fit alerts." << std::endl;
    return;
  }

  for (const auto &message : report.errors)
    std::cout << "ERROR: " << message << std::endl;
  for (const auto &message : report.warnings)
    std::cout << "WARNING: " << message << std::endl;
  for (const auto &message : report.fitAlerts)
    std::cout << "FIT ALERT: " << message << std::endl;
}

inline const char *asymmetryAxisTitle() { return "Asymmetry A"; }

inline const char *asymmetryRatioTitle() { return "A_{MC}/A_{Data}"; }

inline const char *normalizedAsymmetryRatioTitle() { return "r_{A}(#alpha)"; }

inline const char *asymmetryPtCorrectionTitle() { return "c_{A}(p_{T}^{avg})"; }

inline const char *asymmetryEtaCorrectionTitle() { return "C_{A}(#eta)"; }

inline std::string sanitizeSubdir(std::string subdir) {
  for (char &ch : subdir) {
    if (ch == '/')
      ch = '_';
  }
  return subdir;
}

inline std::vector<std::string> splitSubdir(const std::string &subdir) {
  std::vector<std::string> parts;
  std::stringstream ss(subdir);
  std::string token;
  while (std::getline(ss, token, '/')) {
    if (!token.empty())
      parts.push_back(token);
  }
  return parts;
}

inline std::string selectionTag(const std::string &subdir) {
  return sanitizeSubdir(subdir.empty() ? defaultSubdir() : subdir);
}

inline std::string formatCentLabel(const std::string &subdir) {
  const auto parts = splitSubdir(subdir);
  if (parts.empty())
    return "";
  double low = 0.;
  double high = 0.;
  if (std::sscanf(parts[0].c_str(), "hibin_%lf_%lf", &low, &high) == 2) {
    return Form("Centrality %.0f-%.0f", low, high);
  }
  return parts[0];
}

inline std::string formatEtaSelectionLabel(const std::string &subdir) {
  const auto parts = splitSubdir(subdir);
  if (parts.size() < 2)
    return "";
  double low = 0.;
  double high = 0.;
  if (std::sscanf(parts[1].c_str(), "eta_%lf_%lf", &low, &high) == 2) {
    return Form("Selection: %.3f < #eta < %.3f", low, high);
  }
  return parts[1];
}

inline void applyStyle() {
  static bool initialized = false;
  if (initialized)
    return;
  initialized = true;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.12);
}

inline void styleHistAxis(TH1 *h, const char *xTitle, const char *yTitle) {
  if (!h)
    return;
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetXaxis()->SetTitleSize(0.046);
  h->GetYaxis()->SetTitleSize(0.046);
  h->GetXaxis()->SetLabelSize(0.034);
  h->GetYaxis()->SetLabelSize(0.034);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleOffset(1.30);
}

inline void styleLogxAxis(TH1 *h, const char *xTitle, const char *yTitle) {
  if (!h)
    return;
  styleHistAxis(h, xTitle, yTitle);
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNoExponent(kTRUE);
  h->GetXaxis()->SetNdivisions(510);
}

inline void styleLegend(TLegend *legend, double textSize = 0.028,
                        int columns = 1) {
  if (!legend)
    return;
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(textSize);
  legend->SetNColumns(columns);
}

inline int seriesColor(int index) {
  static const int colors[] = {
      kBlue + 1, kRed + 1,    kGreen + 2, kMagenta + 1, kOrange + 7,
      kCyan + 2, kViolet + 1, kTeal + 2,  kPink + 7,    kBlack};
  return colors[index % (sizeof(colors) / sizeof(colors[0]))];
}

inline int seriesMarker(int index) {
  static const int markers[] = {20, 21, 22, 23, 33, 29, 24, 25, 26, 27};
  return markers[index % (sizeof(markers) / sizeof(markers[0]))];
}

inline void drawCmsLikeLabels(const std::string &centLabel,
                              const std::string &line1 = "",
                              const std::string &line2 = "",
                              const std::string &line3 = "") {
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(61);
  latex.SetTextSize(0.05);
  latex.DrawLatex(0.16, 0.93, "CMS");

  latex.SetTextFont(52);
  latex.SetTextSize(0.038);
  latex.DrawLatex(0.27, 0.93, cmsExtraText().c_str());

  latex.SetTextFont(42);
  latex.SetTextSize(0.033);
  latex.DrawLatex(
      0.60, 0.93,
      Form("%s (%s)", cmsLumiText().c_str(), cmsCollisionText().c_str()));

  if (!centLabel.empty())
    latex.DrawLatex(0.18, 0.86, centLabel.c_str());
  if (!line1.empty())
    latex.DrawLatex(0.18, 0.81, line1.c_str());
  if (!line2.empty())
    latex.DrawLatex(0.18, 0.76, line2.c_str());
  if (!line3.empty())
    latex.DrawLatex(0.18, 0.71, line3.c_str());
}

inline TH1D *getH1(TFile *file, const std::string &name,
                   const std::string &subdir = "") {
  if (!file)
    return nullptr;
  return dynamic_cast<TH1D *>(file->Get(selectionPath(subdir, name).c_str()));
}

inline TDirectory *ensureSubdir(TDirectory *root, const std::string &subdir) {
  if (!root || subdir.empty())
    return root;

  TDirectory *current = root;
  std::stringstream ss(subdir);
  std::string token;
  while (std::getline(ss, token, '/')) {
    if (token.empty())
      continue;
    TDirectory *next = current->GetDirectory(token.c_str());
    if (!next)
      next = current->mkdir(token.c_str());
    current = next;
  }
  return current;
}

inline std::string ensureDirectory(const std::string &path) {
  gSystem->mkdir(path.c_str(), true);
  return path;
}

inline TH1D *loadKFactors(TFile *file, const std::string &subdir) {
  if (!file)
    return nullptr;
  if (auto *h = dynamic_cast<TH1D *>(
          file->Get(selectionPath(subdir, "kfactors").c_str())))
    return h;
  return dynamic_cast<TH1D *>(
      file->Get(selectionPath(subdir, "ratio").c_str()));
}

inline TF1 *loadParam(TFile *file, const std::string &subdir,
                      const std::string &name) {
  if (!file)
    return nullptr;
  return dynamic_cast<TF1 *>(file->Get(selectionPath(subdir, name).c_str()));
}

inline bool loadStoredFitRange(TFile *file, const std::string &subdir,
                               double &fitMin, double &fitMax) {
  if (!file)
    return false;

  auto *minParam = dynamic_cast<TParameter<double> *>(
      file->Get(selectionPath(subdir, "pt_fit_min").c_str()));
  auto *maxParam = dynamic_cast<TParameter<double> *>(
      file->Get(selectionPath(subdir, "pt_fit_max").c_str()));
  if (!minParam || !maxParam)
    return false;

  fitMin = minParam->GetVal();
  fitMax = maxParam->GetVal();
  return true;
}

inline std::string formatParamTextRow(double etaMin, double etaMax, int nParams,
                                      double ptMin, double ptMax,
                                      const std::vector<double> &parameters) {
  std::ostringstream output;
  output << Form("  %7.3f %7.3f  %d  %5.0f %5.0f", etaMin, etaMax, nParams,
                 ptMin, ptMax);
  for (double parameter : parameters) {
    output << Form("  %10.6f", parameter);
  }
  output << '\n';
  return output.str();
}

inline std::string formatBinnedTextRow(double etaMin, double etaMax,
                                       double validPtMin, double validPtMax,
                                       int nParams, double formulaPtMin,
                                       double formulaPtMax, double p0,
                                       double p1, double p2) {
  return Form(
      "  %7.3f %7.3f  %8.3f %8.3f  %d  %8.3f %8.3f  %10.6f  %10.6f  %10.6f\n",
      etaMin, etaMax, validPtMin, validPtMax, nParams, formulaPtMin,
      formulaPtMax, p0, p1, p2);
}

} // namespace l2residual
