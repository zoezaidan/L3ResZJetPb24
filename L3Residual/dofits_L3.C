// Fit L3 residual corrections from one or multiple derived inputs.
//
// Workflow:
// 1. Load one or more derived ROOT files.
// 2. For each input, build or reuse an alpha->0 correction independently.
// 3. Combine the prepared pT corrections into one final L3 fit.
// 4. Write L3 and optional combined L2L3 text outputs.

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TObjArray.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "tdrStyle.C"
#include "CMS_lumi.C"

using namespace std;

struct FitWindow {
  double min;
  double max;
};

struct PreparedInput {
  TString path;
  TString label;
  TString token;
  TString jecTag;
  FitWindow fitWindow;
  TH1D* corrHist = nullptr;
  TH1D* kfsrHist = nullptr;
  double refAlphaVal = 0.0;
  double etaMin = 0.0;
  double etaMax = 0.0;
  bool haveEtaRange = false;
};

struct JecRecord {
  double etaMin;
  double etaMax;
  int nPar;
  double ptMin;
  double ptMax;
  std::vector<double> par;
};

static TH1D* drawCleaned(TH1D* h, double ptmin, double ptmax,
                         int marker, int color, bool applyWindowCut = true) {
  static int cloneCounter = 0;
  if (!h) return nullptr;

  TH1D* hc = (TH1D*)h->Clone(Form("hc_%s_%d", h->GetName(), cloneCounter++));
  hc->SetDirectory(nullptr);

  for (int i = 1; i <= hc->GetNbinsX(); ++i) {
    const double ptbinmin = hc->GetBinLowEdge(i);
    const double ptbinmax = hc->GetBinLowEdge(i + 1);
    bool keep = (ptbinmin >= ptmin && ptbinmax <= ptmax);
    if (applyWindowCut && (h->GetBinContent(i) > 1.5 || h->GetBinContent(i) < 0.5)) keep = false;

    const double val = h->GetBinContent(i);
    double err = h->GetBinError(i);
    const double errFloor = 0.002;

    if (val <= 0.0) keep = false;
    if (err <= 0.0) {
      if (applyWindowCut) keep = false;
      else err = errFloor;
    }

    if (!keep) {
      hc->SetBinContent(i, 0.0);
      hc->SetBinError(i, 0.0);
    } else {
      hc->SetBinError(i, std::sqrt(err * err + errFloor * errFloor));
    }
  }

  hc->SetMarkerStyle(marker);
  hc->SetMarkerColor(color);
  hc->SetLineColor(color);
  return hc;
}

static TGraphErrors* cleanGraph(TGraphErrors* g) {
  if (!g) return nullptr;
  for (int i = g->GetN() - 1; i >= 0; --i) {
    if (g->GetY()[i] == 0.0 && g->GetEY()[i] == 0.0) g->RemovePoint(i);
  }
  return g;
}

static std::vector<TString> splitInputFiles(const TString& in) {
  std::vector<TString> out;
  TString s = in;
  s.ReplaceAll("\n", ",");
  s.ReplaceAll("\t", ",");
  s.ReplaceAll(";", ",");
  while (s.Contains(" ")) s.ReplaceAll(" ", "");

  TObjArray* parts = s.Tokenize(",");
  if (!parts) return out;
  for (int i = 0; i < parts->GetEntriesFast(); ++i) {
    TObject* obj = parts->At(i);
    if (!obj) continue;
    TString token = obj->GetName();
    token = token.Strip(TString::kBoth);
    if (token.Length() > 0) out.push_back(token);
  }
  parts->Delete();
  delete parts;
  return out;
}

static TString defaultInputLabel(const TString& path, int indexOneBased) {
  TString base = gSystem->BaseName(path.Data());
  TString low = base;
  low.ToLower();
  if (low.Contains("photonjet") || (low.Contains("photon") && low.Contains("jet"))) return "#gamma+jet";
  if (low.Contains("zjet") || low.Contains("z+jet") || (low.Contains("z") && low.Contains("jet"))) return "Z+jet";
  if (low.Contains("dijet")) return "dijet";
  if (low.Contains("multijet")) return "multijet";
  return Form("Input %d", indexOneBased);
}

static TString inputOutputToken(const TString& label, int indexOneBased) {
  TString token = label;
  token.ToLower();
  token.ReplaceAll("#gamma", "photon");
  token.ReplaceAll("+", "plus");
  token.ReplaceAll("#", "");
  token.ReplaceAll("{", "");
  token.ReplaceAll("}", "");
  token.ReplaceAll("(", "");
  token.ReplaceAll(")", "");
  token.ReplaceAll("/", "_");
  token.ReplaceAll(" ", "_");
  while (token.Contains("__")) token.ReplaceAll("__", "_");
  token = token.Strip(TString::kBoth, '_');
  if (token.Length() == 0) token = Form("input%d", indexOneBased);
  return Form("%02d_%s", indexOneBased, token.Data());
}

static TString inputJecTag(const TString& label, int indexOneBased) {
  TString low = label;
  low.ToLower();
  if (low.Contains("#gamma") || low.Contains("photon")) return "photonjet";
  if (low.Contains("z")) return "zjet";
  return Form("input%d", indexOneBased);
}

static TH1D* makePtFrameFromHist(const TH1* href, const TString& name, const TString& yTitle,
                                 double ymin, double ymax) {
  const TString title = Form(";p_{T}^{#gamma} (GeV);%s", yTitle.Data());
  TH1D* h = nullptr;
  if (!href) {
    h = new TH1D(name, title, 100, 20, 600);
  } else {
    const int nb = href->GetXaxis()->GetNbins();
    const TArrayD* bins = href->GetXaxis()->GetXbins();
    if (bins && bins->GetSize() > 0) h = new TH1D(name, title, nb, bins->GetArray());
    else h = new TH1D(name, title, nb, href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax());
  }
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->GetXaxis()->SetTitleSize(0.038);
  h->GetXaxis()->SetLabelSize(0.032);
  h->GetYaxis()->SetTitleSize(0.038);
  h->GetYaxis()->SetLabelSize(0.032);
  h->GetYaxis()->SetTitleOffset(yTitle.Contains("#frac") ? 1.7 : 1.4);
  return h;
}

static void styleLogxAxis(TH1* h) {
  if (!h) return;
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNdivisions(510);
  h->GetXaxis()->SetNoExponent(kTRUE);
}

static TH1D* cloneHist(TFile* f, const TString& name, const TString& cloneName) {
  if (!f) return nullptr;
  TH1D* h = (TH1D*)f->Get(name);
  if (!h) return nullptr;
  TH1D* hc = (TH1D*)h->Clone(cloneName);
  hc->SetDirectory(nullptr);
  return hc;
}

static TH1D* loadPtCorrectionHist(TFile* f, int refAlphaBin) {
  if (!f) return nullptr;
  if (TH1D* h0 = cloneHist(f, "ratio_vspT_alpha0", "ratio_vspT_alpha0_clone")) return h0;
  if (TH1D* hc = cloneHist(f, "corr_vspT", "corr_vspT_clone")) return hc;
  return cloneHist(f, Form("ratio_vsphotonpt_alpha%d", refAlphaBin), Form("ratio_vsphotonpt_alpha%d_clone", refAlphaBin));
}

static std::map<int, TH1D*> loadRatioVsPtHists(TFile* f) {
  std::map<int, TH1D*> out;
  if (!f) return out;
  for (int alphaBin = 1; alphaBin <= 50; ++alphaBin) {
    TH1D* h = (TH1D*)f->Get(Form("ratio_vsphotonpt_alpha%d", alphaBin));
    if (!h) {
      if (!out.empty()) break;
      continue;
    }
    TH1D* hc = (TH1D*)h->Clone(Form("ratio_vsphotonpt_alpha%d_clone", alphaBin));
    hc->SetDirectory(nullptr);
    out[alphaBin] = hc;
  }
  return out;
}

static FitWindow defaultFitWindowFromHist(const TH1D* h) {
  if (!h) return {0.0, 0.0};
  return {h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinLowEdge(h->GetNbinsX() + 1)};
}

static FitWindow resolveFitWindow(const TString& rangeSpec, const TH1D* h) {
  FitWindow window = defaultFitWindowFromHist(h);
  if (!h) return window;

  TString token = rangeSpec;
  token.ReplaceAll(" ", "");
  if (token.Length() == 0 || !token.Contains("-")) return window;

  TString low = token(0, token.First('-'));
  TString high = token(token.First('-') + 1, token.Length());
  const double lo = low.Atof();
  const double hi = high.Atof();

  if (lo > 0.0) window.min = std::max(window.min, lo);
  if (hi > 0.0) window.max = std::min(window.max, hi);
  if (window.max < window.min) window = defaultFitWindowFromHist(h);
  return window;
}

static std::vector<double> defaultAlphaEdges() {
  return {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
}

static std::vector<double> findAlphaEdges(TFile* f, int nAlphaBins) {
  if (f) {
    TProfile3D* balance3D = (TProfile3D*)f->Get("balance3D_mc");
    if (balance3D) {
      std::vector<double> edges(balance3D->GetZaxis()->GetNbins() + 1);
      for (int i = 1; i <= balance3D->GetZaxis()->GetNbins() + 1; ++i) {
        edges[i - 1] = balance3D->GetZaxis()->GetBinLowEdge(i);
      }
      return edges;
    }
  }
  std::vector<double> edges = defaultAlphaEdges();
  if ((int)edges.size() >= nAlphaBins + 1) {
    edges.resize(nAlphaBins + 1);
    return edges;
  }
  edges.clear();
  for (int i = 0; i <= nAlphaBins; ++i) edges.push_back(0.05 * i);
  return edges;
}

static bool findEtaRange(TFile* f, double& etaMin, double& etaMax) {
  if (!f) return false;
  if (TProfile3D* balance3D = (TProfile3D*)f->Get("balance3D_mc")) {
    etaMin = balance3D->GetYaxis()->GetBinLowEdge(1);
    etaMax = balance3D->GetYaxis()->GetBinLowEdge(balance3D->GetYaxis()->GetNbins() + 1);
    return true;
  }
  if (TH2D* l3resMap = (TH2D*)f->Get("l3resMap")) {
    etaMin = l3resMap->GetYaxis()->GetBinLowEdge(1);
    etaMax = l3resMap->GetYaxis()->GetBinLowEdge(l3resMap->GetYaxis()->GetNbins() + 1);
    return true;
  }
  return false;
}

static bool isSkippableLine(const std::string& line) {
  for (char c : line) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    return (c == '#');
  }
  return true;
}

static std::string trimCopy(const std::string& text) {
  size_t begin = 0;
  while (begin < text.size() && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
  size_t end = text.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
  return text.substr(begin, end - begin);
}

static std::string replaceAllCopy(std::string text, const std::string& from, const std::string& to) {
  if (from.empty()) return text;
  size_t pos = 0;
  while ((pos = text.find(from, pos)) != std::string::npos) {
    text.replace(pos, from.size(), to);
    pos += to.size();
  }
  return text;
}

static std::string extractSimpleJecFormula(const std::string& headerLine) {
  const size_t corrPos = headerLine.find(" Correction ");
  if (corrPos == std::string::npos) return "";
  const size_t jetPtPos = headerLine.find(" JetPt ");
  if (jetPtPos == std::string::npos) return "";
  const size_t exprStart = jetPtPos + std::string(" JetPt ").size();
  if (exprStart >= corrPos) return "";
  return trimCopy(headerLine.substr(exprStart, corrPos - exprStart));
}

static std::string canonicalizeTf1Formula(const std::string& exprText) {
  std::string out = trimCopy(exprText);
  for (int index = 0; index < 20; ++index) {
    TString fromLower = Form("[p%d]", index);
    TString fromUpper = Form("[P%d]", index);
    TString to = Form("[%d]", index);
    out = replaceAllCopy(out, fromLower.Data(), to.Data());
    out = replaceAllCopy(out, fromUpper.Data(), to.Data());
  }
  return out;
}

static std::string buildSimpleJecHeader(const std::string& exprText, const std::string& correctionName) {
  return "{ 1 JetEta 1 JetPt " + trimCopy(exprText) + " Correction " + correctionName + "}";
}

static std::string extractCorrectionName(const std::string& headerLine) {
  const size_t corrPos = headerLine.find(" Correction ");
  if (corrPos == std::string::npos) return "L2Relative";
  const size_t nameStart = corrPos + std::string(" Correction ").size();
  const size_t endPos = headerLine.rfind('}');
  if (endPos == std::string::npos || endPos <= nameStart) return trimCopy(headerLine.substr(nameStart));
  return trimCopy(headerLine.substr(nameStart, endPos - nameStart));
}

static std::string buildCombinedJecHeader(const std::string& l2ExprText,
                                          const std::string& l3ExprText,
                                          const std::string& correctionName) {
  return buildSimpleJecHeader(trimCopy(l2ExprText) + "*" + trimCopy(l3ExprText), correctionName);
}

static const char* l3PtFitExpr() {
  return "1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x)";
}

static std::string l3PtFitExprShifted(int offset) {
  return Form("1./([%d]+[%d]/x+[%d]*log(x)/x+[%d]*(pow(x/[%d],[%d])-1)/(pow(x/[%d],[%d])+1)+[%d]*pow(x,-0.3051)+[%d]*x)",
              offset + 0, offset + 1, offset + 2, offset + 3,
              offset + 4, offset + 5, offset + 4, offset + 5,
              offset + 6, offset + 7);
}

static const char* l3PtFitLabel() {
  return "Run 3 L3 factor (+8 params)";
}

static int l3PtFitNPar() {
  return 8;
}

static int l3PtTxtNPar() {
  return l3PtFitNPar() + 2;
}

static void initL3PtFitParams(TF1* f, double ySeed = 1.0) {
  if (!f) return;
  const double safeSeed = (ySeed > 0.0 ? ySeed : 1.0);
  f->SetParameters(0.9722 / safeSeed, 0.7944, 2.14069, 0.10229, 10.52, 1.5550, -0.72222, -2.25e-06);
  f->SetParLimits(4, 1.0, 100.0);
  f->SetParLimits(5, 0.1, 5.0);
}

static bool readSimpleJecFile(const std::string& path, std::string& headerLine, std::vector<JecRecord>& records) {
  std::ifstream fin(path);
  if (!fin.is_open()) return false;

  records.clear();
  headerLine.clear();

  std::string line;
  while (std::getline(fin, line)) {
    if (isSkippableLine(line)) continue;
    headerLine = line;
    break;
  }
  if (headerLine.empty()) return false;

  while (std::getline(fin, line)) {
    if (isSkippableLine(line)) continue;
    std::istringstream iss(line);
    JecRecord rec;
    iss >> rec.etaMin >> rec.etaMax >> rec.nPar >> rec.ptMin >> rec.ptMax;
    if (!iss || rec.nPar < 2) continue;
    rec.par.resize(rec.nPar - 2, 0.0);
    for (int i = 0; i < rec.nPar - 2; ++i) {
      iss >> rec.par[i];
      if (!iss) {
        rec.par.clear();
        break;
      }
    }
    if (!rec.par.empty()) records.push_back(rec);
  }
  return !records.empty();
}

static std::string formatL3TxtRecord(double etaMin, double etaMax,
                                     double ptMin, double ptMax,
                                     const TF1* fit) {
  std::ostringstream out;
  out << Form("  %6.3f %6.3f  %d  %5.0f %5.0f", etaMin, etaMax, l3PtTxtNPar(), ptMin, ptMax);
  for (int index = 0; index < l3PtFitNPar(); ++index) out << Form("  %10.6f", fit->GetParameter(index));
  out << "\n";
  return out.str();
}

static std::string formatCombinedTxtRecord(const JecRecord& rec, const TF1* fit) {
  std::ostringstream out;
  out << Form("  %6.3f %6.3f  %d  %5.0f %5.0f", rec.etaMin, rec.etaMax, rec.nPar + l3PtFitNPar(), rec.ptMin, rec.ptMax);
  for (size_t index = 0; index < rec.par.size(); ++index) out << Form("  %10.6f", rec.par[index]);
  for (int index = 0; index < l3PtFitNPar(); ++index) out << Form("  %10.6f", fit->GetParameter(index));
  out << "\n";
  return out.str();
}

static void saveRatioOverlay(const std::map<int, TH1D*>& ratioVsPt,
                             const std::vector<double>& alphaEdges,
                             const TString& label,
                             const TString& token,
                             const std::string& pngFolder) {
  if (ratioVsPt.empty()) return;

  TCanvas* c = new TCanvas(Form("cRatioOverlay_%s", token.Data()), Form("cRatioOverlay_%s", token.Data()), 900, 700);
  c->SetLogx();
  c->cd();

  TH1D* href = ratioVsPt.begin()->second;
  TH1D* hFrame = makePtFrameFromHist(href, Form("hFrame_ratioOverlay_%s", token.Data()), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
  hFrame->Draw();
  styleLogxAxis(hFrame);

  const double xMin = hFrame->GetXaxis()->GetBinLowEdge(1);
  const double xMax = hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetNbinsX() + 1);
  TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 1);
  line->Draw("SAME");

  const int colors[] = {kBlack, kBlue, kRed, kGreen + 2, kMagenta + 2, kOrange + 2, kCyan + 2, kViolet + 2, kTeal + 2};
  const int nColors = sizeof(colors) / sizeof(colors[0]);

  TLegend* leg = new TLegend(0.52, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.028);

  int idx = 0;
  for (const auto& kv : ratioVsPt) {
    TH1D* h = kv.second;
    if (!h) continue;
    const int color = colors[idx % nColors];
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(0.8);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->Draw("PE1 SAME");
    const int alphaBin = kv.first;
    const double alphaCut = (alphaBin >= 1 && alphaBin < (int)alphaEdges.size()) ? alphaEdges[alphaBin] : alphaBin;
    leg->AddEntry(h, Form("#alpha < %.2f", alphaCut), "PLE");
    ++idx;
  }

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.18, 0.85, label.Data());
  leg->Draw();
  CMS_lumi(c, 0, 0);
  c->SaveAs(Form("%s/L3Res_%s_ratio_overlay_vspt.png", pngFolder.c_str(), token.Data()));

  delete tex;
  delete leg;
  delete line;
  delete hFrame;
  delete c;
}

static void saveKfsrSummaryPlots(const TString& label,
                                 const TString& token,
                                 double refAlphaVal,
                                 int refAlphaBin,
                                 TH1D* hkFSR,
                                 const std::string& alphaFolder) {
  if (!hkFSR) return;

  TCanvas* ckFSR = new TCanvas(Form("ckFSR_%s", token.Data()), Form("ckFSR_%s", token.Data()), 800, 600);
  ckFSR->SetLogx();
  ckFSR->cd();

  TH1D* hFrame = makePtFrameFromHist(hkFSR, Form("hFramekFSR_%s", token.Data()), "k_{FSR} (#alpha#rightarrow0 intercept)", 0.9, 1.1);
  hFrame->Draw();
  styleLogxAxis(hFrame);

  const double xMin = hFrame->GetXaxis()->GetBinLowEdge(1);
  const double xMax = hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetNbinsX() + 1);
  TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 1);
  line->Draw("SAME");

  hkFSR->SetMarkerStyle(kFullCircle);
  hkFSR->SetMarkerColor(kBlue);
  hkFSR->SetLineColor(kBlue);
  hkFSR->Draw("PE1 SAME");

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.18, 0.30, label.Data());
  tex->DrawLatex(0.18, 0.25, Form("Ref #alpha < %.2f (bin %d)", refAlphaVal, refAlphaBin));

  CMS_lumi(ckFSR, 0, 0);
  ckFSR->SaveAs(Form("%s/L3Res_%s_kFSR_vspT.png", alphaFolder.c_str(), token.Data()));

  delete tex;
  delete line;
  delete hFrame;
  delete ckFSR;
}

static void saveCorrPlot(const TString& label,
                         const TString& token,
                         TH1D* hRef,
                         TH1D* hCorr,
                         const std::string& alphaFolder) {
  if (!hCorr) return;

  TCanvas* c = new TCanvas(Form("cCorr_%s", token.Data()), Form("cCorr_%s", token.Data()), 900, 700);
  c->SetLogx();
  c->cd();

  TH1D* hFrame = makePtFrameFromHist(hCorr, Form("hFrameCorr_%s", token.Data()), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
  hFrame->Draw();
  styleLogxAxis(hFrame);

  const double xMin = hFrame->GetXaxis()->GetBinLowEdge(1);
  const double xMax = hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetNbinsX() + 1);
  TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 1);
  line->Draw("SAME");

  if (hRef) {
    hRef->SetMarkerStyle(kOpenCircle);
    hRef->SetMarkerColor(kBlue);
    hRef->SetLineColor(kBlue);
    hRef->Draw("PE1 SAME");
  }
  hCorr->SetMarkerStyle(kFullCircle);
  hCorr->SetMarkerColor(kRed);
  hCorr->SetLineColor(kRed);
  hCorr->Draw("PE1 SAME");

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.035);
  tex->DrawLatex(0.18, 0.30, label.Data());

  CMS_lumi(c, 0, 0);
  c->SaveAs(Form("%s/L3Res_%s_corr_alpha0_fromkFSR.png", alphaFolder.c_str(), token.Data()));

  delete tex;
  delete line;
  delete hFrame;
  delete c;
}

static TH1D* buildCorrFromRawRatios(TFile* inFile,
                                    const TString& label,
                                    const TString& token,
                                    int refAlphaBin,
                                    double fitAlphaMin,
                                    double fitAlphaMax,
                                    const std::string& pngFolder,
                                    const std::string& alphaFolder,
                                    TFile* diagFile,
                                    TH1D*& hkFSROut,
                                    double& refAlphaValOut) {
  std::map<int, TH1D*> ratioVsPt = loadRatioVsPtHists(inFile);
  if (ratioVsPt.empty() || !ratioVsPt.count(refAlphaBin)) return nullptr;

  const int nAlphaBins = (int)ratioVsPt.size();
  std::vector<double> alphaEdges = findAlphaEdges(inFile, nAlphaBins);
  refAlphaValOut = (refAlphaBin >= 1 && refAlphaBin < (int)alphaEdges.size()) ? alphaEdges[refAlphaBin] : refAlphaBin;

  saveRatioOverlay(ratioVsPt, alphaEdges, label, token, pngFolder);

  TH1D* hRef = ratioVsPt[refAlphaBin];
  std::vector<double> ptEdges(hRef->GetNbinsX() + 1);
  for (int bin = 1; bin <= hRef->GetNbinsX() + 1; ++bin) ptEdges[bin - 1] = hRef->GetBinLowEdge(bin);
  TH1D* hkFSR = new TH1D(Form("kFSR_vsPt_%s", token.Data()),
                         "k_{FSR};p_{T}^{#gamma} (GeV);k_{FSR}",
                         hRef->GetNbinsX(), ptEdges.data());
  hkFSR->SetDirectory(nullptr);

  TCanvas* cAlphaOverlay = new TCanvas(Form("cAlphaOverlay_%s", token.Data()), Form("cAlphaOverlay_%s", token.Data()), 900, 700);
  cAlphaOverlay->cd();
  TH1D* hFrameAlpha = new TH1D(Form("hFrameAlpha_%s", token.Data()), ";#alpha;Normalized Ratio", 100, 0, 0.5);
  hFrameAlpha->SetMinimum(0.97);
  hFrameAlpha->SetMaximum(1.05);
  hFrameAlpha->GetXaxis()->SetTitleSize(0.038);
  hFrameAlpha->GetXaxis()->SetLabelSize(0.032);
  hFrameAlpha->GetYaxis()->SetTitleSize(0.038);
  hFrameAlpha->GetYaxis()->SetLabelSize(0.032);
  hFrameAlpha->GetYaxis()->SetTitleOffset(1.7);
  hFrameAlpha->Draw();
  TLine* lineRefOverlay = new TLine(0, 1.0, 0.5, 1.0);
  lineRefOverlay->SetLineStyle(kDashed);
  lineRefOverlay->SetLineColor(kGray + 1);
  lineRefOverlay->Draw("SAME");

  TLegend* legOverlay = new TLegend(0.55, 0.55, 0.88, 0.88);
  legOverlay->SetBorderSize(0);
  legOverlay->SetFillStyle(0);
  legOverlay->SetTextFont(42);
  legOverlay->SetTextSize(0.025);
  std::vector<TGraph*> overlayGraphs;

  const int ptColors[] = {kBlue, kRed, kGreen + 2, kMagenta + 2, kOrange + 2, kCyan + 2, kViolet + 2, kTeal + 2, kPink + 2, kAzure + 2};
  const int nPtColors = sizeof(ptColors) / sizeof(ptColors[0]);

  for (int ptBin = 1; ptBin <= hRef->GetNbinsX(); ++ptBin) {
    const double ptLo = hRef->GetBinLowEdge(ptBin);
    const double ptHi = hRef->GetBinLowEdge(ptBin + 1);
    const double refRatio = hRef->GetBinContent(ptBin);
    if (refRatio <= 0.0) continue;

    std::vector<double> xAll;
    std::vector<double> yAll;
    std::vector<double> xFit;
    std::vector<double> yFit;

    for (const auto& kv : ratioVsPt) {
      const int alphaBin = kv.first;
      TH1D* h = kv.second;
      const double ratio = h->GetBinContent(ptBin);
      if (ratio <= 0.0) continue;
      const double alphaCenter = 0.5 * (alphaEdges[alphaBin - 1] + alphaEdges[alphaBin]);
      const double norm = ratio / refRatio;
      xAll.push_back(alphaCenter);
      yAll.push_back(norm);
      if (alphaBin != refAlphaBin && alphaCenter >= fitAlphaMin && alphaCenter <= fitAlphaMax) {
        xFit.push_back(alphaCenter);
        yFit.push_back(norm);
      }
    }

    if (xAll.size() < 2 || xFit.size() < 2) continue;

    TGraph* gAll = new TGraph(xAll.size(), xAll.data(), yAll.data());
    gAll->SetName(Form("gAlphaNorm_pt%d", ptBin));
    const int color = ptColors[(ptBin - 1) % nPtColors];
    gAll->SetMarkerStyle(kFullCircle);
    gAll->SetMarkerColor(color);
    gAll->SetLineColor(color);
    gAll->Draw("P SAME");
    overlayGraphs.push_back(gAll);

    TGraph gFit(xFit.size(), xFit.data(), yFit.data());
    TF1 fAlpha(Form("fAlpha_pt%d_col", ptBin), "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
    fAlpha.SetParameters(1.0, 0.0);
    gFit.Fit(&fAlpha, "QNR");

    const double kfsr = fAlpha.GetParameter(0);
    const double kfsrErr = fAlpha.GetParError(0);
    hkFSR->SetBinContent(ptBin, kfsr);
    hkFSR->SetBinError(ptBin, kfsrErr);

    TCanvas* cPt = new TCanvas(Form("cAlphaFit_pt%d_%s", ptBin, token.Data()), Form("cAlphaFit_pt%d_%s", ptBin, token.Data()), 800, 600);
    cPt->cd();
    TH1D* hFramePt = new TH1D(Form("hFrameAlphaFit_pt%d_%s", ptBin, token.Data()), ";#alpha;Normalized Ratio", 100, 0, 0.5);
    hFramePt->SetMinimum(0.97);
    hFramePt->SetMaximum(1.05);
    hFramePt->GetXaxis()->SetTitleSize(0.038);
    hFramePt->GetXaxis()->SetLabelSize(0.032);
    hFramePt->GetYaxis()->SetTitleSize(0.038);
    hFramePt->GetYaxis()->SetLabelSize(0.032);
    hFramePt->GetYaxis()->SetTitleOffset(1.7);
    hFramePt->Draw();
    TLine* l1 = new TLine(0, 1.0, 0.5, 1.0);
    l1->SetLineStyle(kDashed);
    l1->SetLineColor(kGray + 1);
    l1->Draw("SAME");
    gAll->Draw("P SAME");
    TF1* fLine = new TF1(Form("fAlphaLine_pt%d_%s", ptBin, token.Data()), "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
    fLine->SetParameters(fAlpha.GetParameter(0), fAlpha.GetParameter(1));
    fLine->SetLineColor(kRed);
    fLine->SetLineWidth(2);
    fLine->Draw("SAME");
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.18, 0.86, label.Data());
    tex->DrawLatex(0.18, 0.81, Form("p_{T}^{#gamma}: %.0f-%.0f GeV", ptLo, ptHi));
    tex->DrawLatex(0.18, 0.76, Form("Ref #alpha < %.2f (bin %d), excluded", refAlphaValOut, refAlphaBin));
    tex->DrawLatex(0.18, 0.71, Form("k_{FSR} = %.4f #pm %.4f", kfsr, kfsrErr));
    CMS_lumi(cPt, 0, 0);
    cPt->SaveAs(Form("%s/L3Res_%s_kFSR_alphaFit_pt%.0fto%.0f.png", alphaFolder.c_str(), token.Data(), ptLo, ptHi));

    if (diagFile) {
      diagFile->cd();
      gAll->Write(Form("gAlphaNorm_pt%d", ptBin));
      fLine->Write(Form("fAlpha_pt%d_col", ptBin));
    }

    legOverlay->AddEntry(gAll, Form("%.0f-%.0f GeV (k_{FSR}=%.3f)", ptLo, ptHi, kfsr), "P");

    delete tex;
    delete fLine;
    delete l1;
    delete hFramePt;
    delete cPt;
  }

  legOverlay->Draw();
  TLatex* texOverlay = new TLatex();
  texOverlay->SetNDC();
  texOverlay->SetTextFont(42);
  texOverlay->SetTextSize(0.030);
  texOverlay->DrawLatex(0.18, 0.85, label.Data());
  texOverlay->DrawLatex(0.18, 0.80, Form("Fit: %.2f < #alpha < %.2f", fitAlphaMin, fitAlphaMax));
  CMS_lumi(cAlphaOverlay, 0, 0);
  cAlphaOverlay->SaveAs(Form("%s/L3Res_%s_alpha_overlay.png", alphaFolder.c_str(), token.Data()));

  TH1D* hCorrAlpha0 = (TH1D*)hRef->Clone(Form("corr_vspT_%s", token.Data()));
  hCorrAlpha0->SetDirectory(nullptr);
  for (int bin = 1; bin <= hCorrAlpha0->GetNbinsX(); ++bin) {
    const double ref = hRef->GetBinContent(bin);
    const double err = hRef->GetBinError(bin);
    const double kfsr = hkFSR->GetBinContent(bin);
    if (ref > 0.0 && kfsr > 0.0) {
      const double value = ref * kfsr;
      const double rel = (ref > 0.0 && err > 0.0) ? err / ref : 0.0;
      hCorrAlpha0->SetBinContent(bin, value);
      hCorrAlpha0->SetBinError(bin, value * rel);
    } else {
      hCorrAlpha0->SetBinContent(bin, 0.0);
      hCorrAlpha0->SetBinError(bin, 0.0);
    }
  }

  saveKfsrSummaryPlots(label, token, refAlphaValOut, refAlphaBin, hkFSR, alphaFolder);
  saveCorrPlot(label, token, hRef, hCorrAlpha0, alphaFolder);

  if (diagFile) {
    diagFile->cd();
    hkFSR->Write("kFSR_vsPt");
    hCorrAlpha0->Write("corr_vspT");
  }

  hkFSROut = (TH1D*)hkFSR->Clone(Form("kFSR_vsPt_%s_clone", token.Data()));
  hkFSROut->SetDirectory(nullptr);

  delete texOverlay;
  delete legOverlay;
  delete lineRefOverlay;
  delete hFrameAlpha;
  delete cAlphaOverlay;
  delete hkFSR;
  for (auto* g : overlayGraphs) delete g;
  for (auto& kv : ratioVsPt) delete kv.second;

  return hCorrAlpha0;
}

static TH1D* reuseArchivedAlphaDiagnostics(TFile* inFile,
                                           const TString& label,
                                           const TString& token,
                                           int refAlphaBin,
                                           const std::string& alphaFolder,
                                           TFile* diagFile,
                                           TH1D*& hkFSROut,
                                           double& refAlphaValOut) {
  TH1D* hkFSR = cloneHist(inFile, "kFSR_vsPt", Form("kFSR_vsPt_%s_clone", token.Data()));
  TH1D* hAlpha0 = cloneHist(inFile, "ratio_vspT_alpha0", Form("ratio_vspT_alpha0_%s_clone", token.Data()));
  if (!hAlpha0) hAlpha0 = cloneHist(inFile, "corr_vspT", Form("corr_vspT_%s_clone", token.Data()));
  if (!hkFSR || !hAlpha0) {
    delete hkFSR;
    delete hAlpha0;
    return nullptr;
  }

  const std::vector<double> alphaEdges = findAlphaEdges(inFile, 9);
  refAlphaValOut = (refAlphaBin >= 1 && refAlphaBin < (int)alphaEdges.size()) ? alphaEdges[refAlphaBin] : refAlphaBin;

  TCanvas* cAlphaOverlay = new TCanvas(Form("cAlphaOverlayStored_%s", token.Data()), Form("cAlphaOverlayStored_%s", token.Data()), 900, 700);
  cAlphaOverlay->cd();
  TH1D* hFrameAlpha = new TH1D(Form("hFrameAlphaStored_%s", token.Data()), ";#alpha;Normalized Ratio", 100, 0, 0.5);
  hFrameAlpha->SetMinimum(0.97);
  hFrameAlpha->SetMaximum(1.05);
  hFrameAlpha->GetXaxis()->SetTitleSize(0.038);
  hFrameAlpha->GetXaxis()->SetLabelSize(0.032);
  hFrameAlpha->GetYaxis()->SetTitleSize(0.038);
  hFrameAlpha->GetYaxis()->SetLabelSize(0.032);
  hFrameAlpha->GetYaxis()->SetTitleOffset(1.7);
  hFrameAlpha->Draw();
  TLine* lineRefOverlay = new TLine(0, 1.0, 0.5, 1.0);
  lineRefOverlay->SetLineStyle(kDashed);
  lineRefOverlay->SetLineColor(kGray + 1);
  lineRefOverlay->Draw("SAME");

  TLegend* legOverlay = new TLegend(0.55, 0.55, 0.88, 0.88);
  legOverlay->SetBorderSize(0);
  legOverlay->SetFillStyle(0);
  legOverlay->SetTextFont(42);
  legOverlay->SetTextSize(0.025);

  for (int ptBin = 1; ptBin <= hkFSR->GetNbinsX(); ++ptBin) {
    TGraph* gAll = (TGraph*)inFile->Get(Form("gAlphaNorm_pt%d", ptBin));
    TF1* fAlpha = (TF1*)inFile->Get(Form("fAlpha_pt%d_col", ptBin));
    if (!gAll) continue;

    const double ptLo = hkFSR->GetBinLowEdge(ptBin);
    const double ptHi = hkFSR->GetBinLowEdge(ptBin + 1);
    gAll->SetMarkerStyle(kFullCircle);
    gAll->Draw("P SAME");

    TCanvas* cPt = new TCanvas(Form("cStoredAlphaFit_pt%d_%s", ptBin, token.Data()), Form("cStoredAlphaFit_pt%d_%s", ptBin, token.Data()), 800, 600);
    cPt->cd();
    TH1D* hFramePt = new TH1D(Form("hFrameStoredAlphaFit_pt%d_%s", ptBin, token.Data()), ";#alpha;Normalized Ratio", 100, 0, 0.5);
    hFramePt->SetMinimum(0.97);
    hFramePt->SetMaximum(1.05);
    hFramePt->GetXaxis()->SetTitleSize(0.038);
    hFramePt->GetXaxis()->SetLabelSize(0.032);
    hFramePt->GetYaxis()->SetTitleSize(0.038);
    hFramePt->GetYaxis()->SetLabelSize(0.032);
    hFramePt->GetYaxis()->SetTitleOffset(1.7);
    hFramePt->Draw();
    TLine* l1 = new TLine(0, 1.0, 0.5, 1.0);
    l1->SetLineStyle(kDashed);
    l1->SetLineColor(kGray + 1);
    l1->Draw("SAME");
    gAll->Draw("P SAME");
    if (fAlpha) {
      fAlpha->SetLineColor(kRed);
      fAlpha->SetLineWidth(2);
      fAlpha->Draw("SAME");
    }
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    tex->DrawLatex(0.18, 0.86, label.Data());
    tex->DrawLatex(0.18, 0.81, Form("p_{T}^{#gamma}: %.0f-%.0f GeV", ptLo, ptHi));
    tex->DrawLatex(0.18, 0.76, "Stored alpha fit from input file");
    tex->DrawLatex(0.18, 0.71, Form("k_{FSR} = %.4f #pm %.4f", hkFSR->GetBinContent(ptBin), hkFSR->GetBinError(ptBin)));
    CMS_lumi(cPt, 0, 0);
    cPt->SaveAs(Form("%s/L3Res_%s_kFSR_alphaFit_pt%.0fto%.0f.png", alphaFolder.c_str(), token.Data(), ptLo, ptHi));

    if (diagFile) {
      diagFile->cd();
      gAll->Write(Form("gAlphaNorm_pt%d", ptBin));
      if (fAlpha) fAlpha->Write(Form("fAlpha_pt%d_col", ptBin));
    }

    legOverlay->AddEntry(gAll, Form("%.0f-%.0f GeV (k_{FSR}=%.3f)", ptLo, ptHi, hkFSR->GetBinContent(ptBin)), "P");

    delete tex;
    delete l1;
    delete hFramePt;
    delete cPt;
  }

  legOverlay->Draw();
  TLatex* texOverlay = new TLatex();
  texOverlay->SetNDC();
  texOverlay->SetTextFont(42);
  texOverlay->SetTextSize(0.030);
  texOverlay->DrawLatex(0.18, 0.85, label.Data());
  texOverlay->DrawLatex(0.18, 0.80, "Stored alpha-fit diagnostics from input file");
  CMS_lumi(cAlphaOverlay, 0, 0);
  cAlphaOverlay->SaveAs(Form("%s/L3Res_%s_alpha_overlay.png", alphaFolder.c_str(), token.Data()));

  saveKfsrSummaryPlots(label, token, refAlphaValOut, refAlphaBin, hkFSR, alphaFolder);
  saveCorrPlot(label, token, nullptr, hAlpha0, alphaFolder);

  if (diagFile) {
    diagFile->cd();
    hkFSR->Write("kFSR_vsPt");
    hAlpha0->Write("corr_vspT");
  }

  hkFSROut = hkFSR;

  delete texOverlay;
  delete legOverlay;
  delete lineRefOverlay;
  delete hFrameAlpha;
  delete cAlphaOverlay;

  return hAlpha0;
}

static PreparedInput prepareInput(const TString& inputPath,
                                  const TString& inputLabel,
                                  int inputIndex,
                                  const std::string& outfolder,
                                  int refAlphaBin,
                                  double fitAlphaMin,
                                  double fitAlphaMax,
                                  bool saveAlphaExtrap) {
  PreparedInput prepared;
  prepared.path = inputPath;
  prepared.label = inputLabel;
  prepared.token = inputOutputToken(inputLabel, inputIndex + 1);
  prepared.jecTag = inputJecTag(inputLabel, inputIndex + 1);

  const std::string inputOutfolder = outfolder + "/inputs/" + std::string(prepared.token.Data());
  const std::string pngFolder = inputOutfolder + "/pdf";
  const std::string alphaFolder = inputOutfolder + "/alpha_extrap";
  const std::string txtFolder = inputOutfolder + "/textfiles";
  gSystem->mkdir(inputOutfolder.c_str(), kTRUE);
  gSystem->mkdir(pngFolder.c_str(), kTRUE);
  gSystem->mkdir(txtFolder.c_str(), kTRUE);
  if (saveAlphaExtrap) gSystem->mkdir(alphaFolder.c_str(), kTRUE);

  TFile* inFile = TFile::Open(inputPath, "READ");
  if (!inFile || inFile->IsZombie()) {
    cout << "WARNING: cannot open input " << inputPath << endl;
    return prepared;
  }

  findEtaRange(inFile, prepared.etaMin, prepared.etaMax);
  prepared.haveEtaRange = (prepared.etaMax > prepared.etaMin);

  TFile* diagFile = nullptr;
  if (saveAlphaExtrap) {
    diagFile = new TFile(Form("%s/%s_alphaDiagnostics.root", inputOutfolder.c_str(), prepared.token.Data()), "RECREATE");
  }

  if (saveAlphaExtrap) {
    prepared.corrHist = buildCorrFromRawRatios(inFile, prepared.label, prepared.token,
                                               refAlphaBin, fitAlphaMin, fitAlphaMax,
                                               pngFolder, alphaFolder, diagFile,
                                               prepared.kfsrHist, prepared.refAlphaVal);
    if (!prepared.corrHist) {
      prepared.corrHist = reuseArchivedAlphaDiagnostics(inFile, prepared.label, prepared.token,
                                                        refAlphaBin, alphaFolder, diagFile,
                                                        prepared.kfsrHist, prepared.refAlphaVal);
    }
  }

  if (!prepared.corrHist) {
    prepared.corrHist = loadPtCorrectionHist(inFile, refAlphaBin);
    if (prepared.corrHist && prepared.refAlphaVal <= 0.0) {
      const std::vector<double> edges = findAlphaEdges(inFile, refAlphaBin + 1);
      if (refAlphaBin >= 1 && refAlphaBin < (int)edges.size()) prepared.refAlphaVal = edges[refAlphaBin];
    }
  }

  if (diagFile) {
    diagFile->Close();
    delete diagFile;
  }

  inFile->Close();
  delete inFile;
  return prepared;
}

void dofits_L3(TString inFileL3Derived = "L3_derived.root",
               TString inputPtRangesCSV = "",
               string outfilename = "L3Res_photonjet",
               string runLabel = "2024ppRef",
               string lumiLabel = "pp 480.4 pb^{-1}",
               bool saveAlphaExtrap = true,
               int refAlphaBin = 5,
               double fitAlphaMin = 0.0,
               double fitAlphaMax = 0.4,
               string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
               string outBaseDir = "L3Residual",
               TString inputLabelsCSV = "") {

  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  writeExtraText = true;
  extraText = "Preliminary";
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.c_str());

  const std::string outfolder = outBaseDir + "/" + outfilename;
  const std::string pngFolder = outfolder + "/pdf";
  const std::string txtFolder = outfolder + "/textfiles";
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  gSystem->mkdir(pngFolder.c_str(), kTRUE);
  gSystem->mkdir(txtFolder.c_str(), kTRUE);

  std::vector<TString> inputFiles = splitInputFiles(inFileL3Derived);
  if (inputFiles.empty()) inputFiles.push_back(inFileL3Derived);

  std::vector<TString> inputLabels = splitInputFiles(inputLabelsCSV);
  if (inputLabels.empty()) {
    if (inputFiles.size() == 1) {
      inputLabels.push_back("#gamma+jet");
    } else if (inputFiles.size() == 2) {
      inputLabels.push_back("#gamma+jet");
      inputLabels.push_back("Z+jet");
    }
  }
  while (inputLabels.size() < inputFiles.size()) {
    inputLabels.push_back(defaultInputLabel(inputFiles[inputLabels.size()], (int)inputLabels.size() + 1));
  }

  std::vector<TString> fitRangeSpecs = splitInputFiles(inputPtRangesCSV);

  cout << "Using L3 fit model: " << l3PtFitLabel() << " (" << l3PtFitExpr() << ")" << endl;
  cout << "============================================" << endl;
  cout << "L3 Residual Fitting" << endl;
  cout << "Input file(s): " << inFileL3Derived << endl;
  if (inputPtRangesCSV.Length() > 0) cout << "Requested pT fit windows: " << inputPtRangesCSV << endl;
  else cout << "Requested pT fit windows: full histogram range per input" << endl;
  cout << "Reference alpha bin: " << refAlphaBin << endl;
  cout << "Alpha fit range: " << fitAlphaMin << " - " << fitAlphaMax << endl;
  cout << "============================================" << endl;

  if (inputFiles.size() > 1 && saveAlphaExtrap) {
    cout << "NOTE: saveAlphaExtrap=true with multiple inputs." << endl;
    cout << "      Each input will run its own alpha-extrapolation diagnostics in a separate subdirectory." << endl;
    cout << "      Only the final pT fit is combined across inputs." << endl;
  }

  std::vector<PreparedInput> preparedInputs;
  preparedInputs.reserve(inputFiles.size());
  for (size_t i = 0; i < inputFiles.size(); ++i) {
    PreparedInput prepared = prepareInput(inputFiles[i], inputLabels[i], (int)i, outfolder,
                                          refAlphaBin, fitAlphaMin, fitAlphaMax, saveAlphaExtrap);
    if (!prepared.corrHist) {
      cout << "WARNING: skipping input with no valid pT correction histogram: " << inputFiles[i] << endl;
      continue;
    }
    const TString fitRangeSpec = (i < fitRangeSpecs.size() ? fitRangeSpecs[i] : TString(""));
    prepared.fitWindow = resolveFitWindow(fitRangeSpec, prepared.corrHist);
    preparedInputs.push_back(prepared);
  }

  if (preparedInputs.empty()) {
    cout << "ERROR: no valid inputs available for the final pT fit." << endl;
    return;
  }

  double fitRangeMinOut = preparedInputs.front().fitWindow.min;
  double fitRangeMaxOut = preparedInputs.front().fitWindow.max;
  double etaMinOut = preparedInputs.front().etaMin;
  double etaMaxOut = preparedInputs.front().etaMax;
  bool haveEtaRange = preparedInputs.front().haveEtaRange;
  double refAlphaValOut = preparedInputs.front().refAlphaVal;

  for (size_t i = 1; i < preparedInputs.size(); ++i) {
    fitRangeMinOut = std::min(fitRangeMinOut, preparedInputs[i].fitWindow.min);
    fitRangeMaxOut = std::max(fitRangeMaxOut, preparedInputs[i].fitWindow.max);
    if (preparedInputs[i].haveEtaRange) {
      if (!haveEtaRange) {
        etaMinOut = preparedInputs[i].etaMin;
        etaMaxOut = preparedInputs[i].etaMax;
        haveEtaRange = true;
      } else {
        etaMinOut = std::min(etaMinOut, preparedInputs[i].etaMin);
        etaMaxOut = std::max(etaMaxOut, preparedInputs[i].etaMax);
      }
    }
  }

  TFile* outfile = new TFile(Form("%s/%s_fit.root", outfolder.c_str(), outfilename.c_str()), "RECREATE");

  TCanvas* cFinal = new TCanvas("cFinal_combined", "L3Residual final pT fit", 900, 700);
  cFinal->SetLogx();
  cFinal->SetLeftMargin(0.14);
  cFinal->SetBottomMargin(0.12);
  cFinal->SetRightMargin(0.05);
  cFinal->SetTopMargin(0.06);
  cFinal->cd();

  double xMin = preparedInputs.front().corrHist->GetXaxis()->GetBinLowEdge(1);
  double xMax = preparedInputs.front().corrHist->GetXaxis()->GetBinLowEdge(preparedInputs.front().corrHist->GetNbinsX() + 1);
  for (size_t i = 1; i < preparedInputs.size(); ++i) {
    xMin = std::min(xMin, preparedInputs[i].corrHist->GetXaxis()->GetBinLowEdge(1));
    xMax = std::max(xMax, preparedInputs[i].corrHist->GetXaxis()->GetBinLowEdge(preparedInputs[i].corrHist->GetNbinsX() + 1));
  }

  const double yMin = 0.7;
  const double yMax = 1.5;
  TMultiGraph* mg = new TMultiGraph("mg_combined_pt", "mg_combined_pt");
  TLegend* leg = new TLegend(0.50, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.028);

  const int colors[] = {kBlue, kRed, kGreen + 2, kMagenta + 2, kOrange + 2, kCyan + 2, kViolet + 2, kTeal + 2};
  const int nColors = sizeof(colors) / sizeof(colors[0]);
  const int markers[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kFullStar, kFullCircle};
  const int nMarkers = sizeof(markers) / sizeof(markers[0]);

  std::vector<TH1D*> displayHists;
  std::vector<TH1D*> fitHists;
  std::vector<TGraphErrors*> displayGraphs;
  std::vector<TGraphErrors*> fitGraphs;

  for (size_t i = 0; i < preparedInputs.size(); ++i) {
    const int color = colors[i % nColors];
    const int marker = markers[i % nMarkers];
    TH1D* hAll = drawCleaned(preparedInputs[i].corrHist, xMin, xMax, marker, color, false);
    TH1D* hFit = drawCleaned(preparedInputs[i].corrHist, preparedInputs[i].fitWindow.min, preparedInputs[i].fitWindow.max, marker, color, true);
    TGraphErrors* gAll = (hAll ? cleanGraph(new TGraphErrors(hAll)) : nullptr);
    TGraphErrors* gFit = (hFit ? cleanGraph(new TGraphErrors(hFit)) : nullptr);

    displayHists.push_back(hAll);
    fitHists.push_back(hFit);
    displayGraphs.push_back(gAll);
    fitGraphs.push_back(gFit);

    if (gAll && gAll->GetN() > 0) {
      gAll->SetMarkerStyle(marker);
      gAll->SetMarkerColor(color);
      gAll->SetLineColor(color);
      gAll->SetMarkerSize(1.05);
      mg->Add(gAll, "P");
      leg->AddEntry(gAll, Form("%s (%.0f-%.0f)", preparedInputs[i].label.Data(),
                               preparedInputs[i].fitWindow.min, preparedInputs[i].fitWindow.max), "P");
    }

    cout << Form("INFO: input '%s' -> plotted points = %d, in-range points = %d\n",
                 preparedInputs[i].label.Data(), (gAll ? gAll->GetN() : 0), (gFit ? gFit->GetN() : 0));
  }

  mg->SetTitle(";p_{T}^{#gamma} (GeV);#frac{R_{MC}}{R_{Data}}");
  mg->Draw("AP");
  mg->GetXaxis()->SetLimits(xMin, xMax);
  mg->SetMinimum(yMin);
  mg->SetMaximum(yMax);
  if (mg->GetHistogram()) {
    mg->GetHistogram()->SetMinimum(yMin);
    mg->GetHistogram()->SetMaximum(yMax);
    styleLogxAxis(mg->GetHistogram());
  }
  mg->GetXaxis()->SetTitleSize(0.038);
  mg->GetXaxis()->SetLabelSize(0.032);
  mg->GetYaxis()->SetTitleSize(0.038);
  mg->GetYaxis()->SetLabelSize(0.032);
  mg->GetYaxis()->SetTitleOffset(1.7);
  mg->GetYaxis()->SetRangeUser(yMin, yMax);
  mg->Draw("P SAME");

  TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kGray + 1);
  line->Draw("SAME");

  TGraphErrors* gFitTotal = new TGraphErrors();
  int nTotal = 0;
  for (auto* g : fitGraphs) {
    if (!g) continue;
    for (int ip = 0; ip < g->GetN(); ++ip) {
      gFitTotal->SetPoint(nTotal, g->GetX()[ip], g->GetY()[ip]);
      gFitTotal->SetPointError(nTotal, (g->GetEX() ? g->GetEX()[ip] : 0.0), (g->GetEY() ? g->GetEY()[ip] : 0.0));
      ++nTotal;
    }
  }

  if (gFitTotal->GetN() < 2) {
    cout << "ERROR: not enough in-range points for the final pT fit (N=" << gFitTotal->GetN() << ")" << endl;
    outfile->Close();
    delete outfile;
    delete gFitTotal;
    return;
  }

  double ySeed = 0.0;
  for (int ip = 0; ip < gFitTotal->GetN(); ++ip) ySeed += gFitTotal->GetY()[ip];
  ySeed /= gFitTotal->GetN();

  TF1* fref = new TF1("fit_combined_ref", l3PtFitExpr(), fitRangeMinOut, fitRangeMaxOut);
  initL3PtFitParams(fref, ySeed);
  gFitTotal->Fit(fref, "QRN");
  fref->SetLineColor(kRed + 1);
  fref->SetLineWidth(2);
  fref->Draw("SAME");

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.032);
  if (preparedInputs.size() == 1) tex->DrawLatex(0.18, 0.85, Form("Input: %s", preparedInputs.front().label.Data()));
  else tex->DrawLatex(0.18, 0.85, Form("Inputs: %zu", preparedInputs.size()));
  tex->DrawLatex(0.18, 0.80, Form("Ref #alpha < %.2f (bin %d)", refAlphaValOut, refAlphaBin));
  tex->DrawLatex(0.18, 0.75, Form("Fit model: %s", l3PtFitLabel()));
  tex->DrawLatex(0.18, 0.70, Form("Fit range: %.0f-%.0f GeV", fitRangeMinOut, fitRangeMaxOut));
  tex->DrawLatex(0.18, 0.65, Form("#chi^{2}/ndf = %.1f/%d", fref->GetChisquare(), fref->GetNDF()));

  leg->Draw();
  CMS_lumi(cFinal, 0, 0);
  cFinal->SaveAs(Form("%s/L3Res_%s_combined_final.png", pngFolder.c_str(), runLabel.c_str()));

  outfile->cd();
  mg->Write();
  fref->Write();
  for (size_t i = 0; i < preparedInputs.size(); ++i) {
    if (preparedInputs[i].corrHist) preparedInputs[i].corrHist->Write(Form("corr_input_%zu", i));
    if (preparedInputs[i].kfsrHist) preparedInputs[i].kfsrHist->Write(Form("kFSR_input_%zu", i));
  }

  std::string l2Header;
  std::vector<JecRecord> l2recs;
  const bool haveL2 = (!l2ResidualFile.empty() && readSimpleJecFile(l2ResidualFile, l2Header, l2recs));
  const std::string l3Header = buildSimpleJecHeader(l3PtFitExpr(), "L2Relative");
  const std::string localTxt = Form("%s/%s.txt", txtFolder.c_str(), outfilename.c_str());
  std::ofstream outL3(localTxt);
  outL3 << l3Header << "\n";
  if (haveL2) {
    for (const auto& rec : l2recs) outL3 << formatL3TxtRecord(rec.etaMin, rec.etaMax, fitRangeMinOut, fitRangeMaxOut, fref);
  } else {
    const double etaMin = haveEtaRange ? etaMinOut : -5.191;
    const double etaMax = haveEtaRange ? etaMaxOut : 5.191;
    outL3 << formatL3TxtRecord(etaMin, etaMax, fitRangeMinOut, fitRangeMaxOut, fref);
  }
  outL3.close();

  const std::string outTag = (preparedInputs.size() > 1 ? "combined" : std::string(preparedInputs.front().jecTag.Data()));
  const std::string l3TxtJec = Form("%s/L3Residuals_%s_%s_AK4PF.txt", txtFolder.c_str(), runLabel.c_str(), outTag.c_str());
  std::ofstream outL3Jec(l3TxtJec);
  outL3Jec << l3Header << "\n";
  if (haveL2) {
    for (const auto& rec : l2recs) outL3Jec << formatL3TxtRecord(rec.etaMin, rec.etaMax, fitRangeMinOut, fitRangeMaxOut, fref);
  } else {
    const double etaMin = haveEtaRange ? etaMinOut : -5.191;
    const double etaMax = haveEtaRange ? etaMaxOut : 5.191;
    outL3Jec << formatL3TxtRecord(etaMin, etaMax, fitRangeMinOut, fitRangeMaxOut, fref);
  }
  outL3Jec.close();
  cout << "\nWrote L3Residual JEC text: " << l3TxtJec << endl;

  if (haveL2) {
    const std::string l2ExprTextRaw = extractSimpleJecFormula(l2Header);
    const std::string l2ExprText = (!l2ExprTextRaw.empty() ? l2ExprTextRaw : "1./([0]+[1]*log10(0.01*x)+[2]/(x/10.0))");
    const std::string l2ExprCanonical = canonicalizeTf1Formula(l2ExprText);
    const std::string correctionName = extractCorrectionName(l2Header);
    const std::string l2l3Header = buildCombinedJecHeader(l2ExprCanonical, l3PtFitExprShifted((int)l2recs.front().par.size()), correctionName);
    const std::string l2l3Txt = Form("%s/L2L3Residuals_%s_%s_AK4PF.txt", txtFolder.c_str(), runLabel.c_str(), outTag.c_str());
    std::ofstream outL2L3(l2l3Txt);
    outL2L3 << l2l3Header << "\n";
    for (const auto& rec : l2recs) outL2L3 << formatCombinedTxtRecord(rec, fref);
    outL2L3.close();
    cout << "Wrote combined L2L3Residual JEC text: " << l2l3Txt << endl;
  }

  outfile->Close();
  delete outfile;

  cout << "\n============================================" << endl;
  cout << "Output written to: " << outfolder << "/" << outfilename << "_fit.root" << endl;
  cout << "Plots saved in: " << pngFolder << "/" << endl;
  cout << "Final L3Residual text file: " << txtFolder << "/" << outfilename << ".txt" << endl;
  cout << "============================================" << endl;
}