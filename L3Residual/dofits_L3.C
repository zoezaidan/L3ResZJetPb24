// Fit L3 Residual corrections from photon+jet balancing
// Production-ready version with TDR style output
// Follows template from L3Res.C
//
// Usage:
//   root -l -b -q 'dofits_L3.C("input.root", 30, 100, "output", false)'
//
// Author: Based on L3Res.C template

#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>

#include "tdrStyle.C"
#include "CMS_lumi.C"

using namespace std;

TLegend *_leg = nullptr;
string _run = "2024ppRef";
bool fitG = true;
bool fitD = false;
bool doClosure = false;

TH1D* drawCleaned(TH1D* h, string data, double ptmin, double ptmax,
                  int marker, int color, bool applyWindowCut = true) {
  if (!h) return nullptr;
  TH1D* hc = (TH1D*)h->Clone(Form("hc_%s_%s", h->GetName(), _run.c_str()));
  for (int i = 1; i <= hc->GetNbinsX(); ++i) {
    double ptbinmin = hc->GetBinLowEdge(i);
    double ptbinmax = hc->GetBinLowEdge(i+1);
    bool keep = false;
    if (data == "G" && ptbinmin >= ptmin && ptbinmax <= ptmax) keep = true;
    if (data == "D" && ptbinmin >= ptmin && ptbinmax <= ptmax) keep = true;
    if (applyWindowCut && (h->GetBinContent(i) > 1.5 || h->GetBinContent(i) < 0.5)) keep = false;

    const double val = h->GetBinContent(i);
    double err = h->GetBinError(i);
    const double errmin = 0.002; // small floor error used for display and fit stability

    // If bin content is zero or negative, always drop
    if (val <= 0) keep = false;

    // If bin error is zero or negative:
    // - for fit/cleaned mode (applyWindowCut==true) drop the bin (insufficient info)
    // - for display-only mode (applyWindowCut==false) keep the bin but assign a small error
    if (err <= 0) {
      if (applyWindowCut) keep = false;
      else err = errmin;
    }

    if (!keep) {
      hc->SetBinContent(i, 0.);
      hc->SetBinError(i, 0.);
    } else {
      hc->SetBinError(i, sqrt(pow(err, 2) + pow(errmin, 2)));
    }
  }
  hc->SetMarkerStyle(marker);
  hc->SetMarkerColor(color);
  hc->SetLineColor(color);
  return hc;
}

TGraphErrors* cleanGraph(TGraphErrors* g) {
  if (!g) return nullptr;
  for (int i = g->GetN() - 1; i >= 0; --i) {
    if (g->GetY()[i] == 0 && g->GetEY()[i] == 0) g->RemovePoint(i);
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
    auto* obj = parts->At(i);
    if (!obj) continue;
    TString token = obj->GetName();
    token = token.Strip(TString::kBoth);
    if (token.Length() > 0) out.push_back(token);
  }
  parts->Delete();
  delete parts;
  return out;
}

static TH1D* loadPtCorrectionHist(TFile* f, int refAlphaBin) {
  if (!f) return nullptr;
  // Prefer an already-built alpha->0 correction if present
  if (auto* h0 = (TH1D*)f->Get("ratio_vspT_alpha0")) return h0;
  if (auto* hc = (TH1D*)f->Get("corr_vspT")) return hc;
  // Fall back to the nominal-alpha ratio
  return (TH1D*)f->Get(Form("ratio_vsphotonpt_alpha%d", refAlphaBin));
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

static TH1D* makePtFrameFromHist(const TH1* href, const TString& name, const TString& yTitle, double ymin, double ymax) {
  const TString title = Form(";p_{T}^{#gamma} (GeV);%s", yTitle.Data());
  if (!href) {
    TH1D* h = new TH1D(name, title, 100, 20, 600);
    h->SetMinimum(ymin);
    h->SetMaximum(ymax);
    // Slightly reduced sizes for axis titles/labels to improve fit of long fraction titles
    h->GetXaxis()->SetTitleSize(0.038);
    h->GetXaxis()->SetLabelSize(0.032);
    h->GetYaxis()->SetTitleSize(0.038);
    h->GetYaxis()->SetLabelSize(0.032);
    // fraction-like labels need more room
    if (yTitle.Contains("#frac")) h->GetYaxis()->SetTitleOffset(1.7);
    else h->GetYaxis()->SetTitleOffset(1.4);
    return h;
  }
  const int nb = href->GetXaxis()->GetNbins();
  const TArrayD* bins = href->GetXaxis()->GetXbins();
  TH1D* h = nullptr;
  if (bins && bins->GetSize() > 0) {
    h = new TH1D(name, title, nb, bins->GetArray());
  } else {
    h = new TH1D(name, title, nb, href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax());
  }
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  // axis sizing
  // Slightly reduced sizes for axis titles/labels to improve fit of long fraction titles
  h->GetXaxis()->SetTitleSize(0.038);
  h->GetXaxis()->SetLabelSize(0.032);
  h->GetYaxis()->SetTitleSize(0.038);
  h->GetYaxis()->SetLabelSize(0.032);
  if (yTitle.Contains("#frac")) h->GetYaxis()->SetTitleOffset(1.7);
  else h->GetYaxis()->SetTitleOffset(1.4);
  return h;
} 

static void styleLogxAxis(TH1* h) {
  if (!h) return;
  h->GetXaxis()->SetMoreLogLabels(kTRUE);
  h->GetXaxis()->SetNdivisions(510);
  h->GetXaxis()->SetNoExponent(kTRUE);
}

static bool isSkippableLine(const std::string& line) {
  for (char c : line) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    return (c == '#');
  }
  return true;
}

struct JecRecord {
  double etaMin;
  double etaMax;
  int nPar;          // includes 2 range parameters
  double ptMin;
  double ptMax;
  std::vector<double> par; // size = nPar - 2
};

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

// Centralized fit function for the FINAL L3 correction vs pT.
// kFSR is extracted from the linear fit vs alpha; only the final correction is fit vs pT.
static const char* l3PtFitExpr() {
  // Current choice: constant + log(pT).
  // Add higher-order terms here in the future if needed.
  return "[0]+[1]*log10(0.01*x)";
}

static TF1* makeL3PtFitFunc(const char* name, double xmin, double xmax) {
  TF1* f = new TF1(name, l3PtFitExpr(), xmin, xmax);
  f->SetParameters(1.0, 0.0);
  return f;
}

void dofits_L3(TString inFileL3Derived = "L3_derived.root",
               double ptminG = 30.,
               double ptmaxG = 100.,
               string outfilename = "L3Res_photonjet",
               bool closure = false,
               string runLabel = "2024ppRef",
               string lumiLabel = "pp 480.4 pb^{-1}",
               bool plotRawResponses = true,
               bool saveAlphaExtrap = false,
               int refAlphaBin = 5,
               double fitAlphaMin = 0.0,
               double fitAlphaMax = 0.4,
               bool applyKFSRToPtFit = true,
               bool doEtaBinnedAlphaFits = false,
               bool useSingleEtaBin = true,
               int etaBinForL3 = 1,
               bool plotEtaMaps = false,
               bool plotBalanceDistOverlay = true,
               TString mcRawFileForDist = "",
               TString dataRawFileForDist = "",
               bool writeL2L3 = true,
               string l2ResidualFile = "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt",
               string outBaseDir = "L3Residual",
               string jecOutDir = "",
               TString inputLabelsCSV = "",
               bool doCombinedPtFit = false,
               TString inputPtRangesCSV = "",
               bool plotPerAlphaPtFits = false) {

  doClosure = closure;
  _run = runLabel;

  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  writeExtraText = true;
  extraText = "Preliminary";
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.c_str());

  string outfolder = outBaseDir + "/" + outfilename;
  string pngFolder = outfolder + "/pdf"; // keep folder name for compatibility
  string txtFolder = outfolder + "/textfiles";
  string rawFolder = outfolder + "/raw";
  string alphaFolder = outfolder + "/alpha_extrap";
  gSystem->mkdir(outfolder.c_str(), kTRUE);
  gSystem->mkdir(pngFolder.c_str(), kTRUE);
  gSystem->mkdir(txtFolder.c_str(), kTRUE);
  if (jecOutDir.empty()) jecOutDir = outBaseDir + "/jecfiles";
  gSystem->mkdir(jecOutDir.c_str(), kTRUE);
  if (writeL2L3 && l2ResidualFile.empty()) {
    cout << "WARNING: writeL2L3=true but l2ResidualFile is empty; skipping L2L3 output." << endl;
    writeL2L3 = false;
  }
  if (plotRawResponses) gSystem->mkdir(rawFolder.c_str(), kTRUE);
  if (saveAlphaExtrap) gSystem->mkdir(alphaFolder.c_str(), kTRUE);

  // Optional: overlay balance distributions (requires raw analysis ROOT files)
  if (plotBalanceDistOverlay) {
    if (mcRawFileForDist.Length() == 0 || dataRawFileForDist.Length() == 0) {
      cout << "WARNING: plotBalanceDistOverlay=true but mcRawFileForDist/dataRawFileForDist not provided; skipping." << endl;
      plotBalanceDistOverlay = false;
    }
  }

  std::vector<TString> inputFiles = splitInputFiles(inFileL3Derived);
  if (inputFiles.empty()) inputFiles.push_back(inFileL3Derived);
  const bool isMultiInput = (inputFiles.size() > 1);

  if (isMultiInput && saveAlphaExtrap) {
    cout << "NOTE: saveAlphaExtrap=true with multiple inputs." << endl;
    cout << "      This macro does NOT recompute kFSR separately per input." << endl;
    cout << "      For proper combined fits, provide per-input pT corrections as 'corr_vspT' or 'ratio_vspT_alpha0' in each derived ROOT file." << endl;
  }

  std::vector<TString> inputLabels = splitInputFiles(inputLabelsCSV);
  // If labels are not provided or mismatched, fall back to file basenames
  if (inputLabels.size() != inputFiles.size()) {
    inputLabels.clear();
    for (size_t i = 0; i < inputFiles.size(); ++i) inputLabels.push_back(defaultInputLabel(inputFiles[i], (int)i + 1));
  }

  // Optional per-input pT ranges, format like: "60-300,300-1000" (same length as inputs)
  std::vector<TString> inputPtRanges = splitInputFiles(inputPtRangesCSV);

  if (isMultiInput && !doCombinedPtFit) {
    cout << "WARNING: multiple inputs provided but doCombinedPtFit=false; will only use the first input." << endl;
  }

  // Use the first file for determining binning and for optional per-alpha diagnostics
  TFile* inFile = TFile::Open(inputFiles.front(), "READ");
  if (!inFile || inFile->IsZombie()) {
    cout << "ERROR: Cannot open input file: " << inputFiles.front() << endl;
    return;
  }

  cout << "============================================" << endl;
  cout << "L3 Residual Fitting" << endl;
  cout << "Input file(s): " << inFileL3Derived << endl;
  cout << "Photon+jet pT range: " << ptminG << " - " << ptmaxG << " GeV" << endl;
  cout << "Reference alpha bin: " << refAlphaBin << endl;
  cout << "Alpha fit range: " << fitAlphaMin << " - " << fitAlphaMax << endl;
  cout << "Closure test: " << (doClosure ? "YES" : "NO") << endl;
  cout << "============================================" << endl;

  TProfile3D* balance3D_mc = (TProfile3D*)inFile->Get("balance3D_mc");
  TProfile3D* balance3D_data = (TProfile3D*)inFile->Get("balance3D_data");

  TH3D* counts3D_mc = (TH3D*)inFile->Get("counts3D_mc");
  TH3D* counts3D_data = (TH3D*)inFile->Get("counts3D_data");
  if (!counts3D_mc) counts3D_mc = (TH3D*)inFile->Get("counts3D_mc_wide");
  if (!counts3D_mc) counts3D_mc = (TH3D*)inFile->Get("counts3D_mc_narrow");
  if (!counts3D_data) counts3D_data = (TH3D*)inFile->Get("counts3D_data_wide");
  if (!counts3D_data) counts3D_data = (TH3D*)inFile->Get("counts3D_data_narrow");

  if (!balance3D_mc || !balance3D_data) {
    cout << "ERROR: Cannot find balance3D profiles in input file" << endl;
    return;
  }

  int nptbins = balance3D_mc->GetXaxis()->GetNbins();
  int netabins = balance3D_mc->GetYaxis()->GetNbins();
  int nalphabins = balance3D_mc->GetZaxis()->GetNbins();

  cout << "Detected binning: " << nptbins << " pT bins, "
       << netabins << " eta bins, " << nalphabins << " alpha bins" << endl;
  cout << "Alpha bin edges: ";
  for (int a = 1; a <= nalphabins + 1; ++a) cout << balance3D_mc->GetZaxis()->GetBinLowEdge(a) << " ";
  cout << endl;

  if (!counts3D_mc) cout << "WARNING: counts3D_mc not found" << endl;
  if (!counts3D_data) cout << "WARNING: counts3D_data not found" << endl;

  TFile* outfile = new TFile(Form("%s/%s_fit.root", outfolder.c_str(), outfilename.c_str()), "RECREATE");

  // Optional: build alpha-extrapolation (per pT, eta) from the 3D profiles
  // Storage for kFSR extrapolation results (per pT bin, collapsed over eta)
  map<int, double> kFSR_vsPt;       // kFSR = p0 from linear fit at alpha->0
  map<int, double> kFSR_err_vsPt;   // error on p0
  
  // Get reference alpha bin value
  double refAlphaVal = balance3D_mc->GetZaxis()->GetBinLowEdge(refAlphaBin + 1);
  cout << "Reference alpha cut: alpha < " << refAlphaVal << " (bin " << refAlphaBin << ")" << endl;
  
  if (saveAlphaExtrap) {
    cout << "\n=== kFSR Extraction: fitting normalized ratio vs alpha ===" << endl;
    cout << "Normalized ratio = (MC/Data at α) / (MC/Data at ref α)" << endl;
    cout << "Reference alpha bin " << refAlphaBin << " excluded from fit (=1 by construction)" << endl;
    cout << "Fit range: " << fitAlphaMin << " < alpha < " << fitAlphaMax << endl;
    cout << "Saving alpha extrapolation histograms to " << alphaFolder << endl;
    
    TAxis* zaxis = balance3D_mc->GetZaxis();
    int nAlphaBins = zaxis->GetNbins();
    std::vector<double> alphaEdges(nAlphaBins + 1);
    for (int b = 1; b <= nAlphaBins + 1; ++b) alphaEdges[b-1] = zaxis->GetBinLowEdge(b);

    // Colors for multiple pT bins on overlay plots
    int ptColors[] = {kBlue, kRed, kGreen+2, kMagenta+2, kOrange+2, kCyan+2, kViolet+2, kTeal+2, kPink+2, kAzure+2};
    int nPtColors = sizeof(ptColors)/sizeof(ptColors[0]);

    // First compute reference alpha bin MC/Data ratio (for normalization)
    // Store reference ratio for each pT, eta bin
    map<pair<int,int>, double> refRatio;     // (ptbin, etabin) -> MC/Data ratio at ref alpha
    map<pair<int,int>, double> refRatioErr;  // error on ratio
    
    for (int ptbin = 1; ptbin <= nptbins; ++ptbin) {
      for (int etabin = 1; etabin <= netabins; ++etabin) {
        double mc_ref = balance3D_mc->GetBinContent(ptbin, etabin, refAlphaBin);
        double emc_ref = balance3D_mc->GetBinError(ptbin, etabin, refAlphaBin);
        double dt_ref = balance3D_data->GetBinContent(ptbin, etabin, refAlphaBin);
        double edt_ref = balance3D_data->GetBinError(ptbin, etabin, refAlphaBin);
        
        if (dt_ref > 0 && mc_ref > 0) {
          double ratio = mc_ref / dt_ref;
          double relErr2 = 0.0;
          if (mc_ref > 0 && emc_ref > 0) relErr2 += pow(emc_ref / mc_ref, 2);
          if (dt_ref > 0 && edt_ref > 0) relErr2 += pow(edt_ref / dt_ref, 2);
          refRatio[{ptbin, etabin}] = ratio;
          refRatioErr[{ptbin, etabin}] = ratio * sqrt(relErr2);
        } else {
          refRatio[{ptbin, etabin}] = 1.0;
          refRatioErr[{ptbin, etabin}] = 0.0;
        }
      }
    }

    // First loop: per-pT, per-eta alpha fits (detailed)
    // Optional (very verbose). Off by default: L3 is treated as pT-only.
    if (doEtaBinnedAlphaFits) {
    for (int ptbin = 1; ptbin <= nptbins; ++ptbin) {
      for (int etabin = 1; etabin <= netabins; ++etabin) {
        TH1D* hAlphaMc = new TH1D(Form("alpha_mc_pt%d_eta%d", ptbin, etabin), "MC balance vs #alpha;#alpha;Balance", nAlphaBins, alphaEdges.data());
        TH1D* hAlphaDt = new TH1D(Form("alpha_dt_pt%d_eta%d", ptbin, etabin), "Data balance vs #alpha;#alpha;Balance", nAlphaBins, alphaEdges.data());
        TH1D* hAlphaRatio = new TH1D(Form("alpha_ratio_pt%d_eta%d", ptbin, etabin), "Balance Ratio (MC/Data) vs #alpha;#alpha;Balance Ratio", nAlphaBins, alphaEdges.data());
        TH1D* hAlphaRatioNorm = new TH1D(Form("alpha_ratio_norm_pt%d_eta%d", ptbin, etabin), 
                                          "Normalized Ratio vs #alpha;#alpha;Normalized Ratio", nAlphaBins, alphaEdges.data());

        double ref_ratio = refRatio[{ptbin, etabin}];
        double ref_err = refRatioErr[{ptbin, etabin}];
        
        int nValidBins = 0;
        for (int abin = 1; abin <= nAlphaBins; ++abin) {
          double mc = balance3D_mc->GetBinContent(ptbin, etabin, abin);
          double emc = balance3D_mc->GetBinError(ptbin, etabin, abin);
          double dt = balance3D_data->GetBinContent(ptbin, etabin, abin);
          double edt = balance3D_data->GetBinError(ptbin, etabin, abin);

          hAlphaMc->SetBinContent(abin, mc);
          hAlphaMc->SetBinError(abin, emc);
          hAlphaDt->SetBinContent(abin, dt);
          hAlphaDt->SetBinError(abin, edt);

          if (dt > 0 && mc > 0) {
            double ratio = mc / dt;
            double relErr2 = 0.0;
            if (mc > 0 && emc > 0) relErr2 += pow(emc / mc, 2);
            if (dt > 0 && edt > 0) relErr2 += pow(edt / dt, 2);
            hAlphaRatio->SetBinContent(abin, ratio);
            hAlphaRatio->SetBinError(abin, ratio * sqrt(relErr2));
            
            // Normalized ratio = ratio / ref_ratio
            if (ref_ratio > 0) {
              double norm_ratio = ratio / ref_ratio;
              double norm_err = norm_ratio * sqrt(relErr2 + pow(ref_err/ref_ratio, 2));
              hAlphaRatioNorm->SetBinContent(abin, norm_ratio);
              // Errors are fully correlated for this normalized ratio; do not use/show them.
              hAlphaRatioNorm->SetBinError(abin, 0.0);
              
              // Don't count reference bin for fit
              if (abin != refAlphaBin) nValidBins++;
            }
          }
        }

        // Get readable pT and eta bin ranges for labeling
        double ptBinLo = balance3D_mc->GetXaxis()->GetBinLowEdge(ptbin);
        double ptBinHi = balance3D_mc->GetXaxis()->GetBinLowEdge(ptbin+1);
        double etaBinLo = balance3D_mc->GetYaxis()->GetBinLowEdge(etabin);
        double etaBinHi = balance3D_mc->GetYaxis()->GetBinLowEdge(etabin+1);

        // Fit NORMALIZED ratio vs alpha with linear function to extrapolate to alpha->0
        // Exclude reference bin from fit (it's 1 by construction)
        TF1* fAlpha = new TF1(Form("fAlpha_pt%d_eta%d", ptbin, etabin), "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
        fAlpha->SetParameters(1.0, 0.0);
        
        // Create TGraph excluding reference bin for fitting
        // Note: do NOT use y-errors here (fully correlated); treat all points equally.
        std::vector<double> x_fit, y_fit;
        for (int abin = 1; abin <= nAlphaBins; ++abin) {
          if (abin == refAlphaBin) continue; // Skip reference bin
          double alpha_center = zaxis->GetBinCenter(abin);
          double val = hAlphaRatioNorm->GetBinContent(abin);
          if (val > 0 && alpha_center >= fitAlphaMin && alpha_center <= fitAlphaMax) {
            x_fit.push_back(alpha_center);
            y_fit.push_back(val);
          }
        }
        
        TGraph* gFit = nullptr;
        if (x_fit.size() >= 2) {
          gFit = new TGraph(x_fit.size(), x_fit.data(), y_fit.data());
          gFit->Fit(fAlpha, "QNR");
        }

        TCanvas* cAlpha = new TCanvas(Form("cAlpha_pt%.0fto%.0f_eta%.3fto%.3f", ptBinLo, ptBinHi, etaBinLo, etaBinHi), 
                                       Form("cAlpha_pt%.0fto%.0f_eta%.3fto%.3f", ptBinLo, ptBinHi, etaBinLo, etaBinHi), 800, 600);
        cAlpha->cd();
        hAlphaRatioNorm->SetMinimum(0.97);
        hAlphaRatioNorm->SetMaximum(1.05);
        hAlphaRatioNorm->SetMarkerStyle(kFullCircle);
        hAlphaRatioNorm->SetMarkerColor(kBlue);
        hAlphaRatioNorm->SetLineColor(kBlue);
        const TString yNormPt = Form("#frac{(R_{MC}/R_{Data})(#alpha)}{(R_{MC}/R_{Data})(#alpha < %.2f)}", refAlphaVal);
        hAlphaRatioNorm->GetXaxis()->SetTitle("#alpha");
        hAlphaRatioNorm->GetXaxis()->SetTitleSize(0.038);
        hAlphaRatioNorm->GetXaxis()->SetLabelSize(0.032);
        hAlphaRatioNorm->GetYaxis()->SetTitle(yNormPt);
        hAlphaRatioNorm->GetYaxis()->SetTitleSize(0.038);
        hAlphaRatioNorm->GetYaxis()->SetLabelSize(0.032);
        hAlphaRatioNorm->GetYaxis()->SetTitleOffset(1.7);
        // Draw points only (no y-errors)
        hAlphaRatioNorm->Draw("P");
        
        // Mark reference bin differently
        TH1D* hRefPoint = new TH1D(Form("hRefPoint_pt%d_eta%d", ptbin, etabin), "", nAlphaBins, alphaEdges.data());
        hRefPoint->SetBinContent(refAlphaBin, hAlphaRatioNorm->GetBinContent(refAlphaBin));
        hRefPoint->SetBinError(refAlphaBin, hAlphaRatioNorm->GetBinError(refAlphaBin));
        hRefPoint->SetMarkerStyle(kFullSquare);
        hRefPoint->SetMarkerColor(kBlack);
        hRefPoint->SetLineColor(kBlack);
        hRefPoint->Draw("P SAME");
        
        fAlpha->SetLineColor(kRed);
        fAlpha->SetLineWidth(2);
        fAlpha->Draw("SAME");
        
        // Add reference line at 1.0
        TLine* lineRef = new TLine(alphaEdges.front(), 1.0, alphaEdges.back(), 1.0);
        lineRef->SetLineStyle(kDashed);
        lineRef->SetLineColor(kGray+1);
        lineRef->Draw("SAME");
        
        // Legend following dofits.C style
        TLegend* legAlpha = new TLegend(0.55, 0.70, 0.88, 0.88);
        legAlpha->SetBorderSize(0);
        legAlpha->SetFillStyle(0);
        legAlpha->SetTextFont(42);
        legAlpha->SetTextSize(0.030);
        legAlpha->AddEntry(hAlphaRatioNorm, "R (normalized)", "P");
        legAlpha->AddEntry(hRefPoint, Form("Reference (#alpha < %.2f, excluded)", refAlphaVal), "P");
        legAlpha->AddEntry(fAlpha, Form("Linear fit (%.2f < #alpha < %.2f)", fitAlphaMin, fitAlphaMax), "L");
        legAlpha->Draw();
        
        // Labels following dofits.C style
        TLatex* tbin = new TLatex();
        tbin->SetNDC();
        tbin->SetTextFont(42);
        tbin->SetTextSize(0.035);
        tbin->DrawLatex(0.18, 0.85, Form("%.0f < p_{T}^{#gamma} < %.0f GeV", ptBinLo, ptBinHi));
        tbin->DrawLatex(0.18, 0.80, Form("%.3f < |#eta_{jet}| < %.3f", etaBinLo, etaBinHi));
        tbin->DrawLatex(0.18, 0.75, Form("p0 = %.4f", fAlpha->GetParameter(0)));
        tbin->DrawLatex(0.18, 0.70, Form("p1 = %.4f", fAlpha->GetParameter(1)));
        
        CMS_lumi(cAlpha, 0, 0);
        
        string alphaFile = Form("%s/L3Res_%s_pt%.0fto%.0f_eta%.3fto%.3f_alpha.png", alphaFolder.c_str(), _run.c_str(), ptBinLo, ptBinHi, etaBinLo, etaBinHi);
        cAlpha->SaveAs(alphaFile.c_str());

        outfile->cd();
        hAlphaMc->Write();
        hAlphaDt->Write();
        hAlphaRatio->Write();
        hAlphaRatioNorm->Write();
        fAlpha->Write();

        delete lineRef;
        delete legAlpha;
        delete tbin;
        delete cAlpha;
        delete hAlphaMc;
        delete hAlphaDt;
        delete hAlphaRatio;
        delete hAlphaRatioNorm;
        delete hRefPoint;
        if (gFit) delete gFit;
        delete fAlpha;
      }
    }
    }
    
        // Second loop: kFSR extraction per pT bin.
        // L3 is treated as pT-only: choose a single wide |eta| bin (barrel) or collapse over eta.
        // Do NOT draw/use y-error bars for the normalized ratio points (fully correlated).
        cout << "\n=== Extracting kFSR vs photon pT (alpha->0 extrapolation) ===" << endl;
        cout << "kFSR extraction uses " << (useSingleEtaBin ? "single eta bin" : "eta-collapsed")
          << ", etaBinForL3=" << etaBinForL3 << endl;
    
    // Create histogram for kFSR vs pT
    std::vector<double> ptEdges(nptbins + 1);
    for (int b = 1; b <= nptbins + 1; ++b) ptEdges[b-1] = balance3D_mc->GetXaxis()->GetBinLowEdge(b);
    
    TH1D* hkFSR = new TH1D("kFSR_vsPt", "k_{FSR} (Normalized Ratio extrapolated to #alpha#rightarrow0) vs p_{T}^{#gamma};p_{T}^{#gamma} (GeV);k_{FSR}", 
                           nptbins, ptEdges.data());
    
    auto ratioAt = [&](int ptbin, int alphabin) -> double {
      if (ptbin < 1 || ptbin > nptbins || alphabin < 1 || alphabin > nalphabins) return 0.0;
      double mc = 0.0, dt = 0.0;
      if (useSingleEtaBin) {
        int eb = etaBinForL3;
        if (eb < 1) eb = 1;
        if (eb > netabins) eb = netabins;
        mc = balance3D_mc->GetBinContent(ptbin, eb, alphabin);
        dt = balance3D_data->GetBinContent(ptbin, eb, alphabin);
      } else {
        double sum_mc = 0.0, ent_mc = 0.0;
        double sum_dt = 0.0, ent_dt = 0.0;
        for (int eb = 1; eb <= netabins; ++eb) {
          const double v_mc = balance3D_mc->GetBinContent(ptbin, eb, alphabin);
          const double n_mc = balance3D_mc->GetBinEntries(balance3D_mc->GetBin(ptbin, eb, alphabin));
          const double v_dt = balance3D_data->GetBinContent(ptbin, eb, alphabin);
          const double n_dt = balance3D_data->GetBinEntries(balance3D_data->GetBin(ptbin, eb, alphabin));
          if (v_mc > 0 && n_mc > 0 && !TMath::IsNaN(v_mc)) { sum_mc += v_mc * n_mc; ent_mc += n_mc; }
          if (v_dt > 0 && n_dt > 0 && !TMath::IsNaN(v_dt)) { sum_dt += v_dt * n_dt; ent_dt += n_dt; }
        }
        mc = (ent_mc > 0) ? (sum_mc / ent_mc) : 0.0;
        dt = (ent_dt > 0) ? (sum_dt / ent_dt) : 0.0;
      }
      if (dt <= 0 || mc <= 0) return 0.0;
      return mc / dt;
    };

    // Create overlay plot with all pT bins on same canvas
    TCanvas* cAlphaOverlay = new TCanvas("cAlphaOverlay", "Normalized Ratio vs alpha (all pT bins)", 900, 700);
    cAlphaOverlay->cd();
    const TString yNormTitle = Form("#frac{(R_{MC}/R_{Data})(#alpha)}{(R_{MC}/R_{Data})(#alpha < %.2f)}", refAlphaVal);
    TH1D* hFrameAlpha = new TH1D("hFrameAlpha", Form(";#alpha;%s", yNormTitle.Data()), 100, 0, 0.5);
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
    lineRefOverlay->SetLineColor(kGray+1);
    lineRefOverlay->Draw("SAME");
    
    TLegend* legOverlay = new TLegend(0.55, 0.55, 0.88, 0.88);
    legOverlay->SetBorderSize(0);
    legOverlay->SetFillStyle(0);
    legOverlay->SetTextFont(42);
    legOverlay->SetTextSize(0.025);
    
    // Store points for each pT bin
    map<int, TGraph*> gAlphaNorm;
    
    for (int ptbin = 1; ptbin <= nptbins; ++ptbin) {
      double ptBinLo = balance3D_mc->GetXaxis()->GetBinLowEdge(ptbin);
      double ptBinHi = balance3D_mc->GetXaxis()->GetBinLowEdge(ptbin+1);
      const double ref_ratio = ratioAt(ptbin, refAlphaBin);
      if (ref_ratio <= 0) continue;

      std::vector<double> x_all;
      std::vector<double> y_all;
      std::vector<double> x_fit;
      std::vector<double> y_fit;

      for (int abin = 1; abin <= nAlphaBins; ++abin) {
        const double a = zaxis->GetBinCenter(abin);
        const double ratio = ratioAt(ptbin, abin);
        if (ratio <= 0) continue;
        const double y = ratio / ref_ratio;
        x_all.push_back(a);
        y_all.push_back(y);
        if (abin != refAlphaBin && a >= fitAlphaMin && a <= fitAlphaMax) {
          x_fit.push_back(a);
          y_fit.push_back(y);
        }
      }

      if (x_all.size() < 2 || x_fit.size() < 2) continue;

      const int colorIdx = (ptbin - 1) % nPtColors;
      TGraph* gAll = new TGraph(x_all.size(), x_all.data(), y_all.data());
      gAll->SetName(Form("gAlphaNorm_pt%d", ptbin));
      gAll->SetMarkerStyle(kFullCircle);
      gAll->SetMarkerColor(ptColors[colorIdx]);
      gAll->SetLineColor(ptColors[colorIdx]);
      gAll->Draw("P SAME");
      gAlphaNorm[ptbin] = gAll;

      TGraph gFit(x_fit.size(), x_fit.data(), y_fit.data());
      TF1 fAlphaCol(Form("fAlpha_pt%d_col", ptbin), "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
      fAlphaCol.SetParameters(1.0, 0.0);
      gFit.Fit(&fAlphaCol, "QNR");

      const double kfsr = fAlphaCol.GetParameter(0);
      const double kfsr_err = fAlphaCol.GetParError(0);
      kFSR_vsPt[ptbin] = kfsr;
      kFSR_err_vsPt[ptbin] = kfsr_err;
      hkFSR->SetBinContent(ptbin, kfsr);
      hkFSR->SetBinError(ptbin, kfsr_err);
      cout << "  pT [" << ptBinLo << "-" << ptBinHi << "]: kFSR = " << kfsr << " +/- " << kfsr_err << endl;

      // Save a dedicated per-pT plot showing the points and the linear alpha fit.
      {
        TCanvas* cPt = new TCanvas(Form("cAlphaFit_pt%d", ptbin), Form("Normalized ratio vs alpha (ptbin %d)", ptbin), 800, 600);
        cPt->cd();
        const TString yNormTitlePt = Form("#frac{(R_{MC}/R_{Data})(#alpha)}{(R_{MC}/R_{Data})(#alpha < %.2f)}", refAlphaVal);
        TH1D* hFramePt = new TH1D(Form("hFrameAlphaFit_pt%d", ptbin), Form(";#alpha;%s", yNormTitlePt.Data()), 100, 0, 0.5);
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
        l1->SetLineColor(kGray+1);
        l1->Draw("SAME");

        gAll->Draw("P SAME");
        TGraph* gFitPts = new TGraph(x_fit.size(), x_fit.data(), y_fit.data());
        gFitPts->SetMarkerStyle(kOpenCircle);
        gFitPts->SetMarkerColor(kRed+1);
        gFitPts->SetLineColor(kRed+1);
        gFitPts->Draw("P SAME");

        TF1* fLine = new TF1(Form("fAlphaLine_pt%d", ptbin), "[0]+[1]*x", fitAlphaMin, fitAlphaMax);
        fLine->SetParameters(fAlphaCol.GetParameter(0), fAlphaCol.GetParameter(1));
        fLine->SetLineColor(kRed);
        fLine->SetLineWidth(2);
        fLine->Draw("SAME");

        TLatex* t = new TLatex();
        t->SetNDC();
        t->SetTextFont(42);
        t->SetTextSize(0.035);
        t->DrawLatex(0.18, 0.86, Form("p_{T}^{#gamma}: %.0f-%.0f GeV", ptBinLo, ptBinHi));
        t->DrawLatex(0.18, 0.81, Form("Ref #alpha < %.2f (bin %d), excluded", refAlphaVal, refAlphaBin));
        t->DrawLatex(0.18, 0.76, Form("Fit: %.2f < #alpha < %.2f", fitAlphaMin, fitAlphaMax));
        t->DrawLatex(0.18, 0.71, Form("k_{FSR} = %.4f #pm %.4f", kfsr, kfsr_err));
        CMS_lumi(cPt, 0, 0);
        cPt->SaveAs(Form("%s/L3Res_%s_kFSR_alphaFit_pt%.0fto%.0f.png", alphaFolder.c_str(), _run.c_str(), ptBinLo, ptBinHi));

        outfile->cd();
        gFitPts->Write(Form("gAlphaNorm_fitPts_pt%d", ptbin));
        fLine->Write();

        delete t;
        delete fLine;
        delete gFitPts;
        delete l1;
        delete hFramePt;
        delete cPt;
      }

      legOverlay->AddEntry(gAll, Form("%.0f-%.0f GeV (k_{FSR}=%.3f)", ptBinLo, ptBinHi, kfsr), "P");

      outfile->cd();
      gAll->Write();
      fAlphaCol.Write();
    }

    legOverlay->Draw();
    
    TLatex* texOverlay = new TLatex();
    texOverlay->SetNDC();
    texOverlay->SetTextFont(42);
    texOverlay->SetTextSize(0.030);
    texOverlay->DrawLatex(0.18, 0.85, Form("Denominator: (#alpha < %.2f), bin %d (excluded)", refAlphaVal, refAlphaBin));
    texOverlay->DrawLatex(0.18, 0.80, Form("Fit: %.2f < #alpha < %.2f", fitAlphaMin, fitAlphaMax));
    
    CMS_lumi(cAlphaOverlay, 0, 0);
    cAlphaOverlay->SaveAs(Form("%s/L3Res_%s_alpha_overlay.png", alphaFolder.c_str(), _run.c_str()));
    
    outfile->cd();
    
    delete cAlphaOverlay;
    delete hFrameAlpha;
    delete lineRefOverlay;
    delete legOverlay;
    delete texOverlay;
    
    // Create kFSR vs pT diagnostic plot (no functional fit)
    // kFSR is defined per pT bin from the linear alpha-extrapolation.
    // Only the FINAL correction is fit vs pT.
    cout << "\n=== kFSR vs photon pT (diagnostic; no fit) ===" << endl;
    
    TCanvas* ckFSR = new TCanvas("ckFSR", "kFSR vs photon pT", 800, 600);
    ckFSR->SetLogx();
    ckFSR->cd();
    
    TH1D* hFramekFSR = makePtFrameFromHist(hkFSR, "hFramekFSR", "k_{FSR} (\"#alpha#rightarrow0\" intercept)", 0.9, 1.1);
    hFramekFSR->Draw();
    styleLogxAxis(hFramekFSR);

    const double xMinK = hFramekFSR->GetXaxis()->GetBinLowEdge(1);
    const double xMaxK = hFramekFSR->GetXaxis()->GetBinLowEdge(hFramekFSR->GetNbinsX() + 1);
    TLine* linekFSR = new TLine(xMinK, 1.0, xMaxK, 1.0);
    linekFSR->SetLineStyle(kDashed);
    linekFSR->SetLineColor(kGray+1);
    linekFSR->Draw("SAME");
    
    hkFSR->SetMarkerStyle(kFullCircle);
    hkFSR->SetMarkerColor(kBlue);
    hkFSR->SetLineColor(kBlue);
    hkFSR->Draw("PE1 SAME");
    
    TLegend* legkFSR = new TLegend(0.50, 0.65, 0.88, 0.88);
    legkFSR->SetBorderSize(0);
    legkFSR->SetFillStyle(0);
    legkFSR->SetTextFont(42);
    legkFSR->SetTextSize(0.030);
    legkFSR->AddEntry(hkFSR, "k_{FSR} (#alpha#rightarrow0 extrap.)", "PLE");
    legkFSR->Draw();
    
    TLatex* texkFSR = new TLatex();
    texkFSR->SetNDC();
    texkFSR->SetTextFont(42);
    texkFSR->SetTextSize(0.035);
    texkFSR->DrawLatex(0.18, 0.30, Form("Ref #alpha < %.2f (bin %d)", refAlphaVal, refAlphaBin));
    texkFSR->DrawLatex(0.18, 0.25, "No functional fit applied");
    
    CMS_lumi(ckFSR, 0, 0);
    ckFSR->SaveAs(Form("%s/L3Res_%s_kFSR_vspT.png", alphaFolder.c_str(), _run.c_str()));
    
    outfile->cd();
    hkFSR->Write();
    
    delete ckFSR;
    delete hFramekFSR;
    delete linekFSR;
    delete legkFSR;
    delete texkFSR;
  }

  cout << "\n=== Loading pT-dependent response ratios ===" << endl;
  map<int, TH1D*> ratioVsPt;
  for (int alpha = 1; alpha <= nalphabins; ++alpha) {
    TH1D* h = (TH1D*)inFile->Get(Form("ratio_vsphotonpt_alpha%d", alpha));
    if (h) {
      ratioVsPt[alpha] = h;
      cout << "  Loaded ratio_vsphotonpt_alpha" << alpha << endl;
    }
  }
  if (ratioVsPt.empty()) {
    cout << "ERROR: No ratio histograms available" << endl;
    return;
  }

  // If kFSR was extracted in this run, build the alpha->0 correction histogram explicitly:
  //   Corr_alpha0(pT) = Ratio_refAlpha(pT) * kFSR(pT)
  // where Ratio_refAlpha is the MC/Data balance ratio at the reference alpha cut.
  TH1D* hCorrAlpha0 = nullptr;
  if (saveAlphaExtrap && applyKFSRToPtFit && ratioVsPt.count(refAlphaBin)) {
    TH1D* hRef = ratioVsPt[refAlphaBin];
    if (hRef) {
      hCorrAlpha0 = (TH1D*)hRef->Clone("ratio_vspT_alpha0_fromkFSR");
      hCorrAlpha0->SetDirectory(nullptr);
      for (int b = 1; b <= hCorrAlpha0->GetNbinsX(); ++b) {
        const double r = hRef->GetBinContent(b);
        const double er = hRef->GetBinError(b);
        const double k = (kFSR_vsPt.count(b) ? kFSR_vsPt[b] : 1.0);
        const double ek = (kFSR_err_vsPt.count(b) ? kFSR_err_vsPt[b] : 0.0);
        if (r > 0 && k > 0) {
          const double val = r * k;
          double rel2 = 0.0;
          if (er > 0) rel2 += pow(er / r, 2);
          if (ek > 0) rel2 += pow(ek / k, 2);
          hCorrAlpha0->SetBinContent(b, val);
          hCorrAlpha0->SetBinError(b, val * sqrt(rel2));
        } else {
          hCorrAlpha0->SetBinContent(b, 0.0);
          hCorrAlpha0->SetBinError(b, 0.0);
        }
      }

      // Plot the reference-alpha ratio and the kFSR-corrected (alpha->0) ratio.
      TCanvas* cCorr = new TCanvas("cCorrAlpha0", "Ref-alpha ratio and kFSR-corrected alpha->0 ratio", 900, 700);
      cCorr->SetLogx();
      cCorr->cd();
      TH1D* hFrame = makePtFrameFromHist(hRef, "hFrameCorrAlpha0", "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
      hFrame->Draw();
      styleLogxAxis(hFrame);
      const double xMin = hFrame->GetXaxis()->GetBinLowEdge(1);
      const double xMax = hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetNbinsX() + 1);
      TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
      line->SetLineStyle(kDashed);
      line->SetLineColor(kGray+1);
      line->Draw("SAME");

      hRef->SetMarkerStyle(kOpenCircle);
      hRef->SetMarkerColor(kBlue);
      hRef->SetLineColor(kBlue);
      hRef->Draw("PE1 SAME");

      hCorrAlpha0->SetMarkerStyle(kFullCircle);
      hCorrAlpha0->SetMarkerColor(kRed);
      hCorrAlpha0->SetLineColor(kRed);
      hCorrAlpha0->Draw("PE1 SAME");

      TLegend* leg = new TLegend(0.50, 0.72, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.030);
      leg->AddEntry(hRef, Form("Ref #alpha < %.2f (bin %d)", refAlphaVal, refAlphaBin), "PLE");
      leg->AddEntry(hCorrAlpha0, "k_{FSR}(p_{T}) #times ref(#alpha)  (#alpha#rightarrow0)", "PLE");
      leg->Draw();
      CMS_lumi(cCorr, 0, 0);
      cCorr->SaveAs(Form("%s/L3Res_%s_corr_alpha0_fromkFSR.png", alphaFolder.c_str(), _run.c_str()));

      outfile->cd();
      hCorrAlpha0->Write("ratio_vspT_alpha0_fromkFSR");

      delete leg;
      delete line;
      delete hFrame;
      delete cCorr;
    }
  }

  // Overlay unnormalized balance ratios (MC/Data) vs pT for all alpha cuts.
  // This is a diagnostic plot only (no normalization to the reference alpha bin).
  {
    TCanvas* c = new TCanvas("cRatioOverlayPt", "Unnormalized ratio vs pT (all alpha cuts)", 900, 700);
    c->SetLogx();
    c->cd();

    TH1D* href = ratioVsPt.begin()->second;
    TH1D* hFrame = makePtFrameFromHist(href, Form("hFrame_ratioOverlayPt_%s", outfilename.c_str()), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
    hFrame->Draw();
    styleLogxAxis(hFrame);

    const double xMin = hFrame->GetXaxis()->GetBinLowEdge(1);
    const double xMax = hFrame->GetXaxis()->GetBinLowEdge(hFrame->GetNbinsX() + 1);
    TLine* line = new TLine(xMin, 1.0, xMax, 1.0);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray+1);
    line->Draw("SAME");

    const int colors[] = {kBlack, kBlue, kRed, kGreen+2, kMagenta+2, kOrange+2, kCyan+2, kViolet+2, kTeal+2};
    const int nColors = sizeof(colors) / sizeof(colors[0]);

    TLegend* leg = new TLegend(0.52, 0.60, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    TAxis* zax = balance3D_mc->GetZaxis();
    int idx = 0;
    for (const auto& kv : ratioVsPt) {
      const int alphaBin = kv.first;
      TH1D* h = kv.second;
      if (!h) continue;
      const int col = colors[idx % nColors];
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerSize(0.8);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->Draw("PE1 SAME");
      const double alphaCut = zax->GetBinLowEdge(alphaBin + 1);
      leg->AddEntry(h, Form("#alpha < %.2f", alphaCut), "PLE");
      ++idx;
    }

    leg->Draw();
    CMS_lumi(c, 0, 0);
    c->SaveAs(Form("%s/L3Res_%s_ratio_overlay_vspt.png", pngFolder.c_str(), _run.c_str()));

    delete leg;
    delete line;
    delete hFrame;
    delete c;
  }

  // --- Build per-input pT correction(s) and do one combined pT fit (multi-input) ---
  // For multi-input workflows (photon+jet + Z+jet, etc.), provide a comma-separated list
  // of derived ROOT files as the first argument. This macro will fit all pT points together.
  std::vector<TString> inputFiles2 = splitInputFiles(inFileL3Derived);
  if (inputFiles2.empty()) inputFiles2.push_back(inFileL3Derived);

  const bool isMultiInput2 = (inputFiles2.size() > 1);
  const bool doCombined = (doCombinedPtFit && isMultiInput2);

  // Collect pT correction histograms from each input
  std::vector<TH1D*> hCorrInputs;
  std::vector<TString> labelInputs;
  if (!doCombined && saveAlphaExtrap && applyKFSRToPtFit) {
    cout << "\n=== Building alpha->0 correction for final pT fit ===" << endl;
    cout << "Using: Corr_{alpha->0}(pT) = Ratio_{ref}(pT) * kFSR(pT)" << endl;
    cout << "  - Ratio_ref(pT) is MC/Data at the reference alpha cut (ref bin not used in alpha fit)" << endl;
    cout << "  - kFSR(pT) is the intercept from the linear fit of normalized ratio vs alpha" << endl;
  }
  for (size_t i = 0; i < inputFiles2.size(); ++i) {
    if (!doCombined && i > 0) break;
    TFile* fIn = TFile::Open(inputFiles2[i], "READ");
    if (!fIn || fIn->IsZombie()) {
      cout << "WARNING: cannot open input " << inputFiles2[i] << ", skipping" << endl;
      continue;
    }
    TH1D* hIn = loadPtCorrectionHist(fIn, refAlphaBin);
    if (!hIn) {
      cout << "WARNING: could not find pT correction hist in " << inputFiles2[i]
           << " (expected corr_vspT or ratio_vspT_alpha0 or ratio_vsphotonpt_alphaN); skipping" << endl;
      fIn->Close();
      continue;
    }

    TH1D* hClone = (TH1D*)hIn->Clone(Form("corr_input_%zu", i));
    hClone->SetDirectory(nullptr);

    // If this is single-input and kFSR was computed in this run, optionally apply it.
    // For multi-input, provide corr_vspT / ratio_vspT_alpha0 from the corresponding derive step.
    if (!doCombined && applyKFSRToPtFit && saveAlphaExtrap) {
      const std::string hname = hIn->GetName() ? std::string(hIn->GetName()) : std::string();
      const bool alreadyAlpha0 = (hname.find("ratio_vspT_alpha0") != std::string::npos) || (hname.find("corr_vspT") != std::string::npos);
      if (alreadyAlpha0) {
        cout << "NOTE: input hist '" << hname << "' already represents an alpha->0 correction; not applying kFSR again." << endl;
      } else {
        cout << "Applying kFSR(pT) to input hist '" << hname << "' to build alpha->0 correction." << endl;
      for (int b = 1; b <= hClone->GetNbinsX(); ++b) {
        const double r = hClone->GetBinContent(b);
        const double er = hClone->GetBinError(b);
        const double k = (kFSR_vsPt.count(b) ? kFSR_vsPt[b] : 1.0);
        const double ek = (kFSR_err_vsPt.count(b) ? kFSR_err_vsPt[b] : 0.0);
        if (r > 0 && k > 0) {
          const double val = r * k;
          double rel2 = 0.0;
          if (er > 0) rel2 += pow(er / r, 2);
          if (ek > 0) rel2 += pow(ek / k, 2);
          hClone->SetBinContent(b, val);
          hClone->SetBinError(b, val * sqrt(rel2));
        } else {
          hClone->SetBinContent(b, 0.0);
          hClone->SetBinError(b, 0.0);
        }
      }
      }
    }

    hCorrInputs.push_back(hClone);
    TString lbl = (i < inputLabels.size()) ? inputLabels[i] : defaultInputLabel(inputFiles2[i], (int)i + 1);
    labelInputs.push_back(lbl);
    fIn->Close();
  }

  if (hCorrInputs.empty()) {
    cout << "ERROR: no valid inputs for pT fit." << endl;
    return;
  }

  // Keep a handle to the final pT-fit function for optional L2L3 output
  TF1* frefFinal = nullptr;

  // Plot + fit combined pT correction
  {
    TCanvas* cFinal = new TCanvas("cFinal_combined", "L3Residual combined pT fit", 900, 700);
    cFinal->SetLogx();
    // Set margins BEFORE drawing the frame, otherwise axis titles/labels can get clipped
    cFinal->SetLeftMargin(0.14);
    cFinal->SetBottomMargin(0.12);
    cFinal->SetRightMargin(0.05);
    cFinal->SetTopMargin(0.06);
    cFinal->cd();

    const double xMin = hCorrInputs.front()->GetXaxis()->GetBinLowEdge(1);
    const double xMax = hCorrInputs.front()->GetXaxis()->GetBinLowEdge(hCorrInputs.front()->GetNbinsX() + 1);
    const double yMin = 0.7;
    const double yMax = 1.5;

    TMultiGraph* mg = new TMultiGraph("mg_combined_pt", "mg_combined_pt");
    TLegend* leg = new TLegend(0.50, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.028);

    const int colors[] = {kBlue, kRed, kGreen+2, kMagenta+2, kOrange+2, kCyan+2, kViolet+2, kTeal+2};
    const int nColors = sizeof(colors) / sizeof(colors[0]);
    const int markers[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross, kFullStar, kFullCircle};
    const int nMarkers = sizeof(markers) / sizeof(markers[0]);

    std::vector<TH1D*> cleaned;
    std::vector<TGraphErrors*> graphs;
    std::vector<TGraphErrors*> fitGraphs;
    for (size_t i = 0; i < hCorrInputs.size(); ++i) {
      TH1D* h = hCorrInputs[i];
      if (!h) continue;
      const int col = colors[i % nColors];
      const int mkr = markers[i % nMarkers];

      double ptminThis = ptminG;
      double ptmaxThis = ptmaxG;
      if (i < inputPtRanges.size() && inputPtRanges[i].Length() > 0) {
        TString r = inputPtRanges[i];
        r.ReplaceAll(" ", "");
        if (r.Contains("-")) {
          TString a = r(0, r.First('-'));
          TString b = r(r.First('-') + 1, r.Length());
          const double lo = a.Atof();
          const double hi = b.Atof();
          if (lo > 0) ptminThis = std::max(ptminThis, lo);
          if (hi > 0) ptmaxThis = std::min(ptmaxThis, hi);
        }
      }

      // Plot the full distribution for display (no window cut so display bins are preserved)
      TH1D* hAll = drawCleaned(h, "G", xMin, xMax, mkr, col, false);
      if (hAll) cleaned.push_back(hAll);
      TGraphErrors* gAll = (hAll ? cleanGraph(new TGraphErrors(hAll)) : nullptr);
      if (gAll && gAll->GetN() > 0) {
        gAll->SetMarkerStyle(mkr);
        gAll->SetMarkerColor(col);
        gAll->SetLineColor(col);
        gAll->SetMarkerSize(1.05);
        graphs.push_back(gAll);
        mg->Add(gAll, "P");
      }

      // Also prepare the in-fit-range graph used for fitting (do not add to mg)
      TH1D* hIn = drawCleaned(h, "G", ptminThis, ptmaxThis, mkr, col, true);
      TGraphErrors* gIn = (hIn ? cleanGraph(new TGraphErrors(hIn)) : nullptr);
      if (gIn && gIn->GetN() > 0) {
        gIn->SetMarkerStyle(mkr);
        gIn->SetMarkerColor(col);
        gIn->SetLineColor(col);
        gIn->SetMarkerSize(1.05);
        fitGraphs.push_back(gIn);
      }

      // Report counts for diagnostics
      {
        const int nAll = (gAll ? gAll->GetN() : 0);
        const int nIn = (gIn ? gIn->GetN() : 0);
        cout << Form("INFO: input '%s' -> plotted points = %d, in-range points = %d\n", labelInputs[i].Data(), nAll, nIn);
      }

      // Legend: include fit range text in the label (visual only)
      if (gAll && gAll->GetN() > 0) {
        leg->AddEntry(gAll, Form("%s (%.0f-%.0f)", labelInputs[i].Data(), ptminThis, ptmaxThis), "P");
      }

    }

    // Draw axes from the multigraph so titles and log-x behave like other plots
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

    // Draw full distributions (colored points) via mg
    mg->Draw("AP");
    mg->Draw("P SAME");

    // Horizontal reference line at 1
    TLine* line = new TLine(xMin, 1, xMax, 1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray+1);
    line->Draw("SAME");

    gPad->RedrawAxis();

    // Build a concatenated graph containing ONLY the in-fit-range points and fit that
    TF1* fref = makeL3PtFitFunc("fref_combined", 15., 3500.);
    TGraphErrors* gFitTotal = new TGraphErrors();
    int nTot = 0;
    for (auto* g : fitGraphs) {
      if (!g) continue;
      for (int ip = 0; ip < g->GetN(); ++ip) {
        const double x = g->GetX()[ip];
        const double y = g->GetY()[ip];
        const double ex = (g->GetEX() ? g->GetEX()[ip] : 0.0);
        const double ey = (g->GetEY() ? g->GetEY()[ip] : 0.0);
        gFitTotal->SetPoint(nTot, x, y);
        gFitTotal->SetPointError(nTot, ex, ey);
        ++nTot;
      }
    }
    if (gFitTotal->GetN() >= 2) {
      gFitTotal->Fit(fref, "QRN");
      frefFinal = fref;
      fref->SetLineColor(kRed+1);
      fref->SetLineWidth(2);
      fref->Draw("SAME");
    } else {
      cout << "WARNING: Not enough in-range points to perform combined pT fit (N=" << gFitTotal->GetN() << ")" << endl;
      delete gFitTotal;
    }

    // Re-apply y-range and axis styling after the fit (fit can trigger pad autoscaling)
    mg->SetMinimum(yMin);
    mg->SetMaximum(yMax);
    mg->GetYaxis()->SetRangeUser(yMin, yMax);
    if (mg->GetHistogram()) {
      mg->GetHistogram()->SetMinimum(yMin);
      mg->GetHistogram()->SetMaximum(yMax);
      styleLogxAxis(mg->GetHistogram());
    }
    // Force canvas update and redraw axes/ticks so labels appear reliably
    cFinal->Modified();
    cFinal->Update();
    gPad->RedrawAxis();

    leg->Draw();

    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.032);
    tex->DrawLatex(0.18, 0.85, Form("Denominator: (#alpha < %.2f), bin %d", refAlphaVal, refAlphaBin));
    tex->DrawLatex(0.18, 0.80, Form("Inputs: %zu", hCorrInputs.size()));
    tex->DrawLatex(0.18, 0.75, Form("p0=%.5f, p1=%.5f", fref->GetParameter(0), fref->GetParameter(1)));
    tex->DrawLatex(0.18, 0.70, Form("Fit: %s", l3PtFitExpr()));
    tex->DrawLatex(0.18, 0.65, Form("#chi^{2}/ndf = %.1f/%d", fref->GetChisquare(), fref->GetNDF()));
    CMS_lumi(cFinal, 0, 0);
    cFinal->SaveAs(Form("%s/L3Res_%s_combined_final.png", pngFolder.c_str(), _run.c_str()));

    outfile->cd();
    mg->Write("mg_combined_pt");
    fref->Write("fit_combined_ref");
    for (size_t i = 0; i < hCorrInputs.size(); ++i) {
      if (hCorrInputs[i]) hCorrInputs[i]->Write(Form("corr_input_%zu", i));
    }

    // Write final L3Residual text in JEC format.
    // If an L2Residual file is provided, match its eta binning + pT validity ranges so it can be used consistently downstream.
    {
      std::string l2HeaderTmp;
      std::vector<JecRecord> l2recsTmp;
      const bool haveL2 = (!l2ResidualFile.empty() && readSimpleJecFile(l2ResidualFile, l2HeaderTmp, l2recsTmp) && !l2recsTmp.empty());

      const std::string l3TxtOut = Form("%s/%s.txt", txtFolder.c_str(), outfilename.c_str());
      std::ofstream ftxt(l3TxtOut);
      ftxt << "{ 1 JetEta 1 JetPt " << l3PtFitExpr() << " Correction L3Residual}\n";
      if (haveL2) {
        for (const auto& rec : l2recsTmp) {
          ftxt << Form("  %6.3f %6.3f  4  %5.0f %5.0f  %10.6f %10.6f\n",
                       rec.etaMin, rec.etaMax, rec.ptMin, rec.ptMax,
                       fref->GetParameter(0), fref->GetParameter(1));
        }
      } else {
        const double etaMinAll = balance3D_mc->GetYaxis()->GetBinLowEdge(1);
        const double etaMaxAll = balance3D_mc->GetYaxis()->GetBinLowEdge(netabins + 1);
        const double ptLo = hCorrInputs.front()->GetXaxis()->GetBinLowEdge(1);
        const double ptHi = hCorrInputs.front()->GetXaxis()->GetBinLowEdge(hCorrInputs.front()->GetNbinsX() + 1);
        ftxt << Form("  %6.3f %6.3f  4  %5.0f %5.0f  %10.6f %10.6f\n",
                     etaMinAll, etaMaxAll, ptLo, ptHi,
                     fref->GetParameter(0), fref->GetParameter(1));
      }
      ftxt.close();

      gSystem->mkdir(jecOutDir.c_str(), kTRUE);
      const std::string outTag = doCombined ? "combined" : "photonjet";
      const std::string l3TxtJec = Form("%s/L3Residuals_%s_%s_AK4PF.txt", jecOutDir.c_str(), _run.c_str(), outTag.c_str());
      std::ofstream ftxt2(l3TxtJec);
      ftxt2 << "{ 1 JetEta 1 JetPt " << l3PtFitExpr() << " Correction L3Residual}\n";
      if (haveL2) {
        // Use L2 eta structure but pT range is the fit range (ptminG, ptmaxG)
        for (const auto& rec : l2recsTmp) {
          ftxt2 << Form("  %6.3f %6.3f  4  %5.0f %5.0f  %10.6f %10.6f\n",
                        rec.etaMin, rec.etaMax, ptminG, ptmaxG,
                        fref->GetParameter(0), fref->GetParameter(1));
        }
      } else {
        const double etaMinAll = balance3D_mc->GetYaxis()->GetBinLowEdge(1);
        const double etaMaxAll = balance3D_mc->GetYaxis()->GetBinLowEdge(netabins + 1);
        ftxt2 << Form("  %6.3f %6.3f  4  %5.0f %5.0f  %10.6f %10.6f\n",
                      etaMinAll, etaMaxAll, ptminG, ptmaxG,
                      fref->GetParameter(0), fref->GetParameter(1));
      }
      ftxt2.close();

      cout << "\nWrote L3Residual JEC text: " << l3TxtJec << endl;
    }

    // cleanup (ROOT will own objects written to file; delete local-only)
    delete tex;
    delete leg;
    delete line;
    delete cFinal;
    // graphs in mg were heap-allocated; ok to keep until process exit
  }

  // Optional: build combined L2L3Residual by multiplying the current L2Residual file with the fitted L3Residual(pt)
  if (writeL2L3) {
    if (!frefFinal) {
      cout << "WARNING: writeL2L3 requested but no final L3 pT fit is available; skipping L2L3 write" << endl;
    } else {
      std::string l2Header;
      std::vector<JecRecord> l2recs;
      if (!readSimpleJecFile(l2ResidualFile, l2Header, l2recs)) {
        cout << "WARNING: Could not read L2Residual file: " << l2ResidualFile << " (skipping L2L3 write)" << endl;
      } else {
        const std::string l2l3Txt = Form("%s/L2L3Residuals_%s_photonjet_AK4PF.txt", jecOutDir.c_str(), _run.c_str());
        std::ofstream out(l2l3Txt);
        // Preserve the L2 file structure (eta bins, pT validity range, and functional form).
        // We refit the product (L2(pt)*L3(pt)) to the SAME functional form as the L2 file.
        std::string l2l3Header = l2Header;
        if (l2l3Header.find("L2Relative") != std::string::npos) {
          size_t p = l2l3Header.find("L2Relative");
          l2l3Header.replace(p, std::string("L2Relative").size(), "L2L3Residual");
        }
        out << l2l3Header << "\n";

        // Use hardcoded L2 expression (matches the L2 file used in this repo)
        const char* l2Expr = "1./([0]+[1]*log10(0.01*x)+[2]/(x/10.0))";
        // Fit combined shape per eta-bin using sampled points in the L2 pT validity range
        for (const auto& rec : l2recs) {
          const double ptLo = rec.ptMin;
          const double ptHi = rec.ptMax;
          const int nSamples = 30;
          std::vector<double> x(nSamples), y(nSamples);

          // Evaluate L2 using the parameters from the file
          TF1 fL2("fL2", l2Expr, ptLo, ptHi);
          for (size_t ip = 0; ip < rec.par.size() && ip < 3; ++ip) fL2.SetParameter((int)ip, rec.par[ip]);

          for (int ip = 0; ip < nSamples; ++ip) {
            const double t = (ip + 0.5) / nSamples;
            const double pt = ptLo * pow(ptHi / ptLo, t); // log-uniform sampling
            const double c2 = fL2.Eval(pt);
            const double c3 = frefFinal->Eval(pt);
            x[ip] = pt;
            y[ip] = c2 * c3;
          }
          TGraph g(nSamples, x.data(), y.data());
          TF1 fComb("fComb", l2Expr, ptLo, ptHi);
          for (size_t ip = 0; ip < rec.par.size() && ip < 3; ++ip) fComb.SetParameter((int)ip, rec.par[ip]);
          g.Fit(&fComb, "QNR");

          std::ostringstream oss;
          oss.setf(std::ios::fixed);
          oss << Form("  %6.3f %6.3f  %d  %5.0f %5.0f", rec.etaMin, rec.etaMax, rec.nPar, ptLo, ptHi);
          for (int ip = 0; ip < rec.nPar - 2; ++ip) {
            oss << Form("  %10.6f", fComb.GetParameter(ip));
          }
          oss << "\n";
          out << oss.str();
        }
        out.close();
        cout << "Wrote combined L2L3Residual JEC text: " << l2l3Txt << endl;
      }
      }
    }

  cout << "\n=== Creating plots for all alpha bins ===" << endl;

  // Optional: overlay balance distributions (MC vs Data) from raw analysis files
  // Uses TH3D photonjet_balance_dist: (photon pT, alpha, balance)
  if (plotBalanceDistOverlay) {
    TFile* fMcRaw = TFile::Open(mcRawFileForDist, "READ");
    TFile* fDtRaw = TFile::Open(dataRawFileForDist, "READ");
    if (!fMcRaw || fMcRaw->IsZombie() || !fDtRaw || fDtRaw->IsZombie()) {
      cout << "WARNING: Could not open raw files for balance_dist overlay; skipping." << endl;
    } else {
      TH3D* hMcDist = (TH3D*)fMcRaw->Get("hibin_-1.0_0.0/eta_-5.2_5.2/photonjet_balance_dist");
      TH3D* hDtDist = (TH3D*)fDtRaw->Get("hibin_-1.0_0.0/eta_-5.2_5.2/photonjet_balance_dist");
      if (!hMcDist || !hDtDist) {
        cout << "WARNING: photonjet_balance_dist not found in raw files; skipping overlay." << endl;
      } else {
        const std::string distFolder = rawFolder + "/balance_dist_overlay";
        gSystem->mkdir(distFolder.c_str(), kTRUE);

        const int yMaxBin = hMcDist->GetYaxis()->FindBin(refAlphaVal - 1e-6);
        const int yMinBin = 1;

        for (int xbin = 1; xbin <= hMcDist->GetXaxis()->GetNbins(); ++xbin) {
          const double ptLo = hMcDist->GetXaxis()->GetBinLowEdge(xbin);
          const double ptHi = hMcDist->GetXaxis()->GetBinLowEdge(xbin + 1);
          // only plot bins overlapping the configured pT fit range
          if (ptHi < ptminG || ptLo > ptmaxG) continue;

          TH1D* hMc = hMcDist->ProjectionZ(Form("hMcDist_ptbin%d", xbin), xbin, xbin, yMinBin, yMaxBin);
          TH1D* hDt = hDtDist->ProjectionZ(Form("hDtDist_ptbin%d", xbin), xbin, xbin, yMinBin, yMaxBin);
          if (!hMc || !hDt) continue;

          if (hMc->Integral() > 0) hMc->Scale(1.0 / hMc->Integral());
          if (hDt->Integral() > 0) hDt->Scale(1.0 / hDt->Integral());

          TCanvas* c = new TCanvas(Form("cBalDist_ptbin%d", xbin), Form("cBalDist_ptbin%d", xbin), 800, 600);
          c->cd();

          hMc->SetLineColor(kBlue);
          hMc->SetLineWidth(2);
          hDt->SetLineColor(kRed);
          hDt->SetLineWidth(2);

          hMc->SetTitle(Form("Balance distribution (normalized), #alpha < %.2f;Balance;1/N dN/dBalance", refAlphaVal));
          hMc->Draw("HIST");
          hDt->Draw("HIST SAME");

          TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
          leg->SetBorderSize(0);
          leg->SetFillStyle(0);
          leg->SetTextFont(42);
          leg->SetTextSize(0.035);
          leg->AddEntry(hMc, "MC", "L");
          leg->AddEntry(hDt, "Data", "L");
          leg->Draw();

          TLatex* t = new TLatex();
          t->SetNDC();
          t->SetTextFont(42);
          t->SetTextSize(0.035);
          t->DrawLatex(0.18, 0.85, Form("%.0f < p_{T}^{#gamma} < %.0f GeV", ptLo, ptHi));
          CMS_lumi(c, 0, 0);

          c->SaveAs(Form("%s/balance_dist_overlay_pt%.0fto%.0f.png", distFolder.c_str(), ptLo, ptHi));

          delete leg;
          delete t;
          delete c;
          delete hMc;
          delete hDt;
        }
      }
      fMcRaw->Close();
      fDtRaw->Close();
    }
  }

  for (auto& alphaPair : ratioVsPt) {
    int alphaBin = alphaPair.first;
    TH1D* hRatio = alphaPair.second;
    if (!hRatio) continue;

    cout << "\nProcessing alpha bin " << alphaBin << endl;
    double alphaCutVal = balance3D_mc->GetZaxis()->GetBinUpEdge(alphaBin);
    cout << "  Alpha cut: alpha < " << alphaCutVal << endl;

    // Ratio plot
    TH1D* hFrame1 = makePtFrameFromHist(hRatio, Form("hFrame1_a%d", alphaBin), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
    TCanvas* c1 = new TCanvas(Form("c1_a%d", alphaBin), Form("c1_a%d", alphaBin), 800, 600);
    c1->SetLogx();
    c1->cd();
    hFrame1->Draw();
    styleLogxAxis(hFrame1);
    const double xMin1 = hFrame1->GetXaxis()->GetBinLowEdge(1);
    const double xMax1 = hFrame1->GetXaxis()->GetBinLowEdge(hFrame1->GetNbinsX() + 1);
    TLine* line = new TLine();
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray+1);
    line->DrawLine(xMin1, 1, xMax1, 1);
    hRatio->SetMarkerStyle(kFullCircle);
    hRatio->SetMarkerColor(kBlue);
    hRatio->SetLineColor(kBlue);
    hRatio->Draw("PE1 SAME");
    TLegend* leg1 = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.035);
    leg1->AddEntry(hRatio, Form("#gamma+jet (#alpha < %.2f)", alphaCutVal), "PLE");
    leg1->Draw();
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.035);
    // L3 is treated as eta-independent here.
    CMS_lumi(c1, 0, 0);
    c1->SaveAs(Form("%s/L3Res_%s_alpha%d_ratio.png", pngFolder.c_str(), _run.c_str(), alphaBin));


    // Optional raw MC/Data responses (no fits, no ratio)
    if (plotRawResponses) {
      TH1D* hMcRaw = (TH1D*)inFile->Get(Form("balance_vsphotonpt_mc_alpha%d", alphaBin));
      TH1D* hDtRaw = (TH1D*)inFile->Get(Form("balance_vsphotonpt_data_alpha%d", alphaBin));
      if (!hMcRaw || !hDtRaw) {
        // Create a placeholder raw plot so users have a consistent output even when raw histos are missing
        TCanvas* cRaw = new TCanvas(Form("cRaw_a%d_missing", alphaBin), Form("cRaw_a%d_missing", alphaBin), 800, 600);
        TH1D* frameRaw = makePtFrameFromHist(hRatio, Form("hFrameRaw_a%d_missing", alphaBin), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
        frameRaw->Draw();
        styleLogxAxis(frameRaw);
        TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.04); t.DrawLatex(0.18, 0.55, "Raw MC/Data responses not found in input");
        CMS_lumi(cRaw, 0, 0);
        cRaw->SaveAs(Form("%s/L3Res_%s_alpha%d_raw_missing.png", rawFolder.c_str(), _run.c_str(), alphaBin));
        delete frameRaw; delete cRaw;
      } else {
        TCanvas* cRaw = new TCanvas(Form("cRaw_a%d", alphaBin), Form("cRaw_a%d", alphaBin), 800, 600);
        cRaw->SetLogx();
        TH1D* frameRaw = makePtFrameFromHist(hMcRaw, Form("hFrameRaw_a%d", alphaBin), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
        frameRaw->Draw();
        styleLogxAxis(frameRaw);
        frameRaw->GetXaxis()->SetNoExponent(kTRUE);
        hMcRaw->SetMarkerStyle(kFullCircle);
        hMcRaw->SetMarkerColor(kBlue);
        hMcRaw->SetLineColor(kBlue);
        hMcRaw->Draw("PE1 SAME");
        hDtRaw->SetMarkerStyle(kFullCircle);
        hDtRaw->SetMarkerColor(kRed);
        hDtRaw->SetLineColor(kRed);
        hDtRaw->Draw("PE1 SAME");
        TLegend* legRaw = new TLegend(0.55, 0.75, 0.88, 0.88);
        legRaw->SetBorderSize(0);
        legRaw->SetFillStyle(0);
        legRaw->SetTextFont(42);
        legRaw->SetTextSize(0.035);
        legRaw->AddEntry(hMcRaw, "MC", "PE");
        legRaw->AddEntry(hDtRaw, "Data", "PE");
        legRaw->Draw();
        TLatex tRaw; tRaw.SetNDC(); tRaw.SetTextFont(42); tRaw.SetTextSize(0.035); tRaw.DrawLatex(0.18, 0.85, Form("#alpha < %.2f", alphaCutVal));
        CMS_lumi(cRaw, 0, 0);
        cRaw->SaveAs(Form("%s/L3Res_%s_alpha%d_raw.png", rawFolder.c_str(), _run.c_str(), alphaBin));
        delete legRaw;
        // Optional diagnostic: per-alpha pT fits (off by default)
        if (plotPerAlphaPtFits) {
          TH1D* hFrame2 = makePtFrameFromHist(hRatio, Form("hFrame2_a%d", alphaBin), "#frac{R_{MC}}{R_{Data}}", 0.7, 1.5);
          TCanvas* c2 = new TCanvas(Form("c2_a%d", alphaBin), Form("c2_a%d", alphaBin), 800, 600);
          c2->SetLogx();
          c2->cd();
          hFrame2->Draw();
          styleLogxAxis(hFrame2);
          const double xMin2 = hFrame2->GetXaxis()->GetBinLowEdge(1);
          const double xMax2 = hFrame2->GetXaxis()->GetBinLowEdge(hFrame2->GetNbinsX() + 1);

          line->DrawLine(xMin2, 1, xMax2, 1);
          line->SetLineStyle(kDotted);
          line->SetLineColor(kGray);
          line->DrawLine(ptminG, 0.5, ptminG, 1.2);
          line->DrawLine(ptmaxG, 0.5, ptmaxG, 1.2);
          line->SetLineStyle(kDashed);
          line->SetLineColor(kGray+1);

          TH1D* hRatioFull = (TH1D*)hRatio->Clone(Form("hRatioFull_a%d", alphaBin));
          hRatioFull->SetMarkerStyle(kOpenCircle);
          hRatioFull->SetMarkerColor(kGray);
          hRatioFull->SetLineColor(kGray);
          hRatioFull->Draw("PE1 SAME");

          TH1D* hRatioClean = drawCleaned(hRatio, "G", ptminG, ptmaxG, kFullCircle, kBlue);
          hRatioClean->Draw("PE1 SAME");

          TMultiGraph* mg = new TMultiGraph(Form("mg_a%d", alphaBin), "mg");
          if (fitG && hRatioClean) {
            TGraphErrors* gG = cleanGraph(new TGraphErrors(hRatioClean));
            if (gG && gG->GetN() > 0) mg->Add(gG, "P");
          }
          const double fitRangeMin = ptminG;
          const double fitRangeMax = ptmaxG;
          TF1* f0 = new TF1(Form("f0_a%d", alphaBin), "[0]", fitRangeMin, fitRangeMax);
          TF1* f1 = new TF1(Form("f1_a%d", alphaBin), l3PtFitExpr(), fitRangeMin, fitRangeMax);
          f0->SetParameter(0, 1.0);
          mg->Fit(f0, "QRN");
          f1->SetParameters(f0->GetParameter(0), -0.01);
          mg->Fit(f1, "QRN");
          f0->SetLineColor(kMagenta+2);
          f0->SetLineStyle(kDotted);
          f0->Draw("SAME");
          f1->SetLineColor(kRed);
          f1->SetLineWidth(2);
          f1->Draw("SAME");

          TLegend* leg2 = new TLegend(0.50, 0.65, 0.88, 0.88);
          leg2->SetBorderSize(0);
          leg2->SetFillStyle(0);
          leg2->SetTextFont(42);
          leg2->SetTextSize(0.030);
          leg2->AddEntry(hRatioClean, "#gamma+jet (cleaned)", "PLE");
          leg2->AddEntry(f0, Form("f_{0}: const (p0=%.4f)", f0->GetParameter(0)), "L");
          leg2->AddEntry(f1, Form("f_{1}: p0=%.4f, p1=%.4f", f1->GetParameter(0), f1->GetParameter(1)), "L");
          leg2->Draw();

          tex->DrawLatex(0.20, 0.85, Form("Fit range: %.0f-%.0f GeV", ptminG, ptmaxG));
          tex->DrawLatex(0.20, 0.80, Form("(#alpha < %.2f), bin %d", alphaCutVal, alphaBin));
          tex->DrawLatex(0.20, 0.75, Form("|#eta_{jet}| < %.1f", balance3D_mc->GetYaxis()->GetBinLowEdge(netabins + 1)));
          CMS_lumi(c2, 0, 0);
          c2->SaveAs(Form("%s/L3Res_%s_alpha%d_fits.png", pngFolder.c_str(), _run.c_str(), alphaBin));

          // Closure plot for this alpha bin (optional)
          if (doClosure) {
            TH1D* hFrame3 = makePtFrameFromHist(hRatio, Form("hFrame3_a%d", alphaBin), "Corrected balance", 0.7, 1.3);
            TCanvas* c3 = new TCanvas(Form("c3_a%d", alphaBin), Form("c3_a%d", alphaBin), 800, 600);
            c3->SetLogx();
            c3->cd();
            hFrame3->Draw();
            styleLogxAxis(hFrame3);
            const double xMin3 = hFrame3->GetXaxis()->GetBinLowEdge(1);
            const double xMax3 = hFrame3->GetXaxis()->GetBinLowEdge(hFrame3->GetNbinsX() + 1);
            line->SetLineStyle(kDashed);
            line->SetLineColor(kGray+1);
            line->DrawLine(xMin3, 1, xMax3, 1);

            TH1D* hCorrected = (TH1D*)hRatio->Clone(Form("hCorrected_a%d", alphaBin));
            for (int i = 1; i <= hCorrected->GetNbinsX(); ++i) {
              const double pt = hCorrected->GetBinCenter(i);
              const double val = hCorrected->GetBinContent(i);
              const double err = hCorrected->GetBinError(i);
              if (val > 0 && pt >= fitRangeMin && pt <= fitRangeMax) {
                const double corr = f1->Eval(pt);
                hCorrected->SetBinContent(i, val / corr);
                hCorrected->SetBinError(i, err / corr);
              }
            }
            hCorrected->SetMarkerStyle(kFullCircle);
            hCorrected->SetMarkerColor(kBlue);
            hCorrected->SetLineColor(kBlue);
            hCorrected->Draw("PE1 SAME");
            TLegend* leg3 = new TLegend(0.50, 0.75, 0.88, 0.88);
            leg3->SetBorderSize(0);
            leg3->SetFillStyle(0);
            leg3->SetTextFont(42);
            leg3->SetTextSize(0.030);
            leg3->AddEntry(hCorrected, "Corrected", "PLE");
            leg3->AddEntry(f1, "Fit applied", "L");
            leg3->Draw();
            tex->DrawLatex(0.20, 0.85, Form("(#alpha < %.2f), bin %d", alphaCutVal, alphaBin));
            tex->DrawLatex(0.20, 0.80, Form("Fit range: %.0f-%.0f GeV", fitRangeMin, fitRangeMax));
            CMS_lumi(c3, 0, 0);
            c3->SaveAs(Form("%s/L3Res_%s_alpha%d_closure.png", pngFolder.c_str(), _run.c_str(), alphaBin));
            delete leg3;
            delete hCorrected;
            delete hFrame3;
            delete c3;
          }

          outfile->cd();
          if (hRatioClean) hRatioClean->Write(Form("ratio_vspT_cleaned_alpha%d", alphaBin));
          f0->Write();
          f1->Write();

          delete leg2;
          delete f0;
          delete f1;
          delete mg;
          delete hRatioFull;
          delete hRatioClean;
          delete hFrame2;
          delete c2;
        }
      }
    } // end if (plotRawResponses)


    // Eta-dependent maps (counts and l3resmap) are optional.
    // L3 is treated as eta-independent by default.
    if (plotEtaMaps) {

    const double etaMinAll = balance3D_mc->GetYaxis()->GetBinLowEdge(1);
    const double etaMaxAll = balance3D_mc->GetYaxis()->GetBinLowEdge(netabins + 1);

    // Counts MC
    if (counts3D_mc) {
      counts3D_mc->GetZaxis()->SetRange(1, alphaBin);
      TH2D* h2dCounts = (TH2D*)counts3D_mc->Project3D("yx");
      h2dCounts->SetName(Form("counts_pteta_mc_alpha%d", alphaBin));
      double ptMinAll = balance3D_mc->GetXaxis()->GetBinLowEdge(1);
      double ptMaxAll = balance3D_mc->GetXaxis()->GetBinLowEdge(nptbins+1);
      h2dCounts->SetTitle(Form("MC Counts (#alpha < %.2f);p_{T}^{#gamma} (GeV) [%.0f-%.0f];|#eta_{jet}| [%.3f-%.3f]", alphaCutVal, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll));
      TCanvas* cCounts = new TCanvas(Form("cCounts_mc_a%d", alphaBin), Form("cCounts_mc_a%d", alphaBin), 1000, 700);
      cCounts->SetRightMargin(0.15);
      cCounts->SetLeftMargin(0.12);
      cCounts->SetLogx();
      cCounts->cd();
      // draw as color map (no overlaid text) and use log-x for pT readability
      gStyle->SetPaintTextFormat("0.0f");
      h2dCounts->SetMarkerSize(1.4);
      // Improve log-x axis labeling for pT range on 2D maps
      h2dCounts->GetXaxis()->SetMoreLogLabels(kTRUE);
      h2dCounts->GetXaxis()->SetNdivisions(510);
      h2dCounts->GetXaxis()->SetNoExponent(kTRUE);
      h2dCounts->Draw("TEXTCOLZ");
      CMS_lumi(cCounts, 0, 0);
      cCounts->SaveAs(Form("%s/L3Res_%s_alpha%d_counts_mc_pt%.0fto%.0f_eta%.3fto%.3f.png", pngFolder.c_str(), _run.c_str(), alphaBin, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll));
      outfile->cd();
      h2dCounts->Write();
      delete cCounts;
      counts3D_mc->GetZaxis()->SetRange(1, nalphabins);
    }

    // Counts Data
    if (counts3D_data) {
      counts3D_data->GetZaxis()->SetRange(1, alphaBin);
      TH2D* h2dCountsData = (TH2D*)counts3D_data->Project3D("yx");
      h2dCountsData->SetName(Form("counts_pteta_data_alpha%d", alphaBin));
      double ptMinAll = balance3D_data->GetXaxis()->GetBinLowEdge(1);
      double ptMaxAll = balance3D_data->GetXaxis()->GetBinLowEdge(nptbins+1);
      h2dCountsData->SetTitle(Form("Data Counts (#alpha < %.2f);p_{T}^{#gamma} (GeV) [%.0f-%.0f];|#eta_{jet}| [%.3f-%.3f]", alphaCutVal, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll));
      TCanvas* cCountsData = new TCanvas(Form("cCounts_data_a%d", alphaBin), Form("cCounts_data_a%d", alphaBin), 1000, 700);
      cCountsData->SetRightMargin(0.15);
      cCountsData->SetLeftMargin(0.12);
      cCountsData->SetLogx();
      cCountsData->cd();
      gStyle->SetPaintTextFormat("0.0f");
      h2dCountsData->SetMarkerSize(1.4);
      // Improve log-x axis labeling for pT range on 2D maps
      h2dCountsData->GetXaxis()->SetMoreLogLabels(kTRUE);
      h2dCountsData->GetXaxis()->SetNdivisions(510);
      h2dCountsData->GetXaxis()->SetNoExponent(kTRUE);
      h2dCountsData->Draw("TEXTCOLZ");
      CMS_lumi(cCountsData, 0, 0);
      cCountsData->SaveAs(Form("%s/L3Res_%s_alpha%d_counts_data_pt%.0fto%.0f_eta%.3fto%.3f.png", pngFolder.c_str(), _run.c_str(), alphaBin, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll));
      outfile->cd();
      h2dCountsData->Write();
      delete cCountsData;
      counts3D_data->GetZaxis()->SetRange(1, nalphabins);
    }

    // L3 Residual Map (MC/Data ratio) as 2D TEXTCOLZ plot
    {
      double ptMinAll = balance3D_mc->GetXaxis()->GetBinLowEdge(1);
      double ptMaxAll = balance3D_mc->GetXaxis()->GetBinLowEdge(nptbins+1);
      
      // Create 2D histogram for L3Res map
      std::vector<double> ptEdges(nptbins + 1), etaEdges(netabins + 1);
      for (int b = 1; b <= nptbins + 1; ++b) ptEdges[b-1] = balance3D_mc->GetXaxis()->GetBinLowEdge(b);
      for (int b = 1; b <= netabins + 1; ++b) etaEdges[b-1] = balance3D_mc->GetYaxis()->GetBinLowEdge(b);
      
      TH2D* h2dL3Res = new TH2D(Form("l3res_pteta_alpha%d", alphaBin),
                                 Form("L3 Residual (MC/Data) (#alpha < %.2f);p_{T}^{#gamma} (GeV) [%.0f-%.0f];|#eta_{jet}| [%.3f-%.3f]",
                                      alphaCutVal, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll),
                                 nptbins, ptEdges.data(), netabins, etaEdges.data());
      
      // Fill from 3D profiles with cumulative alpha cut
      balance3D_mc->GetZaxis()->SetRange(1, alphaBin);
      balance3D_data->GetZaxis()->SetRange(1, alphaBin);
      
      for (int ptbin = 1; ptbin <= nptbins; ++ptbin) {
        for (int etabin = 1; etabin <= netabins; ++etabin) {
          double sum_mc = 0., entries_mc = 0.;
          double sum_dt = 0., entries_dt = 0.;
          
          for (int abin = 1; abin <= alphaBin; ++abin) {
            double val_mc = balance3D_mc->GetBinContent(ptbin, etabin, abin);
            double ent_mc = balance3D_mc->GetBinEntries(balance3D_mc->GetBin(ptbin, etabin, abin));
            double val_dt = balance3D_data->GetBinContent(ptbin, etabin, abin);
            double ent_dt = balance3D_data->GetBinEntries(balance3D_data->GetBin(ptbin, etabin, abin));
            
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
          double l3res = (dt_val > 0 && mc_val > 0) ? mc_val / dt_val : 1.0;
          
          h2dL3Res->SetBinContent(ptbin, etabin, l3res);
        }
      }
      
      balance3D_mc->GetZaxis()->SetRange(1, nalphabins);
      balance3D_data->GetZaxis()->SetRange(1, nalphabins);
      
      TCanvas* cL3Res = new TCanvas(Form("cL3Res_a%d", alphaBin), Form("cL3Res_a%d", alphaBin), 1000, 700);
      cL3Res->SetRightMargin(0.15);
      cL3Res->SetLeftMargin(0.12);
      cL3Res->SetLogx();
      cL3Res->cd();
      gStyle->SetPaintTextFormat("0.3f");
      // Improve log-x axis labeling for pT range on 2D maps
      h2dL3Res->GetXaxis()->SetMoreLogLabels(kTRUE);
      h2dL3Res->GetXaxis()->SetNdivisions(510);
      h2dL3Res->GetXaxis()->SetNoExponent(kTRUE);
      h2dL3Res->SetMarkerSize(1.4);
      h2dL3Res->SetMinimum(0.7);
      h2dL3Res->SetMaximum(1.5);
      h2dL3Res->Draw("TEXTCOLZ");
      CMS_lumi(cL3Res, 0, 0);
      cL3Res->SaveAs(Form("%s/L3Res_%s_alpha%d_l3resmap_pt%.0fto%.0f_eta%.3fto%.3f.png", pngFolder.c_str(), _run.c_str(), alphaBin, ptMinAll, ptMaxAll, etaMinAll, etaMaxAll));
      outfile->cd();
      h2dL3Res->Write();
      delete cL3Res;
    }

    } // plotEtaMaps

    outfile->cd();
    hRatio->Write(Form("ratio_vspT_alpha%d", alphaBin));
    delete c1;
  }

  outfile->Close();
  inFile->Close();

  cout << "\n============================================" << endl;
  cout << "Output written to: " << outfolder << "/" << outfilename << "_fit.root" << endl;
  cout << "Plots saved in: " << pngFolder << "/" << endl;
  cout << "Final L3Residual text file: " << txtFolder << "/" << outfilename << ".txt" << endl;
  cout << "============================================" << endl;
}
