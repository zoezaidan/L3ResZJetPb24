// Purpose: create the final L3 and L2L3 text payloads from the active L3 fit.
//
// Step 1. Read the fit ROOT file produced by L3Res.C.
// Step 2. Refit the shared pTref graph directly with the full Run3-style L3
//         correction model, without converting the fit to JetPt.
// Step 3. Read the requested L2Residual text file to define the row layout
//         and the existing L2 correction term.
// Step 4. Clone each L2Residual row and append the same exported L3
//         function to produce the standalone L3 and combined L2L3 text files.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "CMS_lumi.C"
#include "L3ResCommon.h"
#include "tdrStyle.C"

using namespace std;

static void useExportPlotStyle(const string &lumiLabel) {
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  writeExtraText = true;
  extraText = "Preliminary";
  cmsTextSize = 0.78;
  lumi_sqrtS = Form("%s, #sqrt{s} = 5.36 TeV", lumiLabel.c_str());
}

static void drawCmsStamp(TCanvas &canvas) {
  canvas.cd();
  CMS_lumi(&canvas, 0, 0);
  canvas.Modified();
  canvas.Update();
}

static void styleExportFrame(TH1D *frame, double yTitleOffset = 1.35) {
  if (!frame)
    return;
  frame->GetXaxis()->SetTitleSize(0.046);
  frame->GetYaxis()->SetTitleSize(0.046);
  frame->GetXaxis()->SetLabelSize(0.034);
  frame->GetYaxis()->SetLabelSize(0.034);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(yTitleOffset);
}

struct JecRecord {
  double etaMin = 0.0;
  double etaMax = 0.0;
  int nPar = 0;
  double ptMin = 0.0;
  double ptMax = 0.0;
  vector<double> parameters;
};

static string trimCopy(const string &text) {
  size_t begin = 0;
  while (begin < text.size() &&
         isspace(static_cast<unsigned char>(text[begin])))
    ++begin;
  size_t end = text.size();
  while (end > begin && isspace(static_cast<unsigned char>(text[end - 1])))
    --end;
  return text.substr(begin, end - begin);
}

static bool digitsOnly(const string &text) {
  if (text.empty())
    return false;
  for (char ch : text) {
    if (!isdigit(static_cast<unsigned char>(ch)))
      return false;
  }
  return true;
}

static string rewriteParameterSlots(const string &expression, int offset) {
  string rewritten;
  rewritten.reserve(expression.size() + 16);

  for (size_t pos = 0; pos < expression.size();) {
    if (expression[pos] != '[') {
      rewritten.push_back(expression[pos]);
      ++pos;
      continue;
    }

    const size_t closePos = expression.find(']', pos + 1);
    if (closePos == string::npos) {
      rewritten.append(expression.substr(pos));
      break;
    }

    const string slot =
        trimCopy(expression.substr(pos + 1, closePos - pos - 1));
    string digits = slot;
    if (!slot.empty() && (slot[0] == 'p' || slot[0] == 'P'))
      digits = slot.substr(1);

    if (digitsOnly(digits)) {
      rewritten += "[" + to_string(stoi(digits) + offset) + "]";
    } else {
      rewritten.append(expression.substr(pos, closePos - pos + 1));
    }
    pos = closePos + 1;
  }

  return rewritten;
}

static void drawExportModelBlock(double x, double y, const TF1 *fit) {
  if (!fit)
    return;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.024);
  latex.DrawLatex(x, y,
                  "C_{L3}(x)=1./(p_{0}+p_{1}/x+p_{2}log(x)/"
                  "x+p_{3}T(x)+p_{6}x^{-0.3051}+p_{7}x)");
  latex.DrawLatex(x, y - 0.040,
                  "T(x)=((x/p_{4})^{p_{5}}-1)/((x/p_{4})^{p_{5}}+1)");
  latex.DrawLatex(x, y - 0.085,
                  Form("p_{0}=%.4f, p_{1}=%.4f, p_{2}=%.4f, p_{3}=%.4f",
                       fit->GetParameter(0), fit->GetParameter(1),
                       fit->GetParameter(2), fit->GetParameter(3)));
  latex.DrawLatex(x, y - 0.125,
                  Form("p_{4}=%.4f, p_{5}=%.4f, p_{6}=%.4f, p_{7}=%.6f",
                       fit->GetParameter(4), fit->GetParameter(5),
                       fit->GetParameter(6), fit->GetParameter(7)));
}

static string extractSimpleJecFormula(const string &headerLine) {
  const size_t correctionPos = headerLine.find(" Correction ");
  const size_t jetPtPos = headerLine.find(" JetPt ");
  if (correctionPos == string::npos || jetPtPos == string::npos)
    return "";
  const size_t expressionStart = jetPtPos + string(" JetPt ").size();
  if (expressionStart >= correctionPos)
    return "";
  return trimCopy(
      headerLine.substr(expressionStart, correctionPos - expressionStart));
}

static bool readJecFile(const string &path, string &headerLine,
                        vector<JecRecord> &records) {
  ifstream input(path.c_str());
  if (!input.is_open())
    return false;

  headerLine.clear();
  records.clear();

  string line;
  while (getline(input, line)) {
    const string trimmed = trimCopy(line);
    if (trimmed.empty() || trimmed[0] == '#')
      continue;
    headerLine = trimmed;
    break;
  }
  if (headerLine.empty())
    return false;

  while (getline(input, line)) {
    const string trimmed = trimCopy(line);
    if (trimmed.empty() || trimmed[0] == '#')
      continue;

    istringstream parser(trimmed);
    JecRecord record;
    parser >> record.etaMin >> record.etaMax >> record.nPar >> record.ptMin >>
        record.ptMax;
    if (!parser || record.nPar < 2)
      continue;

    record.parameters.resize(record.nPar - 2, 0.0);
    for (int index = 0; index < record.nPar - 2; ++index)
      parser >> record.parameters[index];
    if (parser)
      records.push_back(record);
  }

  return !records.empty();
}

static pair<double, double> graphXBounds(const TGraphErrors *graph) {
  if (!graph || graph->GetN() <= 0)
    return make_pair(0.0, 0.0);

  double xmin = graph->GetX()[0];
  double xmax = graph->GetX()[0];
  for (int point = 1; point < graph->GetN(); ++point) {
    xmin = std::min(xmin, graph->GetX()[point]);
    xmax = std::max(xmax, graph->GetX()[point]);
  }
  return make_pair(xmin, xmax);
}

static pair<double, double> graphYBounds(const TGraphErrors *graph) {
  if (!graph || graph->GetN() <= 0)
    return make_pair(0.7, 1.3);

  double ymin = graph->GetY()[0];
  double ymax = graph->GetY()[0];
  for (int point = 1; point < graph->GetN(); ++point) {
    ymin = std::min(ymin, graph->GetY()[point]);
    ymax = std::max(ymax, graph->GetY()[point]);
  }
  return make_pair(ymin, ymax);
}

static const char *l3ExportFitExpr() {
  return "1./([0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/"
         "[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x)";
}

static void initL3ExportFitParams(TF1 *fit, double ySeed = 1.0) {
  if (!fit)
    return;

  const double safeSeed = (ySeed > 0.0 ? ySeed : 1.0);
  fit->SetParameters(0.9722 / safeSeed, 0.7944, 2.14069, 0.10229, 10.52, 1.5550,
                     -0.72222, -2.25e-06);
  fit->SetParLimits(0, 0.5, 2.0);
  fit->SetParLimits(4, 1.0, 100.0);
  fit->SetParLimits(5, 0.1, 5.0);
}

static TF1 *cloneFullExportFit(const TF1 *sourceFit, const char *name,
                               double fitMin, double fitMax) {
  if (!sourceFit || sourceFit->GetNpar() != 8)
    return nullptr;

  TF1 *clone = new TF1(name, l3ExportFitExpr(), fitMin, fitMax);
  for (int parameter = 0; parameter < clone->GetNpar(); ++parameter) {
    clone->SetParameter(parameter, sourceFit->GetParameter(parameter));
  }
  return clone;
}

static string formatL3Record(double etaMin, double etaMax, double ptMin,
                             double ptMax, const TF1 *fit) {
  ostringstream output;
  const int nParameters = fit ? fit->GetNpar() : 0;
  output << Form("  %6.3f %6.3f  %d  %5.0f %5.0f", etaMin, etaMax,
                 nParameters + 2, ptMin, ptMax);
  for (int parameter = 0; parameter < nParameters; ++parameter)
    output << Form("  %10.6f", fit->GetParameter(parameter));
  output << "\n";
  return output.str();
}

static string formatCombinedRecord(const JecRecord &record, double ptMin,
                                   double ptMax, const TF1 *fit) {
  ostringstream output;
  const int nParameters = fit ? fit->GetNpar() : 0;
  output << Form("  %6.3f %6.3f  %d  %5.0f %5.0f", record.etaMin, record.etaMax,
                 record.nPar + nParameters, ptMin, ptMax);
  for (size_t index = 0; index < record.parameters.size(); ++index)
    output << Form("  %10.6f", record.parameters[index]);
  for (int parameter = 0; parameter < nParameters; ++parameter)
    output << Form("  %10.6f", fit->GetParameter(parameter));
  output << "\n";
  return output.str();
}

void createL2L3ResTextFile(
    TString fitRootFile = "L3Residual/L3Res_photonjet/L3Res_photonjet_fit.root",
    string l2ResidualFile =
        "fillhistograms/jecfiles/Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt") {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *fitFile = TFile::Open(fitRootFile, "READ");
  if (!fitFile || fitFile->IsZombie()) {
    cout << "ERROR: Cannot open fit ROOT file '" << fitRootFile << "'." << endl;
    return;
  }

  TGraphErrors *combinedGraph =
      (TGraphErrors *)fitFile->Get("g_ptref_combined");
  TF1 *ptRefFit = (TF1 *)fitFile->Get("f_ptref_final");
  TNamed *runLabelObj = (TNamed *)fitFile->Get("run_label");
  TNamed *outputNameObj = (TNamed *)fitFile->Get("output_name");
  TNamed *outputTagObj = (TNamed *)fitFile->Get("output_tag");
  TNamed *lumiLabelObj = (TNamed *)fitFile->Get("lumi_label");
  TParameter<double> *fitMinObj =
      (TParameter<double> *)fitFile->Get("pt_fit_min");
  TParameter<double> *fitMaxObj =
      (TParameter<double> *)fitFile->Get("pt_fit_max");

  if (!combinedGraph || !ptRefFit || !runLabelObj || !outputNameObj ||
      !outputTagObj || !fitMinObj || !fitMaxObj) {
    cout << "ERROR: Fit file '" << fitRootFile
         << "' is missing the required L3Res.C metadata." << endl;
    fitFile->Close();
    delete fitFile;
    return;
  }

  const string outfolder = gSystem->DirName(fitRootFile.Data());
  const string pdfFolder = outfolder + "/pdf";
  const string txtFolder = outfolder + "/textfiles";
  gSystem->mkdir(pdfFolder.c_str(), kTRUE);
  gSystem->mkdir(txtFolder.c_str(), kTRUE);

  const string outfilename = outputNameObj->GetTitle();
  const string runLabel = runLabelObj->GetTitle();
  const string outTag = outputTagObj->GetTitle();
  const string lumiLabel =
      lumiLabelObj ? lumiLabelObj->GetTitle() : "pp 480.4 pb^{-1}";
  const double fitMin = fitMinObj->GetVal();
  const double fitMax = fitMaxObj->GetVal();

  useExportPlotStyle(lumiLabel);

  string l2Header;
  vector<JecRecord> l2Records;
  if (!readJecFile(l2ResidualFile, l2Header, l2Records)) {
    cout << "ERROR: Cannot read the requested L2Residual text file '"
         << l2ResidualFile << "'." << endl;
    fitFile->Close();
    delete fitFile;
    return;
  }

  const string l2CorrectionExpr = extractSimpleJecFormula(l2Header);
  if (l2CorrectionExpr.empty()) {
    cout << "ERROR: Could not parse the L2Residual header expression from '"
         << l2ResidualFile << "'." << endl;
    fitFile->Close();
    delete fitFile;
    return;
  }

  TFile *outputRoot = new TFile(
      Form("%s/%s_textfits.root", outfolder.c_str(), outfilename.c_str()),
      "RECREATE");

  TF1 *exportPtRefFit =
      cloneFullExportFit(ptRefFit, "f_ptref_export", fitMin, fitMax);
  if (!exportPtRefFit) {
    exportPtRefFit =
        new TF1("f_ptref_export", l3ExportFitExpr(), fitMin, fitMax);
    initL3ExportFitParams(exportPtRefFit, graphMeanY(combinedGraph));
    int fitStatus = (int)combinedGraph->Fit(exportPtRefFit, "QRN");
    fitStatus = (int)combinedGraph->Fit(exportPtRefFit, "QRN");
    if (fitStatus != 0) {
      cout << "WARNING: The direct pTref export fit returned status "
           << fitStatus << "." << endl;
    }
  }

  TCanvas ptRefCanvas("cPtRefExport", "cPtRefExport", 900, 700);
  ptRefCanvas.SetLogx();
  const pair<double, double> yBounds = graphYBounds(combinedGraph);
  TH1D *ptRefFrame =
      new TH1D("hPtRefExportFrame",
               ";p_{T}^{ref} (GeV);L3 residual factor = B^{Data}/B^{MC}", 100,
               fitMin, fitMax);
  ptRefFrame->SetMinimum(std::max(0.45, yBounds.first - 0.08));
  ptRefFrame->SetMaximum(std::min(1.8, yBounds.second + 0.30));
  ptRefFrame->GetXaxis()->SetMoreLogLabels(kTRUE);
  ptRefFrame->GetXaxis()->SetNoExponent(kTRUE);
  styleExportFrame(ptRefFrame);
  ptRefFrame->Draw();

  TLine unityPtRef(fitMin, 1.0, fitMax, 1.0);
  unityPtRef.SetLineStyle(kDashed);
  unityPtRef.SetLineColor(kGray + 1);
  unityPtRef.Draw("SAME");

  combinedGraph->SetMarkerStyle(kFullCircle);
  combinedGraph->SetMarkerColor(kBlack);
  combinedGraph->SetLineColor(kBlack);
  combinedGraph->Draw("PE SAME");

  exportPtRefFit->SetLineColor(kRed + 1);
  exportPtRefFit->SetLineWidth(2);
  exportPtRefFit->Draw("SAME");

  TLatex ptRefLabel;
  ptRefLabel.SetNDC();
  ptRefLabel.SetTextFont(42);
  ptRefLabel.SetTextSize(0.032);
  ptRefLabel.DrawLatex(0.18, 0.86, Form("Run label: %s", runLabel.c_str()));
  ptRefLabel.DrawLatex(0.18, 0.81,
                       Form("Fit range: %.0f-%.0f GeV", fitMin, fitMax));
  ptRefLabel.DrawLatex(0.18, 0.76,
                       "Direct global p_{T}^{ref} fit to R_{Data/MC}");
  ptRefLabel.DrawLatex(0.18, 0.71,
                       Form("#chi^{2}/ndf = %.1f/%d",
                            exportPtRefFit->GetChisquare(),
                            exportPtRefFit->GetNDF()));
  ptRefLabel.DrawLatex(0.18, 0.66, "Full L3 model");
  drawExportModelBlock(0.18, 0.61, exportPtRefFit);
  drawCmsStamp(ptRefCanvas);
  ptRefCanvas.SaveAs(Form("%s/L3Res_%s_ptref_export_fit.png", pdfFolder.c_str(),
                          runLabel.c_str()));

  outputRoot->cd();
  combinedGraph->Write("g_ptref_combined");
  exportPtRefFit->Write("f_ptref_export");

  const int l2ParameterCount = (int)l2Records.front().parameters.size();
  const string l3Header = string("{ 1 JetEta 1 JetPt ") +
                          rewriteParameterSlots(l3ExportFitExpr(), 0) +
                          " Correction L3Residual}";
  const string combinedHeader =
      string("{ 1 JetEta 1 JetPt ") +
      rewriteParameterSlots(l2CorrectionExpr, 0) + "*" +
      rewriteParameterSlots(l3ExportFitExpr(), l2ParameterCount) +
      " Correction L2L3Residual}";

  ofstream localL3File(
      Form("%s/%s.txt", txtFolder.c_str(), outfilename.c_str()));
  ofstream exportL3File(Form("%s/L3Residuals_%s_%s_AK4PF.txt",
                             txtFolder.c_str(), runLabel.c_str(),
                             outTag.c_str()));
  ofstream exportCombinedFile(Form("%s/L2L3Residuals_%s_%s_AK4PF.txt",
                                   txtFolder.c_str(), runLabel.c_str(),
                                   outTag.c_str()));
  localL3File << l3Header << "\n";
  exportL3File << l3Header << "\n";
  exportCombinedFile << combinedHeader << "\n";

  for (const JecRecord &record : l2Records) {
    const double validPtMin = record.ptMin;
    const double validPtMax = record.ptMax;
    const string l3Record = formatL3Record(
        record.etaMin, record.etaMax, validPtMin, validPtMax, exportPtRefFit);
    localL3File << l3Record;
    exportL3File << l3Record;
    exportCombinedFile << formatCombinedRecord(record, validPtMin, validPtMax,
                                               exportPtRefFit);
  }

  localL3File.close();
  exportL3File.close();
  exportCombinedFile.close();

  delete ptRefFrame;

  outputRoot->Close();
  fitFile->Close();

  cout << "============================================" << endl;
  cout << "Wrote L3 text file: " << txtFolder << "/" << outfilename << ".txt"
       << endl;
  cout << "Wrote export L3 text file: " << txtFolder << "/L3Residuals_"
       << runLabel << "_" << outTag << "_AK4PF.txt" << endl;
  cout << "Wrote combined L2L3 text file: " << txtFolder << "/L2L3Residuals_"
       << runLabel << "_" << outTag << "_AK4PF.txt" << endl;
  cout << "Used the full direct pTref L3 model for export; no JetPt refit was "
          "applied."
       << endl;
  cout << "Export plots saved in: " << pdfFolder << endl;
  cout << "============================================" << endl;
}
