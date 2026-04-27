// Create .txt file with the L2res corrections

#include "../fillhistograms/histograms.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

int nparams = 5;

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

struct CorrBin {
  double lowPt;
  double highPt;
  TH1D *hist;
};
} // namespace

void doTxt(TString inFileName =
               "../RERECORESULTS/kfactor_rerunall_combined_allpts.root",
           TString outFileName = "L2residual.txt") {

  TFile *inFile = new TFile(inFileName, "READ");
  // TFile *inFile2 = new TFile("", "READ");

  string header =
      "{2 JetEta JetPt 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}";

  ofstream txtfile;
  txtfile.open(outFileName);

  std::cout << "Producing txt file " << outFileName << " with header "
            << header.c_str() << std::endl;
  txtfile << header.c_str() << std::endl;

  std::vector<CorrBin> corrBins;
  for (int ptbin = 1; ptbin <= histograms::nptforjec; ++ptbin) {
    const std::string lbl = ptBinLabel(ptbin);
    TH1D *h =
        dynamic_cast<TH1D *>(inFile->Get(Form("corrections_%s", lbl.c_str())));
    if (!h)
      continue;
    corrBins.push_back(
        {histograms::ptforjec[ptbin - 1], histograms::ptforjec[ptbin], h});
  }

  if (corrBins.empty()) {
    std::cerr << "No corrections_* histograms found in " << inFileName
              << std::endl;
    txtfile.close();
    return;
  }

  std::sort(
      corrBins.begin(), corrBins.end(),
      [](const CorrBin &a, const CorrBin &b) { return a.lowPt < b.lowPt; });

  // Alternatively: take the last bin of pt from a different file
  /*   for (int ptbin = 2; ptbin < 4; ++ptbin) {
    facts[ptbin] =
  (TH1D*)inFile->Get(Form("corrs_%dto%d",pts[ptbin-1],pts[ptbin]));
  }
  facts[4] = (TH1D*)inFile2->Get(Form("corrs_%dto%d",pts[4],pts[5])); */

  const int nEtaBins = corrBins.front().hist->GetXaxis()->GetNbins();

  for (int etabin = nEtaBins; etabin >= 1; --etabin) {
    for (const auto &cb : corrBins) {
      auto factors = cb.hist;

      // For case of different eta bin widths
      float loweta = factors->GetXaxis()->GetBinLowEdge(etabin + 1);
      float higheta = factors->GetXaxis()->GetBinLowEdge(etabin);

      if (factors->GetXaxis()->GetBinLowEdge(etabin) < 2.9) {
        float lowpt = 0.001;
        if (cb.lowPt > 25)
          lowpt = cb.lowPt;

        if (higheta != 0)
          txtfile << -loweta << " " << -higheta << " " << lowpt << " "
                  << cb.highPt << " " << nparams << " " << cb.lowPt << " "
                  << cb.highPt << " "
                  << factors->GetBinContent(factors->FindBin(higheta + 0.001))
                  << " 0 1" << std::endl;

        else
          txtfile << -loweta << " " << higheta << " " << lowpt << " "
                  << cb.highPt << " " << nparams << " " << cb.lowPt << " "
                  << cb.highPt << " "
                  << factors->GetBinContent(factors->FindBin(higheta + 0.001))
                  << " 0 1" << std::endl;
      }
    }
  }

  for (int etabin = 1; etabin <= nEtaBins;
       ++etabin) { // These are the narrow eta bins
    for (const auto &cb : corrBins) {
      auto factors = cb.hist;

      float loweta = factors->GetXaxis()->GetBinLowEdge(etabin);
      float higheta = factors->GetXaxis()->GetBinLowEdge(etabin + 1);

      if (factors->GetXaxis()->GetBinLowEdge(etabin) < 2.9) {
        float lowpt = 0.001;
        if (cb.lowPt > 25)
          lowpt = cb.lowPt;

        txtfile << loweta << " " << higheta << " " << lowpt << " " << cb.highPt
                << " " << nparams << " " << cb.lowPt << " " << cb.highPt << " "
                << factors->GetBinContent(factors->FindBin(loweta + 0.001))
                << " 0 1" << std::endl;
      }
    }
  }

  txtfile.close();
}
