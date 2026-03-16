// Create .txt file with the L2res corrections

int nparams = 5;

void doTxt(TString inFileName = "../RERECORESULTS/kfactor_rerunall_combined_allpts.root", TString outFileName = "L2residual.txt" ) {

  // int pts[] = {40, 55, 80, 120, 170, 1000}; // Binning in files
  int pts[] = {15, 25, 80, 120, 1000}; // Binning in files
 
  // TFile *inFile = new TFile("../fits_combinedbins/kfactor_rerunall_combined_allpts.root", "READ");
   TFile *inFile = new TFile(inFileName, "READ");
  // TFile *inFile2 = new TFile("", "READ");
  
  string header = "{2 JetEta JetPt 1 JetPt ([0]+[1]*log(x))*[2] Correction L2Relative}";
   
  ofstream txtfile;
  txtfile.open(outFileName);

  std::cout << "Producing txt file " << outFileName << " with header " << header.c_str() << std::endl;
  txtfile << header.c_str() << std::endl;

  map<int, TH1D*> facts;
   for (int ptbin = 2; ptbin <= 4; ++ptbin) {
     //       facts[ptbin] = (TH1D*)inFile->Get(Form("corrs_%dto%d",pts[ptbin-1],pts[ptbin]));
     facts[ptbin] = (TH1D*)inFile->Get(Form("corrections_%dto%d",pts[ptbin-1],pts[ptbin]));
     }

  // Alternatively: take the last bin of pt from a different file
  /*   for (int ptbin = 2; ptbin < 4; ++ptbin) {
    facts[ptbin] = (TH1D*)inFile->Get(Form("corrs_%dto%d",pts[ptbin-1],pts[ptbin]));
  }
  facts[4] = (TH1D*)inFile2->Get(Form("corrs_%dto%d",pts[4],pts[5])); */

  for (int etabin = facts[3]->GetXaxis()->GetNbins(); etabin >= 1; --etabin) {
        for (int ptbin = 2; ptbin <= 4; ++ptbin) {
	  auto factors = facts[ptbin];

	  // For case of different eta bin widths
	  float loweta = factors->GetXaxis()->GetBinLowEdge(etabin+1);
	  float higheta = factors->GetXaxis()->GetBinLowEdge(etabin);

	  if (factors->GetXaxis()->GetBinLowEdge(etabin) < 2.9) {
	    if (ptbin == 4 and loweta > 2.5) break;    // pT 120- range stops at |eta| = 2.5
	    //    if (ptbin == 5 and loweta > 2.0) break;    // pT 170-1000 range stops at |eta| = 1.93

	    float lowpt = 0.001;
	    if (pts[ptbin-1] > 25) lowpt = pts[ptbin-1]; 
	 
	    if (higheta != 0) txtfile << -loweta << " " << -higheta << " " << lowpt << " " << pts[ptbin] << " " << nparams << " " << pts[ptbin-1] << " " << pts[ptbin] << " " << factors->GetBinContent(factors->FindBin(higheta+0.001)) << " 0 1"  << std::endl;
	    
	    else txtfile << -loweta << " " << higheta << " " << lowpt << " " << pts[ptbin] << " " << nparams<< " " << pts[ptbin-1] << " " << pts[ptbin] << " " << factors->GetBinContent(factors->FindBin(higheta+0.001)) << " 0 1"  << std::endl;
	 
	  }     
      }
  }
  
  for (int etabin = 1; etabin <= facts[3]->GetXaxis()->GetNbins(); ++etabin) {   // These are the narrow eta bins
         for (int ptbin = 2; ptbin <= 4; ++ptbin) {
	  auto factors = facts[ptbin];

	  float loweta = factors->GetXaxis()->GetBinLowEdge(etabin);
	  float higheta = factors->GetXaxis()->GetBinLowEdge(etabin+1);
	  
       if (factors->GetXaxis()->GetBinLowEdge(etabin) < 2.9) {
	  if (ptbin == 4 and higheta > 2.5) break;    // pT 120- range stops at |eta| = 2.5
	  // if (ptbin == 5 and higheta > 2.0) break;    // pT 170-1000 range stops at |eta| = 1.93

	  float lowpt = 0.001;
	  if (pts[ptbin-1] > 25) lowpt = pts[ptbin-1]; 
	    
	  txtfile << loweta << " " << higheta << " " << lowpt << " " << pts[ptbin] << " " << nparams<< " " << pts[ptbin-1] << " " << pts[ptbin] << " " << factors->GetBinContent(factors->FindBin(loweta+0.001)) << " 0 1"  << std::endl; 
       
	}
     }
 }

  txtfile.close();

}
