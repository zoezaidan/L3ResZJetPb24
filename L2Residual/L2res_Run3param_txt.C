// Print L2residual correction txt file when using pT-dependent parametrization in each bin of pt
// Takes as input a file with fit functions, and also the file containing k factors (eta bins read from the histogram used to store the factors)
// flag "loglin" can be used to switch between loglin and constant fit

int nparams = 6;

void L2res_Run3param_txt(TString inFileName = "ptfits/testing_pt_dep.root", TString kfactFile = "L2fits_AK4CHSjetveto_all/correctionfile.root", TString outFileName = "L2residual-param.txt", bool loglin = true, bool kfsr = false) {

   TFile *inFile = new TFile(inFileName, "READ");
   TFile *inKfacts = new TFile(kfactFile, "READ");

   auto kfacts = (TH1D*)inKfacts->Get("ratio");
   
   string header = "{1 JetEta 1 JetPt 1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))*[3] Correction L2Relative}";
   
   ofstream txtfile;
   txtfile.open(outFileName);
   
   std::cout << "Producing txt file " << outFileName << " with header " << header.c_str() << std::endl;
   txtfile << header.c_str() << std::endl;

   map<int, TH1D*> facts;
 
   float ptmin = 30; float ptmax = 1000;
   for (int etabin = kfacts->GetXaxis()->GetNbins(); etabin >= 1; --etabin) {

     float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin);

     auto params =  (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "const"),  etabin));
     //     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < 2.9) {
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < 2.5) {
	    
       if (higheta != 0) txtfile << -loweta << " " << -higheta << " " << nparams << " " <<  ptmin << " " << ptmax << " " << params->GetParameter(0) << " " << params->GetParameter(1) << " " <<  params->GetParameter(2) <<  " " << kfacts->GetBinContent(etabin) << std::endl;
       else txtfile << -loweta << " " << higheta << " "  << nparams<< " "  << " " <<  ptmin << " " << ptmax << " " << params->GetParameter(0) << " " << params->GetParameter(1) << " " <<  params->GetParameter(2) <<  " " << (kfsr ? kfacts->GetBinContent(etabin) : 0)<< std::endl;
	 
     }     
   }
  
  for (int etabin = 1; etabin <= kfacts->GetXaxis()->GetNbins(); ++etabin) {
    float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin);
    float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);

    auto params =  (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "const"),  etabin));

    if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < 2.5) {
      txtfile << loweta << " " << higheta << " "  << nparams<< " "  << " " <<  ptmin << " " << ptmax << " " << params->GetParameter(0) << " " <<  params->GetParameter(1) << " " <<  params->GetParameter(2) <<  " " << (kfsr ? kfacts->GetBinContent(etabin) : 0) << std::endl;
    }
  }

  txtfile.close();

}
