// Print uncertainty component for JEC uncertainties (up/down variations, stats from alpha-exptrapolations, parametrization)
// Takes as inputs files with fit functions, and also the file containing k factors (eta bins read from the histogram used to store the factors)

/// 
float pts[] = {9.0, 11.0, 13.5, 16.5, 19.5, 22.5, 26.0, 30.0, 34.5, 40.0, 46.0, 52.5, 60.0, 69.0, 79.0, 90.5, 105.5, 123.5, 143.0, 163.5, 185.0, 208.0, 232.5, 258.5, 286.0, 331.0, 396.0, 468.5, 549.5, 639.0, 738.0, 847.5, 968.5, 1102.0, 1249.5, 1412.0, 1590.5, 1787.0, 2003.0, 2241.0, 2503.0, 2790.5, 3107.0, 3455.0, 3837.0, 4257.0, 4719.0, 5226.5, 5784.0, 6538.0};
int npts =  sizeof(pts)/sizeof(int);

float ptmin = 30; float ptmax = 1000;
float maxeta = 2.5;

bool debug = false;

TString inFileName = "../withJERSF/ptfits/fits_pt_dep.root";
TString kfactFile = "../withJERSF/L2fits/correctionfile.root";
TString inFileNameUp = "../ptfits-jerup3per/fits_pt_dep_jerup3per.root";
TString kfactFileUp = "../L2fits-jerup3per/correctionfile-jersfup3per.root";
TString inFileNameDown = "../ptfits-jerdown3per/fits_pt_dep_jerdown3per.root";
TString kfactFileDown = "../L2fits-jerdown3per/correctionfile-jersfdown3per.root";

bool loglin = false; // Used log-linear parametrization for the corrections?
string header = "{1 JetEta 1 JetPt \"\" Correction JECSource}";

////

void symm_uncertainty(TH1D*, TH1D*, TH1D*, TFile*, TFile*, TFile*, bool, TString); // up-down vatiations, symmetrize as the average shift
void stat_uncertainty(TH1D*, TString); // Stat uncertainty from the k factors
void param_uncertainty(TFile*, TH1D*, bool, TString); // Uncertainty from parametrization, loglin/run3 vs. linear

void L2res_all_uncertainties( ) {

  // Nominal corrections
  TFile *inFile = new TFile(inFileName, "READ");
  TFile *inKfacts = new TFile(kfactFile, "READ");

  // These are configured to be the up-down variations from JER SF
  TFile *inFileUp = new TFile(inFileNameUp, "READ");
  TFile *inKfactsUp = new TFile(kfactFileUp, "READ");
  TFile *inFileDown = new TFile(inFileNameDown, "READ");
  TFile *inKfactsDown = new TFile(kfactFileDown, "READ");

   auto kfacts = (TH1D*)inKfacts->Get("ratio");
   auto kfactsUp = (TH1D*)inKfactsUp->Get("ratio");
   auto kfactsDown = (TH1D*)inKfactsDown->Get("ratio");


   stat_uncertainty(kfacts, "Spring23Prompt23PbPb_RelativeStat_unc.txt");
   symm_uncertainty(kfacts, kfactsUp, kfactsDown, inFile, inFileUp, inFileDown, false, "Spring23Prompt23PbPb_RelativeJER_unc.txt");
   param_uncertainty(inFile, kfacts, false, "Spring23Prompt23PbPb_RelativePt_unc.txt");
}

void stat_uncertainty(TH1D* kfacts, TString outFileName) {
   
   ofstream txtfile;
   txtfile.open(outFileName);
   
   std::cout << "Producing txt file " << outFileName << " with header " << header.c_str() << std::endl;
   txtfile << header.c_str() << std::endl;
 
   for (int etabin = kfacts->GetXaxis()->GetNbins(); etabin >= 1; --etabin) {

     float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin);

     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {
       if (higheta!=0) txtfile << -loweta << " " << -higheta << " " << npts*3;
       else txtfile << -loweta << " " << higheta << " " << npts*3;
	 
       for (int i = 0; i < npts; i++) {
	 float unc = kfacts->GetBinError(etabin);
	 // cout << " " << pts[i] << " " << unc << " " << unc;
	 txtfile << " " << pts[i] << " " << unc << " " << unc;
       }
       txtfile << endl;
     }     
   }
  
  for (int etabin = 1; etabin <= kfacts->GetXaxis()->GetNbins(); ++etabin) {
    float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin);
    float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {
       //  cout << loweta << " " << higheta << " " << npts*3;
       txtfile << loweta << " " << higheta << " " << npts*3;
	 
       for (int i = 0; i < npts; i++) {
	 float unc = kfacts->GetBinError(etabin);
	 // cout << " " << pts[i] << " " << unc << " " << unc;
	 txtfile << " " << pts[i] << " " << unc << " " << unc;
       }
       txtfile << endl;
     }
  }
  txtfile.close();
}

void param_uncertainty(TFile* inFile, TH1D* kfacts, bool loglin, TString outFileName) {

   ofstream txtfile;
   txtfile.open(outFileName);
   
   std::cout << "Producing txt file " << outFileName << " with header " << header.c_str() << std::endl;
   txtfile << header.c_str() << std::endl;

   float unc = 0;

   for (int etabin = kfacts->GetXaxis()->GetNbins(); etabin >= 1; --etabin) {

     float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin);

     auto params_eval =  (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     auto params_const =  (TF1*)inFile->Get(Form("const_eta%d",  etabin));
     
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {

       if (higheta!=0) txtfile << -loweta << " " << -higheta << " " << npts*3;
       else txtfile << -loweta << " " << higheta << " " << npts*3;
	 
       for (int i = 0; i < npts; i++) {

	 unc = kfacts->GetBinContent(etabin)*(params_eval->Eval(pts[i])-params_const->Eval(pts[i]))*0.5;
	 // cout << " " << pts[i] << " " << unc << " " << unc;
	 txtfile << " " << pts[i] << " " << unc << " " << unc;
	    
       }
       txtfile << endl;
     }     
   }
  
  for (int etabin = 1; etabin <= kfacts->GetXaxis()->GetNbins(); ++etabin) {
    float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin);
    float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);

    auto params_eval =  (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
    auto params_const =  (TF1*)inFile->Get(Form("const_eta%d",  etabin));
     
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {

       //  cout << loweta << " " << higheta << " " << npts*3;
       txtfile << loweta << " " << higheta << " " << npts*3;
	 
       for (int i = 0; i < npts; i++) {
	 unc = kfacts->GetBinContent(etabin)*(params_eval->Eval(pts[i])-params_const->Eval(pts[i]))*0.5;
	 txtfile << " " << pts[i] << " " << unc << " " << unc;
       }
       txtfile << endl;
     }     
  }
  txtfile.close();
}


void symm_uncertainty(TH1D* kfacts, TH1D* kfactsUp, TH1D* kfactsDown, TFile* inFile, TFile* inFileUp, TFile* inFileDown, bool loglin, TString outFileName) {
   
   ofstream txtfile;
   if (!debug) {
     txtfile.open(outFileName); 
     std::cout << "Producing txt file " << outFileName << " with header " << header.c_str() << std::endl;
     txtfile << header.c_str() << std::endl;
   }
   
   for (int etabin = kfacts->GetXaxis()->GetNbins(); etabin >= 1; --etabin) {

     float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin);

     auto params_nom = (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     auto params_up = (TF1*)inFileUp->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     auto params_down = (TF1*)inFileDown->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {

       if (debug)  cout << -loweta << " " << -higheta << " " << npts*3;

       if (higheta!=0) txtfile << -loweta << " " << -higheta << " " << npts*3;
       else txtfile << -loweta << " " << higheta << " " << npts*3;
	 
       for (int i = 0; i < npts; i++) {
	 float nom, up, down;

	 nom = kfacts->GetBinContent(etabin)* (params_nom->Eval(pts[i]));
	 up = kfactsUp->GetBinContent(etabin)* (params_up->Eval(pts[i]));
	 down = kfactsDown->GetBinContent(etabin)* (params_down->Eval(pts[i]));
 
	 float unc = 0.5*(abs(nom-up) + abs(nom-down) );
	 //float unc = max(abs(nom-up),abs(nom-down) ); // Alternative symmetrization 
	   
	 if (debug) cout << " " << pts[i] << " " << unc << " " << unc;
	 else txtfile << " " << pts[i] << " " << unc << " " << unc;
	    
       }
       if (!debug) txtfile << endl;
       else cout << endl;
     }

   }
   
   for (int etabin = 1; etabin <= kfacts->GetXaxis()->GetNbins(); ++etabin) {
     
     float loweta = kfacts->GetXaxis()->GetBinLowEdge(etabin);
     float higheta = kfacts->GetXaxis()->GetBinLowEdge(etabin+1);
     
     auto params_nom = (TF1*)inFile->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     auto params_up = (TF1*)inFileUp->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     auto params_down = (TF1*)inFileDown->Get(Form("%s_eta%d", (loglin ? "loglin" : "run3"), etabin));
     
     if (kfacts->GetXaxis()->GetBinLowEdge(etabin) < maxeta) {
	 
       if (debug)  cout << loweta << " " << higheta << " " << npts*3;

       else txtfile << loweta << " " << higheta << " " << npts*3;
       
       for (int i = 0; i < npts; i++) {
	 float nom, up, down;
	 
	 nom = kfacts->GetBinContent(etabin)* (params_nom->Eval(pts[i]));
	 up = kfactsUp->GetBinContent(etabin)* (params_up->Eval(pts[i]));
	 down = kfactsDown->GetBinContent(etabin)* (params_down->Eval(pts[i]));
	 
	 float unc = 0.5*(abs(nom-up) + abs(nom-down) );
	 //float unc = max(abs(nom-up),abs(nom-down) ); // Alternative symmetrization 
	 
	 if (debug) cout << " " << pts[i] << " " << unc << " " << unc;
	 else txtfile << " " << pts[i] << " " << unc << " " << unc;
	 
       }
       if (!debug) txtfile << endl;
       else cout << endl;
     }     
   }
   
   txtfile.close();

}
