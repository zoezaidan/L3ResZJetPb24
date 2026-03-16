#include "../fillhistograms/histograms.h"

// Do the fits as function of alpha

// https://root-forum.cern.ch/t/getrms-on-th1d-class/41881
double GetRootMeanSquare( const TH1* h1) { 
   double stats[TH1::kNstat]; 
   h1->GetStats(stats); 
   return (stats[0] > 0) ? std::sqrt(stats[3]/stats[0]) : 0;
}

//void JERSF_RMS(string inFileName = "HIJEC_results/rerunall_combinedbins/MC_AK4_jetid.root", string inFileNameZB = "HIJEC_results/rerunall_combinedbins/zerobiasall_jetid_l2corr.root", string inFileNameDT = "HIJEC_results/rerunall_combinedbins/HP_AK4_PFTRIG_jetid_l2corr.root") {
//  void JERSF_RMS(string inFileName = "HIJEC_results/rerunall_combinedbins/MC_AK4_PFTRIG_jetid_l2corr_forjer.root", string inFileNameZB = "HIJEC_results/rerunall_combinedbins/zerobiasall_jetid_l2corr_forjer.root", string inFileNameDT = "HIJEC_results/rerunall_combinedbins/HP_AK4_PFTRIG_jetid_l2corr_forjer.root") {
// latest with old reco
//void JERSF_RMS(string inFileName = "HIJEC_results/rerunall_combinedbins/MC_AK4_PFTRIG_jetid_l2corr_forjer_wideeta.root", string inFileNameZB = "HIJEC_results/rerunall_combinedbins/zerobiasall_jetid_l2corr_forjer_wideeta.root", string inFileNameDT = "HIJEC_results/rerunall_combinedbins/HP_AK4_PFTRIG_jetid_l2corr_forjer_wideeta.root") {

// RERECO with old corrections
void JERSF_RMS(TString outFileName = "JERSF_sigmas_RMS.root", string inFileName = "/home/laura/Code/jec/HIJEC_rereco_results/RERECOMC_AK4_PFTRIG_jetid_l2jecclosure_forjer.root", string inFileNameZB = "/home/laura/Code/jec/HIJEC_rereco_results/RERECO_ZB_AK4_PFTRIG_jetid_l2closure_forjer.root", string inFileNameDT = "/home/laura/Code/jec/HIJEC_rereco_results/RERECOHP_AK4_PFTRIG_jetid_l2jecclosure_forjer.root") {


  float cut = 0.985; // cut for trunctuating the |A| histograms
  
  //  int pts[] = {40, 55, 80, 120, 170, 1000}; // Temporary
  int pts[] = {15, 25, 80, 120, 1000}; // Temporary

  // float etabins[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}; // Temporary
   float etabins[] = {0, 1.3, 2.5, 3.0}; // Temporary
  
  TFile *inFile = new TFile(inFileName.c_str(), "READ");
  TFile *inFileDT = new TFile(inFileNameDT.c_str(), "READ");
  TFile *inFileZB = new TFile(inFileNameZB.c_str(), "READ"); 

  // Get 3D monsters: 
  TH3D* asymmMC = (TH3D*)inFile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/asymmdist3D");
  TH3D* asymmDT = (TH3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/asymmdist3D");
  TH3D* asymmZB = (TH3D*)inFileZB->Get("hibin_-1.0_0.0/eta_-5.2_5.2/asymmdist3D");

  map<int, TH3D*> asymmMCs;
  map<int, TH3D*> asymmDTs, asymmZBs; // _a10 to a55?
  int alphabins[] = {10, 15, 20, 25, 30, 35, 40, 45};
 
  // Outputs
  map<int, map<int, TH1D*>> sigmasMCmap, sigmasDTmap;
  map<int, map<int, TGraphErrors*>> sigmasGraphMCmap, sigmasGraphDTmap;

  TH1D* sigmasMC = new TH1D("sigmasMC","; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);
  TH1D* sigmasDT = new TH1D("sigmasDT","; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);

  for (int ptbin = 2; ptbin <= asymmDT->GetXaxis()->GetNbins(); ++ptbin) {
    for (int etabin = 1; etabin <= asymmDT->GetYaxis()->GetNbins(); ++etabin) {
      sigmasMCmap[ptbin][etabin] = new TH1D(Form("sigmasMCpt%deta%d",ptbin,etabin),"; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);
      sigmasDTmap[ptbin][etabin] = new TH1D(Form("sigmasDTpt%deta%d",ptbin,etabin),"; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);
    }
  }

  TFile *outfile = new TFile(outFileName,"RECREATE");
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  
  int nalphas= 7;
  for (int i = 0; i < nalphas; ++i) {
    cout << "alpha: " << alphabins[i] << " " << Form("hibin_-1.0_0.0/eta_-5.2_5.2/absasymmdist3D_a%d",alphabins[i]) <<  " file " << inFileNameDT.c_str() << endl;
    asymmMCs[i] = (TH3D*)inFile->Get(Form("hibin_-1.0_0.0/eta_-5.2_5.2/absasymmdist3D_a%d",alphabins[i]));
    asymmDTs[i] = (TH3D*)inFileDT->Get(Form("hibin_-1.0_0.0/eta_-5.2_5.2/absasymmdist3D_a%d",alphabins[i]));
    asymmZBs[i] = (TH3D*)inFileZB->Get(Form("hibin_-1.0_0.0/eta_-5.2_5.2/absasymmdist3D_a%d",alphabins[i]));


    
    // if (i != 0) break; // debug

   // trunct. RMS: 98.5% of events in the core
    int alphabin = i + 1;// binning in array vs. histogram
    for (int ptbin = 2; ptbin <= asymmDTs[i]->GetXaxis()->GetNbins(); ++ptbin) {
      for (int etabin = 1; etabin <= asymmDTs[i]->GetYaxis()->GetNbins(); ++etabin) {

	//	cout << asymmDT->GetBinContent(ptbin,etabin,20) << " " <<  asymmDT->GetBinError(ptbin,etabin,20) << endl;
  
	// 	TH1D* asMC = (TH1D*)asymmMC->ProjectionZ("MC",ptbin,ptbin,etabin,etabin);
	//TH1D* asDT = (TH1D*)asymmDT->ProjectionZ("DT",ptbin,ptbin,etabin,etabin);
	TH1D* asDT(0); 
	// Take HP instead of ZB; TODO: be more flexible for 2024
	if ( ptbin > 3 ) asDT = (TH1D*)asymmDTs[i]->ProjectionZ(Form("DT_%d%d%d",ptbin,etabin,alphabin),ptbin,ptbin,etabin,etabin);
	else asDT = (TH1D*)asymmZBs[i]->ProjectionZ(Form("DT_%d%d%d",ptbin,etabin,alphabin),ptbin,ptbin,etabin,etabin);

	TH1D* asMC = (TH1D*)asymmMCs[i]->ProjectionZ(Form("MC_%d%d%d",ptbin,etabin,alphabin),ptbin,ptbin,etabin,etabin);
	asDT->GetXaxis()->SetRangeUser(0,1.4);
	asMC->GetXaxis()->SetRangeUser(0,1.4);

	
        asDT->Scale(1./asDT->Integral(),"");
	asMC->Scale(1./asMC->Integral(),"");
	
	asDT->SetTitle("");
	asDT->GetXaxis()->SetTitle("|A|");
	cout << "Integrate: " <<  asDT->Integral() << " get entries: " << asDT->GetEntries() << endl;
	cout << "Scan histrogram" << endl;


	//// DATA
	float einbins = 0, fracdt = 0, fracmc = 0, linepointdt = 0, linepointmc = 0;
	float all =  asDT->Integral();
	if (all != all) { cout << "No entries in data, skipping this bin" << endl; continue; }
	
	TH1D* trunctAsDT = (TH1D*)asDT->Clone(Form("trunct%d",etabin));
	trunctAsDT->Reset();
	//	for (int i = 1; i < asDT->GetXaxis()->GetNbins(); ++i) {

	for (int j = 1;  einbins/all < cut ; ++j) {
	  float einbinsold = einbins;
	  einbins += asDT->GetBinContent(j);
	  // Debug: 
	  cout << " in bins " << einbins << " binc " << asDT->GetBinContent(j) << " ratio " <<   einbins/all << endl;

	  if (einbins/all < cut) {
	    trunctAsDT->SetBinContent(j, asDT->GetBinContent(j));
	  }
	  else {
	    fracdt = (cut*all-einbinsold)/ asDT->GetBinContent(j);
	    //    cout << "FRAC is: " << fracdt << " in bin then " << fracdt*asDT->GetBinContent(j) << " " << asDT->GetBinContent(j) << endl; // << " " << cut << " " << all << " " << einbins << endl;
	    trunctAsDT->SetBinContent(j, fracdt*asDT->GetBinContent(j));
	    linepointdt = fracdt*asDT->GetBinWidth(j) + asDT->GetBinLowEdge(j);
	  }
	}
	double RMSDT =  GetRootMeanSquare(trunctAsDT);
	//	sigmasDT->SetBinContent(alphabin,RMSDT);
	sigmasDTmap[ptbin][etabin]->SetBinContent(alphabin,RMSDT);
	cout << "RMS in trunct data: " << RMSDT << endl;

	//// MC
	einbins = 0;
        all =  asMC->Integral();
	TH1D* trunctAsMC = (TH1D*)asMC->Clone(Form("trunct%d",etabin));
	trunctAsMC->Reset();

	for (int j = 1;  einbins/all < cut ; ++j) { // TODO: while loop
	  float einbinsold = einbins;
	  einbins += asMC->GetBinContent(j);
	  // Debug: 
	  //	   cout << " in bins " << einbins << " binc " << asMC->GetBinContent(j) << " ratio " <<   einbins/all << endl;

	  if (einbins/all < cut) trunctAsMC->SetBinContent(j, asMC->GetBinContent(j));
	  else {
	    fracmc = (cut*all-einbinsold)/ asMC->GetBinContent(j);
	    //	    cout << "FRAC is: " << fracmc << " in bin then " << fracmc*asMC->GetBinContent(j) << " " << asMC->GetBinContent(j) << endl; // << " " << cut << " " << all << " " << einbins << endl;
	    trunctAsMC->SetBinContent(j, fracmc*asMC->GetBinContent(j));
	    linepointmc = fracmc*asMC->GetBinWidth(j) + asMC->GetBinLowEdge(j);
	  }
	}
	double RMSMC = GetRootMeanSquare(trunctAsMC);
	//	sigmasMC->SetBinContent(alphabin,RMSMC);
	sigmasMCmap[ptbin][etabin]->SetBinContent(alphabin,RMSMC);
	cout << "RMS in trunct MC: " << RMSMC << endl;


	////////// Draw data
	auto linedt  = new TLine(linepointdt, 0, linepointdt, 0.3); // TODO: max
	asDT->SetStats(0);
	asDT->Draw(""); 
	linedt->Draw();
	//	trunctAsDT->Draw("same");

	///////// Draw MC
	auto linemc  = new TLine(linepointmc, 0, linepointmc, 0.3);
	linemc->SetLineStyle(kDashed);
	asMC->SetLineColor(kRed);
	asMC->Draw("same");
	linemc->Draw();
	trunctAsMC->SetLineColor(kRed);
	//	trunctAsMC->Draw("same");
	
	//	asDT->Fit("gaus");
	auto leg = new TLegend(0.6, 0.43, 0.85, 0.53);
	leg->SetTextSize(0.032);
	leg->SetBorderSize(0);
	leg->AddEntry(asDT,Form("Data"),"l");
	leg->AddEntry(asMC,Form("MC"),"l");
	leg->Draw();
	
	auto leg2 = new TLegend(0.6, 0.12, 0.85, 0.22);
	leg2->SetTextSize(0.032);
	leg2->SetBorderSize(0);
	leg2->AddEntry(linedt,Form("%.1f%% data",cut*100),"l");
	leg2->AddEntry(linemc,Form("%.1f%% mc",cut*100),"l");
	leg2->Draw();

	auto txt = new TLatex();
	txt->SetNDC();
	txt->SetTextSize(0.03);
	txt->DrawLatex( 0.62, 0.4, Form("%d < p_{T,avg} < %d",pts[ptbin-1],pts[ptbin]));
	txt->DrawLatex( 0.62, 0.35, Form("%.1f < |#eta| < %.1f, #alpha < %.2f",etabins[etabin-1],etabins[etabin],0.01*alphabins[i]));
	txt->DrawLatex( 0.62, 0.3, Form("RMS DT: %g",RMSDT));
	txt->DrawLatex( 0.62, 0.25, Form("RMS MC: %g",RMSMC));
       
	c1->SetLogy();
	
        c1->Print(Form("jersfhists/absasymm_RMS_ptbin_%d_etabin_%d_a%d.png",ptbin,etabin,i));

      } // etabin

    } // ptbin

  }


  for (int ptbin = 2; ptbin <= asymmDT->GetXaxis()->GetNbins(); ++ptbin) {
    for (int etabin = 1; etabin <= asymmDT->GetYaxis()->GetNbins(); ++etabin) {
      sigmasMCmap[ptbin][etabin]->Write(Form("sigmasMCpt%deta%d",ptbin,etabin));
      sigmasDTmap[ptbin][etabin]->Write(Form("sigmasDTpt%deta%d",ptbin,etabin));

      Double_t results[histograms::nalphavaluesgraphjer];
      Double_t errors[histograms::nalphavaluesgraphjer];
      Double_t errorsx[histograms::nalphavaluesgraphjer];

      Double_t resultsdt[histograms::nalphavaluesgraphjer];
      Double_t errorsdt[histograms::nalphavaluesgraphjer];
      
      for (int abin = 0; abin < histograms::nalphavaluesgraphjer; ++abin) {
	
        results[abin] =  sigmasMCmap[ptbin][etabin]->GetBinContent(abin+1);
	errors[abin] =   0; //sigmasMCmap[ptbin][etabin]->GetBinError(abin+1);

	resultsdt[abin] =  sigmasDTmap[ptbin][etabin]->GetBinContent(abin+1);
	errorsdt[abin] =   0;//  sigmasDTmap[ptbin][etabin]->GetBinError(abin+1);
       }

      sigmasGraphMCmap[ptbin][etabin] = new TGraphErrors(histograms::nalphavaluesgraphjer, histograms::alphavaluesgraphjer, results, errorsx, errors); // TODO: do not include reference?
      sigmasGraphMCmap[ptbin][etabin]->Write(Form("gsigmasMCpt%deta%d",ptbin,etabin));

      sigmasGraphDTmap[ptbin][etabin] = new TGraphErrors(histograms::nalphavaluesgraphjer, histograms::alphavaluesgraphjer, resultsdt, errorsx, errorsdt); // TODO: do not include reference?
      sigmasGraphDTmap[ptbin][etabin]->Write(Form("gsigmasDTpt%deta%d",ptbin,etabin));
    }
  }

  outfile->Close();
}

