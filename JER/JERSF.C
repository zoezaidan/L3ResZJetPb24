// This version basically takes the plots without alpha cut


void JERSF(string inFileName = "HIJEC_results/debugged_plus_JER/pbpbreco_MC_JER.root", string inFileNameDT = "HIJEC_results/debugged_plus_JER/pbpbreco_DATA_closureandJER.root") {

  int pts[] = {40, 55, 80, 120, 170, 1000}; // Temporary

  float etabins[] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0}; // Temporary
  
  TFile *inFile = new TFile(inFileName.c_str(), "READ");
  TFile *inFileDT = new TFile(inFileNameDT.c_str(), "READ"); 

  // Get 3D monster
  TH3D* asymmMC = (TH3D*)inFile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/asymmdist3D");
  TH3D* asymmDT = (TH3D*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/asymmdist3D");

  //responseprofile->Draw();
  //cout << responseprofile->GetXaxis()->GetNbins() << endl;
  //cout << responseprofile->GetYaxis()->GetNbins() << endl;
  TCanvas *c1 = new TCanvas("c1","c1",800,600);

  // TODO: histos or vectors to save the sigmas


  // Need to fit these histograms and get the width; is it iterative?
  // trunct. RMS: 98.5% of events in the core
  for (int ptbin = 3; ptbin <= asymmDT->GetXaxis()->GetNbins(); ++ptbin) {
      for (int etabin = 1; etabin <= asymmDT->GetYaxis()->GetNbins(); ++etabin) {
	//	cout << asymmDT->GetBinContent(ptbin,etabin,20) << " " <<  asymmDT->GetBinError(ptbin,etabin,20) << endl;
  
	TH1D* asMC = (TH1D*)asymmMC->ProjectionZ("MC",ptbin,ptbin,etabin,etabin);
	TH1D* asDT = (TH1D*)asymmDT->ProjectionZ("DT",ptbin,ptbin,etabin,etabin);
	
	//	asDT->Scale(1./asDT->Integral(),"width"); // For some reason using "scale" leads to histograms not plotted with error bars after the first iteration
	asDT->SetTitle("");
	asDT->GetXaxis()->SetTitle("A");
	

	//	cout << asymmDT->GetBinContent(ptbin,etabin,20) << " " <<  asymmDT->GetBinError(ptbin,etabin,20) << endl;
	//      cout << asDT->GetBinContent(20) << " " <<  asDT->GetBinError(20) << endl;
	asDT->Draw("E0");
	asDT->Fit("gaus");

	TF1 *fit = asDT->GetFunction("gaus");
	auto txt = new TLatex();
	txt->SetNDC();
	txt->SetTextSize(0.03);
	txt->DrawLatex( 0.2, 0.35, Form("%.1f < |#eta| < %.1f",etabins[etabin-1],etabins[etabin]));
	txt->DrawLatex( 0.2, 0.4, Form("%d < p_{T,avg} < %d",pts[ptbin-1],pts[ptbin]));
	txt->DrawLatex( 0.2, 0.45, Form("#sigma = %.5f",fit->GetParameter(2)));
       
	//asymmDT->Draw();
	// Iterate?
	//	asDT->Fit("gaus"); // fit gaussian twice? Plot separately MC and DT? for that use the abs version?
	
	c1->SetLogy();

	// Put the labels in

	//	TF1 *fit = asDT->GetFunction("gaus");
     	//  cout << fit->GetParameter(1) << endl;    //Mean, sigma is 2
       

	c1->Print(Form("jersfhists_L2/asymm_DT_ptbin_%d_etabin_%d.pdf",ptbin,etabin));
      }
      }


  for (int ptbin = 3; ptbin <= asymmDT->GetXaxis()->GetNbins(); ++ptbin) {
      for (int etabin = 1; etabin <= asymmDT->GetYaxis()->GetNbins(); ++etabin) {
	//	cout << asymmDT->GetBinContent(ptbin,etabin,20) << " " <<  asymmDT->GetBinError(ptbin,etabin,20) << endl;
  
	TH1D* asMC = (TH1D*)asymmMC->ProjectionZ("MC",ptbin,ptbin,etabin,etabin);
	
	// asMC->Scale(1./asMC->Integral(),"width");
	// Technically we'd need the lumi normalization?
	
	asMC->SetTitle("");

	asMC->GetXaxis()->SetTitle("A");

	//	cout << asymmDT->GetBinContent(ptbin,etabin,20) << " " <<  asymmDT->GetBinError(ptbin,etabin,20) << endl;
	//      cout << asDT->GetBinContent(20) << " " <<  asDT->GetBinError(20) << endl;
	asMC->Draw("");
	asMC->Fit("gaus");
	
	c1->SetLogy();

	// Put the labels in

	TF1 *fit = asMC->GetFunction("gaus");
     	//  cout << fit->GetParameter(1) << endl;    //Mean, sigma is 2

	auto txt = new TLatex();
	txt->SetNDC();
	txt->SetTextSize(0.03);
	txt->DrawLatex( 0.2, 0.35, Form("%.1f < |#eta| < %.1f",etabins[etabin-1],etabins[etabin]));
	txt->DrawLatex( 0.2, 0.4, Form("%d < p_{T,avg} < %d",pts[ptbin-1],pts[ptbin]));
	txt->DrawLatex( 0.2, 0.45, Form("#sigma = %.5f",fit->GetParameter(2)));
       

	c1->Print(Form("jersfhists_L2/asymm_MC_ptbin_%d_etabin_%d.pdf",ptbin,etabin));
      }
  }
  
  

  // SAVE
}
