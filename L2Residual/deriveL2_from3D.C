// Derive responses from dijet asymmetries using the 3D profile

#include "../fillhistograms/histograms.h"

// Inputs root files are
// - file with MC
// - file with DT (run for all separate data files you have)

void deriveL2_from3D(TString inFileName = "HIJEC_results/rerunall/MC_AK4_jetid.root",TString inFileNameDT = "HIJEC_results/rerunall/zerobiasall_jetid.root", TString outfilename = "RERUNALL/L2residuals_pbpbreco_rerunall_zerobias_jetid.root", int alphabin = 5, bool useabs = true, bool usewideabs = false) {


  // Open file
  TFile *inFile = new TFile(inFileName, "READ"); 
  TFile *inFileDT = new TFile(inFileNameDT, "READ");

  //// These are bins to be processed
  vector<string> etabins = {"eta_-5.2_5.2"};
  int i = 0;

  map<string, TProfile*> asymm3d, data3d, mc3d;
  map<string, TH1D*> nom, denom, responses;
  map<int, TH1D*> respETA;
  // TODO: should one save responses is like 3d thing? means: how to save responses from the different alphas, actually.

  
  TFile *outfile = new TFile(outfilename,"RECREATE");
  
  TH1D* vsalpha = new TH1D("testvsalpha","  ; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);

  // Make these as histogram template, then a map, fill histograms in a map?
  
  TH1D* vsalpha_nom = new TH1D("vsalpha_nom","  ; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);
  TH1D* vsalpha_denom = new TH1D("vsalpha_denom","  ; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);

  TH1D* vsalpha_nom_data = new TH1D("vsalpha_nom_data","  ; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);
  TH1D* vsalpha_denom_data = new TH1D("vsalpha_denom_data","  ; ;",  histograms::nalphavalues, &histograms::alphavalues[0]);


  // Here need to switch between eta bins
  
  TH1D* vseta_nom(0);
  TH1D* vseta_denom(0);

  TH1D* vseta_nom_data(0); 
  TH1D* vseta_denom_data(0); 

  TH1D *aerrormc(0), *aerrordt(0);


  
  if (useabs) {  
    vseta_nom = new TH1D("vseta_nom","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);
    vseta_denom = new TH1D("vseta_denom","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);

    vseta_nom_data = new TH1D("vseta_nom_data","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);
    vseta_denom_data = new TH1D("vseta_denom_data","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);

    aerrormc = new TH1D("aerrormc","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);
    aerrordt = new TH1D("aerrrodt","  ; ;",  histograms::nwabsetas, &histograms::wabsetarange[0]);

    
    asymm3d[etabins[i].c_str()] = (TProfile*)inFile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3Dabseta"); // Bins in order: pT, eta, alpha
    data3d[etabins[i].c_str()] = (TProfile*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3Dabseta"); // Bins in order: pT, eta, alpha
  }
  else if (usewideabs) {  
    vseta_nom = new TH1D("vseta_nom","  ; ;",  histograms::ndwabsetas, &histograms::dwabsetarange[0]);
    vseta_denom = new TH1D("vseta_denom","  ; ;",  histograms::ndwabsetas, &histograms::dwabsetarange[0]);

    vseta_nom_data = new TH1D("vseta_nom_data","  ; ;",  histograms::ndwabsetas, &histograms::dwabsetarange[0]);
    vseta_denom_data = new TH1D("vseta_denom_data","  ; ;",  histograms::ndwabsetas, &histograms::dwabsetarange[0]);

    asymm3d[etabins[i].c_str()] = (TProfile*)inFile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3Dabsetawide"); // Bins in order: pT, eta, alpha
    data3d[etabins[i].c_str()] = (TProfile*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3Dabsetawide"); // Bins in order: pT, eta, alpha
  }
  else {

    vseta_nom = new TH1D("vseta_nom","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);
    vseta_denom = new TH1D("vseta_denom","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);

    vseta_nom_data = new TH1D("vseta_nom_data","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);
    vseta_denom_data = new TH1D("vseta_denom_data","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);

    aerrormc = new TH1D("aerrormc","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);
    aerrordt = new TH1D("aerrordt","  ; ;",  histograms::nwetas, &histograms::wetarange[0]);
     
    asymm3d[etabins[i].c_str()] = (TProfile*)inFile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3D"); // Bins in order: pT, eta, alpha
    data3d[etabins[i].c_str()] = (TProfile*)inFileDT->Get("hibin_-1.0_0.0/eta_-5.2_5.2/dijetasymmetry3D"); // Bins in order: pT, eta, alpha

  }

  cout << etabins[i] <<  " nptbins  " << asymm3d[etabins[i].c_str()]->GetXaxis()->GetNbins() << endl;
  cout << etabins[i] <<  " nptbins  " << data3d[etabins[i].c_str()]->GetXaxis()->GetNbins() << endl;
  //data3d[etabins[i].c_str()]->Draw();

///////////////// Responses against eta in bin of alpha cut, pt

  for (int ptbin = 1; ptbin <= asymm3d[etabins[i].c_str()]->GetXaxis()->GetNbins(); ++ptbin) {

    cout << "Getting corrections as function of eta" << endl;
    cout << "pT bin edges: " << asymm3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin) << " " << asymm3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin+1) << endl;
    string ptstr = Form("%.0fto%.0f",asymm3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin),asymm3d[etabins[i].c_str()]->GetXaxis()->GetBinLowEdge(ptbin+1));
    cout << ptstr.c_str() << endl;
    //  data3d[etabins[i].c_str()]->Draw("same");

    cout << "alpha bin: " << alphabin << endl;
    cout << " alpha bin edges: " << asymm3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphabin) << " " << asymm3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphabin+1) << endl;
    string alphastr = Form("alpha%.1f",asymm3d[etabins[i].c_str()]->GetZaxis()->GetBinLowEdge(alphabin+1));
    
  // This is for MC
  for (int etabin = 1; etabin <= asymm3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
      
    // cout <<  asymm3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin) << endl;
      
      double val = asymm3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin);
      double err = asymm3d[etabins[i].c_str()]->GetBinError(ptbin,etabin,alphabin);

      //      cout << "Check bins " << vseta_nom->GetBinLowEdge(etabin-1) << endl;
      //  cout << "Check bins " << asymm3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin-1) << endl;
      
      vseta_nom->SetBinContent(etabin, 1+val);
      vseta_nom->SetBinError(etabin, err);

      vseta_denom->SetBinContent(etabin, 1-val);
      vseta_denom->SetBinError(etabin, err);

      //aerrormc->SetBinContent(etabin,err);
      aerrormc->SetBinContent(etabin,(TMath::IsNaN(err) ? 0. : err)*2/((1-val)*(1-val)));

    }

  TH1D* resp_eta = (TH1D*)vseta_nom->Clone("resp_eta");
  resp_eta->Divide(vseta_denom);
  for (int bin = 1; bin < resp_eta->GetXaxis()->GetNbins(); ++bin) {
    float error = aerrormc->GetBinContent(bin);
    resp_eta->SetBinError(bin, error);
      
    }
  
  resp_eta->Draw();
   
  // This is for Data
    for (int etabin = 1; etabin <= data3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
      
      //  cout << data3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin) << endl;
    
      double val = data3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin);
      double err = data3d[etabins[i].c_str()]->GetBinError(ptbin,etabin,alphabin);
     
      vseta_nom_data->SetBinContent(etabin, 1+val);
      vseta_nom_data->SetBinError(etabin, (TMath::IsNaN(err) ? 0. : err));

      vseta_denom_data->SetBinContent(etabin, 1-val);
      vseta_denom_data->SetBinError(etabin, (TMath::IsNaN(err) ? 0. : err));

      //cout << "ERROR " << err << " val " << val << " " << (TMath::IsNaN(err) ? 0. : err) << endl;

      aerrordt->SetBinContent(etabin,(TMath::IsNaN(err) ? 0. : err)*2/((1-val)*(1-val)));
      //aerrordt->SetBinContent(etabin,err);
    }

    // First part of the correction comes from the ratio of these
    
    TH1D* resp_eta_data = (TH1D*)vseta_nom_data->Clone("resp_eta_data");
    resp_eta_data->Divide(vseta_denom_data);

    for (int bin = 1; bin < resp_eta_data->GetXaxis()->GetNbins(); ++bin) {
      float error = aerrordt->GetBinContent(bin);
      resp_eta_data->SetBinError(bin, error);
    }
    
    resp_eta_data->SetLineColor(kRed);
    resp_eta_data->Draw("same");
  
    responses[Form("mc_pt%s_%s",ptstr.c_str(),alphastr.c_str())] = resp_eta; 
    responses[Form("dt_pt%s_%s",ptstr.c_str(),alphastr.c_str())] = resp_eta_data; 

    respETA[ptbin] = (TH1D*)resp_eta->Clone("ratio"); 
    respETA[ptbin]->Divide(resp_eta_data);

    responses[Form("mc_pt%s_%s",ptstr.c_str(),alphastr.c_str())]->Write(Form("mc_pt%s_%s",ptstr.c_str(),alphastr.c_str()));
    responses[Form("dt_pt%s_%s",ptstr.c_str(),alphastr.c_str())]->Write(Form("dt_pt%s_%s",ptstr.c_str(),alphastr.c_str()));
    //   responses[Form("ratio_pt%s_%s",ptstr.c_str(),alphastr.c_str())]->Write(Form("ratio_pt%s_%s",ptstr.c_str(),alphastr.c_str()));   // This is the reference for the fit
    respETA[ptbin]->Write(Form("ratio_pt%s_%s",ptstr.c_str(),alphastr.c_str()));   // This is the reference for the fit
 
  }

  // INPUTS FOR THE ISR/FSR CORRECTIONS; need to get the responses in bins of pt, eta
  cout << "number of alpha bins: " << asymm3d[etabins[i].c_str()]->GetZaxis()->GetNbins() << endl;

   
  // loop over bin in pt
  for (int ptbin = 1; ptbin <= asymm3d[etabins[i].c_str()]->GetXaxis()->GetNbins(); ++ptbin) {
  // loop over bin in eta
     cout << "NEW PT BIN " << ptbin << endl;
     for (int etabin = 1; etabin <= data3d[etabins[i].c_str()]->GetYaxis()->GetNbins(); ++etabin) {
       cout << "NEW ETA BIN " << data3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin) << " " << data3d[etabins[i].c_str()]->GetYaxis()->GetBinLowEdge(etabin+1) <<  endl;
  
       
       for (int alphabin = 1; alphabin <= asymm3d[etabins[i].c_str()]->GetZaxis()->GetNbins(); ++alphabin) {
	 double val = asymm3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin);
	 double err = asymm3d[etabins[i].c_str()]->GetBinError(ptbin,etabin,alphabin);
	 
	 vsalpha_nom->SetBinContent(alphabin, 1+val); 
	 vsalpha_nom->SetBinError(alphabin, err);

	 vsalpha_denom->SetBinContent(alphabin, 1-val);
	 vsalpha_denom->SetBinError(alphabin, err); 

	 double valdt = data3d[etabins[i].c_str()]->GetBinContent(ptbin,etabin,alphabin);
	 double errdt = data3d[etabins[i].c_str()]->GetBinError(ptbin,etabin,alphabin);
	 
	 vsalpha_nom_data->SetBinContent(alphabin, 1+valdt); 
	 vsalpha_nom_data->SetBinError(alphabin, errdt);

	 vsalpha_denom_data->SetBinContent(alphabin, 1-valdt);
	 vsalpha_denom_data->SetBinError(alphabin, errdt);

	 aerrormc->SetBinContent(etabin,(TMath::IsNaN(err) ? 0. : err)*2/((1-val)*(1-val)));
	 aerrordt->SetBinContent(etabin,(TMath::IsNaN(errdt) ? 0. : errdt)*2/((1-valdt)*(1-valdt)));

	 //cout <<  alphabin << " " << val << endl;

      }
 
     vsalpha_nom->Divide(vsalpha_denom);
     for (int bin = 1; bin < vsalpha_nom->GetXaxis()->GetNbins(); ++bin) {
      float error = aerrormc->GetBinContent(bin);
      vsalpha_nom->SetBinError(bin, error);
     }
  
     vsalpha_nom->Write(Form("Respvsa_nom_mc_%d_%d",ptbin,etabin));

     vsalpha_nom_data->Divide(vsalpha_denom_data);

     for (int bin = 1; bin < vsalpha_nom_data->GetXaxis()->GetNbins(); ++bin) {
       float error = aerrordt->GetBinContent(bin);
       vsalpha_nom_data->SetBinError(bin, error);
     }
     
     vsalpha_nom_data->Write(Form("Respvsa_denom_data_%d_%d",ptbin,etabin));
     
     vsalpha_nom->Divide(vsalpha_nom_data); // This is the MC/Data responses in bin of alpha
     vsalpha_nom->Write(Form("Respvsa_%d_%d",ptbin,etabin));
     
     TH1D* vsalpha_norm = (TH1D*)vsalpha_nom->Clone(Form("vsalpha_norm_%d_%d",ptbin,etabin));
     for (int bin = 1; bin <= vsalpha_nom->GetXaxis()->GetNbins(); ++bin) {
       
         double val =  vsalpha_nom->GetBinContent(bin);
      	 double err = vsalpha_nom->GetBinError(bin);
	 
	 double norm =  respETA[ptbin]->GetBinContent(etabin);
	 double normerr =  respETA[ptbin]->GetBinError(etabin);
	 
	 vsalpha_norm->SetBinContent(bin,val/norm);
	 vsalpha_norm->SetBinError(bin,err/norm); //

	 cout << "in alphabin " << bin << " norm with  " << norm << " lowedge " <<   vsalpha_norm->GetBinLowEdge(bin) << endl;
     }
     vsalpha_norm->Write(Form("Respvsa_norm_%d_%d",ptbin,etabin)); // This is a histogram that will be eventually fitted if looking purely at the radiation corrections
       
     }

   }
  outfile->Close();
   
}
