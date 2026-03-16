// Do fits for pt-parametization of L2residuals

#include "../fillhistograms/histograms.h"

// 2023 bins
const float pts[] = { 30, 40, 80, 92, 120, 1000};

int nptbins = 5;
string ptbins[] = {"30to40", "40to80", "80to92", "92to120", "120to1000"};

void fit_pt_param(TString inzb = "L2residuals_pbpbreco_rereco_zb_jetid.root", TString inHP =  "L2residuals_pbpbreco_rereco_hp_jetid.root", float fitmin = 25., float fitmax = 1000, string outfilename = "testing_pt_dep", string outfolder = "ptfits", bool doabseta = true) {

  // input file is from 3D derivation
  gStyle->SetOptStat(0);

  TFile *inFilezb = new TFile(inzb, "READ");
  TFile *inFile = new TFile(inHP, "READ");

  gSystem->mkdir(outfolder.c_str());
  TFile *outfile = new TFile(Form("%s/%s.root",outfolder.c_str(),outfilename.c_str()),"RECREATE");
  
 
  // This is for correcton factors - hiso against eta
  auto factors = (TH1D*)inFile->Get(Form("ratio_pt%s_alpha0.3",ptbins[0].c_str()));
  factors->Reset();

  int colours[] = {209, 226, 213, 51, 206, 209};

  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  c3->SetLogx();

  // Fro mresponses
  map<string, TH1D*> histos;
  map<int, TH1D*> histosvspt;

  // Picking up histograms for different pT bins from ZB/HP
  for (int i = 0; i < 2; ++i)   histos[ptbins[i]] = (TH1D*)inFilezb->Get(Form("ratio_pt%s_alpha0.3",ptbins[i].c_str()));
  for (int i = 2; i < nptbins; ++i)   histos[ptbins[i]] = (TH1D*)inFile->Get(Form("ratio_pt%s_alpha0.3",ptbins[i].c_str()));

  auto txt = new TLatex();
  txt->SetTextSize(0.03);
  txt->SetNDC();

  TF1 * f1 = new TF1("f1","[0] + [1]*log(x)",fitmin,fitmax);
  TF1 * f2 = new TF1("f2","pol0",fitmin,fitmax);
  TF1 * f3 = new TF1("f3","1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))",fitmin,fitmax); // Run3 parametrization
  
  // Histograms vs. pT
  for (int ebin = 1; ebin < 12; ++ebin) {   // 12 corresponds to 2.5
    auto leg2 = new TLegend(0.12,0.10,0.35,0.36); //  x, y, x, y
    leg2->SetTextSize(0.03);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    
    histosvspt[ebin] = new TH1D(Form("ebin_%d",ebin),"",nptbins,&pts[0]); 

    // etabins for
    for (int ptbin = 1; ptbin <= nptbins; ++ptbin) {
      histosvspt[ebin]->SetBinContent(ptbin,histos[ptbins[ptbin-1].c_str()]->GetBinContent(ebin));
      histosvspt[ebin]->SetBinError(ptbin,histos[ptbins[ptbin-1].c_str()]->GetBinError(ebin));
    }

    histosvspt[ebin]->SetMaximum(1.1);
    histosvspt[ebin]->SetMinimum(0.9);
    histosvspt[ebin]->GetXaxis()->SetTitle("p_{T,avg} (GeV)");
    histosvspt[ebin]->Draw();
 
    f1->SetParameters(1,0);
    histosvspt[ebin]->Fit("f1","R");
  
    f1->Draw("same");

    f2->SetParameters(1);
    f2->SetLineStyle(kDashed);
    histosvspt[ebin]->Fit("f2","R");
    f2->Draw("same");

    f3->SetLineStyle(kDashed);
    f3->SetLineColor(kBlue+1);
    histosvspt[ebin]->Fit("f3","R");
    f3->Draw("same");

    leg2->AddEntry(f1, Form("a + b*log(p_{T}), a = %.3f #pm %.3f, b = %.3f #pm %.3f",f1->GetParameter(0),f1->GetParError(0),f1->GetParameter(1),f1->GetParError(1)));
    leg2->AddEntry(f1, Form("#chi2/ndof = %.3f / %d",  f1->GetChisquare(), f1->GetNDF()));
    leg2->AddEntry(f2, Form("p0 = %.3f #pm %.3f",f2->GetParameter(0), f2->GetParError(0)));
    leg2->AddEntry(f2, Form("#chi2/ndof = %.3f / %d",  f2->GetChisquare(), f2->GetNDF()));
    leg2->AddEntry(f3, "1./(a+b*log10(0.01*p_{T})+c/(p_{T}/10.))");
        leg2->AddEntry(f3, Form(" a = %.3f #pm %.3f, b = %.3f #pm %.3f, c = %.3f #pm %.3f",f3->GetParameter(0), f3->GetParError(0),f3->GetParameter(1), f3->GetParError(1),f3->GetParameter(2), f3->GetParError(2)));
    leg2->AddEntry(f3, Form("#chi2/ndof = %.3f / %d",  f3->GetChisquare(), f3->GetNDF()));
    
    leg2->Draw();

    txt->DrawLatex(0.5,0.8,Form("%.3f < |#eta| < %.3f ",histos[ptbins[0]]->GetBinLowEdge(ebin),histos[ptbins[0]]->GetBinLowEdge(ebin+1)));

    f1->Write(Form("loglin_eta%d",ebin));
    f2->Write(Form("const_eta%d",ebin));
    f3->Write(Form("run3_eta%d",ebin));
    c3->Print(Form("%s/ebin_%d.png",outfolder.c_str(),ebin));
  }
  
  outfile->Close();
  
}

