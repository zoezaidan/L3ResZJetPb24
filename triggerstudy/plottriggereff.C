// Plot trigger efficiencies as function of leading jet pT, or pT average
// Method is the reference triggers, using the ZB data (2023 pp reference with PbPb settings)

void plottriggereff(string ZBfilename = "/home/laura/Code/jec/HIJEC_rereco_results/RERECO_ZB_ALL_PFTRIG_jetid.root") {

  auto ZBfile = new TFile(ZBfilename.c_str(),"READ");

  // Jet trigger thresholds
  int pthrs[] = {40, 60};
  int nthrs = 2;
  int cols[] = {1, 2};
  
  map<int,TH1D*> leading, ptavg;

  // Reference histos passing ZB
  leading[0] = (TH1D*)ZBfile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/HLTZB");
  ptavg[0] = (TH1D*)ZBfile->Get("hibin_-1.0_0.0/eta_-5.2_5.2/HLTZB_ptav");

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  gStyle->SetOptStat(0);

  auto l0 = new TLine(0,1,1000.,1);
  auto l = new TLine(0,0.95,1000.,0.95);
  auto l2 = new TLine(80,0.,80.,1.1);
  l->SetLineStyle(kDashed);
  l2->SetLineStyle(kDashed);
  for (int i = 1; i < nthrs; ++i) {

    cout << pthrs[i] << endl;

    leading[pthrs[i]] = (TH1D*)ZBfile->Get(Form("hibin_-1.0_0.0/eta_-5.2_5.2/HLT%d",pthrs[i]));
    ptavg[pthrs[i]] = (TH1D*)ZBfile->Get(Form("hibin_-1.0_0.0/eta_-5.2_5.2/HLT%d_ptav",pthrs[i]));

    leading[pthrs[i]]->Divide(leading[pthrs[i]],leading[0],1,1,"b");
    ptavg[pthrs[i]]->Divide(ptavg[pthrs[i]],ptavg[0],1,1,"b");

    leading[pthrs[i]]->SetLineColor(cols[i]);
    ptavg[pthrs[i]]->SetLineColor(cols[i]);

    leading[pthrs[i]]->SetMaximum(1.1);
    leading[pthrs[i]]->SetMinimum(0.);
    leading[pthrs[i]]->SetTitle("");


    ptavg[pthrs[i]]->SetMaximum(1.1);
    ptavg[pthrs[i]]->SetMinimum(0.);
    ptavg[pthrs[i]]->GetXaxis()->SetRangeUser(15, 1000);
    ptavg[pthrs[i]]->SetTitle("");

    c1->cd();
    leading[pthrs[i]]->Draw("same");
    l0->Draw("same");
    leading[pthrs[i]]->Draw("same");
    l->Draw("same");
    l2->Draw("same");
    
    c2->cd();
    ptavg[pthrs[i]]->Draw("same");
    l0->Draw("same");
    ptavg[pthrs[i]]->Draw("same");
    l->Draw("same");
    l2->Draw("same");
  }

  c1->SetLogx();
  c1->Print(Form("triggerturnon_leading_pf_jetid.png"));
  
  c2->SetLogx();
  c2->Print(Form("triggerturnon_ptavg_pf_jetid.png"));

}
