
void plotresponses( string input = "L2residuals_pbpbreco_rereco_zb_jetid.root", string plottag = "ZB", bool closure = false) {

  float h1 = -1.0, h2 = 0.0, eta1 = -5.2, eta2 = 5.2;
  
  float rmax = 1.12, rmin = 0.87; // For initial checks
  //float rmax = 1.3, rmin = 0.85; // For initial checks
  // float rmax = 1.07, rmin = 0.93; // For closure to get better zoom

  float etalimit[] = {3.0, 3.0, 2.964-0.08, 2.5, 1.93+0.12};// for shady box, do manually before figuring out something smart
  bool drawbox = false;
  
  gStyle->SetOptStat(0);
 
  TFile *file = new TFile(input.c_str(),"READ");

  float alpha = 0.3;
  int pts[] = {15, 30, 40, 80, 92, 120, 1000};
  //int pts[] = {15, 30, 80, 120, 1000};

  for (int i = 0; i < 6; ++i) {
    int pt1 = pts[i];
    int pt2= pts[i+1];
    string namelabel = Form("pbpb_%dto%dalpha%.1f_%s",pt1,pt2,alpha,plottag.c_str());
    
    auto mc = (TH1D*)file->Get(Form("mc_pt%dto%d_alpha%.1f",pt1,pt2,alpha));
    auto data = (TH1D*)file->Get(Form("dt_pt%dto%d_alpha%.1f",pt1,pt2,alpha));

    if (closure) {
     mc->GetXaxis()->SetRangeUser(0,2.5);
     data->GetXaxis()->SetRangeUser(0,2.5);
    }
    
    auto leg = new TLegend(0.67,0.7,0.85,0.8); //  x, y, x, y
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(data, "Data");
    leg->AddEntry(mc, "MC");

    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    mc->GetXaxis()->SetTitle("|#eta_{probe}|");
    mc->GetYaxis()->SetTitle("Response");
    auto rp = new TRatioPlot(mc, data);
    
    rp->GetLowYaxis()->SetNdivisions(505);
   
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("MC/Data");
    rp->GetLowerRefXaxis()->SetTitle("probe jet #eta");
    
    rp->GetLowerRefGraph()->SetMaximum(rmax);
    rp->GetLowerRefGraph()->SetMinimum(rmin);

    rp->GetUpperRefYaxis()->SetRangeUser(0.89, 1.59);
  
    leg->Draw("same");

    auto line = new TLine();
    line->DrawLine(0,1.02,2.5,1.02);

    auto box = new TBox((etalimit[i]/5.19), 0.09, 0.9, 0.93);
    box->SetLineColor(kRed);
    box->SetFillColorAlpha(kBlack, 0.2);
    if (drawbox) box->Draw("same");

    auto txt = new TLatex();
    txt->SetTextSize(0.03);
    txt->DrawLatex( 0.2, 0.8, Form("%d < p_{T,avg} < %d",pt1,pt2));
    txt->DrawLatex( 0.2, 0.85, Form("#alpha < %.1f",alpha));
  
    c1->Print(Form("response_eta_%.1f_%.1f_%s.pdf",eta1,eta2,namelabel.c_str()));
 
  }
}
