#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
#include "tnp_weight.h"
using namespace std;


// - Trigger:
//   * idx = 0:  nominal
//   * idx = 1..100: toy variations, stat. only
//   * idx = -1: syst variation, "new_MAX", +1 sigma
//   * idx = -2: syst variation, "new_MAX", -1 sigma
//   * idx = -10: binned
// - MuID, STA:
//   * only one SF (for systematic uncertainty only)
//double tnp_w
//      idx = 200,   MuID,
//       idx = 300,   STA

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");

void addInQuad ( TH1D* h0=0, TH1D* h1=0) ;
void getDevRatio ( TH1D* h0=0, TH1D* h1=0) ;
void getMaxTH1D ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0) ;
void addInQuadFive ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0);
void addInQuadSix ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0);
void getTnpScale(int state =1, int Nsamples=100) { 
  TH1::SetDefaultSumw2();

  TFile* f1 = new TFile(Form("efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnpWeight0_tnpIdx0.root",state) );
  TH1D* eff1 = (TH1D*)f1->Get("hptEffPP");
  TH1D* eff2 = (TH1D*)f1->Get("hptEffAA");
  TH1D* eff3 = (TH1D*)f1->Get("hrapEffPP");
  TH1D* eff4 = (TH1D*)f1->Get("hrapEffAA");
  TH1D* eff5 = (TH1D*)f1->Get("hcentintEffPP");
  TH1D* eff6 = (TH1D*)f1->Get("hcentintEffAA");
  TH1D* eff7 = (TH1D*)f1->Get("hcentEffAA");
  
  TFile* f2 = new TFile(Form("efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnpWeight1_tnpIdx0.root",state) );
  TH1D* tnp1 = (TH1D*)f2->Get("hptEffPP");
  TH1D* tnp2 = (TH1D*)f2->Get("hptEffAA");
  TH1D* tnp3 = (TH1D*)f2->Get("hrapEffPP");
  TH1D* tnp4 = (TH1D*)f2->Get("hrapEffAA");
  TH1D* tnp5 = (TH1D*)f2->Get("hcentintEffPP");
  TH1D* tnp6 = (TH1D*)f2->Get("hcentintEffAA");
  TH1D* tnp7 = (TH1D*)f2->Get("hcentEffAA");

  tnp1->Divide(eff1);
  tnp2->Divide(eff2);
  tnp3->Divide(eff3);
  tnp4->Divide(eff4);
  tnp5->Divide(eff5);
  tnp6->Divide(eff6);
  tnp7->Divide(eff7);

  stripErr(tnp1);
  stripErr(tnp2);
  stripErr(tnp3);
  stripErr(tnp4);
  stripErr(tnp5);
  stripErr(tnp6);
  stripErr(tnp7);


  TH1D* sfptAA1 = new TH1D("sfptAA1","",30,0,20);
  for ( int i =1 ; i<= sfptAA1->GetNbinsX() ; i++) { 
    float pt = sfptAA1->GetBinCenter(i);
    if ( pt < 4 ) continue;
    sfptAA1->SetBinContent( i,  tnp_weight_trg_pbpb( pt, 0.2, 0 ) );
  }
  TH1D* sfptAA2 = new TH1D("sfptAA2","",30,0,20);
  for ( int i =1 ; i<= sfptAA2->GetNbinsX() ; i++) { 
    float pt = sfptAA2->GetBinCenter(i);
    if ( pt < 4 ) continue;
    sfptAA2->SetBinContent( i,  tnp_weight_trg_pbpb( pt, 2.2, 0 ) );
  }

  TH1D* sfptPP1 = new TH1D("sfptPP1","",30,0,20);
  for ( int i =1 ; i<= sfptPP1->GetNbinsX() ; i++) { 
    float pt = sfptPP1->GetBinCenter(i);
    if ( pt < 4 ) continue;
    sfptPP1->SetBinContent( i,  tnp_weight_trg_pp( pt, 0.2, 0 ) );
  }
  TH1D* sfptPP2 = new TH1D("sfptPP2","",30,0,20);
  for ( int i =1 ; i<= sfptPP2->GetNbinsX() ; i++) { 
    float pt = sfptPP2->GetBinCenter(i);
    if ( pt < 4 ) continue;
    sfptPP2->SetBinContent( i,  tnp_weight_trg_pp( pt, 2.2, 0 ) );
  }

  
  TCanvas* c_eff_singleMu =  new TCanvas("cSingleMu","",800,400);
  c_eff_singleMu->Divide(2,1);
  c_eff_singleMu->cd(1);
  handsomeTH1(sfptPP1,1) ;
  handsomeTH1(sfptAA1,1) ;
  sfptPP1->SetMarkerStyle(24);
  sfptPP1->SetAxisRange(0.7,1.3,"Y");
  sfptPP1->Draw("p");
  sfptAA1->Draw("same p");
  c_eff_singleMu->cd(2);
  handsomeTH1(sfptPP2,1) ;
  handsomeTH1(sfptAA2,1) ;
  sfptPP2->SetMarkerStyle(24);
  sfptPP2->SetAxisRange(0.7,1.3,"Y");
  sfptPP2->Draw("p");
  sfptAA2->Draw("same p");



  TCanvas* c_eff_pt =  new TCanvas("c_eff_pt","",400,400);
  TH1D* hptEffAA;
  TH1D* hptEffPP;
  c_eff_pt->cd();
  tnp2->SetAxisRange(0.7,1.3,"Y");
  tnp2->SetYTitle("efficiency");
  tnp2->Draw("p");
  tnp1->SetAxisRange(0.7,1.3,"Y");
  tnp1->SetYTitle("efficiency");
  tnp1->SetMarkerStyle(24);
  tnp1->Draw("same p");
  TLegend* leg2 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg2,"");
  leg2->AddEntry(tnp2, "PbPb (0-100%)");
  leg2->AddEntry(tnp1, "pp");
  leg2->Draw();
  drawText(Form("#Upsilon(%dS),  p_{T}^{#mu} > 4GeV/c",state),0.25,0.87,1,15);
  jumSun(0,1,30,1); 
  //  c_eff_pt->SaveAs("tnpCorrection_pt.pdf");
 
  // Efficiency Rap
  TCanvas* c_eff_rap =  new TCanvas("c_eff_rap","",400,400);
  c_eff_rap->cd();
  tnp4 ->SetAxisRange(0.7,1.3,"Y");
  tnp4->SetYTitle("efficiency");
  tnp4->Draw("p");
  tnp3->SetAxisRange(0.7,1.3,"Y");
  tnp3->SetYTitle("efficiency");
  tnp3->SetMarkerStyle(24);
  tnp3->Draw("same p");
  TLegend* leg3 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg3,"");
  leg3->AddEntry(tnp4, "PbPb (0-100%)");
  leg3->AddEntry(tnp3, "pp");
  leg3->Draw();
  jumSun(0,1,30,1);
  //  c_eff_rap->SaveAs("tnpCorrection_rap.pdf");

  // Centrality Efficiency
  TCanvas* c_eff_cent =  new TCanvas("c_eff_cent","",400,400);
  tnp6 -> SetAxisRange(0.7,1.3,"Y");
  tnp6 -> SetYTitle("efficiency");
  tnp5 ->SetMarkerStyle(24);
  tnp5 ->SetLineStyle(2);
  tnp7 ->SetAxisRange(0,1.2,"Y");
  tnp7 ->SetYTitle("efficiency");
  tnp7 ->Draw("p");
  tnp6 -> Draw("same hist");
  tnp5 ->Draw("same hist");

  TLegend* leg4 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg4,"|y| < 2.4");
  leg4->AddEntry(tnp5, "pp","l");
  leg4->AddEntry(tnp6, "PbPb (0-100%)","l");
  leg4->AddEntry(tnp7, "PbPb","pl");
  leg4->Draw();
  drawText(Form("#Upsilon(%dS),  p_{T}^{#mu} > 4GeV/c",state),0.25,0.87,1,15);
  jumSun(0,1,200,1);
  
  
}
