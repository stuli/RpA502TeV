#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");
void mergeSixInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, TH1D*h6=0, int state=1, TString title="");
void mergeFiveInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, int state=1, TString title="");
void mergeFourInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, int state=1);
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void mergeTwoInQuadCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);
void subtractTwo( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void subtractTwoCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);

void mergeSystematicUnc_new_RFB_new2(int state = 1) { 
  
  TH1::SetDefaultSumw2();

  cout << "Make sure to make the soft link to acceptance and efficiency directory before running this macro!!" << endl;
  cout << "  ln -s ../acceptance" << endl;
  cout << "  ln -s ../efficiency" << endl;
  
  
  TH1D* hptPP1[10];
  TH1D* hptPP2[10];
  TH1D* hptPA1[10];
  TH1D* hptPA2[10];
  TH1D* hrapPP1[10];
  TH1D* hrapPP2[10];
  TH1D* hrapPA1[10];
  TH1D* hrapPA2[10];
  TH1D* hintPA[10];
  TH1D* hintPP[10];

  TH1D* hptPP[10];
  TH1D* hptPA[10];
  TH1D* hrapPP[10];
  TH1D* hrapPA[10];
  TH1D* hptRPA[10];
  TH1D* hrapRPA[10];
  
  TH1D* hHFRFB1[10];
  TH1D* hHFRFB2[10];
  TH1D* hHFRFB3[10];
  TH1D* hHFRFB4[10];
  TH1D* hNtracksRFB1[10];
  TH1D* hNtracksRFB2[10];
  TH1D* hNtracksRFB3[10];
  TH1D* hNtracksRFB4[10];
  TH1D* hrapPAF[10];
  TH1D* hrapPAB[10];
  TH1D* hrapPAFnt[10];
  TH1D* hrapPABnt[10];
  TH1D* hrapPA_rfb[10];

  TH1D* hHFRFB0[10];
  TH1D* hHFRFB[10];
  TH1D* hNtracksRFB0[10];
  TH1D* hNtracksRFB[10];

  TH1D* hptRPA1[10];
  TH1D* hptRPA2[10];
  TH1D* hrapRPA1[10];
  TH1D* hrapRPA2[10];
  TH1D* hintRPA[10]; 


  // 1 : efficiency
  TFile* f1 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_SysRFB_%dS.root",state) );
  TFile* f1_1 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");
  
  hHFRFB[1] = (TH1D*)f1_1->Get(Form("hhf%dS_yInt_Rfbdiff",state)); hHFRFB[1]->Reset();
  hHFRFB1[1]= (TH1D*)f1->Get("EffSysIntRFB");
  for(int i=1;i<=hHFRFB[1]->GetNbinsX();i++){
    hHFRFB[1]->SetBinContent(i,hHFRFB1[1]->GetBinContent(1));
  }

  hNtracksRFB[1] = (TH1D*)f1_1->Get(Form("hnt%dS_yInt_Rfbdiff",state)); hNtracksRFB[1]->Reset();
  for(int i=1;i<=hNtracksRFB[1]->GetNbinsX();i++){
    hNtracksRFB[1]->SetBinContent(i,hHFRFB1[1]->GetBinContent(1));
  }
  // 2 : acceptance
  TFile* f2 = new TFile(Form("../Acceptance/20180328/sys_acceptance_ups%dS_20180328.root",state));
  TH1D* hrapPA_rfbint = (TH1D*) f2->Get("hrapSysAccPA2bin");
 
  hHFRFB[2] = (TH1D*) hHFRFB[1]->Clone(Form("hhf%dS_yInt_RFB_acc",state)); hHFRFB[2]->Reset();
  hNtracksRFB[2] = (TH1D*) hNtracksRFB[1]->Clone(Form("hhf%dS_yInt_RFB_acc",state)); hNtracksRFB[2]->Reset();
  double sys1 = hrapPA_rfbint->GetBinContent(1);
  double sys2 = hrapPA_rfbint->GetBinContent(2);
  double sysf = TMath::Abs(1-(1+sys1)/(1+sys2));
  for(int i=1;i<=hHFRFB[2]->GetNbinsX();i++){
    hHFRFB[2]->SetBinContent(i,sysf);
  }
  for(int i=1;i<=hNtracksRFB[2]->GetNbinsX();i++){
    hNtracksRFB[2]->SetBinContent(i,sysf);
  }
  
  // 3 : signal PDF
  TFile* f3 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%ds.root",state));
  TH1D* hsignt = (TH1D*) f3->Get("hNtracksSysSigRFB000to193");
  hNtracksRFB[3] = (TH1D*) hsignt->Clone("hNtracksSysSigRFB000to193_sig");
  TH1D* hsighf = (TH1D*) f3->Get("hHFSysSigRFB000to193");
  hHFRFB[3] = (TH1D*) hsighf->Clone("hHFSysSigRFB000to193_sig");
    
  // 4 : signal Parameter
  TFile* f4 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%ds_ParamFixingOnly.root",state));
  TH1D* hsignt_ = (TH1D*) f4->Get("hNtracksSysSigRFB000to193");
  hNtracksRFB[4] = (TH1D*) hsignt_->Clone("hNtracksSysSigRFB000to193_sig");
  TH1D* hsighf_ = (TH1D*) f4->Get("hHFSysSigRFB000to193");
  hHFRFB[4] = (TH1D*) hsighf_->Clone("hHFSysSigRFB000to193_sig");
    
  // 4 : background PDF 
  TFile* f5 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");

  TH1D* hbkghf = (TH1D*)f5->Get(Form("hhf%dS_yInt_Rfbdiff",state));
  hHFRFB[5] = (TH1D*) hbkghf->Clone(Form("hhf%dS_yInt_Rfbdiff_bkg",state));
  TH1D* hbkgnt = (TH1D*)f5->Get(Form("hnt%dS_yInt_Rfbdiff",state));
  hNtracksRFB[5] = (TH1D*) hbkgnt->Clone(Form("hnt%dS_yInt_Rfbdiff_bkg",state));
  

  // Merge uncertainties for cross-section 
  hHFRFB[0] = (TH1D*) hHFRFB[1]->Clone("hHFmerged"); hHFRFB[0]->Reset();
  hNtracksRFB[0] = (TH1D*) hNtracksRFB[1]->Clone("hNtracksmerged"); hNtracksRFB[0]->Reset();
  hHFRFB[0]->SetTitle("R_{FB} in HF bins for |y_{CM}|<1.93");
  hNtracksRFB[0]->SetTitle("R_{FB} in Ntracks bins for |y_{CM}|<1.93");
  
  mergeFiveInQuad(hHFRFB[0],hHFRFB[1],hHFRFB[2],hHFRFB[3],hHFRFB[4],hHFRFB[5], state);
  mergeFiveInQuad(hNtracksRFB[0],hNtracksRFB[1],hNtracksRFB[2],hNtracksRFB[3],hNtracksRFB[4],hNtracksRFB[5],state);

  TFile* fout = new TFile(Form("mergedSys_ups%ds_rfb_new.root",state),"recreate" );
  hHFRFB[0]->Write();
  hNtracksRFB[0]->Write();
  fout->Close();
}


void mergeSixInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, int state, TString title) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);
  TH1D* hsigMerged = (TH1D*)h3->Clone("sigMerged");
  mergeTwoInQuad( hsigMerged, h3, h5);

  h0->SetAxisRange(0,0.7,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2);   h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(hsigMerged,4); hsigMerged->SetLineWidth(2); hsigMerged->DrawCopy("hist same");
  handsomeTH1(h4,6);         h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h6,8);         h6->SetLineWidth(2); h6->DrawCopy("hist same");
  

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  easyLeg(leg1,title.Data());
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(hsigMerged,"Signal PDF","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->AddEntry(h6,"TAA Uncertainty","l");
  leg1->Draw();
  c0->SaveAs(Form("pdfFiles/%s_ups%ds.pdf", h0->GetName(),state ) );
  // 6 : TAA uncertainty
  // 5 : CB+Gaus PDF  
  // 4 : background PDF
  // 3 : signal PDF
  // 2 : acceptance
  // 1 : efficiency

  
}

void mergeFiveInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D *h3, TH1D* h4, TH1D* h5, int state, TString title) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5);
    h0->SetBinContent( i, a0);
    cout << "h3 : " << i << " = " << h3->GetBinContent(i) << endl;
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);
  

  gStyle->SetOptStat(0);

  if(state!=3) h0->SetAxisRange(-0.01,0.18,"Y");
  else if(state==3) h0->SetAxisRange(-0.01,0.4,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2);   h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,4);         h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,6);         h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,11);         h5->SetLineWidth(2); h5->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  easyLeg(leg1,title.Data());
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal PDF","l");
  leg1->AddEntry(h4,"Signal Parameter","l");
  leg1->AddEntry(h5,"Background PDF","l");
  leg1->Draw();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.11);
  c0->SaveAs(Form("pdfFiles/%s_ups%ds.pdf", h0->GetName(),state ) );
  // 5 : TAA uncertainty
  // 4 : CB+Gaus PDF  
  // 3 : background PDF
  // 2 : acceptance
  // 1 : efficiency



}

void mergeFourInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, int state) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 );
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  gStyle->SetOptStat(0);

  if(state!=3) h0->SetAxisRange(0,0.18,"Y");
  else if(state==3) h0->SetAxisRange(0,0.4,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); 
  handsomeTH1(h1,        2); h1->SetLineWidth(2); 
  handsomeTH1(h2,        3); h2->SetLineWidth(2); 
  handsomeTH1(h3,        4); h3->SetLineWidth(2); 
  handsomeTH1(h4,        5); h4->SetLineWidth(2); 
  
  h0->SetLineColor(kBlack);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue+1);
  h3->SetLineColor(kViolet-3);
  h4->SetLineColor(kGreen+2);
  
  h0->DrawCopy("hist");
  h1->DrawCopy("hist same");
  h2->DrawCopy("hist same");
  h3->DrawCopy("hist same");
  h4->DrawCopy("hist same");

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal PDF","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->Draw();
  c0->SaveAs(Form("pdfFiles/%s_ups%ds.pdf", h0->GetName(),state ) );
  // 4 : background PDF
  // 3 : signal PDF
  // 2 : acceptance
  // 1 : efficiency




}

void mergeTwoInQuad( TH1D* h0, TH1D* h1, TH1D* h2) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2);
    h0->SetBinContent( i, a0);
  } 
}

void subtractTwo( TH1D* h0, TH1D* h1, TH1D* h2) {
  if ( ( h0->GetNbinsX() != h1->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
      float a1 = h1->GetBinContent(i);
      float a2 = h2->GetBinContent(i);
      float a0 = TMath::Abs((1. + a1) / ( 1. + a2) - 1); 
      h0->SetBinContent( i, a0);
    } 
  }
}
void mergeTwoInQuadCent( TH1D* h0, TH1D* hAA, TH1D* hPP) {
  if ( (hPP->GetNbinsX() != 1 ) )  {
    cout << "Number of hPP bins are not 1!! ERROR" << endl;
  }
  else if ( ( h0->GetNbinsX() != hAA->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else  {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
      float a1 = hAA->GetBinContent(i);
      float a2 = hPP->GetBinContent(1);
      float a0 = sqrt( a1*a1 + a2*a2); 
      h0->SetBinContent( i, a0);
    }
  }
}

void subtractTwoCent( TH1D* h0, TH1D* hAA, TH1D* hPP) {
  if ( (hPP->GetNbinsX() != 1 ) )  {
    cout << "Number of hPP bins are not 1!! ERROR" << endl;
  }
  else if ( ( h0->GetNbinsX() != hAA->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else  {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
      float a1 = hAA->GetBinContent(i);
      float a2 = hPP->GetBinContent(1);
      float a0 = (1. + a1) / ( 1. + a2) - 1;
      h0->SetBinContent( i, a0);
    }
  }
}
