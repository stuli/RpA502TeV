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
void mergeFiveInQuad_pp( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, int state=1, TString title="");
void mergeFourInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, int state=1);
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void mergeTwoInQuadCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);
void subtractTwo( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void subtractTwoCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);

void mergeSystematicUnc_rap_2D_3Sbin(int state = 1) { 
  
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
  TH1D* hptPAdw[10];
  TH1D* hrapPP[10];
  TH1D* hrapPA[10];
  TH1D* hptRPA[10];
  TH1D* hrapRPA[10];
  TH1D* hHFRFB[10];
  TH1D* hNtracksRFB[10];


  TH1D* hptRPA1[10];
  TH1D* hptRPA2[10];
  TH1D* hrapRPA1[10];
  TH1D* hrapRPA2[10];
  TH1D* hintRPA[10]; 


  // 1 : efficiency
  TFile* f1 = new TFile(Form("../Efficiency/RootFiles/EffNomCor_Sys3Sbins2DRpA_%dS.root",state) );
  hrapRPA1[1] = (TH1D*)f1->Get("EffSysRapLowpT"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[1] = (TH1D*)f1->Get("EffSysRapHighpT");// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();

  // 2 : acceptance
  TFile* f2 = new TFile(Form("../Acceptance/20190221/sys_acceptance_ups%dS_20190221.root",state));
  hrapPP1[2] = (TH1D*)f2->Get("hrapSysAccPP2BinPt1");
  hrapPP2[2] = (TH1D*)f2->Get("hrapSysAccPP2BinPt2");
  hrapPA1[2] = (TH1D*)f2->Get("hrapSysAccPA2BinPt1");
  hrapPA2[2] = (TH1D*)f2->Get("hrapSysAccPA2BinPt1");
  
  hrapRPA1[2] = (TH1D*)f2->Get("hrapSysAccPA2BinPt1"); hrapRPA1[2]->Reset();// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[2] = (TH1D*)f2->Get("hrapSysAccPA2BinPt2"); hrapRPA2[2]->Reset();// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();

  subtractTwo(hrapRPA1[2], hrapPP1[2], hrapPA1[2]);
  subtractTwo(hrapRPA2[2], hrapPP2[2], hrapPA2[2]);


  // 3 : signal PDF
  TFile* f3 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%dsCombined.root",state));
  hrapPP1[3] = (TH1D*)f3->Get("hySysSigPPLowPt_in3Sbins");  
  hrapPP2[3] = (TH1D*)f3->Get("hySysSigPPHighPt_in3Sbins"); 
  hrapPA1[3] = (TH1D*)f3->Get("hySysSigPALowPt_in3Sbins");  
  hrapPA2[3] = (TH1D*)f3->Get("hySysSigPPHighPt_in3Sbins"); 

  hrapRPA1[3] = (TH1D*)f3->Get("hySysSigRpALowPt_in3Sbins"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[3] = (TH1D*)f3->Get("hySysSigRpAHighPt_in3Sbins"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();

  // 4 : signal parameter
  TFile* f4 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%ds_ParamFixingOnly.root",state));
  hrapPP1[4] = (TH1D*)f4->Get("hySysSigPPLowPt_in3Sbins");  
  hrapPP2[4] = (TH1D*)f4->Get("hySysSigPPHighPt_in3Sbins"); 
  hrapPA1[4] = (TH1D*)f4->Get("hySysSigPALowPt_in3Sbins");  
  hrapPA2[4] = (TH1D*)f4->Get("hySysSigPAHighPt_in3Sbins"); 

  hrapRPA1[4] = (TH1D*)f4->Get("hySysSigRpALowPt_in3Sbins"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[4] = (TH1D*)f4->Get("hySysSigRpAHighPt_in3Sbins"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  

  // 5 : background PDF 
  TFile* f5 = new TFile("Graham_BkgVariation/BkgPdfSystematics2dCondensed.root");
  hrapPP1[5] = (TH1D*)f5->Get(Form("hy%dS_pt06_PPDiff",state));
  hrapPP2[5] = (TH1D*)f5->Get(Form("hy%dS_pt630_PPDiff",state));
  hrapPA1[5] = (TH1D*)f5->Get(Form("hy%dS_pt06_PADiff",state));
  hrapPA2[5] = (TH1D*)f5->Get(Form("hy%dS_pt630_PADiff",state));

  hrapRPA1[5] = (TH1D*)f5->Get(Form("hy%dS_pt06_RpADiff",state));// (TH1D*)hrapPA[4]->Clone("hrapRPA_4");   hrapRPA[4]->Reset();
  hrapRPA2[5] = (TH1D*)f5->Get(Form("hy%dS_pt630_RpADiff",state));// (TH1D*)hrapPA[4]->Clone("hrapRPA_4");   hrapRPA[4]->Reset();
  

  // Merge uncertainties for cross-section 
  /*hrapPP1[0] = (TH1D*)hrapPP1[1]->Clone("hrapPP_merged1"); hrapPP1[0]->Reset();
  hrapPP2[0] = (TH1D*)hrapPP2[1]->Clone("hrapPP_merged2"); hrapPP2[0]->Reset();
  hrapPA1[0] = (TH1D*)hrapPA1[1]->Clone("hrapPA_merged1"); hrapPA1[0]->Reset();
  hrapPA2[0] = (TH1D*)hrapPA2[1]->Clone("hrapPA_merged2"); hrapPA2[0]->Reset();
*/
  // Merge uncertainties for RPA
  hrapRPA1[0] = (TH1D*)hrapRPA1[1]->Clone("hrapRPA_merged1"); hrapRPA1[0]->Reset();  
  hrapRPA2[0] = (TH1D*)hrapRPA2[1]->Clone("hrapRPA_merged2"); hrapRPA2[0]->Reset();  


  hrapRPA1[0]->SetTitle("R_{pPb} in y_{CM} bins for p_{T} < 6 GeV");
  hrapRPA2[0]->SetTitle("R_{pPb} in y_{CM} bins for 6 < p_{T} < 30 GeV");
  

  mergeFiveInQuad( hrapRPA1[0], hrapRPA1[1], hrapRPA1[2], hrapRPA1[3], hrapRPA1[4],hrapRPA1[5], state);
  mergeFiveInQuad( hrapRPA2[0], hrapRPA2[1], hrapRPA2[2], hrapRPA2[3], hrapRPA2[4],hrapRPA2[5], state);
  
//  mergeFourInQuad( hrapRPA1[0], hrapRPA1[1], hrapRPA1[2], hrapRPA1[3], hrapRPA1[5], state);
  cout << "asdasdas" << endl;
//  mergeFourInQuad( hrapRPA2[0], hrapRPA2[1], hrapRPA2[2], hrapRPA2[3], hrapRPA2[5], state);

  TFile* fout = new TFile(Form("rap2D_3Sbin_mergedSys_ups%ds.root",state),"recreate" );

  hrapRPA1[0]->Write();
  hrapRPA2[0]->Write();
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
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);
  

  gStyle->SetOptStat(0);
  
  if(state!=3) h0->SetAxisRange(-0.01,0.195,"Y");
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
  c0->SaveAs(Form("pdfFiles/rap2D_3Sbin_%s_ups%ds.pdf", h0->GetName(),state ) );
  // 5 : Bkg PDF
  // 4 : Signal Parameter 
  // 3 : signal PDF
  // 2 : acceptance
  // 1 : efficiency
}

void mergeFiveInQuad_pp( TH1D* h0, TH1D* h1, TH1D* h2, TH1D *h3, TH1D* h4, TH1D* h5, int state, TString title) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);
  

  gStyle->SetOptStat(0);
  
  if(state!=3) h0->SetAxisRange(-0.01,0.195,"Y");
  else if(state==3) h0->SetAxisRange(-0.01,0.4,"Y");
  h0->SetYTitle("Relative Uncertainty");
  h0->GetXaxis()->SetLimits(0,1.93); 
  h0->GetXaxis()->SetRangeUser(0,1.93); 

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
  // 5 : Bkg PDF
  // 4 : Signal Parameter 
  // 3 : signal PDF
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
  
  if(state!=3) h0->SetAxisRange(0,0.195,"Y");
  else if(state==3) h0->SetAxisRange(0,0.4,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(3); 
  handsomeTH1(h1,        2); h1->SetLineWidth(2); 
  handsomeTH1(h2,        3); h2->SetLineWidth(2); 
  handsomeTH1(h3,        4); h3->SetLineWidth(2); 
  handsomeTH1(h4,        5); h4->SetLineWidth(2); 
  
  h0->SetLineColor(kBlack);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue+1);
  h3->SetLineColor(kViolet-3);
  h4->SetLineColor(kGreen+2);
  
  h0->SetTitleSize(0.03);
  gStyle->SetTitleSize(0.1);
  gPad->SetLeftMargin(0.15);
  h0->GetYaxis()->SetTitleOffset(1.7);

  h0->DrawCopy("hist");
  h1->DrawCopy("hist same");
  h2->DrawCopy("hist same");
  h3->DrawCopy("hist same");
  h4->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.6,0.7, 0.9,0.9,NULL,"brNDC");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal PDF","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->Draw();
  c0->SaveAs(Form("pdfFiles/rap2D_3Sbin_%s_ups%ds.pdf", h0->GetName(),state ) );
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
