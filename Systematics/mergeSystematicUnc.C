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

void mergeSystematicUnc(int state = 1) { 
  
  TH1::SetDefaultSumw2();

  cout << "Make sure to make the soft link to acceptance and efficiency directory before running this macro!!" << endl;
  cout << "  ln -s ../acceptance" << endl;
  cout << "  ln -s ../efficiency" << endl;
  
  
  TH1D* hptPP[10];
  TH1D* hptPA[10];
  TH1D* hrapPP[10];
  TH1D* hrapPA[10];
  TH1D* hintPA[10];
  TH1D* hintPP[10];
  
  TH1D* hptRPA[10];
  TH1D* hrapRPA[10];
  TH1D* hintRPA[10]; 


  // 1 : efficiency
  TFile* f1 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20171114.root",state) );
  //TFile* f1 = new TFile(Form("../efficiency/sys_efficiency_ups%d.root",state) );
  hptPP[1] = (TH1D*)f1->Get("hptSysPP");   hptPP[1]->Reset();
  hptPA[1] = (TH1D*)f1->Get("hptSysPA");   hptPA[1]->Reset();
  hrapPP[1] = (TH1D*)f1->Get("hrapSysPP"); hrapPP[1]->Reset();
  hrapPA[1] = (TH1D*)f1->Get("hrapSysPA"); hrapPA[1]->Reset();
  hintPP[1] = (TH1D*)f1->Get("hIntSysPP"); hintPP[1]->Reset();
  hintPA[1] = (TH1D*)f1->Get("hIntSysPA"); hintPA[1]->Reset();
  
  hptRPA[1] = (TH1D*)f1->Get("hptSysPP"); hptRPA[1]->Reset();//  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hrapRPA[1] = (TH1D*)f1->Get("hrapSysPP"); hrapRPA[1]->Reset();// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hintRPA[1] = (TH1D*)f1->Get("hIntSysPP"); hintRPA[1]->Reset();//(TH1D*)hintPA[1]->Clone("hintRPA_1");    hintRPA[1]->Reset();
  /*  
  mergeTwoInQuad( hptRPA[1], hptPA[1], hptPP[1] );
  mergeTwoInQuad( hrapRPA[1], hrapPA[1], hrapPP[1] );
  mergeTwoInQuad( hintRPA[1], hintPA[1], hintPP[1] );
  */


  // 2 : acceptance
  TFile* f2 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20171114.root",state));
  hptPP[2] = (TH1D*)f2->Get("hptSysPP");
  hptPA[2] = (TH1D*)f2->Get("hptSysPA");
  hrapPP[2] = (TH1D*)f2->Get("hrapSysPP");
  hrapPA[2] = (TH1D*)f2->Get("hrapSysPA");
  hintPP[2] = (TH1D*)f2->Get("hIntSysPP");
  hintPA[2] = (TH1D*)f2->Get("hIntSysPA");
  
  hptRPA[2] = (TH1D*)hptPA[2]->Clone("hptRPA_2");   hptRPA[2]->Reset();
  hrapRPA[2] = (TH1D*)hrapPA[2]->Clone("hrapRPA_2");   hrapRPA[2]->Reset();
  hintRPA[2] = (TH1D*)hintPA[2]->Clone("hintRPA_2");    hintRPA[2]->Reset();
  

  // 3 : signal PDF
  TFile* f3 = new TFile(Form("Jared_SignalShapeVariation/HistoSystematicErrorSignal%ds.root",state));
  hptPP[3] = (TH1D*)f3->Get("hSignalErrptPP");
  hptPA[3] = (TH1D*)f3->Get("hSignalErrptPA");
  hrapPP[3] = (TH1D*)f3->Get("hSignalErryPP");
  hrapPA[3] = (TH1D*)f3->Get("hSignalErryPA");
//  hintPP[3] = (TH1D*)f3->Get("hintSigPPSys");
//  hintPA[3] = (TH1D*)f3->Get("hintSigPASys");
  
  hptRPA[3] = (TH1D*)f3->Get("hSignalErrptRpA");// (TH1D*)hptPA[3]->Clone("hptRPA_3");   hptRPA[3]->Reset();
  hrapRPA[3] = (TH1D*)f3->Get("hSignalErryRpA");//(TH1D*)hrapPA[3]->Clone("hrapRPA_3");   hrapRPA[3]->Reset();
//  hintRPA[3] = (TH1D*)f3->Get("hintSigRpPbSys");//(TH1D*)hintPA[3]->Clone("hintRPA_3");    hintRPA[3]->Reset();
 /* 
  mergeTwoInQuad( hptRPA[3], hptPA[3], hptPP[3] );
  mergeTwoInQuad( hrapRPA[3], hrapPA[3], hrapPP[3] );
  mergeTwoInQuad( hintRPA[3], hintPA[3], hintPP[3] );
  */

  // 4 : background PDF 
  TFile* f4 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");
  hptPP[4] = (TH1D*)f4->Get(Form("hpt%dS_PPDiff",state));
  hptPA[4] = (TH1D*)f4->Get(Form("hpt%dS_PADiff",state));
  hrapPP[4] = (TH1D*)f4->Get(Form("hy%dS_PPDiff",state));
  hrapPA[4] = (TH1D*)f4->Get(Form("hy%dS_PADiff",state));
  hintPP[4] = (TH1D*)f4->Get(Form("hInt%dS_PPDiff",state));
  hintPA[4] = (TH1D*)f4->Get(Form("hInt%dS_PADiff",state));
  
  hptRPA[4] = (TH1D*)f4->Get(Form("hpt%dS_RpADiff",state));//(TH1D*)hptPA[4]->Clone("hptRPA_4");   hptRPA[4]->Reset();
  hrapRPA[4] = (TH1D*)f4->Get(Form("hy%dS_RpADiff",state));// (TH1D*)hrapPA[4]->Clone("hrapRPA_4");   hrapRPA[4]->Reset();
  hintRPA[4] = (TH1D*)f4->Get(Form("hInt%dS_RpADiff",state));//(TH1D*)hintPA[4]->Clone("hintRPA_4");    hintRPA[4]->Reset();
 /* 
  mergeTwoInQuad( hptRPA[4], hptPA[4], hptPP[4] );
  mergeTwoInQuad( hrapRPA[4], hrapPA[4], hrapPP[4] );
  mergeTwoInQuad( hintRPA[4], hintPA[4], hintPP[4] );
  */

  // Merge uncertainties for cross-section 
  hptPP[0] = (TH1D*)hptPP[1]->Clone("hptPP_merged"); hptPP[0]->Reset();
  hptPA[0] = (TH1D*)hptPA[1]->Clone("hptPA_merged"); hptPA[0]->Reset();
  hrapPP[0] = (TH1D*)hrapPP[1]->Clone("hrapPP_merged"); hrapPP[0]->Reset();
  hrapPA[0] = (TH1D*)hrapPA[1]->Clone("hrapPA_merged"); hrapPA[0]->Reset();
  hintPA[0] = (TH1D*)hintPA[1]->Clone("hintPA_merged"); hintPA[0]->Reset();
  hintPP[0] = (TH1D*)hintPP[1]->Clone("hintPP_merged"); hintPP[0]->Reset();

  // Merge uncertainties for RPA
  hptRPA[0] = (TH1D*)hptRPA[1]->Clone("hptRPA_merged"); hptRPA[0]->Reset();
  hrapRPA[0] = (TH1D*)hrapRPA[1]->Clone("hrapRPA_merged"); hrapRPA[0]->Reset();  
  hintRPA[0] = (TH1D*)hintRPA[1]->Clone("hintRPA_merged"); hintRPA[0]->Reset();

  hintPP[3] = (TH1D*)hintPP[2]->Clone("hintPP_3");hintPP[3]->Reset();
  hintPA[3] = (TH1D*)hintPA[2]->Clone("hintPA_3");hintPA[3]->Reset();
  hintRPA[3] = (TH1D*)hintRPA[2]->Clone("hintRPA_3");hintRPA[3]->Reset();

  mergeFourInQuad( hptPP[0], hptPP[1], hptPP[2], hptPP[3],hptPP[4],state);
  mergeFourInQuad( hrapPP[0], hrapPP[1], hrapPP[2], hrapPP[3], hrapPP[4],state);
  mergeFourInQuad( hptPA[0], hptPA[1], hptPA[2], hptPA[3], hptPA[4],state);
  mergeFourInQuad( hrapPA[0], hrapPA[1], hrapPA[2], hrapPA[3], hrapPA[4],state);
  mergeFourInQuad( hintPA[0], hintPA[1], hintPA[2], hintPA[3], hintPA[4],state);
  mergeFourInQuad( hintPP[0], hintPP[1], hintPP[2], hintPP[3], hintPP[4],state);

  mergeFourInQuad( hptRPA[0], hptRPA[1], hptRPA[2], hptRPA[3], hptRPA[4],state);
  mergeFourInQuad( hrapRPA[0], hrapRPA[1], hrapRPA[2], hrapRPA[3], hrapRPA[4],state);
  mergeFourInQuad( hintRPA[0], hintRPA[1], hintRPA[2], hintRPA[3], hintRPA[4],state);
/*  
  TCanvas* c1= new TCanvas("c1","",800,800);
  c1->Divide(2,4);
  c1->cd(1);
  hptPP[0]->Draw();
  c1->cd(2);
  hptAA[0]->Draw();
  c1->cd(3);
  hrapPP[0]->Draw();
  c1->cd(4);
  hrapAA[0]->Draw();
  c1->cd(5);
  hintPP[0]->Draw();
  c1->cd(6);
  hintAA[0]->Draw();
  c1->cd(7);
  hcentAA[0]->Draw();

  TCanvas* c2= new TCanvas("c2","",800,800);
  c2->Divide(2,2);
  c2->cd(1);
  hptRAA[0]->Draw();
  c2->cd(2);
  hrapRAA[0]->Draw();
  c2->cd(3);
  hcentRAA[0]->Draw();
  c2->cd(4);
  hintRAA[0]->Draw();

*/
  TFile* fout = new TFile(Form("mergedSys_ups%ds.root",state),"recreate" );
  hptPP[0]->Write();
  hptPA[0]->Write();
  hrapPP[0]->Write();
  hrapPA[0]->Write();
  hintPP[0]->Write();
  hintPA[0]->Write();
  //hrapCrossPA[0]->Write();

  hptRPA[0]->Write();
  hrapRPA[0]->Write();
  hintRPA[0]->Write();
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
  

  h0->SetAxisRange(-0.5,1.1,"Y");
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
  leg1->AddEntry(h3,"Background PDF","l");
  leg1->AddEntry(h4,"Signal PDF","l");
  leg1->AddEntry(h5,"TAA Uncertainty","l");
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

  h0->SetAxisRange(0,0.7,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2);   h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  

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
      float a0 = (1. + a1) / ( 1. + a2) - 1; 
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
