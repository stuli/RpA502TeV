#include "../../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../cutsAndBin.h"
#include "../../multiTreeUtil.h"
using namespace std;
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void mergeTwoInQuadCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);
void subtractTwo( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void subtractTwoCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);

void bkgMerge(int state=1)
{
  TH1D* hptPP[10];
  TH1D* hptAA[10];
  TH1D* hrapPP[10];
  TH1D* hrapAA[10];
  TH1D* hcentAA[10];
  TH1D* hintAA[10];
  TH1D* hintPP[10];

  TFile* f1 = new TFile(Form("bkglinear/sys_bkgPDFVariaion_linear_%ds.root",state));
  hptPP[0] = (TH1D*)f1->Get("hptPP"); hptPP[0]->SetName("hptPP1");
  hptAA[0] = (TH1D*)f1->Get("hptAA"); hptAA[0]->SetName("hptAA1");
  hrapPP[0] = (TH1D*)f1->Get("hrapPP"); hrapPP[0]->SetName("hrapPP1");
  hrapAA[0] = (TH1D*)f1->Get("hrapAA"); hrapAA[0]->SetName("hrapAA1");
  hcentAA[0]= (TH1D*)f1->Get("hcentAA"); hcentAA[0]->SetName("hcentAA1");
  hintAA[0] = (TH1D*)f1->Get("hIntAA"); hintAA[0]->SetName("hintAA1");
  hintPP[0] = (TH1D*)f1->Get("hIntPP"); hintPP[0]->SetName("hintPP1");

  TFile* f2 = new TFile(Form("4thorder/sys_bkgPDFVariaion_4th_%ds.root",state));
  //TFile* f2 = new TFile(Form("../ToyMC/sys_toyMC_bkg_%ds.root",state));
  hptPP[1] = (TH1D*)f2->Get("hptPP"); hptPP[1]->SetName("hptPP2");
  hptAA[1] = (TH1D*)f2->Get("hptAA"); hptAA[1]->SetName("hptAA2");
  hrapPP[1] = (TH1D*)f2->Get("hrapPP"); hrapPP[1]->SetName("hrapPP2");
  hrapAA[1] = (TH1D*)f2->Get("hrapAA"); hrapAA[1]->SetName("hrapAA2");
  hcentAA[1]= (TH1D*)f2->Get("hcentAA"); hcentAA[1]->SetName("hcentAA2");
  hintAA[1] = (TH1D*)f2->Get("hIntAA"); hintAA[1]->SetName("hintAA2");
  hintPP[1] = (TH1D*)f2->Get("hIntPP"); hintPP[1]->SetName("hintPP2");

  hptPP[2] = (TH1D*)hptPP[1]->Clone("hptPP");   hptPP[2]->Reset();
  hptAA[2] = (TH1D*)hptAA[1]->Clone("hptAA");   hptAA[2]->Reset();
  hrapPP[2] = (TH1D*)hrapPP[1]->Clone("hrapPP");   hrapPP[2]->Reset();
  hrapAA[2] = (TH1D*)hrapAA[1]->Clone("hrapAA");   hrapAA[2]->Reset();
  hcentAA[2] = (TH1D*)hcentAA[1]->Clone("hcentAA");   hcentAA[2]->Reset();
  hintPP[2] = (TH1D*)hintPP[1]->Clone("hIntPP");   hintPP[2]->Reset();
  hintAA[2] = (TH1D*)hintAA[1]->Clone("hIntAA");   hintAA[2]->Reset();

  mergeTwoInQuad( hptPP[2], hptPP[0], hptPP[1] );
  mergeTwoInQuad( hptAA[2], hptAA[0], hptAA[1] );
  mergeTwoInQuad( hrapPP[2], hrapPP[0], hrapPP[1] );
  mergeTwoInQuad( hrapAA[2], hrapAA[0], hrapAA[1] );
  mergeTwoInQuad( hintPP[2], hintPP[0], hintPP[1] );
  mergeTwoInQuad( hintAA[2], hintAA[0], hintAA[1] );
  mergeTwoInQuad( hcentAA[2], hcentAA[0], hcentAA[1] );

  TFile* fout = new TFile(Form("mergedBkg_ups%ds.root",state),"recreate" );
  
  hptPP[2]->Write();
  hptAA[2]->Write();
  hrapPP[2]->Write();
  hrapAA[2]->Write();
  hintPP[2]->Write();
  hintAA[2]->Write();
  hcentAA[2]->Write();
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
