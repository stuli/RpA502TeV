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
void mergeFourInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0);
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void mergeTwoInQuadCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);
void subtractTwo( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void subtractTwoCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);

void mergeSystematicUnc_asym(int state = 2, TString Asym="Hi") { 
  
  TH1::SetDefaultSumw2();

  cout << "Make sure to make the soft link to acceptance and efficiency directory before running this macro!!" << endl;
  cout << "  ln -s ../acceptance" << endl;
  cout << "  ln -s ../efficiency" << endl;
  
  
  TH1D* hptPP[10];
  TH1D* hptAA[10];
  TH1D* hrapPP[10];
  TH1D* hrapAA[10];
  TH1D* hcentAA[10];
  TH1D* hintAA[10];
  TH1D* hintPP[10];

  TH1D* hptRAA[10];
  TH1D* hrapRAA[10];
  TH1D* hcentRAA[10];
  TH1D* hintRAA[10];
 

  // for bkg deviation choice 
  TH1D* hbkgptPP[5];
  TH1D* hbkgrapPP[5];
  TH1D* hbkgintPP[5];

  TH1D* hbkgptAA[5];
  TH1D* hbkgrapAA[5];
  TH1D* hbkgintAA[5];
  TH1D* hbkgcentAA[5];
  
  TH1D* hbkgptRAA[5];
  TH1D* hbkgrapRAA[5];
  TH1D* hbkgintRAA[5];
  TH1D* hbkgcentRAA[5];



  // 1 : efficiency
  TFile* f1 = new TFile(Form("../efficiency/sys_efficiency_ups%d.root",state) );
  hptPP[1] = (TH1D*)f1->Get("hptEffPPSys");
  hptAA[1] = (TH1D*)f1->Get("hptEffAASys");
  hrapPP[1] = (TH1D*)f1->Get("hrapEffPPSys");
  hrapAA[1] = (TH1D*)f1->Get("hrapEffAASys");
  hintPP[1] = (TH1D*)f1->Get("hcentintEffPPSys");
  hintAA[1] = (TH1D*)f1->Get("hcentintEffAASys");
  hcentAA[1]= (TH1D*)f1->Get("hcentEffAASys");

  hptRAA[1] = (TH1D*)hptAA[1]->Clone("hptRAA_1");   hptRAA[1]->Reset();
  hrapRAA[1] = (TH1D*)hrapAA[1]->Clone("hrapRAA_1");   hrapRAA[1]->Reset();
  hcentRAA[1] = (TH1D*)hcentAA[1]->Clone("hcentRAA_1");    hcentRAA[1]->Reset();
  hintRAA[1] = (TH1D*)hintAA[1]->Clone("hintRAA_1");    hintRAA[1]->Reset();
  
  mergeTwoInQuad( hptRAA[1], hptAA[1], hptPP[1] );
  mergeTwoInQuad( hrapRAA[1], hrapAA[1], hrapPP[1] );
  mergeTwoInQuadCent( hcentRAA[1], hcentAA[1], hintPP[1] );
  mergeTwoInQuad( hintRAA[1], hintAA[1], hintPP[1] );
  


  // 2 : acceptance
  //TFile* f2 = new TFile(Form("../acceptance/sys_acceptance_ups%d.root",state));
  TFile* f2 = new TFile(Form("../compareDataMc/sys_acceptance_ups%dS_draft.root",state));
  hptPP[2] = (TH1D*)f2->Get("hptSysPP");
  hptAA[2] = (TH1D*)f2->Get("hptSysAA");
  hrapPP[2] = (TH1D*)f2->Get("hrapSysPP");
  hrapAA[2] = (TH1D*)f2->Get("hrapSysAA");
  hcentAA[2] = (TH1D*)hcentAA[1]->Clone("hcentSysAA");
  hcentAA[2]->Reset();
  //hcentAA[2]= (TH1D*)f2->Get("hcentSysAA");
  hintAA[2] = (TH1D*)f2->Get("hcentSysAA_int");
  hintPP[2] = (TH1D*)f2->Get("hcentSysPP");
  
  for ( int ibin = 1 ; ibin <= hcentAA[2]->GetNbinsX() ; ibin++ ) {
    hcentAA[2]->SetBinContent(ibin, 0.000000001);
  }

  
  hptRAA[2] = (TH1D*)hptAA[2]->Clone("hptRAA_2");   hptRAA[2]->Reset();
  hrapRAA[2] = (TH1D*)hrapAA[2]->Clone("hrapRAA_2");   hrapRAA[2]->Reset();
  hcentRAA[2] = (TH1D*)hcentAA[2]->Clone("hcentRAA_2");    hcentRAA[2]->Reset();
  hintRAA[2] = (TH1D*)hintAA[2]->Clone("hintRAA_2");    hintRAA[2]->Reset();

  mergeTwoInQuad( hptRAA[2], hptAA[2], hptPP[2] );
  mergeTwoInQuad( hrapRAA[2], hrapAA[2], hrapPP[2] );
  mergeTwoInQuadCent( hcentRAA[2], hcentAA[2], hintPP[2] );
  mergeTwoInQuad( hintRAA[2], hintAA[2], hintPP[2] );

  // 3 : signal PDF
  TFile* f3 = new TFile(Form("SignalVariation/sys_signalPDFVariaion_%ds.root",state));
  //TFile* f3 = new TFile(Form("SignalVariation_keep/sys_signalPDFVariaion_%ds.root",state));
  hptPP[3] = (TH1D*)f3->Get("hptPP"); hptPP[3]->SetName("hptPPsig");
  hptAA[3] = (TH1D*)f3->Get("hptAA"); hptAA[3]->SetName("hptAAsig");
  hrapPP[3] = (TH1D*)f3->Get("hrapPP"); hrapPP[3]->SetName("hrapPPsig");
  hrapAA[3] = (TH1D*)f3->Get("hrapAA"); hrapAA[3]->SetName("hrapAAsig");
  hcentAA[3]= (TH1D*)f3->Get("hcentAA"); hcentAA[3]->SetName("hcentAAsig");
  hintAA[3] = (TH1D*)f3->Get("hIntAA"); hintAA[3]->SetName("hintAAsig");
  hintPP[3] = (TH1D*)f3->Get("hIntPP"); hintPP[3]->SetName("hintPPsig");

  hptRAA[3] = (TH1D*)hptAA[3]->Clone("hptRAA_3");   hptRAA[3]->Reset();
  hrapRAA[3] = (TH1D*)hrapAA[3]->Clone("hrapRAA_3");   hrapRAA[3]->Reset();
  hcentRAA[3] = (TH1D*)hcentAA[3]->Clone("hcentRAA_3");    hcentRAA[3]->Reset();
  hintRAA[3] = (TH1D*)hintAA[3]->Clone("hintRAA_3");    hintRAA[3]->Reset();

  mergeTwoInQuad( hptRAA[3], hptAA[3], hptPP[3] );
  mergeTwoInQuad( hrapRAA[3], hrapAA[3], hrapPP[3] );
  mergeTwoInQuadCent( hcentRAA[3], hcentAA[3], hintPP[3] );
  mergeTwoInQuad( hintRAA[3], hintAA[3], hintPP[3] );

  // 4 : background PDF 
  TFile* f4 = new TFile(Form("BkgVariation/4thorder/sys_bkgPDFVariaion_4th_%ds.root",state));
  TFile* f4_1 = new TFile(Form("BkgVariation/bkglinear/sys_bkgPDFVariaion_linear_%ds.root",state));
  hbkgptPP[0] = (TH1D*)f4->Get("hptPP"); hbkgptPP[0]->SetName("hbkgptPPbkg");
  hbkgptAA[0] = (TH1D*)f4->Get("hptAA"); hbkgptAA[0]->SetName("hbkgptAAbkg");
  hbkgrapPP[0] = (TH1D*)f4->Get("hrapPP"); hbkgrapPP[0]->SetName("hbkgrapPPbkg");
  hbkgrapAA[0] = (TH1D*)f4->Get("hrapAA"); hbkgrapAA[0]->SetName("hbkgrapAAbkg");
  hbkgcentAA[0]= (TH1D*)f4->Get("hcentAA"); hbkgcentAA[0]->SetName("hbkgcentAAbkg");
  hbkgintAA[0] = (TH1D*)f4->Get("hIntAA"); hbkgintAA[0]->SetName("hbkgintAAbkg");
  hbkgintPP[0] = (TH1D*)f4->Get("hIntPP"); hbkgintPP[0]->SetName("hbkgintPPbkg");
  
  hbkgptPP[1] = (TH1D*)f4_1->Get("hptPP"); hbkgptPP[1]->SetName("hbkgptPPbkg1");
  hbkgptAA[1] = (TH1D*)f4_1->Get("hptAA"); hbkgptAA[1]->SetName("hbkgptAAbkg1");
  hbkgrapPP[1] = (TH1D*)f4_1->Get("hrapPP"); hbkgrapPP[1]->SetName("hbkgrapPPbkg1");
  hbkgrapAA[1] = (TH1D*)f4_1->Get("hrapAA"); hbkgrapAA[1]->SetName("hbkgrapAAbkg1");
  hbkgcentAA[1]= (TH1D*)f4_1->Get("hcentAA"); hbkgcentAA[1]->SetName("hbkgcentAAbkg1");
  hbkgintAA[1] = (TH1D*)f4_1->Get("hIntAA"); hbkgintAA[1]->SetName("hbkgintAAbkg1");
  hbkgintPP[1] = (TH1D*)f4_1->Get("hIntPP"); hbkgintPP[1]->SetName("hbkgintPPbkg1");
  
  hptPP[4] = (TH1D*)f4->Get("hptPP"); hptPP[4]->SetName("hptPPsig");
  hptAA[4] = (TH1D*)f4->Get("hptAA"); hptAA[4]->SetName("hptAAsig");
  hrapPP[4] = (TH1D*)f4->Get("hrapPP"); hrapPP[4]->SetName("hrapPPsig");
  hrapAA[4] = (TH1D*)f4->Get("hrapAA"); hrapAA[4]->SetName("hrapAAsig");
  hcentAA[4]= (TH1D*)f4->Get("hcentAA"); hcentAA[4]->SetName("hcentAAsig");
  hintAA[4] = (TH1D*)f4->Get("hIntAA"); hintAA[4]->SetName("hintAAsig");
  hintPP[4] = (TH1D*)f4->Get("hIntPP"); hintPP[4]->SetName("hintPPsig");

  hptRAA[4] = (TH1D*)hptAA[4]->Clone("hptRAA_4");   hptRAA[4]->Reset(); hptPP[4]->Reset(); hptAA[4]->Reset();
  hrapRAA[4] = (TH1D*)hrapAA[4]->Clone("hrapRAA_4");   hrapRAA[4]->Reset(); hrapPP[4]->Reset(); hrapAA[4]->Reset();
  hcentRAA[4] = (TH1D*)hcentAA[4]->Clone("hcentRAA_4");    hcentRAA[4]->Reset(); hcentAA[4]->Reset();
  hintRAA[4] = (TH1D*)hintAA[4]->Clone("hintRAA_4");    hintRAA[4]->Reset(); hintAA[4]->Reset(); hintPP[4]->Reset();

  hbkgptRAA[0] = (TH1D*)hptAA[4]->Clone("hbkgptRAA");   hbkgptRAA[0]->Reset();
  hbkgrapRAA[0] = (TH1D*)hrapAA[4]->Clone("hbkgrapRAA");   hbkgrapRAA[0]->Reset();
  hbkgcentRAA[0] = (TH1D*)hcentAA[4]->Clone("hbkgcentRAA");    hbkgcentRAA[0]->Reset();
  hbkgintRAA[0] = (TH1D*)hintAA[4]->Clone("hbkgintRAA");    hbkgintRAA[0]->Reset();

  hbkgptRAA[1] = (TH1D*)hptAA[4]->Clone("hbkgptRAA");   hbkgptRAA[1]->Reset();
  hbkgrapRAA[1] = (TH1D*)hrapAA[4]->Clone("hbkgrapRAA");   hbkgrapRAA[1]->Reset();
  hbkgcentRAA[1] = (TH1D*)hcentAA[4]->Clone("hbkgcentRAA");    hbkgcentRAA[1]->Reset();
  hbkgintRAA[1] = (TH1D*)hintAA[4]->Clone("hbkgintRAA");    hbkgintRAA[1]->Reset();
/*
  subtractTwo( hbkgptRAA[0], hbkgptAA[0], hbkgptPP[0] );
  subtractTwo( hbkgrapRAA[0], hbkgrapAA[0], hbkgrapPP[0] );
  subtractTwoCent( hbkgcentRAA[0], hbkgcentAA[0], hbkgintPP[0] );
  subtractTwo( hbkgintRAA[0], hbkgintAA[0], hbkgintPP[0] );

  subtractTwo( hbkgptRAA[1], hbkgptAA[1], hbkgptPP[1] );
  subtractTwo( hbkgrapRAA[1], hbkgrapAA[1], hbkgrapPP[1] );
  subtractTwoCent( hbkgcentRAA[1], hbkgcentAA[1], hbkgintPP[1] );
  subtractTwo( hbkgintRAA[1], hbkgintAA[1], hbkgintPP[1] );
*/  
/*
  for(int i=1;i<=hbkgptRAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgptRAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgptRAA[1]->GetBinContent(i))) hptRAA[4] -> SetBinContent(i,hbkgptRAA[1]->GetBinContent(i));
    else hptRAA[4] -> SetBinContent(i,hbkgptRAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgrapRAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgrapRAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgrapRAA[1]->GetBinContent(i))) hrapRAA[4] -> SetBinContent(i,hbkgrapRAA[1]->GetBinContent(i));
    else hrapRAA[4] -> SetBinContent(i,hbkgrapRAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgcentRAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgcentRAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgcentRAA[1]->GetBinContent(i))) hcentRAA[4] -> SetBinContent(i,hbkgcentRAA[1]->GetBinContent(i));
    else hcentRAA[4] -> SetBinContent(i,hbkgcentRAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgintRAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgintRAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgintRAA[1]->GetBinContent(i))) hintRAA[4] -> SetBinContent(i,hbkgintRAA[1]->GetBinContent(i));
    else hintRAA[4] -> SetBinContent(i,hbkgintRAA[0]->GetBinContent(i));
  }
*/

  for(int i=1;i<=hbkgptPP[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgptPP[0]->GetBinContent(i)) <= TMath::Abs(hbkgptPP[1]->GetBinContent(i))) hptPP[4] -> SetBinContent(i,hbkgptPP[1]->GetBinContent(i));
    else hptPP[4] -> SetBinContent(i,hbkgptPP[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgptAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgptAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgptAA[1]->GetBinContent(i))) hptAA[4] -> SetBinContent(i,hbkgptAA[1]->GetBinContent(i));
    else hptAA[4] -> SetBinContent(i,hbkgptAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgrapPP[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgrapPP[0]->GetBinContent(i)) <= TMath::Abs(hbkgrapPP[1]->GetBinContent(i))) hrapPP[4] -> SetBinContent(i,hbkgrapPP[1]->GetBinContent(i));
    else hrapPP[4] -> SetBinContent(i,hbkgrapPP[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgrapAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgrapAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgrapAA[1]->GetBinContent(i))) hrapAA[4] -> SetBinContent(i,hbkgrapAA[1]->GetBinContent(i));
    else hrapAA[4] -> SetBinContent(i,hbkgrapAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgintPP[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgintPP[0]->GetBinContent(i)) <= TMath::Abs(hbkgintPP[1]->GetBinContent(i))) hintPP[4] -> SetBinContent(i,hbkgintPP[1]->GetBinContent(i));
    else hintPP[4] -> SetBinContent(i,hbkgintPP[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgintAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgintAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgintAA[1]->GetBinContent(i))) hintAA[4] -> SetBinContent(i,hbkgintAA[1]->GetBinContent(i));
    else hintAA[4] -> SetBinContent(i,hbkgintAA[0]->GetBinContent(i));
  }
  for(int i=1;i<=hbkgcentAA[0]->GetNbinsX();i++)
  {
    if(TMath::Abs(hbkgcentAA[0]->GetBinContent(i)) <= TMath::Abs(hbkgcentAA[1]->GetBinContent(i))) hcentAA[4] -> SetBinContent(i,hbkgcentAA[1]->GetBinContent(i));
    else hcentAA[4] -> SetBinContent(i,hbkgcentAA[0]->GetBinContent(i));
  }

  mergeTwoInQuad( hptRAA[4], hptAA[4], hptPP[4] );
  mergeTwoInQuad( hrapRAA[4], hrapAA[4], hrapPP[4] );
  mergeTwoInQuadCent( hcentRAA[4], hcentAA[4], hintPP[4] );
  mergeTwoInQuad( hintRAA[4], hintAA[4], hintPP[4] );
  

  // 5 : CB+Gaus PDF 
  TFile* f5 = new TFile(Form("CBGaus_Variation/sys_CBGausVariaion_%ds.root",state));
  hptPP[5] = (TH1D*)f5->Get("hptPP"); hptPP[5]->SetName("hptPPCBGaus");
  hptAA[5] = (TH1D*)f5->Get("hptAA"); hptAA[5]->SetName("hptAACBGaus");
  hrapPP[5] = (TH1D*)f5->Get("hrapPP"); hrapPP[5]->SetName("hrapPPCBGaus");
  hrapAA[5] = (TH1D*)f5->Get("hrapAA"); hrapAA[5]->SetName("hrapAACBGaus");
  hcentAA[5]= (TH1D*)f5->Get("hcentAA"); hcentAA[5]->SetName("hcentAACBGaus");
  hintAA[5] = (TH1D*)f5->Get("hIntAA"); hintAA[5]->SetName("hintAACBGaus");
  hintPP[5] = (TH1D*)f5->Get("hIntPP"); hintPP[5]->SetName("hintPPCBGaus");
  
  hptRAA[5] = (TH1D*)hptAA[5]->Clone("hptRAA_5");   hptRAA[5]->Reset();
  hrapRAA[5] = (TH1D*)hrapAA[5]->Clone("hrapRAA_5");   hrapRAA[5]->Reset();
  hcentRAA[5] = (TH1D*)hcentAA[5]->Clone("hcentRAA_5");    hcentRAA[5]->Reset();
  hintRAA[5] = (TH1D*)hintAA[5]->Clone("hintRAA_5");    hintRAA[5]->Reset();

  mergeTwoInQuad( hptRAA[5], hptAA[5], hptPP[5] );
  mergeTwoInQuad( hrapRAA[5], hrapAA[5], hrapPP[5] );
  mergeTwoInQuadCent( hcentRAA[5], hcentAA[5], hintPP[5] );
  mergeTwoInQuad( hintRAA[5], hintAA[5], hintPP[5] );


  // 6 : TAA uncertainty 
  TFile* f6 = new TFile(Form("TAA_UNC/sys_TAA_%ds_%s.root",state,Asym.Data()));
  hptPP[6] = (TH1D*)f6->Get("hptPP"); hptPP[6]->SetName("hptPPTAA");
  hptAA[6] = (TH1D*)f6->Get("hptAA"); hptAA[6]->SetName("hptAATAA");
  hrapPP[6] = (TH1D*)f6->Get("hrapPP"); hrapPP[6]->SetName("hrapPPTAA");
  hrapAA[6] = (TH1D*)f6->Get("hrapAA"); hrapAA[6]->SetName("hrapAATAA");
  hcentAA[6]= (TH1D*)f6->Get("hcentAA"); hcentAA[6]->SetName("hcentAATAA");
  hintAA[6] = (TH1D*)f6->Get("hIntAA"); hintAA[6]->SetName("hintAATAA");
  hintPP[6] = (TH1D*)f6->Get("hIntPP"); hintPP[6]->SetName("hintPPTAA");
  
  hptRAA[6] = (TH1D*)hptAA[6]->Clone("hptRAA_6");   hptRAA[6]->Reset();
  hrapRAA[6] = (TH1D*)hrapAA[6]->Clone("hrapRAA_6");   hrapRAA[6]->Reset();
  hcentRAA[6] = (TH1D*)hcentAA[6]->Clone("hcentRAA_6");    hcentRAA[6]->Reset();
  hintRAA[6] = (TH1D*)hintAA[6]->Clone("hintRAA_6");    hintRAA[6]->Reset();

  mergeTwoInQuad( hptRAA[6], hptAA[6], hptPP[6] );
  mergeTwoInQuad( hrapRAA[6], hrapAA[6], hrapPP[6] );
  mergeTwoInQuadCent( hcentRAA[6], hcentAA[6], hintPP[6] );
  mergeTwoInQuad( hintRAA[6], hintAA[6], hintPP[6] );


  // Merge uncertainties for cross-section 
  hptPP[0] = (TH1D*)hptPP[1]->Clone("hptPP_merged"); hptPP[0]->Reset();
  hptAA[0] = (TH1D*)hptAA[1]->Clone("hptAA_merged"); hptAA[0]->Reset();
  hrapPP[0] = (TH1D*)hrapPP[1]->Clone("hrapPP_merged"); hrapPP[0]->Reset();
  hrapAA[0] = (TH1D*)hrapAA[1]->Clone("hrapAA_merged"); hrapAA[0]->Reset();
  hcentAA[0] = (TH1D*)hcentAA[1]->Clone("hcentAA_merged"); hcentAA[0]->Reset();
  hintAA[0] = (TH1D*)hintAA[1]->Clone("hintAA_merged"); hintAA[0]->Reset();
  hintPP[0] = (TH1D*)hintPP[1]->Clone("hintPP_merged"); hintPP[0]->Reset();
  hcentAA[0]->SetXTitle("Centrality x 2 (%)");


  // Merge uncertainties for RAA
  hptRAA[0] = (TH1D*)hptRAA[1]->Clone("hptRAA_merged"); hptRAA[0]->Reset();
  hrapRAA[0] = (TH1D*)hrapRAA[1]->Clone("hrapRAA_merged"); hrapRAA[0]->Reset();  
  hcentRAA[0] = (TH1D*)hcentRAA[1]->Clone("hcentRAA_merged"); hcentRAA[0]->Reset();
  hcentRAA[0]->SetXTitle("Centrality x 2 (%)");
  hintRAA[0] = (TH1D*)hintRAA[1]->Clone("hintRAA_merged"); hintRAA[0]->Reset();

  /*
  mergeSixInQuad( hptPP[0], hptPP[1], hptPP[2], hptPP[3], hptPP[4], hptPP[5], hptPP[6], state, "pp Unc. vs p_{T}" );
  mergeSixInQuad( hrapPP[0], hrapPP[1], hrapPP[2], hrapPP[3], hrapPP[4], hrapPP[5], hrapPP[6], state, "pp Unc. vs y" );
  mergeSixInQuad( hptAA[0], hptAA[1], hptAA[2], hptAA[3], hptAA[4], hptAA[5], hptAA[6], state, "PbPb Unc. vs p_{T}" );
  mergeSixInQuad( hrapAA[0], hrapAA[1], hrapAA[2], hrapAA[3], hrapAA[4], hrapAA[5], hrapAA[6] , state, "PbPb Unc. vs y");
  mergeSixInQuad( hcentAA[0], hcentAA[1], hcentAA[2], hcentAA[3], hcentAA[4], hcentAA[5], hcentAA[6] , state, "PbPb Unc. vs Centrality");
  mergeSixInQuad( hintAA[0], hintAA[1], hintAA[2], hintAA[3], hintAA[4], hintAA[5], hintAA[6] , state);
  mergeSixInQuad( hintPP[0], hintPP[1], hintPP[2], hintPP[3], hintPP[4], hintPP[5], hintPP[6] , state);

  mergeSixInQuad( hptRAA[0], hptRAA[1], hptRAA[2], hptRAA[3], hptRAA[4], hptRAA[5], hptRAA[6] , state);
  mergeSixInQuad( hrapRAA[0], hrapRAA[1], hrapRAA[2], hrapRAA[3], hrapRAA[4], hrapRAA[5], hrapRAA[6] , state);
  mergeSixInQuad( hcentRAA[0], hcentRAA[1], hcentRAA[2], hcentRAA[3], hcentRAA[4], hcentRAA[5], hcentRAA[6] , state);
  mergeSixInQuad( hintRAA[0], hintRAA[1], hintRAA[2], hintRAA[3], hintRAA[4], hintRAA[5], hintRAA[6] , state);
  */
  
  mergeFiveInQuad( hptPP[0], hptPP[1], hptPP[2], hptPP[4], hptPP[5], hptPP[6] , state);
  mergeFiveInQuad( hrapPP[0], hrapPP[1], hrapPP[2], hrapPP[4], hrapPP[5], hrapPP[6] , state);
  mergeFiveInQuad( hptAA[0], hptAA[1], hptAA[2], hptAA[4], hptAA[5], hptAA[6] , state);
  mergeFiveInQuad( hrapAA[0], hrapAA[1], hrapAA[2], hrapAA[4], hrapAA[5], hrapAA[6] , state);
  mergeFiveInQuad( hcentAA[0], hcentAA[1], hcentAA[2], hcentAA[4], hcentAA[5], hcentAA[6] , state);
  mergeFiveInQuad( hintAA[0], hintAA[1], hintAA[2], hintAA[4], hintAA[5], hintAA[6] , state);
  mergeFiveInQuad( hintPP[0], hintPP[1], hintPP[2], hintPP[4], hintPP[5], hintPP[6] , state);

  mergeFiveInQuad( hptRAA[0], hptRAA[1], hptRAA[2], hptRAA[4], hptRAA[5], hptRAA[6] , state);
  mergeFiveInQuad( hrapRAA[0], hrapRAA[1], hrapRAA[2], hrapRAA[4], hrapRAA[5], hrapRAA[6] , state);
  mergeFiveInQuad( hcentRAA[0], hcentRAA[1], hcentRAA[2], hcentRAA[4], hcentRAA[5], hcentRAA[6] , state);
  mergeFiveInQuad( hintRAA[0], hintRAA[1], hintRAA[2], hintRAA[4], hintRAA[5], hintRAA[6] , state);
  


  
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


  TFile* fout = new TFile(Form("mergedSys_ups%ds_asym%s.root",state,Asym.Data()),"recreate" );
  hptPP[0]->Write();
  hptAA[0]->Write();
  hrapPP[0]->Write();
  hrapAA[0]->Write();
  hintPP[0]->Write();
  hintAA[0]->Write();
  hcentAA[0]->Write();

  hptRAA[0]->Write();
  hrapRAA[0]->Write();
  hcentRAA[0]->Write();
  hintRAA[0]->Write();
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

void mergeFourInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 );
    h0->SetBinContent( i, a0);
  } 
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
