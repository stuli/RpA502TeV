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

void mergeSystematicUnc_new(int state = 1) { 
  
  TH1::SetDefaultSumw2();

  cout << "Make sure to make the soft link to acceptance and efficiency directory before running this macro!!" << endl;
  cout << "  ln -s ../acceptance" << endl;
  cout << "  ln -s ../efficiency" << endl;
  
cout << "asdsadsadsa" << endl;
  
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
  TFile* f1 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_Sys2DRpA_%dS.root",state) );
  TFile* f1_1 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffCor_SyspPbXSMergedBkwdBin_%dS.root",state) );
  TFile* f1_2 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_SysRpA_%dS.root",state) );
  hptPP1[1] = (TH1D*)f1->Get("EffSysPtRapPos");   
  hptPP2[1] = (TH1D*)f1->Get("EffSysPtRapPos");   
  hptPA1[1] = (TH1D*)f1->Get("EffSysPtRapPos");   
  hptPA2[1] = (TH1D*)f1->Get("EffSysPtRapPos");   
  hrapPP1[1] = (TH1D*)f1->Get("EffSysRapLowpT");
  hrapPP2[1] = (TH1D*)f1->Get("EffSysRapLowpT");
  hrapPA1[1] = (TH1D*)f1->Get("EffSysRapLowpT"); 
  hrapPA2[1] = (TH1D*)f1->Get("EffSysRapLowpT"); 
  hintPP[1] = (TH1D*)f1_1->Get("EffSysInt"); 
  hintPA[1] = (TH1D*)f1_1->Get("EffSysInt"); 

  hptRPA1[1] = (TH1D*)f1->Get("EffSysPtRapNeg"); 
  hptRPA2[1] = (TH1D*)f1->Get("EffSysPtRapPos"); //  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hrapRPA1[1] = (TH1D*)f1->Get("EffSysRapLowpT"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[1] = (TH1D*)f1->Get("EffSysRapHighpT");// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hintRPA[1] = (TH1D*)f1_2->Get("EffSysInt"); //(TH1D*)hintPA[1]->Clone("hintRPA_1");    hintRPA[1]->Reset();

  hptPA[1] = (TH1D*) f1_1->Get("EffSysPt");
  hptPAdw[1] = (TH1D*) hptPA[1]->Clone("EffSysPtdw");
  hrapPA[1] = (TH1D*) f1_1->Get("EffSysRap");
  hrapRPA[1] = (TH1D*) f1_2->Get("EffSysRap");
  hptRPA[1] = (TH1D*) f1_2->Get("EffSysPt");

  /*  
  mergeTwoInQuad( hptRPA[1], hptPA[1], hptPP[1] );
  mergeTwoInQuad( hrapRPA[1], hrapPA[1], hrapPP[1] );
  mergeTwoInQuad( hintRPA[1], hintPA[1], hintPP[1] );
  */

  // 2 : acceptance
  TFile* f2 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20180219.root",state));
  TFile* f2_1 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20171121.root",state));
  hptPP1[2] = (TH1D*)f2->Get("hptSysAccPPRap1");   
  hptPP2[2] = (TH1D*)f2->Get("hptSysAccPPRap2");   
  hptPA1[2] = (TH1D*)f2->Get("hptSysAccPARap1");   
  hptPA2[2] = (TH1D*)f2->Get("hptSysAccPARap2");   
  hrapPP1[2] = (TH1D*)f2->Get("hrapSysAccPPPt1");
  hrapPP2[2] = (TH1D*)f2->Get("hrapSysAccPPPt2");
  hrapPA1[2] = (TH1D*)f2->Get("hrapSysAccPAPt1");
  hrapPA2[2] = (TH1D*)f2->Get("hrapSysAccPAPt2");
  hintPP[2] = (TH1D*)f2->Get("hIntSysAccPP");
  hintPA[2] = (TH1D*)f2->Get("hIntSysAccPA");
  
  hptRPA1[2] = (TH1D*)f2->Get("hptSysAccPPRap1"); hptRPA1[2]->Reset();//  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hptRPA2[2] = (TH1D*)f2->Get("hptSysAccPPRap2"); hptRPA2[2]->Reset();//  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hrapRPA1[2] = (TH1D*)f2->Get("hrapSysAccPPPt1"); hrapRPA1[2]->Reset();// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[2] = (TH1D*)f2->Get("hrapSysAccPPPt2"); hrapRPA2[2]->Reset();// (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hintRPA[2] = (TH1D*)f2->Get("hIntSysAccPA"); hintRPA[2]->Reset();//(TH1D*)hintPA[1]->Clone("hintRPA_1");    hintRPA[1]->Reset();

  hptPP[2] = (TH1D*) f2_1->Get("hptSysPP");
  hptPA[2] = (TH1D*) f2_1->Get("hptSysXsPA");
  hptPAdw[2] = (TH1D*) f2->Get("hptSysAccCross");
  TH1D* hrapPA_cal = (TH1D*) f2_1->Get("hrapSysXsPA");
  hrapPP[2] = (TH1D*) f2_1->Get("hrapSysPP");
  TH1D* hrapRPA_acc = (TH1D*)hrapPA_cal->Clone("hrapRPA_acc"); hrapRPA_acc->Reset();
  subtractTwo(hrapRPA_acc,hrapPP[2],hrapPA_cal);
  hrapRPA[2] = (TH1D*) hrapRPA[1]->Clone("hrapSysAccRpA");hrapRPA[2]->Reset();
  
  cout << "hrapPA_cal Nbins : " << hrapPA_cal->GetNbinsX() << endl;
  cout << " hrapPP[2] Nbins : " <<  hrapPP[2]->GetNbinsX() << endl;
  cout << " hrapRPA[2] Nbins : " <<  hrapRPA[2]->GetNbinsX() << endl;
  cout << " hrapRPA_acc Nbins : " <<  hrapRPA_acc->GetNbinsX() << endl;
  cout << "hrapRPA[1] : " << hrapRPA[1]->GetBinContent(1) << endl;

  for(int i=1; i<=hrapRPA[2]->GetNbinsX(); i++){
    hrapRPA[2]->SetBinContent(i, hrapRPA_acc->GetBinContent(i+1));
    cout << " hrapRPA[2] i : " << hrapRPA[2]->GetBinContent(i) << endl;
  }
  hrapPA[2] = (TH1D*) f2->Get("hrapSysAccCross"); 
  hptRPA[2] = (TH1D*) f2_1->Get("hptSysPP"); hptRPA[2]->Reset();

  subtractTwo(hptRPA[2], hptPP[2], hptPA[2]);
  subtractTwo(hptRPA1[2], hptPP1[2], hptPA1[2]);
  subtractTwo(hptRPA2[2], hptPP2[2], hptPA2[2]);
  subtractTwo(hrapRPA1[2], hrapPP1[2], hrapPA1[2]);
  subtractTwo(hrapRPA2[2], hrapPP2[2], hrapPA2[2]);
  subtractTwo(hintRPA[2], hintPP[2], hintPA[2]);


  // 3 : signal PDF
  TFile* f3 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%ds.root",state));
  hptPP1[3] = (TH1D*)f3->Get("hptSysSigPPBackwardY");
  hptPP2[3] = (TH1D*)f3->Get("hptSysSigPPForwardY"); 
  hptPA1[3] = (TH1D*)f3->Get("hptSysSigPABackwardY");
  hptPA2[3] = (TH1D*)f3->Get("hptSysSigPAForwardY"); 
  hrapPP1[3] = (TH1D*)f3->Get("hySysSigPPLowPt");  
  hrapPP2[3] = (TH1D*)f3->Get("hySysSigPPHighPt"); 
  hrapPA1[3] = (TH1D*)f3->Get("hySysSigPALowPt");  
  hrapPA2[3] = (TH1D*)f3->Get("hySysSigPAHighPt"); 
  hintPP[3] = (TH1D*)f3->Get("hintSysSigPP"); 
  hintPA[3] = (TH1D*)f3->Get("hintSysSigPA");
  hptPA[3] = (TH1D*)f3->Get("hptSysSigPA"); 
  hptPAdw[3] = (TH1D*)f3->Get("hptSysSigPA_y287to193"); 
  hrapPA[3] = (TH1D*)f3->Get("hrapSysSigCross");
  
  hptRPA1[3] = (TH1D*)f3->Get("hptSysSigRpABackwardY"); //  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hptRPA2[3] = (TH1D*)f3->Get("hptSysSigRpAForwardY"); //  (TH1D*)hptPA[1]->Clone("hptRPA_1");   hptRPA[1]->Reset();
  hrapRPA1[3] = (TH1D*)f3->Get("hySysSigRpALowPt"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hrapRPA2[3] = (TH1D*)f3->Get("hySysSigRpAHighPt"); // (TH1D*)hrapPA[1]->Clone("hrapRPA_1");   hrapRPA[1]->Reset();
  hintRPA[3] = (TH1D*)f3->Get("hintSysSigRpA"); //(TH1D*)hintPA[1]->Clone("hintRPA_1");    hintRPA[1]->Reset();
  hptRPA[3] = (TH1D*)f3->Get("hptSysSigRpA");
  hrapRPA[3] = (TH1D*)f3->Get("hySysSigRpA");

 /* TH2F* hSignalErryptPP = (TH2F*) f3_rap -> Get("hSignalErryptPP");
  TH2F* hSignalErryptPA = (TH2F*) f3_rap -> Get("hSignalErryptPA");
  TH2F* hSignalErryptRpA = (TH2F*) f3_rap -> Get("hSignalErryptRpA");
  
  */
  
  
  // 4 : background PDF 
  TFile* f4 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");
  hptPP1[4] = (TH1D*)f4->Get(Form("hpt%dS_ym_PPDiff",state));
  hptPP2[4] = (TH1D*)f4->Get(Form("hpt%dS_yp_PPDiff",state));
  hptPA1[4] = (TH1D*)f4->Get(Form("hpt%dS_ym_PADiff",state));
  hptPA2[4] = (TH1D*)f4->Get(Form("hpt%dS_yp_PADiff",state));
  hrapPP1[4] = (TH1D*)f4->Get(Form("hy%dS_pt06_PPDiff",state));
  hrapPP2[4] = (TH1D*)f4->Get(Form("hy%dS_pt630_PPDiff",state));
  hrapPA1[4] = (TH1D*)f4->Get(Form("hy%dS_pt06_PADiff",state));
  hrapPA2[4] = (TH1D*)f4->Get(Form("hy%dS_pt630_PADiff",state));
  hintPP[4] = (TH1D*)f4->Get(Form("hInt%dS_PPDiff",state));
  hintPA[4] = (TH1D*)f4->Get(Form("hInt%dS_PADiff",state));
  hptPA[4] = (TH1D*)f4->Get(Form("hpt%dS_PADiff",state));
  hptPAdw[4] = (TH1D*)f4->Get(Form("hpt%dS_cr_PADiff",state));
  hrapPA[4] = (TH1D*)f4->Get(Form("hy%dS_cr_PADiff",state));


  hptRPA1[4] = (TH1D*)f4->Get(Form("hpt%dS_ym_RpADiff",state));//(TH1D*)hptPA[4]->Clone("hptRPA_4");   hptRPA[4]->Reset();
  hptRPA2[4] = (TH1D*)f4->Get(Form("hpt%dS_yp_RpADiff",state));//(TH1D*)hptPA[4]->Clone("hptRPA_4");   hptRPA[4]->Reset();
  hrapRPA1[4] = (TH1D*)f4->Get(Form("hy%dS_pt06_RpADiff",state));// (TH1D*)hrapPA[4]->Clone("hrapRPA_4");   hrapRPA[4]->Reset();
  hrapRPA2[4] = (TH1D*)f4->Get(Form("hy%dS_pt630_RpADiff",state));// (TH1D*)hrapPA[4]->Clone("hrapRPA_4");   hrapRPA[4]->Reset();
  hintRPA[4] = (TH1D*)f4->Get(Form("hInt%dS_RpADiff",state));//(TH1D*)hintPA[4]->Clone("hintRPA_4");    hintRPA[4]->Reset();
  hptRPA[4] = (TH1D*)f4->Get(Form("hpt%dS_RpADiff",state));
  hrapRPA[4] = (TH1D*)f4->Get(Form("hy%dS_RpADiff",state));
  
  // Merge uncertainties for cross-section 
  hptPP1[0] = (TH1D*)hptPP1[1]->Clone("hptPP_merged1"); hptPP1[0]->Reset();
  hptPP2[0] = (TH1D*)hptPP2[1]->Clone("hptPP_merged2"); hptPP2[0]->Reset();
  hptPA1[0] = (TH1D*)hptPA1[1]->Clone("hptPA_merged1"); hptPA1[0]->Reset();
  hptPA2[0] = (TH1D*)hptPA2[1]->Clone("hptPA_merged2"); hptPA2[0]->Reset();
  hrapPP1[0] = (TH1D*)hrapPP1[1]->Clone("hrapPP_merged1"); hrapPP1[0]->Reset();
  hrapPP2[0] = (TH1D*)hrapPP2[1]->Clone("hrapPP_merged2"); hrapPP2[0]->Reset();
  hrapPA1[0] = (TH1D*)hrapPA1[1]->Clone("hrapPA_merged1"); hrapPA1[0]->Reset();
  hrapPA2[0] = (TH1D*)hrapPA2[1]->Clone("hrapPA_merged2"); hrapPA2[0]->Reset();
  hintPA[0] = (TH1D*)hintPA[1]->Clone("hintPA_merged"); hintPA[0]->Reset();
  hintPP[0] = (TH1D*)hintPP[1]->Clone("hintPP_merged"); hintPP[0]->Reset();
  hptPA[0] = (TH1D*)hptPA[1]->Clone("hptPA_merged"); hptPA[0]->Reset();
  hptPAdw[0] = (TH1D*)hptPAdw[1]->Clone("hptPAdw_merged"); hptPAdw[0]->Reset();
  hrapPA[0] = (TH1D*)hrapPA[1]->Clone("hrapPA_merged"); hrapPA[0]->Reset();

  // Merge uncertainties for RPA
  hptRPA1[0] = (TH1D*)hptRPA1[1]->Clone("hptRPA_merged1"); hptRPA1[0]->Reset();
  hptRPA2[0] = (TH1D*)hptRPA2[1]->Clone("hptRPA_merged2"); hptRPA2[0]->Reset();
  hrapRPA1[0] = (TH1D*)hrapRPA1[1]->Clone("hrapRPA_merged1"); hrapRPA1[0]->Reset();  
  hrapRPA2[0] = (TH1D*)hrapRPA2[1]->Clone("hrapRPA_merged2"); hrapRPA2[0]->Reset();  
  hintRPA[0] = (TH1D*)hintRPA[1]->Clone("hintRPA_merged"); hintRPA[0]->Reset();
  hptRPA[0] = (TH1D*)hptRPA[1]->Clone("hptRPA_merged"); hptRPA[0]->Reset();
  hrapRPA[0] = (TH1D*)hrapRPA[1]->Clone("hrapRPA_merged"); hrapRPA[0]->Reset();  


  hptPA1[0]->SetTitle("pPb in p_{T} bins for -1.93 < y_{CM} < 0");
  hptPA2[0]->SetTitle("pPb in p_{T} bins for 0 < y_{CM} < 1.93");
  hrapPA1[0]->SetTitle("pPb in y_{CM} bins for p_{T} < 6 GeV");
  hrapPA2[0]->SetTitle("pPb in y_{CM} bins for 6 < p_{T} < 30 GeV");
  hptPA[0]->SetTitle("pPb in p_{T} bins for -2.4 < y_{CM} < 1.93");
  hptPAdw[0]->SetTitle("pPb in p_{T} bins for -2.87 < y_{CM} <1.93");
  hrapPA[0]->SetTitle("pPb in y_{CM} bins for p_{T} < 30 GeV");
  hptRPA1[0]->SetTitle("R_{pPb} in p_{T} bins for -1.93 < y_{CM} < 0");
  hptRPA2[0]->SetTitle("R_{pPb} in p_{T} bins for 0 < y_{CM} < 1.93");
  hrapRPA1[0]->SetTitle("R_{pPb} in y_{CM} bins for p_{T} < 6 GeV");
  hrapRPA2[0]->SetTitle("R_{pPb} in y_{CM} bins for 6 < p_{T} < 30 GeV");
  hptRPA[0]->SetTitle("R_{pPb} in p_{T} bins for -1.93 < y_{CM} < 1.93");
  hrapRPA[0]->SetTitle("R_{pPb} in y_{CM} bins for p_{T} < 30 GeV");


  mergeFourInQuad( hptPP1[0], hptPP1[1], hptPP1[2], hptPP1[3],hptPP1[4],state);
  mergeFourInQuad( hptPP2[0], hptPP2[1], hptPP2[2], hptPP2[3],hptPP2[4],state);
  mergeFourInQuad( hrapPP1[0], hrapPP1[1], hrapPP1[2], hrapPP1[3], hrapPP1[4],state);
  mergeFourInQuad( hrapPP2[0], hrapPP2[1], hrapPP2[2], hrapPP2[3], hrapPP2[4],state);
  mergeFourInQuad( hptPA1[0], hptPA1[1], hptPA1[2], hptPA1[3], hptPA1[4],state);
  mergeFourInQuad( hptPA2[0], hptPA2[1], hptPA2[2], hptPA2[3], hptPA2[4],state);
  mergeFourInQuad( hrapPA1[0], hrapPA1[1], hrapPA1[2], hrapPA1[3], hrapPA1[4],state);
  mergeFourInQuad( hrapPA2[0], hrapPA2[1], hrapPA2[2], hrapPA2[3], hrapPA2[4],state);
  mergeFourInQuad( hintPA[0], hintPA[1], hintPA[2], hintPA[3], hintPA[4],state);
  mergeFourInQuad( hintPP[0], hintPP[1], hintPP[2], hintPP[3], hintPP[4],state);
  mergeFourInQuad( hrapPA[0], hrapPA[1], hrapPA[2], hrapPA[3], hrapPA[4],state);
  mergeFourInQuad( hptPA[0], hptPA[1], hptPA[2], hptPA[3], hptPA[4],state);
  mergeFourInQuad( hptPAdw[0], hptPAdw[1], hptPAdw[2], hptPAdw[3], hptPAdw[4],state);

  mergeFourInQuad( hptRPA1[0], hptRPA1[1], hptRPA1[2], hptRPA1[3], hptRPA1[4],state);
  mergeFourInQuad( hptRPA2[0], hptRPA2[1], hptRPA2[2], hptRPA2[3], hptRPA2[4],state);
  mergeFourInQuad( hrapRPA1[0], hrapRPA1[1], hrapRPA1[2], hrapRPA1[3], hrapRPA1[4],state);
  mergeFourInQuad( hrapRPA2[0], hrapRPA2[1], hrapRPA2[2], hrapRPA2[3], hrapRPA2[4],state);
  mergeFourInQuad( hintRPA[0], hintRPA[1], hintRPA[2], hintRPA[3], hintRPA[4],state);
  mergeFourInQuad( hptRPA[0], hptRPA[1], hptRPA[2], hptRPA[3], hptRPA[4],state);
  for(int i=0;i<=4;i++){
  hrapRPA[i]->GetXaxis()->SetRangeUser(-1.93,1.93);
  hrapRPA[i]->GetXaxis()->SetLimits(-1.93,1.93);
  }
  
  mergeFourInQuad( hrapRPA[0], hrapRPA[1], hrapRPA[2], hrapRPA[3], hrapRPA[4],state);

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
  hptPP1[0]->Write();
  hptPP2[0]->Write();
  hptPA1[0]->Write();
  hptPA2[0]->Write();
  hrapPP1[0]->Write();
  hrapPP2[0]->Write();
  hrapPA1[0]->Write();
  hrapPA2[0]->Write();
  hptPA[0]->Write();
  hptPAdw[0]->Write();
  hrapPA[0]->Write();
  //hrapCrossPA[0]->Write();

  hptRPA[0]->Write();
  hrapRPA[0]->Write();
  hptRPA1[0]->Write();
  hptRPA2[0]->Write();
  hrapRPA1[0]->Write();
  hrapRPA2[0]->Write();
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
