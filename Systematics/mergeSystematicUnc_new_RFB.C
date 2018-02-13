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

void mergeSystematicUnc_new_RFB(int state = 1) { 
  
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
  TH1D* hNtracksRFB0[10];



  TH1D* hptRPA1[10];
  TH1D* hptRPA2[10];
  TH1D* hrapRPA1[10];
  TH1D* hrapRPA2[10];
  TH1D* hintRPA[10]; 


  // 1 : efficiency
  TFile* f1 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/EffNomCor_SysRFB_%dS.root",state) );
  TFile* f1_1 = new TFile(Form("../CrossChecks/efficiency_Santona/ForAnalysisNote/EffCor_SyspPbXS_%dS.root",state) );
  TFile* f1_2 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");
  

  hrapPA_rfb[1]=(TH1D*)f1->Get("EffSysRapRFB");  
  hrapPA_rfb[2]=(TH1D*)f1->Get("EffSysIntRFB");
  for(int i=0;i<10;i++){
    hHFRFB1[i]=(TH1D*)f1_2->Get(Form("hhf%dS_yInt_Rfbdiff",state)); hHFRFB1[i]->Reset();
    hHFRFB2[i]=(TH1D*)hHFRFB1[i] -> Clone(Form("hhf2_%dS_%d",state,i+1)); hHFRFB2[i]->Reset();
    hHFRFB3[i]=(TH1D*)hHFRFB1[i] -> Clone(Form("hhf3_%dS_%d",state,i+1)); hHFRFB3[i]->Reset();
    hHFRFB4[i]=(TH1D*)hHFRFB1[i] -> Clone(Form("hhf4_%dS_%d",state,i+1)); hHFRFB4[i]->Reset();
    hHFRFB0[i]=(TH1D*)hHFRFB1[i] -> Clone(Form("hhf0_%dS_%d",state,i+1)); hHFRFB0[i]->Reset();
    hNtracksRFB1[i]=(TH1D*)f1_2->Get(Form("hnt%dS_yInt_Rfbdiff",state)); hNtracksRFB1[i]->Reset();
    hNtracksRFB2[i]=(TH1D*)hNtracksRFB1[i] -> Clone(Form("hnt2_%dS_%d",state,i+1)); hNtracksRFB2[i]->Reset();
    hNtracksRFB3[i]=(TH1D*)hNtracksRFB1[i] -> Clone(Form("hnt3_%dS_%d",state,i+1)); hNtracksRFB3[i]->Reset();
    hNtracksRFB4[i]=(TH1D*)hNtracksRFB1[i] -> Clone(Form("hnt4_%dS_%d",state,i+1)); hNtracksRFB4[i]->Reset();
    hNtracksRFB0[i]=(TH1D*)hNtracksRFB1[i] -> Clone(Form("hnt0_%dS_%d",state,i+1)); hNtracksRFB0[i]->Reset();
  }
  if(state==1){
    for(int i=1; i<=hHFRFB1[1]->GetNbinsX();i++){
      hHFRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hHFRFB2[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(2));
      hHFRFB3[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(3));
      hHFRFB4[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(4));
      hHFRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
      hNtracksRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hNtracksRFB2[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(2));
      hNtracksRFB3[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(3));
      hNtracksRFB4[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(4));
      hNtracksRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
    }
  }
  else if(state==2){
    for(int i=1; i<=hHFRFB1[1]->GetNbinsX();i++){
      hHFRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hHFRFB2[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(2));
      hHFRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
      hNtracksRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hNtracksRFB2[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(2));
      hNtracksRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
    }
  }
  else if(state==3){
    for(int i=1; i<=hHFRFB1[1]->GetNbinsX();i++){
      hHFRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hHFRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
      hNtracksRFB1[1]->SetBinContent(i,hrapPA_rfb[1]->GetBinContent(1));
      hNtracksRFB0[1]->SetBinContent(i,hrapPA_rfb[2]->GetBinContent(1));
    }
  }


  // 2 : acceptance
  TFile* f2 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20180213.root",state));
  TFile* f2_1 = new TFile(Form("../Acceptance/sys_acceptance_ups%dS_20171121.root",state));
  hrapPA_rfb[3] = (TH1D*) f2->Get("hrapSysAccCross");
  hrapPA_rfb[4] = (TH1D*) f2->Get("hrapSysAccCross");
  TH1D* hrapPA_rfbint = (TH1D*) f2->Get("hrapSysAccPA2bin");
  for(int i=0;i<10;i++){
    hrapPAF[i] = (TH1D*)f1_2->Get("hhf1S_yInt_Rfbdiff"); hrapPAF[i]->Reset();
    hrapPAB[i] = (TH1D*)f1_2->Get("hhf1S_yInt_Rfbdiff"); hrapPAB[i]->Reset();
    hrapPAFnt[i] = (TH1D*)f1_2->Get("hnt1S_yInt_Rfbdiff"); hrapPAFnt[i]->Reset();
    hrapPABnt[i] = (TH1D*)f1_2->Get("hnt1S_yInt_Rfbdiff"); hrapPABnt[i]->Reset();
  }
  if(state==1){
    for(int i=1;i<=hrapPAF[1]->GetNbinsX();i++){
      hrapPAF[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPAB[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAF[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(8));
      hrapPAB[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(5));
      hrapPAF[3]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(9));
      hrapPAB[3]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(4));
      hrapPAF[4]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(10));
      hrapPAB[4]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(3));
      hrapPAF[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPAB[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    for(int i=1;i<=hrapPAFnt[1]->GetNbinsX();i++){
      hrapPAFnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPABnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAFnt[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(8));
      hrapPABnt[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(5));
      hrapPAFnt[3]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(9));
      hrapPABnt[3]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(4));
      hrapPAFnt[4]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(10));
      hrapPABnt[4]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(3));
      hrapPAFnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPABnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    subtractTwo(hHFRFB1[2],hrapPAB[1],hrapPAF[1]);
    subtractTwo(hHFRFB2[2],hrapPAB[2],hrapPAF[2]);
    subtractTwo(hHFRFB3[2],hrapPAB[3],hrapPAF[3]);
    subtractTwo(hHFRFB4[2],hrapPAB[4],hrapPAF[4]);
    subtractTwo(hHFRFB0[2],hrapPAB[5],hrapPAF[5]);
    subtractTwo(hNtracksRFB1[2],hrapPABnt[1],hrapPAFnt[1]);
    subtractTwo(hNtracksRFB2[2],hrapPABnt[2],hrapPAFnt[2]);
    subtractTwo(hNtracksRFB3[2],hrapPABnt[3],hrapPAFnt[3]);
    subtractTwo(hNtracksRFB4[2],hrapPABnt[4],hrapPAFnt[4]);
    subtractTwo(hNtracksRFB0[2],hrapPABnt[5],hrapPAFnt[5]);
  }
  else if(state==2){
    for(int i=1;i<=hrapPAF[1]->GetNbinsX();i++){
      hrapPAF[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPAB[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAF[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(8));
      hrapPAB[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(5));
      hrapPAF[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPAB[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    for(int i=1;i<=hrapPAFnt[1]->GetNbinsX();i++){
      hrapPAFnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPABnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAFnt[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(8));
      hrapPABnt[2]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(5));
      hrapPAFnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPABnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    subtractTwo(hHFRFB1[2],hrapPAB[1],hrapPAF[1]);
    subtractTwo(hHFRFB2[2],hrapPAB[2],hrapPAF[2]);
    subtractTwo(hHFRFB0[2],hrapPAB[5],hrapPAF[5]);
    subtractTwo(hNtracksRFB1[2],hrapPABnt[1],hrapPAFnt[1]);
    subtractTwo(hNtracksRFB2[2],hrapPABnt[2],hrapPAFnt[2]);
    subtractTwo(hNtracksRFB0[2],hrapPABnt[5],hrapPAFnt[5]);
  }
  else if(state==3){
    for(int i=1;i<=hrapPAF[1]->GetNbinsX();i++){
      hrapPAF[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPAB[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAF[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPAB[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    for(int i=1;i<=hrapPAFnt[1]->GetNbinsX();i++){
      hrapPAFnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(7));
      hrapPABnt[1]->SetBinContent(i,hrapPA_rfb[3]->GetBinContent(6));
      hrapPAFnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(2));
      hrapPABnt[5]->SetBinContent(i,hrapPA_rfbint->GetBinContent(1));
    }
    subtractTwo(hHFRFB1[2],hrapPAB[1],hrapPAF[1]);
    subtractTwo(hHFRFB0[2],hrapPAB[5],hrapPAF[5]);
    subtractTwo(hNtracksRFB1[2],hrapPABnt[1],hrapPAFnt[1]);
    subtractTwo(hNtracksRFB0[2],hrapPABnt[5],hrapPAFnt[5]);
  }

  // 3 : signal PDF
  TFile* f3 = new TFile(Form("Jared_SignalShapeVariation/ErrorEstimates/SysSig%ds.root",state));
  if(state==1){
    hHFRFB1[3] = (TH1D*)f3->Get("hHFSysSigRFB000to040");
    hHFRFB2[3] = (TH1D*)f3->Get("hHFSysSigRFB040to080");
    hHFRFB3[3] = (TH1D*)f3->Get("hHFSysSigRFB080to120");
    hHFRFB4[3] = (TH1D*)f3->Get("hHFSysSigRFB120to193");
    hHFRFB0[3] = (TH1D*)f3->Get("hHFSysSigRFB000to193");
    hNtracksRFB1[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to040");
    hNtracksRFB2[3] = (TH1D*)f3->Get("hNtracksSysSigRFB040to080");
    hNtracksRFB3[3] = (TH1D*)f3->Get("hNtracksSysSigRFB080to120");
    hNtracksRFB4[3] = (TH1D*)f3->Get("hNtracksSysSigRFB120to193");
    hNtracksRFB0[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to193");
  }
  else if(state==2){
    hHFRFB1[3] = (TH1D*)f3->Get("hHFSysSigRFB000to080");
    hHFRFB2[3] = (TH1D*)f3->Get("hHFSysSigRFB080to193");
    hHFRFB0[3] = (TH1D*)f3->Get("hHFSysSigRFB000to193");
    hNtracksRFB1[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to080");
    hNtracksRFB2[3] = (TH1D*)f3->Get("hNtracksSysSigRFB080to193");
    hNtracksRFB0[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to193");
  }
  else if(state==3){
    hHFRFB1[3] = (TH1D*)f3->Get("hHFSysSigRFB000to193");
    hHFRFB0[3] = (TH1D*)f3->Get("hHFSysSigRFB000to193");
    hNtracksRFB1[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to193");
    hNtracksRFB0[3] = (TH1D*)f3->Get("hNtracksSysSigRFB000to193");
  }
  
  
  
  // 4 : background PDF 
  TFile* f4 = new TFile("Graham_BkgVariation/BkgPdfSystematics.root");
  if(state==1){
    hHFRFB1[4] = (TH1D*)f4->Get("hhf1S_y1_yRfbdiff");
    hHFRFB2[4] = (TH1D*)f4->Get("hhf1S_y2_yRfbdiff");
    hHFRFB3[4] = (TH1D*)f4->Get("hhf1S_y3_yRfbdiff");
    hHFRFB4[4] = (TH1D*)f4->Get("hnt1S_y4_yRfbdiff");
    hHFRFB0[4] = (TH1D*)f4->Get("hhf1S_yInt_Rfbdiff");
    hNtracksRFB1[4] = (TH1D*)f4->Get("hnt1S_y1_yRfbdiff");
    hNtracksRFB2[4] = (TH1D*)f4->Get("hnt1S_y2_yRfbdiff");
    hNtracksRFB3[4] = (TH1D*)f4->Get("hnt1S_y3_yRfbdiff");
    hNtracksRFB4[4] = (TH1D*)f4->Get("hnt1S_y4_yRfbdiff");
    hNtracksRFB0[4] = (TH1D*)f4->Get("hnt1S_yInt_Rfbdiff");
  }
  else if(state==2){
    hHFRFB1[4] = (TH1D*)f4->Get("hhf2S_y1_yRfbdiff");
    hHFRFB2[4] = (TH1D*)f4->Get("hhf2S_y2_yRfbdiff");
    hHFRFB0[4] = (TH1D*)f4->Get("hhf2S_yInt_Rfbdiff");
    hNtracksRFB1[4] = (TH1D*)f4->Get("hnt2S_y1_yRfbdiff");
    hNtracksRFB2[4] = (TH1D*)f4->Get("hnt2S_y2_yRfbdiff");
    hNtracksRFB0[4] = (TH1D*)f4->Get("hnt2S_yInt_Rfbdiff");
  }
  else if(state==3){
    hHFRFB1[4] = (TH1D*)f4->Get("hhf3S_y1_yRfbdiff");
    hHFRFB0[4] = (TH1D*)f4->Get("hhf3S_yInt_Rfbdiff");
    hNtracksRFB1[4] = (TH1D*)f4->Get("hnt3S_y1_yRfbdiff");
    hNtracksRFB0[4] = (TH1D*)f4->Get("hnt3S_yInt_Rfbdiff");
  }

  // Merge uncertainties for cross-section 
  if(state==1){
    mergeFourInQuad(hHFRFB1[0],hHFRFB1[1],hHFRFB1[2],hHFRFB1[3],hHFRFB1[4],state);
    mergeFourInQuad(hHFRFB2[0],hHFRFB2[1],hHFRFB2[2],hHFRFB2[3],hHFRFB2[4],state);
    mergeFourInQuad(hHFRFB3[0],hHFRFB3[1],hHFRFB3[2],hHFRFB3[3],hHFRFB3[4],state);
    mergeFourInQuad(hHFRFB4[0],hHFRFB4[1],hHFRFB4[2],hHFRFB4[3],hHFRFB4[4],state);
    mergeFourInQuad(hHFRFB0[0],hHFRFB0[1],hHFRFB0[2],hHFRFB0[3],hHFRFB0[4],state);
    mergeFourInQuad(hNtracksRFB1[0],hNtracksRFB1[1],hNtracksRFB1[2],hNtracksRFB1[3],hNtracksRFB1[4],state);
    mergeFourInQuad(hNtracksRFB2[0],hNtracksRFB2[1],hNtracksRFB2[2],hNtracksRFB2[3],hNtracksRFB2[4],state);
    mergeFourInQuad(hNtracksRFB3[0],hNtracksRFB3[1],hNtracksRFB3[2],hNtracksRFB3[3],hNtracksRFB3[4],state);
    mergeFourInQuad(hNtracksRFB4[0],hNtracksRFB4[1],hNtracksRFB4[2],hNtracksRFB4[3],hNtracksRFB4[4],state);
    mergeFourInQuad(hNtracksRFB0[0],hNtracksRFB0[1],hNtracksRFB0[2],hNtracksRFB0[3],hNtracksRFB0[4],state);
  }
  else if(state==2){
    mergeFourInQuad(hHFRFB1[0],hHFRFB1[1],hHFRFB1[2],hHFRFB1[3],hHFRFB1[4],state);
    mergeFourInQuad(hHFRFB2[0],hHFRFB2[1],hHFRFB2[2],hHFRFB2[3],hHFRFB2[4],state);
    mergeFourInQuad(hHFRFB0[0],hHFRFB0[1],hHFRFB0[2],hHFRFB0[3],hHFRFB0[4],state);
    mergeFourInQuad(hNtracksRFB1[0],hNtracksRFB1[1],hNtracksRFB1[2],hNtracksRFB1[3],hNtracksRFB1[4],state);
    mergeFourInQuad(hNtracksRFB2[0],hNtracksRFB2[1],hNtracksRFB2[2],hNtracksRFB2[3],hNtracksRFB2[4],state);
    mergeFourInQuad(hNtracksRFB0[0],hNtracksRFB0[1],hNtracksRFB0[2],hNtracksRFB0[3],hNtracksRFB0[4],state);
  }
  else if(state==3){
    mergeFourInQuad(hHFRFB1[0],hHFRFB1[1],hHFRFB1[2],hHFRFB1[3],hHFRFB1[4],state);
    mergeFourInQuad(hHFRFB0[0],hHFRFB0[1],hHFRFB0[2],hHFRFB0[3],hHFRFB0[4],state);
    mergeFourInQuad(hNtracksRFB1[0],hNtracksRFB1[1],hNtracksRFB1[2],hNtracksRFB1[3],hNtracksRFB1[4],state);
    mergeFourInQuad(hNtracksRFB0[0],hNtracksRFB0[1],hNtracksRFB0[2],hNtracksRFB0[3],hNtracksRFB0[4],state);
  }

  TFile* fout = new TFile(Form("mergedSys_ups%ds_rfb.root",state),"recreate" );
  hHFRFB1[0]->SetName("hHFmerged1");
  hHFRFB2[0]->SetName("hHFmerged2");
  hHFRFB3[0]->SetName("hHFmerged3");
  hHFRFB4[0]->SetName("hHFmerged4");
  hHFRFB0[0]->SetName("hHFmergedint");
  hNtracksRFB1[0]->SetName("hNtracksmerged1");
  hNtracksRFB2[0]->SetName("hNtracksmerged2");
  hNtracksRFB3[0]->SetName("hNtracksmerged3");
  hNtracksRFB4[0]->SetName("hNtracksmerged4");
  hNtracksRFB0[0]->SetName("hNtracksmergedint");
  
  if(state==1){
    hHFRFB1[0]->Write();
    hHFRFB2[0]->Write();
    hHFRFB3[0]->Write();
    hHFRFB4[0]->Write();
    hHFRFB0[0]->Write();
    hNtracksRFB1[0]->Write();
    hNtracksRFB2[0]->Write();
    hNtracksRFB3[0]->Write();
    hNtracksRFB4[0]->Write();
    hNtracksRFB0[0]->Write();
  }
  else if(state==2){
    hHFRFB1[0]->Write();
    hHFRFB2[0]->Write();
    hHFRFB0[0]->Write();
    hNtracksRFB1[0]->Write();
    hNtracksRFB2[0]->Write();
    hNtracksRFB0[0]->Write();
  }
  else if(state==3){
    hHFRFB1[0]->Write();
    hHFRFB0[0]->Write();
    hNtracksRFB1[0]->Write();
    hNtracksRFB0[0]->Write();
  }
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
//  c0->SaveAs(Form("pdfFiles/%s_ups%ds.pdf", h0->GetName(),state ) );
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
