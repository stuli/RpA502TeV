#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
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
void getMaxTH1D_eight ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0,TH1D* h5=0,TH1D* h6=0, TH1D* h7=0,TH1D* h8=0) ;
void getMaxTH1D_six ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0,TH1D* h5=0,TH1D* h6=0) ;
void getMaxTH1D_four ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0) ;
void addInQuadFive ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0);
void addInQuadSix ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0);
void addInQuadNine ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0, TH1D* h7=0, TH1D* h8=0, TH1D* h9=0);
void addInQuadEight ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0, TH1D* h7=0, TH1D* h8=0);
void getEffSys(int state =1, int Nsamples=100) { 
  TH1::SetDefaultSumw2();

  TFile* f1 = new TFile(Form("efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
  TH1D* eff1 = (TH1D*)f1->Get("hptEffPP");
  TH1D* eff2 = (TH1D*)f1->Get("hptEffAA");
  TH1D* eff3 = (TH1D*)f1->Get("hrapEffPP");
  TH1D* eff4 = (TH1D*)f1->Get("hrapEffAA");
  TH1D* eff5 = (TH1D*)f1->Get("hcentintEffPP");
  TH1D* eff6 = (TH1D*)f1->Get("hcentintEffAA");
  TH1D* eff7 = (TH1D*)f1->Get("hcentEffAA");

  TH1D* eff1sys = (TH1D*)eff1->Clone("hptEffPPSys"); eff1sys->Reset();
  TH1D* eff2sys = (TH1D*)eff2->Clone("hptEffAASys"); eff2sys->Reset();
  TH1D* eff3sys = (TH1D*)eff3->Clone("hrapEffPPSys"); eff3sys->Reset();
  TH1D* eff4sys = (TH1D*)eff4->Clone("hrapEffAASys"); eff4sys->Reset();
  TH1D* eff5sys = (TH1D*)eff5->Clone("hcentintEffPPSys"); eff5sys->Reset();
  TH1D* eff6sys = (TH1D*)eff6->Clone("hcentintEffAASys"); eff6sys->Reset();
  TH1D* eff7sys = (TH1D*)eff7->Clone("hcentEffAASys"); eff7sys->Reset();
  // stat. fluc. Trig
  TH1D* eff1statTrg = (TH1D*)eff1sys->Clone("hptEffPPStatTrg"); 
  TH1D* eff2statTrg = (TH1D*)eff2sys->Clone("hptEffAAStatTrg"); 
  TH1D* eff3statTrg = (TH1D*)eff3sys->Clone("hrapEffPPStatTrg");
  TH1D* eff4statTrg = (TH1D*)eff4sys->Clone("hrapEffAAStatTrg");
  TH1D* eff5statTrg = (TH1D*)eff5sys->Clone("hcentintEffPPStatTrg"); 
  TH1D* eff6statTrg = (TH1D*)eff6sys->Clone("hcentintEffAAStatTrg"); 
  TH1D* eff7statTrg = (TH1D*)eff7sys->Clone("hcentEffAAStatTrg");
  // stat. fluc. MuId
  TH1D* eff1statMuId = (TH1D*)eff1sys->Clone("hptEffPPStatMuId"); 
  TH1D* eff2statMuId = (TH1D*)eff2sys->Clone("hptEffAAStatMuId"); 
  TH1D* eff3statMuId = (TH1D*)eff3sys->Clone("hrapEffPPStatMuId");
  TH1D* eff4statMuId = (TH1D*)eff4sys->Clone("hrapEffAAStatMuId");
  TH1D* eff5statMuId = (TH1D*)eff5sys->Clone("hcentintEffPPStatMuId"); 
  TH1D* eff6statMuId = (TH1D*)eff6sys->Clone("hcentintEffAAStatMuId"); 
  TH1D* eff7statMuId = (TH1D*)eff7sys->Clone("hcentEffAAStatMuId");
  // stat. fluc. STA
  TH1D* eff1statSta = (TH1D*)eff1sys->Clone("hptEffPPStatSta"); 
  TH1D* eff2statSta = (TH1D*)eff2sys->Clone("hptEffAAStatSta"); 
  TH1D* eff3statSta = (TH1D*)eff3sys->Clone("hrapEffPPStatSta");
  TH1D* eff4statSta = (TH1D*)eff4sys->Clone("hrapEffAAStatSta");
  TH1D* eff5statSta = (TH1D*)eff5sys->Clone("hcentintEffPPStatSta"); 
  TH1D* eff6statSta = (TH1D*)eff6sys->Clone("hcentintEffAAStatSta"); 
  TH1D* eff7statSta = (TH1D*)eff7sys->Clone("hcentEffAAStatSta");



  ///////////////////// sys. trig variation
  TH1D* eff1TrgVar = (TH1D*)eff1sys->Clone("hptEffPPTrgVar"); 
  TH1D* eff2TrgVar = (TH1D*)eff2sys->Clone("hptEffAATrgVar"); 
  TH1D* eff3TrgVar = (TH1D*)eff3sys->Clone("hrapEffPPTrgVar");
  TH1D* eff4TrgVar = (TH1D*)eff4sys->Clone("hrapEffAATrgVar");
  TH1D* eff5TrgVar = (TH1D*)eff5sys->Clone("hcentintEffPPTrgVar"); 
  TH1D* eff6TrgVar = (TH1D*)eff6sys->Clone("hcentintEffAATrgVar"); 
  TH1D* eff7TrgVar = (TH1D*)eff7sys->Clone("hcentEffAATrgVar");
  // sys.var_1
  TFile* ftrg_1 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId-1_trkId0_muId-100_staId-100.root",state));
  TH1D* relTrgVar1_1 = (TH1D*)ftrg_1->Get("hptEffPP");
  TH1D* relTrgVar2_1 = (TH1D*)ftrg_1->Get("hptEffAA");
  TH1D* relTrgVar3_1 = (TH1D*)ftrg_1->Get("hrapEffPP");
  TH1D* relTrgVar4_1 = (TH1D*)ftrg_1->Get("hrapEffAA");
  TH1D* relTrgVar5_1 = (TH1D*)ftrg_1->Get("hcentintEffPP");
  TH1D* relTrgVar6_1 = (TH1D*)ftrg_1->Get("hcentintEffAA");
  TH1D* relTrgVar7_1 = (TH1D*)ftrg_1->Get("hcentEffAA");
  relTrgVar1_1->Add( eff1, -1 );     relTrgVar1_1->Divide( eff1);
  relTrgVar2_1->Add( eff2, -1 );     relTrgVar2_1->Divide( eff2);
  relTrgVar3_1->Add( eff3, -1 );     relTrgVar3_1->Divide( eff3);
  relTrgVar4_1->Add( eff4, -1 );     relTrgVar4_1->Divide( eff4);
  relTrgVar5_1->Add( eff5, -1 );     relTrgVar5_1->Divide( eff5);
  relTrgVar6_1->Add( eff6, -1 );     relTrgVar6_1->Divide( eff6);
  relTrgVar7_1->Add( eff7, -1 );     relTrgVar7_1->Divide( eff7);
  // sys.var_2
  TFile* ftrg_2 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId-2_trkId0_muId-100_staId-100.root",state));
  TH1D* relTrgVar1_2 = (TH1D*)ftrg_2->Get("hptEffPP");
  TH1D* relTrgVar2_2 = (TH1D*)ftrg_2->Get("hptEffAA");
  TH1D* relTrgVar3_2 = (TH1D*)ftrg_2->Get("hrapEffPP");
  TH1D* relTrgVar4_2 = (TH1D*)ftrg_2->Get("hrapEffAA");
  TH1D* relTrgVar5_2 = (TH1D*)ftrg_2->Get("hcentintEffPP");
  TH1D* relTrgVar6_2 = (TH1D*)ftrg_2->Get("hcentintEffAA");
  TH1D* relTrgVar7_2 = (TH1D*)ftrg_2->Get("hcentEffAA");
  relTrgVar1_2->Add( eff1, -1 );     relTrgVar1_2->Divide( eff1);
  relTrgVar2_2->Add( eff2, -1 );     relTrgVar2_2->Divide( eff2);
  relTrgVar3_2->Add( eff3, -1 );     relTrgVar3_2->Divide( eff3);
  relTrgVar4_2->Add( eff4, -1 );     relTrgVar4_2->Divide( eff4);
  relTrgVar5_2->Add( eff5, -1 );     relTrgVar5_2->Divide( eff5);
  relTrgVar6_2->Add( eff6, -1 );     relTrgVar6_2->Divide( eff6);
  relTrgVar7_2->Add( eff7, -1 );     relTrgVar7_2->Divide( eff7);
  // get Max deviation of each bin:
  getMaxTH1D(eff1TrgVar, relTrgVar1_1, relTrgVar1_2);
  getMaxTH1D(eff2TrgVar, relTrgVar2_1, relTrgVar2_2);
  getMaxTH1D(eff3TrgVar, relTrgVar3_1, relTrgVar3_2);
  getMaxTH1D(eff4TrgVar, relTrgVar4_1, relTrgVar4_2);
  getMaxTH1D(eff5TrgVar, relTrgVar5_1, relTrgVar5_2);
  getMaxTH1D(eff6TrgVar, relTrgVar6_1, relTrgVar6_2);
  getMaxTH1D(eff7TrgVar, relTrgVar7_1, relTrgVar7_2);



  ///////////////////// sys. muid variation
  TH1D* eff1MuIdVar = (TH1D*)eff1sys->Clone("hptEffPPMuIdVar"); 
  TH1D* eff2MuIdVar = (TH1D*)eff2sys->Clone("hptEffAAMuIdVar"); 
  TH1D* eff3MuIdVar = (TH1D*)eff3sys->Clone("hrapEffPPMuIdVar");
  TH1D* eff4MuIdVar = (TH1D*)eff4sys->Clone("hrapEffAAMuIdVar");
  TH1D* eff5MuIdVar = (TH1D*)eff5sys->Clone("hcentintEffPPMuIdVar"); 
  TH1D* eff6MuIdVar = (TH1D*)eff6sys->Clone("hcentintEffAAMuIdVar"); 
  TH1D* eff7MuIdVar = (TH1D*)eff7sys->Clone("hcentEffAAMuIdVar");
  // sys.var_1
  TFile* fmuId_1 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-1_staId-100.root",state));
  TH1D* relMuIdVar1_1 = (TH1D*)fmuId_1->Get("hptEffPP");
  TH1D* relMuIdVar2_1 = (TH1D*)fmuId_1->Get("hptEffAA");
  TH1D* relMuIdVar3_1 = (TH1D*)fmuId_1->Get("hrapEffPP");
  TH1D* relMuIdVar4_1 = (TH1D*)fmuId_1->Get("hrapEffAA");
  TH1D* relMuIdVar5_1 = (TH1D*)fmuId_1->Get("hcentintEffPP");
  TH1D* relMuIdVar6_1 = (TH1D*)fmuId_1->Get("hcentintEffAA");
  TH1D* relMuIdVar7_1 = (TH1D*)fmuId_1->Get("hcentEffAA");
  relMuIdVar1_1->Add( eff1, -1 );     relMuIdVar1_1->Divide( eff1);
  relMuIdVar2_1->Add( eff2, -1 );     relMuIdVar2_1->Divide( eff2);
  relMuIdVar3_1->Add( eff3, -1 );     relMuIdVar3_1->Divide( eff3);
  relMuIdVar4_1->Add( eff4, -1 );     relMuIdVar4_1->Divide( eff4);
  relMuIdVar5_1->Add( eff5, -1 );     relMuIdVar5_1->Divide( eff5);
  relMuIdVar6_1->Add( eff6, -1 );     relMuIdVar6_1->Divide( eff6);
  relMuIdVar7_1->Add( eff7, -1 );     relMuIdVar7_1->Divide( eff7);
  // sys.var_2
  TFile* fmuId_2 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-2_staId-100.root",state));
  TH1D* relMuIdVar1_2 = (TH1D*)fmuId_2->Get("hptEffPP");
  TH1D* relMuIdVar2_2 = (TH1D*)fmuId_2->Get("hptEffAA");
  TH1D* relMuIdVar3_2 = (TH1D*)fmuId_2->Get("hrapEffPP");
  TH1D* relMuIdVar4_2 = (TH1D*)fmuId_2->Get("hrapEffAA");
  TH1D* relMuIdVar5_2 = (TH1D*)fmuId_2->Get("hcentintEffPP");
  TH1D* relMuIdVar6_2 = (TH1D*)fmuId_2->Get("hcentintEffAA");
  TH1D* relMuIdVar7_2 = (TH1D*)fmuId_2->Get("hcentEffAA");
  relMuIdVar1_2->Add( eff1, -1 );     relMuIdVar1_2->Divide( eff1);
  relMuIdVar2_2->Add( eff2, -1 );     relMuIdVar2_2->Divide( eff2);
  relMuIdVar3_2->Add( eff3, -1 );     relMuIdVar3_2->Divide( eff3);
  relMuIdVar4_2->Add( eff4, -1 );     relMuIdVar4_2->Divide( eff4);
  relMuIdVar5_2->Add( eff5, -1 );     relMuIdVar5_2->Divide( eff5);
  relMuIdVar6_2->Add( eff6, -1 );     relMuIdVar6_2->Divide( eff6);
  relMuIdVar7_2->Add( eff7, -1 );     relMuIdVar7_2->Divide( eff7);
  // get Max deviation of each bin:
  getMaxTH1D(eff1MuIdVar, relMuIdVar1_1, relMuIdVar1_2);
  getMaxTH1D(eff2MuIdVar, relMuIdVar2_1, relMuIdVar2_2);
  getMaxTH1D(eff3MuIdVar, relMuIdVar3_1, relMuIdVar3_2);
  getMaxTH1D(eff4MuIdVar, relMuIdVar4_1, relMuIdVar4_2);
  getMaxTH1D(eff5MuIdVar, relMuIdVar5_1, relMuIdVar5_2);
  getMaxTH1D(eff6MuIdVar, relMuIdVar6_1, relMuIdVar6_2);
  getMaxTH1D(eff7MuIdVar, relMuIdVar7_1, relMuIdVar7_2);

  ///////////////////// sys. sta variation
  TH1D* eff1StaVar = (TH1D*)eff1sys->Clone("hptEffPPStaVar"); 
  TH1D* eff2StaVar = (TH1D*)eff2sys->Clone("hptEffAAStaVar"); 
  TH1D* eff3StaVar = (TH1D*)eff3sys->Clone("hrapEffPPStaVar");
  TH1D* eff4StaVar = (TH1D*)eff4sys->Clone("hrapEffAAStaVar");
  TH1D* eff5StaVar = (TH1D*)eff5sys->Clone("hcentintEffPPStaVar"); 
  TH1D* eff6StaVar = (TH1D*)eff6sys->Clone("hcentintEffAAStaVar"); 
  TH1D* eff7StaVar = (TH1D*)eff7sys->Clone("hcentEffAAStaVar");
  // sys.var_1
  TFile* fsta_1 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-1.root",state));
  TH1D* relStaVar1_1 = (TH1D*)fsta_1->Get("hptEffPP");
  TH1D* relStaVar2_1 = (TH1D*)fsta_1->Get("hptEffAA");
  TH1D* relStaVar3_1 = (TH1D*)fsta_1->Get("hrapEffPP");
  TH1D* relStaVar4_1 = (TH1D*)fsta_1->Get("hrapEffAA");
  TH1D* relStaVar5_1 = (TH1D*)fsta_1->Get("hcentintEffPP");
  TH1D* relStaVar6_1 = (TH1D*)fsta_1->Get("hcentintEffAA");
  TH1D* relStaVar7_1 = (TH1D*)fsta_1->Get("hcentEffAA");
  relStaVar1_1->Add( eff1, -1 );     relStaVar1_1->Divide( eff1);
  relStaVar2_1->Add( eff2, -1 );     relStaVar2_1->Divide( eff2);
  relStaVar3_1->Add( eff3, -1 );     relStaVar3_1->Divide( eff3);
  relStaVar4_1->Add( eff4, -1 );     relStaVar4_1->Divide( eff4);
  relStaVar5_1->Add( eff5, -1 );     relStaVar5_1->Divide( eff5);
  relStaVar6_1->Add( eff6, -1 );     relStaVar6_1->Divide( eff6);
  relStaVar7_1->Add( eff7, -1 );     relStaVar7_1->Divide( eff7);
  // sys.var_2
  TFile* fsta_2 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-2.root",state));
  TH1D* relStaVar1_2 = (TH1D*)fsta_2->Get("hptEffPP");
  TH1D* relStaVar2_2 = (TH1D*)fsta_2->Get("hptEffAA");
  TH1D* relStaVar3_2 = (TH1D*)fsta_2->Get("hrapEffPP");
  TH1D* relStaVar4_2 = (TH1D*)fsta_2->Get("hrapEffAA");
  TH1D* relStaVar5_2 = (TH1D*)fsta_2->Get("hcentintEffPP");
  TH1D* relStaVar6_2 = (TH1D*)fsta_2->Get("hcentintEffAA");
  TH1D* relStaVar7_2 = (TH1D*)fsta_2->Get("hcentEffAA");
  relStaVar1_2->Add( eff1, -1 );     relStaVar1_2->Divide( eff1);
  relStaVar2_2->Add( eff2, -1 );     relStaVar2_2->Divide( eff2);
  relStaVar3_2->Add( eff3, -1 );     relStaVar3_2->Divide( eff3);
  relStaVar4_2->Add( eff4, -1 );     relStaVar4_2->Divide( eff4);
  relStaVar5_2->Add( eff5, -1 );     relStaVar5_2->Divide( eff5);
  relStaVar6_2->Add( eff6, -1 );     relStaVar6_2->Divide( eff6);
  relStaVar7_2->Add( eff7, -1 );     relStaVar7_2->Divide( eff7);
  // get Max deviation of each bin:
  getMaxTH1D(eff1StaVar, relStaVar1_1, relStaVar1_2);
  getMaxTH1D(eff2StaVar, relStaVar2_1, relStaVar2_2);
  getMaxTH1D(eff3StaVar, relStaVar3_1, relStaVar3_2);
  getMaxTH1D(eff4StaVar, relStaVar4_1, relStaVar4_2);
  getMaxTH1D(eff5StaVar, relStaVar5_1, relStaVar5_2);
  getMaxTH1D(eff6StaVar, relStaVar6_1, relStaVar6_2);
  getMaxTH1D(eff7StaVar, relStaVar7_1, relStaVar7_2);


  ///////////////////// sys. trk variation
  TH1D* eff1TrkVar = (TH1D*)eff1sys->Clone("hptEffPPTrkVar"); 
  TH1D* eff2TrkVar = (TH1D*)eff2sys->Clone("hptEffAATrkVar"); 
  TH1D* eff3TrkVar = (TH1D*)eff3sys->Clone("hrapEffPPTrkVar");
  TH1D* eff4TrkVar = (TH1D*)eff4sys->Clone("hrapEffAATrkVar");
  TH1D* eff5TrkVar = (TH1D*)eff5sys->Clone("hcentintEffPPTrkVar"); 
  TH1D* eff6TrkVar = (TH1D*)eff6sys->Clone("hcentintEffAATrkVar"); 
  TH1D* eff7TrkVar = (TH1D*)eff7sys->Clone("hcentEffAATrkVar");
  // sys.var_1
  TFile* ftrk_1 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId-1_muId-100_staId-100.root",state));
  TH1D* relTrkVar1_1 = (TH1D*)ftrk_1->Get("hptEffPP");
  TH1D* relTrkVar2_1 = (TH1D*)ftrk_1->Get("hptEffAA");
  TH1D* relTrkVar3_1 = (TH1D*)ftrk_1->Get("hrapEffPP");
  TH1D* relTrkVar4_1 = (TH1D*)ftrk_1->Get("hrapEffAA");
  TH1D* relTrkVar5_1 = (TH1D*)ftrk_1->Get("hcentintEffPP");
  TH1D* relTrkVar6_1 = (TH1D*)ftrk_1->Get("hcentintEffAA");
  TH1D* relTrkVar7_1 = (TH1D*)ftrk_1->Get("hcentEffAA");
  relTrkVar1_1->Add( eff1, -1 );     relTrkVar1_1->Divide( eff1);
  relTrkVar2_1->Add( eff2, -1 );     relTrkVar2_1->Divide( eff2);
  relTrkVar3_1->Add( eff3, -1 );     relTrkVar3_1->Divide( eff3);
  relTrkVar4_1->Add( eff4, -1 );     relTrkVar4_1->Divide( eff4);
  relTrkVar5_1->Add( eff5, -1 );     relTrkVar5_1->Divide( eff5);
  relTrkVar6_1->Add( eff6, -1 );     relTrkVar6_1->Divide( eff6);
  relTrkVar7_1->Add( eff7, -1 );     relTrkVar7_1->Divide( eff7);
  // sys.var_2
  TFile* ftrk_2 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId-2_muId-100_staId-100.root",state));
  TH1D* relTrkVar1_2 = (TH1D*)ftrk_2->Get("hptEffPP");
  TH1D* relTrkVar2_2 = (TH1D*)ftrk_2->Get("hptEffAA");
  TH1D* relTrkVar3_2 = (TH1D*)ftrk_2->Get("hrapEffPP");
  TH1D* relTrkVar4_2 = (TH1D*)ftrk_2->Get("hrapEffAA");
  TH1D* relTrkVar5_2 = (TH1D*)ftrk_2->Get("hcentintEffPP");
  TH1D* relTrkVar6_2 = (TH1D*)ftrk_2->Get("hcentintEffAA");
  TH1D* relTrkVar7_2 = (TH1D*)ftrk_2->Get("hcentEffAA");
  relTrkVar1_2->Add( eff1, -1 );     relTrkVar1_2->Divide( eff1);
  relTrkVar2_2->Add( eff2, -1 );     relTrkVar2_2->Divide( eff2);
  relTrkVar3_2->Add( eff3, -1 );     relTrkVar3_2->Divide( eff3);
  relTrkVar4_2->Add( eff4, -1 );     relTrkVar4_2->Divide( eff4);
  relTrkVar5_2->Add( eff5, -1 );     relTrkVar5_2->Divide( eff5);
  relTrkVar6_2->Add( eff6, -1 );     relTrkVar6_2->Divide( eff6);
  relTrkVar7_2->Add( eff7, -1 );     relTrkVar7_2->Divide( eff7);
  // get Max deviation of each bin:
  getMaxTH1D(eff1TrkVar, relTrkVar1_1, relTrkVar1_2);
  getMaxTH1D(eff2TrkVar, relTrkVar2_1, relTrkVar2_2);
  getMaxTH1D(eff3TrkVar, relTrkVar3_1, relTrkVar3_2);
  getMaxTH1D(eff4TrkVar, relTrkVar4_1, relTrkVar4_2);
  getMaxTH1D(eff5TrkVar, relTrkVar5_1, relTrkVar5_2);
  getMaxTH1D(eff6TrkVar, relTrkVar6_1, relTrkVar6_2);
  getMaxTH1D(eff7TrkVar, relTrkVar7_1, relTrkVar7_2);

  // SF from binnined table:
  TFile* fid_10 = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId-10_trkId0_muId-100_staId-100.root",state ) );
  TH1D* eff1binned = (TH1D*)fid_10->Get("hptEffPP");
  TH1D* eff2binned = (TH1D*)fid_10->Get("hptEffAA");
  TH1D* eff3binned = (TH1D*)fid_10->Get("hrapEffPP");
  TH1D* eff4binned = (TH1D*)fid_10->Get("hrapEffAA");
  TH1D* eff5binned = (TH1D*)fid_10->Get("hcentintEffPP");
  TH1D* eff6binned = (TH1D*)fid_10->Get("hcentintEffAA");
  TH1D* eff7binned = (TH1D*)fid_10->Get("hcentEffAA");
  eff1binned ->SetName("hptRelSysBinnedPP");
  eff2binned ->SetName("hptRelSysBinnedAA");
  eff3binned ->SetName("hrapRelSysBinnedPP");
  eff4binned ->SetName("hrapRelSysBinnedAA");
  eff5binned ->SetName("hcentintRelSysBinnedPP");
  eff6binned ->SetName("hcentintRelSysBinnedAA");
  eff7binned ->SetName("hcentRelSysBinnedAA");
  eff1binned->Add( eff1, -1 );     eff1binned->Divide( eff1);
  eff2binned->Add( eff2, -1 );     eff2binned->Divide( eff2);
  eff3binned->Add( eff3, -1 );     eff3binned->Divide( eff3);
  eff4binned->Add( eff4, -1 );     eff4binned->Divide( eff4);
  eff5binned->Add( eff5, -1 );     eff5binned->Divide( eff5);
  eff6binned->Add( eff6, -1 );     eff6binned->Divide( eff6);
  eff7binned->Add( eff7, -1 );     eff7binned->Divide( eff7);
  
  
  // pT reweight (nothing to do with TNP
  
  TH1D* eff1_Cup; 
  TH1D* eff2_Cup; 
  TH1D* eff3_Cup; 
  TH1D* eff4_Cup; 
  TH1D* eff5_Cup; 
  TH1D* eff6_Cup; 
  TH1D* eff7_Cup; 

  TH1D* eff1_Dup; 
  TH1D* eff2_Dup; 
  TH1D* eff3_Dup; 
  TH1D* eff4_Dup; 
  TH1D* eff5_Dup; 
  TH1D* eff6_Dup; 
  TH1D* eff7_Dup; 

  TH1D* eff1_Cdo; 
  TH1D* eff2_Cdo; 
  TH1D* eff3_Cdo; 
  TH1D* eff4_Cdo; 
  TH1D* eff5_Cdo; 
  TH1D* eff6_Cdo; 
  TH1D* eff7_Cdo; 

  TH1D* eff1_Ddo; 
  TH1D* eff2_Ddo; 
  TH1D* eff3_Ddo; 
  TH1D* eff4_Ddo; 
  TH1D* eff5_Ddo; 
  TH1D* eff6_Ddo; 
  TH1D* eff7_Ddo; 


  TFile* fAup = new TFile(Form("efficiencyTable/VarAup_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
  TFile* fAdo = new TFile(Form("efficiencyTable/VarAdo_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
  TFile* fBup = new TFile(Form("efficiencyTable/VarBup_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
  TFile* fBdo = new TFile(Form("efficiencyTable/VarBdo_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
  TFile* fCup;
  TFile* fDup;
  TFile* fCdo;
  TFile* fDdo;
  if(state == 1)
  {
    fCup = new TFile(Form("efficiencyTable/VarCup_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
    fCdo = new TFile(Form("efficiencyTable/VarCdo_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
    fDup = new TFile(Form("efficiencyTable/VarDup_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
    fDdo = new TFile(Form("efficiencyTable/VarDdo_efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state) );
    eff1_Cup = (TH1D*)fCup->Get("hptEffPP");
    eff2_Cup = (TH1D*)fCup->Get("hptEffAA");
    eff3_Cup = (TH1D*)fCup->Get("hrapEffPP");
    eff4_Cup = (TH1D*)fCup->Get("hrapEffAA");
    eff5_Cup = (TH1D*)fCup->Get("hcentintEffPP");
    eff6_Cup = (TH1D*)fCup->Get("hcentintEffAA");
    eff7_Cup = (TH1D*)fCup->Get("hcentEffAA");

    eff1_Dup = (TH1D*)fDup->Get("hptEffPP");
    eff2_Dup = (TH1D*)fDup->Get("hptEffAA");
    eff3_Dup = (TH1D*)fDup->Get("hrapEffPP");
    eff4_Dup = (TH1D*)fDup->Get("hrapEffAA");
    eff5_Dup = (TH1D*)fDup->Get("hcentintEffPP");
    eff6_Dup = (TH1D*)fDup->Get("hcentintEffAA");
    eff7_Dup = (TH1D*)fDup->Get("hcentEffAA");

    eff1_Cdo = (TH1D*)fCdo->Get("hptEffPP");
    eff2_Cdo = (TH1D*)fCdo->Get("hptEffAA");
    eff3_Cdo = (TH1D*)fCdo->Get("hrapEffPP");
    eff4_Cdo = (TH1D*)fCdo->Get("hrapEffAA");
    eff5_Cdo = (TH1D*)fCdo->Get("hcentintEffPP");
    eff6_Cdo = (TH1D*)fCdo->Get("hcentintEffAA");
    eff7_Cdo = (TH1D*)fCdo->Get("hcentEffAA");

    eff1_Ddo = (TH1D*)fDdo->Get("hptEffPP");
    eff2_Ddo = (TH1D*)fDdo->Get("hptEffAA");
    eff3_Ddo = (TH1D*)fDdo->Get("hrapEffPP");
    eff4_Ddo = (TH1D*)fDdo->Get("hrapEffAA");
    eff5_Ddo = (TH1D*)fDdo->Get("hcentintEffPP");
    eff6_Ddo = (TH1D*)fDdo->Get("hcentintEffAA");
    eff7_Ddo = (TH1D*)fDdo->Get("hcentEffAA");
  
    eff1_Cup->Add( eff1, -1 );     eff1_Cup->Divide( eff1);
    eff2_Cup->Add( eff2, -1 );     eff2_Cup->Divide( eff2);
    eff3_Cup->Add( eff3, -1 );     eff3_Cup->Divide( eff3);
    eff4_Cup->Add( eff4, -1 );     eff4_Cup->Divide( eff4);
    eff5_Cup->Add( eff5, -1 );     eff5_Cup->Divide( eff5);
    eff6_Cup->Add( eff6, -1 );     eff6_Cup->Divide( eff6);
    eff7_Cup->Add( eff7, -1 );     eff7_Cup->Divide( eff7);

    eff1_Dup->Add( eff1, -1 );     eff1_Dup->Divide( eff1);
    eff2_Dup->Add( eff2, -1 );     eff2_Dup->Divide( eff2);
    eff3_Dup->Add( eff3, -1 );     eff3_Dup->Divide( eff3);
    eff4_Dup->Add( eff4, -1 );     eff4_Dup->Divide( eff4);
    eff5_Dup->Add( eff5, -1 );     eff5_Dup->Divide( eff5);
    eff6_Dup->Add( eff6, -1 );     eff6_Dup->Divide( eff6);
    eff7_Dup->Add( eff7, -1 );     eff7_Dup->Divide( eff7);

    eff1_Cdo->Add( eff1, -1 );     eff1_Cdo->Divide( eff1);
    eff2_Cdo->Add( eff2, -1 );     eff2_Cdo->Divide( eff2);
    eff3_Cdo->Add( eff3, -1 );     eff3_Cdo->Divide( eff3);
    eff4_Cdo->Add( eff4, -1 );     eff4_Cdo->Divide( eff4);
    eff5_Cdo->Add( eff5, -1 );     eff5_Cdo->Divide( eff5);
    eff6_Cdo->Add( eff6, -1 );     eff6_Cdo->Divide( eff6);
    eff7_Cdo->Add( eff7, -1 );     eff7_Cdo->Divide( eff7);

    eff1_Ddo->Add( eff1, -1 );     eff1_Ddo->Divide( eff1);
    eff2_Ddo->Add( eff2, -1 );     eff2_Ddo->Divide( eff2);
    eff3_Ddo->Add( eff3, -1 );     eff3_Ddo->Divide( eff3);
    eff4_Ddo->Add( eff4, -1 );     eff4_Ddo->Divide( eff4);
    eff5_Ddo->Add( eff5, -1 );     eff5_Ddo->Divide( eff5);
    eff6_Ddo->Add( eff6, -1 );     eff6_Ddo->Divide( eff6);
    eff7_Ddo->Add( eff7, -1 );     eff7_Ddo->Divide( eff7);

  }

  TH1D* eff1_Aup = (TH1D*)fAup->Get("hptEffPP");
  TH1D* eff2_Aup = (TH1D*)fAup->Get("hptEffAA");
  TH1D* eff3_Aup = (TH1D*)fAup->Get("hrapEffPP");
  TH1D* eff4_Aup = (TH1D*)fAup->Get("hrapEffAA");
  TH1D* eff5_Aup = (TH1D*)fAup->Get("hcentintEffPP");
  TH1D* eff6_Aup = (TH1D*)fAup->Get("hcentintEffAA");
  TH1D* eff7_Aup = (TH1D*)fAup->Get("hcentEffAA");

  TH1D* eff1_Bup = (TH1D*)fBup->Get("hptEffPP");
  TH1D* eff2_Bup = (TH1D*)fBup->Get("hptEffAA");
  TH1D* eff3_Bup = (TH1D*)fBup->Get("hrapEffPP");
  TH1D* eff4_Bup = (TH1D*)fBup->Get("hrapEffAA");
  TH1D* eff5_Bup = (TH1D*)fBup->Get("hcentintEffPP");
  TH1D* eff6_Bup = (TH1D*)fBup->Get("hcentintEffAA");
  TH1D* eff7_Bup = (TH1D*)fBup->Get("hcentEffAA");

  TH1D* eff1_Ado = (TH1D*)fAdo->Get("hptEffPP");
  TH1D* eff2_Ado = (TH1D*)fAdo->Get("hptEffAA");
  TH1D* eff3_Ado = (TH1D*)fAdo->Get("hrapEffPP");
  TH1D* eff4_Ado = (TH1D*)fAdo->Get("hrapEffAA");
  TH1D* eff5_Ado = (TH1D*)fAdo->Get("hcentintEffPP");
  TH1D* eff6_Ado = (TH1D*)fAdo->Get("hcentintEffAA");
  TH1D* eff7_Ado = (TH1D*)fAdo->Get("hcentEffAA");

  TH1D* eff1_Bdo = (TH1D*)fBdo->Get("hptEffPP");
  TH1D* eff2_Bdo = (TH1D*)fBdo->Get("hptEffAA");
  TH1D* eff3_Bdo = (TH1D*)fBdo->Get("hrapEffPP");
  TH1D* eff4_Bdo = (TH1D*)fBdo->Get("hrapEffAA");
  TH1D* eff5_Bdo = (TH1D*)fBdo->Get("hcentintEffPP");
  TH1D* eff6_Bdo = (TH1D*)fBdo->Get("hcentintEffAA");
  TH1D* eff7_Bdo = (TH1D*)fBdo->Get("hcentEffAA");


  eff1_Aup->Add( eff1, -1 );     eff1_Aup->Divide( eff1);
  eff2_Aup->Add( eff2, -1 );     eff2_Aup->Divide( eff2);
  eff3_Aup->Add( eff3, -1 );     eff3_Aup->Divide( eff3);
  eff4_Aup->Add( eff4, -1 );     eff4_Aup->Divide( eff4);
  eff5_Aup->Add( eff5, -1 );     eff5_Aup->Divide( eff5);
  eff6_Aup->Add( eff6, -1 );     eff6_Aup->Divide( eff6);
  eff7_Aup->Add( eff7, -1 );     eff7_Aup->Divide( eff7);

  eff1_Bup->Add( eff1, -1 );     eff1_Bup->Divide( eff1);
  eff2_Bup->Add( eff2, -1 );     eff2_Bup->Divide( eff2);
  eff3_Bup->Add( eff3, -1 );     eff3_Bup->Divide( eff3);
  eff4_Bup->Add( eff4, -1 );     eff4_Bup->Divide( eff4);
  eff5_Bup->Add( eff5, -1 );     eff5_Bup->Divide( eff5);
  eff6_Bup->Add( eff6, -1 );     eff6_Bup->Divide( eff6);
  eff7_Bup->Add( eff7, -1 );     eff7_Bup->Divide( eff7);

  eff1_Ado->Add( eff1, -1 );     eff1_Ado->Divide( eff1);
  eff2_Ado->Add( eff2, -1 );     eff2_Ado->Divide( eff2);
  eff3_Ado->Add( eff3, -1 );     eff3_Ado->Divide( eff3);
  eff4_Ado->Add( eff4, -1 );     eff4_Ado->Divide( eff4);
  eff5_Ado->Add( eff5, -1 );     eff5_Ado->Divide( eff5);
  eff6_Ado->Add( eff6, -1 );     eff6_Ado->Divide( eff6);
  eff7_Ado->Add( eff7, -1 );     eff7_Ado->Divide( eff7);

  eff1_Bdo->Add( eff1, -1 );     eff1_Bdo->Divide( eff1);
  eff2_Bdo->Add( eff2, -1 );     eff2_Bdo->Divide( eff2);
  eff3_Bdo->Add( eff3, -1 );     eff3_Bdo->Divide( eff3);
  eff4_Bdo->Add( eff4, -1 );     eff4_Bdo->Divide( eff4);
  eff5_Bdo->Add( eff5, -1 );     eff5_Bdo->Divide( eff5);
  eff6_Bdo->Add( eff6, -1 );     eff6_Bdo->Divide( eff6);
  eff7_Bdo->Add( eff7, -1 );     eff7_Bdo->Divide( eff7);

  TH1D* eff1ptw = (TH1D*)eff1_Aup->Clone("hptRelSysPtwPP");   eff1ptw->Reset();
  TH1D* eff2ptw = (TH1D*)eff2_Aup->Clone("hptRelSysPtwAA");   eff2ptw->Reset();
  TH1D* eff3ptw = (TH1D*)eff3_Aup->Clone("hrapRelSysPtwPP");   eff3ptw->Reset();
  TH1D* eff4ptw = (TH1D*)eff4_Aup->Clone("hrapRelSysPtwAA");   eff4ptw->Reset();
  TH1D* eff5ptw = (TH1D*)eff5_Aup->Clone("hcentintRelSysPtwPP");   eff5ptw->Reset();
  TH1D* eff6ptw = (TH1D*)eff6_Aup->Clone("hcentintRelSysPtwAA");   eff6ptw->Reset();
  TH1D* eff7ptw = (TH1D*)eff7_Aup->Clone("hcentRelSysPtwAA");   eff7ptw->Reset();
  
  eff1ptw ->SetName("hptRelSysPtwPP");
  eff2ptw ->SetName("hptRelSysPtwAA");
  eff3ptw ->SetName("hrapRelSysPtwPP");
  eff4ptw ->SetName("hrapRelSysPtwAA");
  eff5ptw ->SetName("hcentintRelSysPtwPP");
  eff6ptw ->SetName("hcentintRelSysPtwAA");
  eff7ptw ->SetName("hcentRelSysPtwAA");

  cout << "eff1_Aup : " << eff1_Aup->GetBinContent(1) << endl;
  if(state == 1)
  {
    getMaxTH1D_eight(eff1ptw, eff1_Aup, eff1_Ado, eff1_Bup, eff1_Bdo, eff1_Cup, eff1_Cdo, eff1_Dup, eff1_Ddo);
    getMaxTH1D_eight(eff2ptw, eff2_Aup, eff2_Ado, eff2_Bup, eff2_Bdo, eff2_Cup, eff2_Cdo, eff2_Dup, eff2_Ddo);
    getMaxTH1D_eight(eff3ptw, eff3_Aup, eff3_Ado, eff3_Bup, eff3_Bdo, eff3_Cup, eff3_Cdo, eff3_Dup, eff3_Ddo);
    getMaxTH1D_eight(eff4ptw, eff4_Aup, eff4_Ado, eff4_Bup, eff4_Bdo, eff4_Cup, eff4_Cdo, eff4_Dup, eff4_Ddo);
    getMaxTH1D_eight(eff5ptw, eff5_Aup, eff5_Ado, eff5_Bup, eff5_Bdo, eff5_Cup, eff5_Cdo, eff5_Dup, eff5_Ddo);
    getMaxTH1D_eight(eff6ptw, eff6_Aup, eff6_Ado, eff6_Bup, eff6_Bdo, eff6_Cup, eff6_Cdo, eff6_Dup, eff6_Ddo);
    getMaxTH1D_eight(eff7ptw, eff7_Aup, eff7_Ado, eff7_Bup, eff7_Bdo, eff7_Cup, eff7_Cdo, eff7_Dup, eff7_Ddo);
  }
  else if(state != 1)
  {
    getMaxTH1D_four(eff1ptw, eff1_Bup, eff1_Bdo, eff1_Aup, eff1_Ado);
    getMaxTH1D_four(eff2ptw, eff2_Bup, eff2_Bdo, eff2_Aup, eff2_Ado);
    getMaxTH1D_four(eff3ptw, eff3_Bup, eff3_Bdo, eff3_Aup, eff3_Ado);
    getMaxTH1D_four(eff4ptw, eff4_Bup, eff4_Bdo, eff4_Aup, eff4_Ado);
    getMaxTH1D_four(eff5ptw, eff5_Bup, eff5_Bdo, eff5_Aup, eff5_Ado);
    getMaxTH1D_four(eff6ptw, eff6_Bup, eff6_Bdo, eff6_Aup, eff6_Ado);
    getMaxTH1D_four(eff7ptw, eff7_Bup, eff7_Bdo, eff7_Aup, eff7_Ado);
  }


  cout << "ptw output : " << eff1ptw->GetBinContent(1) << endl;
  
  // stat trig
  for ( int idx=1 ; idx<= Nsamples ; idx++) {
    if ( idx%50 == 0) { cout <<"Reading "  << idx<<"th file..." << endl;} 
    TFile* fid = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId%d_trkId0_muId-100_staId-100.root",state,idx) );
    TH1D* relativeDev1 = (TH1D*)fid->Get("hptEffPP");
    TH1D* relativeDev2 = (TH1D*)fid->Get("hptEffAA");
    TH1D* relativeDev3 = (TH1D*)fid->Get("hrapEffPP");
    TH1D* relativeDev4 = (TH1D*)fid->Get("hrapEffAA");
    TH1D* relativeDev5 = (TH1D*)fid->Get("hcentintEffPP");
    TH1D* relativeDev6 = (TH1D*)fid->Get("hcentintEffAA");
    TH1D* relativeDev7 = (TH1D*)fid->Get("hcentEffAA");
    
    cout << " id : " << idx << endl;
    cout << "relativeDev1 : " << relativeDev1->GetBinContent(1) << endl;
    cout << "efficiency : " << eff1->GetBinContent(1) << endl;
    cout << endl;
    relativeDev1->Add( eff1, -1 );     relativeDev1->Divide( eff1);
    relativeDev2->Add( eff2, -1 );     relativeDev2->Divide( eff2);
    relativeDev3->Add( eff3, -1 );     relativeDev3->Divide( eff3);
    relativeDev4->Add( eff4, -1 );     relativeDev4->Divide( eff4);
    relativeDev5->Add( eff5, -1 );     relativeDev5->Divide( eff5);
    relativeDev6->Add( eff6, -1 );     relativeDev6->Divide( eff6);
    relativeDev7->Add( eff7, -1 );     relativeDev7->Divide( eff7);


    addInQuad( eff1statTrg, relativeDev1 ) ;
    addInQuad( eff2statTrg, relativeDev2 ) ;
    addInQuad( eff3statTrg, relativeDev3 ) ;
    addInQuad( eff4statTrg, relativeDev4 ) ;
    addInQuad( eff5statTrg, relativeDev5 ) ;
    addInQuad( eff6statTrg, relativeDev6 ) ;
    addInQuad( eff7statTrg, relativeDev7 ) ;
    fid->Close();
  }
  eff1statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff2statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff3statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff4statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff5statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff6statTrg->Scale( 1./sqrt( float(Nsamples) ) );
  eff7statTrg->Scale( 1./sqrt( float(Nsamples) ) );

  cout << "stat Tgr :!!! " << eff1statTrg->GetBinContent(1) << endl; 


  // stat muid
  for ( int idx=1 ; idx<= Nsamples ; idx++) {
    if ( idx%50 == 0) { cout <<"Reading "  << idx<<"th file..." << endl;}
    TFile* fid = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId%d_staId-100.root",state,idx) );
    TH1D* relativeDev1 = (TH1D*)fid->Get("hptEffPP");
    TH1D* relativeDev2 = (TH1D*)fid->Get("hptEffAA");
    TH1D* relativeDev3 = (TH1D*)fid->Get("hrapEffPP");
    TH1D* relativeDev4 = (TH1D*)fid->Get("hrapEffAA");
    TH1D* relativeDev5 = (TH1D*)fid->Get("hcentintEffPP");
    TH1D* relativeDev6 = (TH1D*)fid->Get("hcentintEffAA");
    TH1D* relativeDev7 = (TH1D*)fid->Get("hcentEffAA");
    
    relativeDev1->Add( eff1, -1 );     relativeDev1->Divide( eff1);
    relativeDev2->Add( eff2, -1 );     relativeDev2->Divide( eff2);
    relativeDev3->Add( eff3, -1 );     relativeDev3->Divide( eff3);
    relativeDev4->Add( eff4, -1 );     relativeDev4->Divide( eff4);
    relativeDev5->Add( eff5, -1 );     relativeDev5->Divide( eff5);
    relativeDev6->Add( eff6, -1 );     relativeDev6->Divide( eff6);
    relativeDev7->Add( eff7, -1 );     relativeDev7->Divide( eff7);

    addInQuad( eff1statMuId, relativeDev1 ) ;
    addInQuad( eff2statMuId, relativeDev2 ) ;
    addInQuad( eff3statMuId, relativeDev3 ) ;
    addInQuad( eff4statMuId, relativeDev4 ) ;
    addInQuad( eff5statMuId, relativeDev5 ) ;
    addInQuad( eff6statMuId, relativeDev6 ) ;
    addInQuad( eff7statMuId, relativeDev7 ) ;
    fid->Close();
  }
  eff1statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff2statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff3statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff4statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff5statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff6statMuId->Scale( 1./sqrt( float(Nsamples) ) );
  eff7statMuId->Scale( 1./sqrt( float(Nsamples) ) );

  // stat sta
  for ( int idx=1 ; idx<= Nsamples ; idx++) {
    if ( idx%50 == 0) { cout <<"Reading "  << idx<<"th file..." << endl;}
    TFile* fid = new TFile(Form("efficiencyTableSys/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId%d.root",state,idx) );
    TH1D* relativeDev1 = (TH1D*)fid->Get("hptEffPP");
    TH1D* relativeDev2 = (TH1D*)fid->Get("hptEffAA");
    TH1D* relativeDev3 = (TH1D*)fid->Get("hrapEffPP");
    TH1D* relativeDev4 = (TH1D*)fid->Get("hrapEffAA");
    TH1D* relativeDev5 = (TH1D*)fid->Get("hcentintEffPP");
    TH1D* relativeDev6 = (TH1D*)fid->Get("hcentintEffAA");
    TH1D* relativeDev7 = (TH1D*)fid->Get("hcentEffAA");
    
    relativeDev1->Add( eff1, -1 );     relativeDev1->Divide( eff1);
    relativeDev2->Add( eff2, -1 );     relativeDev2->Divide( eff2);
    relativeDev3->Add( eff3, -1 );     relativeDev3->Divide( eff3);
    relativeDev4->Add( eff4, -1 );     relativeDev4->Divide( eff4);
    relativeDev5->Add( eff5, -1 );     relativeDev5->Divide( eff5);
    relativeDev6->Add( eff6, -1 );     relativeDev6->Divide( eff6);
    relativeDev7->Add( eff7, -1 );     relativeDev7->Divide( eff7);

    addInQuad( eff1statSta, relativeDev1 ) ;
    addInQuad( eff2statSta, relativeDev2 ) ;
    addInQuad( eff3statSta, relativeDev3 ) ;
    addInQuad( eff4statSta, relativeDev4 ) ;
    addInQuad( eff5statSta, relativeDev5 ) ;
    addInQuad( eff6statSta, relativeDev6 ) ;
    addInQuad( eff7statSta, relativeDev7 ) ;
    fid->Close();
  }
  eff1statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff2statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff3statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff4statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff5statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff6statSta->Scale( 1./sqrt( float(Nsamples) ) );
  eff7statSta->Scale( 1./sqrt( float(Nsamples) ) );
  

  cout << " eff3TrgVar : " << eff3TrgVar->GetBinContent(1) << endl;
  cout << " eff3MuIdVar : " << eff3MuIdVar->GetBinContent(1) << endl;
  cout << " eff3StaVar : " << eff3StaVar->GetBinContent(1) << endl;
  cout << " eff3TrkVar : " << eff3TrkVar->GetBinContent(1) << endl;
  cout << " eff3binned : " << eff3binned->GetBinContent(1) << endl;
  cout << " eff3ptw : " << eff1ptw->GetBinContent(1) << endl;
  cout << " eff3statTrg : " << eff3statTrg->GetBinContent(1) << endl;
  cout << " eff3statMuId : " << eff3statMuId->GetBinContent(1) << endl;
  cout << " eff3statsta : " << eff3statSta->GetBinContent(1) << endl;
  // Merge them! 
  addInQuadNine( eff1sys,  eff1TrgVar, eff1MuIdVar, eff1StaVar, eff1TrkVar, eff1binned, eff1ptw, eff1statTrg, eff1statMuId, eff1statSta);
  addInQuadNine( eff2sys,  eff2TrgVar, eff2MuIdVar, eff2StaVar, eff2TrkVar, eff2binned, eff2ptw, eff2statTrg, eff2statMuId, eff2statSta);
  addInQuadNine( eff3sys,  eff3TrgVar, eff3MuIdVar, eff3StaVar, eff3TrkVar, eff3binned, eff3ptw, eff3statTrg, eff3statMuId, eff3statSta);
  addInQuadNine( eff4sys,  eff4TrgVar, eff4MuIdVar, eff4StaVar, eff4TrkVar, eff4binned, eff4ptw, eff4statTrg, eff4statMuId, eff4statSta);
  addInQuadNine( eff5sys,  eff5TrgVar, eff5MuIdVar, eff5StaVar, eff5TrkVar, eff5binned, eff5ptw, eff5statTrg, eff5statMuId, eff5statSta);
  addInQuadNine( eff6sys,  eff6TrgVar, eff6MuIdVar, eff6StaVar, eff6TrkVar, eff6binned, eff6ptw, eff6statTrg, eff6statMuId, eff6statSta);
  addInQuadNine( eff7sys,  eff7TrgVar, eff7MuIdVar, eff7StaVar, eff7TrkVar, eff7binned, eff7ptw, eff7statTrg, eff7statMuId, eff7statSta);
  
/*  addInQuadEight( eff1sys,  eff1TrgVar, eff1MuIdVar, eff1StaVar, eff1TrkVar, eff1binned, eff1statTrg, eff1statMuId, eff1statSta);
  addInQuadEight( eff2sys,  eff2TrgVar, eff2MuIdVar, eff2StaVar, eff2TrkVar, eff2binned, eff2statTrg, eff2statMuId, eff2statSta);
  addInQuadEight( eff3sys,  eff3TrgVar, eff3MuIdVar, eff3StaVar, eff3TrkVar, eff3binned, eff3statTrg, eff3statMuId, eff3statSta);
  addInQuadEight( eff4sys,  eff4TrgVar, eff4MuIdVar, eff4StaVar, eff4TrkVar, eff4binned, eff4statTrg, eff4statMuId, eff4statSta);
  addInQuadEight( eff5sys,  eff5TrgVar, eff5MuIdVar, eff5StaVar, eff5TrkVar, eff5binned, eff5statTrg, eff5statMuId, eff5statSta);
  addInQuadEight( eff6sys,  eff6TrgVar, eff6MuIdVar, eff6StaVar, eff6TrkVar, eff6binned, eff6statTrg, eff6statMuId, eff6statSta);
  addInQuadEight( eff7sys,  eff7TrgVar, eff7MuIdVar, eff7StaVar, eff7TrkVar, eff7binned, eff7statTrg, eff7statMuId, eff7statSta);
  */
  TFile* fout = new TFile(Form("sys_efficiency_ups%d.root",state),"recreate");
  eff1sys->Write(); // eff1sysVar->Write(); eff1binned->Write(); eff1sta->Write(); eff1muid->Write(); eff1stat->Write(); eff1ptw->Write();
  eff2sys->Write(); // eff2sysVar->Write(); eff2binned->Write(); eff2sta->Write(); eff2muid->Write(); eff2stat->Write(); eff2ptw->Write();
  eff3sys->Write(); // eff3sysVar->Write(); eff3binned->Write(); eff3sta->Write(); eff3muid->Write(); eff3stat->Write(); eff3ptw->Write();
  eff4sys->Write(); // eff4sysVar->Write(); eff4binned->Write(); eff4sta->Write(); eff4muid->Write(); eff4stat->Write(); eff4ptw->Write();
  eff5sys->Write();  //eff5sysVar->Write(); eff5binned->Write(); eff5sta->Write(); eff5muid->Write(); eff5stat->Write(); eff5ptw->Write();
  eff6sys->Write();  //eff6sysVar->Write(); eff6binned->Write(); eff6sta->Write(); eff6muid->Write(); eff6stat->Write(); eff6ptw->Write();
  eff7sys->Write();  //eff7sysVar->Write(); eff7binned->Write(); eff7sta->Write(); eff7muid->Write(); eff7stat->Write(); eff7ptw->Write();
  
  //Print the results for AN table
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double *yBin;
  if ( state == 1 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
    nCentBins = nCentBins1s;  centBin = centBin1s;
  }
  else if ( state == 2 ) {
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
    nYBins = nYBins2S;  yBin = yBin2S;
  }
  else if ( state == 3 ) {
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
    nYBins = nYBins3S;  yBin = yBin3S;
  }
  
  //  pt dependence
  //  TH1D* eff1 = (TH1D*)f1->Get("hptEffPP");
  //  TH1D* eff2 = (TH1D*)f1->Get("hptEffAA");


  for ( int ii = 1 ; ii<= nPtBins ; ii++)   {
    if ( state == 1 ) {
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << "$ \\GeVc &" <<  int(eff1sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<  int
	(eff2sys->GetBinContent(ii) *1000) / 10. << "\\% & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << "$ \\GeVc & & & " <<  int(eff1sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<	int(eff2sys->GetBinContent(ii) *1000) / 10. << "\\% & &   \\\\ " << endl;
    }
    if ( state == 3 ) {
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << " $ \\GeVc & & & & & " <<  int(eff1sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<  int(eff2sys->GetBinContent(ii) *1000) / 10. << "\\% \\\\ " << endl;
    }
  }

  // rap dependence 
  //  TH1D* eff3 = (TH1D*)f1->Get("hrapEffPP");
  //  TH1D* eff4 = (TH1D*)f1->Get("hrapEffAA");

  for ( int ii = 1 ; ii<= nYBins ; ii++)   {
    if ( state == 1 ) {
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ &" <<  int(eff3sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<  int(eff4sys->GetBinContent(ii) *1000) / 10. << "\\% & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ & & & " <<  int(eff3sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<  int(eff4sys->GetBinContent(ii) *1000) / 10. << "\\% & &   \\\\ " << endl;
    }
    if ( state == 3 ) {
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ & & & & & " <<  int(eff3sys->GetBinContent(ii) *1000) / 10. << "\\% & " <<  int
	(eff4sys->GetBinContent(ii) *1000) / 10. << "\\% \\\\ " << endl;
    }
  }

  // centrality dependence 
  //  TH1D* eff5 = (TH1D*)f1->Get("hcentintEffPP");
  //  TH1D* eff6 = (TH1D*)f1->Get("hcentintEffAA");
  //  TH1D* eff7 = (TH1D*)f1->Get("hcentEffAA");

  // Print the results for the table in for AN
  for ( int ii = 1 ; ii<= nCentBins ; ii++)   {
    if ( state == 1 ) {
      cout << "$" << centBin[ii-1]/2 << "\\% - " << centBin[ii]/2 << "\\% $ &   & " <<  int(eff7sys->GetBinContent(ii) *1000) / 10. << "\\% & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {
      cout << "$" << centBin[ii-1]/2 << "\\% - " << centBin[ii]/2 << "\\% $ & & &   & " <<  int(eff7sys->GetBinContent(ii) *1000) / 10. << "\\% &&   \\\\ " << endl;
    }
    if ( state == 3 ) {
      cout << "$" << centBin[ii-1]/2 << "\\% - " << centBin[ii]/2 << "\\% $ & & & & &   & " <<  int(eff7sys->GetBinContent(ii) *1000) / 10. << 
	"\\% \\\\ " << endl;
    }
  }

  // Integrated bin:
  cout << " pp = " << int(eff5sys->GetBinContent(1)*1000)/10. << "\\%,  PbPb = " << int(eff6sys->GetBinContent(1)*1000)/10. << "\\%" <<endl;




}

void addInQuad ( TH1D* h0, TH1D* h1) { 
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x0 =  h0->GetBinContent(i);
    float x1 =  h1->GetBinContent(i);
    h0->SetBinContent(i,  sqrt ( x0*x0 + x1*x1 ) );
  }
}
void addInQuadFive ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5) { 
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  h1->GetBinContent(i);
    float x2 =  h2->GetBinContent(i);
    float x3 =  h3->GetBinContent(i);
    float x4 =  h4->GetBinContent(i);
    float x5 =  h5->GetBinContent(i);
    h0->SetBinContent(i,  sqrt ( x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 ) );
  }
}
void addInQuadSix ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6) { 
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  h1->GetBinContent(i);
    float x2 =  h2->GetBinContent(i);
    float x3 =  h3->GetBinContent(i);
    float x4 =  h4->GetBinContent(i);
    float x5 =  h5->GetBinContent(i);
    float x6 =  h6->GetBinContent(i);
    h0->SetBinContent(i,  sqrt ( x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6) );
  }
}

void addInQuadNine ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, TH1D* h9) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h7->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h8->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h9->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  h1->GetBinContent(i);
    float x2 =  h2->GetBinContent(i);
    float x3 =  h3->GetBinContent(i);
    float x4 =  h4->GetBinContent(i);
    float x5 =  h5->GetBinContent(i);
    float x6 =  h6->GetBinContent(i);
    float x7 =  h7->GetBinContent(i);
    float x8 =  h8->GetBinContent(i);
    float x9 =  h9->GetBinContent(i);
    h0->SetBinContent(i,  sqrt ( x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8 + x9*x9) );
  }
}

void addInQuadEight ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h7->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h8->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  h1->GetBinContent(i);
    float x2 =  h2->GetBinContent(i);
    float x3 =  h3->GetBinContent(i);
    float x4 =  h4->GetBinContent(i);
    float x5 =  h5->GetBinContent(i);
    float x6 =  h6->GetBinContent(i);
    float x7 =  h7->GetBinContent(i);
    float x8 =  h8->GetBinContent(i);
    h0->SetBinContent(i,  sqrt ( x1*x1 + x2*x2 + x3*x3 + x4*x4 + x5*x5 + x6*x6 + x7*x7 + x8*x8) );
  }
}

void getMaxTH1D_eight ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h7->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h8->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    float x3 =  fabs(h3->GetBinContent(i));
    float x4 =  fabs(h4->GetBinContent(i));
    float x5 =  fabs(h5->GetBinContent(i));
    float x6 =  fabs(h6->GetBinContent(i));
    float x7 =  fabs(h7->GetBinContent(i));
    float x8 =  fabs(h8->GetBinContent(i));
    float x = 0;
    if ( x1 > x2 ) x = x1;
    else x = x2;
    if ( x <= x3) x = x3;
    if ( x <= x4) x = x4;
    if ( x <= x5) x = x5;
    if ( x <= x6) x = x6;
    if ( x <= x7) x = x7;
    if ( x <= x8) x = x8;
    h0->SetBinContent(i, x );
  }
}


void getMaxTH1D_six ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h5->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h6->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    float x3 =  fabs(h3->GetBinContent(i));
    float x4 =  fabs(h4->GetBinContent(i));
    float x5 =  fabs(h5->GetBinContent(i));
    float x6 =  fabs(h6->GetBinContent(i));
    float x = 0;
    if ( x1 > x2 ) x = x1;
    else x = x2;
    if ( x <= x3) x = x3;
    if ( x <= x4) x = x4;
    if ( x <= x5) x = x5;
    if ( x <= x6) x = x6;
    h0->SetBinContent(i, x );
  }
}

void getMaxTH1D_four ( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h3->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h4->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    float x3 =  fabs(h3->GetBinContent(i));
    float x4 =  fabs(h4->GetBinContent(i));
    float x = 0;
    if ( x1 > x2 ) x = x1;
    else x = x2;
    if ( x <= x3) x = x3;
    if ( x <= x4) x = x4;
    h0->SetBinContent(i, x );
  }
}

void getMaxTH1D ( TH1D* h0, TH1D* h1, TH1D* h2) {
  if ( h0->GetNbinsX() != h1->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;
  if ( h0->GetNbinsX() != h2->GetNbinsX() )   cout << " Bin numbers are not consistent!" << endl;

  for ( int i=1 ;  i<=h0->GetNbinsX(); i++) {
    float x1 =  fabs(h1->GetBinContent(i));
    float x2 =  fabs(h2->GetBinContent(i));
    if ( x1 > x2 ) 
      h0->SetBinContent(i, x1 );
    else
      h0->SetBinContent(i, x2 );
  }
}

