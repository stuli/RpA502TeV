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
void addInQuadFive ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0);
void addInQuadSix ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0);
void addInQuadNine ( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D* h5=0, TH1D* h6=0, TH1D* h7=0, TH1D* h8=0, TH1D* h9=0);
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
  getMaxTH1D(eff1TrgVar, relTrgVar1_2, relTrgVar1_2);
  getMaxTH1D(eff2TrgVar, relTrgVar2_2, relTrgVar2_2);
  getMaxTH1D(eff3TrgVar, relTrgVar3_2, relTrgVar3_2);
  getMaxTH1D(eff4TrgVar, relTrgVar4_2, relTrgVar4_2);
  getMaxTH1D(eff5TrgVar, relTrgVar5_2, relTrgVar5_2);
  getMaxTH1D(eff6TrgVar, relTrgVar6_2, relTrgVar6_2);
  getMaxTH1D(eff7TrgVar, relTrgVar7_2, relTrgVar7_2);



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
  getMaxTH1D(eff1MuIdVar, relMuIdVar1_2, relMuIdVar1_2);
  getMaxTH1D(eff2MuIdVar, relMuIdVar2_2, relMuIdVar2_2);
  getMaxTH1D(eff3MuIdVar, relMuIdVar3_2, relMuIdVar3_2);
  getMaxTH1D(eff4MuIdVar, relMuIdVar4_2, relMuIdVar4_2);
  getMaxTH1D(eff5MuIdVar, relMuIdVar5_2, relMuIdVar5_2);
  getMaxTH1D(eff6MuIdVar, relMuIdVar6_2, relMuIdVar6_2);
  getMaxTH1D(eff7MuIdVar, relMuIdVar7_2, relMuIdVar7_2);

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
  getMaxTH1D(eff1StaVar, relStaVar1_2, relStaVar1_2);
  getMaxTH1D(eff2StaVar, relStaVar2_2, relStaVar2_2);
  getMaxTH1D(eff3StaVar, relStaVar3_2, relStaVar3_2);
  getMaxTH1D(eff4StaVar, relStaVar4_2, relStaVar4_2);
  getMaxTH1D(eff5StaVar, relStaVar5_2, relStaVar5_2);
  getMaxTH1D(eff6StaVar, relStaVar6_2, relStaVar6_2);
  getMaxTH1D(eff7StaVar, relStaVar7_2, relStaVar7_2);


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
  getMaxTH1D(eff1TrkVar, relTrkVar1_2, relTrkVar1_2);
  getMaxTH1D(eff2TrkVar, relTrkVar2_2, relTrkVar2_2);
  getMaxTH1D(eff3TrkVar, relTrkVar3_2, relTrkVar3_2);
  getMaxTH1D(eff4TrkVar, relTrkVar4_2, relTrkVar4_2);
  getMaxTH1D(eff5TrkVar, relTrkVar5_2, relTrkVar5_2);
  getMaxTH1D(eff6TrkVar, relTrkVar6_2, relTrkVar6_2);
  getMaxTH1D(eff7TrkVar, relTrkVar7_2, relTrkVar7_2);

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
  TFile* fid_ptw = new TFile(Form("efficiencyTable/efficiency_ups%ds_useDataPtWeight0_tnp_trgId0_trkId0_muId-100_staId-100.root",state )) ;
  TH1D* eff1ptw = (TH1D*)fid_ptw->Get("hptEffPP");
  TH1D* eff2ptw = (TH1D*)fid_ptw->Get("hptEffAA");
  TH1D* eff3ptw = (TH1D*)fid_ptw->Get("hrapEffPP");
  TH1D* eff4ptw = (TH1D*)fid_ptw->Get("hrapEffAA");
  TH1D* eff5ptw = (TH1D*)fid_ptw->Get("hcentintEffPP");
  TH1D* eff6ptw = (TH1D*)fid_ptw->Get("hcentintEffAA");
  TH1D* eff7ptw = (TH1D*)fid_ptw->Get("hcentEffAA");
  eff1ptw ->SetName("hptRelSysPtwPP");
  eff2ptw ->SetName("hptRelSysPtwAA");
  eff3ptw ->SetName("hrapRelSysPtwPP");
  eff4ptw ->SetName("hrapRelSysPtwAA");
  eff5ptw ->SetName("hcentintRelSysPtwPP");
  eff6ptw ->SetName("hcentintRelSysPtwAA");
  eff7ptw ->SetName("hcentRelSysPtwAA");
  eff1ptw->Add( eff1, -1 );     eff1ptw->Divide( eff1);
  eff2ptw->Add( eff2, -1 );     eff2ptw->Divide( eff2);
  eff3ptw->Add( eff3, -1 );     eff3ptw->Divide( eff3);
  eff4ptw->Add( eff4, -1 );     eff4ptw->Divide( eff4);
  eff5ptw->Add( eff5, -1 );     eff5ptw->Divide( eff5);
  eff6ptw->Add( eff6, -1 );     eff6ptw->Divide( eff6);
  eff7ptw->Add( eff7, -1 );     eff7ptw->Divide( eff7);

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
  

  // Merge them! 
  addInQuadNine( eff1sys,  eff1TrgVar, eff1MuIdVar, eff1StaVar, eff1TrkVar, eff1binned, eff1ptw, eff1statTrg, eff1statMuId, eff1statSta);
  addInQuadNine( eff2sys,  eff2TrgVar, eff2MuIdVar, eff2StaVar, eff2TrkVar, eff2binned, eff2ptw, eff2statTrg, eff2statMuId, eff2statSta);
  addInQuadNine( eff3sys,  eff3TrgVar, eff3MuIdVar, eff3StaVar, eff3TrkVar, eff3binned, eff3ptw, eff3statTrg, eff3statMuId, eff3statSta);
  addInQuadNine( eff4sys,  eff4TrgVar, eff4MuIdVar, eff4StaVar, eff4TrkVar, eff4binned, eff4ptw, eff4statTrg, eff4statMuId, eff4statSta);
  addInQuadNine( eff5sys,  eff5TrgVar, eff5MuIdVar, eff5StaVar, eff5TrkVar, eff5binned, eff5ptw, eff5statTrg, eff5statMuId, eff5statSta);
  addInQuadNine( eff6sys,  eff6TrgVar, eff6MuIdVar, eff6StaVar, eff6TrkVar, eff6binned, eff6ptw, eff6statTrg, eff6statMuId, eff6statSta);
  addInQuadNine( eff7sys,  eff7TrgVar, eff7MuIdVar, eff7StaVar, eff7TrkVar, eff7binned, eff7ptw, eff7statTrg, eff7statMuId, eff7statSta);
  
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
