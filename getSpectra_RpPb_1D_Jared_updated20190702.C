#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "cutsAndBin_Santona.h"
#include "multiTreeUtil.h"
using namespace std;

TString ResultDir  = "nominalFits";


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0) ;

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 

void getSpectra_RpPb_1D_Jared_updated20190702(int state = 1) {  

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  int nPtBins=0;
  int nYBins=0;
  int nYBins_cr=0;
  double* ptBin;
  double* yBin;
  double* yBin_cr;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;
  //  double nColl1s[nCentBins] = {1819,1432,1005,606,349,186,90.7,40.1,7.67};
  //  double nPart1s[nCentBins] = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};



  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; nYBins_cr = nYBins1S_cr;   yBin_cr = yBin1S_cr;
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S;  nYBins_cr = nYBins2S_cr;   yBin_cr = yBin2S_cr;
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S;   nYBins_cr = nYBins3S_cr;   yBin_cr = yBin3S_cr;
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  

  // signals :
  TH1D* hRPAraw_rap;   // w/o efficiency correction
  TH1D* hRPAraw_int;   // w/o efficiency correction
  TH1D* hRPA_rap;   // w efficiency correction
  TH1D* hRPA_int;   // w efficiency correction
  TH1D* hRPAraw_pt;   // w/o efficiency correction
  TH1D* hRPA_pt;   // w efficiency correction

  TH1D* hrapAccPA;
  TH1D* hrapAccPA_cross;
  TH1D* hrapAccPP;
  TH1D* hrapAccPP_cross;
  TH1D* hptAccPA;
  TH1D* hptAccPAdw;
  TH1D* hptAccPP;
  TH1D* hintAccPA;
  TH1D* hintAccPP;

  TH1D* hrapEffPA;
  TH1D* hrapEffPA_cross;
  TH1D* hrapEffPP;
  TH1D* hrapEffPP_cross;
  TH1D* hptEffPA;
  TH1D* hptEffPP;
  TH1D* hintEffPA;
  TH1D* hintEffPP;

  TH1D* hrapEff;
  TH1D* hptEff;
  TH1D* hintEff;
  TH1D* hEffPA_cross_rap;
  TH1D* hEffPA_cross_pt;
  TH1D* hEffPP_cross_rap;
  TH1D* hEffPP_cross_pt;
  TH1D* hEffPA_cross_int;
  TH1D* hEffPP_cross_int;

  TH1D* hrapEffPA_gen;
  TH1D* hrapEffPA_gen_cross;
  TH1D* hrapEffPP_gen;
  TH1D* hrapEffPP_gen_cross;
  TH1D* hptEffPA_gen;
  TH1D* hptEffPP_gen;
  TH1D* hintEffPA_gen;
  TH1D* hintEffPP_gen;

  //TFile* infaccdw = new TFile(Form("Acceptance/20180328/acceptance_wgt_%dS_20180328_2Dplot.root",state),"read");
  //TFile* infaccdw = new TFile(Form("Acceptance/20180724/acceptance_wgt_%dS_20180724_2Dplot.root",state),"read");
  //TFile* infaccdw = new TFile(Form("Acceptance/20181217/acceptance_wgt_%dS_20181217_2Dplot.root",state),"read");
  TFile* infaccdw = new TFile(Form("Corrections/Acceptance/20190221/acceptance_wgt_%dS_20190221_2Dplot.root",state),"read");
  TH1D* hrapAccPA_  = (TH1D*)infaccdw->Get(Form("hrapAccPA_%dS",state));
  TH1D* hrapAccPP_  = (TH1D*)infaccdw->Get(Form("hrapAccPP_%dS",state));
  hrapAccPA = (TH1D*)infaccdw->Get(Form("hrapAccPA_%dS",state));
  hrapAccPP = (TH1D*)infaccdw->Get(Form("hrapAccPP_%dS",state));
/*  for(int i=1;i<=hrapAccPA->GetNbinsX();i++){
    if(i>1) {hrapAccPA->SetBinContent(i,hrapAccPA_->GetBinContent(i-1)); hrapAccPP->SetBinContent(i,hrapAccPP_->GetBinContent(i-1));}
  }
*/
  hptAccPA  = (TH1D*) infaccdw->Get(Form("hptAccPA_%dS",state));
  hptAccPAdw  = (TH1D*) infaccdw->Get(Form("hptAccCross_%dS",state));
  hptAccPP  = (TH1D*) infaccdw->Get(Form("hptAccPP_%dS",state));
  hrapAccPA_cross = (TH1D*) infaccdw->Get(Form("hrapAccCross_%dS",state));
  hrapAccPP_cross = (TH1D*) infaccdw->Get(Form("hrapAccPP_%dS",state));
  hintAccPA  = (TH1D*) infaccdw->Get(Form("hIntAccPA_%dS",state));
  hintAccPP  = (TH1D*) infaccdw->Get(Form("hIntAccPP_%dS",state));

  cout << "hrapAccPA : " << hrapAccPA->GetBinContent(2) << endl;
  cout << "hrapAccPP : " << hrapAccPP->GetBinContent(2) << endl;

  //TFile* infeff_pPb = new TFile(Form("Efficiency_rootfiles/pPb/Eff_pPb_%dS_8_22_NewPtReweights.root",state),"read");
  //TFile* infeff_pPb = new TFile(Form("Efficiency_rootfiles/pPb/Eff_pPb_%dS_11_20_NewRpABin.root",state),"read");
  TFile* infeff = new TFile(Form("Corrections/Efficiency/RootFiles/EffNomCor_SysRpA_%dS.root",state),"read");
  TFile* infeff_crossdw = new TFile(Form("Corrections/Efficiency/RootFiles/EffCor_SyspPbXSAsymm_%dS.root",state),"read");
  TFile* infeff_cross = new TFile(Form("Corrections/Efficiency/RootFiles/EffCor_SyspPbXSSymm_%dS.root",state),"read");
  TFile* infeffPP_cross = new TFile(Form("Corrections/Efficiency/RootFiles/EffCor_SysPPXS_%dS.root",state),"read");
/*  hrapEff  = (TH1D*)infeff->Get("RecoEventsRap");
  hrapEffPA_gen  = (TH1D*)infeff_pPb->Get("GenEventsRap");
  hptEffPA  = (TH1D*) infeff_pPb->Get("RecoEventsPtRpA");
  hptEffPA_gen  = (TH1D*) infeff_pPb->Get("GenEventsPtRpA");
  hintEffPA  = (TH1D*) infeff_pPb->Get("RecoEventsIntRpA");
  hintEffPA_gen  = (TH1D*) infeff_pPb->Get("GenEventsIntRpA");
  hrapEffPP  = (TH1D*)infeff_pp->Get("RecoEventsRap");
  hrapEffPP_gen  = (TH1D*)infeff_pp->Get("GenEventsRap");
  hptEffPP  = (TH1D*) infeff_pp->Get("RecoEventsPt");
  hptEffPP_gen  = (TH1D*) infeff_pp->Get("GenEventsPt");
  hintEffPP  = (TH1D*) infeff_pp->Get("RecoEventsInt");
  hintEffPP_gen  = (TH1D*) infeff_pp->Get("GenEventsInt");

  hrapEffPA->Divide(hrapEffPA_gen);
  hrapEffPP->Divide(hrapEffPP_gen);
  hptEffPA->Divide(hptEffPA_gen);
  hptEffPP->Divide(hptEffPP_gen);
  hintEffPA->Divide(hintEffPA_gen);
  hintEffPP->Divide(hintEffPP_gen);
*/
  hrapEff = (TH1D*)infeff->Get("EffNomRatRap");
  hptEff = (TH1D*)infeff->Get("EffNomRatPt");
  hintEff = (TH1D*)infeff->Get("EffNomRatInt");
  hEffPA_cross_rap = (TH1D*)infeff_crossdw->Get("EffNomRap");
  hEffPA_cross_pt = (TH1D*)infeff_cross->Get("EffNomPt");
  TH1D* hEffPA_cross_ptdw = (TH1D*)infeff_crossdw->Get("EffNomPt");
  hEffPA_cross_int = (TH1D*)infeff_cross->Get("EffNomInt");
  hEffPP_cross_rap = (TH1D*)infeffPP_cross->Get("EffNomRap");
  hEffPP_cross_pt = (TH1D*)infeffPP_cross->Get("EffNomPt");
  hEffPP_cross_int = (TH1D*)infeffPP_cross->Get("EffNomInt");

  stripErrorBars(hrapAccPA_cross);
  stripErrorBars(hrapAccPA);
  stripErrorBars(hrapAccPP_cross);
  stripErrorBars(hrapAccPP);

  stripErrorBars(hptAccPA);
  stripErrorBars(hptAccPAdw);
  stripErrorBars(hptAccPP);
  stripErrorBars(hrapEff);
  stripErrorBars(hintAccPA);
  stripErrorBars(hintAccPP);
  stripErrorBars(hptEff);
  stripErrorBars(hintEff);
  stripErrorBars(hEffPA_cross_rap);
  stripErrorBars(hEffPA_cross_pt);
  stripErrorBars(hEffPA_cross_ptdw);
  stripErrorBars(hEffPA_cross_int);
  stripErrorBars(hEffPP_cross_rap);
  stripErrorBars(hEffPP_cross_pt);
  stripErrorBars(hEffPP_cross_int);
  
  TH1D* hrapSigPP = (TH1D*) hrapAccPP -> Clone("hrapPP");
  TH1D* hrapSigPA = (TH1D*) hrapAccPA -> Clone("hrapPA");
  TH1D* hrapSigPA_cross = (TH1D*) hrapAccPA_cross->Clone("hrapPA_Cross");
  TH1D* hrapSigPP_cross = (TH1D*) hrapAccPP_cross->Clone("hrapPP_Cross");
  TH1D* hptSigPP = (TH1D*) hptAccPP -> Clone("hptPP");
  TH1D* hptSigPA = (TH1D*) hptAccPA -> Clone("hptPA");
  TH1D* hptSigPA_dw = (TH1D*) hptAccPA -> Clone("hptPAdw");
  TH1D* hintSigPP = (TH1D*) hintAccPP -> Clone("hintPP");
  TH1D* hintSigPA = (TH1D*) hintAccPA -> Clone("hintPA");
  hptSigPP->Reset();
  hptSigPA->Reset();
  hptSigPA_dw->Reset();
  hrapSigPA_cross->Reset();
  hrapSigPP_cross->Reset();
  hrapSigPP->Reset();
  hrapSigPA->Reset();
  hintSigPP->Reset();
  hintSigPA->Reset();


  //***yCM***
  //valErr yieldPP;
  TCanvas* c_rap =  new TCanvas("c_rap","",400,400);
  for ( int irap = 1 ; irap<= nYBins ; irap++) {
    valErr yieldPP;
    if(irap <= (nYBins/2)) {
      yieldPP = getYield(state, kPPDATA, 0,30, TMath::Abs(yBin[irap]), TMath::Abs(yBin[irap-1]), 0,200,0,100);
    }
    else if(irap > (nYBins/2)) {
      cout << "irap : " << irap << endl;
      yieldPP = getYield(state, kPPDATA, 0,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    }
    valErr yieldPA = getYield(state, kPADATA, 0,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    hrapSigPP->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP->SetBinError( irap, yieldPP.err/2 ) ;
    hrapSigPP_cross->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP_cross->SetBinError( irap, yieldPP.err/2 ) ;
    hrapSigPA->SetBinContent( irap, yieldPA.val ) ;
    hrapSigPA->SetBinError( irap, yieldPA.err ) ;
  }

  for(int irap = 1; irap<=nYBins_cr; irap++){
    cout << "nYBins_cr : " << nYBins_cr<< endl;
    cout << "yBin_cr[irap-1] to yBin_cr[irap] : " << yBin_cr[irap-1] << " - " << yBin_cr[irap] << endl;
    valErr yieldPA = getYield(state, kPADATA, 0,30, yBin_cr[irap-1],yBin_cr[irap],0,200,0,100);
    cout << "yieldPA.val : " << yieldPA.val << endl;
    cout << "yieldPA.err : " << yieldPA.err << endl;
    hrapSigPA_cross->SetBinContent(irap, yieldPA.val);
    hrapSigPA_cross->SetBinError(irap, yieldPA.err);
    //hrapSigPP_cross->SetBinContent(irap, yieldPP.val/2);
    //hrapSigPP_cross->SetBinError(irap, yieldPP.err/2);
  }


  //yCM Yield 
  hrapSigPP->SetAxisRange(10,1e5,"Y");
  hrapSigPA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hrapSigPP,2);
  handsomeTH1(hrapSigPA,2);
  hrapSigPP->SetMarkerStyle(24);
 
  //***pt***
  TCanvas* c_pt =  new TCanvas("c_pt","",400,400);
  for ( int ipt = 1 ; ipt<= nPtBins ; ipt++) {
    valErr yieldPP = getYield(state, kPPDATA, ptBin[ipt-1],ptBin[ipt], 0.00,1.93, 0,200,0,100);
    valErr yieldPA = getYield(state, kPADATA, ptBin[ipt-1],ptBin[ipt], -1.93,1.93, 0,200,0,100);
    valErr yieldPAdw = getYield(state, kPADATA, ptBin[ipt-1],ptBin[ipt], -2.87,1.93, 0,200,0,100);
    hptSigPA->SetBinContent( ipt, yieldPA.val ) ;
    hptSigPA->SetBinError( ipt, yieldPA.err ) ;
    hptSigPA_dw->SetBinContent( ipt, yieldPAdw.val ) ;
    hptSigPA_dw->SetBinError( ipt, yieldPAdw.err ) ;
    hptSigPP->SetBinContent( ipt, yieldPP.val ) ;
    hptSigPP->SetBinError( ipt, yieldPP.err ) ;
  }

  valErr yieldPP = getYield(state, kPPDATA, 0,30, 0.00,1.93, 0,200,0,100);
  valErr yieldPA = getYield(state, kPADATA, 0,30, -1.93,1.93, 0,200,0,100);
  hintSigPA->SetBinContent( 1, yieldPA.val ) ;
  hintSigPA->SetBinError( 1, yieldPA.err ) ;
  hintSigPP->SetBinContent( 1, yieldPP.val ) ;
  hintSigPP->SetBinError( 1, yieldPP.err ) ;

  //pt Yield
  hptSigPP->SetAxisRange(10,1e5,"Y");
  hptSigPA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hptSigPP,2);
  handsomeTH1(hptSigPA,2);
  hptSigPP->SetMarkerStyle(24);
  double Factor_Scale = lumi_pp*1000/(lumi_pa*208);

 //*****Draw rap***** 
  c_rap->cd(); 
  hrapSigPP->Draw();
//  gPad->SetLogy();
  hrapSigPA->Draw("same");

  hRPAraw_rap = (TH1D*)hrapSigPA->Clone("rpa_vs_rap_raw");
  hRPAraw_rap->Divide(hrapSigPP);

  TH1D* hrel_Acc_rap = (TH1D*) hrapAccPA -> Clone("hrel_Acc_rap");
  TH1D* hrel_Eff_rap = (TH1D*) hrapEff -> Clone("hrel_Eff_rap");

  cout << "hrapSigPA : " << hrapSigPA->GetBinContent(1) << endl;
  cout << "hrapSigPP : " << hrapSigPP->GetBinContent(1) << endl;
  cout << "rap first bin : " << hRPAraw_rap->GetBinContent(1) << endl;

  //Cros sec
  TH1D* hrap_cross_pA = (TH1D*)hrapSigPA_cross->Clone("rpa_vs_rap_cross");
  hrap_cross_pA->Divide(hrapAccPA_cross);
  cout << " :: " << hrapAccPA_cross->GetNbinsX() << endl;
  hrap_cross_pA->Divide(hEffPA_cross_rap);
  hrap_cross_pA->Scale(1./(1000.*lumi_pa));
  TH1ScaleByWidth(hrap_cross_pA);
  //  
  TH1D* hrap_cross_pp = (TH1D*)hrapSigPP_cross->Clone("rpa_vs_rap_cross");
  hrap_cross_pp->Divide(hrapAccPP_cross);
  cout << " :: " << hrapAccPP_cross->GetNbinsX() << endl;
  hrap_cross_pp->Divide(hEffPP_cross_rap);
  hrap_cross_pp->Scale(1./(1000.*lumi_pp));
  TH1ScaleByWidth(hrap_cross_pp);
  //

  hrel_Acc_rap->Divide(hrapAccPP);
  hRPAraw_rap->Divide(hrel_Acc_rap);
  cout << "rap first bin : " << hRPAraw_rap->GetBinContent(1) << endl;
  hRPAraw_rap->Multiply(hrel_Eff_rap);
  cout << "rap first bin : " << hRPAraw_rap->GetBinContent(1) << endl;
  hRPAraw_rap->Scale(Factor_Scale);
  cout << "rap first bin : " << hRPAraw_rap->GetBinContent(1) << endl;


 //*****Draw pt***** 
  c_pt->cd(); 
  //hptSigPP->Draw();
  //gPad->SetLogy();
  //hptSigPA->Draw("same");

  hRPAraw_pt = (TH1D*)hptSigPA->Clone("raa_vs_pt_raw");
  hRPAraw_pt->Divide(hptSigPP);

  TH1D* hrel_Acc_pt = (TH1D*) hptAccPA -> Clone("hrel_Acc_pt");
  TH1D* hrel_Acc_ptdw = (TH1D*) hptAccPAdw -> Clone("hrel_Acc_ptdw");
  TH1D* hrel_Acc_int = (TH1D*) hintAccPA -> Clone("hrel_Acc_int");
  TH1D* hrel_Eff_pt = (TH1D*) hptEff -> Clone("hrel_Eff_pt");
  TH1D* hrel_Eff_int = (TH1D*) hintEff -> Clone("hrel_Eff_int");

  //Cros sec
  TH1D* hpt_cross_pA = (TH1D*)hptSigPA->Clone("rpa_vs_pt_cross");
  TH1D* hpt_cross_pAdw = (TH1D*)hptSigPA_dw->Clone("rpa_vs_pt_crossdw");
  hpt_cross_pA->Divide(hrel_Acc_pt);
  hpt_cross_pA->Divide(hEffPA_cross_pt);
  hpt_cross_pA->Scale(1./(1000.*lumi_pa));
  TH1ScaleByWidth(hpt_cross_pA);
  
  hpt_cross_pAdw->Divide(hrel_Acc_ptdw);
  hpt_cross_pAdw->Divide(hEffPA_cross_ptdw);
  hpt_cross_pAdw->Scale(1./(1000.*lumi_pa));
  TH1ScaleByWidth(hpt_cross_pAdw);
  //
  TH1D* hint_cross_pA = (TH1D*)hintSigPA->Clone("rpa_vs_int_cross");
  hint_cross_pA->Divide(hintAccPA);
  hint_cross_pA->Divide(hEffPA_cross_int);
  hint_cross_pA->Scale(1./(1000.*lumi_pa));
  TH1ScaleByWidth(hint_cross_pA);
  //
  TH1D* hpt_cross_pp = (TH1D*)hptSigPP->Clone("rpa_vs_pt_cross");
  hpt_cross_pp->Divide(hrel_Acc_pt);
  hpt_cross_pp->Divide(hEffPP_cross_pt);
  hpt_cross_pp->Scale(1./(1000.*lumi_pp));
  TH1ScaleByWidth(hpt_cross_pp);
  //
  TH1D* hint_cross_pp = (TH1D*)hintSigPP->Clone("rpa_vs_int_cross");
  hint_cross_pp->Divide(hintAccPP);
  hint_cross_pp->Divide(hEffPP_cross_int);
  hint_cross_pp->Scale(1./(1000.*lumi_pp));
  TH1ScaleByWidth(hint_cross_pp);
  //
  hrel_Acc_pt->Divide(hptAccPP);

  hRPAraw_pt->Divide(hrel_Acc_pt);
  hRPAraw_pt->Multiply(hrel_Eff_pt);
  hRPAraw_pt->Scale(Factor_Scale);
  hrel_Eff_pt->Draw();

  cout << "Acc rap first bin : " << hrapAccPP_cross->GetBinContent(1) << endl;
  cout << "Acc rap int bin : " << hintAccPP->GetBinContent(1) << endl; 
  cout << "Eff rap first bin : " << hEffPP_cross_rap->GetBinContent(1) << endl;
  cout << "Eff rap int bin : " << hEffPP_cross_int->GetBinContent(1) << endl;  
  cout << "XS rap first bin : " << hrap_cross_pp->GetBinContent(1) << endl;
  cout << "XS rap int bin : " << hint_cross_pp->GetBinContent(1) << endl;

  //***** Save ******

  //***** rap *****
  TCanvas* cRapRPA =  new TCanvas("cRpA_rap","",400,400);
  hRPA_rap = (TH1D*)hRPAraw_rap->Clone("rpa_vs_rap");

  hRPA_rap->SetMarkerStyle(24);
  hRPA_rap->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_rap->GetYaxis()->SetTitle("R_{pPb} corrected");
  cRapRPA->cd();
  hRPA_rap->Draw("pe");
  jumSun(yMin,1,yMax,1,1,1);

  TGraphErrors *gRPA_rap = new TGraphErrors(hRPAraw_rap);
  gRPA_rap->SetName("gRPA_rap");
  gRPA_rap->SetTitle("rpa_vs_rap");

  TGraphErrors *gCross_rap = new TGraphErrors(hrap_cross_pA);
  gCross_rap->SetName("gCross_rap");

  TGraphErrors *gCrossPP_rap = new TGraphErrors(hrap_cross_pp);
  gCrossPP_rap->SetName("gCrossPP_rap");


  //******pt*******
  TCanvas* cPtRPA =  new TCanvas("cRpA_pt","",400,400);
  hRPA_pt = (TH1D*)hRPAraw_pt->Clone("rpa_vs_pt");

  hRPA_pt->SetMarkerStyle(24);
  hRPA_pt->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_pt->GetYaxis()->SetTitle("R_{pPb} uncorrected");
  cPtRPA->cd();
  hRPA_pt->Draw("pe");
  jumSun(ptMin,1,ptMax,1,1,1);

  TGraphErrors *gRPA_pt = new TGraphErrors(hRPAraw_pt);
  gRPA_pt->SetName("gRPA_pt");
  gRPA_pt->SetTitle("rpa_vs_pt");
  
  TGraphErrors *gCross_pt = new TGraphErrors(hpt_cross_pA);
  gCross_pt->SetName("gCross_pt");
  TGraphErrors *gCross_ptdw = new TGraphErrors(hpt_cross_pAdw);
  gCross_ptdw->SetName("gCross_ptdw");
  TGraphErrors *gCross_int = new TGraphErrors(hint_cross_pA);
  gCross_int->SetName("gCross_int");
  TGraphErrors *gCrossPP_pt = new TGraphErrors(hpt_cross_pp);
  gCrossPP_pt->SetName("gCrossPP_pt");
  TGraphErrors *gCrossPP_int = new TGraphErrors(hint_cross_pp);
  gCrossPP_int->SetName("gCrossPP_int");

  TFile *wf = new TFile(Form("finalResults/Ups_%d_1D.root",state),"recreate");
  gRPA_rap->Write();
  gRPA_pt->Write();
  gCross_rap->Write();
  gCross_pt->Write();
  gCross_ptdw->Write();
  gCross_int->Write();
  gCrossPP_rap->Write();
  gCrossPP_pt->Write();
  gCrossPP_int->Write();
  wf->Close();
}

  valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh,   float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TFile* inf = new TFile(Form("/home/jared/Documents/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("NominalFitResult/jaredFit/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/afs/cern.ch/work/j/jaebeom/private/Analysis/RpA502TeV/NominalFitResult/jaredFit/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //Santona
  //TFile* inf = new TFile(Form("/home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/WorkingPlots/FitResults/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //Jaebeom
  //TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/samba/UpsilonAnalysis/fitResultFiles/mcFit_MuPt4_2016_11_04/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

void stripErrorBars( TH1* h, double defaultErr  ) {
  
  for ( int i=1;  i<= h->GetNbinsX() ; i++) {
    h->SetBinError( i, defaultErr);
  }
}
double getScale(int fTAA, double* TAA, double* centBin, int nCentBins)
{
  double flumi_;
  if(fTAA == nCentBins+1) flumi_ = 368;
  else if(centBin[fTAA-1]>=60 && centBin[fTAA-1]<120 && fTAA !=nCentBins+1) flumi_ = 464;
  else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 464;
  //else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 334.82249848;
  else if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1) flumi_ = 368;
  double nMBColl = NumberOfMBColl;
  if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1)  nMBColl = NumberOfMBColl;
  else if(centBin[fTAA-1]>=60 && fTAA !=nCentBins+1) nMBColl = NumberOfMBColl1;
  double scaleFactor;
  if(fTAA == nCentBins+1) scaleFactor = 28000000000000./(nMBColl*TAA[fTAA-1]*1000);
  else scaleFactor = 28000000000000./(nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200.*TAA[fTAA-1]*1000);
  
  if(fTAA!=nCentBins+1){
  cout << endl;
  cout << "icent : " << centBin[fTAA-1] << " - " << centBin[fTAA] << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "nMBColl*(centBin[fTAA]-centBin[fTAA-1])/100. : " << nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200. << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }

  else if(fTAA==nCentBins+1){
  cout << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }
  return scaleFactor;
}

