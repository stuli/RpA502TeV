#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "cutsAndBin.h"
#include "multiTreeUtil.h"
using namespace std;

TString ResultDir  = "nominalFits";


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float hflow=0, float hfhigh=400, int ntlow=0, int nhigh=400) ;

void getSpectra_RFB_ntracks(int state = 1 ) {  

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  int nPtBins=0;
  int nYBins=0;
  double* ptBin;
  double* yBin;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;
  double* rap_rfbmin;
  double* rap_rfbmax;
  //  double nColl1s[nCentBins] = {1819,1432,1005,606,349,186,90.7,40.1,7.67};
  //  double nPart1s[nCentBins] = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};
  double* HFBin;
  double* NtracksBin;
  int nHFBins=0;
  int nNtracksBins=0;
  int nYrange=0;

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S_cr;    yBin = yBin1S; 
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
    HFBin = HFBin1s;
    NtracksBin = NtracksBin1s;
    nHFBins=nHFBins1s;
    nNtracksBins=nNtracksBins1s;
    rap_rfbmin = rap_rfbmin1s;
    rap_rfbmax = rap_rfbmax1s;
    nYrange = nYrange1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S_cr;    yBin = yBin2S; 
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
    HFBin = HFBin2s;
    NtracksBin = NtracksBin1s;
    nHFBins=nHFBins2s;
    nNtracksBins=nNtracksBins2s;
    rap_rfbmin = rap_rfbmin2s;
    rap_rfbmax = rap_rfbmax2s;
    nYrange = nYrange2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S_cr;    yBin = yBin3S; 
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
    HFBin = HFBin3s;
    NtracksBin = NtracksBin1s;
    nHFBins=nHFBins3s;
    nNtracksBins=nNtracksBins3s;
    rap_rfbmin = rap_rfbmin3s;
    rap_rfbmax = rap_rfbmax3s;
    nYrange = nYrange3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  
  double HFMin = HFBin[0]; double HFMax = HFBin[nHFBins];
  double NtracksMin = NtracksBin[0]; double NtracksMax = NtracksBin[nNtracksBins];


  TH1D* hptSigPAF[10];
  TH1D* hptSigPAB[10]; 
  for(int i=0;i<10;i++){
   hptSigPAF[i] = new TH1D(Form("hptPAF%d",i),";Ntracks;",nNtracksBins,NtracksBin);
   hptSigPAB[i] = new TH1D(Form("hptPAB%d",i),";Ntracks;",nNtracksBins,NtracksBin);
  }
  
  // signals :
  TH1D* hRPAraw_pt[10];   // w/o efficiency correction
  TH1D* hRPA_pt[10];   // w efficiency correction

  //***Pt***
  TCanvas* c1 =  new TCanvas("c1","",400,400);
  for ( int iyr = 0; iyr<nYrange; iyr++){
    for ( int ipt = 1 ; ipt<=nNtracksBins ; ipt++) {
      valErr yieldPP;
      if(rap_rfbmin[iyr]==0) yieldPP = getYield(state, kPADATA, 0,30, -rap_rfbmax[iyr], rap_rfbmin[iyr], 0,200, 0,400,NtracksBin[ipt-1],NtracksBin[ipt]);
      else yieldPP = getYield(state, kPADATA, 0,30, -rap_rfbmax[iyr], -rap_rfbmin[iyr], 0,200, 0,400,NtracksBin[ipt-1],NtracksBin[ipt]);
      valErr yieldPA = getYield(state, kPADATA, 0,30,  rap_rfbmin[iyr],  rap_rfbmax[iyr], 0,200, 0,400,NtracksBin[ipt-1],NtracksBin[ipt]);
      hptSigPAF[iyr]->SetBinContent( ipt, yieldPA.val ) ;
      hptSigPAF[iyr]->SetBinError( ipt, yieldPA.err ) ;
      hptSigPAB[iyr]->SetBinContent( ipt, yieldPP.val ) ;
      hptSigPAB[iyr]->SetBinError( ipt, yieldPP.err ) ;
      hptSigPAF[iyr]->SetAxisRange(10,1e5,"Y");
      hptSigPAB[iyr]->SetAxisRange(10,1e5,"Y");
      handsomeTH1(hptSigPAF[iyr],2);
      handsomeTH1(hptSigPAB[iyr],2);
    }
    hRPAraw_pt[iyr] = (TH1D*)hptSigPAF[iyr]->Clone(Form("raa_vs_pt_raw_%d",iyr+1));
    hRPAraw_pt[iyr]->Divide(hptSigPAB[iyr]);
  }

  //~*~*~*~* Corrections ~*~*~*~*~*
  TFile *f_acc = new TFile(Form("Acceptance/acceptance_wgt_%dS_20171121.root",state),"read");
  TFile *f_eff = new TFile(Form("CrossChecks/efficiency_Santona/ForAnalysisNote/EffCor_SyspPbXS_%dS.root",state),"read");
  TH1D* hAcc = (TH1D*)f_acc->Get(Form("hrapAccXsPA%dS",state));
  TH1D* hEff = (TH1D*)f_eff->Get("EffNomRap");
  double corr_b;
  double corr_f;
  for(int iyr=0;iyr<nYrange;iyr++){
   corr_b = hEff->GetBinContent(nYrange+2-iyr);
   corr_f = hEff->GetBinContent(nYrange+2+1+iyr);
   hRPAraw_pt[iyr]->Scale(corr_b/corr_f);
   corr_b = hAcc->GetBinContent(nYrange+2-iyr);
   corr_f = hAcc->GetBinContent(nYrange+2+1+iyr);
   hRPAraw_pt[iyr]->Scale(corr_b/corr_f);
  }

  for(int ipt =1 ; ipt<=nNtracksBins; ipt++){
    valErr yieldIntB = getYield(state, kPADATA, 0,30, -1.93, 0.00, 0,200, 0,400,NtracksBin[ipt-1],NtracksBin[ipt]);
    valErr yieldIntF = getYield(state, kPADATA, 0,30, 0.00, 1.93, 0,200, 0,400,NtracksBin[ipt-1],NtracksBin[ipt]);
    hptSigPAB[9] -> SetBinContent(ipt, yieldIntB.val);
    hptSigPAF[9] -> SetBinContent(ipt, yieldIntF.val);
    hptSigPAB[9] -> SetBinError(ipt, yieldIntB.err);
    hptSigPAF[9] -> SetBinError(ipt, yieldIntF.err);
  }
  hptSigPAF[9]->Divide(hptSigPAB[9]);
  
  TFile* fEff_int = new TFile(Form("CrossChecks/efficiency_Santona/ForAnalysisNote/EffNomCor_Sys2DRpA_%dS.root",state),"read");
  TH1D* hEff_intF = (TH1D*) fEff_int -> Get("EffNomRatIntRapPos");
  TH1D* hEff_intB = (TH1D*) fEff_int -> Get("EffNomRatIntRapNeg");
  corr_b = hEff_intB->GetBinContent(1); 
  corr_f = hEff_intB->GetBinContent(1); 
  hRPAraw_pt[9] = (TH1D*)hptSigPAF[9]->Clone("raa_vs_pt_raw_9");
  hRPAraw_pt[9]->Scale(corr_b/corr_f);

  TGraphErrors *gRPA_pt[10];
  TFile *wf = new TFile(Form("finalResults/Ups_%d_RFB_Ntracks.root",state),"recreate");
  for(int iyr=0; iyr<nYrange;iyr++){
   gRPA_pt[iyr] = new TGraphErrors(hRPAraw_pt[iyr]);
   gRPA_pt[iyr]->SetName(Form("gRFB_%d",iyr+1));
   gRPA_pt[iyr]->SetTitle(Form("gRFB_vs_Ntracks_%d",iyr+1));
   gRPA_pt[iyr]->Write();
  }
  gRPA_pt[9] = new TGraphErrors(hRPAraw_pt[9]);
  gRPA_pt[9]->SetName("gRFB_int");
  gRPA_pt[9]->SetTitle("gRFB_vs_Ntrackl_int");
  gRPA_pt[9]->Write(); 
  wf->Close();
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh, float hflow, float hfhigh, int ntlow, int nthigh) {
  TString kineLabel = getKineLabel_eventAct (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, hflow, hfhigh, ntlow, nthigh) ;
  TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/NominalFitResult/jaredFit/NtracksFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/samba/UpsilonAnalysis/fitResultFiles/mcFit_MuPt4_2016_11_04/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

void stripErrorBars( TH1* h, double defaultErr  ) {
  
  for ( int i=0;  i<= h->GetNbinsX() ; i++) {
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

