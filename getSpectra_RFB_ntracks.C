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

  TH1D* hptSigPAF = new TH1D("hNTPAF",";N_{tracks};",nNtracksBins,NtracksBin);
  TH1D* hptSigPAB = new TH1D("hNTPAB",";N_{tracks};",nNtracksBins,NtracksBin);
  

  // signals :
  TH1D* hRFB_HF; 

  //***Pt***
  TCanvas* c1 =  new TCanvas("c1","",400,400);
  for ( int ipt = 1 ; ipt<=nNtracksBins ; ipt++) {
    valErr yieldPP = getYield(state, kPADATA, 0,30, -1.93, 0.00, 0,200, 0,120, NtracksBin[ipt-1],NtracksBin[ipt]);
    valErr yieldPA = getYield(state, kPADATA, 0,30,  0.00, 1.93, 0,200, 0,120, NtracksBin[ipt-1],NtracksBin[ipt]);
    hptSigPAF->SetBinContent( ipt, yieldPA.val ) ;
    hptSigPAF->SetBinError( ipt, yieldPA.err ) ;
    hptSigPAB->SetBinContent( ipt, yieldPP.val ) ;
    hptSigPAB->SetBinError( ipt, yieldPP.err ) ;
    hptSigPAF->SetAxisRange(10,1e5,"Y");
    hptSigPAB->SetAxisRange(10,1e5,"Y");
    handsomeTH1(hptSigPAF,2);
    handsomeTH1(hptSigPAB,2);
  }
  hRFB_HF = (TH1D*)hptSigPAF->Clone("rfb_nt");
  hRFB_HF->Divide(hptSigPAB);



  //~*~*~*~* Corrections ~*~*~*~*~*
  TFile *f_acc = new TFile(Form("Acceptance/20180328/acceptance_wgt_%dS_20180328_2Dplot.root",state),"read");
  TFile *f_eff = new TFile(Form("CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_SysRFB_%dS.root",state),"read");

  TH1D* hAcc_div = (TH1D*) hptSigPAF->Clone("hAcc_div"); hAcc_div->Reset();
  TH1D* hEff_div = (TH1D*) hptSigPAF->Clone("hAcc_div"); hEff_div->Reset();
  TH1D* hAcc_r = (TH1D*) f_acc->Get(Form("hrapAccPA2bin_%dS",state));
  TH1D* hEff_r = (TH1D*) f_eff->Get("EffNomRatIntRFB");
  double corr_acc = hAcc_r->GetBinContent(1)/hAcc_r->GetBinContent(2);
  double corr_eff = hEff_r->GetBinContent(1);
  for(int i=1;i<=hAcc_div->GetNbinsX();i++){
    hAcc_div->SetBinContent(i,corr_acc);
    hEff_div->SetBinContent(i,corr_eff);
  }
  hRFB_HF->Multiply(hAcc_div);
  hRFB_HF->Multiply(hEff_div);


  TGraphErrors *gRFB_HF;
  TFile *wf = new TFile(Form("finalResults/Ups_%d_RFB_Ntracks.root",state),"recreate");
  gRFB_HF = new TGraphErrors(hRFB_HF);
  gRFB_HF->SetName("gRFB_NT");
  gRFB_HF->SetTitle("gRFB_vs_NT");
  gRFB_HF->Write();
  wf->Close();
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh, float hflow, float hfhigh, int ntlow, int nthigh) {
  TString kineLabel = getKineLabel_eventAct (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, hflow, hfhigh, ntlow, nthigh) ;
  TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/NominalFitResult/jaredFit/NominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
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

