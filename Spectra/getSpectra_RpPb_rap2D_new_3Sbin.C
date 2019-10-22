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
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0) ;

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 

void getSpectra_RpPb_rap2D_new_3Sbin(int state = 1) {  

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
    nYBins = nYBins3S_2D;    yBin = yBin3S_2D; nYBins_cr = nYBins1S_cr;   yBin_cr = yBin1S_cr;
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins3S_2D;    yBin = yBin3S_2D;  nYBins_cr = nYBins2S_cr;   yBin_cr = yBin2S_cr;
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S_2D;    yBin = yBin3S_2D;   nYBins_cr = nYBins3S_cr;   yBin_cr = yBin3S_cr;
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  

  // signals :
  TH1D* hRPAraw_rap_low;   // w/o efficiency correction
  TH1D* hRPAraw_rap_high;   // w/o efficiency correction
  TH1D* hRPA_rap_low;   // w efficiency correction
  TH1D* hRPA_rap_high;   // w efficiency correction

  TH1D* hrapAccPA_low;
  TH1D* hrapAccPP_low;
  TH1D* hrapAccPA_high;
  TH1D* hrapAccPP_high;

  TH1D* hrapEff_low;
  TH1D* hrapEff_high;


  TFile* infacc = new TFile(Form("Acceptance/20180724/acceptance_wgt_%dS_20180724_2Dplot.root",state),"read");
  hrapAccPP_low  = (TH1D*)infacc->Get(Form("hrapAccPP2BinPt1_%dS",state));
  hrapAccPA_low  = (TH1D*)infacc->Get(Form("hrapAccPA2BinPt1_%dS",state));
  hrapAccPP_high  = (TH1D*)infacc->Get(Form("hrapAccPP2BinPt2_%dS",state));
  hrapAccPA_high  = (TH1D*)infacc->Get(Form("hrapAccPA2BinPt2_%dS",state));

  TFile* infeff = new TFile(Form("Efficiency/RootFiles/EffNomCor_Sys3Sbins2DRpA_%dS.root",state),"read");
  hrapEff_low = (TH1D*)infeff->Get("EffNomRatRapLowpT");
  hrapEff_high = (TH1D*)infeff->Get("EffNomRatRapHighpT");

  stripErrorBars(hrapAccPP_low);
  stripErrorBars(hrapAccPA_low);
  stripErrorBars(hrapAccPP_high);
  stripErrorBars(hrapAccPA_high);
  stripErrorBars(hrapEff_low);
  stripErrorBars(hrapEff_high);
  
  TH1D* hrapSigPP_low = (TH1D*) hrapAccPP_low -> Clone("hrapPP_low");
  TH1D* hrapSigPA_low = (TH1D*) hrapAccPA_low -> Clone("hrapPA_low");
  TH1D* hrapSigPP_high = (TH1D*) hrapAccPP_high -> Clone("hrapPP_high");
  TH1D* hrapSigPA_high = (TH1D*) hrapAccPA_high -> Clone("hrapPA_high");

  hrapSigPP_low->Reset();
  hrapSigPA_low->Reset();
  hrapSigPP_high->Reset();
  hrapSigPA_high->Reset();


  //***yCM***
  valErr yieldPP;
  TCanvas* c_rap =  new TCanvas("c_rap","",400,400);

  for ( int irap = 1 ; irap<= nYBins ; irap++) {
    valErr yieldPP;
    if(irap <= (nYBins/2)) {
      yieldPP = getYield(state, kPPDATA, 0,6, TMath::Abs(yBin[irap]), TMath::Abs(yBin[irap-1]), 0,200,0,100);
    }
    else if(irap > (nYBins/2)) {
      cout << "irap : " << irap << endl;
      yieldPP = getYield(state, kPPDATA, 0,6, yBin[irap-1], yBin[irap], 0,200,0,100);
    }
    valErr yieldPA = getYield(state, kPADATA, 0,6, yBin[irap-1], yBin[irap], 0,200,0,100);
    hrapSigPP_low->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP_low->SetBinError( irap, yieldPP.err/2 ) ;
    hrapSigPA_low->SetBinContent( irap, yieldPA.val ) ;
    hrapSigPA_low->SetBinError( irap, yieldPA.err ) ;
  }

  for ( int irap = 1 ; irap<= nYBins ; irap++) {
    valErr yieldPP;
    if(irap <= (nYBins/2)) {
      yieldPP = getYield(state, kPPDATA, 6,30, TMath::Abs(yBin[irap]), TMath::Abs(yBin[irap-1]), 0,200,0,100);
    }
    else if(irap > (nYBins/2)) {
      cout << "irap : " << irap << endl;
      yieldPP = getYield(state, kPPDATA, 6,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    }
    valErr yieldPA = getYield(state, kPADATA, 6,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    hrapSigPP_high->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP_high->SetBinError( irap, yieldPP.err/2 ) ;
    hrapSigPA_high->SetBinContent( irap, yieldPA.val ) ;
    hrapSigPA_high->SetBinError( irap, yieldPA.err ) ;
  }

  //yCM Yield 
  hrapSigPP_low->SetAxisRange(10,1e5,"Y");
  hrapSigPA_low->SetAxisRange(10,1e5,"Y");
  hrapSigPP_high->SetAxisRange(10,1e5,"Y");
  hrapSigPA_high->SetAxisRange(10,1e5,"Y");
  handsomeTH1(hrapSigPP_low,2);
  handsomeTH1(hrapSigPA_low,2);
  handsomeTH1(hrapSigPP_high,2);
  handsomeTH1(hrapSigPA_high,2);
  hrapSigPP_low->SetMarkerStyle(24);
  hrapSigPP_high->SetMarkerStyle(24);
  hrapSigPA_low->SetMarkerStyle(24);
  hrapSigPA_high->SetMarkerStyle(24);
  
  double Factor_Scale = lumi_pp*1000/(lumi_pa*208);

  hRPAraw_rap_low = (TH1D*)hrapSigPA_low->Clone("rpa_vs_rap_raw_low");
  hRPAraw_rap_low -> Divide(hrapSigPP_low);
  hRPAraw_rap_high = (TH1D*)hrapSigPA_high->Clone("rpa_vs_rap_raw_high");
  hRPAraw_rap_high -> Divide(hrapSigPP_high);

  TH1D* hrel_Acc_rap_low = (TH1D*) hrapAccPA_low -> Clone("hrel_Acc_rap_low");
  TH1D* hrel_Eff_rap_low = (TH1D*) hrapEff_low -> Clone("hrel_Eff_rap_low");
  TH1D* hrel_Acc_rap_high = (TH1D*) hrapAccPA_high -> Clone("hrel_Acc_rap_high");
  TH1D* hrel_Eff_rap_high = (TH1D*) hrapEff_high -> Clone("hrel_Eff_rap_high");

  hrel_Acc_rap_low->Divide(hrapAccPP_low);
  hRPAraw_rap_low->Divide(hrel_Acc_rap_low);
  hRPAraw_rap_low->Multiply(hrel_Eff_rap_low);
  hRPAraw_rap_low->Scale(Factor_Scale);

  hrel_Acc_rap_high->Divide(hrapAccPP_high);
  hRPAraw_rap_high->Divide(hrel_Acc_rap_high);
  hRPAraw_rap_high->Multiply(hrel_Eff_rap_high);
  hRPAraw_rap_high->Scale(Factor_Scale);


  //***** Save ******

  //***** rap *****
  TCanvas* cRapRPA_low =  new TCanvas("cRpA_rap_low","",400,400);
  hRPA_rap_low = (TH1D*)hRPAraw_rap_low -> Clone("rpa_vs_rap_low");
  TCanvas* cRapRPA_high =  new TCanvas("cRpA_rap_high","",400,400);
  hRPA_rap_high = (TH1D*)hRPAraw_rap_high -> Clone("rpa_vs_rap_high");

  hRPA_rap_low -> SetMarkerStyle(24);
  hRPA_rap_low -> GetYaxis()->SetRangeUser(0,1.5);
  hRPA_rap_low -> GetYaxis()->SetTitle("R_{pPb} lowPt corrected");

  hRPA_rap_high -> SetMarkerStyle(24);
  hRPA_rap_high -> GetYaxis()->SetRangeUser(0,1.5);
  hRPA_rap_high -> GetYaxis()->SetTitle("R_{pPb} highPt corrected");

  cRapRPA_low->cd();
  hRPA_rap_low -> Draw("pe");
  jumSun(yMin,1,yMax,1,1,1);

  TGraphErrors *gRPA_rap_low = new TGraphErrors(hRPAraw_rap_low);
  gRPA_rap_low->SetName("gRPA_rap_low");
  gRPA_rap_low->SetTitle("rpa_vs_rap_low");

  cRapRPA_high->cd();
  hRPA_rap_high -> Draw("pe");
  jumSun(yMin,1,yMax,1,1,1);

  TGraphErrors *gRPA_rap_high = new TGraphErrors(hRPAraw_rap_high);
  gRPA_rap_high->SetName("gRPA_rap_high");
  gRPA_rap_high->SetTitle("rpa_vs_rap_high");


  TFile *wf = new TFile(Form("finalResults/Ups_%d_RPA_2D_rap_3Sbin.root",state),"recreate");
  gRPA_rap_low->Write();
  gRPA_rap_high->Write();
  wf->Close();
}

  valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh,   float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TFile* inf = new TFile(Form("NominalFitResult/jaredFit/NominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
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

