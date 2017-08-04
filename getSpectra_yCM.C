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

void getSpectra_yCM(int state = 3 ) {  

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
  //  double nColl1s[nCentBins] = {1819,1432,1005,606,349,186,90.7,40.1,7.67};
  //  double nPart1s[nCentBins] = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};



  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; 
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S; 
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S; 
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  


  TH1D* hrapSigPP = new TH1D("hrapPP",";y_{CM};",nYBins,yBin);
  TH1D* hrapSigPA = new TH1D("hrapPA",";y_{CM};",nYBins,yBin);
  

  // signals :
  TH1D* hRPAraw_rap;   // w/o efficiency correction
  TH1D* hRPA_rap;   // w efficiency correction

  //***Pt***
  TCanvas* c1 =  new TCanvas("c1","",400,400);
  for ( int irap = 1 ; irap<= nYBins ; irap++) {
    valErr yieldPP = getYield(state, kPPDATA, 0,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    valErr yieldPA = getYield(state, kPADATA, 0,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    hrapSigPA->SetBinContent( irap, yieldPA.val ) ;
    hrapSigPA->SetBinError( irap, yieldPA.err ) ;
    hrapSigPP->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP->SetBinError( irap, yieldPP.err/TMath::Sqrt(2) ) ;
  }

  //*****Pt Yield***** 
  hrapSigPP->SetAxisRange(10,1e5,"Y");
  hrapSigPA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hrapSigPP,2);
  handsomeTH1(hrapSigPA,2);
  hrapSigPP->SetMarkerStyle(24);
 
  c1->cd(); 
  hrapSigPP->Draw();
  gPad->SetLogy();
  hrapSigPA->Draw("same");

  hRPAraw_rap = (TH1D*)hrapSigPA->Clone("raa_vs_rap_raw");
  hRPAraw_rap->Divide(hrapSigPP);
  
  double Factor_Scale = lumi_pp*1000/(lumi_pa*208);
  
  hRPAraw_rap->Scale(Factor_Scale);

  TCanvas* cRapRPA =  new TCanvas("cRpA_rap","",400,400);
  hRPA_rap = (TH1D*)hRPAraw_rap->Clone("rpa_vs_rap");

  hRPA_rap->SetMarkerStyle(24);
  hRPA_rap->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_rap->GetYaxis()->SetTitle("R_{pPb} uncorrected");
  cRapRPA->cd();
  hRPA_rap->Draw("pe");
  jumSun(yMin,1,yMax,1,1,1);

  TGraphErrors *gRPA_rap = new TGraphErrors(hRPAraw_rap);
  gRPA_rap->SetName("gRPA_rap");
  gRPA_rap->SetTitle("rpa_vs_rap");

  TFile *wf = new TFile(Form("finalResults/Ups_%d_RPA_rap.root",state),"recreate");
  gRPA_rap->Write();
  wf->Close();
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh,   float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
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

