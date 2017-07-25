#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>

#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>

using namespace std;
using namespace RooFit;

int draw_bkgparam_rap(TString szAA = "PP", int states =1)
{
  /////////////////////////////////////////////////////////
  //// set style
  /////////////////////////////////////////////////////////
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12) ; 
  
  /////////////////////////////////////////////////////////
  //// binning setting
  /////////////////////////////////////////////////////////
  
  double tmpArr1s[7] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  double tmpArr2s[4] = {0.0, 0.8, 1.6, 2.4};
  double tmpArr3s[3] = {0.0, 1.2, 2.4};

  int tmpBin;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin = 6;
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin = 3;
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin = 2;
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl; return 0;
  }
  
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double binArr[nArrNum]; // array
  cout << "nBin = " << nBin << endl;
  for (int ib =0; ib < nArrNum; ib ++ ) {
    if (states ==1) { binArr[ib] = tmpArr1s[ib]; }
    else if (states ==2) { binArr[ib] = tmpArr2s[ib]; }
    else if (states ==3) { binArr[ib] = tmpArr3s[ib]; }
    cout << ib <<"th bin = " << binArr[ib] << endl;
  }

  /////////////////////////////////////////////////////////
  //// Open RooDataFile
  /////////////////////////////////////////////////////////
 
  //file and ws
  TFile *fileIn[nBin];
  RooWorkspace* ws[nBin];
  // parameters 
  double lambda[nBin];
  double lambdaErr[nBin];
  double mu[nBin];
  double muErr[nBin];
  double sigma[nBin];
  double sigmaErr[nBin];
  
  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    if (szAA == "PP" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else if (szAA == "AA" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_%s_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
    //cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    //if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    //ws[ib]->Print();
    //// get parameters
    //cout << ws[ib]->var("#lambda")->getVal() << endl;
    lambda[ib]=ws[ib]->var("#lambda")->getVal();
    lambdaErr[ib]=ws[ib]->var("#lambda")->getError();
    mu[ib]=ws[ib]->var("#mu")->getVal();
    muErr[ib]=ws[ib]->var("#mu")->getError();
    sigma[ib]=ws[ib]->var("#sigma")->getVal();
    sigmaErr[ib]=ws[ib]->var("#sigma")->getError();
    //cout << ib << "th lambda = " << lambda[ib] << endl;
    //cout << ib << "th mu = " << mu[ib] << endl;
    //cout << ib << "th sigma = " << sigma[ib] << endl;
  }
 
  //// histogram
  TH1D* h1_lambda = new TH1D("h1_lambda","h1_lambda;|y|;#lambda",nBin,binArr); 
  TH1D* h1_mu = new TH1D("h1_mu","h1_mu;|y|;#mu",nBin,binArr); 
  TH1D* h1_sigma = new TH1D("h1_sigma","h1_sigma;|y|;#sigma",nBin,binArr); 
  
  for (int ib =0; ib < nBin; ib ++ ) {
    h1_lambda->SetBinContent(ib+1,lambda[ib]);   
    h1_lambda->SetBinError(ib+1,lambdaErr[ib]);   
    h1_mu->SetBinContent(ib+1,mu[ib]);   
    h1_mu->SetBinError(ib+1,muErr[ib]);   
    h1_sigma->SetBinContent(ib+1,sigma[ib]);   
    h1_sigma->SetBinError(ib+1,sigmaErr[ib]);   
  }
  
  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);

  TCanvas* c_lambda = new TCanvas("c_lambda","c_lambda",600,600);
  c_lambda->cd();
  h1_lambda->GetXaxis()->CenterTitle(1);
  h1_lambda->GetYaxis()->CenterTitle(1);
  h1_lambda->GetYaxis()->SetRangeUser(0,35);
  h1_lambda->SetMarkerColor(kRed);
  h1_lambda->SetLineColor(kRed);
  h1_lambda->SetMarkerStyle(kFullCircle);
  h1_lambda->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_lambda->SaveAs(Form("bkgparam/rap_lambda_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_mu = new TCanvas("c_mu","c_mu",600,600);
  c_mu->cd();
  h1_mu->GetXaxis()->CenterTitle(1);
  h1_mu->GetYaxis()->CenterTitle(1);
  h1_mu->GetYaxis()->SetRangeUser(0,12);
  h1_mu->SetMarkerColor(kRed);
  h1_mu->SetLineColor(kRed);
  h1_mu->SetMarkerStyle(kFullCircle);
  h1_mu->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_mu->SaveAs(Form("bkgparam/rap_mu_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_sigma = new TCanvas("c_sigma","c_sigma",600,600);
  c_sigma->cd();
  h1_sigma->GetXaxis()->CenterTitle(1);
  h1_sigma->GetYaxis()->CenterTitle(1);
  h1_sigma->GetYaxis()->SetRangeUser(0,2.5);
  h1_sigma->SetMarkerColor(kRed);
  h1_sigma->SetLineColor(kRed);
  h1_sigma->SetMarkerStyle(kFullCircle);
  h1_sigma->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_sigma->SaveAs(Form("bkgparam/rap_sigma_%s_%dS.pdf",szAA.Data(),states));

  return 0;

}

