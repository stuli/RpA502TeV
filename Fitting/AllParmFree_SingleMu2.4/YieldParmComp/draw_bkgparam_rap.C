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
#include "../../../cutsAndBin.h"

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
  int nPtBins=0;
  int nYBins=0;
  double* ptBin;
  double* yBin;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;
  if ( states == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; 
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( states == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S; 
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( states == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S; 
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  int tmpBin = nYBins;
  
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double *binArr = yBin; // array
  cout << "nBin = " << nBin << endl;

  for (int ib =0; ib < nArrNum; ib ++ ) {
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
    if (szAA == "PP" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else if (szAA == "PA" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
    //cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    //if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    //// get parameters
    //cout << ws[ib]->var("alpha1s")->getVal() << endl;
    lambda[ib]=ws[ib]->var("#lambda")->getVal();
    lambdaErr[ib]=ws[ib]->var("#lambda")->getError();
    if ( (states==1 && ib < 3) || (states==2 && ib < 2) || (states==3 && ib < 1) ) {
      mu[ib]=ws[ib]->var("#mu")->getVal();
      muErr[ib]=ws[ib]->var("#mu")->getError();
      sigma[ib]=ws[ib]->var("#sigma")->getVal();
      sigmaErr[ib]=ws[ib]->var("#sigma")->getError();
    } else {
      mu[ib] = 0;
      muErr[ib] = 0;
      sigma[ib] = 0;
      sigmaErr[ib] = 0;
    }
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

