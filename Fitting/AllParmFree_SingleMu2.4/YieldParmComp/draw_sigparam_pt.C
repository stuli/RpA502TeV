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

int draw_sigparam_pt(TString szAA = "PP", int states =1)
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
  
  int tmpBin = nPtBins;
  
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double *binArr = ptBin; // array
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
  double alpha1s[nBin];
  double alpha1sErr[nBin];
  double n1s[nBin];
  double n1sErr[nBin];
  double sigma1s_1[nBin];
  double sigma1s_1Err[nBin];
  double f1s[nBin];
  double f1sErr[nBin];
  double x1s[nBin];
  double x1sErr[nBin];
  
  
  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    if (szAA == "PP" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt%.1f-%.1f_y-1.93-1.93_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else if (szAA == "PA" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt%.1f-%.1f_y-1.93-1.93_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
    //cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    //if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    //// get parameters
    //cout << ws[ib]->var("alpha1s")->getVal() << endl;
    alpha1s[ib]=ws[ib]->var("alpha1s_1")->getVal();
    alpha1sErr[ib]=ws[ib]->var("alpha1s_1")->getError();
    n1s[ib]=ws[ib]->var("n1s_1")->getVal();
    n1sErr[ib]=ws[ib]->var("n1s_1")->getError();
    sigma1s_1[ib]=ws[ib]->var("sigma1s_1")->getVal();
    sigma1s_1Err[ib]=ws[ib]->var("sigma1s_1")->getError();
    f1s[ib]=ws[ib]->var("f1s")->getVal();
    f1sErr[ib]=ws[ib]->var("f1s")->getError();
    x1s[ib]=ws[ib]->var("x1s")->getVal();
    x1sErr[ib]=ws[ib]->var("x1s")->getError();
  }
 
  //// histogram
  TH1D* h1_alpha1s = new TH1D("h1_alpha1s","h1_alpha1s;p_{T} (GeV/c);alpha1s",nBin,binArr); 
  TH1D* h1_n1s = new TH1D("h1_n1s","h1_n1s;p_{T} (GeV/c);n1s",nBin,binArr); 
  TH1D* h1_sigma1s_1 = new TH1D("h1_sigma1s_1","h1_sigma1s_1;p_{T} (GeV/c);sigma1s_1",nBin,binArr); 
  TH1D* h1_f1s = new TH1D("h1_f1s","h1_f1s;p_{T} (GeV/c);f1s",nBin,binArr); 
  TH1D* h1_x1s = new TH1D("h1_x1s","h1_x1s;p_{T} (GeV/c);x1s",nBin,binArr); 
  
  for (int ib =0; ib < nBin; ib ++ ) {
    h1_alpha1s->SetBinContent(ib+1,alpha1s[ib]);   
    h1_alpha1s->SetBinError(ib+1,alpha1sErr[ib]);   
    h1_n1s->SetBinContent(ib+1,n1s[ib]);   
    h1_n1s->SetBinError(ib+1,n1sErr[ib]);   
    h1_sigma1s_1->SetBinContent(ib+1,sigma1s_1[ib]);   
    h1_sigma1s_1->SetBinError(ib+1,sigma1s_1Err[ib]);   
    h1_f1s->SetBinContent(ib+1,f1s[ib]);   
    h1_f1s->SetBinError(ib+1,f1sErr[ib]);   
    h1_x1s->SetBinContent(ib+1,x1s[ib]);   
    h1_x1s->SetBinError(ib+1,x1sErr[ib]);   
  }

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
  
  TCanvas* c_alpha1s = new TCanvas("c_alpha1s","c_alpha1s",600,600);
  c_alpha1s->cd();
  h1_alpha1s->GetXaxis()->CenterTitle(1);
  h1_alpha1s->GetYaxis()->CenterTitle(1);
  h1_alpha1s->GetYaxis()->SetRangeUser(0,3);
  h1_alpha1s->SetMarkerColor(kRed);
  h1_alpha1s->SetLineColor(kRed);
  h1_alpha1s->SetMarkerStyle(kFullCircle);
  h1_alpha1s->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_alpha1s->SaveAs(Form("sigparam/pt_alpha1s_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_n1s = new TCanvas("c_n1s","c_n1s",600,600);
  c_n1s->cd();
  h1_n1s->GetXaxis()->CenterTitle(1);
  h1_n1s->GetYaxis()->CenterTitle(1);
  h1_n1s->GetYaxis()->SetRangeUser(0,6);
  h1_n1s->SetMarkerColor(kRed);
  h1_n1s->SetLineColor(kRed);
  h1_n1s->SetMarkerStyle(kFullCircle);
  h1_n1s->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_n1s->SaveAs(Form("sigparam/pt_n1s_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_sigma1s_1 = new TCanvas("c_sigma1s_1","c_sigma1s_1",600,600);
  c_sigma1s_1->cd();
  h1_sigma1s_1->GetXaxis()->CenterTitle(1);
  h1_sigma1s_1->GetYaxis()->CenterTitle(1);
  h1_sigma1s_1->GetYaxis()->SetRangeUser(0,0.2);
  h1_sigma1s_1->SetMarkerColor(kRed);
  h1_sigma1s_1->SetLineColor(kRed);
  h1_sigma1s_1->SetMarkerStyle(kFullCircle);
  h1_sigma1s_1->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_sigma1s_1->SaveAs(Form("sigparam/pt_sigma1s_1_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_f1s = new TCanvas("c_f1s","c_f1s",600,600);
  c_f1s->cd();
  h1_f1s->GetXaxis()->CenterTitle(1);
  h1_f1s->GetYaxis()->CenterTitle(1);
  h1_f1s->GetYaxis()->SetRangeUser(0,1.0);
  h1_f1s->SetMarkerColor(kRed);
  h1_f1s->SetLineColor(kRed);
  h1_f1s->SetMarkerStyle(kFullCircle);
  h1_f1s->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_f1s->SaveAs(Form("sigparam/pt_f1s_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_x1s = new TCanvas("c_x1s","c_x1s",600,600);
  c_x1s->cd();
  h1_x1s->GetXaxis()->CenterTitle(1);
  h1_x1s->GetYaxis()->CenterTitle(1);
  h1_x1s->GetYaxis()->SetRangeUser(0,4);
  h1_x1s->SetMarkerColor(kRed);
  h1_x1s->SetLineColor(kRed);
  h1_x1s->SetMarkerStyle(kFullCircle);
  h1_x1s->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_x1s->SaveAs(Form("sigparam/pt_x1s_%s_%dS.pdf",szAA.Data(),states));

  return 0;

}

