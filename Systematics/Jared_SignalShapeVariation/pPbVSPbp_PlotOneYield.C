#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "../../HeaderFiles/PsetCollection.h"
#include "../../HeaderFiles/CMS_lumi.C"
#include "../../HeaderFiles/tdrstyle.C"
#include "../../HeaderFiles/StyleSetting.h"


void pPbVSPbp_PlotOneYield(int collId=kPADATA, int whichUpsilon=1, int whichParam=1) {

  float scale = whichUpsilon;
  if (collId==kPPDATA) scale = scale*8;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetErrorX(0);

  //arrays of upper and lower limits. The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_sigma,err_mu,m_lambda,mass, nSig1s,nSig2s,nSig3s,nBkg}
  double paramsupper[13] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0, 9.56, 1000000,360000,260000,5000000};
  double paramslower[13] = {0.02, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.36, 0,-20,-50,0};
  double plotupper[13] = {0.6, 5, 5, 5.2, 2, 25, 25, 25, 9.6, 1600*scale,800*scale,800*scale,8000*scale};
  double plotlower[13] = {-0.1, 0, 0, -0.2, -1, 0, 0, 0, 9.35, 0,-20,-50,0};

  //Name of histogram:
  TString paramNames[13] = {"sigma1s","x1s","alpha","n1s","f1s","errsigma","errmu","lambda","mass", "nSig1s","nSig2s","nSig3s","nBkg"};

  TString ptHistoName = Form("h%spt",paramNames[whichParam].Data());
  TString yHistoName = Form("h%sy",paramNames[whichParam].Data());
  cout << ptHistoName << endl;
  cout << yHistoName << endl;

  //set up canvases
  TCanvas *cParam1 = new TCanvas("cParam1","cParam1",4,45,400,400);
  TCanvas *cParam2 = new TCanvas("cParam2","cParam2",4,45,400,400);

  //Extract histograms
  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";
  TString HistoFileName = Form("FittedParamsHistos_%s_%is_run%i.root",strId.Data(),whichUpsilon,3);
  cout << HistoFileName << endl;
  HistoFile = TFile::Open(HistoFileName,"READ");
  TH1F* hParampt = (TH1F*)HistoFile->Get(ptHistoName);
  TH1F* hParamy = (TH1F*)HistoFile->Get(yHistoName);

  TString NewHistoFileName = Form("FittedParamsHistos_%s_%is_run%i.root",strId.Data(),whichUpsilon,1);
  cout << NewHistoFileName << endl;
  NewHistoFile = TFile::Open(NewHistoFileName,"READ");
  TH1F* NewhParampt = (TH1F*)NewHistoFile->Get(ptHistoName);
  TH1F* NewhParamy = (TH1F*)NewHistoFile->Get(yHistoName);

  TString NewerHistoFileName = Form("FittedParamsHistos_%s_%is_run%i.root",strId.Data(),whichUpsilon,2);
  cout << NewerHistoFileName << endl;
  NewerHistoFile = TFile::Open(NewerHistoFileName,"READ");
  TH1F* NewerhParampt = (TH1F*)NewerHistoFile->Get(ptHistoName);
  TH1F* NewerhParamy = (TH1F*)NewerHistoFile->Get(yHistoName);

  hParampt->SetMarkerStyle(2);
  NewhParampt->SetMarkerStyle(2);
  NewerhParampt->SetMarkerStyle(2);

  hParamy->SetMarkerStyle(2);
  NewhParamy->SetMarkerStyle(2);
  NewerhParamy->SetMarkerStyle(2);

  //draw pt plots
  cParam1->cd();
  hParampt->SetXTitle("pT");
  hParampt->GetYaxis()->SetRangeUser(plotlower[whichParam],plotupper[whichParam]);
  hParampt->Draw();
  hParampt->SetMarkerColor(3);
  hParampt->SetLineColor(3);
  NewhParampt->Draw("same");
  NewhParampt->SetMarkerColor(6);
  NewhParampt->SetLineColor(6);
  NewerhParampt->Draw("same");
  NewerhParampt->SetMarkerColor(4);
  NewerhParampt->SetLineColor(4);
  TLine *nBkgptlineLow = new TLine(0,paramslower[whichParam],30,paramslower[whichParam]);
  TLine *nBkgptlineHigh = new TLine(0,paramsupper[whichParam],30,paramsupper[whichParam]);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();
  nBkgptlineHigh->SetLineColor(kRed);
  nBkgptlineHigh->Draw();
  TLegend *nBkgLegend1 = new TLegend(0.4,0.7,0.6,0.9);
  nBkgLegend1->AddEntry(hParampt,"combined","p");
  nBkgLegend1->AddEntry(NewhParampt,"run 1","p");
  nBkgLegend1->AddEntry(NewerhParampt,"run 2","p");
  nBkgLegend1->Draw("same");

  cParam2->cd();
  hParamy->SetXTitle("y");
  hParamy->GetYaxis()->SetRangeUser(plotlower[whichParam],plotupper[whichParam]);
  hParamy->Draw();
  hParamy->SetMarkerColor(3);
  hParamy->SetLineColor(3);
  NewhParamy->Draw("same");
  NewhParamy->SetMarkerColor(6);
  NewhParamy->SetLineColor(6);
  NewerhParamy->Draw("same");
  NewerhParamy->SetMarkerColor(4);
  NewerhParamy->SetLineColor(4);
  TLine *nBkgylineLow = new TLine(-1.93,paramslower[whichParam],1.93,paramslower[whichParam]);
  TLine *nBkgylineHigh = new TLine(-1.93,paramsupper[whichParam],1.93,paramsupper[whichParam]);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();
  nBkgylineHigh->SetLineColor(kRed);
  nBkgylineHigh->Draw();
  TLegend *nBkgLegend = new TLegend(0.4,0.7,0.6,0.9);
  nBkgLegend->AddEntry(hParamy,"combined","p");
  nBkgLegend->AddEntry(NewhParamy,"run 1","p");
  nBkgLegend->AddEntry(NewerhParamy,"run 2","p");
  nBkgLegend->Draw("same");

  //save plots
  cParam1->SaveAs(Form("ParameterPlots/%sfitted_%s_%isPtbins_pPbVSPbp.png",strId.Data(),paramNames[whichParam].Data(),whichUpsilon));
  cParam2->SaveAs(Form("ParameterPlots/%sfitted_%s_%isYbins_pPbVSPbp.png",strId.Data(),paramNames[whichParam].Data(),whichUpsilon));

}
