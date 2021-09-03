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


void PlotChiSquare(int collId=kPPDATA, int whichUpsilon=1) {

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

  //set up canvases
  TCanvas *cParam1 = new TCanvas("cParam1","cParam1",4,45,400,400);
  TCanvas *cParam2 = new TCanvas("cParam2","cParam2",4,45,400,400);

  TString ptHistoName = "hchisqpt";
  TString yHistoName = "hchisqy";

  //Extract histograms
  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";
  TString HistoFileName = Form("ChisqHisto_%s_%is_fixedn%i.root",strId.Data(),whichUpsilon,3);
  cout << HistoFileName << endl;
  HistoFile = TFile::Open(HistoFileName,"READ");
  TH1F* hParampt = (TH1F*)HistoFile->Get(ptHistoName);
  TH1F* hParamy = (TH1F*)HistoFile->Get(yHistoName);

  TString NewHistoFileName = Form("ChisqHisto_%s_%is_fixedn%i.root",strId.Data(),whichUpsilon,2);
  cout << NewHistoFileName << endl;
  NewHistoFile = TFile::Open(NewHistoFileName,"READ");
  TH1F* NewhParampt = (TH1F*)NewHistoFile->Get(ptHistoName);
  TH1F* NewhParamy = (TH1F*)NewHistoFile->Get(yHistoName);

  TString NewerHistoFileName = Form("ChisqHisto_%s_%is_fixedn%i.root",strId.Data(),whichUpsilon,4);
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
  hParampt->GetYaxis()->SetRangeUser(0,4);
  hParampt->Draw("p");
  hParampt->SetMarkerColor(3);
  hParampt->SetLineColor(3);
  NewhParampt->Draw("same p");
  NewhParampt->SetMarkerColor(6);
  NewhParampt->SetLineColor(6);
  NewerhParampt->Draw("same p");
  NewerhParampt->SetMarkerColor(4);
  NewerhParampt->SetLineColor(4);
  TLegend *nBkgLegend1 = new TLegend(0.4,0.7,0.6,0.9);
  nBkgLegend1->AddEntry(hParampt,"n=2","p");
  nBkgLegend1->AddEntry(NewhParampt,"n=3","p");
  nBkgLegend1->AddEntry(NewerhParampt,"n=4","p");
  nBkgLegend1->Draw("same");

  cParam2->cd();
  hParamy->SetXTitle("y");
  hParamy->GetYaxis()->SetRangeUser(0,4);
  hParamy->Draw("p");
  hParamy->SetMarkerColor(3);
  hParamy->SetLineColor(3);
  NewhParamy->Draw("same p");
  NewhParamy->SetMarkerColor(6);
  NewhParamy->SetLineColor(6);
  NewerhParamy->Draw("same p");
  NewerhParamy->SetMarkerColor(4);
  NewerhParamy->SetLineColor(4);
  TLegend *nBkgLegend = new TLegend(0.4,0.7,0.6,0.9);
  nBkgLegend->AddEntry(hParamy,"n=2","p");
  nBkgLegend->AddEntry(NewhParamy,"n=3","p");
  nBkgLegend->AddEntry(NewerhParamy,"n=4","p");
  nBkgLegend->Draw("same");

  //save plots
  TString paramName = "chisq";
  cParam1->SaveAs(Form("ParameterPlots/%sfitted_%s_%isPtbins_compare.png",strId.Data(),paramName.Data(),whichUpsilon));
  cParam2->SaveAs(Form("ParameterPlots/%sfitted_%s_%isYbins_compare.png",strId.Data(),paramName.Data(),whichUpsilon));

}
