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

using namespace std;
void PlotHFYields(int whichUpsilon=3) {

  int collId = kPADATA;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = 0.0;
  float yHigh = 1.93;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;
  float hfLow = 0;
  float hfHigh = 400;
  float ntracksLow = 0;
  float ntracksHigh = 400;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //if (whichUpsilon==2) yLow=0.0;
  //if (whichUpsilon==3) yLow=0.0;

  float HFBins[4] = {0,12,19,120};

  const int numHFBins = sizeof(HFBins)/sizeof(float)-1;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetMarkerStyle(2);
  //gStyle->SetErrorX(0);

  //set up canvases
  TCanvas *cyield = new TCanvas("cyield","cyield",4,45,500,400);
  //cmass->Divide(2,1);

  TH1F* hyieldF = new TH1F("hyieldF",Form("Fitted Upsilon Yield"),numHFBins,HFBins);
  TH1F* hyieldB = new TH1F("hyieldB",Form("Fitted Upsilon(%iS) Yield",whichUpsilon),numHFBins,HFBins);

//HF Forward loop
  for (int i = 0; i<numHFBins; i++) {

    hfLow = HFBins[i];
    hfHigh = HFBins[i+1];

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, (int)ntracksLow, (int)ntracksHigh );
    TString FileName = Form("/home/jared/Desktop/CMS_Research/RpA502TeV/SignalFitting/FitsWithSimICs/nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float yield1 = PAnomws->var(Form("nSig%is",1))->getVal();  
    float yield1err = PAnomws->var(Form("nSig%is",1))->getError();
    float yield2 = PAnomws->var(Form("nSig%is",2))->getVal();  
    float yield2err = PAnomws->var(Form("nSig%is",2))->getError();
    float yield3 = PAnomws->var(Form("nSig%is",3))->getVal();  
    float yield3err = PAnomws->var(Form("nSig%is",3))->getError();
    float yield = yield3;  
    float yielderr = yield3err;
    PANomFile->Close("R");

    //fill histograms
    hyieldF->SetBinContent(i+1, yield);
    hyieldF->SetBinError(i+1, yielderr);
  }

//HF Backward loop
  yLow = -1.93;
  yHigh = 0.0;
  for (int i = 0; i<numHFBins; i++) {

    hfLow = HFBins[i];
    hfHigh = HFBins[i+1];

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, (int)ntracksLow, (int)ntracksHigh );
    TString FileName = Form("/home/jared/Desktop/CMS_Research/RpA502TeV/SignalFitting/FitsWithSimICs/nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float yield = PAnomws->var(Form("nSig%is",whichUpsilon))->getVal();  
    float yielderr = PAnomws->var(Form("nSig%is",whichUpsilon))->getError();
    PANomFile->Close("R");

    //fill histograms
    hyieldB->SetBinContent(i+1, yield);
    hyieldB->SetBinError(i+1, yielderr);
  }

  //draw plot
  cyield->cd();
  hyieldF->SetMinimum(0);
  hyieldF->SetXTitle("HF");
  hyieldF->GetYaxis()->SetRangeUser(0,300);
  //hyieldF->SetMarkerStyle(4);
  //hyieldF->SetMarkerSize(1);
  hyieldF->SetMarkerColor(kBlue);
  hyieldF->Draw();
  hyieldB->SetLineColor(3);
  hyieldB->Draw("same");

  TLegend* yieldLegend = new TLegend(0.7,0.1,0.9,0.3);
  yieldLegend->SetTextSize(16);
  yieldLegend->SetTextFont(43);
  yieldLegend->AddEntry("hyieldF","Forward","le");
  yieldLegend->AddEntry("hyieldB","Backward","le");
  yieldLegend->Draw("same");

  //save plots
  TString BinsString = Form("_Bins%i",(int)HFBins[0]);
  for(int j = 1; j<=numHFBins; j++) BinsString = BinsString + Form("-%i",(int)HFBins[j]);
  TString PicName = Form("HFYields_%is",whichUpsilon);
  PicName = PicName + BinsString + Form(".png");
  cyield->SaveAs(PicName);
}
