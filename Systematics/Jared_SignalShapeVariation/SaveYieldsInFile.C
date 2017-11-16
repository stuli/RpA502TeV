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


void SaveYieldsInFile() {

  int collId = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  TString outFileName = "PPNominalYields1s.root";

  //choose a set of bins
  int whichUpsilon = 1;

  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[3] = {-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");

  //set up canvases
  TCanvas *cnSig1s = new TCanvas("cnSig1s","cnSig1s",4,45,800,400);
  cnSig1s->Divide(2,1);
  TCanvas *cnSig2s = new TCanvas("cnSig2s","cnSig2s",4,45,800,400);
  cnSig2s->Divide(2,1);
  TCanvas *cnSig3s = new TCanvas("cnSig3s","cnSig3s",4,45,800,400);
  cnSig3s->Divide(2,1);

  //declare histograms
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs pt",numptbins,ptbins);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs pt",numptbins,ptbins);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs pt",numptbins,ptbins);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs pt",numptbins,ptbins);

  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs y",numybins,ybins);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs y",numybins,ybins);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs y",numybins,ybins);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  RooWorkspace *ws;

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    float ptLow = ptbins[ipt];
    float ptHigh = ptbins[ipt+1];
    float yLow = 0.0;
    float yHigh = 1.93;

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    ws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float tempnSig1s = ws->var("nSig1s")->getVal();  
    float tempnSig1serr = ws->var("nSig1s")->getError();
    float tempnSig2s = ws->var("nSig2s")->getVal();  
    float tempnSig2serr = ws->var("nSig2s")->getError();
    float tempnSig3s = ws->var("nSig3s")->getVal();  
    float tempnSig3serr = ws->var("nSig3s")->getError();

    //fill histograms
    hnSig1spt->SetBinContent(ipt+1, tempnSig1s);
    hnSig1spt->SetBinError  (ipt+1, tempnSig1serr);
    hnSig2spt->SetBinContent(ipt+1, tempnSig2s);
    hnSig2spt->SetBinError  (ipt+1, tempnSig2serr);
    hnSig3spt->SetBinContent(ipt+1, tempnSig3s);
    hnSig3spt->SetBinError  (ipt+1, tempnSig3serr);
  }

  //draw pt plots
  cnSig1s->cd(1);
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600);
  hnSig1spt->Draw();
  hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,30,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd(1);
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->GetYaxis()->SetRangeUser(0,1600);
  hnSig2spt->Draw();
  hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,0,30,0);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd(1);
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->GetYaxis()->SetRangeUser(0,1600);
  hnSig3spt->Draw();
  hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,0,30,0);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  
  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    float ptLow = 0;
    float ptHigh = 30;
    float yLowPA = ybins[iy];
    float yHighPA = ybins[iy+1];
    if (yLowPA<0) {
      yLow = TMath::Abs(yHighPA);
      yHigh = TMath::Abs(yLowPA);
    }
    else {
      yLow = yLowPA;
      yHigh = yHighPA;
    }
    

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    yws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float tempnSig1s = yws->var("nSig1s")->getVal();  
    float tempnSig1serr = yws->var("nSig1s")->getError();
    float tempnSig2s = yws->var("nSig2s")->getVal();  
    float tempnSig2serr = yws->var("nSig2s")->getError();
    float tempnSig3s = yws->var("nSig3s")->getVal();  
    float tempnSig3serr = yws->var("nSig3s")->getError();

    //fill histograms
    hnSig1sy->SetBinContent(iy+1, tempnSig1s);
    hnSig1sy->SetBinError  (iy+1, tempnSig1serr);
    hnSig2sy->SetBinContent(iy+1, tempnSig2s);
    hnSig2sy->SetBinError  (iy+1, tempnSig2serr);
    hnSig3sy->SetBinContent(iy+1, tempnSig3s);
    hnSig3sy->SetBinError  (iy+1, tempnSig3serr);
  }

  //draw rapidity plots
  cnSig1s->cd(2);
  hnSig1sy->SetXTitle("y");
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600);
  hnSig1sy->Draw();
  hnSig1sy->Fit("pol1");
  TLine *nSig1sylineLow = new TLine(-1.93,0,1.93,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();

  cnSig2s->cd(2);
  hnSig2sy->SetXTitle("y");
  hnSig2sy->GetYaxis()->SetRangeUser(0,1600);
  hnSig2sy->Draw();
  hnSig2sy->Fit("pol1");
  TLine *nSig2sylineLow = new TLine(-1.93,0,1.93,0);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();

  cnSig3s->cd(2);
  hnSig3sy->SetXTitle("y");
  hnSig3sy->GetYaxis()->SetRangeUser(0,1600);
  hnSig3sy->Draw();
  hnSig3sy->Fit("pol1");
  TLine *nSig3sylineLow = new TLine(-1.93,0,1.93,0);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();

  //save plots
  TFile outFile(outFileName, "RECREATE");
  hnSig1spt->Write();
  hnSig2spt->Write();
  hnSig3spt->Write();
  hnSig1sy->Write();
  hnSig2sy->Write();
  hnSig3sy->Write();
  outFile.Close();
}
