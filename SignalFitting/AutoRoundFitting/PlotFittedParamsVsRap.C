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
#include "RoundsHeader.h"

void PlotFittedParamsVsRap(int whichUpsilon=1, int collId=kPADATA, int whichRound=R3b) {

  TString directory;
  if (whichRound==R0) directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithSimICs_AllParamFree_March15/";
  else if (whichRound==R1a) directory = "RoundFits_R1a/";
  else if (whichRound==R1b) directory = "RoundFits_R1b/";
  else if (whichRound==R2a) directory = "RoundFits_R2a/";
  else if (whichRound==R2b) directory = "RoundFits_R2b/";
  else if (whichRound==R3a) directory = "RoundFits_R3a/";
  else if (whichRound==R3b) directory = "RoundFits_R3b/";

  float scale = whichUpsilon;
  if (collId==kPPDATA) scale = scale*8;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double plotlimupper[8] = {0.3, 1.2, 5.5, 5.5, 1.2, 15.0, 15.0, 25.0};
  double plotlimlower[8] = {0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0};

  //choose a set of bins
  float ybins0[2] = {-1.93,1.93};
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ybins3[3] = {-1.93,0.0,1.93};
  int numybinstmp;
  float* ybinsptr;

  if (whichUpsilon==0) {
    ybinsptr = &ybins0[0];
    numybinstmp = sizeof(ybins0)/sizeof(float)-1;
  }
  else if (whichUpsilon==1) {
    ybinsptr = &ybins1[0];
    numybinstmp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ybinsptr = &ybins2[0];
    numybinstmp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ybinsptr = &ybins3[0];
    numybinstmp = sizeof(ybins3)/sizeof(float)-1;
  }

  const int numybins = numybinstmp;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");

  //set up canvases
  TCanvas *c1 = new TCanvas("c1","c1",4,45,1000,800);
  c1->Divide(2,2);

  //declare histograms
  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybinsptr);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybinsptr);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybinsptr);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybinsptr);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float ptLow, ptHigh, yLow, yHigh;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;


  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    ptLow = 0;
    ptHigh = 30;
    yLow = *(ybinsptr+iy);
    yHigh = *(ybinsptr+iy+1);
    if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(*(ybinsptr+iy+1));
      yHigh = TMath::Abs(*(ybinsptr+iy));
    }

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *yws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float ywidth = yHigh-yLow;
    if (collId==kPPDATA) ywidth = 2*ywidth;
    tempalpha = yws->var("alpha1s_1")->getVal();  
    tempalphaerr = yws->var("alpha1s_1")->getError();
    tempf1s = yws->var("f1s")->getVal();  
    tempf1serr = yws->var("f1s")->getError();
    tempn1s = yws->var("n1s_1")->getVal();  
    tempn1serr = yws->var("n1s_1")->getError();
    tempx1s = yws->var("x1s")->getVal();  
    tempx1serr = yws->var("x1s")->getError();

    //fill histograms
    halphay->SetBinContent(iy+1, tempalpha);
    halphay->SetBinError  (iy+1, tempalphaerr);
    hf1sy->SetBinContent(iy+1, tempf1s);
    hf1sy->SetBinError  (iy+1, tempf1serr);
    hn1sy->SetBinContent(iy+1, tempn1s);
    hn1sy->SetBinError  (iy+1, tempn1serr);
    hx1sy->SetBinContent(iy+1, tempx1s);
    hx1sy->SetBinError  (iy+1, tempx1serr);
    delete yws;
    NomFile->Close("R");
    delete NomFile;
  }

  //draw rapidity plots
  c1->cd(1);
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphay->Draw();
  //halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(-1.93,paramslower[2],1.93,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(-1.93,paramsupper[2],1.93,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  c1->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1sy->Draw();
  //hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(-1.93,paramslower[4],1.93,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(-1.93,paramsupper[4],1.93,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  c1->cd(3);
  hn1sy->SetXTitle("y");
  hn1sy->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1sy->Draw();
  //hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(-1.93,paramslower[3],1.93,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(-1.93,paramsupper[3],1.93,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  c1->cd(4);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1sy->Draw();
  //hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(-1.93,paramslower[1],1.93,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(-1.93,paramsupper[1],1.93,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

  c1->SaveAs(Form("SignalParamPlots_%s.png",roundLabel[whichRound].Data()));
}
