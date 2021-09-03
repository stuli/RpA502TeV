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


void MakeParamsHistosYBins(int whichUpsilon=1, int whichPtRange=0, int collId=kPADATA) {

  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  if (whichUpsilon==1) {
    if (whichPtRange==0) {
      float ybins[11] = {-2.87,-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    }
    else {
      float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    }
  }
  else if (whichUpsilon==2) {
    if (whichPtRange==0) {
      float ybins[7] = {-2.87,-2.4,-1.93,-0.8,0.0,0.8,1.93};
    }
    else {
      float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
    }
  }
  else if (whichUpsilon==3) {
    if (whichPtRange==0) {
      float ybins[5] = {-2.87,-2.4,-1.93,0.0,1.93};
    }
    else {
      float ybins[3] = {-1.93,0.0,1.93};
    }
  }

  const int numybins = sizeof(ybins)/sizeof(float)-1;

  TString kineLabel, NomFileName, outfilename;
  TFile* NomFile;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow, yHigh;
  int iystart = 0;
  int iyend = numybins;
  if (whichPtRange==0) {
    ptLow = 0;
    ptHigh = 30;
  }
  else if (whichPtRange==1) {
    ptLow = 0;
    ptHigh = 6;
  }
  else if (whichPtRange==2) {
    ptLow = 6;
    ptHigh = 30;
  }
  if (collId==kPADATA) {
    iyend = numybins;
    outfilename = Form("FittedParamsHistos_PA_%is_YBins_PtRange%i.root",whichUpsilon,whichPtRange);
  }
  else if (collId==kPPDATA) {
    if (whichPtRange==0) {
      iystart = 1;
      iyend = numybins/2+1;
    }
    else {
      iystart = 0;
      iyend = numybins/2;
    }
    outfilename = Form("FittedParamsHistos_PP_%is_YBins_PtRange%i.root",whichUpsilon,whichPtRange);
  }

  TFile outFile(outfilename, "RECREATE");

  //declare histograms
  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybins);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybins);
  TH1F* hmassy = new TH1F("hmassy","mass vs y",numybins,ybins);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybins);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs y",numybins,ybins);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybins);
  TH1F* herrmuy = new TH1F("herrmuy","errmu vs y",numybins,ybins);
  TH1F* herrsigmay = new TH1F("herrsigmay","errsigma vs y",numybins,ybins);
  TH1F* hlambday = new TH1F("hlambday","errlambda vs y",numybins,ybins);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs y",numybins,ybins);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs y",numybins,ybins);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs y",numybins,ybins);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs y",numybins,ybins);

  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //y loop
  for (int iy = iystart; iy<iyend; iy++) {
    yLow = ybins[iy];
    yHigh = ybins[iy+1];
    if (collId==kPPDATA) {
      yLow = TMath::Abs(ybins[iy+1]);
      yHigh = TMath::Abs(ybins[iy]);
    }

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("OfficialNominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
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
    tempmass = yws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = yws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = yws->var("n1s_1")->getVal();  
    tempn1serr = yws->var("n1s_1")->getError();
    tempsigma1s = yws->var("sigma1s_1")->getVal();  
    tempsigma1serr = yws->var("sigma1s_1")->getError();
    tempx1s = yws->var("x1s")->getVal();  
    tempx1serr = yws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = yws->var("#mu")->getVal();  
      temperrmuerr = yws->var("#mu")->getError();
      temperrsigma = yws->var("#sigma")->getVal();  
      temperrsigmaerr = yws->var("#sigma")->getError();
    }
    templambda = yws->var("#lambda")->getVal();  
    templambdaerr = yws->var("#lambda")->getError();
    tempnSig1s = yws->var("nSig1s")->getVal();//ywidth;  
    tempnSig1serr = yws->var("nSig1s")->getError();
    tempnSig2s = yws->var("nSig2s")->getVal();//ywidth;  
    tempnSig2serr = yws->var("nSig2s")->getError();
    tempnSig3s = yws->var("nSig3s")->getVal();//ywidth;  
    tempnSig3serr = yws->var("nSig3s")->getError();
    tempnBkg = yws->var("nBkg")->getVal();//ywidth;  
    tempnBkgerr = yws->var("nBkg")->getError();

    //fill histograms
    halphay->SetBinContent(iy+1, tempalpha);
    halphay->SetBinError  (iy+1, tempalphaerr);
    hf1sy->SetBinContent(iy+1, tempf1s);
    hf1sy->SetBinError  (iy+1, tempf1serr);
    hmassy->SetBinContent(iy+1, tempmass);
    hmassy->SetBinError  (iy+1, tempmasserr);
    hn1sy->SetBinContent(iy+1, tempn1s);
    hn1sy->SetBinError  (iy+1, tempn1serr);
    hsigma1sy->SetBinContent(iy+1, tempsigma1s);
    hsigma1sy->SetBinError  (iy+1, tempsigma1serr);
    hx1sy->SetBinContent(iy+1, tempx1s);
    hx1sy->SetBinError  (iy+1, tempx1serr);
    if (ptLow<5) {
      herrmuy->SetBinContent(iy+1, temperrmu);
      herrmuy->SetBinError  (iy+1, temperrmuerr);
      herrsigmay->SetBinContent(iy+1, temperrsigma);
      herrsigmay->SetBinError  (iy+1, temperrsigmaerr);
    }
    hlambday->SetBinContent(iy+1, templambda);
    hlambday->SetBinError  (iy+1, templambdaerr);
    hnSig1sy->SetBinContent(iy+1, tempnSig1s);
    hnSig1sy->SetBinError  (iy+1, tempnSig1serr);
    hnSig2sy->SetBinContent(iy+1, tempnSig2s);
    hnSig2sy->SetBinError  (iy+1, tempnSig2serr);
    hnSig3sy->SetBinContent(iy+1, tempnSig3s);
    hnSig3sy->SetBinError  (iy+1, tempnSig3serr);
    hnBkgy->SetBinContent(iy+1, tempnBkg);
    hnBkgy->SetBinError  (iy+1, tempnBkgerr);
  }

  //save histograms
  outFile.cd();

  halphay->Write();
  hf1sy->Write();
  hmassy->Write();
  hn1sy->Write();
  hsigma1sy->Write();
  hx1sy->Write();
  herrmuy->Write();
  herrsigmay->Write();
  hlambday->Write();
  hnSig1sy->Write();
  hnSig2sy->Write();
  hnSig3sy->Write();
  hnBkgy->Write();

  outFile.Close();

}
