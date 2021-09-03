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


void CompareFittedParams(int whichUpsilon=1, int collId=kPADATA) {

  TString newdirectory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";
  TString olddirectory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";
  float scale = whichUpsilon;
  TString yTitle = "y";
  if (collId==kPPDATA) {
    scale = scale*8;
    yTitle = "|y|";
  }
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  double plotlimupper[8] = {0.3, 1.2, 5.5, 5.5, 1.2, 16.0, 16.0, 26.0};
  double plotlimlower[8] = {0.0, -0.2, 0.5, 0.5, -0.2, -1.0, -1.0, -1.0};

  //Name of histogram:
  TString paramNames[13] = {"sigma1s","x1s","alpha","n1s","f1s","errmu","errsigma","lambda","mass", "nSig1s","nSig2s","nSig3s","nBkg"};

  //choose a set of bins
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[4] = {-2.87,-1.93,0.0,1.93};
  float ybinsPP1[5] = {0.0,0.4,0.8,1.2,1.93};
  float ybinsPP2[3] = {0.0,0.8,1.93};
  float ybinsPP3[2] = {0.0,1.93};

  float* ptbinsptr;
  float* ybinsptr;
  int numptbinstemp, numybinstemp;

if (collId==kPADATA) {
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    ybinsptr = &ybins1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
    numybinstemp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    ybinsptr = &ybins2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
    numybinstemp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    ybinsptr = &ybins3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
    numybinstemp = sizeof(ybins3)/sizeof(float)-1;
  }
}
else if (collId==kPPDATA) {
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    ybinsptr = &ybinsPP1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
    numybinstemp = sizeof(ybinsPP1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    ybinsptr = &ybinsPP2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
    numybinstemp = sizeof(ybinsPP2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    ybinsptr = &ybinsPP3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
    numybinstemp = sizeof(ybinsPP3)/sizeof(float)-1;
  }
}
  const int numptbins = numptbinstemp;
  const int numybins = numybinstemp;
  const int numtot = numptbins + numybins;

  float ymin = *ybinsptr;
  float ymax = *(ybinsptr+numybins);

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetMarkerStyle(7);
  gStyle->SetMarkerColor(4);

  //set up canvases
  TCanvas *calpha = new TCanvas("calpha","calpha",4,45,400,400);
  TCanvas *cf1s = new TCanvas("cf1s","cf1s",4,45,400,400);
  TCanvas *cmass = new TCanvas("cmass","cmass",4,45,400,400);
  TCanvas *cn1s = new TCanvas("cn1s","cn1s",4,45,400,400);
  TCanvas *csigma1s = new TCanvas("csigma1s","csigma1s",4,45,400,400);
  TCanvas *cx1s = new TCanvas("cx1s","cx1s",4,45,400,400);
  TCanvas *cerrmu = new TCanvas("cerrmu","cerrmu",4,45,400,400);
  TCanvas *cerrsigma = new TCanvas("cerrsigma","cerrsigma",4,45,400,400);
  TCanvas *clambda = new TCanvas("clambda","clambda",4,45,400,400);
  TCanvas *cnSig1s = new TCanvas("cnSig1s","cnSig1s",4,45,400,400);
  TCanvas *cnSig2s = new TCanvas("cnSig2s","cnSig2s",4,45,400,400);
  TCanvas *cnSig3s = new TCanvas("cnSig3s","cnSig3s",4,45,400,400);
  TCanvas *cnBkg = new TCanvas("cnBkg","cnBkg",4,45,400,400);

  TCanvas *calpha2 = new TCanvas("calpha2","calpha2",4,45,400,400);
  TCanvas *cf1s2 = new TCanvas("cf1s2","cf1s2",4,45,400,400);
  TCanvas *cmass2 = new TCanvas("cmass2","cmass2",4,45,400,400);
  TCanvas *cn1s2 = new TCanvas("cn1s2","cn1s2",4,45,400,400);
  TCanvas *csigma1s2 = new TCanvas("csigma1s2","csigma1s2",4,45,400,400);
  TCanvas *cx1s2 = new TCanvas("cx1s2","cx1s2",4,45,400,400);
  TCanvas *cerrmu2 = new TCanvas("cerrmu2","cerrmu2",4,45,400,400);
  TCanvas *cerrsigma2 = new TCanvas("cerrsigma2","cerrsigma2",4,45,400,400);
  TCanvas *clambda2 = new TCanvas("clambda2","clambda2",4,45,400,400);
  TCanvas *cnSig1s2 = new TCanvas("cnSig1s2","cnSig1s2",4,45,400,400);
  TCanvas *cnSig2s2 = new TCanvas("cnSig2s2","cnSig2s2",4,45,400,400);
  TCanvas *cnSig3s2 = new TCanvas("cnSig3s2","cnSig3s2",4,45,400,400);
  TCanvas *cnBkg2 = new TCanvas("cnBkg2","cnBkg2",4,45,400,400);


  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs pt",numptbins,ptbinsptr);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs pt",numptbins,ptbinsptr);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs pt",numptbins,ptbinsptr);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs pt",numptbins,ptbinsptr);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs pt",numptbins,ptbinsptr);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs pt",numptbins,ptbinsptr);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs pt",numptbins,ptbinsptr);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs pt",numptbins,ptbinsptr);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs pt",numptbins,ptbinsptr);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs pt",numptbins,ptbinsptr);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs pt",numptbins,ptbinsptr);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs pt",numptbins,ptbinsptr);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs pt",numptbins,ptbinsptr);

  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybinsptr);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybinsptr);
  TH1F* hmassy = new TH1F("hmassy","mass vs y",numybins,ybinsptr);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybinsptr);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs y",numybins,ybinsptr);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybinsptr);
  TH1F* herrmuy = new TH1F("herrmuy","errmu vs y",numybins,ybinsptr);
  TH1F* herrsigmay = new TH1F("herrsigmay","errsigma vs y",numybins,ybinsptr);
  TH1F* hlambday = new TH1F("hlambday","errlambda vs y",numybins,ybinsptr);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs y",numybins,ybinsptr);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs y",numybins,ybinsptr);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs y",numybins,ybinsptr);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs y",numybins,ybinsptr);

  TH1F* oldhalphapt = new TH1F("oldhalphapt","alpha vs pt",numptbins,ptbinsptr);
  TH1F* oldhf1spt = new TH1F("oldhf1spt","f1s vs pt",numptbins,ptbinsptr);
  TH1F* oldhmasspt = new TH1F("oldhmasspt","mass vs pt",numptbins,ptbinsptr);
  TH1F* oldhn1spt = new TH1F("oldhn1spt","n1s vs pt",numptbins,ptbinsptr);
  TH1F* oldhsigma1spt = new TH1F("oldhsigma1spt","sigma vs pt",numptbins,ptbinsptr);
  TH1F* oldhx1spt = new TH1F("oldhx1spt","x1s vs pt",numptbins,ptbinsptr);
  TH1F* oldherrmupt = new TH1F("oldherrmupt","errmu vs pt",numptbins,ptbinsptr);
  TH1F* oldherrsigmapt = new TH1F("oldherrsigmapt","errsigma vs pt",numptbins,ptbinsptr);
  TH1F* oldhlambdapt = new TH1F("oldhlambdapt","errlambda vs pt",numptbins,ptbinsptr);
  TH1F* oldhnSig1spt = new TH1F("oldhnSig1spt","nSig1s vs pt",numptbins,ptbinsptr);
  TH1F* oldhnSig2spt = new TH1F("oldhnSig2spt","nSig2s vs pt",numptbins,ptbinsptr);
  TH1F* oldhnSig3spt = new TH1F("oldhnSig3spt","nSig3s vs pt",numptbins,ptbinsptr);
  TH1F* oldhnBkgpt = new TH1F("oldhnBkgpt","nBkg vs pt",numptbins,ptbinsptr);

  TH1F* oldhalphay = new TH1F("oldhalphay","alpha vs y",numybins,ybinsptr);
  TH1F* oldhf1sy = new TH1F("oldhf1sy","f1s vs y",numybins,ybinsptr);
  TH1F* oldhmassy = new TH1F("oldhmassy","mass vs y",numybins,ybinsptr);
  TH1F* oldhn1sy = new TH1F("oldhn1sy","n1s vs y",numybins,ybinsptr);
  TH1F* oldhsigma1sy = new TH1F("oldhsigma1sy","sigma vs y",numybins,ybinsptr);
  TH1F* oldhx1sy = new TH1F("oldhx1sy","x1s vs y",numybins,ybinsptr);
  TH1F* oldherrmuy = new TH1F("oldherrmuy","errmu vs y",numybins,ybinsptr);
  TH1F* oldherrsigmay = new TH1F("oldherrsigmay","errsigma vs y",numybins,ybinsptr);
  TH1F* oldhlambday = new TH1F("oldhlambday","errlambda vs y",numybins,ybinsptr);
  TH1F* oldhnSig1sy = new TH1F("oldhnSig1sy","nSig1s vs y",numybins,ybinsptr);
  TH1F* oldhnSig2sy = new TH1F("oldhnSig2sy","nSig2s vs y",numybins,ybinsptr);
  TH1F* oldhnSig3sy = new TH1F("oldhnSig3sy","nSig3s vs y",numybins,ybinsptr);
  TH1F* oldhnBkgy = new TH1F("oldhnBkgy","nBkg vs y",numybins,ybinsptr);

  TString kineLabel, NomFileName, OldFileName;
  TFile* NomFile;
  TFile* OldFile;
  float ptLow, ptHigh, yLow, yHigh;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = *(ptbinsptr+ipt);
    ptHigh = *(ptbinsptr+ipt+1);
    if (collId==kPADATA) {
      yLow = -1.93;
      yHigh = 1.93;
    }
    else if (collId==kPPDATA) {
      yLow = 0.0;
      yHigh = 1.93;
    }

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",newdirectory.Data(),kineLabel.Data());
    cout << NomFileName << endl;
    NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");  

    //extract parameter values
    float ywidth = yHigh-yLow;
    if (collId==kPPDATA) ywidth = 2*ywidth;
    tempalpha = ws->var("alpha1s_1")->getVal();  
    tempalphaerr = ws->var("alpha1s_1")->getError();
    tempf1s = ws->var("f1s")->getVal();  
    tempf1serr = ws->var("f1s")->getError();
    tempmass = ws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = ws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = ws->var("n1s_1")->getVal();  
    tempn1serr = ws->var("n1s_1")->getError();
    tempsigma1s = ws->var("sigma1s_1")->getVal();  
    tempsigma1serr = ws->var("sigma1s_1")->getError();
    tempx1s = ws->var("x1s")->getVal();  
    tempx1serr = ws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = ws->var("#mu")->getVal();  
      temperrmuerr = ws->var("#mu")->getError();
      temperrsigma = ws->var("#sigma")->getVal();  
      temperrsigmaerr = ws->var("#sigma")->getError();
    }
    else {
      temperrmu = 0;
      temperrmuerr = 0;
      temperrsigma = 0;  
      temperrsigmaerr = 0;
    }
    templambda = ws->var("#lambda")->getVal();  
    templambdaerr = ws->var("#lambda")->getError();
    tempnSig1s = ws->var("nSig1s")->getVal();//ywidth;  
    tempnSig1serr = ws->var("nSig1s")->getError();
    tempnSig2s = ws->var("nSig2s")->getVal();//ywidth;  
    tempnSig2serr = ws->var("nSig2s")->getError();
    tempnSig3s = ws->var("nSig3s")->getVal();//ywidth;  
    tempnSig3serr = ws->var("nSig3s")->getError();
    tempnBkg = ws->var("nBkg")->getVal();//ywidth;  
    tempnBkgerr = ws->var("nBkg")->getError();

    //fill histograms
    halphapt->SetBinContent(ipt+1, tempalpha);
    halphapt->SetBinError  (ipt+1, tempalphaerr);
    hf1spt->SetBinContent(ipt+1, tempf1s);
    hf1spt->SetBinError  (ipt+1, tempf1serr);
    hmasspt->SetBinContent(ipt+1, tempmass);
    hmasspt->SetBinError  (ipt+1, tempmasserr);
    hn1spt->SetBinContent(ipt+1, tempn1s);
    hn1spt->SetBinError  (ipt+1, tempn1serr);
    hsigma1spt->SetBinContent(ipt+1, tempsigma1s);
    hsigma1spt->SetBinError  (ipt+1, tempsigma1serr);
    hx1spt->SetBinContent(ipt+1, tempx1s);
    hx1spt->SetBinError  (ipt+1, tempx1serr);
    if (ptLow<5) {
      herrmupt->SetBinContent(ipt+1, temperrmu);
      herrmupt->SetBinError  (ipt+1, temperrmuerr);
      herrsigmapt->SetBinContent(ipt+1, temperrsigma);
      herrsigmapt->SetBinError  (ipt+1, temperrsigmaerr);
    }
    hlambdapt->SetBinContent(ipt+1, templambda);
    hlambdapt->SetBinError  (ipt+1, templambdaerr);
    hnSig1spt->SetBinContent(ipt+1, tempnSig1s);
    hnSig1spt->SetBinError  (ipt+1, tempnSig1serr);
    hnSig2spt->SetBinContent(ipt+1, tempnSig2s);
    hnSig2spt->SetBinError  (ipt+1, tempnSig2serr);
    hnSig3spt->SetBinContent(ipt+1, tempnSig3s);
    hnSig3spt->SetBinError  (ipt+1, tempnSig3serr);
    hnBkgpt->SetBinContent(ipt+1, tempnBkg);
    hnBkgpt->SetBinError  (ipt+1, tempnBkgerr);

    delete ws;
    delete NomFile;

    OldFileName = Form("%snomfitresults_upsilon_%s.root",olddirectory.Data(),kineLabel.Data());
    cout << OldFileName << endl;
    OldFile = TFile::Open(OldFileName,"READ");
    ws = (RooWorkspace*)OldFile->Get("workspace");  

    //extract parameter values
    if (collId==kPPDATA) ywidth = 2*ywidth;
    tempalpha = ws->var("alpha1s_1")->getVal();  
    tempalphaerr = ws->var("alpha1s_1")->getError();
    tempf1s = ws->var("f1s")->getVal();  
    tempf1serr = ws->var("f1s")->getError();
    tempmass = ws->var("m_{#Upsilon(1S)}")->getVal();  
    tempmasserr = ws->var("m_{#Upsilon(1S)}")->getError();
    tempn1s = ws->var("n1s_1")->getVal();  
    tempn1serr = ws->var("n1s_1")->getError();
    tempsigma1s = ws->var("sigma1s_1")->getVal();  
    tempsigma1serr = ws->var("sigma1s_1")->getError();
    tempx1s = ws->var("x1s")->getVal();  
    tempx1serr = ws->var("x1s")->getError();
    if (ptLow<5) {
      temperrmu = ws->var("#mu")->getVal();  
      temperrmuerr = ws->var("#mu")->getError();
      temperrsigma = ws->var("#sigma")->getVal();  
      temperrsigmaerr = ws->var("#sigma")->getError();
    }
    else {
      temperrmu = 0;
      temperrmuerr = 0;
      temperrsigma = 0;  
      temperrsigmaerr = 0;
    }
    templambda = ws->var("#lambda")->getVal();  
    templambdaerr = ws->var("#lambda")->getError();
    tempnSig1s = ws->var("nSig1s")->getVal();//ywidth;  
    tempnSig1serr = ws->var("nSig1s")->getError();
    tempnSig2s = ws->var("nSig2s")->getVal();//ywidth;  
    tempnSig2serr = ws->var("nSig2s")->getError();
    tempnSig3s = ws->var("nSig3s")->getVal();//ywidth;  
    tempnSig3serr = ws->var("nSig3s")->getError();
    tempnBkg = ws->var("nBkg")->getVal();//ywidth;  
    tempnBkgerr = ws->var("nBkg")->getError();

    //fill histograms
    oldhalphapt->SetBinContent(ipt+1, tempalpha);
    oldhalphapt->SetBinError  (ipt+1, tempalphaerr);
    oldhf1spt->SetBinContent(ipt+1, tempf1s);
    oldhf1spt->SetBinError  (ipt+1, tempf1serr);
    oldhmasspt->SetBinContent(ipt+1, tempmass);
    oldhmasspt->SetBinError  (ipt+1, tempmasserr);
    oldhn1spt->SetBinContent(ipt+1, tempn1s);
    oldhn1spt->SetBinError  (ipt+1, tempn1serr);
    oldhsigma1spt->SetBinContent(ipt+1, tempsigma1s);
    oldhsigma1spt->SetBinError  (ipt+1, tempsigma1serr);
    oldhx1spt->SetBinContent(ipt+1, tempx1s);
    oldhx1spt->SetBinError  (ipt+1, tempx1serr);
    if (ptLow<5) {
      oldherrmupt->SetBinContent(ipt+1, temperrmu);
      oldherrmupt->SetBinError  (ipt+1, temperrmuerr);
      oldherrsigmapt->SetBinContent(ipt+1, temperrsigma);
      oldherrsigmapt->SetBinError  (ipt+1, temperrsigmaerr);
    }
    oldhlambdapt->SetBinContent(ipt+1, templambda);
    oldhlambdapt->SetBinError  (ipt+1, templambdaerr);
    oldhnSig1spt->SetBinContent(ipt+1, tempnSig1s);
    oldhnSig1spt->SetBinError  (ipt+1, tempnSig1serr);
    oldhnSig2spt->SetBinContent(ipt+1, tempnSig2s);
    oldhnSig2spt->SetBinError  (ipt+1, tempnSig2serr);
    oldhnSig3spt->SetBinContent(ipt+1, tempnSig3s);
    oldhnSig3spt->SetBinError  (ipt+1, tempnSig3serr);
    oldhnBkgpt->SetBinContent(ipt+1, tempnBkg);
    oldhnBkgpt->SetBinError  (ipt+1, tempnBkgerr);

    delete ws;
    delete OldFile;
  }

  //draw pt plots
  calpha->cd();
  halphapt->SetXTitle("pT");
  halphapt->SetYTitle(paramNames[2]);
  halphapt->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphapt->SetLineColor(kGreen);
  halphapt->SetMarkerColor(kGreen);
  halphapt->Draw();
  oldhalphapt->SetLineColor(kBlue);
  oldhalphapt->SetMarkerColor(kBlue);
  oldhalphapt->Draw("same");
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],30,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],30,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();
  TLegend* lalphapt = new TLegend(0.6,0.7,0.9,0.9);
  lalphapt->AddEntry(oldhalphapt,"Old","l");
  lalphapt->AddEntry(halphapt,"New","l");
  lalphapt->Draw("same");

  cf1s->cd();
  hf1spt->SetXTitle("pT");
  hf1spt->SetYTitle(paramNames[4]);
  hf1spt->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1spt->SetLineColor(kGreen);
  hf1spt->SetMarkerColor(kGreen);
  hf1spt->Draw();
  oldhf1spt->SetLineColor(kBlue);
  oldhf1spt->SetMarkerColor(kBlue);
  oldhf1spt->Draw("same");
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],30,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],30,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();
  TLegend* lf1spt = new TLegend(0.6,0.7,0.9,0.9);
  lf1spt->AddEntry(oldhf1spt,"Old","l");
  lf1spt->AddEntry(hf1spt,"New","l");
  lf1spt->Draw("same");

  cmass->cd();
  hmasspt->SetXTitle("pT");
  hmasspt->SetYTitle(paramNames[8]);
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->SetLineColor(kGreen);
  hmasspt->SetMarkerColor(kGreen);
  hmasspt->Draw();
  oldhmasspt->SetLineColor(kBlue);
  oldhmasspt->SetMarkerColor(kBlue);
  oldhmasspt->Draw("same");
  //hmasspt->Fit("pol1");
  TLine *massptlinepdg = new TLine(0,9.46,30,9.46);
  massptlinepdg->SetLineColor(3);
  massptlinepdg->Draw();
  TLine *massptlineLow = new TLine(0,9.36,30,9.36);
  massptlineLow->SetLineColor(kRed);
  massptlineLow->Draw();
  TLine *massptlineHigh = new TLine(0,9.56,30,9.56);
  massptlineHigh->SetLineColor(kRed);
  massptlineHigh->Draw();
  TLegend* lmasspt = new TLegend(0.6,0.7,0.9,0.9);
  lmasspt->AddEntry(oldhmasspt,"Old","l");
  lmasspt->AddEntry(hmasspt,"New","l");
  lmasspt->Draw("same");

  cn1s->cd();
  hn1spt->SetXTitle("pT");
  hn1spt->SetYTitle(paramNames[3]);
  hn1spt->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1spt->SetLineColor(kGreen);
  hn1spt->SetMarkerColor(kGreen);
  hn1spt->Draw();
  oldhn1spt->SetLineColor(kBlue);
  oldhn1spt->SetMarkerColor(kBlue);
  oldhn1spt->Draw("same");
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],30,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],30,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();
  TLegend* ln1spt = new TLegend(0.6,0.7,0.9,0.9);
  ln1spt->AddEntry(oldhn1spt,"Old","l");
  ln1spt->AddEntry(hn1spt,"New","l");
  ln1spt->Draw("same");

  csigma1s->cd();
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->SetYTitle(paramNames[0]);
  hsigma1spt->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1spt->SetLineColor(kGreen);
  hsigma1spt->SetMarkerColor(kGreen);
  hsigma1spt->Draw();
  oldhsigma1spt->SetLineColor(kBlue);
  oldhsigma1spt->SetMarkerColor(kBlue);
  oldhsigma1spt->Draw("same");
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],30,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],30,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();
  TLegend* lsigmapt = new TLegend(0.6,0.7,0.9,0.9);
  lsigmapt->AddEntry(oldhsigma1spt,"Old","l");
  lsigmapt->AddEntry(hsigma1spt,"New","l");
  lsigmapt->Draw("same");

  cx1s->cd();
  hx1spt->SetXTitle("pT");
  hx1spt->SetYTitle(paramNames[1]);
  hx1spt->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1spt->SetLineColor(kGreen);
  hx1spt->SetMarkerColor(kGreen);
  hx1spt->Draw();
  oldhx1spt->SetLineColor(kBlue);
  oldhx1spt->SetMarkerColor(kBlue);
  oldhx1spt->Draw("same");
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],30,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],30,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();
  TLegend* lx1spt = new TLegend(0.6,0.7,0.9,0.9);
  lx1spt->AddEntry(oldhx1spt,"Old","l");
  lx1spt->AddEntry(hx1spt,"New","l");
  lx1spt->Draw("same");

  cerrmu->cd();
  herrmupt->SetXTitle("pT");
  herrmupt->SetYTitle(paramNames[6]);
  herrmupt->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmupt->SetLineColor(kGreen);
  herrmupt->SetMarkerColor(kGreen);
  herrmupt->Draw();
  oldherrmupt->SetLineColor(kBlue);
  oldherrmupt->SetMarkerColor(kBlue);
  oldherrmupt->Draw("same");
  TLine *errmuptlineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  errmuptlineLow->SetLineColor(kRed);
  errmuptlineLow->Draw();
  TLine *errmuptlineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  errmuptlineHigh->SetLineColor(kRed);
  errmuptlineHigh->Draw();
  TLegend* lerrmupt = new TLegend(0.6,0.7,0.9,0.9);
  lerrmupt->AddEntry(oldherrmupt,"Old","l");
  lerrmupt->AddEntry(herrmupt,"New","l");
  lerrmupt->Draw("same");

  cerrsigma->cd();
  herrsigmapt->SetXTitle("pT");
  herrsigmapt->SetYTitle(paramNames[5]);
  herrsigmapt->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmapt->SetLineColor(kGreen);
  herrsigmapt->SetMarkerColor(kGreen);
  herrsigmapt->Draw();
  oldherrsigmapt->SetLineColor(kBlue);
  oldherrsigmapt->SetMarkerColor(kBlue);
  oldherrsigmapt->Draw("same");
  TLine *errsigmaptlineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  errsigmaptlineLow->SetLineColor(kRed);
  errsigmaptlineLow->Draw();
  TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  errsigmaptlineHigh->SetLineColor(kRed);
  errsigmaptlineHigh->Draw();
  TLegend* lerrsigmapt = new TLegend(0.6,0.7,0.9,0.9);
  lerrsigmapt->AddEntry(oldherrsigmapt,"Old","l");
  lerrsigmapt->AddEntry(herrsigmapt,"New","l");
  lerrsigmapt->Draw("same");

  clambda->cd();
  hlambdapt->SetXTitle("pT");
  hlambdapt->SetYTitle(paramNames[7]);
  hlambdapt->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdapt->SetLineColor(kGreen);
  hlambdapt->SetMarkerColor(kGreen);
  hlambdapt->Draw();
  oldhlambdapt->SetLineColor(kBlue);
  oldhlambdapt->SetMarkerColor(kBlue);
  oldhlambdapt->Draw("same");
  TLine *lambdaptlineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  lambdaptlineLow->SetLineColor(kRed);
  lambdaptlineLow->Draw();
  TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  lambdaptlineHigh->SetLineColor(kRed);
  lambdaptlineHigh->Draw();
  TLegend* llambdapt = new TLegend(0.6,0.7,0.9,0.9);
  llambdapt->AddEntry(oldhlambdapt,"Old","l");
  llambdapt->AddEntry(hlambdapt,"New","l");
  llambdapt->Draw("same");

  cnSig1s->cd();
  //gPad->SetLogy();
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->SetYTitle(paramNames[9]);
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->SetLineColor(kGreen);
  hnSig1spt->SetMarkerColor(kGreen);
  hnSig1spt->Draw();
  oldhnSig1spt->SetLineColor(kBlue);
  oldhnSig1spt->SetMarkerColor(kBlue);
  oldhnSig1spt->Draw("same");
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,30,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();
  TLegend* lnSig1spt = new TLegend(0.6,0.7,0.9,0.9);
  lnSig1spt->AddEntry(oldhnSig1spt,"Old","l");
  lnSig1spt->AddEntry(hnSig1spt,"New","l");
  lnSig1spt->Draw("same");

  cnSig2s->cd();
  //gPad->SetLogy();
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->SetYTitle(paramNames[10]);
  hnSig2spt->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2spt->SetLineColor(kGreen);
  hnSig2spt->SetMarkerColor(kGreen);
  hnSig2spt->Draw();
  oldhnSig2spt->SetLineColor(kBlue);
  oldhnSig2spt->SetMarkerColor(kBlue);
  oldhnSig2spt->Draw("same");
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,-20,30,-20);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();
  TLegend* lnSig2spt = new TLegend(0.6,0.7,0.9,0.9);
  lnSig2spt->AddEntry(oldhnSig2spt,"Old","l");
  lnSig2spt->AddEntry(hnSig2spt,"New","l");
  lnSig2spt->Draw("same");

  cnSig3s->cd();
  //gPad->SetLogy();
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->SetYTitle(paramNames[11]);
  hnSig3spt->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3spt->SetLineColor(kGreen);
  hnSig3spt->SetMarkerColor(kGreen);
  hnSig3spt->Draw();
  oldhnSig3spt->SetLineColor(kBlue);
  oldhnSig3spt->SetMarkerColor(kBlue);
  oldhnSig3spt->Draw("same");
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,-50,30,-50);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();
  TLegend* lnSig3spt = new TLegend(0.6,0.7,0.9,0.9);
  lnSig3spt->AddEntry(oldhnSig3spt,"Old","l");
  lnSig3spt->AddEntry(hnSig3spt,"New","l");
  lnSig3spt->Draw("same");

  cnBkg->cd();
  //gPad->SetLogy();
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->SetYTitle(paramNames[12]);
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgpt->SetLineColor(kGreen);
  hnBkgpt->SetMarkerColor(kGreen);
  hnBkgpt->Draw();
  oldhnBkgpt->SetLineColor(kBlue);
  oldhnBkgpt->SetMarkerColor(kBlue);
  oldhnBkgpt->Draw("same");
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,30,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();
  TLegend* lnBkgpt = new TLegend(0.6,0.7,0.9,0.9);
  lnBkgpt->AddEntry(oldhnBkgpt,"Old","l");
  lnBkgpt->AddEntry(hnBkgpt,"New","l");
  lnBkgpt->Draw("same");

  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    ptLow = 0;
    ptHigh = 30;
    yLow = *(ybinsptr+iy);
    yHigh = *(ybinsptr+iy+1);
    /*if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(ybins[iy+1]);
      yHigh = TMath::Abs(ybins[iy]);
    }*/

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",newdirectory.Data(),kineLabel.Data());
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

    delete yws;
    delete NomFile;

    OldFileName = Form("%snomfitresults_upsilon_%s.root",olddirectory.Data(),kineLabel.Data());
    cout << OldFileName << endl;
    OldFile = TFile::Open(OldFileName,"READ");
    yws = (RooWorkspace*)OldFile->Get("workspace");  

    //extract parameter values
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
    oldhalphay->SetBinContent(iy+1, tempalpha);
    oldhalphay->SetBinError  (iy+1, tempalphaerr);
    oldhf1sy->SetBinContent(iy+1, tempf1s);
    oldhf1sy->SetBinError  (iy+1, tempf1serr);
    oldhmassy->SetBinContent(iy+1, tempmass);
    oldhmassy->SetBinError  (iy+1, tempmasserr);
    oldhn1sy->SetBinContent(iy+1, tempn1s);
    oldhn1sy->SetBinError  (iy+1, tempn1serr);
    oldhsigma1sy->SetBinContent(iy+1, tempsigma1s);
    oldhsigma1sy->SetBinError  (iy+1, tempsigma1serr);
    oldhx1sy->SetBinContent(iy+1, tempx1s);
    oldhx1sy->SetBinError  (iy+1, tempx1serr);
    if (ptLow<5) {
      oldherrmuy->SetBinContent(iy+1, temperrmu);
      oldherrmuy->SetBinError  (iy+1, temperrmuerr);
      oldherrsigmay->SetBinContent(iy+1, temperrsigma);
      oldherrsigmay->SetBinError  (iy+1, temperrsigmaerr);
    }
    oldhlambday->SetBinContent(iy+1, templambda);
    oldhlambday->SetBinError  (iy+1, templambdaerr);
    oldhnSig1sy->SetBinContent(iy+1, tempnSig1s);
    oldhnSig1sy->SetBinError  (iy+1, tempnSig1serr);
    oldhnSig2sy->SetBinContent(iy+1, tempnSig2s);
    oldhnSig2sy->SetBinError  (iy+1, tempnSig2serr);
    oldhnSig3sy->SetBinContent(iy+1, tempnSig3s);
    oldhnSig3sy->SetBinError  (iy+1, tempnSig3serr);
    oldhnBkgy->SetBinContent(iy+1, tempnBkg);
    oldhnBkgy->SetBinError  (iy+1, tempnBkgerr);

    delete yws;
    delete OldFile;
  }

  //draw rapidity plots
  calpha2->cd();
  halphay->SetXTitle(yTitle);
  halphay->SetYTitle(paramNames[2]);
  halphay->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphay->SetLineColor(kGreen);
  halphay->SetMarkerColor(kGreen);
  halphay->Draw();
  oldhalphay->SetLineColor(kBlue);
  oldhalphay->SetMarkerColor(kBlue);
  oldhalphay->Draw("same");
  //halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(ymin,paramslower[2],ymax,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(ymin,paramsupper[2],ymax,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();
  TLegend* lalphay = new TLegend(0.6,0.7,0.9,0.9);
  lalphay->AddEntry(oldhalphay,"Old","l");
  lalphay->AddEntry(halphay,"New","l");
  lalphay->Draw("same");

  cf1s2->cd();
  hf1sy->SetXTitle(yTitle);
  hf1sy->SetYTitle(paramNames[4]);
  hf1sy->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1sy->SetLineColor(kGreen);
  hf1sy->SetMarkerColor(kGreen);
  hf1sy->Draw();
  oldhf1sy->SetLineColor(kBlue);
  oldhf1sy->SetMarkerColor(kBlue);
  oldhf1sy->Draw("same");
  //hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(ymin,paramslower[4],ymax,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(ymin,paramsupper[4],ymax,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();
  TLegend* lf1sy = new TLegend(0.6,0.7,0.9,0.9);
  lf1sy->AddEntry(oldhf1sy,"Old","l");
  lf1sy->AddEntry(hf1sy,"New","l");
  lf1sy->Draw("same");

  cmass2->cd();
  hmassy->SetXTitle(yTitle);
  hmassy->SetYTitle(paramNames[8]);
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->SetLineColor(kGreen);
  hmassy->SetMarkerColor(kGreen);
  hmassy->Draw();
  oldhmassy->SetLineColor(kBlue);
  oldhmassy->SetMarkerColor(kBlue);
  oldhmassy->Draw("same");
  //hmassy->Fit("pol1");
  TLine *massylinepdg = new TLine(ymin,9.46,ymax,9.46);
  massylinepdg->SetLineColor(3);
  massylinepdg->Draw();
  TLine *massylineLow = new TLine(ymin,9.36,ymax,9.36);
  massylineLow->SetLineColor(kRed);
  massylineLow->Draw();
  TLine *massylineHigh = new TLine(ymin,9.56,ymax,9.56);
  massylineHigh->SetLineColor(kRed);
  massylineHigh->Draw();
  TLegend* lmassy = new TLegend(0.6,0.7,0.9,0.9);
  lmassy->AddEntry(oldhmassy,"Old","l");
  lmassy->AddEntry(hmassy,"New","l");
  lmassy->Draw("same");

  cn1s2->cd();
  hn1sy->SetXTitle(yTitle);
  hn1sy->SetYTitle(paramNames[3]);
  hn1sy->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1sy->SetLineColor(kGreen);
  hn1sy->SetMarkerColor(kGreen);
  hn1sy->Draw();
  oldhn1sy->SetLineColor(kBlue);
  oldhn1sy->SetMarkerColor(kBlue);
  oldhn1sy->Draw("same");
  //hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(ymin,paramslower[3],ymax,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(ymin,paramsupper[3],ymax,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();
  TLegend* ln1sy = new TLegend(0.6,0.7,0.9,0.9);
  ln1sy->AddEntry(oldhn1sy,"Old","l");
  ln1sy->AddEntry(hn1sy,"New","l");
  ln1sy->Draw("same");

  csigma1s2->cd();
  hsigma1sy->SetXTitle(yTitle);
  hsigma1sy->SetYTitle(paramNames[0]);
  hsigma1sy->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1sy->SetLineColor(kGreen);
  hsigma1sy->SetMarkerColor(kGreen);
  hsigma1sy->Draw();
  oldhsigma1sy->SetLineColor(kBlue);
  oldhsigma1sy->SetMarkerColor(kBlue);
  oldhsigma1sy->Draw("same");
  //hsigma1sy->Fit("pol1");
  TLine *sigmaylineLow = new TLine(ymin,paramslower[0],ymax,paramslower[0]);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(ymin,paramsupper[0],ymax,paramsupper[0]);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();
  TLegend* lsigmay = new TLegend(0.6,0.7,0.9,0.9);
  lsigmay->AddEntry(oldhsigma1sy,"Old","l");
  lsigmay->AddEntry(hsigma1sy,"New","l");
  lsigmay->Draw("same");

  cx1s2->cd();
  hx1sy->SetXTitle(yTitle);
  hx1sy->SetYTitle(paramNames[1]);
  hx1sy->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1sy->SetLineColor(kGreen);
  hx1sy->SetMarkerColor(kGreen);
  hx1sy->Draw();
  oldhx1sy->SetLineColor(kBlue);
  oldhx1sy->SetMarkerColor(kBlue);
  oldhx1sy->Draw("same");
  //hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(ymin,paramslower[1],ymax,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(ymin,paramsupper[1],ymax,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();
  TLegend* lx1sy = new TLegend(0.6,0.7,0.9,0.9);
  lx1sy->AddEntry(oldhx1sy,"Old","l");
  lx1sy->AddEntry(hx1sy,"New","l");
  lx1sy->Draw("same");

  cerrmu2->cd();
  herrmuy->SetXTitle(yTitle);
  herrmuy->SetYTitle(paramNames[6]);
  herrmuy->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmuy->SetLineColor(kGreen);
  herrmuy->SetMarkerColor(kGreen);
  herrmuy->Draw();
  oldherrmuy->SetLineColor(kBlue);
  oldherrmuy->SetMarkerColor(kBlue);
  oldherrmuy->Draw("same");
  TLine *errmuylineLow = new TLine(ymin,paramslower[6],ymax,paramslower[6]);
  errmuylineLow->SetLineColor(kRed);
  errmuylineLow->Draw();
  TLine *errmuylineHigh = new TLine(ymin,paramsupper[6],ymax,paramsupper[6]);
  errmuylineHigh->SetLineColor(kRed);
  errmuylineHigh->Draw();
  TLegend* lerrmuy = new TLegend(0.6,0.7,0.9,0.9);
  lerrmuy->AddEntry(oldherrmuy,"Old","l");
  lerrmuy->AddEntry(herrmuy,"New","l");
  lerrmuy->Draw("same");

  cerrsigma2->cd();
  herrsigmay->SetXTitle(yTitle);
  herrsigmay->SetYTitle(paramNames[5]);
  herrsigmay->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmay->SetLineColor(kGreen);
  herrsigmay->SetMarkerColor(kGreen);
  herrsigmay->Draw();
  oldherrsigmay->SetLineColor(kBlue);
  oldherrsigmay->SetMarkerColor(kBlue);
  oldherrsigmay->Draw("same");
  TLine *errsigmaylineLow = new TLine(ymin,paramslower[5],ymax,paramslower[5]);
  errsigmaylineLow->SetLineColor(kRed);
  errsigmaylineLow->Draw();
  TLine *errsigmaylineHigh = new TLine(ymin,paramsupper[5],ymax,paramsupper[5]);
  errsigmaylineHigh->SetLineColor(kRed);
  errsigmaylineHigh->Draw();
  TLegend* lerrsigmay = new TLegend(0.6,0.7,0.9,0.9);
  lerrsigmay->AddEntry(oldherrsigmay,"Old","l");
  lerrsigmay->AddEntry(herrsigmay,"New","l");
  lerrsigmay->Draw("same");

  clambda2->cd();
  hlambday->SetXTitle(yTitle);
  hlambday->SetYTitle(paramNames[7]);
  hlambday->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambday->SetLineColor(kGreen);
  hlambday->SetMarkerColor(kGreen);
  hlambday->Draw();
  oldhlambday->SetLineColor(kBlue);
  oldhlambday->SetMarkerColor(kBlue);
  oldhlambday->Draw("same");
  TLine *lambdaylineLow = new TLine(ymin,paramslower[7],ymax,paramslower[7]);
  lambdaylineLow->SetLineColor(kRed);
  lambdaylineLow->Draw();
  TLine *lambdaylineHigh = new TLine(ymin,paramsupper[7],ymax,paramsupper[7]);
  lambdaylineHigh->SetLineColor(kRed);
  lambdaylineHigh->Draw();
  TLegend* llambday = new TLegend(0.6,0.7,0.9,0.9);
  llambday->AddEntry(oldhlambday,"Old","l");
  llambday->AddEntry(hlambday,"New","l");
  llambday->Draw("same");

  cnSig1s2->cd();
  //gPad->SetLogy();
  hnSig1sy->SetXTitle(yTitle);
  hnSig1sy->SetYTitle(paramNames[9]);
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sy->SetLineColor(kGreen);
  hnSig1sy->SetMarkerColor(kGreen);
  hnSig1sy->Draw();
  oldhnSig1sy->SetLineColor(kBlue);
  oldhnSig1sy->SetMarkerColor(kBlue);
  oldhnSig1sy->Draw("same");
  //hnSig1sy->Fit("pol1");
  TLine *nSig1sylineLow = new TLine(ymin,0,ymax,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();
  TLegend* lnSig1sy = new TLegend(0.6,0.7,0.9,0.9);
  lnSig1sy->AddEntry(oldhnSig1sy,"Old","l");
  lnSig1sy->AddEntry(hnSig1sy,"New","l");
  lnSig1sy->Draw("same");

  cnSig2s2->cd();
  //gPad->SetLogy();
  hnSig2sy->SetXTitle(yTitle);
  hnSig2sy->SetYTitle(paramNames[10]);
  hnSig2sy->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2sy->SetLineColor(kGreen);
  hnSig2sy->SetMarkerColor(kGreen);
  hnSig2sy->Draw();
  oldhnSig2sy->SetLineColor(kBlue);
  oldhnSig2sy->SetMarkerColor(kBlue);
  oldhnSig2sy->Draw("same");
  //hnSig2sy->Fit("pol1");
  TLine *nSig2sylineLow = new TLine(ymin,-20,ymax,-20);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();
  TLegend* lnSig2sy = new TLegend(0.6,0.7,0.9,0.9);
  lnSig2sy->AddEntry(oldhnSig2sy,"Old","l");
  lnSig2sy->AddEntry(hnSig2sy,"New","l");
  lnSig2sy->Draw("same");

  cnSig3s2->cd();
  //gPad->SetLogy();
  hnSig3sy->SetXTitle(yTitle);
  hnSig3sy->SetYTitle(paramNames[11]);
  hnSig3sy->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3sy->SetLineColor(kGreen);
  hnSig3sy->SetMarkerColor(kGreen);
  hnSig3sy->Draw();
  oldhnSig3sy->SetLineColor(kBlue);
  oldhnSig3sy->SetMarkerColor(kBlue);
  oldhnSig3sy->Draw("same");
  //hnSig3sy->Fit("pol1");
  TLine *nSig3sylineLow = new TLine(ymin,-50,ymax,-50);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();
  TLegend* lnSig3sy = new TLegend(0.6,0.7,0.9,0.9);
  lnSig3sy->AddEntry(oldhnSig3sy,"Old","l");
  lnSig3sy->AddEntry(hnSig3sy,"New","l");
  lnSig3sy->Draw("same");

  cnBkg2->cd();
  //gPad->SetLogy();
  hnBkgy->SetXTitle(yTitle);
  hnBkgy->SetYTitle(paramNames[12]);
  hnBkgy->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgy->SetLineColor(kGreen);
  hnBkgy->SetMarkerColor(kGreen);
  hnBkgy->Draw();
  oldhnBkgy->SetLineColor(kBlue);
  oldhnBkgy->SetMarkerColor(kBlue);
  oldhnBkgy->Draw("same");
  //hnBkgy->Fit("pol1");
  TLine *nBkgylineLow = new TLine(ymin,0,ymax,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();
  TLegend* lnBkgy = new TLegend(0.6,0.7,0.9,0.9);
  lnBkgy->AddEntry(oldhnBkgy,"Old","l");
  lnBkgy->AddEntry(hnBkgy,"New","l");
  lnBkgy->Draw("same");

  TString strId;
  if (collId==kPADATA) strId = "ComparePA";
  else if (collId==kPPDATA) strId = "ComparePP";

  //save plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins.png",strId.Data(),whichUpsilon));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins.png",strId.Data(),whichUpsilon));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins.png",strId.Data(),whichUpsilon));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins.png",strId.Data(),whichUpsilon));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins.png",strId.Data(),whichUpsilon));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins.png",strId.Data(),whichUpsilon));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins.png",strId.Data(),whichUpsilon));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins.png",strId.Data(),whichUpsilon));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins.png",strId.Data(),whichUpsilon));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins.png",strId.Data(),whichUpsilon));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins.png",strId.Data(),whichUpsilon));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins.png",strId.Data(),whichUpsilon));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins.png",strId.Data(),whichUpsilon));

  calpha2->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins2.png",strId.Data(),whichUpsilon));
  cf1s2->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins2.png",strId.Data(),whichUpsilon));
  cmass2->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins2.png",strId.Data(),whichUpsilon));
  cn1s2->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins2.png",strId.Data(),whichUpsilon));
  csigma1s2->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins2.png",strId.Data(),whichUpsilon));
  cx1s2->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins2.png",strId.Data(),whichUpsilon));
  cerrmu2->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins2.png",strId.Data(),whichUpsilon));
  cerrsigma2->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins2.png",strId.Data(),whichUpsilon));
  clambda2->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins2.png",strId.Data(),whichUpsilon));
  cnSig1s2->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins2.png",strId.Data(),whichUpsilon));
  cnSig2s2->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins2.png",strId.Data(),whichUpsilon));
  cnSig3s2->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins2.png",strId.Data(),whichUpsilon));
  cnBkg2->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins2.png",strId.Data(),whichUpsilon));

  //save pdf plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins.pdf",strId.Data(),whichUpsilon));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins.pdf",strId.Data(),whichUpsilon));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins.pdf",strId.Data(),whichUpsilon));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins.pdf",strId.Data(),whichUpsilon));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins.pdf",strId.Data(),whichUpsilon));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins.pdf",strId.Data(),whichUpsilon));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins.pdf",strId.Data(),whichUpsilon));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins.pdf",strId.Data(),whichUpsilon));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins.pdf",strId.Data(),whichUpsilon));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins.pdf",strId.Data(),whichUpsilon));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins.pdf",strId.Data(),whichUpsilon));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins.pdf",strId.Data(),whichUpsilon));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins.pdf",strId.Data(),whichUpsilon));

  calpha2->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins2.pdf",strId.Data(),whichUpsilon));
  cf1s2->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cmass2->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins2.pdf",strId.Data(),whichUpsilon));
  cn1s2->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins2.pdf",strId.Data(),whichUpsilon));
  csigma1s2->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cx1s2->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cerrmu2->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins2.pdf",strId.Data(),whichUpsilon));
  cerrsigma2->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins2.pdf",strId.Data(),whichUpsilon));
  clambda2->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins2.pdf",strId.Data(),whichUpsilon));
  cnSig1s2->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cnSig2s2->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cnSig3s2->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins2.pdf",strId.Data(),whichUpsilon));
  cnBkg2->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins2.pdf",strId.Data(),whichUpsilon));
/*
  calpha->Print("h1.pdf(","pdf");
  cf1s->Print("h1.pdf","pdf");
  cmass->Print("h1.pdf","pdf");
  cn1s->Print("h1.pdf","pdf");
  csigma1s->Print("h1.pdf","pdf");
  cx1s->Print("h1.pdf","pdf");
  cerrmu->Print("h1.pdf","pdf");
  cerrsigma->Print("h1.pdf","pdf");
  clambda->Print("h1.pdf","pdf");
  cnSig1s->Print("h1.pdf","pdf");
  cnSig2s->Print("h1.pdf","pdf");
  cnSig3s->Print("h1.pdf","pdf");
  cnBkg->Print("h1.pdf)","pdf");
*/

  delete calpha;
  delete cf1s;
  delete cmass;
  delete cn1s;
  delete csigma1s;
  delete cx1s;
  delete cerrmu;
  delete cerrsigma;
  delete clambda;
  delete cnSig1s;
  delete cnSig2s;
  delete cnSig3s;
  delete cnBkg;

  delete calpha2;
  delete cf1s2;
  delete cmass2;
  delete cn1s2;
  delete csigma1s2;
  delete cx1s2;
  delete cerrmu2;
  delete cerrsigma2;
  delete clambda2;
  delete cnSig1s2;
  delete cnSig2s2;
  delete cnSig3s2;
  delete cnBkg2;

  delete halphapt;
  delete hf1spt;
  delete hmasspt;
  delete hn1spt;
  delete hsigma1spt;
  delete hx1spt;
  delete herrmupt;
  delete herrsigmapt;
  delete hlambdapt;
  delete hnSig1spt;
  delete hnSig2spt;
  delete hnSig3spt;
  delete hnBkgpt;

  delete halphay;
  delete hf1sy;
  delete hmassy;
  delete hn1sy;
  delete hsigma1sy;
  delete hx1sy;
  delete herrmuy;
  delete herrsigmay;
  delete hlambday;
  delete hnSig1sy;
  delete hnSig2sy;
  delete hnSig3sy;
  delete hnBkgy;
}
