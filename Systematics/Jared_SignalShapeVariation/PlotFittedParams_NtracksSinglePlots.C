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


void PlotFittedParams_NtracksSinglePlots(int whichUpsilon=1) {

  TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";
  float scale = 1.0;
  if (whichUpsilon==3) scale = 2.0;
  int collId=kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.25, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  double plotlimupper[8] = {0.3, 1.2, 5.5, 5.5, 1.2, 16.0, 16.0, 26.0};
  double plotlimlower[8] = {0.0, -0.2, 0.5, 0.5, -0.2, -1.0, -1.0, -1.0};

  //choose a set of bins
  float ybins[3] = {-1.93,0.0,1.93};
  float hfbins3[3] = {0,12,120};
  float ntbins3[3] = {0,40,400};
  float hfbins1[5] = {0,12,19,27,120};
  float ntbins1[5] = {0,40,62,88,400};

  float* hfbinsptr;
  float* ntbinsptr;
  int numhfbinstemp, numntbinstemp;
  if (whichUpsilon==3) {
    hfbinsptr = &hfbins3[0];
    ntbinsptr = &ntbins3[0];
    numhfbinstemp = sizeof(hfbins3)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins3)/sizeof(float)-1;
  }
  else {
    hfbinsptr = &hfbins1[0];
    ntbinsptr = &ntbins1[0];
    numhfbinstemp = sizeof(hfbins1)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins1)/sizeof(float)-1;
  }

  const int numybins = 2;
  const int numhfbins = numhfbinstemp;
  const int numntbins = numntbinstemp;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetMarkerStyle(7);
  gStyle->SetMarkerColor(4);

  float ptLow = 0;
  float ptHigh = 30;
  float yLow, yHigh;
  float hfLow = 0;
  float hfHigh = 120;
  int ntLow, ntHigh;

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

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs ntracks",numntbins,ntbinsptr);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs ntracks",numntbins,ntbinsptr);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs ntracks",numntbins,ntbinsptr);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs ntracks",numntbins,ntbinsptr);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs ntracks",numntbins,ntbinsptr);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs ntracks",numntbins,ntbinsptr);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs ntracks",numntbins,ntbinsptr);

  TH1F* halphay = new TH1F("halphay","alpha vs ntracks",numntbins,ntbinsptr);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hmassy = new TH1F("hmassy","mass vs ntracks",numntbins,ntbinsptr);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs ntracks",numntbins,ntbinsptr);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs ntracks",numntbins,ntbinsptr);
  TH1F* herrmuy = new TH1F("herrmuy","errmu vs ntracks",numntbins,ntbinsptr);
  TH1F* herrsigmay = new TH1F("herrsigmay","errsigma vs ntracks",numntbins,ntbinsptr);
  TH1F* hlambday = new TH1F("hlambday","errlambda vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs ntracks",numntbins,ntbinsptr);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs ntracks",numntbins,ntbinsptr);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  yLow = -1.93;
  yHigh = 0.0;

  //1S backward loop
  for (int ipt = 0; ipt<numntbins; ipt++) {

    ntLow = *(ntbinsptr+ipt);
    ntHigh = *(ntbinsptr+ipt+1);

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntLow, ntHigh );
    NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
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

  }

  //draw pt plots
  calpha->cd();
  halphapt->SetXTitle("Ntracks");
  halphapt->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphapt->Draw();
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],400,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],400,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd();
  hf1spt->SetXTitle("Ntracks");
  hf1spt->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1spt->Draw();
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],400,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],400,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd();
  hmasspt->SetXTitle("Ntracks");
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  //hmasspt->Fit("pol1");
  TLine *massptlinepdg = new TLine(0,9.46,400,9.46);
  massptlinepdg->SetLineColor(3);
  massptlinepdg->Draw();
  TLine *massptlineLow = new TLine(0,9.36,400,9.36);
  massptlineLow->SetLineColor(kRed);
  massptlineLow->Draw();
  TLine *massptlineHigh = new TLine(0,9.56,400,9.56);
  massptlineHigh->SetLineColor(kRed);
  massptlineHigh->Draw();

  cn1s->cd();
  hn1spt->SetXTitle("Ntracks");
  hn1spt->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1spt->Draw();
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],400,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],400,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd();
  hsigma1spt->SetXTitle("Ntracks");
  hsigma1spt->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1spt->Draw();
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],400,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],400,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd();
  hx1spt->SetXTitle("Ntracks");
  hx1spt->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1spt->Draw();
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],400,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],400,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd();
  herrmupt->SetXTitle("Ntracks");
  herrmupt->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmupt->Draw();
  TLine *errmuptlineLow = new TLine(0,paramslower[6],400,paramslower[6]);
  errmuptlineLow->SetLineColor(kRed);
  errmuptlineLow->Draw();
  TLine *errmuptlineHigh = new TLine(0,paramsupper[6],400,paramsupper[6]);
  errmuptlineHigh->SetLineColor(kRed);
  errmuptlineHigh->Draw();

  cerrsigma->cd();
  herrsigmapt->SetXTitle("Ntracks");
  herrsigmapt->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmapt->Draw();
  TLine *errsigmaptlineLow = new TLine(0,paramslower[5],400,paramslower[5]);
  errsigmaptlineLow->SetLineColor(kRed);
  errsigmaptlineLow->Draw();
  TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],400,paramsupper[5]);
  errsigmaptlineHigh->SetLineColor(kRed);
  errsigmaptlineHigh->Draw();

  clambda->cd();
  hlambdapt->SetXTitle("Ntracks");
  hlambdapt->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdapt->Draw();
  TLine *lambdaptlineLow = new TLine(0,paramslower[7],400,paramslower[7]);
  lambdaptlineLow->SetLineColor(kRed);
  lambdaptlineLow->Draw();
  TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],400,paramsupper[7]);
  lambdaptlineHigh->SetLineColor(kRed);
  lambdaptlineHigh->Draw();

  cnSig1s->cd();
  //gPad->SetLogy();
  hnSig1spt->SetXTitle("Ntracks");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,400,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd();
  //gPad->SetLogy();
  hnSig2spt->SetXTitle("Ntracks");
  hnSig2spt->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2spt->Draw();
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,-20,400,-20);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd();
  //gPad->SetLogy();
  hnSig3spt->SetXTitle("Ntracks");
  hnSig3spt->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3spt->Draw();
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,-50,400,-50);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd();
  //gPad->SetLogy();
  hnBkgpt->SetXTitle("Ntracks");
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgpt->Draw();
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,400,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  yLow = 0.0;
  yHigh = 1.93;

  //1S forward loop
  for (int iy = 0; iy<numhfbins; iy++) {

    ntLow = *(ntbinsptr+iy);
    ntHigh = *(ntbinsptr+iy+1);

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntLow, ntHigh );
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

  //draw rapidity plots
  calpha->cd();
  halphay->SetLineColor(3);
  halphay->SetLineColor(3);
  halphay->Draw("same");
  TLegend* alphaLegend = new TLegend(0.6,0.7,0.9,0.9);
  alphaLegend->SetTextSize(16);
  alphaLegend->SetTextFont(43);
  alphaLegend->AddEntry("halphapt","Backward","lep");
  alphaLegend->AddEntry("halphay","Forward","lep");
  alphaLegend->Draw("same");

  cf1s->cd();
  hf1sy->SetLineColor(3);
  hf1sy->SetLineColor(3);
  hf1sy->Draw("same");
  TLegend* f1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  f1sLegend->SetTextSize(16);
  f1sLegend->SetTextFont(43);
  f1sLegend->AddEntry("hf1spt","Backward","lep");
  f1sLegend->AddEntry("hf1sy","Forward","lep");
  f1sLegend->Draw("same");

  cmass->cd();
  hmassy->SetLineColor(3);
  hmassy->SetLineColor(3);
  hmassy->Draw("same");
  TLegend* nmassLegend = new TLegend(0.6,0.7,0.9,0.9);
  nmassLegend->SetTextSize(16);
  nmassLegend->SetTextFont(43);
  nmassLegend->AddEntry("hmasspt","Backward","lep");
  nmassLegend->AddEntry("hmassy","Forward","lep");
  nmassLegend->Draw("same");

  cn1s->cd();
  hn1sy->SetLineColor(3);
  hn1sy->SetLineColor(3);
  hn1sy->Draw("same");
  TLegend* n1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  n1sLegend->SetTextSize(16);
  n1sLegend->SetTextFont(43);
  n1sLegend->AddEntry("hn1spt","Backward","lep");
  n1sLegend->AddEntry("hn1sy","Forward","lep");
  n1sLegend->Draw("same");

  csigma1s->cd();
  hsigma1sy->SetLineColor(3);
  hsigma1sy->SetLineColor(3);
  hsigma1sy->Draw("same");
  TLegend* sigma1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  sigma1sLegend->SetTextSize(16);
  sigma1sLegend->SetTextFont(43);
  sigma1sLegend->AddEntry("hsigma1spt","Backward","lep");
  sigma1sLegend->AddEntry("hsigma1sy","Forward","lep");
  sigma1sLegend->Draw("same");

  cx1s->cd();
  hx1sy->SetLineColor(3);
  hx1sy->SetLineColor(3);
  hx1sy->Draw("same");
  TLegend* x1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  x1sLegend->SetTextSize(16);
  x1sLegend->SetTextFont(43);
  x1sLegend->AddEntry("hx1spt","Backward","lep");
  x1sLegend->AddEntry("hx1sy","Forward","lep");
  x1sLegend->Draw("same");

  cerrmu->cd();
  herrmuy->SetLineColor(3);
  herrmuy->SetLineColor(3);
  herrmuy->Draw("same");
  TLegend* errmuLegend = new TLegend(0.6,0.7,0.9,0.9);
  errmuLegend->SetTextSize(16);
  errmuLegend->SetTextFont(43);
  errmuLegend->AddEntry("herrmupt","Backward","lep");
  errmuLegend->AddEntry("herrmuy","Forward","lep");
  errmuLegend->Draw("same");

  cerrsigma->cd();
  herrsigmay->SetLineColor(3);
  herrsigmay->SetLineColor(3);
  herrsigmay->Draw("same");
  TLegend* errsigmaLegend = new TLegend(0.6,0.7,0.9,0.9);
  errsigmaLegend->SetTextSize(16);
  errsigmaLegend->SetTextFont(43);
  errsigmaLegend->AddEntry("herrsigmapt","Backward","lep");
  errsigmaLegend->AddEntry("herrsigmay","Forward","lep");
  errsigmaLegend->Draw("same");

  clambda->cd();
  hlambday->SetLineColor(3);
  hlambday->SetLineColor(3);
  hlambday->Draw("same");
  TLegend* lambdaLegend = new TLegend(0.6,0.7,0.9,0.9);
  lambdaLegend->SetTextSize(16);
  lambdaLegend->SetTextFont(43);
  lambdaLegend->AddEntry("hlambdapt","Backward","lep");
  lambdaLegend->AddEntry("hlambday","Forward","lep");
  lambdaLegend->Draw("same");

  cnSig1s->cd();
  hnSig1sy->SetLineColor(3);
  hnSig1sy->SetLineColor(3);
  hnSig1sy->Draw("same");
  TLegend* nSig1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig1sLegend->SetTextSize(16);
  nSig1sLegend->SetTextFont(43);
  nSig1sLegend->AddEntry("hnSig1spt","Backward","lep");
  nSig1sLegend->AddEntry("hnSig1sy","Forward","lep");
  nSig1sLegend->Draw("same");

  cnSig2s->cd();
  hnSig2sy->SetLineColor(3);
  hnSig2sy->SetLineColor(3);
  hnSig2sy->Draw("same");
  TLegend* nSig2sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig2sLegend->SetTextSize(16);
  nSig2sLegend->SetTextFont(43);
  nSig2sLegend->AddEntry("hnSig2spt","Backward","lep");
  nSig2sLegend->AddEntry("hnSig2sy","Forward","lep");
  nSig2sLegend->Draw("same");

  cnSig3s->cd();
  hnSig3sy->SetLineColor(3);
  hnSig3sy->SetLineColor(3);
  hnSig3sy->Draw("same");
  TLegend* nSig3sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig3sLegend->SetTextSize(16);
  nSig3sLegend->SetTextFont(43);
  nSig3sLegend->AddEntry("hnSig3spt","Backward","lep");
  nSig3sLegend->AddEntry("hnSig3sy","Forward","lep");
  nSig3sLegend->Draw("same");

  cnBkg->cd();
  hnBkgy->SetLineColor(3);
  hnBkgy->SetLineColor(3);
  hnBkgy->Draw("same");
  TLegend* nBkgLegend = new TLegend(0.6,0.7,0.9,0.9);
  nBkgLegend->SetTextSize(16);
  nBkgLegend->SetTextFont(43);
  nBkgLegend->AddEntry("hnBkgpt","Backward","lep");
  nBkgLegend->AddEntry("hnBkgy","Forward","lep");
  nBkgLegend->Draw("same");

  //save plots
  calpha->SaveAs(Form("ParameterPlots/PAfitted_alpha_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cf1s->SaveAs(Form("ParameterPlots/PAfitted_f1s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cmass->SaveAs(Form("ParameterPlots/PAfitted_mass_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cn1s->SaveAs(Form("ParameterPlots/PAfitted_n1s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  csigma1s->SaveAs(Form("ParameterPlots/PAfitted_sigma1s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cx1s->SaveAs(Form("ParameterPlots/PAfitted_x1s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cerrmu->SaveAs(Form("ParameterPlots/PAfitted_errmu_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cerrsigma->SaveAs(Form("ParameterPlots/PAfitted_errsigma_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  clambda->SaveAs(Form("ParameterPlots/PAfitted_lambda_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cnSig1s->SaveAs(Form("ParameterPlots/PAfitted_nSig1s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cnSig2s->SaveAs(Form("ParameterPlots/PAfitted_nSig2s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cnSig3s->SaveAs(Form("ParameterPlots/PAfitted_nSig3s_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));
  cnBkg->SaveAs(Form("ParameterPlots/PAfitted_nBkg_%isbins_Ntracks_y%.2f-%.2f.png",whichUpsilon,yLow,yHigh));

  //save pdf plots
  calpha->SaveAs(Form("ParameterPlots/PAfitted_alpha_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cf1s->SaveAs(Form("ParameterPlots/PAfitted_f1s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cmass->SaveAs(Form("ParameterPlots/PAfitted_mass_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cn1s->SaveAs(Form("ParameterPlots/PAfitted_n1s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  csigma1s->SaveAs(Form("ParameterPlots/PAfitted_sigma1s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cx1s->SaveAs(Form("ParameterPlots/PAfitted_x1s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cerrmu->SaveAs(Form("ParameterPlots/PAfitted_errmu_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cerrsigma->SaveAs(Form("ParameterPlots/PAfitted_errsigma_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  clambda->SaveAs(Form("ParameterPlots/PAfitted_lambda_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cnSig1s->SaveAs(Form("ParameterPlots/PAfitted_nSig1s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cnSig2s->SaveAs(Form("ParameterPlots/PAfitted_nSig2s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cnSig3s->SaveAs(Form("ParameterPlots/PAfitted_nSig3s_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
  cnBkg->SaveAs(Form("ParameterPlots/PAfitted_nBkg_%isbins_Ntracks_y%.2f-%.2f.pdf",whichUpsilon,yLow,yHigh));
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
}
