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


void PlotFittedParams_pt0to6to30SinglePlots(int whichUpsilon=1, int collId=kPADATA) {

  TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";
  float scale = whichUpsilon*0.6;
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
  double paramsupper[8] = {0.25, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  double plotlimupper[8] = {0.3, 1.2, 5.5, 5.5, 1.2, 16.0, 16.0, 26.0};
  double plotlimlower[8] = {0.0, -0.2, 0.5, 0.5, -0.2, -1.0, -1.0, -1.0};

  //choose a set of bins
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[3] = {-1.93,0.0,1.93};
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

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs y",numybins,ybinsptr);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs y",numybins,ybinsptr);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs y",numybins,ybinsptr);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs y",numybins,ybinsptr);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs y",numybins,ybinsptr);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs y",numybins,ybinsptr);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs y",numybins,ybinsptr);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs y",numybins,ybinsptr);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs y",numybins,ybinsptr);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs y",numybins,ybinsptr);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs y",numybins,ybinsptr);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs y",numybins,ybinsptr);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs y",numybins,ybinsptr);

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

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float ptLow, ptHigh, yLow, yHigh;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //1S pt loop
  for (int ipt = 0; ipt<numybins; ipt++) {

    ptLow = 0;
    ptHigh = 6;
    yLow = *(ybinsptr+ipt);
    yHigh = *(ybinsptr+ipt+1);
    /*if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(ybins[ipt+1]);
      yHigh = TMath::Abs(ybins[ipt]);
    }*/

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
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

    delete NomFile;
  }

  //draw pt plots
  calpha->cd();
  halphapt->SetXTitle(yTitle);
  halphapt->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphapt->Draw();
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(ymin,paramslower[2],ymax,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(ymin,paramsupper[2],ymax,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd();
  hf1spt->SetXTitle(yTitle);
  hf1spt->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1spt->Draw();
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(ymin,paramslower[4],ymax,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(ymin,paramsupper[4],ymax,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd();
  hmasspt->SetXTitle(yTitle);
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  //hmasspt->Fit("pol1");
  TLine *massptlinepdg = new TLine(ymin,9.46,ymax,9.46);
  massptlinepdg->SetLineColor(3);
  massptlinepdg->Draw();
  TLine *massptlineLow = new TLine(ymin,9.36,ymax,9.36);
  massptlineLow->SetLineColor(kRed);
  massptlineLow->Draw();
  TLine *massptlineHigh = new TLine(ymin,9.56,ymax,9.56);
  massptlineHigh->SetLineColor(kRed);
  massptlineHigh->Draw();

  cn1s->cd();
  hn1spt->SetXTitle(yTitle);
  hn1spt->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1spt->Draw();
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(ymin,paramslower[3],ymax,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(ymin,paramsupper[3],ymax,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd();
  hsigma1spt->SetXTitle(yTitle);
  hsigma1spt->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1spt->Draw();
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(ymin,paramslower[0],ymax,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(ymin,paramsupper[0],ymax,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd();
  hx1spt->SetXTitle(yTitle);
  hx1spt->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1spt->Draw();
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(ymin,paramslower[1],ymax,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(ymin,paramsupper[1],ymax,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd();
  herrmupt->SetXTitle(yTitle);
  herrmupt->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmupt->Draw();
  TLine *errmuptlineLow = new TLine(ymin,paramslower[6],ymax,paramslower[6]);
  errmuptlineLow->SetLineColor(kRed);
  errmuptlineLow->Draw();
  TLine *errmuptlineHigh = new TLine(ymin,paramsupper[6],ymax,paramsupper[6]);
  errmuptlineHigh->SetLineColor(kRed);
  errmuptlineHigh->Draw();

  cerrsigma->cd();
  herrsigmapt->SetXTitle(yTitle);
  herrsigmapt->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmapt->Draw();
  TLine *errsigmaptlineLow = new TLine(ymin,paramslower[5],ymax,paramslower[5]);
  errsigmaptlineLow->SetLineColor(kRed);
  errsigmaptlineLow->Draw();
  TLine *errsigmaptlineHigh = new TLine(ymin,paramsupper[5],ymax,paramsupper[5]);
  errsigmaptlineHigh->SetLineColor(kRed);
  errsigmaptlineHigh->Draw();

  clambda->cd();
  hlambdapt->SetXTitle(yTitle);
  hlambdapt->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdapt->Draw();
  TLine *lambdaptlineLow = new TLine(ymin,paramslower[7],ymax,paramslower[7]);
  lambdaptlineLow->SetLineColor(kRed);
  lambdaptlineLow->Draw();
  TLine *lambdaptlineHigh = new TLine(ymin,paramsupper[7],ymax,paramsupper[7]);
  lambdaptlineHigh->SetLineColor(kRed);
  lambdaptlineHigh->Draw();

  cnSig1s->cd();
  //gPad->SetLogy();
  hnSig1spt->SetXTitle(yTitle);
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(ymin,0,ymax,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd();
  //gPad->SetLogy();
  hnSig2spt->SetXTitle(yTitle);
  hnSig2spt->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2spt->Draw();
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(ymin,-20,ymax,-20);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd();
  //gPad->SetLogy();
  hnSig3spt->SetXTitle(yTitle);
  hnSig3spt->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3spt->Draw();
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(ymin,-50,ymax,-50);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd();
  //gPad->SetLogy();
  hnBkgpt->SetXTitle(yTitle);
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgpt->Draw();
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(ymin,0,ymax,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    ptLow = 6;
    ptHigh = 30;
    yLow = *(ybinsptr+iy);
    yHigh = *(ybinsptr+iy+1);
    /*if (collId==kPPDATA && yLow<0) {
      yLow = TMath::Abs(ybins[iy+1]);
      yHigh = TMath::Abs(ybins[iy]);
    }*/

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

    delete NomFile;
  }


  //draw rapidity plots
  calpha->cd();
  halphay->SetLineColor(3);
  halphay->SetMarkerColor(3);
  halphay->Draw("same");
  TLegend* alphaLegend = new TLegend(0.6,0.7,0.9,0.9);
  alphaLegend->SetTextSize(16);
  alphaLegend->SetTextFont(43);
  alphaLegend->AddEntry("halphapt","Low Pt","lep");
  alphaLegend->AddEntry("halphay","High Pt","lep");
  alphaLegend->Draw("same");

  cf1s->cd();
  hf1sy->SetLineColor(3);
  hf1sy->SetMarkerColor(3);
  hf1sy->Draw("same");
  TLegend* f1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  f1sLegend->SetTextSize(16);
  f1sLegend->SetTextFont(43);
  f1sLegend->AddEntry("hf1spt","Low Pt","lep");
  f1sLegend->AddEntry("hf1sy","High Pt","lep");
  f1sLegend->Draw("same");

  cmass->cd();
  hmassy->SetLineColor(3);
  hmassy->SetMarkerColor(3);
  hmassy->Draw("same");
  TLegend* nmassLegend = new TLegend(0.6,0.7,0.9,0.9);
  nmassLegend->SetTextSize(16);
  nmassLegend->SetTextFont(43);
  nmassLegend->AddEntry("hmasspt","Low Pt","lep");
  nmassLegend->AddEntry("hmassy","High Pt","lep");
  nmassLegend->Draw("same");

  cn1s->cd();
  hn1sy->SetLineColor(3);
  hn1sy->SetMarkerColor(3);
  hn1sy->Draw("same");
  TLegend* n1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  n1sLegend->SetTextSize(16);
  n1sLegend->SetTextFont(43);
  n1sLegend->AddEntry("hn1spt","Low Pt","lep");
  n1sLegend->AddEntry("hn1sy","High Pt","lep");
  n1sLegend->Draw("same");

  csigma1s->cd();
  hsigma1sy->SetLineColor(3);
  hsigma1sy->SetMarkerColor(3);
  hsigma1sy->Draw("same");
  TLegend* sigma1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  sigma1sLegend->SetTextSize(16);
  sigma1sLegend->SetTextFont(43);
  sigma1sLegend->AddEntry("hsigma1spt","Low Pt","lep");
  sigma1sLegend->AddEntry("hsigma1sy","High Pt","lep");
  sigma1sLegend->Draw("same");

  cx1s->cd();
  hx1sy->SetLineColor(3);
  hx1sy->SetMarkerColor(3);
  hx1sy->Draw("same");
  TLegend* x1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  x1sLegend->SetTextSize(16);
  x1sLegend->SetTextFont(43);
  x1sLegend->AddEntry("hx1spt","Low Pt","lep");
  x1sLegend->AddEntry("hx1sy","High Pt","lep");
  x1sLegend->Draw("same");

  cerrmu->cd();
  herrmuy->SetLineColor(3);
  herrmuy->SetMarkerColor(3);
  herrmuy->Draw("same");
  TLegend* errmuLegend = new TLegend(0.6,0.7,0.9,0.9);
  errmuLegend->SetTextSize(16);
  errmuLegend->SetTextFont(43);
  errmuLegend->AddEntry("herrmupt","Low Pt","lep");
  errmuLegend->AddEntry("herrmuy","High Pt","lep");
  errmuLegend->Draw("same");

  cerrsigma->cd();
  herrsigmay->SetLineColor(3);
  herrsigmay->SetMarkerColor(3);
  herrsigmay->Draw("same");
  TLegend* errsigmaLegend = new TLegend(0.6,0.7,0.9,0.9);
  errsigmaLegend->SetTextSize(16);
  errsigmaLegend->SetTextFont(43);
  errsigmaLegend->AddEntry("herrsigmapt","Low Pt","lep");
  errsigmaLegend->AddEntry("herrsigmay","High Pt","lep");
  errsigmaLegend->Draw("same");

  clambda->cd();
  hlambday->SetLineColor(3);
  hlambday->SetMarkerColor(3);
  hlambday->Draw("same");
  TLegend* lambdaLegend = new TLegend(0.6,0.7,0.9,0.9);
  lambdaLegend->SetTextSize(16);
  lambdaLegend->SetTextFont(43);
  lambdaLegend->AddEntry("hlambdapt","Low Pt","lep");
  lambdaLegend->AddEntry("hlambday","High Pt","lep");
  lambdaLegend->Draw("same");

  cnSig1s->cd();
  hnSig1sy->SetLineColor(3);
  hnSig1sy->SetMarkerColor(3);
  hnSig1sy->Draw("same");
  TLegend* nSig1sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig1sLegend->SetTextSize(16);
  nSig1sLegend->SetTextFont(43);
  nSig1sLegend->AddEntry("hnSig1spt","Low Pt","lep");
  nSig1sLegend->AddEntry("hnSig1sy","High Pt","lep");
  nSig1sLegend->Draw("same");

  cnSig2s->cd();
  hnSig2sy->SetLineColor(3);
  hnSig2sy->SetMarkerColor(3);
  hnSig2sy->Draw("same");
  TLegend* nSig2sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig2sLegend->SetTextSize(16);
  nSig2sLegend->SetTextFont(43);
  nSig2sLegend->AddEntry("hnSig2spt","Low Pt","lep");
  nSig2sLegend->AddEntry("hnSig2sy","High Pt","lep");
  nSig2sLegend->Draw("same");

  cnSig3s->cd();
  hnSig3sy->SetLineColor(3);
  hnSig3sy->SetMarkerColor(3);
  hnSig3sy->Draw("same");
  TLegend* nSig3sLegend = new TLegend(0.6,0.7,0.9,0.9);
  nSig3sLegend->SetTextSize(16);
  nSig3sLegend->SetTextFont(43);
  nSig3sLegend->AddEntry("hnSig3spt","Low Pt","lep");
  nSig3sLegend->AddEntry("hnSig3sy","High Pt","lep");
  nSig3sLegend->Draw("same");

  cnBkg->cd();
  hnBkgy->SetLineColor(3);
  hnBkgy->SetMarkerColor(3);
  hnBkgy->Draw("same");
  TLegend* nBkgLegend = new TLegend(0.6,0.7,0.9,0.9);
  nBkgLegend->SetTextSize(16);
  nBkgLegend->SetTextFont(43);
  nBkgLegend->AddEntry("hnBkgpt","Low Pt","lep");
  nBkgLegend->AddEntry("hnBkgy","High Pt","lep");
  nBkgLegend->Draw("same");

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

  //save plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins_pt0to6to30.png",strId.Data(),whichUpsilon));

  //save pdf plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_%isbins_pt0to6to30.pdf",strId.Data(),whichUpsilon));
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
