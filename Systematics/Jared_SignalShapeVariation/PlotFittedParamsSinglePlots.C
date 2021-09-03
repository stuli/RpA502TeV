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


void PlotFittedParamsSinglePlots(int whichUpsilon=1, int collId=kPADATA) {

  TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";
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
  double paramsupper[8] = {0.25, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
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

  TString kineLabel, NomFileName;
  TFile* NomFile;
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
  halphapt->SetXTitle("pT");
  halphapt->SetYTitle(paramNames[2]);
  halphapt->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphapt->Draw();
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],30,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],30,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd();
  hf1spt->SetXTitle("pT");
  hf1spt->SetYTitle(paramNames[4]);
  hf1spt->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1spt->Draw();
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],30,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],30,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd();
  hmasspt->SetXTitle("pT");
  hmasspt->SetYTitle(paramNames[8]);
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
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

  cn1s->cd();
  hn1spt->SetXTitle("pT");
  hn1spt->SetYTitle(paramNames[3]);
  hn1spt->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1spt->Draw();
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],30,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],30,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd();
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->SetYTitle(paramNames[0]);
  hsigma1spt->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1spt->Draw();
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],30,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],30,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd();
  hx1spt->SetXTitle("pT");
  hx1spt->SetYTitle(paramNames[1]);
  hx1spt->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1spt->Draw();
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],30,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],30,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd();
  herrmupt->SetXTitle("pT");
  herrmupt->SetYTitle(paramNames[6]);
  herrmupt->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmupt->Draw();
  TLine *errmuptlineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  errmuptlineLow->SetLineColor(kRed);
  errmuptlineLow->Draw();
  TLine *errmuptlineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  errmuptlineHigh->SetLineColor(kRed);
  errmuptlineHigh->Draw();

  cerrsigma->cd();
  herrsigmapt->SetXTitle("pT");
  herrsigmapt->SetYTitle(paramNames[5]);
  herrsigmapt->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmapt->Draw();
  TLine *errsigmaptlineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  errsigmaptlineLow->SetLineColor(kRed);
  errsigmaptlineLow->Draw();
  TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  errsigmaptlineHigh->SetLineColor(kRed);
  errsigmaptlineHigh->Draw();

  clambda->cd();
  hlambdapt->SetXTitle("pT");
  hlambdapt->SetYTitle(paramNames[7]);
  hlambdapt->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambdapt->Draw();
  TLine *lambdaptlineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  lambdaptlineLow->SetLineColor(kRed);
  lambdaptlineLow->Draw();
  TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  lambdaptlineHigh->SetLineColor(kRed);
  lambdaptlineHigh->Draw();

  cnSig1s->cd();
  //gPad->SetLogy();
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->SetYTitle(paramNames[9]);
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,30,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd();
  //gPad->SetLogy();
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->SetYTitle(paramNames[10]);
  hnSig2spt->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2spt->Draw();
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,-20,30,-20);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd();
  //gPad->SetLogy();
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->SetYTitle(paramNames[11]);
  hnSig3spt->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3spt->Draw();
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,-50,30,-50);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd();
  //gPad->SetLogy();
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->SetYTitle(paramNames[12]);
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgpt->Draw();
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,30,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


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
  calpha2->cd();
  halphay->SetXTitle(yTitle);
  halphay->SetYTitle(paramNames[2]);
  halphay->GetYaxis()->SetRangeUser(plotlimlower[2],plotlimupper[2]);
  halphay->Draw();
  //halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(ymin,paramslower[2],ymax,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(ymin,paramsupper[2],ymax,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  cf1s2->cd();
  hf1sy->SetXTitle(yTitle);
  hf1sy->SetYTitle(paramNames[4]);
  hf1sy->GetYaxis()->SetRangeUser(plotlimlower[4],plotlimupper[4]);
  hf1sy->Draw();
  //hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(ymin,paramslower[4],ymax,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(ymin,paramsupper[4],ymax,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  cmass2->cd();
  hmassy->SetXTitle(yTitle);
  hmassy->SetYTitle(paramNames[8]);
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
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

  cn1s2->cd();
  hn1sy->SetXTitle(yTitle);
  hn1sy->SetYTitle(paramNames[3]);
  hn1sy->GetYaxis()->SetRangeUser(plotlimlower[3],plotlimupper[3]);
  hn1sy->Draw();
  //hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(ymin,paramslower[3],ymax,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(ymin,paramsupper[3],ymax,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  csigma1s2->cd();
  hsigma1sy->SetXTitle(yTitle);
  hsigma1sy->SetYTitle(paramNames[0]);
  hsigma1sy->GetYaxis()->SetRangeUser(plotlimlower[0],plotlimupper[0]);
  hsigma1sy->Draw();
  //hsigma1sy->Fit("pol1");
  TLine *sigmaylineLow = new TLine(ymin,paramslower[0],ymax,paramslower[0]);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(ymin,paramsupper[0],ymax,paramsupper[0]);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();

  cx1s2->cd();
  hx1sy->SetXTitle(yTitle);
  hx1sy->SetYTitle(paramNames[1]);
  hx1sy->GetYaxis()->SetRangeUser(plotlimlower[1],plotlimupper[1]);
  hx1sy->Draw();
  //hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(ymin,paramslower[1],ymax,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(ymin,paramsupper[1],ymax,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

  cerrmu2->cd();
  herrmuy->SetXTitle(yTitle);
  herrmuy->SetYTitle(paramNames[6]);
  herrmuy->GetYaxis()->SetRangeUser(plotlimlower[6],plotlimupper[6]);
  herrmuy->Draw();
  TLine *errmuylineLow = new TLine(ymin,paramslower[6],ymax,paramslower[6]);
  errmuylineLow->SetLineColor(kRed);
  errmuylineLow->Draw();
  TLine *errmuylineHigh = new TLine(ymin,paramsupper[6],ymax,paramsupper[6]);
  errmuylineHigh->SetLineColor(kRed);
  errmuylineHigh->Draw();

  cerrsigma2->cd();
  herrsigmay->SetXTitle(yTitle);
  herrsigmay->SetYTitle(paramNames[5]);
  herrsigmay->GetYaxis()->SetRangeUser(plotlimlower[5],plotlimupper[5]);
  herrsigmay->Draw();
  TLine *errsigmaylineLow = new TLine(ymin,paramslower[5],ymax,paramslower[5]);
  errsigmaylineLow->SetLineColor(kRed);
  errsigmaylineLow->Draw();
  TLine *errsigmaylineHigh = new TLine(ymin,paramsupper[5],ymax,paramsupper[5]);
  errsigmaylineHigh->SetLineColor(kRed);
  errsigmaylineHigh->Draw();

  clambda2->cd();
  hlambday->SetXTitle(yTitle);
  hlambday->SetYTitle(paramNames[7]);
  hlambday->GetYaxis()->SetRangeUser(plotlimlower[7],plotlimupper[7]);
  hlambday->Draw();
  TLine *lambdaylineLow = new TLine(ymin,paramslower[7],ymax,paramslower[7]);
  lambdaylineLow->SetLineColor(kRed);
  lambdaylineLow->Draw();
  TLine *lambdaylineHigh = new TLine(ymin,paramsupper[7],ymax,paramsupper[7]);
  lambdaylineHigh->SetLineColor(kRed);
  lambdaylineHigh->Draw();

  cnSig1s2->cd();
  //gPad->SetLogy();
  hnSig1sy->SetXTitle(yTitle);
  hnSig1sy->SetYTitle(paramNames[9]);
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sy->Draw();
  //hnSig1sy->Fit("pol1");
  TLine *nSig1sylineLow = new TLine(ymin,0,ymax,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();

  cnSig2s2->cd();
  //gPad->SetLogy();
  hnSig2sy->SetXTitle(yTitle);
  hnSig2sy->SetYTitle(paramNames[10]);
  hnSig2sy->GetYaxis()->SetRangeUser(-20,600*scale);
  hnSig2sy->Draw();
  //hnSig2sy->Fit("pol1");
  TLine *nSig2sylineLow = new TLine(ymin,-20,ymax,-20);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();

  cnSig3s2->cd();
  //gPad->SetLogy();
  hnSig3sy->SetXTitle(yTitle);
  hnSig3sy->SetYTitle(paramNames[11]);
  hnSig3sy->GetYaxis()->SetRangeUser(-50,300*scale);
  hnSig3sy->Draw();
  //hnSig3sy->Fit("pol1");
  TLine *nSig3sylineLow = new TLine(ymin,-50,ymax,-50);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();

  cnBkg2->cd();
  //gPad->SetLogy();
  hnBkgy->SetXTitle(yTitle);
  hnBkgy->SetYTitle(paramNames[12]);
  hnBkgy->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgy->Draw();
  //hnBkgy->Fit("pol1");
  TLine *nBkgylineLow = new TLine(ymin,0,ymax,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

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
