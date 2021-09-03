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


void PlotFittedParams_Ntracks(int whichUpsilon=1, int j=0) {

  float scale = whichUpsilon*0.3;
  int collId=kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.3, 3.0, 3.321, 5.0, 1.0, 25.0, 25.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  //choose a set of bins
  if (whichUpsilon==1) {
    float ybins[5] = {0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[3] = {0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[2] = {0.0,1.93};
  }
  float hfbins[5] = {0,15,22,30,120};
  float ntbins[5] = {0,40,65,90,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntbins = sizeof(ntbins)/sizeof(int)-1;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");

  float ptLow = 0;
  float ptHigh = 30;
  float yLow, yHigh;
  float hfLow = 0;
  float hfHigh = 400;
  int ntLow, ntHigh;

  //set up canvases
  TCanvas *calpha = new TCanvas("calpha","calpha",4,45,800,400);
  calpha->Divide(2,1);
  TCanvas *cf1s = new TCanvas("cf1s","cf1s",4,45,800,400);
  cf1s->Divide(2,1);
  TCanvas *cmass = new TCanvas("cmass","cmass",4,45,900,400);
  cmass->Divide(2,1);
  TCanvas *cn1s = new TCanvas("cn1s","cn1s",4,45,800,400);
  cn1s->Divide(2,1);
  TCanvas *csigma1s = new TCanvas("csigma1s","csigma1s",4,45,800,400);
  csigma1s->Divide(2,1);
  TCanvas *cx1s = new TCanvas("cx1s","cx1s",4,45,800,400);
  cx1s->Divide(2,1);
  TCanvas *cerrmu = new TCanvas("cerrmu","cerrmu",4,45,800,400);
  cerrmu->Divide(2,1);
  TCanvas *cerrsigma = new TCanvas("cerrsigma","cerrsigma",4,45,800,400);
  cerrsigma->Divide(2,1);
  TCanvas *clambda = new TCanvas("clambda","clambda",4,45,800,400);
  clambda->Divide(2,1);
  TCanvas *cnSig1s = new TCanvas("cnSig1s","cnSig1s",4,45,800,400);
  cnSig1s->Divide(2,1);
  TCanvas *cnSig2s = new TCanvas("cnSig2s","cnSig2s",4,45,800,400);
  cnSig2s->Divide(2,1);
  TCanvas *cnSig3s = new TCanvas("cnSig3s","cnSig3s",4,45,800,400);
  cnSig3s->Divide(2,1);
  TCanvas *cnBkg = new TCanvas("cnBkg","cnBkg",4,45,900,400);
  cnBkg->Divide(2,1);

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs ntracks",numntbins,ntbins);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs ntracks",numntbins,ntbins);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs ntracks",numntbins,ntbins);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs ntracks",numntbins,ntbins);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs ntracks",numntbins,ntbins);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs ntracks",numntbins,ntbins);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs ntracks",numntbins,ntbins);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs ntracks",numntbins,ntbins);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs ntracks",numntbins,ntbins);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs ntracks",numntbins,ntbins);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs ntracks",numntbins,ntbins);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs ntracks",numntbins,ntbins);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs ntracks",numntbins,ntbins);

  TH1F* halphay = new TH1F("halphay","alpha vs ntracks",numntbins,ntbins);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs ntracks",numntbins,ntbins);
  TH1F* hmassy = new TH1F("hmassy","mass vs ntracks",numntbins,ntbins);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs ntracks",numntbins,ntbins);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs ntracks",numntbins,ntbins);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs ntracks",numntbins,ntbins);
  TH1F* herrmuy = new TH1F("herrmuy","errmu vs ntracks",numntbins,ntbins);
  TH1F* herrsigmay = new TH1F("herrsigmay","errsigma vs ntracks",numntbins,ntbins);
  TH1F* hlambday = new TH1F("hlambday","errlambda vs ntracks",numntbins,ntbins);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs ntracks",numntbins,ntbins);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs ntracks",numntbins,ntbins);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs ntracks",numntbins,ntbins);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs ntracks",numntbins,ntbins);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  yLow = -ybins[j+1];
  if (ybins[j]>0) {
    yHigh = -ybins[j];
  }
  else yHigh = 0;

  //1S backward loop
  for (int ipt = 0; ipt<numntbins; ipt++) {

    ntLow = ntbins[ipt];
    ntHigh = ntbins[ipt+1];

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntLow, ntHigh );
    NomFileName = Form("HFNtracksFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
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
  calpha->cd(1);
  halphapt->SetXTitle("Ntracks");
  halphapt->GetYaxis()->SetRangeUser(0,5);
  halphapt->Draw();
  //halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],400,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],400,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd(1);
  hf1spt->SetXTitle("Ntracks");
  hf1spt->GetYaxis()->SetRangeUser(-1,2);
  hf1spt->Draw();
  //hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],400,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],400,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd(1);
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

  cn1s->cd(1);
  hn1spt->SetXTitle("Ntracks");
  hn1spt->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1spt->Draw();
  //hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],400,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],400,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd(1);
  hsigma1spt->SetXTitle("Ntracks");
  hsigma1spt->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1spt->Draw();
  //hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],400,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],400,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd(1);
  hx1spt->SetXTitle("Ntracks");
  hx1spt->GetYaxis()->SetRangeUser(0.0,8);
  hx1spt->Draw();
  //hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],400,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],400,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd(1);
  herrmupt->SetXTitle("Ntracks");
  herrmupt->GetYaxis()->SetRangeUser(0.0,25);
  herrmupt->Draw();
  //TLine *errmuptlineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  //errmuptlineLow->SetLineColor(kRed);
  //errmuptlineLow->Draw();
  //TLine *errmuptlineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  //errmuptlineHigh->SetLineColor(kRed);
  //errmuptlineHigh->Draw();

  cerrsigma->cd(1);
  herrsigmapt->SetXTitle("Ntracks");
  herrsigmapt->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmapt->Draw();
  //TLine *errsigmaptlineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  //errsigmaptlineLow->SetLineColor(kRed);
  //errsigmaptlineLow->Draw();
  //TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  //errsigmaptlineHigh->SetLineColor(kRed);
  //errsigmaptlineHigh->Draw();

  clambda->cd(1);
  hlambdapt->SetXTitle("Ntracks");
  hlambdapt->GetYaxis()->SetRangeUser(0.0,25);
  hlambdapt->Draw();
  //TLine *lambdaptlineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  //lambdaptlineLow->SetLineColor(kRed);
  //lambdaptlineLow->Draw();
  //TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  //lambdaptlineHigh->SetLineColor(kRed);
  //lambdaptlineHigh->Draw();

  cnSig1s->cd(1);
  //gPad->SetLogy();
  hnSig1spt->SetXTitle("Ntracks");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  //hnSig1spt->Fit("pol1");
  TLine *nSig1sptlineLow = new TLine(0,0,400,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd(1);
  //gPad->SetLogy();
  hnSig2spt->SetXTitle("Ntracks");
  hnSig2spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig2spt->Draw();
  //hnSig2spt->Fit("pol1");
  TLine *nSig2sptlineLow = new TLine(0,0,400,0);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd(1);
  //gPad->SetLogy();
  hnSig3spt->SetXTitle("Ntracks");
  hnSig3spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig3spt->Draw();
  //hnSig3spt->Fit("pol1");
  TLine *nSig3sptlineLow = new TLine(0,0,400,0);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd(1);
  //gPad->SetLogy();
  hnBkgpt->SetXTitle("Ntracks");
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgpt->Draw();
  //hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,400,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  yLow = ybins[j];
  yHigh = ybins[j+1];

  //1S forward loop
  for (int iy = 0; iy<numntbins; iy++) {

    ntLow = ntbins[iy];
    ntHigh = ntbins[iy+1];

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntLow, ntHigh );
    NomFileName = Form("HFNtracksFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
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
  calpha->cd(2);
  halphay->SetXTitle("Ntracks");
  halphay->GetYaxis()->SetRangeUser(0,5);
  halphay->Draw();
  //halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(0,paramslower[2],400,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(0,paramsupper[2],400,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  cf1s->cd(2);
  hf1sy->SetXTitle("Ntracks");
  hf1sy->GetYaxis()->SetRangeUser(-1,2);
  hf1sy->Draw();
  //hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(0,paramslower[4],400,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(0,paramsupper[4],400,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  cmass->cd(2);
  hmassy->SetXTitle("Ntracks");
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
  //hmassy->Fit("pol1");
  TLine *massylinepdg = new TLine(0,9.46,400,9.46);
  massylinepdg->SetLineColor(3);
  massylinepdg->Draw();
  TLine *massylineLow = new TLine(0,9.36,400,9.36);
  massylineLow->SetLineColor(kRed);
  massylineLow->Draw();
  TLine *massylineHigh = new TLine(0,9.56,400,9.56);
  massylineHigh->SetLineColor(kRed);
  massylineHigh->Draw();

  cn1s->cd(2);
  hn1sy->SetXTitle("Ntracks");
  hn1sy->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1sy->Draw();
  //hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(0,paramslower[3],400,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(0,paramsupper[3],400,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  csigma1s->cd(2);
  hsigma1sy->SetXTitle("Ntracks");
  hsigma1sy->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1sy->Draw();
  //hsigma1sy->Fit("pol1");
  TLine *sigmaylineLow = new TLine(0,paramslower[0],400,paramslower[0]);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(0,paramsupper[0],400,paramsupper[0]);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();

  cx1s->cd(2);
  hx1sy->SetXTitle("Ntracks");
  hx1sy->GetYaxis()->SetRangeUser(0.0,8);
  hx1sy->Draw();
  //hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(0,paramslower[1],400,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(0,paramsupper[1],400,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

  cerrmu->cd(2);
  herrmuy->SetXTitle("Ntracks");
  herrmuy->GetYaxis()->SetRangeUser(0.0,25);
  herrmuy->Draw();
  //TLine *errmuylineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  //errmuylineLow->SetLineColor(kRed);
  //errmuylineLow->Draw();
  //TLine *errmuylineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  //errmuylineHigh->SetLineColor(kRed);
  //errmuylineHigh->Draw();

  cerrsigma->cd(2);
  herrsigmay->SetXTitle("Ntracks");
  herrsigmay->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmay->Draw();
  //TLine *errsigmaylineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  //errsigmaylineLow->SetLineColor(kRed);
  //errsigmaylineLow->Draw();
  //TLine *errsigmaylineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  //errsigmaylineHigh->SetLineColor(kRed);
  //errsigmaylineHigh->Draw();

  clambda->cd(2);
  hlambday->SetXTitle("Ntracks");
  hlambday->GetYaxis()->SetRangeUser(0.0,25);
  hlambday->Draw();
  //TLine *lambdaylineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  //lambdaylineLow->SetLineColor(kRed);
  //lambdaylineLow->Draw();
  //TLine *lambdaylineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  //lambdaylineHigh->SetLineColor(kRed);
  //lambdaylineHigh->Draw();

  cnSig1s->cd(2);
  //gPad->SetLogy();
  hnSig1sy->SetXTitle("Ntracks");
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sy->Draw();
  //hnSig1sy->Fit("pol1");
  TLine *nSig1sylineLow = new TLine(0,0,400,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();

  cnSig2s->cd(2);
  //gPad->SetLogy();
  hnSig2sy->SetXTitle("Ntracks");
  hnSig2sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig2sy->Draw();
  //hnSig2sy->Fit("pol1");
  TLine *nSig2sylineLow = new TLine(0,0,400,0);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();

  cnSig3s->cd(2);
  //gPad->SetLogy();
  hnSig3sy->SetXTitle("Ntracks");
  hnSig3sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig3sy->Draw();
  //hnSig3sy->Fit("pol1");
  TLine *nSig3sylineLow = new TLine(0,0,400,0);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();

  cnBkg->cd(2);
  //gPad->SetLogy();
  hnBkgy->SetXTitle("Ntracks");
  hnBkgy->GetYaxis()->SetRangeUser(10,10000*scale);
  hnBkgy->Draw();
  //hnBkgy->Fit("pol1");
  TLine *nBkgylineLow = new TLine(0,0,400,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();

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
