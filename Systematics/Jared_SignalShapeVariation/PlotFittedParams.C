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


void PlotFittedParams() {

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

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
  TCanvas *cnSig1s = new TCanvas("cnSig1s","cnSig1s",4,45,800,400);
  cnSig1s->Divide(2,1);
  TCanvas *cnSig2s = new TCanvas("cnSig2s","cnSig2s",4,45,800,400);
  cnSig2s->Divide(2,1);
  TCanvas *cnSig3s = new TCanvas("cnSig3s","cnSig3s",4,45,800,400);
  cnSig3s->Divide(2,1);
  TCanvas *cnBkg = new TCanvas("cnBkg","cnBkg",4,45,800,400);
  cnBkg->Divide(2,1);

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs pt",numptbins,ptbins);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs pt",numptbins,ptbins);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs pt",numptbins,ptbins);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs pt",numptbins,ptbins);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs pt",numptbins,ptbins);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs pt",numptbins,ptbins);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs pt",numptbins,ptbins);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs pt",numptbins,ptbins);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs pt",numptbins,ptbins);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs pt",numptbins,ptbins);

  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybins);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybins);
  TH1F* hmassy = new TH1F("hmassy","mass vs y",numybins,ybins);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybins);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs y",numybins,ybins);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybins);
  TH1F* hnSig1sy = new TH1F("hnSig1sy","nSig1s vs y",numybins,ybins);
  TH1F* hnSig2sy = new TH1F("hnSig2sy","nSig2s vs y",numybins,ybins);
  TH1F* hnSig3sy = new TH1F("hnSig3sy","nSig3s vs y",numybins,ybins);
  TH1F* hnBkgy = new TH1F("hnBkgy","nBkg vs y",numybins,ybins);

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    float ptLow = ptbins[ipt];
    float ptHigh = ptbins[ipt+1];
    float yLow = -1.93;
    float yHigh = 1.93;

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *ws = NomFile->Get("workspace");  

    //extract parameter values
    float tempalpha = ws->var("alpha1s_1")->getVal();  
    float tempalphaerr = ws->var("alpha1s_1")->getError();
    float tempf1s = ws->var("f1s")->getVal();  
    float tempf1serr = ws->var("f1s")->getError();
    float tempmass = ws->var("m_{#Upsilon(1S)}")->getVal();  
    float tempmasserr = ws->var("m_{#Upsilon(1S)}")->getError();
    float tempn1s = ws->var("n1s_1")->getVal();  
    float tempn1serr = ws->var("n1s_1")->getError();
    float tempsigma1s = ws->var("sigma1s_1")->getVal();  
    float tempsigma1serr = ws->var("sigma1s_1")->getError();
    float tempx1s = ws->var("x1s")->getVal();  
    float tempx1serr = ws->var("x1s")->getError();
    float tempnSig1s = ws->var("nSig1s")->getVal();  
    float tempnSig1serr = ws->var("nSig1s")->getError();
    float tempnSig2s = ws->var("nSig2s")->getVal();  
    float tempnSig2serr = ws->var("nSig2s")->getError();
    float tempnSig3s = ws->var("nSig3s")->getVal();  
    float tempnSig3serr = ws->var("nSig3s")->getError();
    float tempnBkg = ws->var("nBkg")->getVal();  
    float tempnBkgerr = ws->var("nBkg")->getError();

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
  halphapt->SetXTitle("pT");
  halphapt->GetYaxis()->SetRangeUser(0,5);
  halphapt->Draw();
  halphapt->Fit("pol1");
  TLine *alphaptlineLow = new TLine(0,1.0,30,1.0);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,3.321,30,3.321);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd(1);
  hf1spt->SetXTitle("pT");
  hf1spt->GetYaxis()->SetRangeUser(-1,2);
  hf1spt->Draw();
  hf1spt->Fit("pol1");
  TLine *f1sptlineLow = new TLine(0,0,30,0);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,1,30,1);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd(1);
  hmasspt->SetXTitle("pT");
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  hmasspt->Fit("pol1");
  TLine *massptlineLow = new TLine(0,9.36,30,9.36);
  massptlineLow->SetLineColor(kRed);
  massptlineLow->Draw();
  TLine *massptlineHigh = new TLine(0,9.56,30,9.56);
  massptlineHigh->SetLineColor(kRed);
  massptlineHigh->Draw();

  cn1s->cd(1);
  hn1spt->SetXTitle("pT");
  hn1spt->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1spt->Draw();
  hn1spt->Fit("pol1");
  TLine *n1sptlineLow = new TLine(0,1.416,30,1.416);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,3.357,30,3.357);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd(1);
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1spt->Draw();
  hsigma1spt->Fit("pol1");
  TLine *sigmaptlineLow = new TLine(0,0.02,30,0.02);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,0.3,30,0.3);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd(1);
  hx1spt->SetXTitle("pT");
  hx1spt->GetYaxis()->SetRangeUser(-0.5,4);
  hx1spt->Draw();
  hx1spt->Fit("pol1");
  TLine *x1sptlineLow = new TLine(0,0,30,0);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,1,30,1);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

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

  cnBkg->cd(1);
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->GetYaxis()->SetRangeUser(0,10000);
  hnBkgpt->Draw();
  hnBkgpt->Fit("pol1");
  TLine *nBkgptlineLow = new TLine(0,0,30,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    float ptLow = 0;
    float ptHigh = 30;
    float yLow = ybins[iy];
    float yHigh = ybins[iy+1];

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *yws = NomFile->Get("workspace");  

    //extract parameter values
    float tempalpha = yws->var("alpha1s_1")->getVal();  
    float tempalphaerr = yws->var("alpha1s_1")->getError();
    float tempf1s = yws->var("f1s")->getVal();  
    float tempf1serr = yws->var("f1s")->getError();
    float tempmass = yws->var("m_{#Upsilon(1S)}")->getVal();  
    float tempmasserr = yws->var("m_{#Upsilon(1S)}")->getError();
    float tempn1s = yws->var("n1s_1")->getVal();  
    float tempn1serr = yws->var("n1s_1")->getError();
    float tempsigma1s = yws->var("sigma1s_1")->getVal();  
    float tempsigma1serr = yws->var("sigma1s_1")->getError();
    float tempx1s = yws->var("x1s")->getVal();  
    float tempx1serr = yws->var("x1s")->getError();
    float tempnSig1s = yws->var("nSig1s")->getVal();  
    float tempnSig1serr = yws->var("nSig1s")->getError();
    float tempnSig2s = yws->var("nSig2s")->getVal();  
    float tempnSig2serr = yws->var("nSig2s")->getError();
    float tempnSig3s = yws->var("nSig3s")->getVal();  
    float tempnSig3serr = yws->var("nSig3s")->getError();
    float tempnBkg = yws->var("nBkg")->getVal();  
    float tempnBkgerr = yws->var("nBkg")->getError();

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
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(0,5);
  halphay->Draw();
  halphay->Fit("pol1");
  TLine *alphaylineLow = new TLine(-1.93,1.0,1.93,1.0);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(-1.93,3.321,1.93,3.321);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  cf1s->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(-1,2);
  hf1sy->Draw();
  hf1sy->Fit("pol1");
  TLine *f1sylineLow = new TLine(-1.93,0,1.93,0);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(-1.93,1,1.93,1);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  cmass->cd(2);
  hmassy->SetXTitle("y");
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
  hmassy->Fit("pol1");
  TLine *massylineLow = new TLine(-1.93,9.36,1.93,9.36);
  massylineLow->SetLineColor(kRed);
  massylineLow->Draw();
  TLine *massylineHigh = new TLine(-1.93,9.56,1.93,9.56);
  massylineHigh->SetLineColor(kRed);
  massylineHigh->Draw();

  cn1s->cd(2);
  hn1sy->SetXTitle("y");
  hn1sy->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1sy->Draw();
  hn1sy->Fit("pol1");
  TLine *n1sylineLow = new TLine(-1.93,1.416,1.93,1.416);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(-1.93,3.357,1.93,3.357);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  csigma1s->cd(2);
  hsigma1sy->SetXTitle("y");
  hsigma1sy->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1sy->Draw();
  hsigma1sy->Fit("pol1");
  TLine *sigmaylineLow = new TLine(-1.93,0.02,1.93,0.02);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(-1.93,0.3,1.93,0.3);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();

  cx1s->cd(2);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(-0.5,4);
  hx1sy->Draw();
  hx1sy->Fit("pol1");
  TLine *x1sylineLow = new TLine(-1.93,0,1.93,0);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(-1.93,1,1.93,1);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

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

  cnBkg->cd(2);
  hnBkgy->SetXTitle("y");
  hnBkgy->GetYaxis()->SetRangeUser(0,10000);
  hnBkgy->Draw();
  hnBkgy->Fit("pol1");
  TLine *nBkgylineLow = new TLine(-1.93,0,1.93,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();

  //save plots
  calpha->SaveAs(Form("fitted_alpha_%isbins.png",whichUpsilon));
  cf1s->SaveAs(Form("fitted_f1s_%isbins.png",whichUpsilon));
  cmass->SaveAs(Form("fitted_mass_%isbins.png",whichUpsilon));
  cn1s->SaveAs(Form("fitted_n1s_%isbins.png",whichUpsilon));
  csigma1s->SaveAs(Form("fitted_sigma1s_%isbins.png",whichUpsilon));
  cx1s->SaveAs(Form("fitted_x1s_%isbins.png",whichUpsilon));
  cnSig1s->SaveAs(Form("fitted_nSig1s_%isbins.png",whichUpsilon));
  cnSig2s->SaveAs(Form("fitted_nSig2s_%isbins.png",whichUpsilon));
  cnSig3s->SaveAs(Form("fitted_nSig3s_%isbins.png",whichUpsilon));
  cnBkg->SaveAs(Form("fitted_nBkg_%isbins.png",whichUpsilon));

}
