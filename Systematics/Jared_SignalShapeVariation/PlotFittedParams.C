#include <iostream>
#include "../HeaderFiles/rootFitHeaders.h"
#include "../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "../HeaderFiles/PsetCollection.h"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/StyleSetting.h"


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
    float ybins[9] = {-2.4,-1.67,-1.27,-0.87,-0.47,-0.07,0.33,0.73,1.46};
    float ybinsCM[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[5] = {-2.4,-1.27,-0.47,0.33,1.46};
    float ybinsCM[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[3] = {-2.4,-0.47,1.46};
    float ybinsCM[3] = {-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;

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

  //declare histograms
  TH1F* halphapt = new TH1F("halphapt","alpha vs pt",numptbins,ptbins);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs pt",numptbins,ptbins);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs pt",numptbins,ptbins);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs pt",numptbins,ptbins);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs pt",numptbins,ptbins);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs pt",numptbins,ptbins);
  TH1F* halphay = new TH1F("halphay","alpha vs y",numybins,ybinsCM);
  TH1F* hf1sy = new TH1F("hf1sy","f1s vs y",numybins,ybinsCM);
  TH1F* hmassy = new TH1F("hmassy","mass vs y",numybins,ybinsCM);
  TH1F* hn1sy = new TH1F("hn1sy","n1s vs y",numybins,ybinsCM);
  TH1F* hsigma1sy = new TH1F("hsigma1sy","sigma vs y",numybins,ybinsCM);
  TH1F* hx1sy = new TH1F("hx1sy","x1s vs y",numybins,ybinsCM);

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    float ptLow = ptbins[ipt];
    float ptHigh = ptbins[ipt+1];
    if (collId==kPADATA) {
      float yLow = -2.4;
      float yHigh = 1.46;
    }
    else if (collId==kPPDATA) {
      float yLow = -1.93;
      float yHigh = 1.93;
    }

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

  }

  //draw pt plots
  calpha->cd(1);
  halphapt->SetXTitle("pT");
  halphapt->GetYaxis()->SetRangeUser(0,5);
  halphapt->Draw();
  halphapt->Fit("pol1");
  //TLine *alphaptline = new TLine(x,y,x,y);
  //line->SetLineColor(kRed);
  //line->Draw();
  cf1s->cd(1);
  hf1spt->SetXTitle("pT");
  hf1spt->GetYaxis()->SetRangeUser(-1,1.2);
  hf1spt->Draw();
  hf1spt->Fit("pol1");
  cmass->cd(1);
  hmasspt->SetXTitle("pT");
  hmasspt->GetYaxis()->SetRangeUser(9.435,9.49);
  hmasspt->Draw();
  hmasspt->Fit("pol1");
  cn1s->cd(1);
  hn1spt->SetXTitle("pT");
  hn1spt->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1spt->Draw();
  hn1spt->Fit("pol1");
  csigma1s->cd(1);
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1spt->Draw();
  hsigma1spt->Fit("pol1");
  cx1s->cd(1);
  hx1spt->SetXTitle("pT");
  hx1spt->GetYaxis()->SetRangeUser(-0.5,4);
  hx1spt->Draw();
  hx1spt->Fit("pol1");


  //1S y loop
  for (int iy = 0; iy<numybins; iy++) {
    float ptLow = 0;
    float ptHigh = 30;
    if (collId==kPADATA) {
      float yLow = ybins[iy];
      float yHigh = ybins[iy+1];
    }
    else if (collId==kPPDATA) {
      float yLow = ybinsCM[iy];
      float yHigh = ybinsCM[iy+1];
    }
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
  }

  //draw rapidity plots
  calpha->cd(2);
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(0,5);
  halphay->Draw();
  halphay->Fit("pol1");
  cf1s->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(-1,1.2);
  hf1sy->Draw();
  hf1sy->Fit("pol1");
  cmass->cd(2);
  hmassy->SetXTitle("y");
  hmassy->GetYaxis()->SetRangeUser(9.435,9.49);
  hmassy->Draw();
  hmassy->Fit("pol1");
  cn1s->cd(2);
  hn1sy->SetXTitle("y");
  hn1sy->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1sy->Draw();
  hn1sy->Fit("pol1");
  csigma1s->cd(2);
  hsigma1sy->SetXTitle("y");
  hsigma1sy->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1sy->Draw();
  hsigma1sy->Fit("pol1");
  cx1s->cd(2);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(-0.5,4);
  hx1sy->Draw();
  hx1sy->Fit("pol1");

  //save plots
  calpha->SaveAs(Form("fitted_alpha_%isbins.png",whichUpsilon));
  cf1s->SaveAs(Form("fitted_f1s_%isbins.png",whichUpsilon));
  cmass->SaveAs(Form("fitted_mass_%isbins.png",whichUpsilon));
  cn1s->SaveAs(Form("fitted_n1s_%isbins.png",whichUpsilon));
  csigma1s->SaveAs(Form("fitted_sigma1s_%isbins.png",whichUpsilon));
  cx1s->SaveAs(Form("fitted_x1s_%isbins.png",whichUpsilon));

}
