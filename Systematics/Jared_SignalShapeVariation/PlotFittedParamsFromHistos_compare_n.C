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


void PlotFittedParamsFromHistos_compare_n() {

  int whichUpsilon = 1;
  int collId = kPPDATA;
  float scale = whichUpsilon;
  if (collId==kPPDATA) scale = scale*8;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.2, 5.0, 3.321, 5.0, 1.0, 25.0, 25.0, 25.0};
  double paramslower[8] = {0.02, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetErrorX(0);
  gStyle->SetMarkerStyle(2);

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

  //Extract histograms
  TString HistoFileName = Form("FittedParamsHistos_PP_%is_fixedn%i.root",whichUpsilon,2);
  cout << HistoFileName << endl;
  HistoFile = TFile::Open(HistoFileName,"READ");

  TH1F* halphapt = (TH1F*)HistoFile->Get("halphapt");
  TH1F* hf1spt = (TH1F*)HistoFile->Get("hf1spt");
  TH1F* hmasspt = (TH1F*)HistoFile->Get("hmasspt");
  TH1F* hn1spt = (TH1F*)HistoFile->Get("hn1spt");
  TH1F* hsigma1spt = (TH1F*)HistoFile->Get("hsigma1spt");
  TH1F* hx1spt = (TH1F*)HistoFile->Get("hx1spt");
  TH1F* herrmupt = (TH1F*)HistoFile->Get("herrmupt");
  TH1F* herrsigmapt = (TH1F*)HistoFile->Get("herrsigmapt");
  TH1F* hlambdapt = (TH1F*)HistoFile->Get("hlambdapt");
  TH1F* hnSig1spt = (TH1F*)HistoFile->Get("hnSig1spt");
  TH1F* hnSig2spt = (TH1F*)HistoFile->Get("hnSig2spt");
  TH1F* hnSig3spt = (TH1F*)HistoFile->Get("hnSig3spt");
  TH1F* hnBkgpt = (TH1F*)HistoFile->Get("hnBkgpt");

  TH1F* halphay = (TH1F*)HistoFile->Get("halphay");
  TH1F* hf1sy = (TH1F*)HistoFile->Get("hf1sy");
  TH1F* hmassy = (TH1F*)HistoFile->Get("hmassy");
  TH1F* hn1sy = (TH1F*)HistoFile->Get("hn1sy");
  TH1F* hsigma1sy = (TH1F*)HistoFile->Get("hsigma1sy");
  TH1F* hx1sy = (TH1F*)HistoFile->Get("hx1sy");
  TH1F* herrmuy = (TH1F*)HistoFile->Get("herrmuy");
  TH1F* herrsigmay = (TH1F*)HistoFile->Get("herrsigmay");
  TH1F* hlambday = (TH1F*)HistoFile->Get("hlambday");
  TH1F* hnSig1sy = (TH1F*)HistoFile->Get("hnSig1sy");
  TH1F* hnSig2sy = (TH1F*)HistoFile->Get("hnSig2sy");
  TH1F* hnSig3sy = (TH1F*)HistoFile->Get("hnSig3sy");
  TH1F* hnBkgy = (TH1F*)HistoFile->Get("hnBkgy");

  TString NewHistoFileName = Form("FittedParamsHistos_PP_%is_fixedn%i.root",whichUpsilon,3);
  cout << NewHistoFileName << endl;
  NewHistoFile = TFile::Open(NewHistoFileName,"READ");

  TH1F* Newhalphapt = (TH1F*)NewHistoFile->Get("halphapt");
  TH1F* Newhf1spt = (TH1F*)NewHistoFile->Get("hf1spt");
  TH1F* Newhmasspt = (TH1F*)NewHistoFile->Get("hmasspt");
  TH1F* Newhn1spt = (TH1F*)NewHistoFile->Get("hn1spt");
  TH1F* Newhsigma1spt = (TH1F*)NewHistoFile->Get("hsigma1spt");
  TH1F* Newhx1spt = (TH1F*)NewHistoFile->Get("hx1spt");
  TH1F* Newherrmupt = (TH1F*)NewHistoFile->Get("herrmupt");
  TH1F* Newherrsigmapt = (TH1F*)NewHistoFile->Get("herrsigmapt");
  TH1F* Newhlambdapt = (TH1F*)NewHistoFile->Get("hlambdapt");
  TH1F* NewhnSig1spt = (TH1F*)NewHistoFile->Get("hnSig1spt");
  TH1F* NewhnSig2spt = (TH1F*)NewHistoFile->Get("hnSig2spt");
  TH1F* NewhnSig3spt = (TH1F*)NewHistoFile->Get("hnSig3spt");
  TH1F* NewhnBkgpt = (TH1F*)NewHistoFile->Get("hnBkgpt");

  TH1F* Newhalphay = (TH1F*)NewHistoFile->Get("halphay");
  TH1F* Newhf1sy = (TH1F*)NewHistoFile->Get("hf1sy");
  TH1F* Newhmassy = (TH1F*)NewHistoFile->Get("hmassy");
  TH1F* Newhn1sy = (TH1F*)NewHistoFile->Get("hn1sy");
  TH1F* Newhsigma1sy = (TH1F*)NewHistoFile->Get("hsigma1sy");
  TH1F* Newhx1sy = (TH1F*)NewHistoFile->Get("hx1sy");
  TH1F* Newherrmuy = (TH1F*)NewHistoFile->Get("herrmuy");
  TH1F* Newherrsigmay = (TH1F*)NewHistoFile->Get("herrsigmay");
  TH1F* Newhlambday = (TH1F*)NewHistoFile->Get("hlambday");
  TH1F* NewhnSig1sy = (TH1F*)NewHistoFile->Get("hnSig1sy");
  TH1F* NewhnSig2sy = (TH1F*)NewHistoFile->Get("hnSig2sy");
  TH1F* NewhnSig3sy = (TH1F*)NewHistoFile->Get("hnSig3sy");
  TH1F* NewhnBkgy = (TH1F*)NewHistoFile->Get("hnBkgy");

  TString NewerHistoFileName = Form("FittedParamsHistos_PP_%is_fixedn%i.root",whichUpsilon,4);
  cout << NewerHistoFileName << endl;
  NewerHistoFile = TFile::Open(NewerHistoFileName,"READ");

  TH1F* Newerhalphapt = (TH1F*)NewerHistoFile->Get("halphapt");
  TH1F* Newerhf1spt = (TH1F*)NewerHistoFile->Get("hf1spt");
  TH1F* Newerhmasspt = (TH1F*)NewerHistoFile->Get("hmasspt");
  TH1F* Newerhn1spt = (TH1F*)NewerHistoFile->Get("hn1spt");
  TH1F* Newerhsigma1spt = (TH1F*)NewerHistoFile->Get("hsigma1spt");
  TH1F* Newerhx1spt = (TH1F*)NewerHistoFile->Get("hx1spt");
  TH1F* Newerherrmupt = (TH1F*)NewerHistoFile->Get("herrmupt");
  TH1F* Newerherrsigmapt = (TH1F*)NewerHistoFile->Get("herrsigmapt");
  TH1F* Newerhlambdapt = (TH1F*)NewerHistoFile->Get("hlambdapt");
  TH1F* NewerhnSig1spt = (TH1F*)NewerHistoFile->Get("hnSig1spt");
  TH1F* NewerhnSig2spt = (TH1F*)NewerHistoFile->Get("hnSig2spt");
  TH1F* NewerhnSig3spt = (TH1F*)NewerHistoFile->Get("hnSig3spt");
  TH1F* NewerhnBkgpt = (TH1F*)NewerHistoFile->Get("hnBkgpt");

  TH1F* Newerhalphay = (TH1F*)NewerHistoFile->Get("halphay");
  TH1F* Newerhf1sy = (TH1F*)NewerHistoFile->Get("hf1sy");
  TH1F* Newerhmassy = (TH1F*)NewerHistoFile->Get("hmassy");
  TH1F* Newerhn1sy = (TH1F*)NewerHistoFile->Get("hn1sy");
  TH1F* Newerhsigma1sy = (TH1F*)NewerHistoFile->Get("hsigma1sy");
  TH1F* Newerhx1sy = (TH1F*)NewerHistoFile->Get("hx1sy");
  TH1F* Newerherrmuy = (TH1F*)NewerHistoFile->Get("herrmuy");
  TH1F* Newerherrsigmay = (TH1F*)NewerHistoFile->Get("herrsigmay");
  TH1F* Newerhlambday = (TH1F*)NewerHistoFile->Get("hlambday");
  TH1F* NewerhnSig1sy = (TH1F*)NewerHistoFile->Get("hnSig1sy");
  TH1F* NewerhnSig2sy = (TH1F*)NewerHistoFile->Get("hnSig2sy");
  TH1F* NewerhnSig3sy = (TH1F*)NewerHistoFile->Get("hnSig3sy");
  TH1F* NewerhnBkgy = (TH1F*)NewerHistoFile->Get("hnBkgy");

  halphapt->SetMarkerStyle(2);
  hf1spt->SetMarkerStyle(2);
  hmasspt->SetMarkerStyle(2);
  hn1spt->SetMarkerStyle(2);
  hsigma1spt->SetMarkerStyle(2);
  hx1spt->SetMarkerStyle(2);
  herrmupt->SetMarkerStyle(2);
  herrsigmapt->SetMarkerStyle(2);
  hlambdapt->SetMarkerStyle(2);
  hnSig1spt->SetMarkerStyle(2);
  hnSig2spt->SetMarkerStyle(2);
  hnSig3spt->SetMarkerStyle(2);
  hnBkgpt->SetMarkerStyle(2);

  Newhalphapt->SetMarkerStyle(2);
  Newhf1spt->SetMarkerStyle(2);
  Newhmasspt->SetMarkerStyle(2);
  Newhn1spt->SetMarkerStyle(2);
  Newhsigma1spt->SetMarkerStyle(2);
  Newhx1spt->SetMarkerStyle(2);
  Newherrmupt->SetMarkerStyle(2);
  Newherrsigmapt->SetMarkerStyle(2);
  Newhlambdapt->SetMarkerStyle(2);
  NewhnSig1spt->SetMarkerStyle(2);
  NewhnSig2spt->SetMarkerStyle(2);
  NewhnSig3spt->SetMarkerStyle(2);
  NewhnBkgpt->SetMarkerStyle(2);

  Newerhalphapt->SetMarkerStyle(2);
  Newerhf1spt->SetMarkerStyle(2);
  Newerhmasspt->SetMarkerStyle(2);
  Newerhn1spt->SetMarkerStyle(2);
  Newerhsigma1spt->SetMarkerStyle(2);
  Newerhx1spt->SetMarkerStyle(2);
  Newerherrmupt->SetMarkerStyle(2);
  Newerherrsigmapt->SetMarkerStyle(2);
  Newerhlambdapt->SetMarkerStyle(2);
  NewerhnSig1spt->SetMarkerStyle(2);
  NewerhnSig2spt->SetMarkerStyle(2);
  NewerhnSig3spt->SetMarkerStyle(2);
  NewerhnBkgpt->SetMarkerStyle(2);

  halphay->SetMarkerStyle(2);
  hf1sy->SetMarkerStyle(2);
  hmassy->SetMarkerStyle(2);
  hn1sy->SetMarkerStyle(2);
  hsigma1sy->SetMarkerStyle(2);
  hx1sy->SetMarkerStyle(2);
  herrmuy->SetMarkerStyle(2);
  herrsigmay->SetMarkerStyle(2);
  hlambday->SetMarkerStyle(2);
  hnSig1sy->SetMarkerStyle(2);
  hnSig2sy->SetMarkerStyle(2);
  hnSig3sy->SetMarkerStyle(2);
  hnBkgy->SetMarkerStyle(2);

  Newhalphay->SetMarkerStyle(2);
  Newhf1sy->SetMarkerStyle(2);
  Newhmassy->SetMarkerStyle(2);
  Newhn1sy->SetMarkerStyle(2);
  Newhsigma1sy->SetMarkerStyle(2);
  Newhx1sy->SetMarkerStyle(2);
  Newherrmuy->SetMarkerStyle(2);
  Newherrsigmay->SetMarkerStyle(2);
  Newhlambday->SetMarkerStyle(2);
  NewhnSig1sy->SetMarkerStyle(2);
  NewhnSig2sy->SetMarkerStyle(2);
  NewhnSig3sy->SetMarkerStyle(2);
  NewhnBkgy->SetMarkerStyle(2);

  Newerhalphay->SetMarkerStyle(2);
  Newerhf1sy->SetMarkerStyle(2);
  Newerhmassy->SetMarkerStyle(2);
  Newerhn1sy->SetMarkerStyle(2);
  Newerhsigma1sy->SetMarkerStyle(2);
  Newerhx1sy->SetMarkerStyle(2);
  Newerherrmuy->SetMarkerStyle(2);
  Newerherrsigmay->SetMarkerStyle(2);
  Newerhlambday->SetMarkerStyle(2);
  NewerhnSig1sy->SetMarkerStyle(2);
  NewerhnSig2sy->SetMarkerStyle(2);
  NewerhnSig3sy->SetMarkerStyle(2);
  NewerhnBkgy->SetMarkerStyle(2);

  //draw pt plots
  calpha->cd(1);
  halphapt->SetXTitle("pT");
  halphapt->GetYaxis()->SetRangeUser(0,5);
  halphapt->Draw();
  halphapt->SetMarkerColor(3);
  halphapt->SetLineColor(3);
  Newhalphapt->Draw("same");
  Newhalphapt->SetMarkerColor(6);
  Newhalphapt->SetLineColor(6);
  Newerhalphapt->Draw("same");
  TLine *alphaptlineLow = new TLine(0,paramslower[2],30,paramslower[2]);
  alphaptlineLow->SetLineColor(kRed);
  alphaptlineLow->Draw();
  TLine *alphaptlineHigh = new TLine(0,paramsupper[2],30,paramsupper[2]);
  alphaptlineHigh->SetLineColor(kRed);
  alphaptlineHigh->Draw();

  cf1s->cd(1);
  hf1spt->SetXTitle("pT");
  hf1spt->GetYaxis()->SetRangeUser(-1,2);
  hf1spt->Draw();
  hf1spt->SetMarkerColor(3);
  hf1spt->SetLineColor(3);
  Newhf1spt->Draw("same");
  Newhf1spt->SetMarkerColor(6);
  Newhf1spt->SetLineColor(6);
  Newerhf1spt->Draw("same");
  TLine *f1sptlineLow = new TLine(0,paramslower[4],30,paramslower[4]);
  f1sptlineLow->SetLineColor(kRed);
  f1sptlineLow->Draw();
  TLine *f1sptlineHigh = new TLine(0,paramsupper[4],30,paramsupper[4]);
  f1sptlineHigh->SetLineColor(kRed);
  f1sptlineHigh->Draw();

  cmass->cd(1);
  hmasspt->SetXTitle("pT");
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  hmasspt->SetMarkerColor(3);
  hmasspt->SetLineColor(3);
  Newhmasspt->Draw("same");
  Newhmasspt->SetMarkerColor(6);
  Newhmasspt->SetLineColor(6);
  Newerhmasspt->Draw("same");
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
  hn1spt->SetMarkerColor(3);
  hn1spt->SetLineColor(3);
  Newhn1spt->Draw("same");
  Newhn1spt->SetMarkerColor(6);
  Newhn1spt->SetLineColor(6);
  Newerhn1spt->Draw("same");
  TLine *n1sptlineLow = new TLine(0,paramslower[3],30,paramslower[3]);
  n1sptlineLow->SetLineColor(kRed);
  n1sptlineLow->Draw();
  TLine *n1sptlineHigh = new TLine(0,paramsupper[3],30,paramsupper[3]);
  n1sptlineHigh->SetLineColor(kRed);
  n1sptlineHigh->Draw();

  csigma1s->cd(1);
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1spt->Draw();
  hsigma1spt->SetMarkerColor(3);
  hsigma1spt->SetLineColor(3);
  Newhsigma1spt->Draw("same");
  Newhsigma1spt->SetMarkerColor(6);
  Newhsigma1spt->SetLineColor(6);
  Newerhsigma1spt->Draw("same");
  TLine *sigmaptlineLow = new TLine(0,paramslower[0],30,paramslower[0]);
  sigmaptlineLow->SetLineColor(kRed);
  sigmaptlineLow->Draw();
  TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],30,paramsupper[0]);
  sigmaptlineHigh->SetLineColor(kRed);
  sigmaptlineHigh->Draw();

  cx1s->cd(1);
  hx1spt->SetXTitle("pT");
  hx1spt->GetYaxis()->SetRangeUser(0.0,8);
  hx1spt->Draw();
  hx1spt->SetMarkerColor(3);
  hx1spt->SetLineColor(3);
  Newhx1spt->Draw("same");
  Newhx1spt->SetMarkerColor(6);
  Newhx1spt->SetLineColor(6);
  Newerhx1spt->Draw("same");
  TLine *x1sptlineLow = new TLine(0,paramslower[1],30,paramslower[1]);
  x1sptlineLow->SetLineColor(kRed);
  x1sptlineLow->Draw();
  TLine *x1sptlineHigh = new TLine(0,paramsupper[1],30,paramsupper[1]);
  x1sptlineHigh->SetLineColor(kRed);
  x1sptlineHigh->Draw();

  cerrmu->cd(1);
  herrmupt->SetXTitle("pT");
  herrmupt->GetYaxis()->SetRangeUser(0.0,25);
  herrmupt->Draw();
  herrmupt->SetMarkerColor(3);
  herrmupt->SetLineColor(3);
  Newherrmupt->Draw("same");
  Newherrmupt->SetMarkerColor(6);
  Newherrmupt->SetLineColor(6);
  Newerherrmupt->Draw("same");
  TLine *errmuptlineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  errmuptlineLow->SetLineColor(kRed);
  errmuptlineLow->Draw();
  TLine *errmuptlineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  errmuptlineHigh->SetLineColor(kRed);
  errmuptlineHigh->Draw();

  cerrsigma->cd(1);
  herrsigmapt->SetXTitle("pT");
  herrsigmapt->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmapt->Draw();
  herrsigmapt->SetMarkerColor(3);
  herrsigmapt->SetLineColor(3);
  Newherrsigmapt->Draw("same");
  Newherrsigmapt->SetMarkerColor(6);
  Newherrsigmapt->SetLineColor(6);
  Newerherrsigmapt->Draw("same");
  TLine *errsigmaptlineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  errsigmaptlineLow->SetLineColor(kRed);
  errsigmaptlineLow->Draw();
  TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  errsigmaptlineHigh->SetLineColor(kRed);
  errsigmaptlineHigh->Draw();

  clambda->cd(1);
  hlambdapt->SetXTitle("pT");
  hlambdapt->GetYaxis()->SetRangeUser(0.0,25);
  hlambdapt->Draw();
  hlambdapt->SetMarkerColor(3);
  hlambdapt->SetLineColor(3);
  Newhlambdapt->Draw("same");
  Newhlambdapt->SetMarkerColor(6);
  Newhlambdapt->SetLineColor(6);
  Newerhlambdapt->Draw("same");
  TLine *lambdaptlineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  lambdaptlineLow->SetLineColor(kRed);
  lambdaptlineLow->Draw();
  TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  lambdaptlineHigh->SetLineColor(kRed);
  lambdaptlineHigh->Draw();

  cnSig1s->cd(1);
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1spt->Draw();
  hnSig1spt->SetMarkerColor(3);
  hnSig1spt->SetLineColor(3);
  NewhnSig1spt->Draw("same");
  NewhnSig1spt->SetMarkerColor(6);
  NewhnSig1spt->SetLineColor(6);
  NewerhnSig1spt->Draw("same");
  TLine *nSig1sptlineLow = new TLine(0,0,30,0);
  nSig1sptlineLow->SetLineColor(kRed);
  nSig1sptlineLow->Draw();

  cnSig2s->cd(1);
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->GetYaxis()->SetRangeUser(0,800*scale);
  hnSig2spt->Draw();
  hnSig2spt->SetMarkerColor(3);
  hnSig2spt->SetLineColor(3);
  NewhnSig2spt->Draw("same");
  NewhnSig2spt->SetMarkerColor(6);
  NewhnSig2spt->SetLineColor(6);
  NewerhnSig2spt->Draw("same");
  TLine *nSig2sptlineLow = new TLine(0,0,30,0);
  nSig2sptlineLow->SetLineColor(kRed);
  nSig2sptlineLow->Draw();

  cnSig3s->cd(1);
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->GetYaxis()->SetRangeUser(0,800*scale);
  hnSig3spt->Draw();
  hnSig3spt->SetMarkerColor(3);
  hnSig3spt->SetLineColor(3);
  NewhnSig3spt->Draw("same");
  NewhnSig3spt->SetMarkerColor(6);
  NewhnSig3spt->SetLineColor(6);
  NewerhnSig3spt->Draw("same");
  TLine *nSig3sptlineLow = new TLine(0,0,30,0);
  nSig3sptlineLow->SetLineColor(kRed);
  nSig3sptlineLow->Draw();

  cnBkg->cd(1);
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->GetYaxis()->SetRangeUser(10,5000*scale);
  hnBkgpt->Draw();
  hnBkgpt->SetMarkerColor(3);
  hnBkgpt->SetLineColor(3);
  NewhnBkgpt->Draw("same");
  NewhnBkgpt->SetMarkerColor(6);
  NewhnBkgpt->SetLineColor(6);
  NewerhnBkgpt->Draw("same");
  TLine *nBkgptlineLow = new TLine(0,0,30,0);
  nBkgptlineLow->SetLineColor(kRed);
  nBkgptlineLow->Draw();


  //draw rapidity plots
  calpha->cd(2);
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(0,5);
  halphay->Draw();
  halphay->SetMarkerColor(3);
  halphay->SetLineColor(3);
  Newhalphay->Draw("same");
  Newhalphay->SetMarkerColor(6);
  Newhalphay->SetLineColor(6);
  Newerhalphay->Draw("same");
  TLine *alphaylineLow = new TLine(-1.93,paramslower[2],1.93,paramslower[2]);
  alphaylineLow->SetLineColor(kRed);
  alphaylineLow->Draw();
  TLine *alphaylineHigh = new TLine(-1.93,paramsupper[2],1.93,paramsupper[2]);
  alphaylineHigh->SetLineColor(kRed);
  alphaylineHigh->Draw();

  cf1s->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(-1,2);
  hf1sy->Draw();
  hf1sy->SetMarkerColor(3);
  hf1sy->SetLineColor(3);
  Newhf1sy->Draw("same");
  Newhf1sy->SetMarkerColor(6);
  Newhf1sy->SetLineColor(6);
  Newerhf1sy->Draw("same");
  TLine *f1sylineLow = new TLine(-1.93,paramslower[4],1.93,paramslower[4]);
  f1sylineLow->SetLineColor(kRed);
  f1sylineLow->Draw();
  TLine *f1sylineHigh = new TLine(-1.93,paramsupper[4],1.93,paramsupper[4]);
  f1sylineHigh->SetLineColor(kRed);
  f1sylineHigh->Draw();

  cmass->cd(2);
  hmassy->SetXTitle("y");
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
  hmassy->SetMarkerColor(3);
  hmassy->SetLineColor(3);
  Newhmassy->Draw("same");
  Newhmassy->SetMarkerColor(6);
  Newhmassy->SetLineColor(6);
  Newerhmassy->Draw("same");
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
  hn1sy->SetMarkerColor(3);
  hn1sy->SetLineColor(3);
  Newhn1sy->Draw("same");
  Newhn1sy->SetMarkerColor(6);
  Newhn1sy->SetLineColor(6);
  Newerhn1sy->Draw("same");
  TLine *n1sylineLow = new TLine(-1.93,paramslower[3],1.93,paramslower[3]);
  n1sylineLow->SetLineColor(kRed);
  n1sylineLow->Draw();
  TLine *n1sylineHigh = new TLine(-1.93,paramsupper[3],1.93,paramsupper[3]);
  n1sylineHigh->SetLineColor(kRed);
  n1sylineHigh->Draw();

  csigma1s->cd(2);
  hsigma1sy->SetXTitle("y");
  hsigma1sy->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1sy->Draw();
  hsigma1sy->SetMarkerColor(3);
  hsigma1sy->SetLineColor(3);
  Newhsigma1sy->Draw("same");
  Newhsigma1sy->SetMarkerColor(6);
  Newhsigma1sy->SetLineColor(6);
  Newerhsigma1sy->Draw("same");
  TLine *sigmaylineLow = new TLine(-1.93,paramslower[0],1.93,paramslower[0]);
  sigmaylineLow->SetLineColor(kRed);
  sigmaylineLow->Draw();
  TLine *sigmaylineHigh = new TLine(-1.93,paramsupper[0],1.93,paramsupper[0]);
  sigmaylineHigh->SetLineColor(kRed);
  sigmaylineHigh->Draw();

  cx1s->cd(2);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(0.0,8);
  hx1sy->Draw();
  hx1sy->SetMarkerColor(3);
  hx1sy->SetLineColor(3);
  Newhx1sy->Draw("same");
  Newhx1sy->SetMarkerColor(6);
  Newhx1sy->SetLineColor(6);
  Newerhx1sy->Draw("same");
  TLine *x1sylineLow = new TLine(-1.93,paramslower[1],1.93,paramslower[1]);
  x1sylineLow->SetLineColor(kRed);
  x1sylineLow->Draw();
  TLine *x1sylineHigh = new TLine(-1.93,paramsupper[1],1.93,paramsupper[1]);
  x1sylineHigh->SetLineColor(kRed);
  x1sylineHigh->Draw();

  cerrmu->cd(2);
  herrmuy->SetXTitle("y");
  herrmuy->GetYaxis()->SetRangeUser(0.0,25);
  herrmuy->Draw();
  herrmuy->SetMarkerColor(3);
  herrmuy->SetLineColor(3);
  Newherrmuy->Draw("same");
  Newherrmuy->SetMarkerColor(6);
  Newherrmuy->SetLineColor(6);
  Newerherrmuy->Draw("same");
  TLine *errmuylineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  errmuylineLow->SetLineColor(kRed);
  errmuylineLow->Draw();
  TLine *errmuylineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  errmuylineHigh->SetLineColor(kRed);
  errmuylineHigh->Draw();

  cerrsigma->cd(2);
  herrsigmay->SetXTitle("y");
  herrsigmay->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmay->Draw();
  herrsigmay->SetMarkerColor(3);
  herrsigmay->SetLineColor(3);
  Newherrsigmay->Draw("same");
  Newherrsigmay->SetMarkerColor(6);
  Newherrsigmay->SetLineColor(6);
  Newerherrsigmay->Draw("same");
  TLine *errsigmaylineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  errsigmaylineLow->SetLineColor(kRed);
  errsigmaylineLow->Draw();
  TLine *errsigmaylineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  errsigmaylineHigh->SetLineColor(kRed);
  errsigmaylineHigh->Draw();

  clambda->cd(2);
  hlambday->SetXTitle("y");
  hlambday->GetYaxis()->SetRangeUser(0.0,25);
  hlambday->Draw();
  hlambday->SetMarkerColor(3);
  hlambday->SetLineColor(3);
  Newhlambday->Draw("same");
  Newhlambday->SetMarkerColor(6);
  Newhlambday->SetLineColor(6);
  Newerhlambday->Draw("same");
  TLine *lambdaylineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  lambdaylineLow->SetLineColor(kRed);
  lambdaylineLow->Draw();
  TLine *lambdaylineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  lambdaylineHigh->SetLineColor(kRed);
  lambdaylineHigh->Draw();

  cnSig1s->cd(2);
  hnSig1sy->SetXTitle("y");
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*scale);
  hnSig1sy->Draw();
  hnSig1sy->SetMarkerColor(3);
  hnSig1sy->SetLineColor(3);
  NewhnSig1sy->Draw("same");
  NewhnSig1sy->SetMarkerColor(6);
  NewhnSig1sy->SetLineColor(6);
  NewerhnSig1sy->Draw("same");
  TLine *nSig1sylineLow = new TLine(-1.93,0,1.93,0);
  nSig1sylineLow->SetLineColor(kRed);
  nSig1sylineLow->Draw();

  cnSig2s->cd(2);
  hnSig2sy->SetXTitle("y");
  hnSig2sy->GetYaxis()->SetRangeUser(0,800*scale);
  hnSig2sy->Draw();
  hnSig2sy->SetMarkerColor(3);
  hnSig2sy->SetLineColor(3);
  NewhnSig2sy->Draw("same");
  NewhnSig2sy->SetMarkerColor(6);
  NewhnSig2sy->SetLineColor(6);
  NewerhnSig2sy->Draw("same");
  TLine *nSig2sylineLow = new TLine(-1.93,0,1.93,0);
  nSig2sylineLow->SetLineColor(kRed);
  nSig2sylineLow->Draw();

  cnSig3s->cd(2);
  hnSig3sy->SetXTitle("y");
  hnSig3sy->GetYaxis()->SetRangeUser(0,800*scale);
  hnSig3sy->Draw();
  hnSig3sy->SetMarkerColor(3);
  hnSig3sy->SetLineColor(3);
  NewhnSig3sy->Draw("same");
  NewhnSig3sy->SetMarkerColor(6);
  NewhnSig3sy->SetLineColor(6);
  NewerhnSig3sy->Draw("same");
  TLine *nSig3sylineLow = new TLine(-1.93,0,1.93,0);
  nSig3sylineLow->SetLineColor(kRed);
  nSig3sylineLow->Draw();
  TLegend *nSig3sLegend = new TLegend(0.4,0.7,0.6,0.9);
  nSig3sLegend->AddEntry(hnSig3sy,"n=2","p");
  nSig3sLegend->AddEntry(NewhnSig3sy,"n=3","p");
  nSig3sLegend->AddEntry(NewerhnSig3sy,"n=4","p");
  nSig3sLegend->Draw("same");

  cnBkg->cd(2);
  hnBkgy->SetXTitle("y");
  hnBkgy->GetYaxis()->SetRangeUser(10,5000*scale);
  hnBkgy->Draw();
  hnBkgy->SetMarkerColor(3);
  hnBkgy->SetLineColor(3);
  NewhnBkgy->Draw("same");
  NewhnBkgy->SetMarkerColor(6);
  NewhnBkgy->SetLineColor(6);
  NewerhnBkgy->Draw("same");
  TLine *nBkgylineLow = new TLine(-1.93,0,1.93,0);
  nBkgylineLow->SetLineColor(kRed);
  nBkgylineLow->Draw();
  TLegend *nBkgLegend = new TLegend(0.4,0.7,0.6,0.9);
  nBkgLegend->AddEntry(hnBkgy,"n=2","p");
  nBkgLegend->AddEntry(NewhnBkgy,"n=3","p");
  nBkgLegend->AddEntry(NewerhnBkgy,"n=4","p");
  nBkgLegend->Draw("same");

  //save plots
  calpha->SaveAs(Form("fitted_alpha_PA_%isbins_change.png",whichUpsilon));
  cf1s->SaveAs(Form("fitted_f1s_PA_%isbins_change.png",whichUpsilon));
  cmass->SaveAs(Form("fitted_mass_PA_%isbins_change.png",whichUpsilon));
  cn1s->SaveAs(Form("fitted_n1s_PA_%isbins_change.png",whichUpsilon));
  csigma1s->SaveAs(Form("fitted_sigma1s_PA_%isbins_change.png",whichUpsilon));
  cx1s->SaveAs(Form("fitted_x1s_PA_%isbins_change.png",whichUpsilon));
  cerrmu->SaveAs(Form("fitted_errmu_PA_%isbins_change.png",whichUpsilon));
  cerrsigma->SaveAs(Form("fitted_errsigma_PA_%isbins_change.png",whichUpsilon));
  clambda->SaveAs(Form("fitted_lambda_PA_%isbins_change.png",whichUpsilon));
  cnSig1s->SaveAs(Form("fitted_nSig1s_PA_%isbins_change.png",whichUpsilon));
  cnSig2s->SaveAs(Form("fitted_nSig2s_PA_%isbins_change.png",whichUpsilon));
  cnSig3s->SaveAs(Form("fitted_nSig3s_PA_%isbins_change.png",whichUpsilon));
  cnBkg->SaveAs(Form("fitted_nBkg_PA_%isbins_change.png",whichUpsilon));

}
