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


void PlotFittedParamsFromHistos() {

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.2, 5.0, 3.321, 5.0, 1.0, 25.0, 25.0, 25.0};
  double paramslower[8] = {0.02, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  //choose a set of bins
  int whichUpsilon = 1;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.04,"Y");
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
  HistoFileName = Form("NominalFitsAsOfJan012018/FittedParamsHistos_PA_%is_new.root",whichUpsilon);
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

  NewHistoFileName = Form("FittedParamsHistos_PA_%is_fixedsigman.root",whichUpsilon);
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

  //draw pt plots
  calpha->cd(1);
  halphapt->SetXTitle("pT");
  halphapt->GetYaxis()->SetRangeUser(0,5);
  halphapt->Draw();
  halphapt->SetLineColor(3);
  Newhalphapt->Draw("same");
  //TLine *alphaptlineLow = new TLine(0,paramslower[2],30,paramslower[2]);
  //alphaptlineLow->SetLineColor(kRed);
  //alphaptlineLow->Draw();
  //TLine *alphaptlineHigh = new TLine(0,paramsupper[2],30,paramsupper[2]);
  //alphaptlineHigh->SetLineColor(kRed);
  //alphaptlineHigh->Draw();

  cf1s->cd(1);
  hf1spt->SetXTitle("pT");
  hf1spt->GetYaxis()->SetRangeUser(-1,2);
  hf1spt->Draw();
  hf1spt->SetLineColor(3);
  Newhf1spt->Draw("same");
  //TLine *f1sptlineLow = new TLine(0,paramslower[4],30,paramslower[4]);
  //f1sptlineLow->SetLineColor(kRed);
  //f1sptlineLow->Draw();
  //TLine *f1sptlineHigh = new TLine(0,paramsupper[4],30,paramsupper[4]);
  //f1sptlineHigh->SetLineColor(kRed);
  //f1sptlineHigh->Draw();

  cmass->cd(1);
  hmasspt->SetXTitle("pT");
  hmasspt->GetYaxis()->SetRangeUser(9.35,9.6);
  hmasspt->Draw();
  hmasspt->SetLineColor(3);
  Newhmasspt->Draw("same");
  //TLine *massptlineLow = new TLine(0,9.36,30,9.36);
  //massptlineLow->SetLineColor(kRed);
  //massptlineLow->Draw();
  //TLine *massptlineHigh = new TLine(0,9.56,30,9.56);
  //massptlineHigh->SetLineColor(kRed);
  //massptlineHigh->Draw();

  cn1s->cd(1);
  hn1spt->SetXTitle("pT");
  hn1spt->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1spt->Draw();
  hn1spt->SetLineColor(3);
  Newhn1spt->Draw("same");
  //TLine *n1sptlineLow = new TLine(0,paramslower[3],30,paramslower[3]);
  //n1sptlineLow->SetLineColor(kRed);
  //n1sptlineLow->Draw();
  //TLine *n1sptlineHigh = new TLine(0,paramsupper[3],30,paramsupper[3]);
  //n1sptlineHigh->SetLineColor(kRed);
  //n1sptlineHigh->Draw();

  csigma1s->cd(1);
  hsigma1spt->SetXTitle("pT");
  hsigma1spt->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1spt->Draw();
  hsigma1spt->SetLineColor(3);
  Newhsigma1spt->Draw("same");
  //TLine *sigmaptlineLow = new TLine(0,paramslower[0],30,paramslower[0]);
  //sigmaptlineLow->SetLineColor(kRed);
  //sigmaptlineLow->Draw();
  //TLine *sigmaptlineHigh = new TLine(0,paramsupper[0],30,paramsupper[0]);
  //sigmaptlineHigh->SetLineColor(kRed);
  //sigmaptlineHigh->Draw();

  cx1s->cd(1);
  hx1spt->SetXTitle("pT");
  hx1spt->GetYaxis()->SetRangeUser(0.0,8);
  hx1spt->Draw();
  hx1spt->SetLineColor(3);
  Newhx1spt->Draw("same");
  //TLine *x1sptlineLow = new TLine(0,paramslower[1],30,paramslower[1]);
  //x1sptlineLow->SetLineColor(kRed);
  //x1sptlineLow->Draw();
  //TLine *x1sptlineHigh = new TLine(0,paramsupper[1],30,paramsupper[1]);
  //x1sptlineHigh->SetLineColor(kRed);
  //x1sptlineHigh->Draw();

  cerrmu->cd(1);
  herrmupt->SetXTitle("pT");
  herrmupt->GetYaxis()->SetRangeUser(0.0,25);
  herrmupt->Draw();
  herrmupt->SetLineColor(3);
  Newherrmupt->Draw("same");
  //TLine *errmuptlineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  //errmuptlineLow->SetLineColor(kRed);
  //errmuptlineLow->Draw();
  //TLine *errmuptlineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  //errmuptlineHigh->SetLineColor(kRed);
  //errmuptlineHigh->Draw();

  cerrsigma->cd(1);
  herrsigmapt->SetXTitle("pT");
  herrsigmapt->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmapt->Draw();
  herrsigmapt->SetLineColor(3);
  Newherrsigmapt->Draw("same");
  //TLine *errsigmaptlineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  //errsigmaptlineLow->SetLineColor(kRed);
  //errsigmaptlineLow->Draw();
  //TLine *errsigmaptlineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  //errsigmaptlineHigh->SetLineColor(kRed);
  //errsigmaptlineHigh->Draw();

  clambda->cd(1);
  hlambdapt->SetXTitle("pT");
  hlambdapt->GetYaxis()->SetRangeUser(0.0,25);
  hlambdapt->Draw();
  hlambdapt->SetLineColor(3);
  Newhlambdapt->Draw("same");
  //TLine *lambdaptlineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  //lambdaptlineLow->SetLineColor(kRed);
  //lambdaptlineLow->Draw();
  //TLine *lambdaptlineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  //lambdaptlineHigh->SetLineColor(kRed);
  //lambdaptlineHigh->Draw();

  cnSig1s->cd(1);
  hnSig1spt->SetXTitle("pT");
  hnSig1spt->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig1spt->Draw();
  hnSig1spt->SetLineColor(3);
  NewhnSig1spt->Draw("same");
  //TLine *nSig1sptlineLow = new TLine(0,0,30,0);
  //nSig1sptlineLow->SetLineColor(kRed);
  //nSig1sptlineLow->Draw();

  cnSig2s->cd(1);
  hnSig2spt->SetXTitle("pT");
  hnSig2spt->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig2spt->Draw();
  hnSig2spt->SetLineColor(3);
  NewhnSig2spt->Draw("same");
  //TLine *nSig2sptlineLow = new TLine(0,0,30,0);
  //nSig2sptlineLow->SetLineColor(kRed);
  //nSig2sptlineLow->Draw();

  cnSig3s->cd(1);
  hnSig3spt->SetXTitle("pT");
  hnSig3spt->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig3spt->Draw();
  hnSig3spt->SetLineColor(3);
  NewhnSig3spt->Draw("same");
  //TLine *nSig3sptlineLow = new TLine(0,0,30,0);
  //nSig3sptlineLow->SetLineColor(kRed);
  //nSig3sptlineLow->Draw();

  cnBkg->cd(1);
  hnBkgpt->SetXTitle("pT");
  hnBkgpt->GetYaxis()->SetRangeUser(10,10000*whichUpsilon);
  hnBkgpt->Draw();
  hnBkgpt->SetLineColor(3);
  NewhnBkgpt->Draw("same");
  //TLine *nBkgptlineLow = new TLine(0,0,30,0);
  //nBkgptlineLow->SetLineColor(kRed);
  //nBkgptlineLow->Draw();


  //draw rapidity plots
  calpha->cd(2);
  halphay->SetXTitle("y");
  halphay->GetYaxis()->SetRangeUser(0,5);
  halphay->Draw();
  halphay->SetLineColor(3);
  Newhalphay->Draw("same");
  //TLine *alphaylineLow = new TLine(-1.93,paramslower[2],1.93,paramslower[2]);
  //alphaylineLow->SetLineColor(kRed);
  //alphaylineLow->Draw();
  //TLine *alphaylineHigh = new TLine(-1.93,paramsupper[2],1.93,paramsupper[2]);
  //alphaylineHigh->SetLineColor(kRed);
  //alphaylineHigh->Draw();

  cf1s->cd(2);
  hf1sy->SetXTitle("y");
  hf1sy->GetYaxis()->SetRangeUser(-1,2);
  hf1sy->Draw();
  hf1sy->SetLineColor(3);
  Newhf1sy->Draw("same");
  //TLine *f1sylineLow = new TLine(-1.93,paramslower[4],1.93,paramslower[4]);
  //f1sylineLow->SetLineColor(kRed);
  //f1sylineLow->Draw();
  //TLine *f1sylineHigh = new TLine(-1.93,paramsupper[4],1.93,paramsupper[4]);
  //f1sylineHigh->SetLineColor(kRed);
  //f1sylineHigh->Draw();

  cmass->cd(2);
  hmassy->SetXTitle("y");
  hmassy->GetYaxis()->SetRangeUser(9.35,9.6);
  hmassy->Draw();
  hmassy->SetLineColor(3);
  Newhmassy->Draw("same");
  //TLine *massylineLow = new TLine(-1.93,9.36,1.93,9.36);
  //massylineLow->SetLineColor(kRed);
  //massylineLow->Draw();
  //TLine *massylineHigh = new TLine(-1.93,9.56,1.93,9.56);
  //massylineHigh->SetLineColor(kRed);
  //massylineHigh->Draw();

  cn1s->cd(2);
  hn1sy->SetXTitle("y");
  hn1sy->GetYaxis()->SetRangeUser(-0.2,5.2);
  hn1sy->Draw();
  hn1sy->SetLineColor(3);
  Newhn1sy->Draw("same");
  //TLine *n1sylineLow = new TLine(-1.93,paramslower[3],1.93,paramslower[3]);
  //n1sylineLow->SetLineColor(kRed);
  //n1sylineLow->Draw();
  //TLine *n1sylineHigh = new TLine(-1.93,paramsupper[3],1.93,paramsupper[3]);
  //n1sylineHigh->SetLineColor(kRed);
  //n1sylineHigh->Draw();

  csigma1s->cd(2);
  hsigma1sy->SetXTitle("y");
  hsigma1sy->GetYaxis()->SetRangeUser(-0.1,0.6);
  hsigma1sy->Draw();
  hsigma1sy->SetLineColor(3);
  Newhsigma1sy->Draw("same");
  //TLine *sigmaylineLow = new TLine(-1.93,paramslower[0],1.93,paramslower[0]);
  //sigmaylineLow->SetLineColor(kRed);
  //sigmaylineLow->Draw();
  //TLine *sigmaylineHigh = new TLine(-1.93,paramsupper[0],1.93,paramsupper[0]);
  //sigmaylineHigh->SetLineColor(kRed);
  //sigmaylineHigh->Draw();

  cx1s->cd(2);
  hx1sy->SetXTitle("y");
  hx1sy->GetYaxis()->SetRangeUser(0.0,8);
  hx1sy->Draw();
  hx1sy->SetLineColor(3);
  Newhx1sy->Draw("same");
  //TLine *x1sylineLow = new TLine(-1.93,paramslower[1],1.93,paramslower[1]);
  //x1sylineLow->SetLineColor(kRed);
  //x1sylineLow->Draw();
  //TLine *x1sylineHigh = new TLine(-1.93,paramsupper[1],1.93,paramsupper[1]);
  //x1sylineHigh->SetLineColor(kRed);
  //x1sylineHigh->Draw();

  cerrmu->cd(2);
  herrmuy->SetXTitle("y");
  herrmuy->GetYaxis()->SetRangeUser(0.0,25);
  herrmuy->Draw();
  herrmuy->SetLineColor(3);
  Newherrmuy->Draw("same");
  //TLine *errmuylineLow = new TLine(0,paramslower[6],30,paramslower[6]);
  //errmuylineLow->SetLineColor(kRed);
  //errmuylineLow->Draw();
  //TLine *errmuylineHigh = new TLine(0,paramsupper[6],30,paramsupper[6]);
  //errmuylineHigh->SetLineColor(kRed);
  //errmuylineHigh->Draw();

  cerrsigma->cd(2);
  herrsigmay->SetXTitle("y");
  herrsigmay->GetYaxis()->SetRangeUser(0.0,25);
  herrsigmay->Draw();
  herrsigmay->SetLineColor(3);
  Newherrsigmay->Draw("same");
  //TLine *errsigmaylineLow = new TLine(0,paramslower[5],30,paramslower[5]);
  //errsigmaylineLow->SetLineColor(kRed);
  //errsigmaylineLow->Draw();
  //TLine *errsigmaylineHigh = new TLine(0,paramsupper[5],30,paramsupper[5]);
  //errsigmaylineHigh->SetLineColor(kRed);
  //errsigmaylineHigh->Draw();

  clambda->cd(2);
  hlambday->SetXTitle("y");
  hlambday->GetYaxis()->SetRangeUser(0.0,25);
  hlambday->Draw();
  hlambday->SetLineColor(3);
  Newhlambday->Draw("same");
  //TLine *lambdaylineLow = new TLine(0,paramslower[7],30,paramslower[7]);
  //lambdaylineLow->SetLineColor(kRed);
  //lambdaylineLow->Draw();
  //TLine *lambdaylineHigh = new TLine(0,paramsupper[7],30,paramsupper[7]);
  //lambdaylineHigh->SetLineColor(kRed);
  //lambdaylineHigh->Draw();

  cnSig1s->cd(2);
  hnSig1sy->SetXTitle("y");
  hnSig1sy->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig1sy->Draw();
  hnSig1sy->SetLineColor(3);
  NewhnSig1sy->Draw("same");
  //TLine *nSig1sylineLow = new TLine(-1.93,0,1.93,0);
  //nSig1sylineLow->SetLineColor(kRed);
  //nSig1sylineLow->Draw();

  cnSig2s->cd(2);
  hnSig2sy->SetXTitle("y");
  hnSig2sy->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig2sy->Draw();
  hnSig2sy->SetLineColor(3);
  NewhnSig2sy->Draw("same");
  //TLine *nSig2sylineLow = new TLine(-1.93,0,1.93,0);
  //nSig2sylineLow->SetLineColor(kRed);
  //nSig2sylineLow->Draw();

  cnSig3s->cd(2);
  hnSig3sy->SetXTitle("y");
  hnSig3sy->GetYaxis()->SetRangeUser(0,1600*whichUpsilon);
  hnSig3sy->Draw();
  hnSig3sy->SetLineColor(3);
  NewhnSig3sy->Draw("same");
  //TLine *nSig3sylineLow = new TLine(-1.93,0,1.93,0);
  //nSig3sylineLow->SetLineColor(kRed);
  //nSig3sylineLow->Draw();

  cnBkg->cd(2);
  hnBkgy->SetXTitle("y");
  hnBkgy->GetYaxis()->SetRangeUser(10,10000*whichUpsilon);
  hnBkgy->Draw();
  hnBkgy->SetLineColor(3);
  NewhnBkgy->Draw("same");
  //TLine *nBkgylineLow = new TLine(-1.93,0,1.93,0);
  //nBkgylineLow->SetLineColor(kRed);
  //nBkgylineLow->Draw();

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
