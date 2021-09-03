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


void PlotFittedParams_intBinSinglePlots(int collId=kPADATA) {

  TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";
  float scale = 6;
  if (collId==kPPDATA) scale = scale*4;
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
  float ptbins[2] = {0,30};

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;

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
  TH1F* halphapt = new TH1F("halphapt","alpha vs pt",numptbins,ptbins);
  TH1F* hf1spt = new TH1F("hf1spt","f1s vs pt",numptbins,ptbins);
  TH1F* hmasspt = new TH1F("hmasspt","mass vs pt",numptbins,ptbins);
  TH1F* hn1spt = new TH1F("hn1spt","n1s vs pt",numptbins,ptbins);
  TH1F* hsigma1spt = new TH1F("hsigma1spt","sigma vs pt",numptbins,ptbins);
  TH1F* hx1spt = new TH1F("hx1spt","x1s vs pt",numptbins,ptbins);
  TH1F* herrmupt = new TH1F("herrmupt","errmu vs pt",numptbins,ptbins);
  TH1F* herrsigmapt = new TH1F("herrsigmapt","errsigma vs pt",numptbins,ptbins);
  TH1F* hlambdapt = new TH1F("hlambdapt","errlambda vs pt",numptbins,ptbins);
  TH1F* hnSig1spt = new TH1F("hnSig1spt","nSig1s vs pt",numptbins,ptbins);
  TH1F* hnSig2spt = new TH1F("hnSig2spt","nSig2s vs pt",numptbins,ptbins);
  TH1F* hnSig3spt = new TH1F("hnSig3spt","nSig3s vs pt",numptbins,ptbins);
  TH1F* hnBkgpt = new TH1F("hnBkgpt","nBkg vs pt",numptbins,ptbins);

  TString kineLabel, NomFileName;
  TFile* NomFile;
  float ptLow, ptHigh, yLow, yHigh;
  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
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

  TString strId;
  if (collId==kPADATA) strId = "PA";
  else if (collId==kPPDATA) strId = "PP";

  //save plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_int.png",strId.Data()));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_int.png",strId.Data()));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_int.png",strId.Data()));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_int.png",strId.Data()));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_int.png",strId.Data()));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_int.png",strId.Data()));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_int.png",strId.Data()));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_int.png",strId.Data()));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_int.png",strId.Data()));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_int.png",strId.Data()));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_int.png",strId.Data()));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_int.png",strId.Data()));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_int.png",strId.Data()));

  //save pdf plots
  calpha->SaveAs(Form("ParameterPlots/%sfitted_alpha_int.pdf",strId.Data()));
  cf1s->SaveAs(Form("ParameterPlots/%sfitted_f1s_int.pdf",strId.Data()));
  cmass->SaveAs(Form("ParameterPlots/%sfitted_mass_int.pdf",strId.Data()));
  cn1s->SaveAs(Form("ParameterPlots/%sfitted_n1s_int.pdf",strId.Data()));
  csigma1s->SaveAs(Form("ParameterPlots/%sfitted_sigma1s_int.pdf",strId.Data()));
  cx1s->SaveAs(Form("ParameterPlots/%sfitted_x1s_int.pdf",strId.Data()));
  cerrmu->SaveAs(Form("ParameterPlots/%sfitted_errmu_int.pdf",strId.Data()));
  cerrsigma->SaveAs(Form("ParameterPlots/%sfitted_errsigma_int.pdf",strId.Data()));
  clambda->SaveAs(Form("ParameterPlots/%sfitted_lambda_int.pdf",strId.Data()));
  cnSig1s->SaveAs(Form("ParameterPlots/%sfitted_nSig1s_int.pdf",strId.Data()));
  cnSig2s->SaveAs(Form("ParameterPlots/%sfitted_nSig2s_int.pdf",strId.Data()));
  cnSig3s->SaveAs(Form("ParameterPlots/%sfitted_nSig3s_int.pdf",strId.Data()));
  cnBkg->SaveAs(Form("ParameterPlots/%sfitted_nBkg_int.pdf",strId.Data()));

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
