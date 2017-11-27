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

using namespace std;
void PlotFittedMass() {

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  int whichUpsilon = 1;

    float ptbins1[7] = {0,2,4,6,9,12,30};
    float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float ptbins2[4] = {0,4,9,30};
    float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
    float ptbins3[3] = {0,6,30};
    float ybins3[3] = {-1.93,0.0,1.93};

  int ptsize, ysize;

  if (whichUpsilon==1) {
    ptsize = 7;
    ysize = 9;
  }
  else if (whichUpsilon==2) {
    ptsize = 4;
    ysize = 5;
  }
  else if (whichUpsilon==3) {
    ptsize = 3;
    ysize = 3;
  }

  const int ptconst = ptsize;
  const int yconst = ysize;
  const int numptbins = ptsize-1;
  const int numybins = ysize-1;
  float ptbins[ptconst];
  float ybins[yconst];
  if (whichUpsilon==1) {
    memcpy(ptbins,ptbins1,sizeof(ptbins));
    memcpy(ybins,ybins1,sizeof(ybins));
  }
  else if (whichUpsilon==2) {
    memcpy(ptbins,ptbins2,sizeof(ptbins));
    memcpy(ybins,ybins2,sizeof(ybins));
  }
  else if (whichUpsilon==3) {
    memcpy(ptbins,ptbins3,sizeof(ptbins));
    memcpy(ybins,ybins3,sizeof(ybins));
  }

  float histmin = 9.42;
  float histmax = 9.5;

  TString kineLabel;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetErrorX(0);

  //set up canvases
  TCanvas *cmass = new TCanvas("cmass","cmass",4,45,500,400);
  //cmass->Divide(2,1);
/*
  float xshift = 0.3;
  float shiftedptbins[ptconst];
  float shiftedptbins2[ptconst];
  float shiftedptbins3[ptconst];
  for (int ishift=0; ishift<ptconst; ishift++) {
    shiftedptbins[ishift] = ptbins[ishift]-xshift;
    shiftedptbins2[ishift] = ptbins[ishift]+xshift;
    shiftedptbins3[ishift] = ptbins[ishift]+2*xshift;
  }

  //declare histograms
  TH1F* hPAnommasspt = new TH1F("hPAnommasspt","Fitted Upsilon(1S) Mass",numptbins,ptbins);
  TH1F* hPAaltmasspt = new TH1F("hPAaltmasspt","mass vs pt",numptbins,shiftedptbins3);
  TH1F* hPPnommasspt = new TH1F("hPPnommasspt","mass vs pt",numptbins,shiftedptbins);
  TH1F* hPPaltmasspt = new TH1F("hPPaltmasspt","mass vs pt",numptbins,shiftedptbins2);

  //pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    float ptLow = ptbins[ipt];
    float ptHigh = ptbins[ipt+1];
    float yLow = -1.93;
    float yHigh = 1.93;

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString FileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float PAnommass = PAnomws->var("m_{#Upsilon(1S)}")->getVal();  
    float PAnommasserr = PAnomws->var("m_{#Upsilon(1S)}")->getError();
    PANomFile->Close("R");

    //import fitted model
    FileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PAAltFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAaltws = (RooWorkspace*)PAAltFile->Get("workspace");

    //extract parameter values
    float PAaltmass = PAaltws->var("m_{#Upsilon(1S)}")->getVal();  
    float PAaltmasserr = PAaltws->var("m_{#Upsilon(1S)}")->getError();
    PAAltFile->Close("R");

    //import fitted model
    kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, 0.0, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    FileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PPNomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PPnomws = (RooWorkspace*)PPNomFile->Get("workspace");

    //extract parameter values
    float PPnommass = PPnomws->var("m_{#Upsilon(1S)}")->getVal();  
    float PPnommasserr = PPnomws->var("m_{#Upsilon(1S)}")->getError();
    PPNomFile->Close("R");

    //import fitted model
    FileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PPAltFile = TFile::Open(FileName,"READ");
    RooWorkspace *PPaltws = (RooWorkspace*)PPAltFile->Get("workspace");

    //extract parameter values
    float PPaltmass = PPaltws->var("m_{#Upsilon(1S)}")->getVal();  
    float PPaltmasserr = PPaltws->var("m_{#Upsilon(1S)}")->getError();
    PPAltFile->Close("R");

    //fill histograms
    hPAnommasspt->SetBinContent(ipt+1, PAnommass);
    hPAnommasspt->SetBinError  (ipt+1, PAnommasserr);
    hPAaltmasspt->SetBinContent(ipt+1, PAaltmass);
    hPAaltmasspt->SetBinError  (ipt+1, PAaltmasserr);
    hPPnommasspt->SetBinContent(ipt+1, PPnommass);
    hPPnommasspt->SetBinError  (ipt+1, PPnommasserr);
    hPPaltmasspt->SetBinContent(ipt+1, PPaltmass);
    hPPaltmasspt->SetBinError  (ipt+1, PPaltmasserr);

    cout << "PAnommasserr = " << PAnommasserr << endl;
    cout << "PAaltmasserr = " << PAaltmasserr << endl;
    cout << "PPnommasserr = " << PPnommasserr << endl;
    cout << "PPaltmasserr = " << PPaltmasserr << endl;
  }

  //draw pt plots
  cmass->cd();
  hPAnommasspt->SetXTitle("pT");
  hPAnommasspt->GetYaxis()->SetRangeUser(histmin,histmax);
  hPAnommasspt->SetMarkerStyle(4);
  hPAnommasspt->SetMarkerSize(1);
  hPAnommasspt->SetMarkerColor(kOrange);
  hPAnommasspt->Draw();

  hPAaltmasspt->SetMarkerStyle(25);
  hPAaltmasspt->SetMarkerSize(1);
  hPAaltmasspt->SetMarkerColor(2);
  hPAaltmasspt->Draw("same");

  hPPnommasspt->SetMarkerStyle(26);
  hPPnommasspt->SetMarkerSize(1);
  hPPnommasspt->SetMarkerColor(3);
  hPPnommasspt->Draw("same");

  hPPaltmasspt->SetMarkerStyle(28);
  hPPaltmasspt->SetMarkerSize(1);
  hPPaltmasspt->SetMarkerColor(4);
  hPPaltmasspt->Draw("same");

  //TLine *massptline = new TLine(0,9.46,30,9.46);
  //massptline->SetLineColor(3);
  //massptline->Draw("same");

  TF1 *pdgmass = new TF1("pdgmass","9.46",0,30);
  pdgmass->Draw("same");
  pdgmass->SetLineColor(7);
  pdgmass->SetLineStyle(3);
  pdgmass->SetLineWidth(0.5);

  for (int ilines = 0; ilines<numptbins; ilines++) {
    float xval = ptbins[ilines];
    TLine *binline = new TLine(xval,9.42,xval,9.45);
    binline->SetLineStyle(2);
    binline->Draw();
    float xmin = ptbins[ilines]+0.1;
    float xmax = ptbins[ilines+1]-0.1;
    TArrow *arbin = new TArrow(xmin,9.43,xmax,9.43,0.02,"<>");
    arbin->Draw();
  }

  TLegend* massptLegend = new TLegend(0.5,0.6,0.9,0.9);
  massptLegend->SetTextSize(16);
  massptLegend->SetTextFont(43);
  massptLegend->AddEntry("hPAnommasspt","pPb nominal","p");
  massptLegend->AddEntry("hPAaltmasspt","pPb alternative","p");
  massptLegend->AddEntry("hPPnommasspt","pp nominal","p");
  massptLegend->AddEntry("hPPaltmasspt","pp alternative","p");
  massptLegend->AddEntry("pdgmass","pdg mass","l");
  massptLegend->Draw("same");

  //save plots
  cmass->SaveAs(Form("fitted_mass_fancy_%isbins_withErrorBars.png",whichUpsilon));

*/

  float xshift = 0.05;
  float shiftedybins[yconst];
  float shiftedybins2[yconst];
  float shiftedybins3[yconst];
  for (int ishift=0; ishift<yconst; ishift++) {
    shiftedybins[ishift] = ybins[ishift]-xshift;
    shiftedybins2[ishift] = ybins[ishift]+xshift;
    shiftedybins3[ishift] = ybins[ishift]+2*xshift;
  }

  TH1F* hPAnommassy = new TH1F("hPAnommassy","Fitted Upsilon(1S) Mass",numybins,ybins);
  TH1F* hPAaltmassy = new TH1F("hPAaltmassy","mass vs y",numybins,shiftedybins3);
  TH1F* hPPnommassy = new TH1F("hPPnommassy","mass vs y",numybins,shiftedybins);
  TH1F* hPPaltmassy = new TH1F("hPPaltmassy","mass vs y",numybins,shiftedybins2);

//1S y loop
  for (int iy = 0; iy<numybins; iy++) {

    float ptLow = 0;
    float ptHigh = 30;
    float yLow = ybins[iy];
    float yHigh = ybins[iy+1];
    float yLowPP, yHighPP;
    if (yLow<0) {
      yLowPP = TMath::Abs(yHigh);
      yHighPP = TMath::Abs(yLow);
    }
    else {
      yLowPP = yLow;
      yHighPP = yHigh;
    }

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString FileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float PAnommass = PAnomws->var("m_{#Upsilon(1S)}")->getVal();  
    float PAnommasserr = PAnomws->var("m_{#Upsilon(1S)}")->getError();
    PANomFile->Close("R");

    //import fitted model
    FileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PAAltFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAaltws = (RooWorkspace*)PAAltFile->Get("workspace");

    //extract parameter values
    float PAaltmass = PAaltws->var("m_{#Upsilon(1S)}")->getVal();  
    float PAaltmasserr = PAaltws->var("m_{#Upsilon(1S)}")->getError();
    PAAltFile->Close("R");

    //import fitted model
    kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    FileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PPNomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PPnomws = (RooWorkspace*)PPNomFile->Get("workspace");

    //extract parameter values
    float PPnommass = PPnomws->var("m_{#Upsilon(1S)}")->getVal();  
    float PPnommasserr = PPnomws->var("m_{#Upsilon(1S)}")->getError();
    PPNomFile->Close("R");

    //import fitted model
    FileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PPAltFile = TFile::Open(FileName,"READ");
    RooWorkspace *PPaltws = (RooWorkspace*)PPAltFile->Get("workspace");

    //extract parameter values
    float PPaltmass = PPaltws->var("m_{#Upsilon(1S)}")->getVal();  
    float PPaltmasserr = PPaltws->var("m_{#Upsilon(1S)}")->getError();
    PPAltFile->Close("R");

    //fill histograms
    hPAnommassy->SetBinContent(iy+1, PAnommass);
    hPAnommassy->SetBinError  (iy+1, PAnommasserr);
    hPAaltmassy->SetBinContent(iy+1, PAaltmass);
    hPAaltmassy->SetBinError  (iy+1, PAaltmasserr);
    hPPnommassy->SetBinContent(iy+1, PPnommass);
    hPPnommassy->SetBinError  (iy+1, PPnommasserr);
    hPPaltmassy->SetBinContent(iy+1, PPaltmass);
    hPPaltmassy->SetBinError  (iy+1, PPaltmasserr);

    cout << "PAnommasserr = " << PAnommasserr << endl;
    cout << "PAaltmasserr = " << PAaltmasserr << endl;
    cout << "PPnommasserr = " << PPnommasserr << endl;
    cout << "PPaltmasserr = " << PPaltmasserr << endl;
  }

  //draw pt plots
  cmass->cd();
  hPAnommassy->SetXTitle("y");
  hPAnommassy->GetYaxis()->SetRangeUser(histmin,histmax);
  hPAnommassy->SetMarkerStyle(4);
  hPAnommassy->SetMarkerSize(1);
  hPAnommassy->SetMarkerColor(kOrange);
  hPAnommassy->Draw();
  hPAaltmassy->SetMarkerStyle(25);
  hPAaltmassy->SetMarkerSize(1);
  hPAaltmassy->SetMarkerColor(2);
  hPAaltmassy->Draw("same");

  hPPnommassy->SetMarkerStyle(26);
  hPPnommassy->SetMarkerSize(1);
  hPPnommassy->SetMarkerColor(3);
  hPPnommassy->Draw("same");

  hPPaltmassy->SetMarkerStyle(28);
  hPPaltmassy->SetMarkerSize(1);
  hPPaltmassy->SetMarkerColor(4);
  hPPaltmassy->Draw("same");

  TF1 *pdgmass = new TF1("pdgmass","9.46",-1.93,1.93);
  pdgmass->Draw("same");
  pdgmass->SetLineColor(7);
  pdgmass->SetLineStyle(3);
  pdgmass->SetLineWidth(0.5);

  for (int ilines = 0; ilines<numybins; ilines++) {
    float xval = ybins[ilines];
    TLine *binline = new TLine(xval,9.42,xval,9.45);
    binline->SetLineStyle(2);
    binline->Draw();
    float xmin = ybins[ilines]+0.04;
    float xmax = ybins[ilines+1]-0.04;
    TArrow *arbin = new TArrow(xmin,9.43,xmax,9.43,0.02,"<>");
    arbin->Draw();
  }

  TLegend* massyLegend = new TLegend(0.5,0.6,0.9,0.9);
  massyLegend->SetTextSize(16);
  massyLegend->SetTextFont(43);
  massyLegend->AddEntry("hPAnommassy","pPb nominal","p");
  massyLegend->AddEntry("hPAaltmassy","pPb alternative","p");
  massyLegend->AddEntry("hPPnommassy","pp nominal","p");
  massyLegend->AddEntry("hPPaltmassy","pp alternative","p");
  massyLegend->AddEntry("pdgmass","pdg mass","l");
  massyLegend->Draw("same");

  //save plots
  cmass->SaveAs(Form("fitted_mass_y_fancy_%isbins.png",whichUpsilon));

}
