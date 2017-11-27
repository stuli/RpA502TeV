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
void PlotFittedYields() {

  int collId = kPPDATA;
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

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetErrorX(0);

  const float histmin = 0;
  //const float histmax = 1200*whichUpsilon;
  const float histmax = 7000*whichUpsilon;

  //set up canvases
  TCanvas *cyield = new TCanvas("cyield","cyield",4,45,500,400);
  //cmass->Divide(2,1);

  /*float xshift = 0.3;
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
    kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
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
  }

  //draw pt plots
  cmass->cd();
  hPAnommasspt->SetXTitle("pT");
  //hPAnommasspt->GetYaxis()->SetRangeUser(histmin,histmax);
  hPAnommasspt->SetMarkerStyle(4);
  hPAnommasspt->SetMarkerSize(1);
  hPAnommasspt->SetMarkerColor(kOrange);
  hPAnommasspt->Draw("hist p");

  hPAaltmasspt->SetMarkerStyle(25);
  hPAaltmasspt->SetMarkerSize(1);
  hPAaltmasspt->SetMarkerColor(2);
  hPAaltmasspt->Draw("same hist p");

  hPPnommasspt->SetMarkerStyle(26);
  hPPnommasspt->SetMarkerSize(1);
  hPPnommasspt->SetMarkerColor(3);
  hPPnommasspt->Draw("same hist p");

  hPPaltmasspt->SetMarkerStyle(28);
  hPPaltmasspt->SetMarkerSize(1);
  hPPaltmasspt->SetMarkerColor(4);
  hPPaltmasspt->Draw("same hist p");

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
  cmass->SaveAs(Form("fitted_mass_fancy_%isbins.png",whichUpsilon));
*/


  float xshift = 0.05;
  float shiftedybins[yconst];
  for (int ishift=0; ishift<yconst; ishift++) {
    shiftedybins[ishift] = ybins[ishift]-xshift;
  }

  TH1F* hPAnomyieldy = new TH1F("hPAnomyieldy","Fitted Upsilon(1S) Yield",numybins,ybins);
  //TH1F* hPAaltyieldy = new TH1F("hPAaltyieldy","yield vs y",numybins,shiftedybins);

//1S y loop
  for (int iy = 0; iy<numybins; iy++) {

    float ptLow = 0;
    float ptHigh = 30;
    float yLow = ybins[iy];
    float yHigh = ybins[iy+1];
    float yLowBad = yLow-0.47;
    float yHighBad = yHigh-0.47;

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString FileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float PAnomyield = PAnomws->var("nSig1s")->getVal();  
    //float PAnomyielderr = PAnomws->var("nSig1s")->getError();
    float deriv = PAnomyield/(yHigh-yLow);
    PANomFile->Close("R");

   /* //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowBad, yHighBad, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    FileName = Form("BadFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PAAltFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAaltws = (RooWorkspace*)PAAltFile->Get("workspace");

    //extract parameter values
    float PAaltyield = PAaltws->var("nSig1s")->getVal();  
    float PAaltyielderr = PAaltws->var("nSig1s")->getError();
    PAAltFile->Close("R");
*/
    //fill histograms
    hPAnomyieldy->SetBinContent(iy+1, deriv);
    //hPAnomyieldy->SetBinError  (iy+1, PAnomyielderr);
    //hPAaltyieldy->SetBinContent(iy+1, PAaltyield);
    //hPAaltyieldy->SetBinError  (iy+1, PAaltyielderr);
  }

  //draw pt plots
  cyield->cd();
  hPAnomyieldy->SetXTitle("y");
  hPAnomyieldy->GetYaxis()->SetRangeUser(histmin,histmax);
  hPAnomyieldy->SetMarkerStyle(4);
  //hPAnomyieldy->SetMarkerSize(1);
  hPAnomyieldy->SetMarkerColor(kBlue);
  hPAnomyieldy->Draw("hist p");
  //hPAaltyieldy->SetMarkerStyle(25);
  //hPAaltyieldy->SetMarkerSize(1);
  //hPAaltyieldy->SetMarkerColor(kRed);
  //hPAaltyieldy->Draw("same hist p");

  /*for (int ilines = 0; ilines<numybins; ilines++) {
    float xval = ybins[ilines];
    TLine *binline = new TLine(xval,9.42,xval,9.45);
    binline->SetLineStyle(2);
    binline->Draw();
    float xmin = ybins[ilines]+0.04;
    float xmax = ybins[ilines+1]-0.04;
    TArrow *arbin = new TArrow(xmin,9.43,xmax,9.43,0.02,"<>");
    arbin->Draw();
  }*/

  TLegend* yieldyLegend = new TLegend(0.4,0.7,0.6,0.9);
  //yieldyLegend->SetTextSize(16);
  //yieldyLegend->SetTextFont(43);
  yieldyLegend->AddEntry("hPAnomyieldy","pPb nominal","p");
  //yieldyLegend->AddEntry("hPAaltyieldy","pPb alternative","p");
  yieldyLegend->Draw("same");

  //save plots
  cyield->SaveAs(Form("PPfitted_yield_y_fancy_%isbins.png",whichUpsilon));
}
