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


void MakeParamsHistosPtBins(int whichUpsilon=1, int whichYRange=0, int collId=kPADATA) {

  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;

  TString kineLabel, NomFileName, outfilename;
  TFile* NomFile;
  float ptLow, ptHigh, yLow, yHigh;

  if (collId==kPADATA) {
    if (whichYRange==0) {
      yLow = -1.93;
      yHigh = 1.93;
    }
    else if (whichYRange==1) {
      yLow = -1.93;
      yHigh = 0.0;
    }
    else if (whichYRange==2) {
      yLow = 0.0;
      yHigh = 1.93;
    }
    outfilename = Form("FittedParamsHistos_PA_%is_PtBins_YRange%i.root",whichUpsilon,whichYRange);
  }
  else if (collId==kPPDATA) {
    yLow = 0.0;
    yHigh = 1.93;
    outfilename = Form("FittedParamsHistos_PP_%is_PtBins_YRange%i.root",whichUpsilon,whichYRange);
  }

  TFile outFile(outfilename, "RECREATE");

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

  float tempalpha, tempalphaerr, tempf1s, tempf1serr, tempmass, tempmasserr, tempn1s, tempn1serr, tempsigma1s, tempsigma1serr, tempx1s, tempx1serr, temperrmu, temperrmuerr, temperrsigma, temperrsigmaerr, templambda, templambdaerr, tempnSig1s, tempnSig1serr, tempnSig2s, tempnSig2serr, tempnSig3s, tempnSig3serr, tempnBkg, tempnBkgerr;

  //1S pt loop
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];

    //import fitted model
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    NomFileName = Form("OfficialNominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
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

  //save histograms
  //TCanvas* c1 = new TCanvas();
  halphapt->Draw();
  //c1->SaveAs(Form("ParameterPlots/halphapt_YRange%i_%iS.pdf",whichYRange,whichUpsilon));
  /*hf1spt->Draw();
  hmasspt->Draw();
  hn1spt->Draw();
  hsigma1spt->Draw();
  hx1spt->Draw();
  herrmupt->Draw();
  herrsigmapt->Draw();
  hlambdapt->Draw();
  hnSig1spt->Draw();
  hnSig2spt->Draw();
  hnSig3spt->Draw();
  hnBkgpt->Draw();

  outFile.cd();

  halphapt->Write();
  hf1spt->Write();
  hmasspt->Write();
  hn1spt->Write();
  hsigma1spt->Write();
  hx1spt->Write();
  herrmupt->Write();
  herrsigmapt->Write();
  hlambdapt->Write();
  hnSig1spt->Write();
  hnSig2spt->Write();
  hnSig3spt->Write();
  hnBkgpt->Write();

  outFile.Close();*/

}
