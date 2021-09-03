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
void PlotNtracksYields(int whichUpsilon=1) {

  int collId = kPADATA;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = -1.93;
  float yHigh = 1.93;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;
  float hfLow = 0;
  float hfHigh = 400;
  float ntracksLow, ntracksHigh;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //if (whichUpsilon==2) yLow=0.8;
  //if (whichUpsilon==3) yLow=0.0;

  float NtracksBins[5] = {0,40,62,88,400};

  const int numNtracksBins = sizeof(NtracksBins)/sizeof(float)-1;

  //set up plot styles
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetMarkerStyle(2);
  //gStyle->SetErrorX(0);

  //set up canvases
  TCanvas *cyield = new TCanvas("cyield","cyield",4,45,500,400);
  //cmass->Divide(2,1);

  TH1F* hyieldF = new TH1F("hyieldF",Form("Fitted Upsilon Yield"),numNtracksBins,NtracksBins);
  //TH1F* hyieldB = new TH1F("hyieldB",Form("Fitted Upsilon(%iS) Yield",whichUpsilon),numNtracksBins,NtracksBins);

//Ntracks loop
  for (int i = 0; i<numNtracksBins; i++) {

    ntracksLow = NtracksBins[i];
    ntracksHigh = NtracksBins[i+1];

    //import fitted model
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, (int)ntracksLow, (int)ntracksHigh );
    TString FileName = Form("HFTestBins/nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << FileName << endl;
    TFile* PANomFile = TFile::Open(FileName,"READ");
    RooWorkspace *PAnomws = (RooWorkspace*)PANomFile->Get("workspace");

    //extract parameter values
    float yield1 = PAnomws->var(Form("nSig%is",1))->getVal();  
    float yield1err = PAnomws->var(Form("nSig%is",1))->getError();
    float yield2 = PAnomws->var(Form("nSig%is",2))->getVal();  
    float yield2err = PAnomws->var(Form("nSig%is",2))->getError();
    float yield3 = PAnomws->var(Form("nSig%is",3))->getVal();  
    float yield3err = PAnomws->var(Form("nSig%is",3))->getError();
    float yield = yield1+yield2+yield3;  
    float yielderr = yield1err+yield2err+yield3err;
    PANomFile->Close("R");

    //fill histograms
    hyieldF->SetBinContent(i+1, yield);
    hyieldF->SetBinError(i+1, yielderr);
  }

  //draw plot
  cyield->cd();
  hyieldF->SetMinimum(0);
  hyieldF->SetXTitle("Ntracks");
  //hyieldF->GetYaxis()->SetRangeUser(histmin,histmax);
  //hyieldF->SetMarkerStyle(4);
  //hyieldF->SetMarkerSize(1);
  //hyieldF->SetMarkerColor(kBlue);
  hyieldF->Draw();
  /*hyieldB->SetLineColor(3);
  hyieldB->Draw("same");

  TLegend* yieldLegend = new TLegend(0.7,0.1,0.9,0.3);
  yieldLegend->SetTextSize(16);
  yieldLegend->SetTextFont(43);
  yieldLegend->AddEntry("hyieldF","Forward","le");
  yieldLegend->AddEntry("hyieldB","Backward","le");
  yieldLegend->Draw("same");
*/
  //save plots
  TString BinsString = Form("_Bins%i",(int)NtracksBins[0]);
  for(int j = 1; j<=numNtracksBins; j++) BinsString = BinsString + Form("-%i",(int)NtracksBins[j]);
  TString PicName = Form("NtracksYields_integrated");
  PicName = PicName + BinsString + Form(".png");
  cyield->SaveAs(PicName);
}
