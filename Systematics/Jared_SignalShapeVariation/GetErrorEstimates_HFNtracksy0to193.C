#include <iostream>
#include <iomanip>
#include <sstream>
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
void GetErrorEstimates_HFNtracksy0to193(int binmode=0) {
//0=hfmode, 1=ntmode

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  whichUpsilon=3;
  if (whichUpsilon==3) {
    float hfbins[3] = {0,12,120};
    int ntracksbins[3] = {0,40,400};
  }
  else {
    float hfbins[5] = {0,15,22,30,120};
    int ntracksbins[5] = {0,40,65,90,400};
  }

  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntbins = sizeof(ntracksbins)/sizeof(float)-1;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -30;
  int max = 30;

  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TCanvas *c1 = new TCanvas("c1","c1",4,45,400,400);

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA F ERR");
  cout << setw(Width1S) << setfill(separator) << Form("PA B ERR");
  cout << setw(Width1S) << setfill(separator) << "RFB ERR";
  cout << endl;

  float yLowB, yHighB, yLowF, yHighF, hfLow, hfHigh;
  int ntracksLow, ntracksHigh;
  float ptLow = 0;
  float ptHigh = 30;
  TString binLabel, kineLabelF, kineLabelB, histFileNameF, histFileNameB, canvasFileName;
  TFile* theFileF;
  TFile* theFileB;

  TString hfntracksbins = "_hfbins.root";
  int numbins = numhfbins;
  if (binmode==1) {
    hfntracksbins = "_ntracksbins.root";
    numbins = numntbins;
  }

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is",whichUpsilon);
  outFileName = outFileName + hfntracksbins;
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple12 = new TNtuple("ntuple3s","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb1sFerr:pPb1sBerr:RFBErr",8);

  //BIN LOOP********************************************************
  for (int ihf = 0; ihf<numbins; ihf++) {
    //Choose the hf bin
    if (binmode==1) {
      hfLow = 0;
      hfHigh = 120;
      ntracksLow = ntracksbins[ihf];
      ntracksHigh = ntracksbins[ihf+1];
    }
    else {
      hfLow = hfbins[ihf];
      hfHigh = hfbins[ihf+1];
      ntracksLow = 0;
      ntracksHigh = 400;
    }
  for (int whichUpsilon = 1; whichUpsilon<3; whichUpsilon++) {

    //Choose the rapidity bin
    yLowB = -1.93;
    yHighB = 0.00;
    yLowF = 0.00;
    yHighF = 1.93;

    //print bin label
    binLabel = Form("hf%.1f-%.1f_y%.2f-%.2f",hfLow,hfHigh,yLowF,yHighF);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabelF = kineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    histFileNameF = Form("PseudoExperimentResultsFeb2018/PseudoExpResults_%s.root",kineLabelF.Data());
    theFileF = new TFile(histFileNameF);
    TNtuple* ntupleResultsF = (TNtuple*)theFileF->Get("ntuple;1");

    kineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    kineLabelB = kineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    histFileNameB = Form("PseudoExperimentResultsFeb2018/PseudoExpResults_%s.root",kineLabelB.Data());
    theFileB = new TFile(histFileNameB);
    TNtuple* ntupleResultsB = (TNtuple*)theFileB->Get("ntuple;1");

    //Extract errors in yields
    if (whichUpsilon==1) TString strBranch = "diff1s";
    else if (whichUpsilon==2) TString strBranch = "diff2s";
    else if (whichUpsilon==3) TString strBranch = "diff3s";
    if (ntupleResultsF->GetMinimum("diff1s")<min) min = -20;
    else min = -10;
    if (ntupleResultsF->GetMaximum("diff1s")>max) max = 20;
    else max = 10;
    cntuple->cd(1);
    TH1F * ntuplehistoF = new TH1F("ntuplehistoF", "pPb %Diff in Yield Forward", 100,min,max);
    ntuplehistoF->SetXTitle("%Diff");
    ntuplehistoF->GetXaxis()->SetTitleSize(0.05);
    ntuplehistoF->GetYaxis()->SetLabelSize(0.05);
    ntuplehistoF->GetXaxis()->SetLabelSize(0.05);
    ntuplehistoF->GetXaxis()->SetRangeUser(min,max);
    if (whichUpsilon==1) ntupleResultsF->Draw("diff1s>>ntuplehistoF");
    else if (whichUpsilon==2) ntupleResultsF->Draw("diff2s>>ntuplehistoF");
    else if (whichUpsilon==3) ntupleResultsF->Draw("diff3s>>ntuplehistoF");
    float temperrF = ntuplehistoF->GetMean();  
    float temprmsF = ntuplehistoF->GetRMS();
    if (TMath::Abs(temprmsF)>TMath::Abs(temperrF)) temperrF = temprmsF;
    //ntuplehistoF->Fit("gaus","Q");
    //float tempmuF = ntuplehistoF->GetFunction("gaus")->GetParameter(1);//fitted mean
    //if (TMath::Abs(tempmuF)>TMath::Abs(temperrF)) temperrF = tempmuF;
    //float tempsigmaF = ntuplehistoF->GetFunction("gaus")->GetParameter(2);//fitted sigma
    //if (TMath::Abs(tempsigmaF)>TMath::Abs(temperrF)) temperrF = tempsigmaF;

    if (ntupleResultsB->GetMinimum("diff1s")<min) min = -20;
    if (ntupleResultsB->GetMaximum("diff1s")>max) max = 20;
    cntuple->cd(2);
    TH1F * ntuplehistoB = new TH1F("ntuplehistoB", "pPb %Diff in Yield Backward", 100,min,max);
    ntuplehistoB->SetXTitle("%Diff");
    ntuplehistoB->GetXaxis()->SetTitleSize(0.05);
    ntuplehistoB->GetYaxis()->SetLabelSize(0.05);
    ntuplehistoB->GetXaxis()->SetLabelSize(0.05);
    ntuplehistoB->GetXaxis()->SetRangeUser(min,max);
    if (whichUpsilon==1) ntupleResultsB->Draw("diff1s>>ntuplehistoB");
    else if (whichUpsilon==2) ntupleResultsB->Draw("diff2s>>ntuplehistoB");
    else if (whichUpsilon==3) ntupleResultsB->Draw("diff3s>>ntuplehistoB");
    float temperrB = ntuplehistoB->GetMean();  
    float temprmsB = ntuplehistoB->GetRMS();
    if (TMath::Abs(temprmsB)>TMath::Abs(temperrB)) temperrB = temprmsB;
    //ntuplehistoB->Fit("gaus","Q");
    //float tempmuB = ntuplehistoB->GetFunction("gaus")->GetParameter(1);//fitted mean
    //if (TMath::Abs(tempmuB)>TMath::Abs(temperrB)) temperrB = tempmuB;
    //float tempsigmaB = ntuplehistoB->GetFunction("gaus")->GetParameter(2);//fitted sigma
    //if (TMath::Abs(tempsigmaB)>TMath::Abs(temperrB)) temperrB = tempsigmaB;

    //RFB
    const int NEvents = 100;
    TH1F* hRFB = new TH1F("hRFB","hRFB",100,min,max);
    TString yieldNom = Form("yield%isNom",whichUpsilon);
    TString yieldAlt = Form("yield%isAlt",whichUpsilon);
    TLeaf *Fyield1sNomLeaf = ntupleResultsF->GetLeaf(yieldNom);
    TLeaf *Byield1sNomLeaf = ntupleResultsB->GetLeaf(yieldNom);
    TLeaf *Fyield1sAltLeaf = ntupleResultsF->GetLeaf(yieldAlt);
    TLeaf *Byield1sAltLeaf = ntupleResultsB->GetLeaf(yieldAlt);
    for (int i = 0; i<NEvents; i++) {
      ntupleResultsF->GetEntry(i);
      float Fyield1sNom = (float)Fyield1sNomLeaf->GetValue();
      float Fyield1sAlt = (float)Fyield1sAltLeaf->GetValue();
      ntupleResultsB->GetEntry(i);
      float Byield1sNom = (float)Byield1sNomLeaf->GetValue();
      float Byield1sAlt = (float)Byield1sAltLeaf->GetValue();
      float RFBNom = Fyield1sNom/Byield1sNom;
      float RFBAlt = Fyield1sAlt/Byield1sAlt;
      float RFBDiff = (RFBAlt-RFBNom)/RFBNom*100;
      hRFB->Fill(RFBDiff);
    }

    cntuple->cd(3);
    hRFB->SetXTitle("%Diff in RFB");
    hRFB->GetXaxis()->SetTitleSize(0.05);
    hRFB->GetYaxis()->SetLabelSize(0.05);
    hRFB->GetXaxis()->SetLabelSize(0.05);
    hRFB->GetXaxis()->SetRangeUser(min,max);
    hRFB->Draw();
    float RFBerr = hRFB->GetMean();
    float RFBrms = hRFB->GetRMS();
    if (TMath::Abs(RFBrms)>TMath::Abs(RFBerr)) RFBerr = RFBrms;
    //hRFB->Fit("gaus","Q");
    //float RFBmu = hRFB->GetFunction("gaus")->GetParameter(1);//fitted mean
    //if (TMath::Abs(RFBmu)>TMath::Abs(RFBerr)) RFBerr = RFBmu;
    //float RFBsigma = hRFB->GetFunction("gaus")->GetParameter(2);//fitted sigma
    //if (TMath::Abs(RFBsigma)>TMath::Abs(RFBerr)) RFBerr = RFBsigma;

    //Save the plots of yields and RFB errors
    canvasFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f_hfsum%.2f-%.2f_ntracks%i-%i.png",whichUpsilon,ptLow,ptHigh,yLowF,yHighF,muPtCut,hfLow,hfHigh,ntracksLow,ntracksHigh);
    cntuple->SaveAs(canvasFileName);
    canvasFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f_hfsum%.2f-%.2f_ntracks%i-%i.pdf",whichUpsilon,ptLow,ptHigh,yLowF,yHighF,muPtCut,hfLow,hfHigh,ntracksLow,ntracksHigh);
    cntuple->SaveAs(canvasFileName);

    //Take absolute values
    temperrF = TMath::Abs(temperrF);
    temperrB = TMath::Abs(temperrB);
    RFBerr = TMath::Abs(RFBerr);

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << temperrF;
    cout << setw(Width1S) << setfill(separator) << temperrB;
    cout << setw(Width1S) << setfill(separator) << RFBerr;
    cout << endl;

    //put errors in ntuple
    int binlow = (int)hfLow;
    int binhigh = (int)hfHigh;
    if (binmode==1) {
      binlow = ntracksLow;
      binhigh = ntracksHigh;
    }
    ntuple12->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBerr);
    theFileF->Close();
    theFileB->Close();

  }//end of rapidity bin loop
  }//end of hf bin loop

  c1->Close();
  cntuple->Close();

  //Write errors in ntuple file
  outFile.cd();
  ntuple12->Write();
  outFile.Close();

}
