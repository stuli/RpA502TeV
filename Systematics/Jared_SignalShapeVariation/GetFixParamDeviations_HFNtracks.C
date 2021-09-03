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
void GetFixParamDeviations_HFNtracks(int whichUpsilon=1,int binmode=0, TString whichParam="alpha") {
//0=hfmode, 1=ntmode

  TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  float ybinsCM[3] = {-1.93,0.0,1.93};
  if (whichUpsilon==3) {
    float hfbins[3] = {0,12,120};
    int ntracksbins[3] = {0,40,400};
  }
  else {
    float hfbins[5] = {0,12,19,27,120};
    int ntracksbins[5] = {0,40,62,88,400};
  }
  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntbins = sizeof(ntracksbins)/sizeof(float)-1;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -30;
  int max = 30;

  TString whichBins = "HF";
  if (binmode==1) whichBins = "NTRACKS";
  cout << whichUpsilon << "S " << whichBins << " BINS IN FORWARD AND BACKWARD Y FOR " << whichParam.Data() << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA F %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA B %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RFB ERR";
  cout << endl;

  float yLowB, yHighB, yLowF, yHighF, hfLow, hfHigh;
  int ntracksLow, ntracksHigh;
  float ptLow = 0;
  float ptHigh = 30;
  TString binLabel, kineLabelF, kineLabelB, histFileNameF, histFileNameB;
  TFile* theFileF;
  TFile* theFileB;

  TString hfntracksbins = "_hfbins.root";
  int numbins = numhfbins;
  if (binmode==1) {
    hfntracksbins = "_ntracksbins.root";
    numbins = numntbins;
  }

  TString outFileName = Form("FixParamErrors/Param%sDev%is",whichParam.Data(),whichUpsilon);
  outFileName = outFileName + hfntracksbins;
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1s = new TNtuple("ntuple1s","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb1sFerr:pPb1sBerr:RFBErr",numybins);
  TNtuple* ntuple2s = new TNtuple("ntuple2s","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb2sFerr:pPb2sBerr:RFBErr",numybins);
  TNtuple* ntuple3s = new TNtuple("ntuple3s","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb3sFerr:pPb3sBerr:RFBErr",numybins);

  //BIN LOOP********************************************************
  for (int ihf = -1; ihf<numbins; ihf++) {
    //Choose the hf bin
    if (ihf<0) {
      hfLow = 0;
      hfHigh = 120;
      ntracksLow = 0;
      ntracksHigh = 400;
    }
    else {
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
    }
  for (int iy = 0; iy<numybins/2; iy++) {

    //Choose the rapidity bin
    yLowB = ybinsCM[numybins/2-iy-1];
    yHighB = ybinsCM[numybins/2-iy];
    yLowF = ybinsCM[iy+numybins/2];
    yHighF = ybinsCM[iy+numybins/2+1];

    //print bin label
    if (binmode==0) binLabel = Form("hf%.1f-%.1f_y%.2f-%.2f",hfLow,hfHigh,yLowF,yHighF);
    else if (binmode==1) binLabel = Form("ntracks%i-%i_y%.2f-%.2f",ntracksLow,ntracksHigh,yLowF,yHighF);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelF = kineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString NomFileNameF = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabelF.Data());
    TFile* NomFileF = TFile::Open(NomFileNameF,"READ");
    RooWorkspace *NomwsF = (RooWorkspace*)NomFileF->Get("workspace");
    NomFileF->Close("R");
    float NomYieldF = NomwsF->var(Form("nSig%is",whichUpsilon))->getVal();
    delete NomwsF;
    delete NomFileF;

    kineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelB = kineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString NomFileNameB = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabelB.Data());
    TFile* NomFileB = TFile::Open(NomFileNameB,"READ");
    RooWorkspace *NomwsB = (RooWorkspace*)NomFileB->Get("workspace");
    NomFileB->Close("R");
    float NomYieldB = NomwsB->var(Form("nSig%is",whichUpsilon))->getVal();
    delete NomwsB;
    delete NomFileB;

    TString AltkineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) AltkineLabelF = AltkineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltFileNameF = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabelF.Data());
    TFile* AltFileF = TFile::Open(AltFileNameF,"READ");
    RooWorkspace *AltwsF = (RooWorkspace*)AltFileF->Get("workspace");
    AltFileF->Close("R");
    float AltYieldF = AltwsF->var(Form("nSig%is",whichUpsilon))->getVal();
    delete AltwsF;
    delete AltFileF;

    TString AltkineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) AltkineLabelB = AltkineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltFileNameB = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabelB.Data());
    TFile* AltFileB = TFile::Open(AltFileNameB,"READ");
    RooWorkspace *AltwsB = (RooWorkspace*)AltFileB->Get("workspace");
    AltFileB->Close("R");
    float AltYieldB = AltwsB->var(Form("nSig%is",whichUpsilon))->getVal();
    delete AltwsB;
    delete AltFileB;

    //Calculate errors
    float temperrF = TMath::Abs((AltYieldF-NomYieldF)/NomYieldF*100);
    float temperrB = TMath::Abs((AltYieldB-NomYieldB)/NomYieldB*100);
    float NomRFB = NomYieldF/NomYieldB;
    float AltRFB = AltYieldF/AltYieldB;
    float RFBerr = TMath::Abs((AltRFB-NomRFB)/NomRFB*100);

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
    if (whichUpsilon==1) ntuple1s->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBerr);
    else if (whichUpsilon==2) ntuple2s->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBerr);
    else if (whichUpsilon==3) ntuple3s->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBerr);

  }//end of rapidity bin loop
  }//end of hf bin loop

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1s->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2s->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3s->Write();
  }
  delete ntuple1s;
  delete ntuple2s;
  delete ntuple3s;

  outFile.Close();

}
