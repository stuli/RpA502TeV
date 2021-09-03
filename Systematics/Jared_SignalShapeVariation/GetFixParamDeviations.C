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
void GetFixParamDeviations(int whichUpsilon = 1, TString whichParam="alpha") {

  TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  bool PPtoo = kTRUE;

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybinsCM[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybinsCM[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybinsCM[4] = {-2.87,-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;
  const int numtot = numptbins + numybins;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "PT AND RAPIDITY BINS FOR " << whichParam.Data() << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPFileName;
  TFile* theFile;

  //TString whichParam = "alpha";//"alpha","n","x","f";
  TString outFileName = Form("FixParamErrors/Param%sDev%is.root",whichParam.Data(),whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple1sy = new TNtuple("ntuple1sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2sy = new TNtuple("ntuple2sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3sy = new TNtuple("ntuple3sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

  //BIN LOOP********************************************************
  for (int ipt = -1; ipt<numtot; ipt++) {

    if (ipt<numptbins){
      yLowCM = -1.93;
      yHighCM = 1.93;
      yLowPP = 0.00;
      yHighPP = 1.93;
      if (ipt<0) {
        ptLow = 0;
        ptHigh = 30;
      }
      else {
        ptLow = ptbins[ipt];
        ptHigh = ptbins[ipt+1];
      }
      binLow = ptLow;
      binHigh = ptHigh;
      binvar = "pt";
    }
    else {
      ptLow = 0;
      ptHigh = 30;
      yLowCM = ybinsCM[ipt-numptbins];
      yHighCM = ybinsCM[ipt-numptbins+1];
      binLow = yLowCM;
      binHigh = yHighCM;
      binvar = "y";
      if (yLowCM<0) {
        yLowPP = TMath::Abs(yHighCM);
        yHighPP = TMath::Abs(yLowCM);
      }
      else {
        yLowPP = yLowCM;
        yHighPP = yHighCM;
      }
    }

    if (yLowCM<-2.5) {
      PPtoo = kFALSE;
      //cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;
    if (whichParam=="f") PPtoo = kFALSE;
    if (whichParam=="MCf") PPtoo = kFALSE;

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    if (ipt<numptbins) binLabel = strbinLow+"<pt<"+strbinHigh;
    else binLabel = strbinLow+"<y<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import fitted yields
    //cout << "Importing fit results" << endl;
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString NomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabel.Data());
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    NomFile->Close("R");
    float PANomYield = Nomws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Nomws;
    delete NomFile;

    TString AltkineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString AltFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabel.Data());
    TFile* AltFile = TFile::Open(AltFileName,"READ");
    RooWorkspace *Altws = (RooWorkspace*)AltFile->Get("workspace");
    AltFile->Close("R");
    float PAAltYield = Altws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Altws;
    delete AltFile;

    if (PPtoo) {
      TString PPkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString PPNomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),PPkineLabel.Data());
      TFile* PPNomFile = TFile::Open(PPNomFileName,"READ");
      RooWorkspace *PPNomws = (RooWorkspace*)PPNomFile->Get("workspace");
      PPNomFile->Close("R");
      float PPNomYield = PPNomws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPNomws;
      delete PPNomFile;

      TString PPAltkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString PPAltFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),PPAltkineLabel.Data());
      TFile* PPAltFile = TFile::Open(PPAltFileName,"READ");
      RooWorkspace *PPAltws = (RooWorkspace*)PPAltFile->Get("workspace");
      PPAltFile->Close("R");
      float PPAltYield = PPAltws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPAltws;
      delete PPAltFile;
    }

    //Calculate errors
    float temperr = TMath::Abs((PAAltYield-PANomYield)/PANomYield*100);
    float PPerr = 0;
    float RpAerr = 0;
    if (PPtoo) {
      PPerr = TMath::Abs((PPAltYield-PPNomYield)/PPNomYield*100);
      float NomRpA = PANomYield/PPNomYield;
      float AltRpA = PAAltYield/PPAltYield;
      RpAerr = TMath::Abs((AltRpA-NomRpA)/NomRpA*100);
    }

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    cout << endl;

    //put errors in ntuple
    if (whichUpsilon==1)
      if (ipt<numptbins)
        ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2)
      if (ipt<numptbins)
        ntuple2spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple2sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==3)
      if (ipt<numptbins)
        ntuple3spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple3sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
  }

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1spt->Write();
    ntuple1sy->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2spt->Write();
    ntuple2sy->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3spt->Write();
    ntuple3sy->Write();
  }
  delete ntuple1spt;
  delete ntuple1sy;
  delete ntuple2spt;
  delete ntuple2sy;
  delete ntuple3spt;
  delete ntuple3sy;

  outFile.Close();

}
