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


void GetErrorFromData() {

  int collIdPA = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a set of bins
  int whichUpsilon = 3;

  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybinsCM[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybinsCM[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybinsCM[3] = {-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;
  const int numtot = numptbins + numybins;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  if (whichUpsilon>1) {
    cout << setw(Width1S) << setfill(separator) << Form("PP R%i1 ERR",whichUpsilon);
    cout << setw(Width1S) << setfill(separator) << Form("PA R%i1 ERR",whichUpsilon);
    cout << setw(Width1S) << setfill(separator) << "DR ERR";
  }
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString kineLabel, PANomFileName, PAAltFileName, PPNomFileName, PPAltFileName;
  TFile* theFile;

  //PT LOOP********************************************************
  for (int ipt = -1; ipt<0; ipt++) {

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

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    TString binLabel = strbinLow+"<pt<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import fitted models
    kineLabel = getKineLabel (collIdPA, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    PANomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    theFile = TFile::Open(PANomFileName,"READ");
    float N1SPANom = workspace->var("nSig1s")->getVal();
    float N2SPANom = workspace->var("nSig2s")->getVal();
    float N3SPANom = workspace->var("nSig3s")->getVal();
    //PANomFile->Close("R");

    PAAltFileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    theFile = TFile::Open(PAAltFileName,"READ");
    float N1SPAAlt = workspace->var("nSig1s")->getVal();
    float N2SPAAlt = workspace->var("nSig2s")->getVal();
    float N3SPAAlt = workspace->var("nSig3s")->getVal();
    //PAAltFile->Close("R");

    kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    PPNomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    theFile = TFile::Open(PPNomFileName,"READ");
    float N1SPPNom = workspace->var("nSig1s")->getVal();
    float N2SPPNom = workspace->var("nSig2s")->getVal();
    float N3SPPNom = workspace->var("nSig3s")->getVal();
    //PPNomFile->Close("R");

    PPAltFileName = Form("altfitresults_upsilon_%s.root",kineLabel.Data());
    theFile = TFile::Open(PPAltFileName,"READ");
    float N1SPPAlt = workspace->var("nSig1s")->getVal();
    float N2SPPAlt = workspace->var("nSig2s")->getVal();
    float N3SPPAlt = workspace->var("nSig3s")->getVal();
    //PPAltFile->Close("R");

    //extract parameter values
    if (whichUpsilon==1) {
      float NthisSPANom = N1SPANom;
      float NthisSPAAlt = N1SPAAlt;
      float NthisSPPNom = N1SPPNom;
      float NthisSPPAlt = N1SPPAlt;
    }
    else if (whichUpsilon==2) {
      float NthisSPANom = N2SPANom;
      float NthisSPAAlt = N2SPAAlt;
      float NthisSPPNom = N2SPPNom;
      float NthisSPPAlt = N2SPPAlt;
    }
    else if (whichUpsilon==3) {
      float NthisSPANom = N3SPANom;
      float NthisSPAAlt = N3SPAAlt;
      float NthisSPPNom = N3SPPNom;
      float NthisSPPAlt = N3SPPAlt;
    }

    //calculate error in yield of this upsilon
    float PAthisSerr = 100*(NthisSPAAlt-NthisSPANom)/N1SPANom;
    float PPthisSerr = 100*(NthisSPPAlt-NthisSPPNom)/N1SPPNom;

    //calculate error in RpA
    float pseudoRpANom = NthisSPANom/NthisSPPNom;
    float pseudoRpAAlt = NthisSPAAlt/NthisSPPAlt;
    float RpAerr = 100*(pseudoRpAAlt-pseudoRpANom)/pseudoRpANom;

    if (whichUpsilon>1) {
      //calculate error in R21 or R31
      float R21PANom = NthisSPANom/N1SPANom;
      float R21PAAlt = NthisSPAAlt/N1SPAAlt;
      float PAR21err = 100*(R21PAAlt-R21PANom)/R21PANom;
      float R21PPNom = NthisSPPNom/N1SPPNom;
      float R21PPAlt = NthisSPPAlt/N1SPPAlt;
      float PPR21err = 100*(R21PPAlt-R21PPNom)/R21PPNom;

      //calculate error in DR
      float DRNom = R21PANom/R21PPNom;
      float DRAlt = R21PAAlt/R21PPAlt;
      float DRerr = 100*(DRAlt-DRNom)/DRNom;
    }

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPthisSerr;
    cout << setw(Width1S) << setfill(separator) << PAthisSerr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    if (whichUpsilon>1) {
      cout << setw(Width1S) << setfill(separator) << PPR21err;
      cout << setw(Width1S) << setfill(separator) << PAR21err;
      cout << setw(Width1S) << setfill(separator) << DRerr;
    }
    cout << endl;

  }

}
