#include <iostream>
#include <iomanip>
#include <sstream>
#include "../HeaderFiles/rootFitHeaders.h"
#include "../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "../HeaderFiles/PsetCollection.h"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/StyleSetting.h"

using namespace std;
void GetErrorOneBin() {

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a bin
  int whichUpsilon = 2;
  float ptLow = 0;
  float ptHigh = 30;
  float yLowCM = 0.8;
  float yHighCM = 1.93;
  float yLow = yLowCM-0.47;
  float yHigh = yHighCM-0.47;
  bool PPdata = kTRUE;



  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  const int min = -10;
  const int max = 10;



  TCanvas *c1 = new TCanvas("c1","c1",4,45,400,400);

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  if (whichUpsilon>1) {
    cout << setw(Width1S) << setfill(separator) << Form("PA R%i1 ERR",whichUpsilon);
  }
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("RpA ERR",whichUpsilon);
  if (whichUpsilon>1) {
    cout << setw(Width1S) << setfill(separator) << Form("PP R%i1 ERR",whichUpsilon);
    cout << setw(Width1S) << setfill(separator) << Form("DR ERR",whichUpsilon);
  }
  cout << endl;


    //import results of pseudo-experiment
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString histFileName = Form("PseudoExpResults_%s.root",kineLabel.Data());
    TFile *histFile = new TFile(histFileName);
    TNtuple* ntupleResults = (TNtuple*)histFile->Get("ntuple;1");

  if (PPdata) {
    TString kineLabelPP = getKineLabel (collIdPP, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString PPFileName = Form("PseudoExpResults_%s.root",kineLabelPP.Data());
    TFile *PPFile = new TFile(PPFileName);
    TNtuple* PPntupleResults = (TNtuple*)PPFile->Get("ntuple;1");
  }

    //extract error values
    TH1F * ntuplehisto = new TH1F("ntuplehisto", "1S %Diff in Yield", 100,min,max);
    if (whichUpsilon==1) ntupleResults->Draw("diff1s>>ntuplehisto");
    else if (whichUpsilon==2) ntupleResults->Draw("diff2s>>ntuplehisto");
    else if (whichUpsilon==3) ntupleResults->Draw("diff3s>>ntuplehisto");
    float temperr = ntuplehisto->GetMean();  
    float temprms = ntuplehisto->GetRMS();
    if (TMath::Abs(temprms)>TMath::Abs(temperr)) temperr = temprms;

  if (PPdata) {
    if (whichUpsilon==1) PPntupleResults->Draw("diff1s>>ntuplehisto");
    else if (whichUpsilon==2) PPntupleResults->Draw("diff2s>>ntuplehisto");
    else if (whichUpsilon==3) PPntupleResults->Draw("diff3s>>ntuplehisto");
    float PPerr = ntuplehisto->GetMean();
    float PPrms = ntuplehisto->GetRMS();
    if (TMath::Abs(PPrms)>TMath::Abs(PPerr)) PPerr = PPrms;
  }

    //Calculate single and double ratios
  TLeaf *yield1sNomLeaf = ntupleResults->GetLeaf("yield1sNom");
  TLeaf *yield1sAltLeaf = ntupleResults->GetLeaf("yield1sAlt");
  if (whichUpsilon==2) {
    TLeaf *yield2sNomLeaf = ntupleResults->GetLeaf("yield2sNom");
    TLeaf *yield2sAltLeaf = ntupleResults->GetLeaf("yield2sAlt");
  }
  else if (whichUpsilon==3) {
    TLeaf *yield2sNomLeaf = ntupleResults->GetLeaf("yield3sNom");
    TLeaf *yield2sAltLeaf = ntupleResults->GetLeaf("yield3sAlt");
  }
  if (PPdata) {
    TLeaf *PPyield1sNomLeaf = PPntupleResults->GetLeaf("yield1sNom");
    TLeaf *PPyield1sAltLeaf = PPntupleResults->GetLeaf("yield1sAlt");
    if (whichUpsilon==2) {
      TLeaf *PPyield2sNomLeaf = PPntupleResults->GetLeaf("yield2sNom");
      TLeaf *PPyield2sAltLeaf = PPntupleResults->GetLeaf("yield2sAlt");
    }
    else if (whichUpsilon==3) {
      TLeaf *PPyield2sNomLeaf = PPntupleResults->GetLeaf("yield3sNom");
      TLeaf *PPyield2sAltLeaf = PPntupleResults->GetLeaf("yield3sAlt");
    }
    TH1F* hratioPP = new TH1F("hratioPP","hratioPP",100,min,max);
    TH1F* hRpA = new TH1F("hRpA","hRpA",100,min,max);
    TH1F* hDoubleRatio = new TH1F("hDoubleRatio","hDoubleRatio",100,min,max);
  }
  const int NEvents = 100;
  TH1F* hratio = new TH1F("hratio","hratio",100,min,max);
  for (int i = 0; i<NEvents; i++) {
    ntupleResults->GetEntry(i);
    float yield1sNom = (float)yield1sNomLeaf->GetValue();
    float yield1sAlt = (float)yield1sAltLeaf->GetValue();
    if (whichUpsilon>1) {
      float yield2sNom = (float)yield2sNomLeaf->GetValue();
      float yield2sAlt = (float)yield2sAltLeaf->GetValue();
      //PA R21
      float r21Nom = yield2sNom/yield1sNom;
      float r21Alt = yield2sAlt/yield1sAlt;
      float r21Diff = (r21Alt-r21Nom)/r21Nom*100;
      hratio->Fill(r21Diff);
    }
    if (PPdata) {
      PPntupleResults->GetEntry(i);
      float PPyield1sNom = (float)PPyield1sNomLeaf->GetValue();
      float PPyield1sAlt = (float)PPyield1sAltLeaf->GetValue();
      //RpA
      float RpANom = yield1sNom/PPyield1sNom;
      float RpAAlt = yield1sAlt/PPyield1sAlt;
      if (whichUpsilon>1) {
        float PPyield2sNom = (float)PPyield2sNomLeaf->GetValue();
        float PPyield2sAlt = (float)PPyield2sAltLeaf->GetValue();
        RpANom = yield2sNom/PPyield2sNom;
        RpAAlt = yield2sAlt/PPyield2sAlt;
        //PP R21
        float r21NomPP = PPyield2sNom/PPyield1sNom;
        float r21AltPP = PPyield2sAlt/PPyield1sAlt;
        float r21DiffPP = (r21AltPP-r21NomPP)/r21NomPP*100;
        hratioPP->Fill(r21DiffPP);
        //double ratio
        float DR21Nom = r21Nom/r21NomPP;
        float DR21Alt = r21Alt/r21AltPP;
        float DR21Diff = (DR21Alt-DR21Nom)/DR21Nom*100;
        hDoubleRatio->Fill(DR21Diff);
      }
      float RpADiff = (RpAAlt-RpANom)/RpANom*100;
      hRpA->Fill(RpADiff);
    }
  }

  hratio->Draw();
  float ratioerr = hratio->GetMean();//This gives the unbinned value, so accuracy has not bin lost by using a 100-bin histogram.
  float ratiorms = hratio->GetRMS();
  if (TMath::Abs(ratiorms)>TMath::Abs(ratioerr)) ratioerr = ratiorms;
  
  if (PPdata) {
    hratioPP->Draw();
    float ratioerrPP = hratioPP->GetMean();
    float ratiormsPP = hratioPP->GetRMS();
    if (TMath::Abs(ratiormsPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiormsPP;

    hRpA->Draw();
    float RpAerr = hRpA->GetMean();
    float RpArms = hRpA->GetRMS();
    if (TMath::Abs(RpArms)>TMath::Abs(RpAerr)) RpAerr = RpArms;

    hDoubleRatio->Draw();
    float DRerr = hDoubleRatio->GetMean();
    float DRrms = hDoubleRatio->GetRMS();
    if (TMath::Abs(DRrms)>TMath::Abs(DRerr)) DRerr = DRrms;
  }

  //print error estimate
  stringstream stream1;
  stream1 << fixed << setprecision(2) << yLowCM;
  string stryLow = stream1.str();
  stringstream stream2;
  stream2 << fixed << setprecision(2) << yHighCM;
  string stryHigh = stream2.str();
  TString binLabel = stryLow+"<y<"+stryHigh;
  cout << setw(binColWidth) << setfill(separator) << binLabel;
  cout << setw(Width1S) << setfill(separator) << temperr;
  if (whichUpsilon>1) {
    cout << setw(Width1S) << setfill(separator) << ratioerr;
  }
  if (PPdata) {
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    if (whichUpsilon>1) {
      cout << setw(Width1S) << setfill(separator) << ratioerrPP;
      cout << setw(Width1S) << setfill(separator) << DRerr;
    }
  }
    cout << endl;

  c1->Close();
}
