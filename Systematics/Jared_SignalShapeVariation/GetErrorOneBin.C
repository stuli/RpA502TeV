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
void GetErrorOneBin() {

  //gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  //choose a bin
  int whichUpsilon = 3;
  float ptLow = 0;
  float ptHigh = 6;
  float yLowCM = -1.93;
  float yHighCM = 1.93;
  float yLowPP = 0.0;
  float yHighPP = yHighCM;
  bool PPdata = kTRUE;


  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  const int min = -30;
  const int max = 30;

  TCanvas *c1 = new TCanvas("c1","c1",4,45,400,400);

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
  TString kineLabel, histFileName, PPFileName;
  TFile* theFile;

  //PT LOOP********************************************************


    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    TString binLabel = strbinLow+"<pt<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    histFileName = Form("PseudoExpResults_%s.root",kineLabel.Data());
    theFile = new TFile(histFileName);
    TNtuple* ntupleResults = (TNtuple*)theFile->Get("ntuple;1");

    kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    PPFileName = Form("PseudoExpResults_%s.root",kineLabel.Data());
    theFile = new TFile(PPFileName);
    TNtuple* PPntupleResults = (TNtuple*)theFile->Get("ntuple;1");

    //extract error values
    TH1F * ntuplehisto = new TH1F("ntuplehisto", "1S %Diff in Yield", 100,min,max);
    if (whichUpsilon==1) ntupleResults->Draw("diff1s>>ntuplehisto");
    else if (whichUpsilon==2) ntupleResults->Draw("diff2s>>ntuplehisto");
    else if (whichUpsilon==3) ntupleResults->Draw("diff3s>>ntuplehisto");
    float temperr = ntuplehisto->GetMean();  
    float temprms = ntuplehisto->GetRMS();
    if (TMath::Abs(temprms)>TMath::Abs(temperr)) temperr = temprms;
    ntuplehisto->Fit("gaus","Q");
    float tempmu = ntuplehisto->GetFunction("gaus")->GetParameter(1);//fitted mean
    if (TMath::Abs(tempmu)>TMath::Abs(temperr)) temperr = tempmu;
    float tempsigma = ntuplehisto->GetFunction("gaus")->GetParameter(2);//fitted sigma
    if (TMath::Abs(tempsigma)>TMath::Abs(temperr)) temperr = tempsigma;

    if (whichUpsilon==1) PPntupleResults->Draw("diff1s>>ntuplehisto");
    else if (whichUpsilon==2) PPntupleResults->Draw("diff2s>>ntuplehisto");
    else if (whichUpsilon==3) PPntupleResults->Draw("diff3s>>ntuplehisto");
    float PPerr = ntuplehisto->GetMean();
    float PPrms = ntuplehisto->GetRMS();
    if (TMath::Abs(PPrms)>TMath::Abs(PPerr)) PPerr = PPrms;
    ntuplehisto->Fit("gaus","Q");
    float PPmu = ntuplehisto->GetFunction("gaus")->GetParameter(1);//fitted mean
    if (TMath::Abs(PPmu)>TMath::Abs(PPerr)) PPerr = PPmu;
    float PPsigma = ntuplehisto->GetFunction("gaus")->GetParameter(2);//fitted sigma
    if (TMath::Abs(PPsigma)>TMath::Abs(PPerr)) PPerr = PPsigma;

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

    if (whichUpsilon>1) {
      hratio->Draw();
      float ratioerr = hratio->GetMean();
      float ratiorms = hratio->GetRMS();
      if (TMath::Abs(ratiorms)>TMath::Abs(ratioerr)) ratioerr = ratiorms;
      hratio->Fit("gaus","Q");
      float ratiomu = hratio->GetFunction("gaus")->GetParameter(1);//fitted mean
      if (TMath::Abs(ratiomu)>TMath::Abs(ratioerr)) ratioerr = ratiomu;
      float ratiosigma = hratio->GetFunction("gaus")->GetParameter(2);//fitted sigma
      if (TMath::Abs(ratiosigma)>TMath::Abs(ratioerr)) ratioerr = ratiosigma;

      hratioPP->Draw();
      float ratioerrPP = hratioPP->GetMean();
      float ratiormsPP = hratioPP->GetRMS();
      if (TMath::Abs(ratiormsPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiormsPP;
      hratioPP->Fit("gaus","Q");
      float ratiomuPP = hratioPP->GetFunction("gaus")->GetParameter(1);//fitted mean
      if (TMath::Abs(ratiomuPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiomuPP;
      float ratiosigmaPP = hratioPP->GetFunction("gaus")->GetParameter(2);//fitted sigma
      if (TMath::Abs(ratiosigmaPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiosigmaPP;
      hDoubleRatio->Draw();
      float DRerr = hDoubleRatio->GetMean();
      float DRrms = hDoubleRatio->GetRMS();
      if (TMath::Abs(DRrms)>TMath::Abs(DRerr)) DRerr = DRrms;
      hDoubleRatio->Fit("gaus","Q");
      float DRmu = hDoubleRatio->GetFunction("gaus")->GetParameter(1);//fitted mean
      if (TMath::Abs(DRmu)>TMath::Abs(DRerr)) DRerr = DRmu;
      float DRsigma = hDoubleRatio->GetFunction("gaus")->GetParameter(2);//fitted sigma
      if (TMath::Abs(DRsigma)>TMath::Abs(DRerr)) DRerr = DRsigma;
    }

    hRpA->Draw();
    float RpAerr = hRpA->GetMean();
    float RpArms = hRpA->GetRMS();
    if (TMath::Abs(RpArms)>TMath::Abs(RpAerr)) RpAerr = RpArms;
    hRpA->Fit("gaus","Q");
    float RpAmu = hRpA->GetFunction("gaus")->GetParameter(1);//fitted mean
    if (TMath::Abs(RpAmu)>TMath::Abs(RpAerr)) RpAerr = RpAmu;
    float RpAsigma = hRpA->GetFunction("gaus")->GetParameter(2);//fitted sigma
    if (TMath::Abs(RpAsigma)>TMath::Abs(RpAerr)) RpAerr = RpAsigma;

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    if (whichUpsilon>1) {
      cout << setw(Width1S) << setfill(separator) << ratioerrPP;
      cout << setw(Width1S) << setfill(separator) << ratioerr;
      cout << setw(Width1S) << setfill(separator) << DRerr;
    }
    cout << endl;


  c1->Close();

  //generate a pretty canvas to present that shows pp yield, pPb yield, and RpA
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
    TH1F * ntuplehistoPP = new TH1F("ntuplehistoPP", "PP 1S %Diff in Yield", 100,min,max);
    TH1F * ntuplehistoPA = new TH1F("ntuplehistoPA", "PA 1S %Diff in Yield", 100,min,max);

  cntuple->cd(1);
  ntuplehistoPP->SetXTitle("%Diff");
  ntuplehistoPP->GetXaxis()->SetTitleSize(0.05);
  ntuplehistoPP->GetYaxis()->SetLabelSize(0.05);
  ntuplehistoPP->GetXaxis()->SetLabelSize(0.05);
  ntuplehistoPP->GetXaxis()->SetRangeUser(min,max);
  PPntupleResults->Draw("diff3s>>ntuplehistoPP");
  ntuplehistoPP->Fit("gaus");

  cntuple->cd(2);
  ntuplehistoPA->SetXTitle("%Diff");
  ntuplehistoPA->GetXaxis()->SetTitleSize(0.05);
  ntuplehistoPA->GetYaxis()->SetLabelSize(0.05);
  ntuplehistoPA->GetXaxis()->SetLabelSize(0.05);
  ntuplehistoPA->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff3s>>ntuplehistoPA");
  ntuplehistoPA->Fit("gaus");

  cntuple->cd(3);
  hRpA->SetXTitle("%Diff");
  hRpA->GetXaxis()->SetTitleSize(0.05);
  hRpA->GetYaxis()->SetLabelSize(0.05);
  hRpA->GetXaxis()->SetLabelSize(0.05);
  hRpA->GetXaxis()->SetRangeUser(min,max);
  hRpA->Draw();
  //hRpA->Fit("gaus");

  cntuple->SaveAs("PseudoExpPlots.png");
}
