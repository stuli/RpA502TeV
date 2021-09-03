#include "MakePseudoExpPlotWithLineHeader.h"

using namespace std;
void MakePseudoExpPlotWithLine() {

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

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << "PP ERR";
  cout << setw(Width1S) << setfill(separator) << "PA ERR";
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPkineLabel, PPFileName;
  TFile* theFile;

  int whichUpsilon = 3;
  ptLow = 0;
  ptHigh = 30;
  yLowCM = -1.93;
  yHighCM = 1.93;
  yLowPP = 0.00;
  yHighPP = 1.93;
  binLow = ptLow;
  binHigh = ptHigh;
  binvar = "pt";

  if (yLowCM<-2.5) {
    PPtoo = kFALSE;
    //cout << "No PP here." << endl;
  }
  else PPtoo = kTRUE;

  float errors[3] = {0};
  float* errorsptr;
  errorsptr = &errors[0];

  float dataerrors[3] = {0.729,0.3408,0.4299};//int bin
  //float dataerrors[3] = {1.77881285548210145e-01,2.54345070570707323e-02,1.52485743165016170e-01};//pt0to2
  //float dataerrors[3] = {1.92089557647705071e-01,1.63206830620765686e-01,2.89303082972764975e-02};//pt2to4
  float* dataerrorsptr;
  dataerrorsptr = &dataerrors[0];

  //print bin label
  stringstream stream1;
  stream1 << fixed << setprecision(2) << binLow;
  string strbinLow = stream1.str();
  stringstream stream2;
  stream2 << fixed << setprecision(2) << binHigh;
  string strbinHigh = stream2.str();
  binLabel = strbinLow+"<pt<"+strbinHigh;
  cout << setw(binColWidth) << setfill(separator) << binLabel;

  //import results of pseudo-experiments
  kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  histFileName = Form("PseudoExperimentResults/PseudoExpResults_%s.root",kineLabel.Data());
  theFile = TFile::Open(histFileName,"READ");
  TNtuple* ntupleResults = (TNtuple*)theFile->Get("ntuple;1");

  if (PPtoo) {
    PPkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    PPFileName = Form("PseudoExperimentResults/PseudoExpResults_%s.root",PPkineLabel.Data());
    theFile = TFile::Open(PPFileName,"READ");
    TNtuple* PPntupleResults = (TNtuple*)theFile->Get("ntuple;1");
    GetYieldError(ntupleResults, PPntupleResults, whichUpsilon, errorsptr, dataerrorsptr, kineLabel);
  }
  else GetPAOnly(ntupleResults, whichUpsilon, errorsptr, dataerrorsptr, kineLabel);

  //Take absolute values
  float PPerr = TMath::Abs(errors[1]);
  float temperr = TMath::Abs(errors[0]);
  float RpAerr = TMath::Abs(errors[2]);

  //print error estimate
  cout << setw(Width1S) << setfill(separator) << PPerr;
  cout << setw(Width1S) << setfill(separator) << temperr;
  cout << setw(Width1S) << setfill(separator) << RpAerr;
  cout << endl;

}
