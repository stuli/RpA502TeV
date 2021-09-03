#include "NewGetErrorsHeader.h"

using namespace std;
void NewGetErrorEstimatespt0to6to30(int whichUpsilon = 1) {

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
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ybins3[3] = {-1.93,0.0,1.93};

  float* ybinsptr;
  int numybinstemp;
  if (whichUpsilon==1) {
    ybinsptr = &ybins1[0];
    numybinstemp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ybinsptr = &ybins2[0];
    numybinstemp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ybinsptr = &ybins3[0];
    numybinstemp = sizeof(ybins3)/sizeof(float)-1;
  }

  const int numybins = numybinstemp;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -10;
  int max = 10;

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPkineLabel, PPFileName;
  TFile* theFile;

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is_pt0to6to30.root",whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1sy = new TNtuple("ntuple1sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple2sy = new TNtuple("ntuple2sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple3sy = new TNtuple("ntuple3sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

  //BIN LOOP********************************************************
  for (int ipt = 0; ipt<2; ipt++) {
    if (ipt==0){
      ptLow = 0;
      ptHigh = 6;
    }
    else {
      ptLow = 6;
      ptHigh = 30;
    }
  for (int iy = 0; iy<numybins; iy++) {

    yLowCM = *(ybinsptr+iy);
    yHighCM = *(ybinsptr+iy+1);
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
  
    if (yLowCM<-2.5) {
      PPtoo = kFALSE;
      //cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;

    float errors[3] = {0};
    float* errorsptr;
    errorsptr = &errors[0];

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    binLabel = strbinLow+"<y<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    histFileName = Form("%sPseudoExpResults_%s.root",inputDir.Data(),kineLabel.Data());
    theFile = TFile::Open(histFileName,"READ");
    TNtuple* ntupleResults = (TNtuple*)theFile->Get("ntuple;1");

    if (PPtoo) {
      PPkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      PPFileName = Form("%sPseudoExpResults_%s.root",inputDir.Data(),PPkineLabel.Data());
      theFile = TFile::Open(PPFileName,"READ");
      TNtuple* PPntupleResults = (TNtuple*)theFile->Get("ntuple;1");
      GetYieldError(ntupleResults, PPntupleResults, whichUpsilon, errorsptr, kineLabel);
    }
    else GetPAOnly(ntupleResults, whichUpsilon, errorsptr, kineLabel);

    //Take absolute values
    float PPerr = TMath::Abs(errors[1]);
    float temperr = TMath::Abs(errors[0]);
    float RpAerr = TMath::Abs(errors[2]);

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    cout << endl;

    //put errors in ntuple
    if (whichUpsilon==1) ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2) ntuple2sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==3) ntuple3sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);

    theFile->Close();
  }// end of y loop
  }// end of pt loop

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1sy->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2sy->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3sy->Write();
  }
  outFile.Close();


}
