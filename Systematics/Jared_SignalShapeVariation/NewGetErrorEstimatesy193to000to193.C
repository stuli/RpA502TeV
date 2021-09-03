#include "NewGetErrorsHeader.h"

using namespace std;
void NewGetErrorEstimatesy193to000to193(int whichUpsilon = 1) {

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
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ptbins2[4] = {0,4,9,30};
  float ptbins3[3] = {0,6,30};

  float* ptbinsptr;
  int numptbinstemp;
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
  }

  const int numptbins = numptbinstemp;

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

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);

  //BIN LOOP********************************************************
  yLowPP = 0.00;
  yHighPP = 1.93;
  for (int iy = 0; iy<2; iy++) {
    if (iy==0) {
      yLowCM = -1.93;
      yHighCM = 0.0;
    }
    else {
      yLowCM = 0.0;
      yHighCM = 1.93;
    }
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = *(ptbinsptr+ipt);
    ptHigh = *(ptbinsptr+ipt+1);
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

    //print bin label
    binLabel = Form("pt%.2f-%.2f_y%.2f-%.2f",ptLow,ptHigh,yLowCM,yHighCM);
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
    if (whichUpsilon==1) ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2) ntuple2spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==3) ntuple3spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);

    theFile->Close();
  }// end of pt loop
  }// end of y loop

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1spt->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2spt->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3spt->Write();
  }
  outFile.Close();


}
