#include "NewGetErrorsHeader.h"

using namespace std;
void NewGetErrorEstimates(int whichUpsilon = 1) {

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
  float ybins1[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[4] = {-2.87,-1.93,0.0,1.93};

  float* ptbinsptr;
  float* ybinsptr;
  int numptbinstemp, numybinstemp;
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    ybinsptr = &ybins1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
    numybinstemp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    ybinsptr = &ybins2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
    numybinstemp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    ybinsptr = &ybins3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
    numybinstemp = sizeof(ybins3)/sizeof(float)-1;
  }

  const int numptbins = numptbinstemp;
  const int numybins = numybinstemp;
  const int numtot = numptbins + numybins;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPkineLabel, PPFileName;
  TFile* theFile;

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is.root",whichUpsilon);
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
        ptLow = *(ptbinsptr+ipt);
        ptHigh = *(ptbinsptr+ipt+1);
      }
      binLow = ptLow;
      binHigh = ptHigh;
      binvar = "pt";
    }
    else {
      ptLow = 0;
      ptHigh = 30;
      yLowCM = *(ybinsptr+ipt-numptbins);
      yHighCM = *(ybinsptr+ipt-numptbins+1);
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
    if (ipt<numptbins) binLabel = strbinLow+"<pt<"+strbinHigh;
    else binLabel = strbinLow+"<y<"+strbinHigh;
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
      GetYieldError(ntupleResults, PPntupleResults, whichUpsilon, errorsptr,kineLabel);
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

    theFile->Close();
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
  outFile.Close();

}
