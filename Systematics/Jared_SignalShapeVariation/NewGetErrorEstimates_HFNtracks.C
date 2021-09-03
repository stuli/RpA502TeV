#include "NewGetErrorsHeader.h"

using namespace std;
void NewGetErrorEstimates_HFNtracks(int whichUpsilon=1,int binmode=0) {
//0=hfmode, 1=ntmode

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
  float hfbins3[3] = {0,12,120};
  int ntbins3[3] = {0,40,400};
  float hfbins1[5] = {0,12,19,27,120};
  int ntbins1[5] = {0,40,62,88,400};

  float* hfbinsptr;
  int* ntbinsptr;
  int numhfbinstemp, numntbinstemp;
  if (whichUpsilon==3) {
    hfbinsptr = &hfbins3[0];
    ntbinsptr = &ntbins3[0];
    numhfbinstemp = sizeof(hfbins3)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins3)/sizeof(float)-1;
  }
  else {
    hfbinsptr = &hfbins1[0];
    ntbinsptr = &ntbins1[0];
    numhfbinstemp = sizeof(hfbins1)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins1)/sizeof(float)-1;
  }

  const int numybins = 2;
  const int numhfbins = numhfbinstemp;
  const int numntbins = numntbinstemp;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -30;
  int max = 30;

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

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is",whichUpsilon);
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
        ntracksLow = *(ntbinsptr+ihf);
        ntracksHigh = *(ntbinsptr+ihf+1);
      }
      else {
        hfLow = *(hfbinsptr+ihf);
        hfHigh = *(hfbinsptr+ihf+1);
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

    float errors[3] = {0};
    float* errorsptr;
    errorsptr = &errors[0];

    //print bin label
    if (binmode==0) binLabel = Form("hf%.1f-%.1f_y%.2f-%.2f",hfLow,hfHigh,yLowF,yHighF);
    else if (binmode==1) binLabel = Form("ntracks%i-%i_y%.2f-%.2f",ntracksLow,ntracksHigh,yLowF,yHighF);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelF = kineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    histFileNameF = Form("%sPseudoExpResults_%s.root",inputDir.Data(),kineLabelF.Data());
    theFileF = TFile::Open(histFileNameF,"READ");
    TNtuple* ntupleResultsF = (TNtuple*)theFileF->Get("ntuple;1");

    kineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelB = kineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    histFileNameB = Form("%sPseudoExpResults_%s.root",inputDir.Data(),kineLabelB.Data());
    theFileB = TFile::Open(histFileNameB,"READ");
    TNtuple* ntupleResultsB = (TNtuple*)theFileB->Get("ntuple;1");

    GetYieldError(ntupleResultsF, ntupleResultsB, whichUpsilon, errorsptr, kineLabelF);

    //Take absolute values
    float temperrB = TMath::Abs(errors[1]);
    float temperrF = TMath::Abs(errors[0]);
    float RFBerr = TMath::Abs(errors[2]);

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
    theFileF->Close();
    theFileB->Close();

  }//end of rapidity bin loop
  }//end of hf bin loop

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) ntuple1s->Write();
  else if (whichUpsilon==2) ntuple2s->Write();
  else if (whichUpsilon==3) ntuple3s->Write();
  outFile.Close();

}
