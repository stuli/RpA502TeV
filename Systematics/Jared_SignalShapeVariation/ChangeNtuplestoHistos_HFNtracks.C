#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"

using namespace std;
void ChangeNtuplestoHistos_HFNtracks(int whichUpsilon=1, int binmode=0) {
  //0=hfmode, 1=ntmode

  TString hfntracksbins = "_hfbins.root";
  if (binmode==1) hfntracksbins = "_ntracksbins.root";

  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is",whichUpsilon);
  filename = filename + hfntracksbins;
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%is",whichUpsilon);
  outfilename = outfilename + hfntracksbins;
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntuplename = Form("ntuple%is;1",whichUpsilon);
  TNtuple* ntuple1s = (TNtuple*)inFile->Get(ntuplename);

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RFBErr;
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRFBErr = "RFBErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[5] = {0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[3] = {0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[2] = {0.0,1.93};
  }
  float hfbins[5] = {0,15,22,30,120};
  int ntracksbins[6] = {0,40,65,90,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numybinsPP = numybins;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  cout << "Number of bins = " << numybins << endl;

  //Declare histograms
  cout << "here1" << endl;
  TH2F* hSignalErryhfRFB = new TH2F("hSignalErryhfRFB","hSignalErryhfRFB",numybins,ybins,numhfbins,hfbins);
  cout << "here2" << endl;

//Fill the histograms
  for (int ihf = 0; ihf<numhfbins; ihf++) {
    for (int iy = 0; iy<numybins; iy++) {
      ntuple1s->GetEntry(iy+ihf*numybins);
      TLeaf *RFBErrLeaf = ntuple1s->GetLeaf(strRFBErr);
      RFBErr = (float)RFBErrLeaf->GetValue();
      hSignalErryhfRFB->SetBinContent(iy+1,ihf+1, TMath::Abs(RFBErr)/100);
    }
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErryhfRFB->Write();
  outFile.Close();
}
