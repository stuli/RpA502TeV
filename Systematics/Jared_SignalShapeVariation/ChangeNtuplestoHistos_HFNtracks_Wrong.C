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
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float ybinsPP[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
    float ybinsPP[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[3] = {-1.93,0.0,1.93};
    float ybinsPP[3] = {-1.93,0.0,1.93};
  }
  float hfbins[5] = {0,20,30,40,120};
  int ntracksbins[6] = {0,20,30,40,120,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numybinsPP = numybins;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  cout << "Number of bins = " << numybins << endl;

  //Declare histograms
  cout << "here1" << endl;
  TH2F* hSignalErryhfPA = new TH2F("hSignalErryhfPA","hSignalErryhfPA",numybins,ybins,numhfbins,hfbins);
  cout << "here2" << endl;
  TH2F* hSignalErryhfPP = new TH2F("hSignalErryhfPP","hSignalErryhfPP",numybinsPP,ybinsPP,numhfbins,hfbins);
  cout << "here3" << endl;
  TH2F* hSignalErryhfRpA = new TH2F("hSignalErryhfRpA","hSignalErryhfRpA",numybins,ybins,numhfbins,hfbins);
  cout << "here4" << endl;

//Fill the histograms
  for (int ihf = 0; ihf<numhfbins; ihf++) {
    for (int iy = 0; iy<numybins; iy++) {
      ntuple1sy->GetEntry(iy+ihf*numybins);
      TLeaf *pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
      pp1sErr = (float)pp1sErrLeaf->GetValue();
      TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
      pPb1sErr = (float)pPb1sErrLeaf->GetValue();
      TLeaf *RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
      RpPbErr = (float)RpPbErrLeaf->GetValue();
      hSignalErryhfPP->SetBinContent(iy+1,ihf+1, TMath::Abs(pp1sErr)/100);
      hSignalErryhfPA->SetBinContent(iy+1,ihf+1, TMath::Abs(pPb1sErr)/100);
      hSignalErryhfRpA->SetBinContent(iy+1,ihf+1, TMath::Abs(RpPbErr)/100);
    }
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErryhfPP->Write();
  hSignalErryhfPA->Write();
  hSignalErryhfRpA->Write();
  outFile.Close();
}
