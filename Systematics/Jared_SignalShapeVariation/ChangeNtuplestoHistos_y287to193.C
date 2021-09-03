#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void ChangeNtuplestoHistos_y287to193(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%is_y-2.87-1.93.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries()-1;
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);
  const int numybins = (int)ntuple1sy->GetEntries();
  const int numybinsPP = numybins;
  cout << "Number of bins = " << numybins << endl;

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[11] = {-2.87,-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[7] = {-2.87,-2.4,-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[5] = {-2.87,-2.4,-1.93,0.0,1.93};
  }

  //Declare histograms
  cout << "here1" << endl;
  TH1F* hSignalErryPA = new TH1F("hSignalErryPA","hSignalErryPA",numybins,ybins);
  cout << "here2" << endl;

//Fill the histograms
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1sy->GetEntry(iy);
    TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    hSignalErryPA->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErryPA->Write();
  outFile.Close();
}
