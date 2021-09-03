#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void ChangeNtuplestoHistos_y240to193(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%is_y-2.40-1.93.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries()-1;
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);
  const int numybins = (int)ntuple1sy->GetEntries()-1;
  const int numybinsPP = numybins;
  cout << "Number of bins = " << numybins << endl;

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[10] = {-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float ybinsPP[10] = {-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[6] = {-2.4,-1.93,-0.8,0.0,0.8,1.93};
    float ybinsPP[6] = {-2.4,-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[4] = {-2.4,-1.93,0.0,1.93};
    float ybinsPP[4] = {-2.4,-1.93,0.0,1.93};
  }

  //Declare histograms
  cout << "here1" << endl;
  TH1F* hSignalErryPP = new TH1F("hSignalErryPP","hSignalErryPP",numybinsPP,ybinsPP);
  cout << "here2" << endl;
  TH1F* hSignalErryPA = new TH1F("hSignalErryPA","hSignalErryPA",numybins,ybins);
  cout << "here3" << endl;
  TH1F* hSignalErryRpA = new TH1F("hSignalErryRpA","hSignalErryRpA",numybinsPP,ybinsPP);
  cout << "here4" << endl;

//Fill the histograms
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1sy->GetEntry(iy+1);
    TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    hSignalErryPA->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
  }
  for (int iyPP = 0; iyPP<numybinsPP; iyPP++) {
    ntuple1sy->GetEntry(iyPP+1);
    TLeaf *pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    hSignalErryPP->SetBinContent(iyPP+1, TMath::Abs(pp1sErr)/100);
    hSignalErryRpA->SetBinContent(iyPP+1, TMath::Abs(RpPbErr)/100);
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErryPP->Write();
  hSignalErryPA->Write();
  hSignalErryRpA->Write();
  outFile.Close();
}
