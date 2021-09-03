#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void ChangeNtuplestoHistosXSBins(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%isXSBins.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%isXSBins.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries()-1;

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
  }

  //Declare histograms
  TH1F* hSignalErrptPA = new TH1F("hSignalErrptPA","hSignalErrptPA",numptbins,ptbins);

//Fill the histograms
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    ntuple1spt->GetEntry(ipt);
    TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf(strpp1sErr);
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf(strpPb1sErr);
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf(strRpPbErr);
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    hSignalErrptPP->SetBinContent(ipt, TMath::Abs(pp1sErr)/100);
    hSignalErrptPA->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
    hSignalErrptRpA->SetBinContent(ipt, TMath::Abs(RpPbErr)/100);
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErrptPA->Write();
  outFile.Close();
}
