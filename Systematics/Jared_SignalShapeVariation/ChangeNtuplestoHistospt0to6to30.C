#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void ChangeNtuplestoHistospt0to6to30(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is_pt0to6to30.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%is_pt0to6to30.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);
  const int numybins = (int)ntuple1sy->GetEntries()/2;
  const int numybinsPP = numybins;
  cout << "Number of bins = " << numybins << endl;

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[10] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float ybinsPP[10] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[6] = {-1.93,-0.8,0.0,0.8,1.93};
    float ybinsPP[6] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[4] = {-1.93,0.0,1.93};
    float ybinsPP[4] = {-1.93,0.0,1.93};
  }
  float ptbins[3] = {0,6,30};
  float numptbins = 2;

  //Declare histograms
  cout << "here1" << endl;
  TH2F* hSignalErryptPP = new TH2F("hSignalErryptPP","hSignalErryptPP",numybins,ybins,numptbins,ptbins);
  cout << "here2" << endl;
  TH2F* hSignalErryptPA = new TH2F("hSignalErryptPA","hSignalErryptPA",numybins,ybins,numptbins,ptbins);
  cout << "here3" << endl;
  TH2F* hSignalErryptRpA = new TH2F("hSignalErryptRpA","hSignalErryptRpA",numybins,ybins,numptbins,ptbins);
  cout << "here4" << endl;

//Fill the histograms
  for (int ipt = 0; ipt<2; ipt++) {
    for (int iy = 0; iy<numybins; iy++) {
      ntuple1sy->GetEntry(iy+ipt*numybins);
      TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
      pPb1sErr = (float)pPb1sErrLeaf->GetValue();
      hSignalErryptPA->SetBinContent(iy+1,ipt+1, TMath::Abs(pPb1sErr)/100);
      TLeaf *pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
      pp1sErr = (float)pp1sErrLeaf->GetValue();
      hSignalErryptPP->SetBinContent(iy+1,ipt+1, TMath::Abs(pp1sErr)/100);
      TLeaf *RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
      RpPbErr = (float)RpPbErrLeaf->GetValue();
      hSignalErryptRpA->SetBinContent(iy+1,ipt+1, TMath::Abs(RpPbErr)/100);
    }
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErryptPP->Write();
  hSignalErryptPA->Write();
  hSignalErryptRpA->Write();
  outFile.Close();
}
