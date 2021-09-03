#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void ChangeNtuplestoHistosy193to000to193(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/HistoSystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries()/2;
  cout << "Number of bins = " << numptbins << endl;

  int binlowlim, binuplim;
  float binlowlimpt, binuplimpt, pp1sErr, pPb1sErr, RpPbErr;
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
  float ybins[3] = {-1.93,0.0,1.93};
  float numybins = 2;

  //Declare histograms
  cout << "here1" << endl;
  TH2F* hSignalErryptPP = new TH2F("hSignalErryptPP","hSignalErryptPP",numybins,ybins,numptbins,ptbins);
  cout << "here2" << endl;
  TH2F* hSignalErryptPA = new TH2F("hSignalErryptPA","hSignalErryptPA",numybins,ybins,numptbins,ptbins);
  cout << "here3" << endl;
  TH2F* hSignalErryptRpA = new TH2F("hSignalErryptRpA","hSignalErryptRpA",numybins,ybins,numptbins,ptbins);
  cout << "here4" << endl;

//Fill the histograms
  for (int iy = 0; iy<2; iy++) {
    for (int ipt = 0; ipt<numptbins; ipt++) {
      ntuple1spt->GetEntry(ipt+iy*numptbins);
      TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf(strpPb1sErr);
      pPb1sErr = (float)pPb1sErrLeaf->GetValue();
      hSignalErryptPA->SetBinContent(iy+1,ipt+1, TMath::Abs(pPb1sErr)/100);
      TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf(strpp1sErr);
      pp1sErr = (float)pp1sErrLeaf->GetValue();
      hSignalErryptPP->SetBinContent(iy+1,ipt+1, TMath::Abs(pp1sErr)/100);
      TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf(strRpPbErr);
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
