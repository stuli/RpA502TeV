#include "../../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../../cutsAndBin.h"
#include <TGraphErrors.h>
using namespace std;

void makeTAATrackingUnc(int state=1)
{
  double pp ;  // placeholder for pp
  double aa ; // placeholder for PbPb
  double finalUnc; 
  
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double *yBin;
  double *TAA_unc;

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
    nCentBins = nCentBins1s;  centBin = centBin1s; TAA_unc = TAA_unc1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
    nYBins = nYBins2S;  yBin = yBin2S; TAA_unc = TAA_unc2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
    nYBins = nYBins3S;  yBin = yBin3S; TAA_unc = TAA_unc3s;
  }
  
  TH1D *hUncPtPP = new TH1D("hptPP",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapPP = new TH1D("hrapPP",";Uncertainty;|y|",nYBins,yBin);
  TH1D *hUncIntPP = new TH1D("hIntPP",";Uncertainty;Centrality",1,0,200);
  
  TH1D *hUncCentAA = new TH1D("hcentAA","; Uncertainty ; Centrality Bins", nCentBins, centBin);
  TH1D *hUncPtAA = new TH1D("hptAA",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapAA = new TH1D("hrapAA",";Uncertainty;|y|",nYBins,yBin);
  TH1D *hUncIntAA = new TH1D("hIntAA",";Uncertainty;Centrality",1,0,200);
  
  TH1D *hUncIntAAoPP = new TH1D("hIntAAoPP",";Uncertainty;Centrality",1,0,200);
  TH1D *hUncCentAAoPP = new TH1D("hcentAAoPP","; Uncertainty ; Centrality Bins", nCentBins, centBin);
  TH1D *hUncPtAAoPP = new TH1D("hptAAoPP",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapAAoPP = new TH1D("hrapAAoPP",";Uncertainty;|y|",nYBins,yBin);


  cout << "state : " << state << endl;

  double tracking_pp = 0.01;
  double tracking_aa = 0.012;

  double trackingUnc = TMath::Sqrt(tracking_pp*tracking_pp+tracking_aa*tracking_aa);

  //Rapidity
  for(int irap=0;irap<nYBins;irap++)
  {
    hUncRapPP->SetBinContent(irap+1, tracking_pp);
    hUncRapAA->SetBinContent(irap+1, tracking_aa);
    pp = hUncRapPP-> GetBinContent(irap+1);
    aa = hUncRapAA-> GetBinContent(irap+1);
    finalUnc = TMath::Sqrt(pp*pp+aa*aa);
    hUncRapAAoPP->SetBinContent(irap+1,finalUnc);
  }

  //Pt
  for(int ipt=0;ipt<nPtBins;ipt++)
  {
    hUncPtPP->SetBinContent(ipt+1, tracking_pp);
    hUncPtAA->SetBinContent(ipt+1, tracking_aa);
    pp = hUncPtPP-> GetBinContent(ipt+1);
    aa = hUncPtAA-> GetBinContent(ipt+1);
    finalUnc = TMath::Sqrt(pp*pp+aa*aa);
    hUncPtAAoPP->SetBinContent(ipt+1,finalUnc);
  }
  
  //Int
  hUncIntPP->SetBinContent(1,tracking_pp);
  hUncIntAA->SetBinContent(1,TMath::Sqrt(TAA_unc[nCentBins]*TAA_unc[nCentBins]+tracking_aa*tracking_aa));
  pp = hUncIntPP->GetBinContent(1);
  aa = hUncIntAA->GetBinContent(1);
  finalUnc = TMath::Sqrt(pp*pp+aa*aa);
  hUncIntAAoPP->SetBinContent(1,finalUnc);

  //Cent
  for(int icent=0;icent<nCentBins;icent++)
  {
    hUncCentAA->SetBinContent(icent+1, TMath::Sqrt(TAA_unc[icent]*TAA_unc[icent]+tracking_aa*tracking_aa));
    pp = hUncIntPP->GetBinContent(1);
    aa = hUncCentAA-> GetBinContent(icent+1);
    finalUnc = TMath::Sqrt(pp*pp+aa*aa);
    hUncCentAAoPP->SetBinContent(icent+1,finalUnc);
  }


  TFile *wf = new TFile(Form("sys_TAA_Tracking_%ds.root",state),"recreate");

  wf->cd();
  hUncRapPP->Write();
  hUncRapAA->Write();
  hUncRapAAoPP->Write();
  hUncPtPP->Write();
  hUncPtAA->Write();
  hUncPtAAoPP->Write();
  hUncIntAA->Write();
  hUncIntPP->Write();
  hUncIntAAoPP->Write();
  hUncCentAA->Write();
  hUncCentAAoPP->Write();
}








