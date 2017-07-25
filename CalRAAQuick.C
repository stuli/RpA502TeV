#include <iostream>
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TColor.h"
#include "cutsAndBin.h"
#include "multiTreeUtil.h"

using namespace std;

valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0);

void CalRAAQuick(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0)
{

  TH1::SetDefaultSumw2();

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;

  TFile *wf = new TFile("Ups_RAA_root","recreate");

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nCentBins = nCentBins1s;  centBin = centBin1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
  }
  
  //kinematic range
  ptLow = 15;
  ptHigh = 30;
  yLow = 0;
  yHigh = 2.4;
  cLow = 0;
  cHigh = 200;

  TH1D* hptEffAA[3];
  TH1D* hptEffPP[3];
  TH1D* hptSigAA[3];
  TH1D* hptSigPP[3];

//  hptEffAA = f->Get("hptEffAA");
//  hptEffPP = f->Get("hptEffPP");
  
  int nBin = 1;
  for(int i=0;i<3;i++){
  hptSigPP[i] = new TH1D(Form("hptSigPP%dS",i+1),"",nBin,ptLow,ptHigh);
  hptSigAA[i] = (TH1D*) hptSigPP[i] -> Clone(Form("hptSigAA%dS",i+1));
  }
//  TFile* f = new TFile(Form("efficiency/efficiency_ups%ds.root",state)); 

  for(int i=0;i<3;i++){
  valErr yieldPP = getYield(i+1, kPPDATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,0,100); 
  valErr yieldAA = getYield(i+1, kAADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,0,100); 
  hptSigPP[i] -> SetBinContent(1,yieldPP.val); 
  hptSigPP[i] -> SetBinError(1,yieldPP.err); 
  hptSigAA[i] -> SetBinContent(1,yieldAA.val); 
  hptSigAA[i] -> SetBinError(1,yieldAA.err); 
  }

  TH1D *hRAA[3];
  for(int i=0;i<3;i++){
  hRAA[i] = (TH1D*) hptSigAA[i] -> Clone(Form("hRAA%dS",i+1));
  hRAA[i] -> Divide(hptSigPP[i]);
  hRAA[i] -> Scale(26000000. / 351);
  hRAA[i] -> Scale(2./(208*208));
  hRAA[i] -> SetAxisRange(0,1.6,"Y");
  hRAA[i] -> SetYTitle("R_{AA} (efficiency uncorrected");
  }
  
  TCanvas *cc = new TCanvas("cc","",700,700);
  cc->cd();
    hRAA[0]->Draw();

  TGraphErrors *gr[3];
  for(int i=0;i<3;i++)
  {
    gr[i] = new TGraphErrors();
    gr[i]->SetName("gr");
    gr[i]->SetPoint(0,22+(i+1)*0.5,hRAA[i]->GetBinContent(1));
    gr[i]->SetPointError(0,0,hRAA[i]->GetBinError(1));
    gr[i]->GetHistogram()->GetXaxis()->SetLimits(14,31);
    gr[i]->GetHistogram()->GetXaxis()->SetRangeUser(15,30);
    gr[i]->GetHistogram()->GetYaxis()->SetRangeUser(0,1.6);
    gr[i]->GetHistogram()->GetYaxis()->SetTitle("R_{AA} (efficiency uncorrected)");
    gr[i]->GetHistogram()->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  }
  
  TLegend *l = new TLegend(0.55,0.65,0.8,0.8);
  

  TCanvas* c1 = new TCanvas("c1","",800,400);
  c1->Divide();
  for(int i=0;i<3;i++){
    gr[i]->SetMarkerStyle(20);}
  gr[0]->SetMarkerColor(1);
  gr[1]->SetMarkerColor(2);
  gr[2]->SetMarkerColor(4);
  gr[0]->Draw();
  gr[1]->Draw("Psame");
  gr[2]->Draw("Psame");

  l->AddEntry(gr[0],"#Upsilon (1S)","P");
  l->AddEntry(gr[1],"#Upsilon (2S)","P");
  l->AddEntry(gr[2],"#Upsilon (3S)","P");
  l->SetFillColor(0);
  l->SetTextSize(12);
  l->SetTextFont(43);
  l->Draw("same");

//  hRAA[0]->SetMarkerStyle(24);
//  hRAA[1]->SetMarkerStyle(24);
//  hRAA[2]->SetMarkerStyle(24);
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
		float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString SignalCB = "Double";
  TFile* inf = new TFile(Form("TEST/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

  

  
