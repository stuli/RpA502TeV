#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "../../HeaderFiles/PsetCollection.h"
#include "../../HeaderFiles/CMS_lumi.C"
#include "../../HeaderFiles/tdrstyle.C"
#include "../../HeaderFiles/StyleSetting.h"


using namespace std;

void Make_yieldDist_plots( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       const int numtrials = 100
			) 
{

  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  //gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  //import histograms
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString histFileName = Form("PseudoExperimentResults_2018_07_18/PseudoExpResults_%s.root",kineLabel.Data());
  cout << histFileName << endl;
  TFile *histFile = new TFile(histFileName);
  TH1F* myhisto1s = (TH1F*)histFile->Get("1SDiff;1");
  TH1F* myhisto2s = (TH1F*)histFile->Get("2SDiff;1");
  TH1F* myhisto3s = (TH1F*)histFile->Get("3SDiff;1");
  TNtuple* ntupleResults = (TNtuple*)histFile->Get("ntuple;1");

  const int min = 5350;
  const int max = 6850;

  //generate canvas with 2 hisograms
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,800,400);
  cntuple->Divide(2,1);
  TH1F * ntuplehisto1s = new TH1F("NomYield", "Nominal Yield", 100,min,max);
  TH1F * ntuplehisto2s = new TH1F("AltYield", "Alternate Yield", 100,min,max);

  cntuple->cd(1);
  ntuplehisto1s->SetXTitle("Nominal Yield");
  ntuplehisto1s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto1s->GetYaxis()->SetLabelSize(0.04);
  ntuplehisto1s->GetXaxis()->SetLabelSize(0.04);
  ntuplehisto1s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("yield1sNom>>NomYield");

  cntuple->cd(2);
  ntuplehisto2s->SetXTitle("Alternate yield");
  ntuplehisto2s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto2s->GetYaxis()->SetLabelSize(0.04);
  ntuplehisto2s->GetXaxis()->SetLabelSize(0.04);
  ntuplehisto2s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("yield1sAlt>>AltYield");

  cntuple->SaveAs("Yields_Nom_and_Alt.png");

  TCanvas* cNomVsAlt =  new TCanvas("cNomVsAlt","Nom vs Alt Yields",550,45,500,400);
  TH2F * hNomVsAlt = new TH2F("hNomVsAlt", "Nominal vs Alternate Yields", 100,min,max, 100, min, max);
  hNomVsAlt->SetMarkerStyle(4);
  hNomVsAlt->SetMarkerSize(1);
  //Syntax: TNtuple->Draw("yvar:xvar>>TH2F");
  hNomVsAlt->GetXaxis()->SetTitle("Alternate Yield");
  hNomVsAlt->GetYaxis()->SetTitle("Nominal Yield");
  hNomVsAlt->GetYaxis()->SetTitleOffset(1.5);
  ntupleResults->Draw("yield1sNom:yield1sAlt>>hNomVsAlt");
  
  TLine* l1 = new TLine(min,min,max,max);
  l1->Draw("same");

  cNomVsAlt->SaveAs("Yields_Nom_vs_Alt.png");
}
