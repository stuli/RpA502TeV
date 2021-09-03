#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"


using namespace std;

void Make_nice_histos() {

  //gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  //import histograms
  const char *inputFile="systematic_histos_Final.root";
  TFile *thefile = new TFile(inputFile);
  TH1F* myhisto1s = (TH1F*)thefile->Get("1SDiff;1");
  TH1F* myhisto2s = (TH1F*)thefile->Get("2SDiff;1");
  TH1F* myhisto3s = (TH1F*)thefile->Get("3SDiff;1");
  TNtuple* ntupleResults = (TNtuple*)thefile->Get("ntuple;1");

  //generate canvas with all 3 hisograms
  TCanvas* cyields =  new TCanvas("canvas3","results",4,45,1100,400);
  cyields->Divide(3,1);

  gStyle->SetStatW(0.4);// Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.2);// Set height of stat-box (fraction of pad size)

  const int min = -70;
  const int max = 70;

  cyields->cd(1);
  myhisto1s->SetXTitle("%Diff");
  myhisto1s->GetXaxis()->SetTitleSize(0.05);
  myhisto1s->GetYaxis()->SetLabelSize(0.05);
  myhisto1s->GetXaxis()->SetLabelSize(0.05);
  myhisto1s->GetXaxis()->SetRangeUser(min,max);
  myhisto1s->Draw();
  myhisto1s->Fit("gaus");

  cyields->cd(2);
  myhisto2s->SetXTitle("%Diff");
  myhisto2s->GetXaxis()->SetTitleSize(0.05);
  myhisto2s->GetYaxis()->SetLabelSize(0.05);
  myhisto2s->GetXaxis()->SetLabelSize(0.05);
  myhisto2s->GetXaxis()->SetRangeUser(min,max);
  myhisto2s->Draw();
  myhisto2s->Fit("gaus");

  cyields->cd(3);
  myhisto3s->SetXTitle("%Diff");
  myhisto3s->GetXaxis()->SetTitleSize(0.05);
  myhisto3s->GetYaxis()->SetLabelSize(0.05);
  myhisto3s->GetXaxis()->SetLabelSize(0.05);
  myhisto3s->GetXaxis()->SetRangeUser(min,max);
  myhisto3s->Draw();
  myhisto3s->Fit("gaus");

  //generate canvas with all 3 hisograms
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TH1F * ntuplehisto1s = new TH1F("ntuplehisto1s", "ntuplehisto1s", 100,-100,100);
  TH1F * ntuplehisto2s = new TH1F("ntuplehisto2s", "ntuplehisto2s", 100,-100,100);
  TH1F * ntuplehisto3s = new TH1F("ntuplehisto3s", "ntuplehisto3s", 100,-100,100);

  cntuple->cd(1);
  ntuplehisto1s->SetXTitle("%Diff");
  ntuplehisto1s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto1s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto1s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto1s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff1s>>ntuplehisto1s","chisqAlt<3 && chisqNom<2");
  ntuplehisto1s->Fit("gaus");

  cntuple->cd(2);
  ntuplehisto2s->SetXTitle("%Diff");
  ntuplehisto2s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto2s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto2s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto2s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff2s>>ntuplehisto2s","chisqAlt<3 && chisqNom<2");
  ntuplehisto2s->Fit("gaus");

  cntuple->cd(3);
  ntuplehisto3s->SetXTitle("%Diff");
  ntuplehisto3s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto3s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto3s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto3s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff3s>>ntuplehisto3s","chisqAlt<3 && chisqNom<2");
  ntuplehisto3s->Fit("gaus");

  TCanvas* cchisq =  new TCanvas("cchisq","chisqResults",550,45,500,400);
  TH2F * hchisq = new TH2F("hchisq", "chisqNom vs diff3s", 100,-70,70, 100, 0, 10);
  hchisq->SetMarkerStyle(4);
  hchisq->SetMarkerSize(1);
  //Syntax: TNtuple->Draw("yvar:xvar>>TH2F");
  ntupleResults->Draw("chisqAlt:diff2s>>hchisq","chisqAlt<3 && chisqNom<2");

  //TCanvas* diff3s =  new TCanvas("diff3s","diff3s",550,45,500,400);
  //TH1F * diff3sntuple = new TH1F("diff3sntuple", "diff3sntuple", 100,min,max);
  //ntupleResults->Draw("diff3s>>diff3sntuple");

  TCanvas* ctest =  new TCanvas("ctest","ctest",4,45,500,400);
  TH1F* chisqAltHist = new TH1F("chisqAltHist","chisqAltHist",100,0.5,3.5);
  ntupleResults->Draw("chisqAlt>>chisqAltHist");

}
