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

void Make_nice_histos( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=2, 
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

  TString chi2cuts = "chisqAlt<2 && chisqNom<2";

  //gStyle->SetStatW(0.4);// Set width of stat-box (fraction of pad size)
  //gStyle->SetStatH(0.2);// Set height of stat-box (fraction of pad size)

  const int min = -10;
  const int max = 10;

  //generate canvas with all 3 hisograms
  /*TCanvas* cyields =  new TCanvas("canvas3","results",4,45,1100,400);
  cyields->Divide(3,1);

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
  myhisto3s->Fit("gaus");*/

  //generate canvas with all 3 hisograms
  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TH1F * ntuplehisto1s = new TH1F("ntuplehisto1s", "1S %Diff in Yield", 100,min,max);
  TH1F * ntuplehisto2s = new TH1F("ntuplehisto2s", "2S %Diff in Yield", 100,min,max);
  TH1F * ntuplehisto3s = new TH1F("ntuplehisto3s", "3S %Diff in Yield", 100,min,max);

  cntuple->cd(1);
  ntuplehisto1s->SetXTitle("%Diff");
  ntuplehisto1s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto1s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto1s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto1s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff1s>>ntuplehisto1s",chi2cuts);
  ntuplehisto1s->Fit("gaus");

  cntuple->cd(2);
  ntuplehisto2s->SetXTitle("%Diff");
  ntuplehisto2s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto2s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto2s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto2s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff2s>>ntuplehisto2s",chi2cuts);
  ntuplehisto2s->Fit("gaus");

  cntuple->cd(3);
  ntuplehisto3s->SetXTitle("%Diff");
  ntuplehisto3s->GetXaxis()->SetTitleSize(0.05);
  ntuplehisto3s->GetYaxis()->SetLabelSize(0.05);
  ntuplehisto3s->GetXaxis()->SetLabelSize(0.05);
  ntuplehisto3s->GetXaxis()->SetRangeUser(min,max);
  ntupleResults->Draw("diff3s>>ntuplehisto3s",chi2cuts);
  ntuplehisto3s->Fit("gaus");

  TCanvas* cchisqNom =  new TCanvas("cchisqNom","chisqNomResults",550,45,500,400);
  TH2F * hchisqNom = new TH2F("hchisqNom", "chisqNom vs diff3s", 100,-70,70, 100, 0, 10);
  hchisqNom->SetMarkerStyle(4);
  hchisqNom->SetMarkerSize(1);
  //Syntax: TNtuple->Draw("yvar:xvar>>TH2F");
  ntupleResults->Draw("chisqNom:diff2s>>hchisqNom",chi2cuts);

  TCanvas* cchisqAlt =  new TCanvas("cchisqAlt","chisqAltResults",550,545,500,400);
  TH2F * hchisqAlt = new TH2F("hchisqAlt", "chisqAlt vs diff3s", 100,-70,70, 100, 0, 10);
  hchisqAlt->SetMarkerStyle(4);
  hchisqAlt->SetMarkerSize(1);
  //Syntax: TNtuple->Draw("yvar:xvar>>TH2F");
  ntupleResults->Draw("chisqAlt:diff2s>>hchisqAlt",chi2cuts);

  //TCanvas* diff3s =  new TCanvas("diff3s","diff3s",550,45,500,400);
  //TH1F * diff3sntuple = new TH1F("diff3sntuple", "diff3sntuple", 100,min,max);
  //ntupleResults->Draw("diff3s>>diff3sntuple");

  TCanvas* ctest =  new TCanvas("ctest","ctest",4,45,500,400);
  TH1F* chisqAltHist = new TH1F("chisqAltHist","chisqAltHist",100,0.5,3.5);
  ntupleResults->Draw("chisqAlt>>chisqAltHist");

  //Calculate single and double ratios
  TLeaf *yield1sNomLeaf = ntupleResults->GetLeaf("yield1sNom");
  TLeaf *yield2sNomLeaf = ntupleResults->GetLeaf("yield2sNom");
  TLeaf *yield3sNomLeaf = ntupleResults->GetLeaf("yield3sNom");
  TLeaf *yield1sAltLeaf = ntupleResults->GetLeaf("yield1sAlt");
  TLeaf *yield2sAltLeaf = ntupleResults->GetLeaf("yield2sAlt");
  TLeaf *yield3sAltLeaf = ntupleResults->GetLeaf("yield3sAlt");
  const int NEvents = (int)ntupleResults->GetEntries();
  TCanvas* cRatios =  new TCanvas("cRatios","cRatios",4,45,800,400);
  cRatios->Divide(2,1);
  TH1F* hR21 = new TH1F("hR21","hR21",100,min,max);
  TH1F* hR31 = new TH1F("hR31","hR31",100,min,max);
  for (int i = 0; i<NEvents; i++) {
    //single ratios
    ntupleResults->GetEntry(i);
    float yield1sNom = (float)yield1sNomLeaf->GetValue();
    float yield2sNom = (float)yield2sNomLeaf->GetValue();
    float yield3sNom = (float)yield3sNomLeaf->GetValue();
    float yield1sAlt = (float)yield1sAltLeaf->GetValue();
    float yield2sAlt = (float)yield2sAltLeaf->GetValue();
    float yield3sAlt = (float)yield3sAltLeaf->GetValue();
    float r21Nom = yield2sNom/yield1sNom;
    float r21Alt = yield2sAlt/yield1sAlt;
    float r31Nom = yield3sNom/yield1sNom;
    float r31Alt = yield3sAlt/yield1sAlt;
    float r21Diff = (r21Alt-r21Nom)/r21Nom*100;
    float r31Diff = (r31Alt-r31Nom)/r31Nom*100;
    cout << r21Diff << "; " << r31Diff << endl;
    hR21->Fill(r21Diff);
    hR31->Fill(r31Diff);
    //double ratios: Need pp data.

  }
  cRatios->cd(1);
  hR21->SetXTitle("%Diff");
  hR21->GetXaxis()->SetTitleSize(0.05);
  hR21->GetYaxis()->SetLabelSize(0.05);
  hR21->GetXaxis()->SetLabelSize(0.05);
  hR21->Draw();
  hR21->Fit("gaus");

  cRatios->cd(2);
  hR31->SetXTitle("%Diff");
  hR31->GetXaxis()->SetTitleSize(0.05);
  hR31->GetYaxis()->SetLabelSize(0.05);
  hR31->GetXaxis()->SetLabelSize(0.05);
  hR31->Draw();
  hR31->Fit("gaus");

}
