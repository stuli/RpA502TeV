#include "rootFitHeaders.h"
#include "commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include "PsetCollection.h"

using namespace std;
using namespace RooFit;
void doFitUpsilon_MC_pT(
       int collId = kPPMCUps1S,  
       float ptLow=0, float ptHigh=6, 
       float yLow=0, float yHigh=1.2,
       int cLow=0, int cHigh=50,
       float muPtCut=4.0
		     ) 
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  

  using namespace RooFit;
  gStyle->SetEndErrorSize(0);
 
  TString SignalCB = "Double";

  float massLow = 0; 
  float massHigh = 100;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && pt1>%.2f && pt2>%.2f ", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) || (collId==kAAMCUps1S) || (collId==kAAMCUps2S) || (collId==kAAMCUps3S) ) {
    kineCut = kineCut + Form(" && cBin>=%d && cBin<%d ",cLow, cHigh);
    cout << endl;
  }
  
  
  TFile* f0;
  
  if ( collId == kPPMCUps1S) {  
    f0 = new TFile("/Users/jaebeom/Desktop/work/CMS/Upsilon_PbPb5TeV_RAA/upsilonRAA5TeV_copy/skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121653_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root");
  }
  else if ( collId == kAAMCUps1S) {  
    f0 = new TFile("/Users/jaebeom/Desktop/work/CMS/Upsilon_PbPb5TeV_RAA/upsilonRAA5TeV_copy/skimmedFiles/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016812170_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root");
    //f0 = new TFile("/Users/jaebeom/Desktop/work/CMS/Upsilon_PbPb5TeV_RAA/upsilonRAA5TeV_copy/skimmedFiles/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168142122_3c54df0419c4813e2d7256dc8952ac699405d027.root");
  }
  /*
    "skimmedfiles/yskimAA_MC_Ups1S00_03_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645164_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 3.10497);
    "skimmedfiles/yskimAA_MC_Ups1S03_06_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645164_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 4.11498);
    "skimmedfiles/yskimAA_MC_Ups1S06_09_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645164_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 2.2579);
    "skimmedfiles/yskimAA_MC_Ups1S09_12_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645165_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 1.2591);
    "skimmedfiles/yskimAA_MC_Ups1S12_15_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645165_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 0.567094);
    "skimmedfiles/yskimAA_MC_Ups1S15_30_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_201645165_96e9f6ebf895348a74ba386ce1eb487d594c0759.root",treeName,"", 0.783399);
  */
  RooDataSet* dataset = (RooDataSet*)f0->Get("dataset");
  
  RooWorkspace *ws = new RooWorkspace(Form("workspace_%s",kineLabel.Data()));
  ws->import(*dataset);
  ws->data("dataset")->Print();
  //  cout << "####################################" << endl;

  RooDataSet *reducedDS_1 = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt1")),*(ws->var("pt2")), *(ws->var("cBin")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight"))) );
  reducedDS_1->SetName("reducedDS1");
  RooDataSet *reducedDS_2  =  new RooDataSet("reducedDS2","A sample",*reducedDS_1->get(),Import(*reducedDS_1),WeightVar(*ws->var("weight")));
  RooDataSet *reducedDS = (RooDataSet*) reducedDS_2->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))));
  reducedDS->SetName("reducedDS");
  cout <<" here " << endl;

 
 ws->import(*reducedDS);
 ws->var("pt")->setRange(massLow, massHigh);
 ws->var("pt")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("pt")->frame(nMassBin); // bins
  //ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"), Layout(0,1,0.95));
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

    
  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
 
  myPlot2->Draw();
  int entries=0;
  for(int i = 1; i<=myPlot2->GetNbinsX();i++)
  {
    entries += myPlot2->GetBinContent(i);
  }
  cout << "Entries : " << entries << endl;
  cout << "Entries : " << myPlot2->GetNbinsX() << endl;
} 

