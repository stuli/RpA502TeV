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
void doFitUpsilon_GausCB_MC(
    int collId = kPPMCUps1S,  
    float ptLow=0, float ptHigh=2.5, 
    float yLow=0, float yHigh=2.4,
    int cLow=0, int cHigh=200,
    float muPtCut=4.0
    ) 
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  using namespace RooFit;
  gStyle->SetEndErrorSize(0);

  TString SignalCB = "Double";

  float massLow = 8.5; 
  float massHigh = 10;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*100;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && pt1>%.2f && pt2>%.2f ", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) || (collId==kAAMCUps1S) || (collId==kAAMCUps2S) || (collId==kAAMCUps3S) ) {
    kineCut = kineCut + Form(" && cBin>=%d && cBin<%d ",cLow, cHigh);
    cout << endl;
  }

  TFile* f0;

  if ( collId == kPPMCUps1S) {  
    f0 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121653_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6_bak.root");
  }
  else if ( collId == kAAMCUps1S) {  
    f0 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/skimmedFiles/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016812170_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6_bak.root");
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

  RooDataSet *reducedDS_1 = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y")), *(ws->var("weight"))) , kineCut.Data() );
  reducedDS_1->SetName("reducedDS1");
  RooDataSet *reducedDS   =  new RooDataSet("reducedDS","A sample",*reducedDS_1->get(),Import(*reducedDS_1),WeightVar(*ws->var("weight")));
  reducedDS->SetName("reducedDS");
  cout <<" here " << endl;

  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  //ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"), Layout(0,1,0.95));
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );

  PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  initPset.SetMCSgl();
  initPset.SetMCBkg();

  //integrated
  //RooRealVar alpha1s("alpha1s","tail shift", 1.50 , 1.45, 1.75);
  //RooRealVar n1s("n1s","power order", 3.56 , 3.3, 3.85);
  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.06, 0.001, 0.5);
  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);
  RooRealVar x1s("x1s","sigma fraction 1s 2nd CB", 0.55, 0., 2.5);
  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1", RooArgSet( sigma1s_1, x1s) );
  //RooRealVar alpha1s("alpha1s","tail shift", 1.80, 1.01, 5.65);
  //RooRealVar n1s("n1s","power order", 1.951, 1.001, 5.7128);
  RooRealVar alpha1s("alpha1s","tail shift", 2.121, 1.0012, 7.221);
  RooRealVar n1s("n1s","power order", 2.125, 1.00112, 6.531);

  RooCBShape* cb1s = new RooCBShape("CB", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s, n1s);
  RooGaussian* gauss = new RooGaussian("gauss","gaussian PDF",*(ws->var("mass")),mean1s,sigma1s_2);

  RooAddPdf* CBGaus = new RooAddPdf("CBGaus","Signal 1S",RooArgList(*cb1s,*gauss), RooArgList(*f1s) );

  RooRealVar *sl1 = new RooRealVar("sl1","sl1",-0.01,-4.,3);
  RooRealVar *cnst1 = new RooRealVar("cnst1","cnst1",-0.1,-3.0,3.0);

  RooGenericPdf *bkg;
  //RooGenericPdf *bkg_1order = new RooGenericPdf("bkg","Background","@0*@1+@2",RooArgList( *(ws->var("mass")), sl1, cnst1) );
  RooChebychev *bkg_1order;
  bkg_1order = new RooChebychev("bkg","Background",*(ws->var("mass")),RooArgList(*sl1,*cnst1));

  bkg = (RooGenericPdf*) bkg_1order;

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",20000,0,10000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",100000,0,10000000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",100, 0,1000000);
  RooRealVar *nBkg= new RooRealVar("nBkg","N(Bkg) signals",300, 0,25000);

  RooAddPdf* model = new RooAddPdf();
  if ( (collId==kPPMCUps1S) || (collId==kAAMCUps1S) ) {
    model = new RooAddPdf("model","1S",RooArgList(*CBGaus, *bkg),RooArgList(*nSig1s, *nBkg));
  }

  ws->import(*model);

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2);

  //  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"),MarkerSize(1.2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb1s)),LineColor(kRed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*gauss)),LineColor(kGreen));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed));

  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.5);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  drawText(getCollID(collId),0.12,0.85,1,11);
  drawText(Form("%.2f < p_{T}^{#mu#mu} < %.2f GeV",ptLow,ptHigh ),0.12,0.80,2,11);
  drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), 0.12,0.75,2,11);
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S) 
  {
    drawText(Form("Cent %d-%d%s",cLow/2,cHigh/2,perc.Data()),0.12,0.7,2,12);
    drawText(Form("(p_{T}^{#mu} > %.2f GeV)", muPtCut ), 0.12,0.65,1,12);
  }
  else drawText(Form("(p_{T}^{#mu} > %.2f GeV)", muPtCut ), 0.12,0.65,1,12);
  //  drawText(Form("Signal Function : %s CB", SignalCB.Data() ), 0.55,0.54,1,14);

  // PULL 

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.25);
  pad2->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();  
  pad2->Draw(); 
  pad2->cd();

  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution"));
  pullFrame->addPlotable(hpull,"PX") ;
  pullFrame->SetTitleSize(2.57);
  pullFrame->GetYaxis()->SetTitleOffset(1.8) ;
  pullFrame->GetYaxis()->SetLabelSize(0.16) ;
  pullFrame->GetYaxis()->SetRange(-10,10) ;
  pullFrame->GetXaxis()->SetTitleOffset(0.7) ;
  pullFrame->GetXaxis()->SetLabelSize(0.1) ;
  pullFrame->GetXaxis()->SetTitleSize(0.13) ;
  pullFrame->Draw() ;

  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin; 
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;

  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  drawText(Form("chi^{2}/ndf : %.3f / %d ",chisq,ndf ),0.15,0.95,1,12);

  TPad *pad3 = new TPad("pad3", "pad3", 0.65, 0.55, 0.85, 0.92);
  pad3->SetBottomMargin(0);
  c1->cd();  
  pad3->Draw(); 
  pad3->cd();

  RooPlot* legFrame = ws->var("mass")->frame(Name("Fit Results"), Title("Fit Results"));

  //// Show floating parameters only! (not observables)
  RooArgList paramList = fitRes2->floatParsFinal();
  paramList.Print("v");
  ws->pdf("model")->paramOn(legFrame,Layout(0,.95, .97),Parameters(paramList));
  legFrame->getAttText()->SetTextAlign(11);
  legFrame->getAttText()->SetTextSize(0.09);

  TPaveText* hh = (TPaveText*)legFrame->findObject(Form("%s_paramBox",ws->pdf("model")->GetName()));
  hh->SetY1(0.01); hh->SetY2(0.95);
  hh->Draw();

  c1->SaveAs(Form("fitresults_upsilon_%sCB_%s.png",SignalCB.Data(),kineLabel.Data()));

  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");

  TFile* outf = new TFile(Form("fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()),"recreate");
  outh->Write();
  c1->Write();
  outf->Close();

  cout << "nSig1s = " << ws->var("nSig1s")->getVal() << endl;
  ///  cout parameters :
  //     void setSignalParMC(float MCn_, float MCalpha_, float MCsigma1S_, float MCm0_, float MCf_, float MCx_)
  cout << "N, alpha, sigma1s, M0, f, X double CB in MC : " << endl;
  //   if ( (muPtCut==4) && (ptLow == 0 ) && (ptHigh == 30 ) && (yLow == 0 ) && (yHigh == 2.4 ) )      {
  ofstream of1(Form("params_results_%s_%s.txt",SignalCB.Data(),kineLabel.Data()));
  if ( (collId == kPPMCUps1S) || (collId == kPPMCUps2S) || (collId == kPPMCUps3S) ) 
  {  
    of1 << "if ( binMatched( "<<muPtCut<<", " << ptLow <<", "<< ptHigh << ", "<< yLow<<", "<< yHigh << ") ) " ; }
  else 
  {
    of1 << "if ( binMatched( "<<muPtCut<<", " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ;
  }
//  of1 << "\nn1s, alpha1s, sigma1s_1, m_{#Upsilon(1S)}, f1s, x1s\n";
  of1 << "  { setSignalParMC( " ;
  of1 << ws->var("n1s")->getVal() << ", " <<  ws->var("alpha1s")->getVal() << ", " <<ws->var("sigma1s_1")->getVal() << ", ";
  of1 << ws->var("m_{#Upsilon(1S)}")->getVal()<< ", " <<  ws->var("f1s")->getVal() <<", "<<  ws->var("x1s")->getVal() << " );} "<< endl;
} 

