
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
void doFitUpsilon_Data( 
       int collId = kAADATA,  
       float ptLow=4, float ptHigh=6, 
       float yLow=0, float yHigh=1.2,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool fixParameters=1
			) 
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  

  using namespace RooFit;
  gStyle->SetEndErrorSize(0);
 
  TString SignalCB = "Double";

  float massLow = 8; 
  float massHigh = 13;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;
  TFile* f1;
  if      ( collId == kPPDATA) f1 = new TFile("/home/deathold/work/CMS/analysis/quarkonia/skimmedFiles/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20164251755_3964bbec2f15f2cf9baa0676644690f40cee27c4.root");
  else if ( collId == kAADATA) f1 = new TFile("/home/deathold/work/CMS/analysis/quarkonia/skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20164272229_95c28a5bdf107c32b9e54843b8c85939ffe1aa23.root");
  else if ( collId == kAADATAPeri) f1 = new TFile("/home/deathold/work/CMS/analysis/quarkonia/skimmedFiles/yskimPbPb_PeripheralPD_Trig-L1DoubleMu0Peripheral_OpSign_EP-OppositeHF_20164272252_95c28a5bdf107c32b9e54843b8c85939ffe1aa23.root");
  else if ( collId == kPPMCUps1S) f1 = new TFile("/home/deathold/work/CMS/analysis/quarkonia/skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20163251233_2b58ba03c4751c9d10cb9d60303271ddd6e1ba3a.root");
  else if ( collId == kAAMCUps1S) f1 = new TFile("/home/deathold/work/CMS/analysis/quarkonia/skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20163251233_2b58ba03c4751c9d10cb9d60303271ddd6e1ba3a.root");
 
  if(collId == kAADATAPeri) collId =2; 
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f)", (float)muPtCut, (float)muPtCut );
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATACentL3) || (collId==kAADATAPeri) )
    kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d) && ( abs(abs(dphiEp2/3.141592)-0.5)>%.3f && abs(abs(dphiEp2/3.141592)-0.5)<%.3f )",cLow, cHigh, dphiEp2Low, dphiEp2High);
  
  
  TTree* tree = (TTree*) f1->Get("mm");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace(Form("workspace_%s",kineLabel.Data()));
  ws->import(*dataset);
  ws->data("dataset")->Print();
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
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
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
          
  PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  initPset.SetMCSgl();

  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma2s_1("sigma2s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma3s_1("sigma3s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma1s_2("sigma1s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma2s_2("sigma2s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma3s_2("sigma3s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  
  RooRealVar alpha1s_1("alpha1s_1","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha2s_1("alpha2s_1","tail shift", 5. , 1.5, 9.2);  
  RooRealVar alpha3s_1("alpha3s_1","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha1s_2("alpha1s_2","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha2s_2("alpha2s_2","tail shift", 5. , 1.5, 9.2);  
  RooRealVar alpha3s_2("alpha3s_2","tail shift", 5. , 1.5, 9.8);
  
  RooRealVar n1s_1("n1s_1","power order", 5. , 1.5, 9.8);
  RooRealVar n2s_1("n2s_1","power order", 5. , 1.5, 10.);
  RooRealVar n3s_1("n3s_1","power order", 4. , 1.5, 9.8);
  RooRealVar n1s_2("n1s_2","power order", 5. , 1.5, 9.8);
  RooRealVar n2s_2("n2s_2","power order", 5. , 1.5, 10.);
  RooRealVar n3s_2("n3s_2","power order", 4. , 1.5, 9.8);
  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);
  RooRealVar *f2s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);
  RooRealVar *f3s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);

  // Fix the parameters 
  if (fixParameters) 
  { 
    if ( initPset.n1s_1 == -1 )
      {
	cout << endl << endl << endl << "#########################  ERROR!!!! ##################" << endl;
	cout << "No Param. set for " << kineLabel << ","<<endl;
	cout << "Fitting macro is stopped!" << endl << endl << endl;
	return;
      }
    else { 
      cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
      cout << endl << "Fixing the parameters..." << endl << endl;
      cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
      cout << "initPset.n1s_1 = " << initPset.n1s_1 << endl;
      n1s_1.setVal(initPset.n1s_1);  n1s_1.setConstant(); 
      n2s_1.setVal(initPset.n2s_1);  n2s_1.setConstant();  
      n3s_1.setVal(initPset.n3s_1);  n3s_1.setConstant();  
      n1s_2.setVal(initPset.n1s_2);  n1s_2.setConstant(); 
      n2s_2.setVal(initPset.n2s_2);  n2s_2.setConstant();  
      n3s_2.setVal(initPset.n3s_2);  n3s_2.setConstant();  
      alpha1s_1.setVal(initPset.alpha1s_1);  alpha1s_1.setConstant();  
      alpha2s_1.setVal(initPset.alpha2s_1);  alpha2s_1.setConstant();  
      alpha3s_1.setVal(initPset.alpha3s_1);  alpha3s_1.setConstant();  
      alpha1s_2.setVal(initPset.alpha1s_2);  alpha1s_2.setConstant();  
      alpha2s_2.setVal(initPset.alpha2s_2);  alpha2s_2.setConstant();  
      alpha3s_2.setVal(initPset.alpha3s_2);  alpha3s_2.setConstant();  
      sigma1s_1.setVal(initPset.sigma1s_1);  sigma1s_1.setConstant();  
      sigma2s_1.setVal(initPset.sigma2s_1);  sigma2s_1.setConstant();  
      sigma3s_1.setVal(initPset.sigma3s_1);  sigma3s_1.setConstant();  
      sigma1s_2.setVal(initPset.sigma1s_2);  sigma1s_2.setConstant();  
      sigma2s_2.setVal(initPset.sigma2s_2);  sigma2s_2.setConstant();  
      sigma3s_2.setVal(initPset.sigma3s_2);  sigma3s_2.setConstant();  
      f1s->setVal(initPset.f1s);  f1s->setConstant();  
      f2s->setVal(initPset.f2s);  f2s->setConstant();  
      f3s->setVal(initPset.f3s);  f3s->setConstant();  
    } 
  }
  
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  cout << " n1s_1.getVal() = " << n1s_1.getVal() << endl;
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_1, alpha3s_1, n3s_1);
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_2, alpha3s_2, n3s_2);

  RooAddPdf*  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  RooAddPdf*  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1s) );
  RooAddPdf*  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*cb3s_2), RooArgList(*f1s) );

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",4000,0,100000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",1000,0,100000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",100, 0,10000);
  
  // background : 
  initPset.SetMCBkg();
  double init_mu = initPset.bkg_mu ;
  double init_sigma = initPset.bkg_sigma ;
  double init_lambda = initPset.bkg_lambda ;

  //  double init_mu_min = init_mu - 5; double init_mu_max = init_mu + 5;
  //  double init_sigma_min = init_sigma - 2.; double init_sigma_max = init_sigma + 2;
  //  double init_lambda_min = init_lambda - 10; double init_lambda_max = init_lambda + 10;
  double init_mu_min = init_mu - 10; double init_mu_max = init_mu + 10;
  double init_sigma_min = init_sigma - 10.; double init_sigma_max = init_sigma + 10;
  double init_lambda_min = init_lambda - 10; double init_lambda_max = init_lambda + 10;
  if(init_mu_min <0) init_mu_min = 0;
  if(init_sigma_min <0) init_sigma_min = 0;
  if(init_lambda_min <0) init_lambda_min = 0;
 
  RooRealVar err_mu("#mu","err_mu",init_mu,  0, 30) ;
  RooRealVar err_sigma("#sigma","err_sigma", init_sigma, 0,30);
  RooRealVar m_lambda("#lambda","m_lambda",  init_lambda, 0,30);
  /* 2-4pT PbPb
  RooRealVar err_mu("#mu","err_mu",init_mu,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", init_sigma, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  5, 0,23);
  //RooRealVar m_lambda("#lambda","m_lambda",  init_lambda, 0,25);
  */ 
  /*
  RooRealVar err_mu("#mu","err_mu",init_mu,  init_mu_min, init_mu_max ) ;
  RooRealVar err_sigma("#sigma","err_sigma", init_sigma, init_sigma_min, init_sigma_max);
  RooRealVar m_lambda("#lambda","m_lambda",  init_lambda, init_lambda_min, init_lambda_max);
  */
  RooGenericPdf *bkg;
  RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
  RooGenericPdf *bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  
  if  (ptLow >= 5)        bkg = bkgHighPt ;
  else bkg = bkgLowPt;

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000);  

  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));

  ws->import(*model);


  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2);
 
  
//  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb1s)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kRed),LineStyle(kDashed));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*bkg)),LineColor(kBlack),LineStyle(kDashed));

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
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
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
  //RooArgSet* paramSet = ws->pdf("model")->getParameters(*reducedDS);
  //paramSet->Print("v"); 
  RooArgList paramList = fitRes2->floatParsFinal();
  paramList.Print("v");
  //ws->pdf("model")->paramOn(legFrame,Layout(0,.95, .97));
  ws->pdf("model")->paramOn(legFrame,Layout(0,.95, .97),Parameters(paramList));
  legFrame->getAttText()->SetTextAlign(11);
  legFrame->getAttText()->SetTextSize(0.09);
  
  TPaveText* hh = (TPaveText*)legFrame->findObject(Form("%s_paramBox",ws->pdf("model")->GetName()));
  hh->SetY1(0.01); hh->SetY2(0.95);
  hh->Draw();
  //legFrame->findObject(Form("%s_paramBox",ws->pdf("model")->GetName()))->Draw();
              
  c1->SaveAs(Form("fitresults_upsilon_%sCB_%s.png",SignalCB.Data(),kineLabel.Data()));
  
  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");
  //  outh->GetXaxis()->SetBinLabel(4,"2S/1S");
  //  outh->GetXaxis()->SetBinLabel(5,"3S/1S");
  
  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  float temp2 = ws->var("nSig2s")->getVal();  
  float temp2err = ws->var("nSig2s")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();  
  
  outh->SetBinContent(1,  temp1 ) ;
  outh->SetBinError  (1,  temp1err ) ;
  outh->SetBinContent(2,  temp2 ) ;
  outh->SetBinError  (2,  temp2err ) ;
  outh->SetBinContent(3,  temp3 ) ;
  outh->SetBinError  (3,  temp3err ) ;

  cout << "1S signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;
  cout << "2S signal    =  " << outh->GetBinContent(2) << " +/- " << outh->GetBinError(2) << endl;
  cout << "3S signal    =  " << outh->GetBinContent(3) << " +/- " << outh->GetBinError(3) << endl;
    
  TFile* outf = new TFile(Form("fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()),"recreate");
  outh->Write();
  c1->Write();
  ws->Write();
  outf->Close();


  ///  cout parameters :
  /*
    cout << "N, alpha, sigma1s, M0, f, X double CB for data " << endl;
    cout << "if ( (muPtCut==(float)"<< muPtCut<<") &&  ( ptLow == (float)"<< ptLow <<" ) && (ptHigh == (float)"<<ptHigh<<" ) && (yLow == (float)"<<yLow<<" ) && (yHigh == (float)"<<yHigh<<" ) )" << endl;
    cout << " {ret.setParMC( " ;
    cout <<  ws->var("n1S")->getVal() << ", " <<  ws->var("alpha1S")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " << endl;
    cout <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;
  */
} 
 
