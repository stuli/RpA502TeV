//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the background shape.

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

#include "RooMyPdf.h"
#include "RooMyPdfPP.h"


using namespace std;
using namespace RooFit;
void FitData( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=0.0, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       int whichModel=0,   // Nominal = 0. Alternative = 1. Chebychev = 2. Power Law = 3.
	   vector<double>* resultVector = nullptr,
	   RooDataSet* pseudoData = nullptr
			) 
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  float eta_low = -2.4;
  float eta_high = 2.4;
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8;
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  TFile* f1;
  TFile* f2;
  float yLowLab;
  float yHighLab;

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString kineCut;
  
  gROOT->ProcessLine(".L RooMyPdf.cxx+");
  gROOT->ProcessLine(".L RooMyPdfPP.cxx+");
  
  TString modelLabel = "";
  
  if (whichModel == 1)
    modelLabel = "alt";
  else if (whichModel == 2)
	modelLabel = "cheb";
  else if (whichModel == 3)
	modelLabel = "pow";
  else
	modelLabel = "nom";

  //Select Data Set
  if (collId==kPADATA) {
    f1 = new TFile("../../../yskimPA1st_OpSign_20177262037_unIdentified.root");
    f2 = new TFile("../../../yskimPA2nd_OpSign_20177262044_unIdentified.root");
    yLowLab = yLow+0.47;
    yHighLab = yHigh+0.47;
    kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }
  else if (collId==kPPDATA) {
    f1 = new TFile("../../../yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20177262158_.root");
    yLowLab = yLow;
    yHighLab = yHigh;
    kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }


  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);

  //import and merge datasets
  RooWorkspace *ws = new RooWorkspace("workspace");
  RooDataSet *dataset;
  if (pseudoData == nullptr)
  {
	dataset = (RooDataSet*)f1->Get("dataset");
    if (collId==kPADATA) {
      RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
      dataset->append(*dataset2);
    }
    ws->import(*dataset);
  }
  cout << "####################################" << endl;
  RooDataSet *reducedDS;
  if (pseudoData == nullptr)
	reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  else {
	reducedDS = pseudoData;
  }
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();
  /*RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  if (collId==kPADATA) {
    RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
    dataset->append(*dataset2);
  }
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();*/

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );

  //SIGNAL:
  double sigma1s_1_init = 0.08;
  double x1s_init = 0.5;
  double alpha1s_1_init = 1.5;
  double n1s_1_init = 2.5;
  double f1s_init = 0.9;
  /*if (1) {
    sigma1s_1_init = 0.3;
    x1s_init = 0.3;
    alpha1s_1_init = 2.6;
    n1s_1_init = 3.0;
    f1s_init = 0.1;
  }*/
  //if (whichModel != 0 || pseudoData != nullptr) {
    TString NomFileName = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    sigma1s_1_init = Nomws->var("sigma1s_1")->getVal();
    x1s_init = Nomws->var("x1s")->getVal();
    alpha1s_1_init = Nomws->var("alpha1s_1")->getVal();
    n1s_1_init = Nomws->var("n1s_1")->getVal();
    f1s_init = Nomws->var("f1s")->getVal();
	
	NomFile->Close();
	delete NomFile;
	delete Nomws;
  //}


  //From Jaebeom's code:
  //PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow+0.47, yHigh+0.47, cLow, cHigh, muPtCut) ; 
  //initPset.SetMCSgl();

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", x1s_init, 0, 1);

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", alpha1s_1_init, 1.0, 3.321);
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", n1s_1_init , 1.416, 3.357);
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );
  
  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", f1s_init, 0, 1);
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(*f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(*f1s) );
  
  //if (whichModel != 0 || pseudoData != 0)
  //{
	sigma1s_1.setConstant(kTRUE);
	x1s->setConstant(kTRUE);
	alpha1s_1.setConstant(kTRUE);
	n1s_1.setConstant(kTRUE);
	f1s->setConstant(kTRUE);
  //}

  // From Jaebeom's code: Set initial parameters
  /*if ( initPset.n1s_1 == -1 )
    {
      cout << endl << endl << endl << "#########################  ERROR!!!! ##################" << endl;
      cout << "No Param. set for " << kineLabel << ","<<endl;
      cout << "Fitting macro is stopped!" << endl << endl << endl;
      return;
    }
  else { 
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << endl << "Setting the initial  parameters..." << endl << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "initPset.n1s_1 = " << initPset.n1s_1 << endl;
    n1s_1.setVal(initPset.n1s_1);  
    cout << "initPset.alpha1s_1 = " << initPset.alpha1s_1 << endl;
    alpha1s_1.setVal(initPset.alpha1s_1);
    cout << "initPset.sigma1s_1 = " << initPset.sigma1s_1 << endl;
    sigma1s_1.setVal(initPset.sigma1s_1);
    cout << "initPset.f1s = " << initPset.f1s << endl;
    f1s->setVal(initPset.f1s); 
    cout << "initPset.x1s = " << initPset.x1s << endl;
    x1s->setVal(initPset.x1s);
  } */

  // Set up crystal ball shapes
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_1, alpha3s_1, n3s_1);

  RooAddPdf* cb1s;
  RooAddPdf* cb2s;
  RooAddPdf* cb3s;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_2, alpha3s_2, n3s_2);
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1s) );
  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*cb3s_2), RooArgList(*f1s) );


  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",-20,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",-50,260000);

  //BACKGROUND
  //From Jaebeom's code:
  /*initPset.SetMCBkg();
  double err_mu_init = initPset.bkg_mu ;
  double err_sigma_init = initPset.bkg_sigma ;
  double m_lambda_init = initPset.bkg_lambda ;
  */
  /*
  double err_mu_init = 8;
  double err_sigma_init = 8;
  double m_lambda_init = 8;
  if (whichModel) {
    TString NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
    cout << NomFileName << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    if (ptLow<5) {
      err_mu_init = Nomws->var("#mu")->getVal();
      err_sigma_init = Nomws->var("#sigma")->getVal();
    }
    m_lambda_init = Nomws->var("#lambda")->getVal();
  }

  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);*/
  
  RooAbsPdf* bkg;
  
  //LOW PT ALTERNATE BACKGROUND
  RooRealVar a1("A1","A1",0.25,0,1);
  RooRealVar a2("A2","A2",0.25,0,1);
  RooRealVar a3("A3","A3",0.25,0,1);
  //RooRealVar a4("A4","A4",0.25,0,1);
  /*RooRealVar a1("A1","A1",1000,0,10000);
  RooRealVar a2("A2","A2",1000,0,10000);
  RooRealVar a3("A3","A3",1000,0,10000);
  RooRealVar a4("A4","A4",1000,0,10000);*/
  /*RooRealVar a5("A5","A5",1000,0,10000);
  a5.setConstant(kTRUE); //using only 4 bins
  RooMyPdf* bkgLowPtAlt = new RooMyPdf("bkgLowPtAlt","Background",*(ws->var("mass")),a1,a2,a3,a4,a5);
  RooMyPdfPP* bkgLowPtAltPP = new RooMyPdfPP("bkgLowPtAltPP","Background",*(ws->var("mass")),a1,a2,a3,a4,a5);*/
  TFile* lowPtBkgFile;
  if (collId == kPADATA)
	lowPtBkgFile = new TFile("altBkgModels.root","READ");
  else if (collId == kPPDATA)
	lowPtBkgFile = new TFile("altBkgModelsPP.root","READ");
  RooWorkspace* bkgWs = (RooWorkspace*)lowPtBkgFile->Get("altBkgWorkspace");
  //RooAbsPdf* bkgLowPtAlt = bkgWs->pdf("bkgLinCom");
  RooAbsPdf* bkgSumExpErf[5];
  for (int i = 1; i <= 4; i++)
  {
	bkgSumExpErf[i] = bkgWs->pdf(Form("bkgSumExpErf_%d",i));
  }
  RooAddPdf* bkgLowPtAlt = new RooAddPdf("bkgLowPtAlt","Background", RooArgList(*bkgSumExpErf[1],*bkgSumExpErf[2],*bkgSumExpErf[3],*bkgSumExpErf[4]), RooArgList(a1,a2,a3/*,a4*/));
  //ws->import(*bkgLowPtAlt);
  
  //delete bkgWs;
  //delete lowPtBkgFile;
  
  //CHEBYCHEV
  RooRealVar ach1("Ach1","Acheb1",0,-1,1);
  RooRealVar ach2("Ach2","Acheb2",-0.1,-1,1);
  RooRealVar ach3("Ach3","Acheb3",0,-1,1);
  RooRealVar ach4("Ach4","Acheb4",0,-1,1);
  RooChebychev* bkgCheb = new RooChebychev("bkgChebychev","Background",*(ws->var("mass")),RooArgList(ach1,ach2,ach3,ach4));
  
  //POWER LAW
  RooRealVar amp("amp","amplitude",200,0,10000);
  RooRealVar m0("m0","m0",1,0,100);
  RooRealVar pow("pow","pow",10,0,100);
  RooRealVar mpow("mpow","mpow",0,0,100);
  RooGenericPdf* bkgPow = new RooGenericPdf("bkgPow","Background","@1*TMath::Power(@0,@4)/TMath::Power(1+@0/@2,@3)",RooArgList(*(ws->var("mass")),amp,m0,pow,mpow));
  
  //NOMINAL BACKGROUND
  double err_mu_init = 8;
  double err_sigma_init = 8;
  double m_lambda_init = 8;
  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);
  RooGenericPdf* bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  RooGenericPdf* bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
  
  if (whichModel == 1)
  {
	bkg = bkgLowPtAlt;
  }
  else if (whichModel == 2)
	bkg = bkgCheb;
  else if (whichModel == 3)
	bkg = bkgPow;
  else
  {
	  if (ptLow >= 5)
		bkg = bkgHighPt;
	  else
		bkg = bkgLowPt;
  }

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000);  

  //Build the model
  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));

  ws->import(*model);

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //Fit the model to the data
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.054);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.075;
  float text_size = 19;
  int text_color = 1;
  if(ptLow==0 && ptHigh!=2.5) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 2.5 && ptHigh==5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 0 && ptHigh==2.5) drawText(Form("p_{T}^{#mu#mu} < %.1f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if (collId==kPPDATA) {
    if(yLow==0) drawText(Form("|y^{#mu#mu}_{cm}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < |y^{#mu#mu}_{cm}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pp
    }
  else drawText(Form("%.2f < y^{#mu#mu}_{cm} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pPb
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S)
  {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }

  TLegend* fitleg = new TLegend(0.76,0.4,0.91,0.7); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg->AddEntry(myPlot2->findObject("Sig1S"),"Signal","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");

  // PULL
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.30);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad1->SetLeftMargin(0.18);
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  pad2->SetLeftMargin(0.18);
  pad2->SetTicks(1,1);
  pad2->cd();
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.43) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.18) ; //19
  pullFrame->GetYaxis()->SetLabelSize(0.113) ; // 113
  pullFrame->GetYaxis()->SetRangeUser(-3.8,3.8) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.20) ; //23
  pullFrame->GetXaxis()->SetTitleSize(0.25) ;  //28
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  //calculate chi-squared
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
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq/ndf << endl;
  
  RooRealVar chisqVar("chisq","chi squared",chisq);
  RooRealVar ndfVar("ndf","ndf",ndf);
  ws->import(chisqVar);
  ws->import(ndfVar);

  //continue beautifying the plot and print out results
  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  pad1->Update();

  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  TString label;
  label="";
  if(collId == kPPDATA) CMS_lumi(pad1, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(pad1, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(pad1, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(pad1, 21 ,33);

  pad1->Update();
  pad2->Update();

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();
  
  
  if (resultVector == nullptr)
  {
  c1->SaveAs(Form("ResultsBkg/") + modelLabel + Form("fitresults_upsilon_%s.png",kineLabel.Data()));
  c1->SaveAs(Form("ResultsBkg/") + modelLabel + Form("fitresults_upsilon_%s.pdf",kineLabel.Data()));
  
  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");

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

	cout << "if ( binMatched( "<<muPtCut<<",  " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ; 
  cout << "  { setSignalParMC( " ;
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("alpha1s_1")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " ;
  cout <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;

  TString outFileName = Form("ResultsBkg/") + modelLabel + Form("fitresults_upsilon_%s.root",kineLabel.Data());
	
  TFile* outf = new TFile(outFileName,"recreate");
  outh->Write();
  c1->Write();
  ws->Write();
  outf->Close();
  }
  else
  {
  (*resultVector)[0] = chisq/ndf;
  (*resultVector)[1] = ws->var("nSig1s")->getVal();
  (*resultVector)[2] = ws->var("nSig2s")->getVal();
  (*resultVector)[3] = ws->var("nSig3s")->getVal();
  /*
  RooRealVar nSig1sOut("nSig1s"," 1S signals",ws->var("nSig1s")->getVal());
  RooRealVar nSig2sOut("nSig2s"," 2S signals",ws->var("nSig2s")->getVal());
  RooRealVar nSig3sOut("nSig3s"," 3S signals",ws->var("nSig3s")->getVal());
  
  mws->import(nSig1sOut);
  mws->import(nSig2sOut);
  mws->import(nSig3sOut);
  mws->import(chisqndf);*/
  
  //delete c1;
  //delete pad1;
  //delete myPlot;
  //delete model;
  //delete myPlot2;
  delete fitRes2;
  //delete fitleg;
  //delete pad2;
  //delete hpull;
  //delete pullFrame;
  //delete l1;
  
  f1->Close();
  delete f1;
  if (collId==kPADATA)
  {
	f2->Close();
	delete f2;
  }
  delete ws;
  
  delete bkgWs;
  delete lowPtBkgFile;
  }
  
  //Clean up
  /*f1->Close();
  delete f1;
  if (collId==kPADATA)
  {
	f2->Close();
	delete f2;
  }
  delete ws;
  
  delete bkgWs;
  delete lowPtBkgFile;*/
} 
 
