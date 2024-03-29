//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

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
using namespace RooFit;
void FitDataWithConstraintsAltBkgNewConstraints( 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-2.87, float yHigh=-1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
	   int whichModel=0,   // Nominal = 0. Alternative = 1. Chebychev = 2. Power Law = 3.
	   vector<double>* resultVector = nullptr,
	   RooDataSet* pseudoData = nullptr,
	   TString passedFileName = "default"
			) 
{

  TString directory = "FitsWithConstraints/";

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

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 4.5, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.5, 0.0, 0.0, 0.0, 0.0};
  
  TString modelLabel = "default";
  if (whichModel==0)
	modelLabel = "nom";
  else if (whichModel==1)
	modelLabel = "alt";
  else if (whichModel==2)
	modelLabel = "cheb";
  else if (whichModel==3)
	modelLabel = "pow";
	

  TString kineCut;

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
    //kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }


  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);

  //import and merge datasets
  RooWorkspace *ws = new RooWorkspace("workspace");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  if (pseudoData == nullptr)
  {
	dataset = (RooDataSet*)f1->Get("dataset");
    if (collId==kPADATA) {
      RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
      dataset->append(*dataset2);
	  delete dataset2;
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
  delete dataset;
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  //import ICs
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("../../../JaredNomFitsConstrained/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << "Passed file name = " << passedFileName << endl;
  if (!passedFileName.Contains("default"))
  {
    cout << "Using passed nom file name" << endl;
	NomFileName = passedFileName;
  }
  else
    cout << "Using default nom file name" << endl;
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");
  cout << "Got nominal workspace" << endl;

  float ups1smass = Nomws->var("m_{#Upsilon(1S)}")->getVal();
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",ups1smass, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
  cout << "Initialized masses and ratios" << endl;

  //SIGNAL:
  double sigma1s_1_init = Nomws->var("sigma1s_1")->getVal();
  if ((sigma1s_1_init>=0.2) || (sigma1s_1_init<=0)) sigma1s_1_init=0.05;
  double x1s_init = Nomws->var("x1s")->getVal();
  if ((x1s_init>=1.0) || (x1s_init<=0)) x1s_init=0.5;
  double alpha1s_1_init = Nomws->var("alpha1s_1")->getVal();
  if ((alpha1s_1_init>=5) || (alpha1s_1_init<=1)) alpha1s_1_init=2.5;
  double n1s_1_init = Nomws->var("n1s_1")->getVal();
  if ((n1s_1_init>=paramsupper[3]) || (n1s_1_init<=paramslower[3])) n1s_1_init=2.5;
  double f1s_init = Nomws->var("f1s")->getVal();
  if ((f1s_init>=1.0) || (f1s_init<=0)) f1s_init=0.5;
  
  cout << "Got nominal initial seeds" << endl;

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, paramslower[0], paramsupper[0]);
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", x1s_init, paramslower[1], paramsupper[1]);

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", alpha1s_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", n1s_1_init, paramslower[3], paramsupper[3]);
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", f1s_init, paramslower[4], paramsupper[4]);
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(*f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(*f1s) );

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

  double nSig1s_init = Nomws->var("nSig1s")->getVal();
  double nSig2s_init = Nomws->var("nSig2s")->getVal();
  double nSig3s_init = Nomws->var("nSig3s")->getVal();
  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",nSig1s_init,0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",nSig2s_init,-20,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",nSig3s_init,-50,260000);
  
  cout << "Set up signal pdfs and variables" << endl;

  //BACKGROUND
  RooAbsPdf* bkg;
  
  //LOW PT ALTERNATE BACKGROUND
  RooRealVar a1("A1","A1",1000,0,1000000);
  RooRealVar a2("A2","A2",1000,0,1000000);
  RooRealVar a3("A3","A3",1000,0,1000000);
  RooRealVar a4("A4","A4",1000,0,1000000);
  /*RooRealVar a5("A5","A5",1000,0,10000);
  a5.setConstant(kTRUE); //using only 4 bins*/
  TFile* lowPtBkgFile;
  if (collId == kPADATA)
	lowPtBkgFile = new TFile("altBkgModels.root","READ");
  else if (collId == kPPDATA)
	lowPtBkgFile = new TFile("altBkgModelsPP.root","READ");
  RooWorkspace* bkgWs = (RooWorkspace*)lowPtBkgFile->Get("altBkgWorkspace");
  RooAbsPdf* bkgSumExpErf[5];
  for (int i = 1; i <= 4; i++)
  {
	bkgSumExpErf[i] = bkgWs->pdf(Form("bkgSumExpErf_%d",i));
  }
  
  //delete bkgWs;
  //delete lowPtBkgFile;
  
  //CHEBYCHEV
  RooRealVar ach1("Ach1","Acheb1",0,-1,1);
  RooRealVar ach2("Ach2","Acheb2",-0.1,-1,1);
  RooRealVar ach3("Ach3","Acheb3",0,-1,1);
  RooRealVar ach4("Ach4","Acheb4",0,-1,1);
  RooChebychev* bkgCheb = new RooChebychev("bkgChebychev","Background",*(ws->var("mass")),RooArgList(ach1,ach2,ach3,ach4));
  
  //POWER LAW
  RooRealVar m0("m0","m0",1,0,100);
  RooRealVar pow("pow","pow",10,0,100);
  RooRealVar mpow("mpow","mpow",0,0,100);
  RooGenericPdf* bkgPow = new RooGenericPdf("bkgPow","Background","TMath::Power(@0,@3)/TMath::Power(1+@0/@1,@2)",RooArgList(*(ws->var("mass")),m0,pow,mpow));
  
  //NOMINAL BACKGROUND
  double err_mu_init = 8;
  double err_sigma_init = 1;
  if (!(ptLow >= 5.0 || (ptLow == 4.0 && ptHigh == 9.0)))
  {
	cout << "Getting nom bkg param seeds from nom" << endl;
	err_mu_init = Nomws->var("#mu")->getVal();
	err_sigma_init = Nomws->var("#sigma")->getVal();
  }
  double m_lambda_init = Nomws->var("#lambda")->getVal();
  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 9) ; //upper limit was 15
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,15);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);
  RooGenericPdf* bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  RooGenericPdf* bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
  
  RooAbsReal* nBkg;
  RooAddPdf* model;
  
  cout << "Made possible bkg functions" << endl;

  // Construct Gaussian constraints on parameters
  float nfix, xfix, alphafix, ffix;
  float ndev, xdev, alphadev, fdev;
  if (collId==kPADATA) {
    nfix = 3.20859; //2.90258;
    xfix = 0.562105; //0.582016;
    alphafix = 4.37737; //1.95596;
    ndev = 1.97516; //1.51273;
    xdev =  0.0544103; //0.284584;
    alphadev = 0.534965; //0.611055;
  }
  else if (collId==kPPDATA) {
    nfix = 1.42638; //2.29348;
    xfix = 0.556545; //0.542762;
    alphafix = 1.97426; //2.48620;
    ndev = 0.0518141; //0.451840;
    xdev = 0.0096801; //0.0585783;
    alphadev = 0.0205957; //0.395293;
  }
  //n1s_1.setVal(nfix);
  //x1s->setVal(xfix);
  //alpha1s_1.setVal(alphafix);
  RooGaussian nconstraint("nconstraint","nconstraint", n1s_1,RooConst(nfix),RooConst(ndev));
  RooGaussian alphaconstraint("alphaconstraint","alphaconstraint", alpha1s_1,RooConst(alphafix),RooConst(alphadev));
  RooGaussian xconstraint("xconstraint","xconstraint", *x1s,RooConst(xfix),RooConst(xdev));

  RooArgSet allConstraints;
  
  if (collId==kPADATA) {
    float yLowpp;
    float yHighpp;
    if (yLow<0) {
      yLowpp = abs(yHigh);
      yHighpp = abs(yLow);
      if (yHigh>0) {
        yLowpp = 0.0;
        yHighpp = 1.93;
      }
      if (yLow<-2.8 && yHigh<0) {
        yLowpp = 1.2;
        yHighpp = 1.93;
      }
    }
    else {
      yLowpp = yLow;
      yHighpp = yHigh;
    }
    float ffix = Nomws->var("f1s")->getVal();
	float fdev = Nomws->var("f1s")->getError();
	TString kineLabelpp = getKineLabel (kPPDATA, ptLow, ptHigh, yLowpp, yHighpp, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString ppFileName = Form("../../../JaredNomFitsConstrained/nomfitresults_upsilon_%s.root",kineLabelpp.Data());
    cout << ppFileName << endl;
    //TFile* ppFile = new TFile(ppFileName,"READ");
	if (yLow < -2 || NomFileName.Contains("hfsum"))
	{
		cout << "Using nom PA for f" << endl;
	}
	else
	{
		TFile* ppFile = new TFile(ppFileName,"READ");
		cout << "Using nom PP for f" << endl;
		RooWorkspace *ppws = (RooWorkspace*)ppFile->Get("workspace");
		ppFile->Close("R");
		ffix = ppws->var("f1s")->getVal();
		fdev = ppws->var("f1s")->getError();
		delete ppws;
		delete ppFile;
	}
    f1s->setVal(ffix);
    RooGaussian fconstraint("fconstraint","fconstraint", *f1s,RooConst(ffix),RooConst(fdev));
    //delete ppws;
    //delete ppFile;
	
	allConstraints = RooArgSet(nconstraint,alphaconstraint,xconstraint,fconstraint);
  }
  else if (collId==kPPDATA) allConstraints = RooArgSet(nconstraint,alphaconstraint,xconstraint);
  /*
  RooArgSet allConstraints;
  if (collId==kPADATA) allConstraints = RooArgSet(nconstraint,alphaconstraint,xconstraint,fconstraint);
  else if (collId==kPPDATA) allConstraints = RooArgSet(nconstraint,alphaconstraint,xconstraint);
  */
  
  cout << "Set f constraint" << endl;
  
  //Build the model
  if (whichModel == 1)
  {
	//bkg = bkgLowPtAlt;
	//nBkg = new RooFormulaVar("nBkg","fraction of component 1 in bkg","@0+@1+@2+@3",RooArgList(a1,a2,a3,a4));
	model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkgSumExpErf[1],*bkgSumExpErf[2],*bkgSumExpErf[3],*bkgSumExpErf[4]),RooArgList(*nSig1s,*nSig2s,*nSig3s,a1,a2,a3,a4));
  }
  else
  {
	if (whichModel == 2)
		bkg = bkgCheb;
	else if (whichModel == 3)
		bkg = bkgPow;
	else
	{
		if (ptLow >= 5.0 || (ptLow == 4.0 && ptHigh == 9.0))
			bkg = bkgHighPt;
		else
			bkg = bkgLowPt;
	}
	nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000);
	model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  }
  
  //RooAddPdf* model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));

  ws->import(*model);

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //Fit the model to the data
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,ExternalConstraints(allConstraints),Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  if (whichModel==1)
	ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkgSumExpErf[1],*bkgSumExpErf[2],*bkgSumExpErf[3],*bkgSumExpErf[4])),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  else
    ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  ws->import(*fitRes2);

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
  cout << "f1s_init = " << f1s_init << endl;
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.075;
  float text_size = 19;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if (collId==kPPDATA) {
    if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  else if (collId==kPADATA) {
    if(yLow==-yHigh) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);

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
  cout << "chisq/dof = " << chisq << "/" << ndf << " = " << chisq/ndf << endl;
  
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
  
  if (resultVector != nullptr)
  {
	  (*resultVector)[0] = chisq/ndf;
	  (*resultVector)[1] = ws->var("nSig1s")->getVal();
	  (*resultVector)[2] = ws->var("nSig2s")->getVal();
	  (*resultVector)[3] = ws->var("nSig3s")->getVal();
  }
  cout << "Made result vector" << endl;

  if (pseudoData == nullptr)
  {
  
  c1->SaveAs(Form("%s%sfitresults_upsilon_%s.png",directory.Data(),modelLabel.Data(),kineLabel.Data()));
  c1->SaveAs(Form("%s%sfitresults_upsilon_%s.pdf",directory.Data(),modelLabel.Data(),kineLabel.Data()));

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

  delete NomFile;
  delete Nomws;

TString outFileName;
  outFileName = Form("%s%sfitresults_upsilon_%s.root",directory.Data(),modelLabel.Data(),kineLabel.Data());

  TFile* outf = new TFile(outFileName,"recreate");
  outh->Write();
  c1->Write();
  ws->Write();
  outf->Close();
  
  if (resultVector != nullptr)
  {
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

  delete outh;
  delete outf;
  
  }
  else
  {
  cout << "Cleaning up fitting objects..." << endl;
  delete f1;
  if (collId==kPADATA) delete f2;
  delete ws;
  delete f1s;
  delete x1s;
  delete cb1s;
  delete cb2s;
  delete cb3s;
  delete nSig1s;
  delete nSig2s;
  delete nSig3s;
  //delete NomFile;
  //delete Nomws;
  if (whichModel != 1)
	delete nBkg;
  delete bkgLowPt;
  delete bkgHighPt;
 // delete bkgCheb;
  //delete bkgPow;
  //delete bkgWs;
  //lowPtBkgFile->Close();
  //delete lowPtBkgFile;
  delete model;
  cout << "here1" << endl;
  delete fitRes2;
  cout << "here2" << endl;
  //delete reducedDS;
  cout << "here3" << endl;
  delete l1;
  cout << "here4" << endl;
  //delete myPlot;
  cout << "here5" << endl;
  delete myPlot2;
  cout << "here6" << endl;
  delete fitleg;
  cout << "here7" << endl;
  delete pad2;
  cout << "here8" << endl;
  delete pad1;
  cout << "here9" << endl;
  delete c1;
  cout << "here10" << endl;
  }
  
  return;
} 
 
