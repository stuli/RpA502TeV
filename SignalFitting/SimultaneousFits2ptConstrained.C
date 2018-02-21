//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "../HeaderFiles/rootFitHeaders.h"
#include "../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "../HeaderFiles/PsetCollection.h"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/StyleSetting.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"


using namespace std;
using namespace RooFit;
void SimultaneousFits2ptConstrained( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0,   // Nominal = 0. Alternative = 1.
       int ICset=1
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

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 3.0, 3.321, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  TString kineCut_A, kineCut_B;

  float ptLow_A, ptHigh_A, ptLow_B, ptHigh_B;
  ptLow_A = 0;
  ptHigh_A = 6;
  ptLow_B = 6;
  ptHigh_B = 30;

  cout << endl << "PERFORMING SIMULTANEOUS FIT IN PT BINS [" << ptLow_A << "," << ptHigh_A << "], [" << ptLow_B << "," << ptHigh_B << "]" << endl << endl;

  //Select Data Set
  if (collId==kPADATA) {
    f1 = new TFile("../yskimPA1st_OpSign_20177262037_unIdentified.root");
    f2 = new TFile("../yskimPA2nd_OpSign_20177262044_unIdentified.root");
    yLowLab = yLow+0.47;
    yHighLab = yHigh+0.47;
    kineCut_A = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_A, ptHigh_A, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_B = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_B, ptHigh_B, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }
  else if (collId==kPPDATA) {
    f1 = new TFile("../yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20177262158_.root");
    yLowLab = yLow;
    yHighLab = yHigh;
    kineCut_A = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_A, ptHigh_A, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_B = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_B, ptHigh_B, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }

  if (muPtCut>0){
    kineCut_A = kineCut_A + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
    kineCut_B = kineCut_B + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
  }

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  if (collId==kPADATA) {
    RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
    dataset->append(*dataset2);
  }
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS_A = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_A.Data() );
  reducedDS_A->SetName("reducedDS_A");
  ws->import(*reducedDS_A);
  RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_B.Data() );
  reducedDS_B->SetName("reducedDS_B");
  ws->import(*reducedDS_B);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  RooCategory tp("tp","tp");
  tp.defineType("A");
  tp.defineType("B");
  tp.defineType("C");

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  RooDataSet* dsABC = new RooDataSet("dsABC","dsABC",RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS_A),Import("B",*reducedDS_B));
  cout << "******** New Combined Dataset ***********" << endl;
  dsABC->Print();
  ws->import(*dsABC);

  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );

  //import info from full bin
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("SimultaneousFits/Normal/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");
  double totSig1 = Nomws->var("nSig1s")->getVal();
  double totSig2 = Nomws->var("nSig2s")->getVal();
  double totSig3 = Nomws->var("nSig3s")->getVal();
  RooRealVar nSig1s_C("nSig1s_C","1S signals",totSig1);
  RooRealVar nSig2s_C("nSig2s_C","1S signals",totSig2);
  RooRealVar nSig3s_C("nSig3s_C","1S signals",totSig3);

  TCanvas* c_A =  new TCanvas("canvas_A","My plots",4,45,550,520);
  c_A->cd();
  RooPlot* myPlot_A = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_A")->plotOn(myPlot_A,Name("dataHist_A"));

  RooRealVar mean1s_A("m_{#Upsilon(1S)}_A","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_A("mean2s_A","m_{#Upsilon(1S)}_A*mRatio21", RooArgSet(mean1s_A,mRatio21) );
  RooFormulaVar mean3s_A("mean3s_A","m_{#Upsilon(1S)}_A*mRatio31", RooArgSet(mean1s_A,mRatio31) );

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",504,45,550,520);
  c_B->cd();
  RooPlot* myPlot_B = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_B")->plotOn(myPlot_B,Name("dataHist_B"));

  RooRealVar mean1s_B("m_{#Upsilon(1S)}_B","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_B("mean2s_B","m_{#Upsilon(1S)}_B*mRatio21", RooArgSet(mean1s_B,mRatio21) );
  RooFormulaVar mean3s_B("mean3s_B","m_{#Upsilon(1S)}_B*mRatio31", RooArgSet(mean1s_B,mRatio31) );

  //SIGNAL:
  double sigma1s_1_init = 0.1;
  double x1s_init = 0.8;
  double alpha1s_1_init = 1.5;
  double n1s_1_init = 2.6;
  double f1s_init = 0.5;
  if (ICset>1 && ICset<4) {
    sigma1s_1_init = 0.3;
    x1s_init = 0.3;
    alpha1s_1_init = 2.6;
    n1s_1_init = 3.0;
    f1s_init = 0.1;
  }

//A:
  RooRealVar    sigma1s_1_A("sigma1s_1_A","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, paramslower[0], paramsupper[0]);
  RooFormulaVar sigma2s_1_A("sigma2s_1_A","@0*@1",RooArgList(sigma1s_1_A,mRatio21) );
  RooFormulaVar sigma3s_1_A("sigma3s_1_A","@0*@1",RooArgList(sigma1s_1_A,mRatio31) );

  RooRealVar *x1s_A = new RooRealVar("x1s_A","sigma ratio ", x1s_init, paramslower[1], paramsupper[1]);

  RooFormulaVar sigma1s_2_A("sigma1s_2_A","@0*@1",RooArgList(sigma1s_1_A, *x1s_A) );
  RooFormulaVar sigma2s_2_A("sigma2s_2_A","@0*@1",RooArgList(sigma1s_2_A,mRatio21) );
  RooFormulaVar sigma3s_2_A("sigma3s_2_A","@0*@1",RooArgList(sigma1s_2_A,mRatio31) );

  RooRealVar alpha1s_1_A("alpha1s_1_A","tail shift", alpha1s_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha2s_1_A("alpha2s_1_A","1.0*@0",RooArgList(alpha1s_1_A) );
  RooFormulaVar alpha3s_1_A("alpha3s_1_A","1.0*@0",RooArgList(alpha1s_1_A) );
  RooFormulaVar alpha1s_2_A("alpha1s_2_A","1.0*@0",RooArgList(alpha1s_1_A) );
  RooFormulaVar alpha2s_2_A("alpha2s_2_A","1.0*@0",RooArgList(alpha1s_1_A) );
  RooFormulaVar alpha3s_2_A("alpha3s_2_A","1.0*@0",RooArgList(alpha1s_1_A) );

  RooRealVar n1s_1_A("n1s_1_A","power order", n1s_1_init , paramslower[3], paramsupper[3]);
  RooFormulaVar n2s_1_A("n2s_1_A","1.0*@0",RooArgList(n1s_1_A) );
  RooFormulaVar n3s_1_A("n3s_1_A","1.0*@0",RooArgList(n1s_1_A) );
  RooFormulaVar n1s_2_A("n1s_2_A","1.0*@0",RooArgList(n1s_1_A) );
  RooFormulaVar n2s_2_A("n2s_2_A","1.0*@0",RooArgList(n1s_1_A) );
  RooFormulaVar n3s_2_A("n3s_2_A","1.0*@0",RooArgList(n1s_1_A) );

  RooRealVar *f1s_A = new RooRealVar("f1s_A","1S CB fraction", f1s_init, paramslower[4], paramsupper[4]);
  RooFormulaVar f2s_A("f2s_A","1.0*@0",RooArgList(*f1s_A) );
  RooFormulaVar f3s_A("f3s_A","1.0*@0",RooArgList(*f1s_A) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1_A = new RooCBShape("cball1s_1_A", "cystal Ball", *(ws->var("mass")), mean1s_A, sigma1s_1_A, alpha1s_1_A, n1s_1_A);
  RooCBShape* cb2s_1_A = new RooCBShape("cball2s_1_A", "cystal Ball", *(ws->var("mass")), mean2s_A, sigma2s_1_A, alpha2s_1_A, n2s_1_A);
  RooCBShape* cb3s_1_A = new RooCBShape("cball3s_1_A", "cystal Ball", *(ws->var("mass")), mean3s_A, sigma3s_1_A, alpha3s_1_A, n3s_1_A);

  RooAddPdf* cb1s_A;
  RooAddPdf* cb2s_A;
  RooAddPdf* cb3s_A;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2_A = new RooCBShape("cball1s_2_A", "cystal Ball", *(ws->var("mass")), mean1s_A, sigma1s_2_A, alpha1s_2_A, n1s_2_A);
  RooCBShape* cb2s_2_A = new RooCBShape("cball2s_2_A", "cystal Ball", *(ws->var("mass")), mean2s_A, sigma2s_2_A, alpha2s_2_A, n2s_2_A);
  RooCBShape* cb3s_2_A = new RooCBShape("cball3s_2_A", "cystal Ball", *(ws->var("mass")), mean3s_A, sigma3s_2_A, alpha3s_2_A, n3s_2_A);
  cb1s_A = new RooAddPdf("cb1s_A","Signal 1S",RooArgList(*cb1s_1_A,*cb1s_2_A), RooArgList(*f1s_A) );
  cb2s_A = new RooAddPdf("cb2s_A","Signal 2S",RooArgList(*cb2s_1_A,*cb2s_2_A), RooArgList(*f1s_A) );
  cb3s_A = new RooAddPdf("cb3s_A","Signal 3S",RooArgList(*cb3s_1_A,*cb3s_2_A), RooArgList(*f1s_A) );

  RooRealVar *nSig1s_A= new RooRealVar("nSig1s_A"," 1S signals",totSig1/3,0,totSig1);
  RooRealVar *nSig2s_A= new RooRealVar("nSig2s_A"," 2S signals",totSig2/3,-20,totSig2);
  RooRealVar *nSig3s_A= new RooRealVar("nSig3s_A"," 3S signals",totSig3/3,-50,totSig3);


//B:
  RooRealVar    sigma1s_1_B("sigma1s_1_B","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, paramslower[0], paramsupper[0]);
  RooFormulaVar sigma2s_1_B("sigma2s_1_B","@0*@1",RooArgList(sigma1s_1_B,mRatio21) );
  RooFormulaVar sigma3s_1_B("sigma3s_1_B","@0*@1",RooArgList(sigma1s_1_B,mRatio31) );

  RooRealVar *x1s_B = new RooRealVar("x1s_B","sigma ratio ", x1s_init, paramslower[1], paramsupper[1]);

  RooFormulaVar sigma1s_2_B("sigma1s_2_B","@0*@1",RooArgList(sigma1s_1_B, *x1s_B) );
  RooFormulaVar sigma2s_2_B("sigma2s_2_B","@0*@1",RooArgList(sigma1s_2_B,mRatio21) );
  RooFormulaVar sigma3s_2_B("sigma3s_2_B","@0*@1",RooArgList(sigma1s_2_B,mRatio31) );

  RooRealVar alpha1s_1_B("alpha1s_1_B","tail shift", alpha1s_1_init, paramslower[2], paramsupper[2]);
  RooFormulaVar alpha2s_1_B("alpha2s_1_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha3s_1_B("alpha3s_1_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha1s_2_B("alpha1s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha2s_2_B("alpha2s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha3s_2_B("alpha3s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );

  RooRealVar n1s_1_B("n1s_1_B","power order", n1s_1_init , paramslower[3], paramsupper[3]);
  RooFormulaVar n2s_1_B("n2s_1_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n3s_1_B("n3s_1_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n1s_2_B("n1s_2_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n2s_2_B("n2s_2_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n3s_2_B("n3s_2_B","1.0*@0",RooArgList(n1s_1_B) );

  RooRealVar *f1s_B = new RooRealVar("f1s_B","1S CB fraction", f1s_init, paramslower[4], paramsupper[4]);
  RooFormulaVar f2s_B("f2s_B","1.0*@0",RooArgList(*f1s_B) );
  RooFormulaVar f3s_B("f3s_B","1.0*@0",RooArgList(*f1s_B) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1_B = new RooCBShape("cball1s_1_B", "cystal Ball", *(ws->var("mass")), mean1s_B, sigma1s_1_B, alpha1s_1_B, n1s_1_B);
  RooCBShape* cb2s_1_B = new RooCBShape("cball2s_1_B", "cystal Ball", *(ws->var("mass")), mean2s_B, sigma2s_1_B, alpha2s_1_B, n2s_1_B);
  RooCBShape* cb3s_1_B = new RooCBShape("cball3s_1_B", "cystal Ball", *(ws->var("mass")), mean3s_B, sigma3s_1_B, alpha3s_1_B, n3s_1_B);

  RooAddPdf* cb1s_B;
  RooAddPdf* cb2s_B;
  RooAddPdf* cb3s_B;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2_B = new RooCBShape("cball1s_2_B", "cystal Ball", *(ws->var("mass")), mean1s_B, sigma1s_2_B, alpha1s_2_B, n1s_2_B);
  RooCBShape* cb2s_2_B = new RooCBShape("cball2s_2_B", "cystal Ball", *(ws->var("mass")), mean2s_B, sigma2s_2_B, alpha2s_2_B, n2s_2_B);
  RooCBShape* cb3s_2_B = new RooCBShape("cball3s_2_B", "cystal Ball", *(ws->var("mass")), mean3s_B, sigma3s_2_B, alpha3s_2_B, n3s_2_B);
  cb1s_B = new RooAddPdf("cb1s_B","Signal 1S",RooArgList(*cb1s_1_B,*cb1s_2_B), RooArgList(*f1s_B) );
  cb2s_B = new RooAddPdf("cb2s_B","Signal 2S",RooArgList(*cb2s_1_B,*cb2s_2_B), RooArgList(*f1s_B) );
  cb3s_B = new RooAddPdf("cb3s_B","Signal 3S",RooArgList(*cb3s_1_B,*cb3s_2_B), RooArgList(*f1s_B) );

  //RooRealVar *nSig1s_B= new RooRealVar("nSig1s_B"," 1S signals",totSig1/3,0,totSig1);
  //RooRealVar *nSig2s_B= new RooRealVar("nSig2s_B"," 2S signals",totSig2/3,-20,totSig2);
  //RooRealVar *nSig3s_B= new RooRealVar("nSig3s_B"," 3S signals",totSig3/3,-50,totSig3);

  //BACKGROUND
  double err_sigma_init = 5;
  double err_mu_init = 8;
  double m_lambda_init = 5;
  if (ICset>2) {
    err_mu_init = 5;
    m_lambda_init = 5;
  }

//A:
  RooRealVar err_mu_A("#mu_A","err_mu", err_mu_init,  paramslower[5], paramsupper[5]) ;
  RooRealVar err_sigma_A("#sigma_A","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
  RooRealVar m_lambda_A("#lambda_A","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);

  RooGenericPdf *bkg_A;
  RooGenericPdf *bkgLowPt_A = new RooGenericPdf("bkgLowPt_A","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_A, err_mu_A, err_sigma_A) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_A = new RooGenericPdf("bkgHighPt_A","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_A));
  if  (ptLow_A >= 5)        bkg_A = bkgHighPt_A ;
  else bkg_A = bkgLowPt_A;

  RooRealVar *nBkg_A = new RooRealVar("nBkg_A","fraction of component 1 in bkg",10000,0,5000000); 

//B:
  RooRealVar err_mu_B("#mu_B","err_mu", err_mu_init,  paramslower[5], paramsupper[5]) ;
  RooRealVar err_sigma_B("#sigma_B","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
  RooRealVar m_lambda_B("#lambda_B","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);

  RooGenericPdf *bkg_B;
  RooGenericPdf *bkgLowPt_B = new RooGenericPdf("bkgLowPt_B","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_B, err_mu_B, err_sigma_B) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_B = new RooGenericPdf("bkgHighPt_B","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_B));
  if  (ptLow_B >= 5)        bkg_B = bkgHighPt_B ;
  else bkg_B = bkgLowPt_B;

  RooRealVar *nBkg_B = new RooRealVar("nBkg_B","fraction of component 1 in bkg",10000,0,5000000);   

//Constraints
  RooFormulaVar *nSig1s_B = new RooFormulaVar("nSig1s_B"," 1S signals","@0-@1",RooArgList(nSig1s_C,*nSig1s_A) );
  RooFormulaVar *nSig2s_B = new RooFormulaVar("nSig2s_B"," 2S signals","@0-@1",RooArgList(nSig2s_C,*nSig2s_A) );
  RooFormulaVar *nSig3s_B = new RooFormulaVar("nSig3s_B"," 3S signals","@0-@1",RooArgList(nSig3s_C,*nSig3s_A) );

  //Build the model
  RooAddPdf* model_A = new RooAddPdf();
  model_A = new RooAddPdf("model_A","1S+2S+3S + Bkg",RooArgList(*cb1s_A, *cb2s_A, *cb3s_A, *bkg_A),RooArgList(*nSig1s_A,*nSig2s_A,*nSig3s_A,*nBkg_A));
  ws->import(*model_A);
  c_A->cd();
  RooPlot* myPlot2_A = (RooPlot*)myPlot_A->Clone();
  ws->data("reducedDS_A")->plotOn(myPlot2_A,Name("dataOS_FIT_A"),MarkerSize(.8));

  RooAddPdf* model_B = new RooAddPdf();
  model_B = new RooAddPdf("model_B","1S+2S+3S + Bkg",RooArgList(*cb1s_B, *cb2s_B, *cb3s_B, *bkg_B),RooArgList(*nSig1s_B,*nSig2s_B,*nSig3s_B,*nBkg_B));
  ws->import(*model_B);
  c_B->cd();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  ws->data("reducedDS_B")->plotOn(myPlot2_B,Name("dataOS_FIT_B"),MarkerSize(.8));

  // Construct simultaneous PDF
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",tp);
  simPdf->addPdf(*(ws->pdf("model_A")),"A") ;
  simPdf->addPdf(*(ws->pdf("model_B")),"B") ;
  ws->import(*simPdf);

  cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("simPdf")->fitTo(*dsABC,Save(), Hesse(kTRUE),Range(massLow,massHigh),Timer(kTRUE),Extended(kTRUE));
  cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;

  c_A->cd();
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("modelHist_A"));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("Sig1S_A"),Components(RooArgSet(*cb1s_A)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_A")->plotOn(myPlot2_A,Components(RooArgSet(*cb2s_A)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_A")->plotOn(myPlot2_A,Components(RooArgSet(*cb3s_A)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_A")->plotOn(myPlot2_A,Name("bkgPDF_A"),Components(RooArgSet(*bkg_A)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_A->SetFillStyle(4000);
  myPlot2_A->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_A->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_A->GetYaxis()->CenterTitle();
  //myPlot2_A->GetYaxis()->SetTitleSize(0.058);
  myPlot2_A->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_A->GetXaxis()->SetLabelSize(0);
  myPlot2_A->GetXaxis()->SetRangeUser(8,14);
  //myPlot2_A->GetXaxis()->SetTitleSize(0);
  myPlot2_A->Draw();

  c_B->cd();
  ws->pdf("model_B")->plotOn(myPlot2_B,Name("modelHist_B"));
  ws->pdf("model_B")->plotOn(myPlot2_B,Name("Sig1S_B"),Components(RooArgSet(*cb1s_B)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_B")->plotOn(myPlot2_B,Components(RooArgSet(*cb2s_B)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_B")->plotOn(myPlot2_B,Components(RooArgSet(*cb3s_B)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_B")->plotOn(myPlot2_B,Name("bkgPDF_B"),Components(RooArgSet(*bkg_B)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_B->SetFillStyle(4000);
  myPlot2_B->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_B->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_B->GetYaxis()->CenterTitle();
  //myPlot2_B->GetYaxis()->SetTitleSize(0.058);
  myPlot2_B->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_B->GetXaxis()->SetLabelSize(0);
  myPlot2_B->GetXaxis()->SetRangeUser(8,14);
  //myPlot2_B->GetXaxis()->SetTitleSize(0);
  myPlot2_B->Draw();

  c_A->Update();
  c_B->Update();

  TString outFileName;
  if (whichModel){
    outFileName = Form("SimultaneousFits/sim_ptSplit3_altfitresults_upsilon_%s.root",kineLabel.Data());
  }
  else {
    outFileName = Form("SimultaneousFits/sim_ptSplit3_nomfitresults_upsilon_%s.root",kineLabel.Data());
  }
  TFile* outf = new TFile(outFileName,"recreate");
  ws->Write();
  outf->Close();

} 
 
