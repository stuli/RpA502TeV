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
#include "RooDataHist.h"
#include "RooCategory.h"


using namespace std;
using namespace RooFit;
void Simultaneous5Fits( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  int ICset = 1;
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

  TString kineCut, kineCut_B, kineCut_C, kineCut_D, kineCut_E;

  ptLow = 0;
  ptHigh = 30;
  float ptLow_B = 0;
  float ptHigh_B = 6;
  float yLow_B = -1.93;
  float yLow_B = 1.93;
  float ptLow_C = 6;
  float ptHigh_C = 30;
  float yLow_C = -1.93;
  float yLow_C = 1.93;
  float ptLow_D = 0;
  float ptHigh_D = 30;
  float yLow_D = -1.93;
  float yLow_D = 0.0;
  float ptLow_E = 0;
  float ptHigh_E = 30;
  float yLow_E = 0.0;
  float yLow_E = 1.93;
  //Select Data Set
  if (collId==kPADATA) {
    f1 = new TFile("../../yskimPA1st_OpSign_20177262037_unIdentified.root");
    f2 = new TFile("../../yskimPA2nd_OpSign_20177262044_unIdentified.root");
    yLowLab = yLow+0.47;
    yHighLab = yHigh+0.47;
    kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_B = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_B, ptHigh_B, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_C = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_C, ptHigh_C, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_D = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_D, ptHigh_D, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_E = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_E, ptHigh_E, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }
  else if (collId==kPPDATA) {
    f1 = new TFile("../../yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20177262158_.root");
    yLowLab = yLow;
    yHighLab = yHigh;
    kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_B = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_B, ptHigh_B, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_C = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_C, ptHigh_C, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_D = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_D, ptHigh_D, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    kineCut_E = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow_E, ptHigh_E, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  }

  if (muPtCut>0){
    kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
    kineCut_B = kineCut_B + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
    kineCut_C = kineCut_C + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
    kineCut_D = kineCut_D + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
    kineCut_E = kineCut_E + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
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
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  RooDataSet *reducedDS_B = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_B.Data() );
  reducedDS_B->SetName("reducedDS_B");
  ws->import(*reducedDS_B);
  RooDataSet *reducedDS_C = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_C.Data() );
  reducedDS_C->SetName("reducedDS_C");
  ws->import(*reducedDS_C);
  RooDataSet *reducedDS_D = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_D.Data() );
  reducedDS_D->SetName("reducedDS_D");
  ws->import(*reducedDS_D);
  RooDataSet *reducedDS_E = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut_E.Data() );
  reducedDS_E->SetName("reducedDS_E");
  ws->import(*reducedDS_E);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  RooCategory tp("tp","tp");
  tp.defineType("A");
  tp.defineType("B");
  tp.defineType("C");
  tp.defineType("D");
  tp.defineType("E");

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  RooDataSet* dsABC = new RooDataSet("dsABC","dsABC",RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("A",*reducedDS),Import("B",*reducedDS_B),Import("C",*reducedDS_C),Import("D",*reducedDS_D),Import("E",*reducedDS_E));
  cout << "******** New Combined Dataset ***********" << endl;
  dsABC->Print();
  ws->import(*dsABC);

  TCanvas* c_A =  new TCanvas("canvas","My plots",4,45,550,520);
  c_A->cd();
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );

  TCanvas* c_B =  new TCanvas("canvas_B","My plots",504,45,550,520);
  c_B->cd();
  RooPlot* myPlot_B = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_B")->plotOn(myPlot_B,Name("dataHist_B"));

  RooRealVar mean1s_B("m_{#Upsilon(1S)}_B","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_B("mean2s_B","m_{#Upsilon(1S)}_B*mRatio21", RooArgSet(mean1s_B,mRatio21) );
  RooFormulaVar mean3s_B("mean3s_B","m_{#Upsilon(1S)}_B*mRatio31", RooArgSet(mean1s_B,mRatio31) );

  TCanvas* c_C =  new TCanvas("canvas_C","My plots",1004,45,550,520);
  c_C->cd();
  RooPlot* myPlot_C = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_C")->plotOn(myPlot_C,Name("dataHist_C"));

  RooRealVar mean1s_C("m_{#Upsilon(1S)}_C","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_C("mean2s_C","m_{#Upsilon(1S)}_C*mRatio21", RooArgSet(mean1s_C,mRatio21) );
  RooFormulaVar mean3s_C("mean3s_C","m_{#Upsilon(1S)}_C*mRatio31", RooArgSet(mean1s_C,mRatio31) );

  TCanvas* c_D =  new TCanvas("canvas_D","My plots",250,245,550,520);
  c_D->cd();
  RooPlot* myPlot_D = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_D")->plotOn(myPlot_D,Name("dataHist_D"));

  RooRealVar mean1s_D("m_{#Upsilon(1S)}_D","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_D("mean2s_D","m_{#Upsilon(1S)}_D*mRatio21", RooArgSet(mean1s_D,mRatio21) );
  RooFormulaVar mean3s_D("mean3s_D","m_{#Upsilon(1S)}_D*mRatio31", RooArgSet(mean1s_D,mRatio31) );

  TCanvas* c_E =  new TCanvas("canvas_E","My plots",250,245,550,520);
  c_E->cd();
  RooPlot* myPlot_E = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS_E")->plotOn(myPlot_E,Name("dataHist_E"));

  RooRealVar mean1s_E("m_{#Upsilon(1S)}_E","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooFormulaVar mean2s_E("mean2s_E","m_{#Upsilon(1S)}_E*mRatio21", RooArgSet(mean1s_E,mRatio21) );
  RooFormulaVar mean3s_E("mean3s_E","m_{#Upsilon(1S)}_E*mRatio31", RooArgSet(mean1s_E,mRatio31) );

  //SIGNAL:
  double sigma1s_1_init = 0.08;
  double x1s_init = 0.5;
  double alpha1s_1_init = 1.5;
  double n1s_1_init = 2.5;
  double f1s_init = 0.9;
  if (ICset>1 && ICset<4) {
    sigma1s_1_init = 0.3;
    x1s_init = 0.3;
    alpha1s_1_init = 2.6;
    n1s_1_init = 3.0;
    f1s_init = 0.1;
  }

//A:
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

  //RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",0,1000000);
  //RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",-20,360000);
  //RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",-50,260000);

//B:
  RooRealVar    sigma1s_1_B("sigma1s_1_B","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  RooFormulaVar sigma2s_1_B("sigma2s_1_B","@0*@1",RooArgList(sigma1s_1_B,mRatio21) );
  RooFormulaVar sigma3s_1_B("sigma3s_1_B","@0*@1",RooArgList(sigma1s_1_B,mRatio31) );

  RooRealVar *x1s_B = new RooRealVar("x1s_B","sigma ratio ", x1s_init, 0, 1);

  RooFormulaVar sigma1s_2_B("sigma1s_2_B","@0*@1",RooArgList(sigma1s_1_B, *x1s_B) );
  RooFormulaVar sigma2s_2_B("sigma2s_2_B","@0*@1",RooArgList(sigma1s_2_B,mRatio21) );
  RooFormulaVar sigma3s_2_B("sigma3s_2_B","@0*@1",RooArgList(sigma1s_2_B,mRatio31) );

  RooRealVar alpha1s_1_B("alpha1s_1_B","tail shift", alpha1s_1_init, 1.0, 3.321);
  RooFormulaVar alpha2s_1_B("alpha2s_1_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha3s_1_B("alpha3s_1_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha1s_2_B("alpha1s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha2s_2_B("alpha2s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );
  RooFormulaVar alpha3s_2_B("alpha3s_2_B","1.0*@0",RooArgList(alpha1s_1_B) );

  RooRealVar n1s_1_B("n1s_1_B","power order", n1s_1_init , 1.416, 3.357);
  RooFormulaVar n2s_1_B("n2s_1_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n3s_1_B("n3s_1_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n1s_2_B("n1s_2_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n2s_2_B("n2s_2_B","1.0*@0",RooArgList(n1s_1_B) );
  RooFormulaVar n3s_2_B("n3s_2_B","1.0*@0",RooArgList(n1s_1_B) );

  RooRealVar *f1s_B = new RooRealVar("f1s_B","1S CB fraction", f1s_init, 0, 1);
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

  RooRealVar *nSig1s_B= new RooRealVar("nSig1s_B"," 1S signals",0,1000000);
  RooRealVar *nSig2s_B= new RooRealVar("nSig2s_B"," 2S signals",-20,360000);
  RooRealVar *nSig3s_B= new RooRealVar("nSig3s_B"," 3S signals",-50,260000);

//C:
  RooRealVar    sigma1s_1_C("sigma1s_1_C","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  RooFormulaVar sigma2s_1_C("sigma2s_1_C","@0*@1",RooArgList(sigma1s_1_C,mRatio21) );
  RooFormulaVar sigma3s_1_C("sigma3s_1_C","@0*@1",RooArgList(sigma1s_1_C,mRatio31) );

  RooRealVar *x1s_C = new RooRealVar("x1s_C","sigma ratio ", x1s_init, 0, 1);

  RooFormulaVar sigma1s_2_C("sigma1s_2_C","@0*@1",RooArgList(sigma1s_1_C, *x1s_C) );
  RooFormulaVar sigma2s_2_C("sigma2s_2_C","@0*@1",RooArgList(sigma1s_2_C,mRatio21) );
  RooFormulaVar sigma3s_2_C("sigma3s_2_C","@0*@1",RooArgList(sigma1s_2_C,mRatio31) );

  RooRealVar alpha1s_1_C("alpha1s_1_C","tail shift", alpha1s_1_init, 1.0, 3.321);
  RooFormulaVar alpha2s_1_C("alpha2s_1_C","1.0*@0",RooArgList(alpha1s_1_C) );
  RooFormulaVar alpha3s_1_C("alpha3s_1_C","1.0*@0",RooArgList(alpha1s_1_C) );
  RooFormulaVar alpha1s_2_C("alpha1s_2_C","1.0*@0",RooArgList(alpha1s_1_C) );
  RooFormulaVar alpha2s_2_C("alpha2s_2_C","1.0*@0",RooArgList(alpha1s_1_C) );
  RooFormulaVar alpha3s_2_C("alpha3s_2_C","1.0*@0",RooArgList(alpha1s_1_C) );

  RooRealVar n1s_1_C("n1s_1_C","power order", n1s_1_init , 1.416, 3.357);
  RooFormulaVar n2s_1_C("n2s_1_C","1.0*@0",RooArgList(n1s_1_C) );
  RooFormulaVar n3s_1_C("n3s_1_C","1.0*@0",RooArgList(n1s_1_C) );
  RooFormulaVar n1s_2_C("n1s_2_C","1.0*@0",RooArgList(n1s_1_C) );
  RooFormulaVar n2s_2_C("n2s_2_C","1.0*@0",RooArgList(n1s_1_C) );
  RooFormulaVar n3s_2_C("n3s_2_C","1.0*@0",RooArgList(n1s_1_C) );

  RooRealVar *f1s_C = new RooRealVar("f1s_C","1S CB fraction", f1s_init, 0, 1);
  RooFormulaVar f2s_C("f2s_C","1.0*@0",RooArgList(*f1s_C) );
  RooFormulaVar f3s_C("f3s_C","1.0*@0",RooArgList(*f1s_C) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1_C = new RooCBShape("cball1s_1_C", "cystal Ball", *(ws->var("mass")), mean1s_C, sigma1s_1_C, alpha1s_1_C, n1s_1_C);
  RooCBShape* cb2s_1_C = new RooCBShape("cball2s_1_C", "cystal Ball", *(ws->var("mass")), mean2s_C, sigma2s_1_C, alpha2s_1_C, n2s_1_C);
  RooCBShape* cb3s_1_C = new RooCBShape("cball3s_1_C", "cystal Ball", *(ws->var("mass")), mean3s_C, sigma3s_1_C, alpha3s_1_C, n3s_1_C);

  RooAddPdf* cb1s_C;
  RooAddPdf* cb2s_C;
  RooAddPdf* cb3s_C;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2_C = new RooCBShape("cball1s_2_C", "cystal Ball", *(ws->var("mass")), mean1s_C, sigma1s_2_C, alpha1s_2_C, n1s_2_C);
  RooCBShape* cb2s_2_C = new RooCBShape("cball2s_2_C", "cystal Ball", *(ws->var("mass")), mean2s_C, sigma2s_2_C, alpha2s_2_C, n2s_2_C);
  RooCBShape* cb3s_2_C = new RooCBShape("cball3s_2_C", "cystal Ball", *(ws->var("mass")), mean3s_C, sigma3s_2_C, alpha3s_2_C, n3s_2_C);
  cb1s_C = new RooAddPdf("cb1s_C","Signal 1S",RooArgList(*cb1s_1_C,*cb1s_2_C), RooArgList(*f1s_C) );
  cb2s_C = new RooAddPdf("cb2s_C","Signal 2S",RooArgList(*cb2s_1_C,*cb2s_2_C), RooArgList(*f1s_C) );
  cb3s_C = new RooAddPdf("cb3s_C","Signal 3S",RooArgList(*cb3s_1_C,*cb3s_2_C), RooArgList(*f1s_C) );

  RooRealVar *nSig1s_C= new RooRealVar("nSig1s_C"," 1S signals",0,1000000);
  RooRealVar *nSig2s_C= new RooRealVar("nSig2s_C"," 2S signals",-20,360000);
  RooRealVar *nSig3s_C= new RooRealVar("nSig3s_C"," 3S signals",-50,260000);

//D:
  RooRealVar    sigma1s_1_D("sigma1s_1_D","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  RooFormulaVar sigma2s_1_D("sigma2s_1_D","@0*@1",RooArgList(sigma1s_1_D,mRatio21) );
  RooFormulaVar sigma3s_1_D("sigma3s_1_D","@0*@1",RooArgList(sigma1s_1_D,mRatio31) );

  RooRealVar *x1s_D = new RooRealVar("x1s_D","sigma ratio ", x1s_init, 0, 1);

  RooFormulaVar sigma1s_2_D("sigma1s_2_D","@0*@1",RooArgList(sigma1s_1_D, *x1s_D) );
  RooFormulaVar sigma2s_2_D("sigma2s_2_D","@0*@1",RooArgList(sigma1s_2_D,mRatio21) );
  RooFormulaVar sigma3s_2_D("sigma3s_2_D","@0*@1",RooArgList(sigma1s_2_D,mRatio31) );

  RooRealVar alpha1s_1_D("alpha1s_1_D","tail shift", alpha1s_1_init, 1.0, 3.321);
  RooFormulaVar alpha2s_1_D("alpha2s_1_D","1.0*@0",RooArgList(alpha1s_1_D) );
  RooFormulaVar alpha3s_1_D("alpha3s_1_D","1.0*@0",RooArgList(alpha1s_1_D) );
  RooFormulaVar alpha1s_2_D("alpha1s_2_D","1.0*@0",RooArgList(alpha1s_1_D) );
  RooFormulaVar alpha2s_2_D("alpha2s_2_D","1.0*@0",RooArgList(alpha1s_1_D) );
  RooFormulaVar alpha3s_2_D("alpha3s_2_D","1.0*@0",RooArgList(alpha1s_1_D) );

  RooRealVar n1s_1_D("n1s_1_D","power order", n1s_1_init , 1.416, 3.357);
  RooFormulaVar n2s_1_D("n2s_1_D","1.0*@0",RooArgList(n1s_1_D) );
  RooFormulaVar n3s_1_D("n3s_1_D","1.0*@0",RooArgList(n1s_1_D) );
  RooFormulaVar n1s_2_D("n1s_2_D","1.0*@0",RooArgList(n1s_1_D) );
  RooFormulaVar n2s_2_D("n2s_2_D","1.0*@0",RooArgList(n1s_1_D) );
  RooFormulaVar n3s_2_D("n3s_2_D","1.0*@0",RooArgList(n1s_1_D) );

  RooRealVar *f1s_D = new RooRealVar("f1s_D","1S CB fraction", f1s_init, 0, 1);
  RooFormulaVar f2s_D("f2s_D","1.0*@0",RooArgList(*f1s_D) );
  RooFormulaVar f3s_D("f3s_D","1.0*@0",RooArgList(*f1s_D) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1_D = new RooCBShape("cball1s_1_D", "cystal Ball", *(ws->var("mass")), mean1s_D, sigma1s_1_D, alpha1s_1_D, n1s_1_D);
  RooCBShape* cb2s_1_D = new RooCBShape("cball2s_1_D", "cystal Ball", *(ws->var("mass")), mean2s_D, sigma2s_1_D, alpha2s_1_D, n2s_1_D);
  RooCBShape* cb3s_1_D = new RooCBShape("cball3s_1_D", "cystal Ball", *(ws->var("mass")), mean3s_D, sigma3s_1_D, alpha3s_1_D, n3s_1_D);

  RooAddPdf* cb1s_D;
  RooAddPdf* cb2s_D;
  RooAddPdf* cb3s_D;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2_D = new RooCBShape("cball1s_2_D", "cystal Ball", *(ws->var("mass")), mean1s_D, sigma1s_2_D, alpha1s_2_D, n1s_2_D);
  RooCBShape* cb2s_2_D = new RooCBShape("cball2s_2_D", "cystal Ball", *(ws->var("mass")), mean2s_D, sigma2s_2_D, alpha2s_2_D, n2s_2_D);
  RooCBShape* cb3s_2_D = new RooCBShape("cball3s_2_D", "cystal Ball", *(ws->var("mass")), mean3s_D, sigma3s_2_D, alpha3s_2_D, n3s_2_D);
  cb1s_D = new RooAddPdf("cb1s_D","Signal 1S",RooArgList(*cb1s_1_D,*cb1s_2_D), RooArgList(*f1s_D) );
  cb2s_D = new RooAddPdf("cb2s_D","Signal 2S",RooArgList(*cb2s_1_D,*cb2s_2_D), RooArgList(*f1s_D) );
  cb3s_D = new RooAddPdf("cb3s_D","Signal 3S",RooArgList(*cb3s_1_D,*cb3s_2_D), RooArgList(*f1s_D) );

  RooRealVar *nSig1s_D= new RooRealVar("nSig1s_D"," 1S signals",0,1000000);
  RooRealVar *nSig2s_D= new RooRealVar("nSig2s_D"," 2S signals",-20,360000);
  RooRealVar *nSig3s_D= new RooRealVar("nSig3s_D"," 3S signals",-50,260000);

//E:
  RooRealVar    sigma1s_1_E("sigma1s_1_E","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  RooFormulaVar sigma2s_1_E("sigma2s_1_E","@0*@1",RooArgList(sigma1s_1_E,mRatio21) );
  RooFormulaVar sigma3s_1_E("sigma3s_1_E","@0*@1",RooArgList(sigma1s_1_E,mRatio31) );

  RooRealVar *x1s_E = new RooRealVar("x1s_E","sigma ratio ", x1s_init, 0, 1);

  RooFormulaVar sigma1s_2_E("sigma1s_2_E","@0*@1",RooArgList(sigma1s_1_E, *x1s_E) );
  RooFormulaVar sigma2s_2_E("sigma2s_2_E","@0*@1",RooArgList(sigma1s_2_E,mRatio21) );
  RooFormulaVar sigma3s_2_E("sigma3s_2_E","@0*@1",RooArgList(sigma1s_2_E,mRatio31) );

  RooRealVar alpha1s_1_E("alpha1s_1_E","tail shift", alpha1s_1_init, 1.0, 3.321);
  RooFormulaVar alpha2s_1_E("alpha2s_1_E","1.0*@0",RooArgList(alpha1s_1_E) );
  RooFormulaVar alpha3s_1_E("alpha3s_1_E","1.0*@0",RooArgList(alpha1s_1_E) );
  RooFormulaVar alpha1s_2_E("alpha1s_2_E","1.0*@0",RooArgList(alpha1s_1_E) );
  RooFormulaVar alpha2s_2_E("alpha2s_2_E","1.0*@0",RooArgList(alpha1s_1_E) );
  RooFormulaVar alpha3s_2_E("alpha3s_2_E","1.0*@0",RooArgList(alpha1s_1_E) );

  RooRealVar n1s_1_E("n1s_1_E","power order", n1s_1_init , 1.416, 3.357);
  RooFormulaVar n2s_1_E("n2s_1_E","1.0*@0",RooArgList(n1s_1_E) );
  RooFormulaVar n3s_1_E("n3s_1_E","1.0*@0",RooArgList(n1s_1_E) );
  RooFormulaVar n1s_2_E("n1s_2_E","1.0*@0",RooArgList(n1s_1_E) );
  RooFormulaVar n2s_2_E("n2s_2_E","1.0*@0",RooArgList(n1s_1_E) );
  RooFormulaVar n3s_2_E("n3s_2_E","1.0*@0",RooArgList(n1s_1_E) );

  RooRealVar *f1s_E = new RooRealVar("f1s_E","1S CB fraction", f1s_init, 0, 1);
  RooFormulaVar f2s_E("f2s_E","1.0*@0",RooArgList(*f1s_E) );
  RooFormulaVar f3s_E("f3s_E","1.0*@0",RooArgList(*f1s_E) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1_E = new RooCBShape("cball1s_1_E", "cystal Ball", *(ws->var("mass")), mean1s_E, sigma1s_1_E, alpha1s_1_E, n1s_1_E);
  RooCBShape* cb2s_1_E = new RooCBShape("cball2s_1_E", "cystal Ball", *(ws->var("mass")), mean2s_E, sigma2s_1_E, alpha2s_1_E, n2s_1_E);
  RooCBShape* cb3s_1_E = new RooCBShape("cball3s_1_E", "cystal Ball", *(ws->var("mass")), mean3s_E, sigma3s_1_E, alpha3s_1_E, n3s_1_E);

  RooAddPdf* cb1s_E;
  RooAddPdf* cb2s_E;
  RooAddPdf* cb3s_E;

  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2_E = new RooCBShape("cball1s_2_E", "cystal Ball", *(ws->var("mass")), mean1s_E, sigma1s_2_E, alpha1s_2_E, n1s_2_E);
  RooCBShape* cb2s_2_E = new RooCBShape("cball2s_2_E", "cystal Ball", *(ws->var("mass")), mean2s_E, sigma2s_2_E, alpha2s_2_E, n2s_2_E);
  RooCBShape* cb3s_2_E = new RooCBShape("cball3s_2_E", "cystal Ball", *(ws->var("mass")), mean3s_E, sigma3s_2_E, alpha3s_2_E, n3s_2_E);
  cb1s_E = new RooAddPdf("cb1s_E","Signal 1S",RooArgList(*cb1s_1_E,*cb1s_2_E), RooArgList(*f1s_E) );
  cb2s_E = new RooAddPdf("cb2s_E","Signal 2S",RooArgList(*cb2s_1_E,*cb2s_2_E), RooArgList(*f1s_E) );
  cb3s_E = new RooAddPdf("cb3s_E","Signal 3S",RooArgList(*cb3s_1_E,*cb3s_2_E), RooArgList(*f1s_E) );

  //RooRealVar *nSig1s_E= new RooRealVar("nSig1s_E"," 1S signals",0,1000000);
  //RooRealVar *nSig2s_E= new RooRealVar("nSig2s_E"," 2S signals",-20,360000);
  //RooRealVar *nSig3s_E= new RooRealVar("nSig3s_E"," 3S signals",-50,260000);

  //BACKGROUND
  double err_sigma_init = 5;
  double err_mu_init = 8;
  double m_lambda_init = 5;
  if (ICset>2) {
    err_mu_init = 5;
    m_lambda_init = 5;
  }

//A:
  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);

  RooGenericPdf *bkg;
  RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  if  (ptLow >= 5)        bkg = bkgHighPt ;
  else bkg = bkgLowPt;

  //RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000); 

//B:
  RooRealVar err_mu_B("#mu_B","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma_B("#sigma_B","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda_B("#lambda_B","m_lambda",  m_lambda_init, 0,25);

  RooGenericPdf *bkg_B;
  RooGenericPdf *bkgLowPt_B = new RooGenericPdf("bkgLowPt_B","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_B, err_mu_B, err_sigma_B) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_B = new RooGenericPdf("bkgHighPt_B","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_B));
  if  (ptLow_B >= 5)        bkg_B = bkgHighPt_B ;
  else bkg_B = bkgLowPt_B;

  RooRealVar *nBkg_B = new RooRealVar("nBkg_B","fraction of component 1 in bkg",10000,0,5000000);  

//C:
  RooRealVar err_mu_C("#mu_C","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma_C("#sigma_C","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda_C("#lambda_C","m_lambda",  m_lambda_init, 0,25);

  RooGenericPdf *bkg_C;
  RooGenericPdf *bkgLowPt_C = new RooGenericPdf("bkgLowPt_C","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_C, err_mu_C, err_sigma_C) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_C = new RooGenericPdf("bkgHighPt_C","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_C));
  if  (ptLow_C >= 5)        bkg_C = bkgHighPt_C ;
  else bkg_C = bkgLowPt_C;

  RooRealVar *nBkg_C = new RooRealVar("nBkg_C","fraction of component 1 in bkg",10000,0,5000000);  

//D:
  RooRealVar err_mu_D("#mu_D","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma_D("#sigma_D","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda_D("#lambda_D","m_lambda",  m_lambda_init, 0,25);

  RooGenericPdf *bkg_D;
  RooGenericPdf *bkgLowPt_D = new RooGenericPdf("bkgLowPt_D","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_D, err_mu_D, err_sigma_D) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_D = new RooGenericPdf("bkgHighPt_D","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_D));
  if  (ptLow_D >= 5)        bkg_D = bkgHighPt_D ;
  else bkg_D = bkgLowPt_D;

  RooRealVar *nBkg_D = new RooRealVar("nBkg_D","fraction of component 1 in bkg",10000,0,5000000);  

//E:
  RooRealVar err_mu_E("#mu_E","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma_E("#sigma_E","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda_E("#lambda_E","m_lambda",  m_lambda_init, 0,25);

  RooGenericPdf *bkg_E;
  RooGenericPdf *bkgLowPt_E = new RooGenericPdf("bkgLowPt_E","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda_E, err_mu_E, err_sigma_E) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt_E = new RooGenericPdf("bkgHighPt_E","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda_E));
  if  (ptLow_E >= 5)        bkg_E = bkgHighPt_E ;
  else bkg_E = bkgLowPt_E;

  //RooRealVar *nBkg_E = new RooRealVar("nBkg_E","fraction of component 1 in bkg",10000,0,5000000);  

//Constraints
  RooFormulaVar *nSig1s = new RooFormulaVar("nSig1s_C"," 1S signals","@0+@1",RooArgList(*nSig1s_B,*nSig1s_C) );
  RooFormulaVar *nSig2s = new RooFormulaVar("nSig2s_C"," 2S signals","@0+@1",RooArgList(*nSig2s_B,*nSig2s_C) );
  RooFormulaVar *nSig3s = new RooFormulaVar("nSig3s_C"," 3S signals","@0+@1",RooArgList(*nSig3s_B,*nSig3s_C) );
  RooFormulaVar *nBkg = new RooFormulaVar("nBkg","fraction of component 1 in bkg","@0+@1",RooArgList(*nBkg_B,*nBkg_C) );

  RooFormulaVar *nSig1s_E = new RooFormulaVar("nSig1s_E"," 1S signals","@0+@1-@2",RooArgList(*nSig1s_B,*nSig1s_C,*nSig1s_D) );
  RooFormulaVar *nSig2s_E = new RooFormulaVar("nSig2s_E"," 2S signals","@0+@1-@2",RooArgList(*nSig2s_B,*nSig2s_C,*nSig2s_D) );
  RooFormulaVar *nSig3s_E = new RooFormulaVar("nSig3s_E"," 3S signals","@0+@1-@2",RooArgList(*nSig3s_B,*nSig3s_C,*nSig3s_D) );
  RooFormulaVar *nBkg_E = new RooFormulaVar("nBkg_E","fraction of component 1 in bkg","@0+@1-@2",RooArgList(*nBkg_B,*nBkg_C,*nBkg_D) );

  //Build the model
  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  ws->import(*model);

  c_A->cd();
  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  RooAddPdf* model_B = new RooAddPdf();
  model_B = new RooAddPdf("model_B","1S+2S+3S + Bkg",RooArgList(*cb1s_B, *cb2s_B, *cb3s_B, *bkg_B),RooArgList(*nSig1s_B,*nSig2s_B,*nSig3s_B,*nBkg_B));

  ws->import(*model_B);

  c_B->cd();
  RooPlot* myPlot2_B = (RooPlot*)myPlot_B->Clone();
  ws->data("reducedDS_B")->plotOn(myPlot2_B,Name("dataOS_FIT_B"),MarkerSize(.8));

  RooAddPdf* model_C = new RooAddPdf();
  model_C = new RooAddPdf("model_C","1S+2S+3S + Bkg",RooArgList(*cb1s_C, *cb2s_C, *cb3s_C, *bkg_C),RooArgList(*nSig1s_C,*nSig2s_C,*nSig3s_C,*nBkg_C));

  ws->import(*model_C);

  c_C->cd();
  RooPlot* myPlot2_C = (RooPlot*)myPlot_C->Clone();
  ws->data("reducedDS_C")->plotOn(myPlot2_C,Name("dataOS_FIT_C"),MarkerSize(.8));

  RooAddPdf* model_D = new RooAddPdf();
  model_D = new RooAddPdf("model_D","1S+2S+3S + Bkg",RooArgList(*cb1s_D, *cb2s_D, *cb3s_D, *bkg_D),RooArgList(*nSig1s_D,*nSig2s_D,*nSig3s_D,*nBkg_D));

  ws->import(*model_D);

  c_D->cd();
  RooPlot* myPlot2_D = (RooPlot*)myPlot_D->Clone();
  ws->data("reducedDS_D")->plotOn(myPlot2_D,Name("dataOS_FIT_D"),MarkerSize(.8));

  RooAddPdf* model_E = new RooAddPdf();
  model_E = new RooAddPdf("model_E","1S+2S+3S + Bkg",RooArgList(*cb1s_E, *cb2s_E, *cb3s_E, *bkg_E),RooArgList(*nSig1s_E,*nSig2s_E,*nSig3s_E,*nBkg_E));

  ws->import(*model_E);

  c_E->cd();
  RooPlot* myPlot2_E = (RooPlot*)myPlot_E->Clone();
  ws->data("reducedDS_E")->plotOn(myPlot2_E,Name("dataOS_FIT_E"),MarkerSize(.8));


  // Construct simultaneous PDF for A and B
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",tp);
  simPdf->addPdf(*(ws->pdf("model")),"A") ;
  simPdf->addPdf(*(ws->pdf("model_B")),"B") ;
  simPdf->addPdf(*(ws->pdf("model_C")),"C") ;
  simPdf->addPdf(*(ws->pdf("model_D")),"D") ;
  simPdf->addPdf(*(ws->pdf("model_E")),"E") ;
  ws->import(*simPdf);

  cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("simPdf")->fitTo(*dsABC,Save(), Hesse(kTRUE),Range(massLow,massHigh),Timer(kTRUE),Extended(kTRUE));
  cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;

  c_A->cd();
  //Fit the model to the data
  //RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow,massHigh),Timer(kTRUE),Extended(kTRUE));
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
  //myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.054);
  //myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  //myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();

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

  c_C->cd();
  ws->pdf("model_C")->plotOn(myPlot2_C,Name("modelHist_C"));
  ws->pdf("model_C")->plotOn(myPlot2_C,Name("Sig1S_C"),Components(RooArgSet(*cb1s_C)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_C")->plotOn(myPlot2_C,Components(RooArgSet(*cb2s_C)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_C")->plotOn(myPlot2_C,Components(RooArgSet(*cb3s_C)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_C")->plotOn(myPlot2_C,Name("bkgPDF_C"),Components(RooArgSet(*bkg_C)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_C->SetFillStyle(4000);
  myPlot2_C->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_C->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_C->GetYaxis()->CenterTitle();
  //myPlot2_C->GetYaxis()->SetTitleSize(0.058);
  myPlot2_C->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_C->GetXaxis()->SetLabelSize(0);
  myPlot2_C->GetXaxis()->SetRangeUser(8,14);
  //myPlot2_C->GetXaxis()->SetTitleSize(0);
  myPlot2_C->Draw();

  c_D->cd();
  ws->pdf("model_D")->plotOn(myPlot2_D,Name("modelHist_D"));
  ws->pdf("model_D")->plotOn(myPlot2_D,Name("Sig1S_D"),Components(RooArgSet(*cb1s_D)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_D")->plotOn(myPlot2_D,Components(RooArgSet(*cb2s_D)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_D")->plotOn(myPlot2_D,Components(RooArgSet(*cb3s_D)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_D")->plotOn(myPlot2_D,Name("bkgPDF_D"),Components(RooArgSet(*bkg_D)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_D->SetFillStyle(4000);
  myPlot2_D->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_D->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_D->GetYaxis()->CenterTitle();
  //myPlot2_D->GetYaxis()->SetTitleSize(0.058);
  myPlot2_D->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_D->GetXaxis()->SetLabelSize(0);
  myPlot2_D->GetXaxis()->SetRangeUser(8,14);
  //myPlot2_D->GetXaxis()->SetTitleSize(0);
  myPlot2_D->Draw();

  c_E->cd();
  ws->pdf("model_E")->plotOn(myPlot2_E,Name("modelHist_E"));
  ws->pdf("model_E")->plotOn(myPlot2_E,Name("Sig1S_E"),Components(RooArgSet(*cb1s_E)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_E")->plotOn(myPlot2_E,Components(RooArgSet(*cb2s_E)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_E")->plotOn(myPlot2_E,Components(RooArgSet(*cb3s_E)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model_E")->plotOn(myPlot2_E,Name("bkgPDF_E"),Components(RooArgSet(*bkg_E)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2_E->SetFillStyle(4000);
  myPlot2_E->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2_E->GetYaxis()->SetTitleOffset(1.43);
  myPlot2_E->GetYaxis()->CenterTitle();
  //myPlot2_E->GetYaxis()->SetTitleSize(0.058);
  myPlot2_E->GetYaxis()->SetLabelSize(0.054);
  //myPlot2_E->GetXaxis()->SetLabelSize(0);
  myPlot2_E->GetXaxis()->SetRangeUser(8,14);
  //myPlot2_E->GetXaxis()->SetTitleSize(0);
  myPlot2_E->Draw();

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString outFileName;
  if (whichModel){
    outFileName = Form("SimultaneousFits/altfitresults_upsilon_%s.root",kineLabel.Data());
  }
  else {
    outFileName = Form("SimultaneousFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
  }
  TFile* outf = new TFile(outFileName,"recreate");

  TNtuple* Range = new TNtuple("Range","Bin","binlowlim:binuplim",1);
  Range->Fill(ptLow,ptHigh);
  TNtuple* Range_B = new TNtuple("Range_B","Bin","binlowlim:binuplim",1);
  Range_B->Fill(ptLow_B,ptHigh_B);
  TNtuple* Range_C = new TNtuple("Range_C","Bin","binlowlim:binuplim",1);
  Range_C->Fill(ptLow_C,ptHigh_C);
  TNtuple* Range_D = new TNtuple("Range_D","Bin","binlowlim:binuplim",1);
  Range_D->Fill(ptLow_D,ptHigh_D);
  TNtuple* Range_E = new TNtuple("Range_E","Bin","binlowlim:binuplim",1);
  Range_E->Fill(ptLow_E,ptHigh_E);

  Range->Write();
  Range_B->Write();
  Range_C->Write();
  Range_D->Write();
  Range_E->Write();
  ws->Write();
  outf->Close();

} 
 
