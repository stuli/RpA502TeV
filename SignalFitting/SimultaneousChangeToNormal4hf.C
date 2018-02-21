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


using namespace std;
using namespace RooFit;
void SimultaneousChangeToNormal4hf( 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-1.93, float yHigh=0.0,//Run 1 has p going in -z direction
       TString ext = "_A",
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  float hfLowGet, hfHighGet;

  bool depBin = kFALSE;
  if (ext=="_D") {
    depBin = kTRUE;
    hfLowGet = 27;
    hfHighGet = 120;
  }
  else if (ext=="_A") {
    hfLowGet = 0;
    hfHighGet = 12;
  }
  else if (ext=="_B") {
    hfLowGet = 12;
    hfHighGet = 19;
  }
  else if (ext=="_C") {
    hfLowGet = 19;
    hfHighGet = 27;
  }

  cout << endl << "EXTRACTING FIT IN HF BIN [" << hfLowGet << "," << hfHighGet << "]" << endl << endl;

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
  //Select Data Set
  if (collId==kPADATA) {
    f1 = new TFile("../yskimPA1st_OpSign_20177262037_unIdentified.root");
    f2 = new TFile("../yskimPA2nd_OpSign_20177262044_unIdentified.root");
    yLowLab = yLow+0.47;
    yHighLab = yHigh+0.47;
  }

  //Define cuts
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",0, 120, 0, 400 );
  TString kineLabelGet = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  kineLabelGet = kineLabelGet + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLowGet, hfHighGet, 0, 400 );
  TString kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow, ptHigh, yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
  TString hfntracksCut = Form("&& (hfpluseta4+hfminuseta4>%.2f) && (hfpluseta4+hfminuseta4<%.2f)",hfLowGet, hfHighGet );

  if (muPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);

  kineCut = kineCut + hfntracksCut;

  //import and merge datasets
  RooDataSet *PAdataset = (RooDataSet*)f1->Get("dataset");
  if (collId==kPADATA) {
    RooDataSet *PAdataset2 = (RooDataSet*)f2->Get("dataset");
    PAdataset->append(*PAdataset2);
  }
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*PAdataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)PAdataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  //import generating models
  cout << "Importing workspace" << endl;
  TString NomFileName = Form("SimultaneousFits/sim_hfSplit4_nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");

  float ups1smass = Nomws->var("m_{#Upsilon(1S)}"+ext)->getVal();
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",ups1smass, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  mean1s.setError(Nomws->var("m_{#Upsilon(1S)}"+ext)->getError());
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );

  //SIGNAL:
  double sigma1s_1_init = Nomws->var("sigma1s_1"+ext)->getVal();
  double x1s_init = Nomws->var("x1s"+ext)->getVal();
  double alpha1s_1_init = Nomws->var("alpha1s_1"+ext)->getVal();
  double n1s_1_init = Nomws->var("n1s_1"+ext)->getVal();
  double f1s_init = Nomws->var("f1s"+ext)->getVal();

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",sigma1s_1_init, 0.02, 0.3);
  sigma1s_1.setError(Nomws->var("sigma1s_1"+ext)->getError());
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", x1s_init, 0, 2.4);
  x1s->setError(Nomws->var("x1s"+ext)->getError());

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", alpha1s_1_init, 1.0, 3.321);
  alpha1s_1.setError(Nomws->var("alpha1s_1"+ext)->getError());
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", n1s_1_init , 1.416, 3.357);
  n1s_1.setError(Nomws->var("n1s_1"+ext)->getError());
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", f1s_init, 0, 1);
  f1s->setError(Nomws->var("f1s"+ext)->getError());
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(*f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(*f1s) );

  // Set up crystal ball shapes
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_1, alpha3s_1, n3s_1);

  RooAddPdf* cb1s;
  RooAddPdf* cb2s;
  RooAddPdf* cb3s;

if (whichModel) {
  //CB+GAUSSIAN
  RooGaussian* gauss1s = new RooGaussian("gauss1s","gaussian PDF",*(ws->var("mass")),mean1s,sigma1s_2);
  RooGaussian* gauss2s = new RooGaussian("gauss2s","gaussian PDF",*(ws->var("mass")),mean2s,sigma2s_2);
  RooGaussian* gauss3s = new RooGaussian("gauss3s","gaussian PDF",*(ws->var("mass")),mean3s,sigma3s_2);
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*gauss1s), RooArgList(*f1s) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*gauss2s), RooArgList(*f1s) );
  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*gauss3s), RooArgList(*f1s) );
}
else {
  //DOUBLE CRYSTAL BALL
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_2, alpha3s_2, n3s_2);
  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1s) );
  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*cb3s_2), RooArgList(*f1s) );
}

    double nSig1s_init;
    double nSig2s_init;
    double nSig3s_init;
  if (depBin) {
    nSig1s_init = Nomws->var("nSig1s_E")->getVal() - Nomws->var("nSig1s_A")->getVal() - Nomws->var("nSig1s_B")->getVal() - Nomws->var("nSig1s_C")->getVal();
    nSig2s_init = Nomws->var("nSig2s_E")->getVal() - Nomws->var("nSig2s_A")->getVal() - Nomws->var("nSig2s_B")->getVal() - Nomws->var("nSig2s_C")->getVal();
    nSig3s_init = Nomws->var("nSig3s_E")->getVal() - Nomws->var("nSig3s_A")->getVal() - Nomws->var("nSig3s_B")->getVal() - Nomws->var("nSig3s_C")->getVal();
  }
  else {
    nSig1s_init = Nomws->var("nSig1s"+ext)->getVal();
    nSig2s_init = Nomws->var("nSig2s"+ext)->getVal();
    nSig3s_init = Nomws->var("nSig3s"+ext)->getVal();
  }
  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",nSig1s_init,0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",nSig2s_init,-20,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",nSig3s_init,-50,260000);

  //Get bogus errors
  if (depBin) {
    double err1A = Nomws->var("nSig1s_A")->getError();
    double err1B = Nomws->var("nSig1s_B")->getError();
    double err2A = Nomws->var("nSig2s_A")->getError();
    double err2B = Nomws->var("nSig2s_B")->getError();
    double err3A = Nomws->var("nSig3s_A")->getError();
    double err3B = Nomws->var("nSig3s_B")->getError();
    nSig1s->setError(0.5*TMath::Sqrt(err1A*err1A + err1B*err1B));
    nSig2s->setError(0.5*TMath::Sqrt(err2A*err2A + err2B*err2B));
    nSig3s->setError(0.5*TMath::Sqrt(err3A*err3A + err3B*err3B));    
  }
  else {
    nSig1s->setError(Nomws->var("nSig1s"+ext)->getError());
    nSig2s->setError(Nomws->var("nSig2s"+ext)->getError());
    nSig3s->setError(Nomws->var("nSig3s"+ext)->getError());
  }

  //BACKGROUND
  double err_mu_init = 8;
  double err_sigma_init = 5;
  if (ptLow < 5) {
    err_mu_init = Nomws->var("#mu"+ext)->getVal();
    err_sigma_init = Nomws->var("#sigma"+ext)->getVal();
  }
  double m_lambda_init = Nomws->var("#lambda"+ext)->getVal();
  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);
  if (ptLow < 5) {
    err_mu.setError(Nomws->var("#mu"+ext)->getError());
    err_sigma.setError(Nomws->var("#sigma"+ext)->getError());
  }
  m_lambda.setError(Nomws->var("#lambda"+ext)->getError());

  RooGenericPdf *bkg;
  RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  
  if  (ptLow >= 5)        bkg = bkgHighPt ;
  else bkg = bkgLowPt;

  double nBkg_init = Nomws->var("nBkg"+ext)->getVal();
  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",nBkg_init,0,5000000);  
  nBkg->setError(Nomws->var("nBkg"+ext)->getError());

  //fix all parameters
  mean1s.setConstant(kTRUE);

  alpha1s_1.setConstant(kTRUE);
  sigma1s_1.setConstant(kTRUE);
  f1s->setConstant(kTRUE);
  x1s->setConstant(kTRUE);
  n1s_1.setConstant(kTRUE);

  nSig1s.setConstant(kTRUE);
  nSig2s.setConstant(kTRUE);
  nSig3s.setConstant(kTRUE);

  err_mu.setConstant(kTRUE);
  err_sigma.setConstant(kTRUE);
  m_lambda.setConstant(kTRUE);
  nBkg.setConstant(kTRUE);

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
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pPb
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

  if (whichModel) {
    c1->SaveAs(Form("SimultaneousFits/Normal/altfitresults_upsilon_%s.png",kineLabelGet.Data()));
    c1->SaveAs(Form("SimultaneousFits/Normal/altfitresults_upsilon_%s.pdf",kineLabelGet.Data()));
  }
  else {
    c1->SaveAs(Form("SimultaneousFits/Normal/nomfitresults_upsilon_%s.png",kineLabelGet.Data()));
    c1->SaveAs(Form("SimultaneousFits/Normal/nomfitresults_upsilon_%s.pdf",kineLabelGet.Data()));
  }
  
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

TString outFileName;
  if (whichModel){
    outFileName = Form("SimultaneousFits/Normal/altfitresults_upsilon_%s.root",kineLabelGet.Data());
  }
  else {
    outFileName = Form("SimultaneousFits/Normal/nomfitresults_upsilon_%s.root",kineLabelGet.Data());
  }
  TFile* outf = new TFile(outFileName,"recreate");
  outh->Write();
  c1->Write();
  ws->Write();
  outf->Close();

} 
 
