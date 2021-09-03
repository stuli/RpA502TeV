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
#include "../../HeaderFiles/CMS_lumi_Jared.C"
#include "../../HeaderFiles/tdrstyle.C"
#include "../../HeaderFiles/StyleSetting.h"


using namespace std;
using namespace RooFit;
void SuperRpARAAPlot_SameBkg_withoutPull_ppbase( 
       int collId = kPPDATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  TGaxis::SetMaxDigits(3);

  int A = 208;
  float PPeff = 0.824;//Efficiency values as of 4/2/2018
  float PAeff = 0.795;//EffIntRpA = 1.036
  float AAeff = 0.729;//From 16-023 integrated bin
  float PPlum = 28000;//28 pb^-1
  float PAlum = 34.6;//34.6 nb^-1
  float AAlum = 0.368;//368 microb^-1
  float RAA1 = 0.3518;
  float RpA1 = 0.786;

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

//*****************************************************************************
//*****************************************************************************
//             MAKE THE PP MODEL
//*****************************************************************************
//*****************************************************************************

  float yLowLab;
  float yHighLab;

  collId = kPPDATA;

  //Select Data Set
  yLow = 0.0;
  yLowLab = yLow;
  yHighLab = yHigh;

  //import the model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, 0, 100, 0.0, 0.5);
  TString NomFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2018_08_02/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");

  RooAbsData* reducedDS = Nomws->data("reducedDS");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*reducedDS);
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->setVal(Nomws->var("mass")->getVal());
  ws->var("mass")->Print();

  //Plot it
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,550);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.19, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  double mean1s_init = Nomws->var("m_{#Upsilon(1S)}")->getVal();
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",mean1s_init, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );

  //SIGNAL:
  double sigma1s_1_init = Nomws->var("sigma1s_1")->getVal();
  double x1s_init = Nomws->var("x1s")->getVal();
  double alpha1s_1_init = Nomws->var("alpha1s_1")->getVal();
  double n1s_1_init = Nomws->var("n1s_1")->getVal();
  double f1s_init = Nomws->var("f1s")->getVal();

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

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",-20,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",-50,260000);
  double nSig1s_init = Nomws->var("nSig1s")->getVal();
  double nSig2s_init = Nomws->var("nSig2s")->getVal();
  double nSig3s_init = Nomws->var("nSig3s")->getVal();
  nSig1s->setVal(nSig1s_init);
  nSig2s->setVal(nSig2s_init);
  nSig3s->setVal(nSig3s_init);

  //BACKGROUND
  double err_mu_init = Nomws->var("#mu")->getVal();
  double err_sigma_init = Nomws->var("#sigma")->getVal();
  double m_lambda_init = Nomws->var("#lambda")->getVal();

  RooRealVar err_mu("#mu","err_mu", err_mu_init,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, 0,25);
  
  RooGenericPdf *bkg;
  RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
  
  if  (ptLow >= 5)        bkg = bkgHighPt ;
  else bkg = bkgLowPt;

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000);
  nBkg_init = Nomws->var("nBkg")->getVal();
  nBkg->setVal(nBkg_init);

  //Build the model
  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));

  ws->import(*model);
  //Plot it
  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"),NormRange("full"));
  //ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2),NormRange("full"));
  //ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2),NormRange("full"));
  //ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2),NormRange("full"));
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2),NormRange("full"));

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.0);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.045);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  TString perc = "%";

  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.075;
  float text_size = 16;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if (collId!=kPADATA) {
    if(yLow==0) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < |y_{CM}^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  else if (collId==kPADATA) {
    if(yLow==-yHigh) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < y_{CM}^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < 2.4"), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  TString perc = "%";
  drawText(Form("Centrality 0-100%s", perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);

  TLegend* fitleg = new TLegend(0.66,0.35,0.81,0.65); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");

  // PULL
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.25);
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
  pullFrame->GetYaxis()->SetTitleOffset(0.3) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.18) ; //19
  pullFrame->GetYaxis()->SetLabelSize(0.113) ; // 113
  pullFrame->GetYaxis()->SetRangeUser(-6.5,6.5) ;
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

  //continue beautifying the plot and print out results
  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  pad1->Update();

  setTDRStyle();
  //writeExtraText = true;
  //extraText = "Preliminary";

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




//*****************************************************************************
//*****************************************************************************
//             MAKE THE AA MODEL
//*****************************************************************************
//*****************************************************************************

  collId = kAADATA,

  //Select Data Set
  yLow = 0.0;
  yLowLab = yLow;
  yHighLab = yHigh;

  //import the model
  cout << "Importing workspace" << endl;
  TString AAkineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, 0, 100, 0.0, 0.5);
  TString AANomFileName = Form("nomfitresults_upsilon_%s.root",AAkineLabel.Data());
  cout << AANomFileName << endl;
  TFile* AANomFile = TFile::Open(AANomFileName,"READ");
  RooWorkspace *AANomws = (RooWorkspace*)AANomFile->Get("workspace");
  AANomFile->Close("R");

  RooAbsData* AAreducedDS = AANomws->data("reducedDS");
  RooWorkspace *AAws = new RooWorkspace("workspace");
  AAws->import(*AAreducedDS);
  AAws->var("mass")->setRange(massLow, massHigh);
  AAws->var("mass")->setVal(AANomws->var("mass")->getVal());
  AAws->var("mass")->Print();

  //Plot it
  TCanvas* AAc1 =  new TCanvas("AAcanvas2","My plots",504,45,550,550);
  AAc1->cd();
  TPad *AApad1 = new TPad("AApad1", "AApad1", 0, 0.16, 0.98, 1.0);
  AApad1->SetTicks(1,1);
  AApad1->Draw(); AApad1->cd();

  RooPlot* AAmyPlot = AAws->var("mass")->frame(nMassBin); // bins
  AAws->data("reducedDS")->plotOn(AAmyPlot,Name("AAdataHist"));
  double AAmean1s_init = AANomws->var("m_{#Upsilon(1S)}")->getVal();
  RooRealVar AAmean1s("AAm_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",AAmean1s_init, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar AAmRatio21("AAmRatio21","AAmRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar AAmRatio31("AAmRatio31","AAmRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar AAmean2s("AAmean2s","AAm_{#Upsilon(1S)}*AAmRatio21", RooArgSet(AAmean1s,AAmRatio21) );
  RooFormulaVar AAmean3s("AAmean3s","AAm_{#Upsilon(1S)}*AAmRatio31", RooArgSet(AAmean1s,AAmRatio31) );

  //SIGNAL:
  double AAsigma1s_1_init = AANomws->var("sigma1s_1")->getVal();
  double AAx1s_init = AANomws->var("x1s")->getVal();
  double AAalpha1s_1_init = AANomws->var("alpha1s_1")->getVal();
  double AAn1s_1_init = AANomws->var("n1s_1")->getVal();
  double AAf1s_init = AANomws->var("f1s")->getVal();

  RooRealVar    AAsigma1s_1("AAsigma1s_1","width/sigma of the signal gaussian mass PDF",AAsigma1s_1_init, 0.02, 0.3);
  RooFormulaVar AAsigma2s_1("AAsigma2s_1","@0*@1",RooArgList(AAsigma1s_1,AAmRatio21) );
  RooFormulaVar AAsigma3s_1("AAsigma3s_1","@0*@1",RooArgList(AAsigma1s_1,AAmRatio31) );

  RooRealVar *AAx1s = new RooRealVar("AAx1s","sigma ratio ", AAx1s_init, 0, 1);

  RooFormulaVar AAsigma1s_2("AAsigma1s_2","@0*@1",RooArgList(AAsigma1s_1, *AAx1s) );
  RooFormulaVar AAsigma2s_2("AAsigma2s_2","@0*@1",RooArgList(AAsigma1s_2,AAmRatio21) );
  RooFormulaVar AAsigma3s_2("AAsigma3s_2","@0*@1",RooArgList(AAsigma1s_2,AAmRatio31) );

  RooRealVar AAalpha1s_1("AAalpha1s_1","tail shift", AAalpha1s_1_init, 1.0, 3.321);
  RooFormulaVar AAalpha2s_1("AAalpha2s_1","1.0*@0",RooArgList(AAalpha1s_1) );
  RooFormulaVar AAalpha3s_1("AAalpha3s_1","1.0*@0",RooArgList(AAalpha1s_1) );
  RooFormulaVar AAalpha1s_2("AAalpha1s_2","1.0*@0",RooArgList(AAalpha1s_1) );
  RooFormulaVar AAalpha2s_2("AAalpha2s_2","1.0*@0",RooArgList(AAalpha1s_1) );
  RooFormulaVar AAalpha3s_2("AAalpha3s_2","1.0*@0",RooArgList(AAalpha1s_1) );

  RooRealVar AAn1s_1("AAn1s_1","power order", AAn1s_1_init , 1.416, 3.357);
  RooFormulaVar AAn2s_1("AAn2s_1","1.0*@0",RooArgList(AAn1s_1) );
  RooFormulaVar AAn3s_1("AAn3s_1","1.0*@0",RooArgList(AAn1s_1) );
  RooFormulaVar AAn1s_2("AAn1s_2","1.0*@0",RooArgList(AAn1s_1) );
  RooFormulaVar AAn2s_2("AAn2s_2","1.0*@0",RooArgList(AAn1s_1) );
  RooFormulaVar AAn3s_2("AAn3s_2","1.0*@0",RooArgList(AAn1s_1) );
  
  RooRealVar *AAf1s = new RooRealVar("AAf1s","1S CB fraction", AAf1s_init, 0, 1);
  RooFormulaVar AAf2s("AAf2s","1.0*@0",RooArgList(*AAf1s) );
  RooFormulaVar AAf3s("AAf3s","1.0*@0",RooArgList(*AAf1s) );

  // Set up crystal ball shapes
  RooCBShape* AAcb1s_1 = new RooCBShape("AAcball1s_1", "cystal Ball", *(AAws->var("mass")), AAmean1s, AAsigma1s_1, AAalpha1s_1, AAn1s_1);
  RooCBShape* AAcb2s_1 = new RooCBShape("AAcball2s_1", "cystal Ball", *(AAws->var("mass")), AAmean2s, AAsigma2s_1, AAalpha2s_1, AAn2s_1);
  RooCBShape* AAcb3s_1 = new RooCBShape("AAcball3s_1", "cystal Ball", *(AAws->var("mass")), AAmean3s, AAsigma3s_1, AAalpha3s_1, AAn3s_1);

  RooAddPdf* AAcb1s;
  RooAddPdf* AAcb2s;
  RooAddPdf* AAcb3s;

if (whichModel) {
  //CB+GAUSSIAN
  RooGaussian* AAgauss1s = new RooGaussian("AAgauss1s","gaussian PDF",*(AAws->var("mass")),AAmean1s,AAsigma1s_2);
  RooGaussian* AAgauss2s = new RooGaussian("AAgauss2s","gaussian PDF",*(AAws->var("mass")),AAmean2s,AAsigma2s_2);
  RooGaussian* AAgauss3s = new RooGaussian("AAgauss3s","gaussian PDF",*(AAws->var("mass")),AAmean3s,AAsigma3s_2);
  AAcb1s = new RooAddPdf("AAcb1s","Signal 1S",RooArgList(*AAcb1s_1,*AAgauss1s), RooArgList(*AAf1s) );
  AAcb2s = new RooAddPdf("AAcb2s","Signal 2S",RooArgList(*AAcb2s_1,*AAgauss2s), RooArgList(*AAf1s) );
  AAcb3s = new RooAddPdf("AAcb3s","Signal 3S",RooArgList(*AAcb3s_1,*AAgauss3s), RooArgList(*AAf1s) );
}
else {
  //DOUBLE CRYSTAL BALL
  RooCBShape* AAcb1s_2 = new RooCBShape("AAcball1s_2", "cystal Ball", *(AAws->var("mass")), AAmean1s, AAsigma1s_2, AAalpha1s_2, AAn1s_2);
  RooCBShape* AAcb2s_2 = new RooCBShape("AAcball2s_2", "cystal Ball", *(AAws->var("mass")), AAmean2s, AAsigma2s_2, AAalpha2s_2, AAn2s_2);
  RooCBShape* AAcb3s_2 = new RooCBShape("AAcball3s_2", "cystal Ball", *(AAws->var("mass")), AAmean3s, AAsigma3s_2, AAalpha3s_2, AAn3s_2);
  AAcb1s = new RooAddPdf("AAcb1s","Signal 1S",RooArgList(*AAcb1s_1,*AAcb1s_2), RooArgList(*AAf1s) );
  AAcb2s = new RooAddPdf("AAcb2s","Signal 2S",RooArgList(*AAcb2s_1,*AAcb2s_2), RooArgList(*AAf1s) );
  AAcb3s = new RooAddPdf("AAcb3s","Signal 3S",RooArgList(*AAcb3s_1,*AAcb3s_2), RooArgList(*AAf1s) );
}

  //float scalefactor = 1;
  float scalefactor = Nomws->var("nSig1s")->getVal()*RAA1/(AANomws->var("nSig1s")->getVal());

  RooRealVar *AAnSig1s= new RooRealVar("AAnSig1s"," 1S signals",0,1000000);
  RooRealVar *AAnSig2s= new RooRealVar("AAnSig2s"," 2S signals",-20,360000);
  RooRealVar *AAnSig3s= new RooRealVar("AAnSig3s"," 3S signals",-50,260000);
  double AAnSig1s_init = AANomws->var("nSig1s")->getVal()*scalefactor;
  double AAnSig2s_init = AANomws->var("nSig2s")->getVal()*scalefactor;
  double AAnSig3s_init = AANomws->var("nSig3s")->getVal()*scalefactor;
  AAnSig1s->setVal(AAnSig1s_init);
  AAnSig2s->setVal(AAnSig2s_init);
  AAnSig3s->setVal(AAnSig3s_init);

  //BACKGROUND
  RooRealVar AAerr_mu("#AAmu","AAerr_mu", err_mu_init,  0, 25) ;
  RooRealVar AAerr_sigma("#AAsigma","AAerr_sigma", err_sigma_init, 0,25);
  RooRealVar AAm_lambda("#AAlambda","AAm_lambda",  m_lambda_init, 0,25);
  
  RooGenericPdf *AAbkg;
  RooGenericPdf *AAbkgLowPt = new RooGenericPdf("AAbkgLowPt","AABackground","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(AAws->var("mass")), AAm_lambda, AAerr_mu, AAerr_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *AAbkgHighPt = new RooGenericPdf("AAbkgHighPt","AABackground","TMath::Exp(-@0/@1)",RooArgList(*(AAws->var("mass")),AAm_lambda));
  
  if  (ptLow >= 5)        AAbkg = AAbkgHighPt ;
  else AAbkg = AAbkgLowPt;

  RooRealVar *AAnBkg = new RooRealVar("AAnBkg","fraction of component 1 in bkg",10000,0,5000000);
  AAnBkg->setVal(nBkg_init);

  //Build the model
  RooAddPdf* AAmodel = new RooAddPdf();
  AAmodel = new RooAddPdf("AAmodel","1S+2S+3S + Bkg",RooArgList(*AAcb1s, *AAcb2s, *AAcb3s, *AAbkg),RooArgList(*AAnSig1s,*AAnSig2s,*AAnSig3s,*AAnBkg));

  AAws->import(*AAmodel);
  //Plot it
  RooPlot* AAmyPlot2 = (RooPlot*)AAmyPlot->Clone();
  AAws->data("reducedDS")->plotOn(AAmyPlot2,Name("AAdataOS_FIT"),MarkerSize(.8));

  AAws->pdf("AAmodel")->plotOn(AAmyPlot2,Name("AAmodelHist"),NormRange("full"));



//*****************************************************************************
//*****************************************************************************
//             MAKE THE PA MODEL
//*****************************************************************************
//*****************************************************************************

  collId = kPADATA;

  //Select Data Set
  yLow = -1.93;
  yLowLab = yLow+0.47;
  yHighLab = yHigh+0.47;

  //import the model
  cout << "Importing workspace" << endl;
  TString PAkineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString PANomFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2018_07_18/nomfitresults_upsilon_%s.root",PAkineLabel.Data());
  cout << PANomFileName << endl;
  TFile* PANomFile = TFile::Open(PANomFileName,"READ");
  RooWorkspace *PANomws = (RooWorkspace*)PANomFile->Get("workspace");
  PANomFile->Close("R");

  RooAbsData* PAreducedDS = PANomws->data("reducedDS");
  RooWorkspace *PAws = new RooWorkspace("workspace");
  PAws->import(*PAreducedDS);
  PAws->var("mass")->setRange(massLow, massHigh);
  PAws->var("mass")->setVal(PANomws->var("mass")->getVal());
  PAws->var("mass")->Print();

  //Plot it
  TCanvas* PAc1 =  new TCanvas("PAcanvas2","My plots",40,45,550,550);
  PAc1->cd();
  TPad *PApad1 = new TPad("PApad1", "PApad1", 0, 0.19, 0.98, 1.0);
  PApad1->SetTicks(1,1);
  PApad1->Draw(); PApad1->cd();

  RooPlot* PAmyPlot = PAws->var("mass")->frame(nMassBin); // bins
  PAws->data("reducedDS")->plotOn(PAmyPlot,Name("PAdataHist"));
  double PAmean1s_init = PANomws->var("m_{#Upsilon(1S)}")->getVal();
  RooRealVar PAmean1s("PAm_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",PAmean1s_init, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar PAmRatio21("PAmRatio21","PAmRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar PAmRatio31("PAmRatio31","PAmRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar PAmean2s("PAmean2s","PAm_{#Upsilon(1S)}*PAmRatio21", RooArgSet(PAmean1s,PAmRatio21) );
  RooFormulaVar PAmean3s("PAmean3s","PAm_{#Upsilon(1S)}*PAmRatio31", RooArgSet(PAmean1s,PAmRatio31) );

  //SIGNAL:
  double PAsigma1s_1_init = PANomws->var("sigma1s_1")->getVal();
  double PAx1s_init = PANomws->var("x1s")->getVal();
  double PAalpha1s_1_init = PANomws->var("alpha1s_1")->getVal();
  double PAn1s_1_init = PANomws->var("n1s_1")->getVal();
  double PAf1s_init = PANomws->var("f1s")->getVal();

  RooRealVar    PAsigma1s_1("PAsigma1s_1","width/sigma of the signal gaussian mass PDF",PAsigma1s_1_init, 0.02, 0.3);
  RooFormulaVar PAsigma2s_1("PAsigma2s_1","@0*@1",RooArgList(PAsigma1s_1,PAmRatio21) );
  RooFormulaVar PAsigma3s_1("PAsigma3s_1","@0*@1",RooArgList(PAsigma1s_1,PAmRatio31) );

  RooRealVar *PAx1s = new RooRealVar("PAx1s","sigma ratio ", PAx1s_init, 0, 1);

  RooFormulaVar PAsigma1s_2("PAsigma1s_2","@0*@1",RooArgList(PAsigma1s_1, *PAx1s) );
  RooFormulaVar PAsigma2s_2("PAsigma2s_2","@0*@1",RooArgList(PAsigma1s_2,PAmRatio21) );
  RooFormulaVar PAsigma3s_2("PAsigma3s_2","@0*@1",RooArgList(PAsigma1s_2,PAmRatio31) );

  RooRealVar PAalpha1s_1("PAalpha1s_1","tail shift", PAalpha1s_1_init, 1.0, 3.321);
  RooFormulaVar PAalpha2s_1("PAalpha2s_1","1.0*@0",RooArgList(PAalpha1s_1) );
  RooFormulaVar PAalpha3s_1("PAalpha3s_1","1.0*@0",RooArgList(PAalpha1s_1) );
  RooFormulaVar PAalpha1s_2("PAalpha1s_2","1.0*@0",RooArgList(PAalpha1s_1) );
  RooFormulaVar PAalpha2s_2("PAalpha2s_2","1.0*@0",RooArgList(PAalpha1s_1) );
  RooFormulaVar PAalpha3s_2("PAalpha3s_2","1.0*@0",RooArgList(PAalpha1s_1) );

  RooRealVar PAn1s_1("PAn1s_1","power order", PAn1s_1_init , 1.416, 3.357);
  RooFormulaVar PAn2s_1("PAn2s_1","1.0*@0",RooArgList(PAn1s_1) );
  RooFormulaVar PAn3s_1("PAn3s_1","1.0*@0",RooArgList(PAn1s_1) );
  RooFormulaVar PAn1s_2("PAn1s_2","1.0*@0",RooArgList(PAn1s_1) );
  RooFormulaVar PAn2s_2("PAn2s_2","1.0*@0",RooArgList(PAn1s_1) );
  RooFormulaVar PAn3s_2("PAn3s_2","1.0*@0",RooArgList(PAn1s_1) );
  
  RooRealVar *PAf1s = new RooRealVar("PAf1s","1S CB fraction", PAf1s_init, 0, 1);
  RooFormulaVar PAf2s("PAf2s","1.0*@0",RooArgList(*PAf1s) );
  RooFormulaVar PAf3s("PAf3s","1.0*@0",RooArgList(*PAf1s) );

  // Set up crystal ball shapes
  RooCBShape* PAcb1s_1 = new RooCBShape("PAcball1s_1", "cystal Ball", *(PAws->var("mass")), PAmean1s, PAsigma1s_1, PAalpha1s_1, PAn1s_1);
  RooCBShape* PAcb2s_1 = new RooCBShape("PAcball2s_1", "cystal Ball", *(PAws->var("mass")), PAmean2s, PAsigma2s_1, PAalpha2s_1, PAn2s_1);
  RooCBShape* PAcb3s_1 = new RooCBShape("PAcball3s_1", "cystal Ball", *(PAws->var("mass")), PAmean3s, PAsigma3s_1, PAalpha3s_1, PAn3s_1);

  RooAddPdf* PAcb1s;
  RooAddPdf* PAcb2s;
  RooAddPdf* PAcb3s;

if (whichModel) {
  //CB+GAUSSIAN
  RooGaussian* PAgauss1s = new RooGaussian("PAgauss1s","gaussian PDF",*(PAws->var("mass")),PAmean1s,PAsigma1s_2);
  RooGaussian* PAgauss2s = new RooGaussian("PAgauss2s","gaussian PDF",*(PAws->var("mass")),PAmean2s,PAsigma2s_2);
  RooGaussian* PAgauss3s = new RooGaussian("PAgauss3s","gaussian PDF",*(PAws->var("mass")),PAmean3s,PAsigma3s_2);
  PAcb1s = new RooAddPdf("PAcb1s","Signal 1S",RooArgList(*PAcb1s_1,*PAgauss1s), RooArgList(*PAf1s) );
  PAcb2s = new RooAddPdf("PAcb2s","Signal 2S",RooArgList(*PAcb2s_1,*PAgauss2s), RooArgList(*PAf1s) );
  PAcb3s = new RooAddPdf("PAcb3s","Signal 3S",RooArgList(*PAcb3s_1,*PAgauss3s), RooArgList(*PAf1s) );
}
else {
  //DOUBLE CRYSTAL BALL
  RooCBShape* PAcb1s_2 = new RooCBShape("PAcball1s_2", "cystal Ball", *(PAws->var("mass")), PAmean1s, PAsigma1s_2, PAalpha1s_2, PAn1s_2);
  RooCBShape* PAcb2s_2 = new RooCBShape("PAcball2s_2", "cystal Ball", *(PAws->var("mass")), PAmean2s, PAsigma2s_2, PAalpha2s_2, PAn2s_2);
  RooCBShape* PAcb3s_2 = new RooCBShape("PAcball3s_2", "cystal Ball", *(PAws->var("mass")), PAmean3s, PAsigma3s_2, PAalpha3s_2, PAn3s_2);
  PAcb1s = new RooAddPdf("PAcb1s","Signal 1S",RooArgList(*PAcb1s_1,*PAcb1s_2), RooArgList(*PAf1s) );
  PAcb2s = new RooAddPdf("PAcb2s","Signal 2S",RooArgList(*PAcb2s_1,*PAcb2s_2), RooArgList(*PAf1s) );
  PAcb3s = new RooAddPdf("PAcb3s","Signal 3S",RooArgList(*PAcb3s_1,*PAcb3s_2), RooArgList(*PAf1s) );
}

  //float PAscalefactor = 1;
  float PAscalefactor = Nomws->var("nSig1s")->getVal()*RpA1/(PANomws->var("nSig1s")->getVal());

  RooRealVar *PAnSig1s= new RooRealVar("PAnSig1s"," 1S signals",0,1000000);
  RooRealVar *PAnSig2s= new RooRealVar("PAnSig2s"," 2S signals",-20,360000);
  RooRealVar *PAnSig3s= new RooRealVar("PAnSig3s"," 3S signals",-50,260000);
  double PAnSig1s_init = PANomws->var("nSig1s")->getVal()*PAscalefactor;
  double PAnSig2s_init = PANomws->var("nSig2s")->getVal()*PAscalefactor;
  double PAnSig3s_init = PANomws->var("nSig3s")->getVal()*PAscalefactor;
  PAnSig1s->setVal(PAnSig1s_init);
  PAnSig2s->setVal(PAnSig2s_init);
  PAnSig3s->setVal(PAnSig3s_init);

  //BACKGROUND
  RooRealVar PAerr_mu("#PAmu","PAerr_mu", err_mu_init,  0, 25) ;
  RooRealVar PAerr_sigma("#PAsigma","PAerr_sigma", err_sigma_init, 0,25);
  RooRealVar PAm_lambda("#PAlambda","PAm_lambda",  m_lambda_init, 0,25);
  
  RooGenericPdf *PAbkg;
  RooGenericPdf *PAbkgLowPt = new RooGenericPdf("PAbkgLowPt","PABackground","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(PAws->var("mass")), PAm_lambda, PAerr_mu, PAerr_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *PAbkgHighPt = new RooGenericPdf("PAbkgHighPt","PABackground","TMath::Exp(-@0/@1)",RooArgList(*(PAws->var("mass")),PAm_lambda));
  
  if  (ptLow >= 5)        PAbkg = PAbkgHighPt ;
  else PAbkg = PAbkgLowPt;

  RooRealVar *PAnBkg = new RooRealVar("PAnBkg","fraction of component 1 in bkg",10000,0,5000000);
  PAnBkg->setVal(nBkg_init);

  //Build the model
  RooAddPdf* PAmodel = new RooAddPdf();
  PAmodel = new RooAddPdf("PAmodel","1S+2S+3S + Bkg",RooArgList(*PAcb1s, *PAcb2s, *PAcb3s, *PAbkg),RooArgList(*PAnSig1s,*PAnSig2s,*PAnSig3s,*PAnBkg));

  PAws->import(*PAmodel);
  //Plot it
  RooPlot* PAmyPlot2 = (RooPlot*)PAmyPlot->Clone();
  PAws->data("reducedDS")->plotOn(PAmyPlot2,Name("PAdataOS_FIT"),MarkerSize(.8));

  PAws->pdf("PAmodel")->plotOn(PAmyPlot2,Name("PAmodelHist"),NormRange("full"));

  PAmyPlot2->Draw()
  PApad1->Update();

  PApad1->Update();
  //PApad2->Update();

  PAc1->cd();
  PApad1->Draw();
  //PApad2->Draw();

  PApad1->Update();
  //PApad2->Update();

//*****************************************************************************
//*****************************************************************************
//             MAKE THE RpA PLOT
//*****************************************************************************
//*****************************************************************************

  //TCanvas* RpAc1 =  new TCanvas("RpAcanvas2","My plots",1004,45,550,520);
  TCanvas* RpAc1 =  new TCanvas("RpAcanvas2","My plots",1004,45,520,570);
  RpAc1->cd();

  //Plot it
  RooPlot* RpAmyPlot2 = (RooPlot*)myPlot2->Clone();
  float AANorm = (AANomws->var("nSig1s")->getVal()*scalefactor) + (AANomws->var("nSig2s")->getVal()*scalefactor) + (AANomws->var("nSig3s")->getVal()*scalefactor) + (Nomws->var("nBkg")->getVal());
  float PANorm = (PANomws->var("nSig1s")->getVal()*PAscalefactor) + (PANomws->var("nSig2s")->getVal()*PAscalefactor) + (PANomws->var("nSig3s")->getVal()*PAscalefactor) + (Nomws->var("nBkg")->getVal());
  PAws->pdf("PAmodel")->plotOn(RpAmyPlot2,Name("PAmodelHistRpA"),LineColor(kGreen+1),LineWidth(2),LineStyle(1),Normalization(PANorm,RooAbsReal::NumEvent),Range(8.8,10.8),NormRange("full"));
  AAws->pdf("AAmodel")->plotOn(RpAmyPlot2,Name("AAmodelHistRpA"),LineColor(kRed),LineWidth(2),LineStyle(2),Normalization(AANorm,RooAbsReal::NumEvent),Range(8.8,10.8),NormRange("full"));

  RpAmyPlot2->Draw();

  //make a pretty plot
  RpAmyPlot2->SetFillStyle(4000);
  RpAmyPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  RpAmyPlot2->GetYaxis()->SetTitle("Events / (0.1 GeV/c^{ 2})");
  RpAmyPlot2->GetYaxis()->SetTitleOffset(1.4);
  RpAmyPlot2->GetYaxis()->CenterTitle();
  RpAmyPlot2->GetYaxis()->SetTitleSize(0.05);
  RpAmyPlot2->GetYaxis()->SetLabelSize(0.05);
  RpAmyPlot2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  RpAmyPlot2->GetXaxis()->SetLabelSize(0.05);
  RpAmyPlot2->GetXaxis()->SetRangeUser(8,14);
  RpAmyPlot2->GetXaxis()->SetTitleSize(0.058);
  RpAmyPlot2->GetXaxis()->CenterTitle();
  RpAmyPlot2->Draw();

  collId = kPADATA;
  yLow = -1.93;
  yHigh = 1.93;
  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.055;
  float text_size = 16;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if (collId==kPPDATA) {
    if(yLow==0) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < |y_{CM}^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  else if (collId==kPADATA) {
    if(yLow==-yHigh) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < y_{CM}^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("|#eta^{#mu}| < 2.4"), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality 0-100%s", perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);

  TLegend* RpAfitleg = new TLegend(0.6,0.3,0.91,0.55); RpAfitleg->SetTextSize(19);
  RpAfitleg->SetTextFont(43);
  RpAfitleg->SetBorderSize(0);
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("dataOS_FIT"),"pp Data","pe");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("modelHist"),"Total Fit","l");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("bkgPDF"),"Background","l");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("PAmodelHistRpA"),"R_{pPb} scaled","l");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("AAmodelHistRpA"),"R_{PbPb} scaled","l");
  RpAfitleg->Draw("same");

  //continue beautifying the plot and print out results
  l1->Draw("same");

  setTDRStyle();
  //writeExtraText = true;
  //extraText = "Preliminary";

  label="";
  CMS_lumi_Jared(RpAc1, 102 ,33);

  RpAc1->cd();

  RpAc1->SaveAs("RpARAAVisualPlot.png");
  RpAc1->SaveAs("RpARAAVisualPlot.pdf");

  cout << "PPnSig1s = " << nSig1s_init << endl;
  cout << "PPnSig2s = " << nSig2s_init << endl;
  cout << "PPnSig3s = " << nSig3s_init << endl;
  cout << "PPnBkg = " << nBkg_init << endl;
  cout << "PAscalefactor = " << PAscalefactor << endl;
  cout << "PAnSig1s = " << PAnSig1s_init << endl;
  cout << "PAnSig2s = " << PAnSig2s_init << endl;
  cout << "PAnSig3s = " << PAnSig3s_init << endl;
  cout << "PAnBkg = " << nBkg_init << endl;
  cout << "AAscalefactor = " << scalefactor << endl;
  cout << "AAnSig1s = " << AAnSig1s_init << endl;
  cout << "AAnSig2s = " << AAnSig2s_init << endl;
  cout << "AAnSig3s = " << AAnSig3s_init << endl;
  cout << "AAnBkg = " << nBkg_init << endl;
} 
 
