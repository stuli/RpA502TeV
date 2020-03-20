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
void draw_RpAVisualPlot( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  TGaxis::SetMaxDigits(3);

  //Useful parameters
  int A = 208;
  float PPeff = 0.824;//Efficiency values as of 4/2/2018
  float PAeff = 0.795;//EffIntRpA = 1.036
  float AAeff = 0.729;//From 16-023 integrated bin
  float PPlum = 28000;//28 pb^-1
  float PAlum = 34.6;//34.6 nb^-1
  float AAlum = 0.368;//368 microb^-1
  float RAA1 = 0.3518;
  float RpA1 = 0.773;//0.786;
  float RpA2 = 0.673;//0.686;
  float RpA3 = 0.514;//0.518;

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
//             MAKE THE PA MODEL
//*****************************************************************************
//*****************************************************************************

  float yLowLab;
  float yHighLab;

  collId = kPADATA;

  //Select Data Set
  yLowLab = yLow;
  yHighLab = yHigh;

  //import the model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, 0, 100, 0.0, 0.5);
  TString NomFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/nomfitresults_upsilon_%s.root",kineLabel.Data());
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
  float nBkg_init = Nomws->var("nBkg")->getVal();
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
  //TString perc = "%";

  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.075;
  float text_size = 20;
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
  drawText(Form("|#eta_{lab}^{#mu}| < 2.4"), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  //TString perc = "%";
  //drawText(Form("Centrality 0-100%s", perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);

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




//*****************************************************************************
//*****************************************************************************
//             MAKE THE PP MODEL
//*****************************************************************************
//*****************************************************************************

  collId = kPPDATA,

  //Select Data Set
  yLow = 0.0;
  yLowLab = yLow;
  yHighLab = yHigh;

  //import the model
  cout << "Importing workspace" << endl;
  //TString PPkineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  //TString PPNomFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2018_07_18/nomfitresults_upsilon_%s.root",PPkineLabel.Data());
  //cout << PPNomFileName << endl;
  TFile* PPNomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *PPNomws = (RooWorkspace*)PPNomFile->Get("workspace");
  PPNomFile->Close("R");

  RooAbsData* PPreducedDS = PPNomws->data("reducedDS");
  RooWorkspace *PPws = new RooWorkspace("workspace");
  PPws->import(*PPreducedDS);
  PPws->var("mass")->setRange(massLow, massHigh);
  PPws->var("mass")->setVal(PPNomws->var("mass")->getVal());
  PPws->var("mass")->Print();

  //Plot it
  TCanvas* PPc1 =  new TCanvas("PPcanvas2","My plots",504,45,550,550);
  PPc1->cd();
  TPad *PPpad1 = new TPad("PPpad1", "PPpad1", 0, 0.16, 0.98, 1.0);
  PPpad1->SetTicks(1,1);
  PPpad1->Draw(); PPpad1->cd();

  RooPlot* PPmyPlot = PPws->var("mass")->frame(nMassBin); // bins
  PPws->data("reducedDS")->plotOn(PPmyPlot,Name("PPdataHist"));
  double PPmean1s_init = PPNomws->var("m_{#Upsilon(1S)}")->getVal();
  RooRealVar PPmean1s("PPm_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",PPmean1s_init, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar PPmRatio21("PPmRatio21","PPmRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar PPmRatio31("PPmRatio31","PPmRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar PPmean2s("PPmean2s","PPm_{#Upsilon(1S)}*PPmRatio21", RooArgSet(PPmean1s,PPmRatio21) );
  RooFormulaVar PPmean3s("PPmean3s","PPm_{#Upsilon(1S)}*PPmRatio31", RooArgSet(PPmean1s,PPmRatio31) );

  //SIGNAL:
  double PPsigma1s_1_init = PPNomws->var("sigma1s_1")->getVal();
  double PPx1s_init = PPNomws->var("x1s")->getVal();
  double PPalpha1s_1_init = PPNomws->var("alpha1s_1")->getVal();
  double PPn1s_1_init = PPNomws->var("n1s_1")->getVal();
  double PPf1s_init = PPNomws->var("f1s")->getVal();

  RooRealVar    PPsigma1s_1("PPsigma1s_1","width/sigma of the signal gaussian mass PDF",PPsigma1s_1_init, 0.02, 0.3);
  RooFormulaVar PPsigma2s_1("PPsigma2s_1","@0*@1",RooArgList(PPsigma1s_1,PPmRatio21) );
  RooFormulaVar PPsigma3s_1("PPsigma3s_1","@0*@1",RooArgList(PPsigma1s_1,PPmRatio31) );

  RooRealVar *PPx1s = new RooRealVar("PPx1s","sigma ratio ", PPx1s_init, 0, 1);

  RooFormulaVar PPsigma1s_2("PPsigma1s_2","@0*@1",RooArgList(PPsigma1s_1, *PPx1s) );
  RooFormulaVar PPsigma2s_2("PPsigma2s_2","@0*@1",RooArgList(PPsigma1s_2,PPmRatio21) );
  RooFormulaVar PPsigma3s_2("PPsigma3s_2","@0*@1",RooArgList(PPsigma1s_2,PPmRatio31) );

  RooRealVar PPalpha1s_1("PPalpha1s_1","tail shift", PPalpha1s_1_init, 1.0, 3.321);
  RooFormulaVar PPalpha2s_1("PPalpha2s_1","1.0*@0",RooArgList(PPalpha1s_1) );
  RooFormulaVar PPalpha3s_1("PPalpha3s_1","1.0*@0",RooArgList(PPalpha1s_1) );
  RooFormulaVar PPalpha1s_2("PPalpha1s_2","1.0*@0",RooArgList(PPalpha1s_1) );
  RooFormulaVar PPalpha2s_2("PPalpha2s_2","1.0*@0",RooArgList(PPalpha1s_1) );
  RooFormulaVar PPalpha3s_2("PPalpha3s_2","1.0*@0",RooArgList(PPalpha1s_1) );

  RooRealVar PPn1s_1("PPn1s_1","power order", PPn1s_1_init , 1.416, 3.357);
  RooFormulaVar PPn2s_1("PPn2s_1","1.0*@0",RooArgList(PPn1s_1) );
  RooFormulaVar PPn3s_1("PPn3s_1","1.0*@0",RooArgList(PPn1s_1) );
  RooFormulaVar PPn1s_2("PPn1s_2","1.0*@0",RooArgList(PPn1s_1) );
  RooFormulaVar PPn2s_2("PPn2s_2","1.0*@0",RooArgList(PPn1s_1) );
  RooFormulaVar PPn3s_2("PPn3s_2","1.0*@0",RooArgList(PPn1s_1) );
  
  RooRealVar *PPf1s = new RooRealVar("PPf1s","1S CB fraction", PPf1s_init, 0, 1);
  RooFormulaVar PPf2s("PPf2s","1.0*@0",RooArgList(*PPf1s) );
  RooFormulaVar PPf3s("PPf3s","1.0*@0",RooArgList(*PPf1s) );

  // Set up crystal ball shapes
  RooCBShape* PPcb1s_1 = new RooCBShape("PPcball1s_1", "cystal Ball", *(PPws->var("mass")), PPmean1s, PPsigma1s_1, PPalpha1s_1, PPn1s_1);
  RooCBShape* PPcb2s_1 = new RooCBShape("PPcball2s_1", "cystal Ball", *(PPws->var("mass")), PPmean2s, PPsigma2s_1, PPalpha2s_1, PPn2s_1);
  RooCBShape* PPcb3s_1 = new RooCBShape("PPcball3s_1", "cystal Ball", *(PPws->var("mass")), PPmean3s, PPsigma3s_1, PPalpha3s_1, PPn3s_1);

  RooAddPdf* PPcb1s;
  RooAddPdf* PPcb2s;
  RooAddPdf* PPcb3s;

if (whichModel) {
  //CB+GAUSSIAN
  RooGaussian* PPgauss1s = new RooGaussian("PPgauss1s","gaussian PDF",*(PPws->var("mass")),PPmean1s,PPsigma1s_2);
  RooGaussian* PPgauss2s = new RooGaussian("PPgauss2s","gaussian PDF",*(PPws->var("mass")),PPmean2s,PPsigma2s_2);
  RooGaussian* PPgauss3s = new RooGaussian("PPgauss3s","gaussian PDF",*(PPws->var("mass")),PPmean3s,PPsigma3s_2);
  PPcb1s = new RooAddPdf("PPcb1s","Signal 1S",RooArgList(*PPcb1s_1,*PPgauss1s), RooArgList(*PPf1s) );
  PPcb2s = new RooAddPdf("PPcb2s","Signal 2S",RooArgList(*PPcb2s_1,*PPgauss2s), RooArgList(*PPf1s) );
  PPcb3s = new RooAddPdf("PPcb3s","Signal 3S",RooArgList(*PPcb3s_1,*PPgauss3s), RooArgList(*PPf1s) );
}
else {
  //DOUBLE CRYSTAL BALL
  RooCBShape* PPcb1s_2 = new RooCBShape("PPcball1s_2", "cystal Ball", *(PPws->var("mass")), PPmean1s, PPsigma1s_2, PPalpha1s_2, PPn1s_2);
  RooCBShape* PPcb2s_2 = new RooCBShape("PPcball2s_2", "cystal Ball", *(PPws->var("mass")), PPmean2s, PPsigma2s_2, PPalpha2s_2, PPn2s_2);
  RooCBShape* PPcb3s_2 = new RooCBShape("PPcball3s_2", "cystal Ball", *(PPws->var("mass")), PPmean3s, PPsigma3s_2, PPalpha3s_2, PPn3s_2);
  PPcb1s = new RooAddPdf("PPcb1s","Signal 1S",RooArgList(*PPcb1s_1,*PPcb1s_2), RooArgList(*PPf1s) );
  PPcb2s = new RooAddPdf("PPcb2s","Signal 2S",RooArgList(*PPcb2s_1,*PPcb2s_2), RooArgList(*PPf1s) );
  PPcb3s = new RooAddPdf("PPcb3s","Signal 3S",RooArgList(*PPcb3s_1,*PPcb3s_2), RooArgList(*PPf1s) );
}

  //Scale the upsilon peaks
  RooRealVar *PPnSig1s= new RooRealVar("PPnSig1s"," 1S signals",0,1000000);
  RooRealVar *PPnSig2s= new RooRealVar("PPnSig2s"," 2S signals",-20,360000);
  RooRealVar *PPnSig3s= new RooRealVar("PPnSig3s"," 3S signals",-50,260000);
  double PPnSig1s_init = PPNomws->var("nSig1s")->getVal()*1/RpA1;
  double PPnSig2s_init = PPNomws->var("nSig2s")->getVal()*1/RpA2;
  double PPnSig3s_init = PPNomws->var("nSig3s")->getVal()*1/RpA3;
  PPnSig1s->setVal(PPnSig1s_init);
  PPnSig2s->setVal(PPnSig2s_init);
  PPnSig3s->setVal(PPnSig3s_init);

  //BACKGROUND
  RooRealVar PPerr_mu("#PPmu","PPerr_mu", err_mu_init,  0, 25) ;
  RooRealVar PPerr_sigma("#PPsigma","PPerr_sigma", err_sigma_init, 0,25);
  RooRealVar PPm_lambda("#PPlambda","PPm_lambda",  m_lambda_init, 0,25);
  
  RooGenericPdf *PPbkg;
  RooGenericPdf *PPbkgLowPt = new RooGenericPdf("PPbkgLowPt","PPBackground","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(PPws->var("mass")), PPm_lambda, PPerr_mu, PPerr_sigma) );

  //THIS IS THE HIGH-PT BACKGROUND FUNCTION
  RooGenericPdf *PPbkgHighPt = new RooGenericPdf("PPbkgHighPt","PPBackground","TMath::Exp(-@0/@1)",RooArgList(*(PPws->var("mass")),PPm_lambda));
  
  if  (ptLow >= 5)        PPbkg = PPbkgHighPt ;
  else PPbkg = PPbkgLowPt;

  RooRealVar *PPnBkg = new RooRealVar("PPnBkg","fraction of component 1 in bkg",10000,0,5000000);
  PPnBkg->setVal(nBkg_init);

  //Build the model
  RooAddPdf* PPmodel = new RooAddPdf();
  PPmodel = new RooAddPdf("PPmodel","1S+2S+3S + Bkg",RooArgList(*PPcb1s, *PPcb2s, *PPcb3s, *PPbkg),RooArgList(*PPnSig1s,*PPnSig2s,*PPnSig3s,*PPnBkg));

  PPws->import(*PPmodel);
  //Plot it
  RooPlot* PPmyPlot2 = (RooPlot*)PPmyPlot->Clone();
  PPws->data("reducedDS")->plotOn(PPmyPlot2,Name("PPdataOS_FIT"),MarkerSize(.8));

  PPws->pdf("PPmodel")->plotOn(PPmyPlot2,Name("PPmodelHist"),NormRange("full"));



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
  float PPNorm = (PPNomws->var("nSig1s")->getVal()*1/RpA1) + (PPNomws->var("nSig2s")->getVal()*1/RpA2) + (PPNomws->var("nSig3s")->getVal()*1/RpA3) + (Nomws->var("nBkg")->getVal());
  PPws->pdf("PPmodel")->plotOn(RpAmyPlot2,Name("PPmodelHistRpA"),LineColor(kRed),LineWidth(2),LineStyle(2),Normalization(PPNorm,RooAbsReal::NumEvent),Range(8.6,10.8),NormRange("full"));

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

  //Write kinematic cuts
  collId = kPADATA;
  yLow = -1.93;
  yHigh = 1.93;
  pos_text_x = 0.43;
  pos_text_y = 0.816;
  pos_y_diff = 0.075;
  text_size = 19;
  text_color = 1;
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
  drawText(Form("|#eta_{lab}^{#mu}| < 2.4"), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);

  //Draw legend
  TLegend* RpAfitleg = new TLegend(0.6,0.3,0.91,0.55); RpAfitleg->SetTextSize(19);
  RpAfitleg->SetTextFont(43);
  RpAfitleg->SetBorderSize(0);
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("dataOS_FIT"),"pPb data","pe");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("modelHist"),"Total fit","l");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("bkgPDF"),"Background","l");
  RpAfitleg->AddEntry(RpAmyPlot2->findObject("PPmodelHistRpA"),"R_{pPb} scaled","l");
  RpAfitleg->Draw("same");

  //l1->Draw("same");

  //Apply some CMS_lumi style
  setTDRStyle();
  writeExtraText = false;
  extraText = "Preliminary";

  label="";
  //CMS_lumi(RpAc1, 101 ,33);
  CMS_lumi(RpAc1, 3 ,33);

  //Save
  RpAc1->cd();

  RpAc1->SaveAs("plots/RpAVisualPlot.png");
  RpAc1->SaveAs("plots/RpAVisualPlot.pdf");

  cout << "PAnSig1s = " << nSig1s_init << endl;
  cout << "PAnSig2s = " << nSig2s_init << endl;
  cout << "PAnSig3s = " << nSig3s_init << endl;
  cout << "PAnBkg = " << nBkg_init << endl;

} 
 
