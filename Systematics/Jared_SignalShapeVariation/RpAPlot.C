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
void RpAPlot( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
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

//*****************************************************************************
//*****************************************************************************
//             MAKE THE PA PLOT
//*****************************************************************************
//*****************************************************************************

  TFile* f1;
  TFile* f2;
  float yLowLab;
  float yHighLab;
  //Select Data Set
  f1 = new TFile("../../yskimPA1st_OpSign_20177262037_unIdentified.root");
  f2 = new TFile("../../yskimPA2nd_OpSign_20177262044_unIdentified.root");
  yLowLab = yLow+0.47;
  yHighLab = yHigh+0.47;

  //import the model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");

  RooAbsData* reducedDS = ws->data("reducedDS");

  RooRealVar* m_lambda = ws->var("#lambda");
  RooRealVar* err_mu = ws->var("#mu");
  RooRealVar* err_sigma = ws->var("#sigma");
  RooRealVar* alpha1s_1 = ws->var("alpha1s_1");
  RooRealVar* f1s = ws->var("f1s");
  RooRealVar* mRatio21 = ws->var("mRatio21");
  RooRealVar* mRatio31 = ws->var("mRatio31");
  RooRealVar* mean1s = ws->var("m_{#Upsilon(1S)}");
  RooRealVar* n1s_1 = ws->var("n1s_1");
  RooRealVar* sigma1s_1 = ws->var("sigma1s_1");
  RooRealVar* x1s = ws->var("x1s");
  RooRealVar* nSig1s = ws->var("nSig1s");
  RooRealVar* nSig2s = ws->var("nSig2s");
  RooRealVar* nSig3s = ws->var("nSig3s");
  RooRealVar* nBkg = ws->var("nBkg");

  RooAbsPdf* cb1s = ws->pdf("cb1s");
  RooAbsPdf* cb2s = ws->pdf("cb2s");
  RooAbsPdf* cb3s = ws->pdf("cb3s");
  RooAbsPdf* bkg = ws->pdf("bkg");
  RooAbsPdf* model = ws->pdf("model");

  //Plot it
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

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
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
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


//*****************************************************************************
//*****************************************************************************
//             MAKE THE PP PLOT
//*****************************************************************************
//*****************************************************************************

  collId = kPPDATA,

  TFile* f3;
  //Select Data Set
  f3 = new TFile("../../yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20177262158_.root");
  yLow = 0.0;
  yLowLab = yLow;
  yHighLab = yHigh;

  //import the model
  cout << "Importing workspace" << endl;
  TString PPkineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString PPNomFileName = Form("nomfitresults_upsilon_%s.root",PPkineLabel.Data());
  cout << PPNomFileName << endl;
  TFile* PPNomFile = TFile::Open(PPNomFileName,"READ");
  RooWorkspace *PPws = (RooWorkspace*)PPNomFile->Get("workspace");
  PPNomFile->Close("R");

  RooAbsData* PPreducedDS = PPws->data("reducedDS");

  RooRealVar* PPm_lambda = PPws->var("#lambda");
  RooRealVar* PPerr_mu = PPws->var("#mu");
  RooRealVar* PPerr_sigma = PPws->var("#sigma");
  RooRealVar* PPalpha1s_1 = PPws->var("alpha1s_1");
  RooRealVar* PPf1s = PPws->var("f1s");
  RooRealVar* PPmRatio21 = PPws->var("mRatio21");
  RooRealVar* PPmRatio31 = PPws->var("mRatio31");
  RooRealVar* PPmean1s = PPws->var("m_{#Upsilon(1S)}");
  RooRealVar* PPn1s_1 = PPws->var("n1s_1");
  RooRealVar* PPsigma1s_1 = PPws->var("sigma1s_1");
  RooRealVar* PPx1s = PPws->var("x1s");
  RooRealVar* PPnSig1s = PPws->var("nSig1s");
  RooRealVar* PPnSig2s = PPws->var("nSig2s");
  RooRealVar* PPnSig3s = PPws->var("nSig3s");
  RooRealVar* PPnBkg = PPws->var("nBkg");

  RooAbsPdf* PPcb1s = PPws->pdf("cb1s");
  RooAbsPdf* PPcb2s = PPws->pdf("cb2s");
  RooAbsPdf* PPcb3s = PPws->pdf("cb3s");
  RooAbsPdf* PPbkg = PPws->pdf("bkg");
  RooAbsPdf* PPmodel = PPws->pdf("model");

  //Plot it
  TCanvas* PPc1 =  new TCanvas("PPcanvas2","My plots",545,45,550,520);
  PPc1->cd();
  TPad *PPpad1 = new TPad("PPpad1", "PPpad1", 0, 0.25, 0.98, 1.0);
  PPpad1->SetTicks(1,1);
  PPpad1->Draw(); PPpad1->cd();

  RooPlot* PPmyPlot = PPws->var("mass")->frame(nMassBin); // bins
  PPws->data("reducedDS")->plotOn(PPmyPlot,Name("PPdataHist"));

  RooPlot* PPmyPlot2 = (RooPlot*)PPmyPlot->Clone();
  PPws->data("reducedDS")->plotOn(PPmyPlot2,Name("PPdataOS_FIT"),MarkerSize(.8));

  PPws->pdf("model")->plotOn(PPmyPlot2,Name("PPmodelHist"));
  ws->pdf("model")->plotOn(PPmyPlot2,Name("PPmodelHistother"));
  PPws->pdf("model")->plotOn(PPmyPlot2,Name("PPSig1S"),Components(RooArgSet(*PPcb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  PPws->pdf("model")->plotOn(PPmyPlot2,Components(RooArgSet(*PPcb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  PPws->pdf("model")->plotOn(PPmyPlot2,Components(RooArgSet(*PPcb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  PPws->pdf("model")->plotOn(PPmyPlot2,Name("PPbkgPDF"),Components(RooArgSet(*PPbkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  PPmyPlot2->SetFillStyle(4000);
  PPmyPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  PPmyPlot2->GetYaxis()->SetTitleOffset(1.43);
  PPmyPlot2->GetYaxis()->CenterTitle();
  PPmyPlot2->GetYaxis()->SetTitleSize(0.058);
  PPmyPlot2->GetYaxis()->SetLabelSize(0.054);
  PPmyPlot2->GetXaxis()->SetLabelSize(0);
  PPmyPlot2->GetXaxis()->SetRangeUser(8,14);
  PPmyPlot2->GetXaxis()->SetTitleSize(0);
  PPmyPlot2->Draw();
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
  else {
    if (collId==kPPDATA) drawText(Form("%.2f < |y^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pp
    else drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pPb
  }
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S)
  {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }

  TLegend* PPfitleg = new TLegend(0.76,0.4,0.91,0.7); PPfitleg->SetTextSize(19);
  PPfitleg->SetTextFont(43);
  PPfitleg->SetBorderSize(0);
  PPfitleg->AddEntry(PPmyPlot2->findObject("PPdataOS_FIT"),"Data","pe");
  PPfitleg->AddEntry(PPmyPlot2->findObject("PPmodelHist"),"Total fit","l");
  PPfitleg->AddEntry(PPmyPlot2->findObject("PPSig1S"),"Signal","l");
  PPfitleg->AddEntry(PPmyPlot2->findObject("PPbkgPDF"),"Background","l");
  PPfitleg->Draw("same");

  // PULL
  TPad *PPpad2 = new TPad("PPpad2", "pad2", 0, 0.05, 0.98, 0.30);
  PPpad2->SetTopMargin(0); // Upper and lower plot are joined
  PPpad2->SetBottomMargin(0.67);
  PPpad1->SetLeftMargin(0.18);
  PPpad1->SetRightMargin(0.02);
  PPpad2->SetRightMargin(0.02);
  PPpad2->SetLeftMargin(0.18);
  PPpad2->SetTicks(1,1);
  PPpad2->cd();
  
  RooHist* PPhpull = PPmyPlot2->pullHist("PPdataHist","PPmodelHist");
  PPhpull->SetMarkerSize(0.8);
  RooPlot* PPpullFrame = PPws->var("mass")->frame(Title("Pull Distribution")) ;
  PPpullFrame->addPlotable(PPhpull,"P") ;
  PPpullFrame->SetTitleSize(0);
  PPpullFrame->GetYaxis()->SetTitleOffset(0.43) ;
  PPpullFrame->GetYaxis()->SetTitle("Pull") ;
  PPpullFrame->GetYaxis()->SetTitleSize(0.18) ; //19
  PPpullFrame->GetYaxis()->SetLabelSize(0.113) ; // 113
  PPpullFrame->GetYaxis()->SetRangeUser(-3.8,3.8) ;
  PPpullFrame->GetYaxis()->CenterTitle();

  PPpullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  PPpullFrame->GetXaxis()->SetTitleOffset(1.05) ;
  PPpullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  PPpullFrame->GetXaxis()->SetLabelSize(0.20) ; //23
  PPpullFrame->GetXaxis()->SetTitleSize(0.25) ;  //28
  PPpullFrame->GetXaxis()->CenterTitle();

  PPpullFrame->GetYaxis()->SetTickSize(0.04);
  PPpullFrame->GetYaxis()->SetNdivisions(404);
  PPpullFrame->GetXaxis()->SetTickSize(0.03);
  PPpullFrame->Draw() ;

  //continue beautifying the plot and print out results
  TLine *PPl1 = new TLine(massLow,0,massHigh,0);
  PPl1->SetLineStyle(9);
  PPl1->Draw("same");
  PPpad1->Update();

  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  TString label;
  label="";
  if(collId == kPPDATA) CMS_lumi(PPpad1, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(PPpad1, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(PPpad1, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(PPpad1, 21 ,33);

  PPpad1->Update();
  PPpad2->Update();

  PPc1->cd();
  PPpad1->Draw();
  PPpad2->Draw();

  PPpad1->Update();
  PPpad2->Update();

//*****************************************************************************
//*****************************************************************************
//             MAKE THE RpA PLOT
//*****************************************************************************
//*****************************************************************************

  /*RooAddPdf* PPoverlay = new RooAddPdf();
  PPoverlay = new RooAddPdf("PPoverlay","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  ws->import(*PPoverlay);
  ws->pdf("PPoverlay")->plotOn(myPlot2,Name("PPoverlayHist"));
*/

  //get the fitted values
  //RooArgList pars(*model.getParameters(RooArgSet(jj) ) );
  //make fitted function
  //RooArgSet prodSet(model); 
  //RooProduct unNormPdf("fitted Function", "fitted Function", prodSet);
  //now save the fitted function as TF1 object


  //TF1* modelfunc = ws->pdf("model")->asTF(*ws->var("mass"),RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  //clone to sample the PL calculation
  //TF1 * modelfunc = (TF1*) tempfunc->Clone();
  //TCanvas * cfunc = new TCanvas();
  //cfunc->cd();
  //modelfunc->Draw();

  //My crazy creation
  /*c1->cd();
  pad1->cd();
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb1s,*cb2s,*cb3s,*bkg)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  pad1->Update();
  pad2->Update();

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();*/

  RooAbsPdf* PAModel = ws->pdf("model");
  RooAbsPdf* PPModel = PPws->pdf("model");
  RooWorkspace *wsRpA = new RooWorkspace("workspace");
  wsRpA->import(*PAModel);
  wsRpA->import(*PPModel);

  TCanvas* RpAc1 =  new TCanvas("RpAcanvas2","My plots",545,545,550,520);
  RpAc1->cd();
  TPad *RpApad1 = new TPad("RpApad1", "RpApad1", 0, 0.25, 0.98, 1.0);
  RpApad1->SetTicks(1,1);
  RpApad1->Draw(); RpApad1->cd();

  RooPlot* RpAmyPlot = RpAws->var("mass")->frame(nMassBin); // bins
  //RpAws->data("reducedDS")->plotOn(RpAmyPlot,Name("PPdataHist"));

  RooPlot* RpAmyPlot2 = (RooPlot*)RpAmyPlot->Clone();
  //RpAws->data("reducedDS")->plotOn(RpAmyPlot2,Name("RpAdataOS_FIT"),MarkerSize(.8));

  RpAws->pdf("PAmodel")->plotOn(RpAmyPlot2,Name("RpAmodelHist"));

} 
 
