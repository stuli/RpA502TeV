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

TString kineLabel, NomFileName;
TString directory = "./";

using namespace std;
using namespace RooFit;
void PlotFitOldWay_noPull( 
       int collId = kPPDATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=0.0, float yHigh=1.93,//Run 1 has p going in -z direction
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  TGaxis::SetMaxDigits(3);

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

  float yLowLab;
  float yHighLab;

  //Select Data Set
  yLowLab = yLow;
  yHighLab = yHigh;

  bool isHFBin;
  if ((hfLow>0) || (hfHigh<120) || (ntracksLow>0) || (ntracksHigh<400)) isHFBin = kTRUE;

  //import the model
  cout << "Importing workspace" << endl;
  kineLabel = Form("%s_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, muPtCut) ;
  if (isHFBin) kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
  NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
  //NomFileName = Form("/home/jared/Documents/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/nomfitresults_upsilon_%s.root",kineLabel.Data());
  //NomFileName = Form("FitsWithConstraints/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace* workspace = (RooWorkspace*)NomFile->Get("workspace");
  workspace->Print();

  setTDRStyle();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",1004,45,520,570);
  c1->cd();

  //Plot it
  c1->cd();

  RooPlot* myPlot = workspace->var("mass")->frame(nMassBin); // bins
  workspace->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  workspace->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  workspace->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  RooAbsPdf* cb1s = workspace->pdf("cb1s");
  RooAbsPdf* cb2s = workspace->pdf("cb2s");
  RooAbsPdf* cb3s = workspace->pdf("cb3s");
  RooAbsPdf* bkg;
  if (ptLow<5) bkg = workspace->pdf("bkgLowPt");
  else bkg = workspace->pdf("bkgHighPt");
  workspace->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(7));
  workspace->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(7));
  workspace->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(7));
  workspace->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitle("Events / (0.1 GeV/c^{ 2})");
  myPlot2->GetYaxis()->SetTitleOffset(1.4);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.05);
  myPlot2->GetYaxis()->SetLabelSize(0.05);
  myPlot2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot2->GetXaxis()->SetLabelSize(0.05);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0.058);
  myPlot2->GetXaxis()->CenterTitle();
  myPlot2->Draw();
  RooFitResult* fitRes2 = (RooFitResult*)workspace->obj("fitresult_model_reducedDS");
  fitRes2->Print("v");
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
    if(yLow==0) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < |y_{CM}^{#mu#mu}| < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  else if (collId==kPADATA) {
    if(yLow==-yHigh) drawText(Form("|y_{CM}^{#mu#mu}| < %.2f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    else drawText(Form("%.2f < y_{CM}^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
    }
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  drawText(Form("|#eta_{lab}^{#mu}| < 2.4"), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);

  TLegend* fitleg = new TLegend(0.6,0.3,0.91,0.55); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"pp data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  //fitleg->AddEntry(myPlot2->findObject("Sig1S"),"Signal","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");

  writeExtraText = false;
  extraText = "Preliminary";

  TString label;
  label="";
  if(collId == kPPDATA) CMS_lumi(c1, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(c1, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(c1, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(c1, 21 ,33);

  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");

  float temp1 = workspace->var("nSig1s")->getVal();  
  float temp1err = workspace->var("nSig1s")->getError();  
  float temp2 = workspace->var("nSig2s")->getVal();  
  float temp2err = workspace->var("nSig2s")->getError();  
  float temp3 = workspace->var("nSig3s")->getVal();  
  float temp3err = workspace->var("nSig3s")->getError();  
  
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
  cout <<  workspace->var("n1s_1")->getVal() << ", " <<  workspace->var("alpha1s_1")->getVal() << ", "<<  workspace->var("sigma1s_1")->getVal() << ", " ;
  cout <<  workspace->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  workspace->var("f1s")->getVal() << ", "<<  workspace->var("x1s")->getVal() << " );} " << endl;

  c1->SaveAs(Form("%snomfitresults_upsilon_%s.png",directory.Data(),kineLabel.Data()));
  c1->SaveAs(Form("%snomfitresults_upsilon_%s.pdf",directory.Data(),kineLabel.Data()));

  delete outh;

  NomFile->Close();
} 
 
