//Author: Jared Jay
//This code is based on the fitting code in FitData.C.
//This code generates pseudo-data from the nominal fit (genModel), and then fits it with the nominal fit and the alternative fit. It then compares the estimated upsilon yields from each fit and plots the percent difference for each of the upsilons 1S, 2S, and 3S. It repeats this process as many times as you want. It also plots the most recent fit, and all the plots are updated after every fit.
  
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
void FitPseudoData_Constrained( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=6, 
       float yLow=-1.93, float yHigh=1.93,
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       const int numtrials = 100 
			)
{
 
  int binmode = 0;
  if (hfHigh-hfLow<120) binmode = 1;
  else if (ntracksHigh-ntracksLow<400) binmode = 2;
 
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8; 
  float massHigh = 14;
  int nMassBin  = (massHigh-massLow)*10;

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  
  //import generating model
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  if (binmode>0) kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
  TString NomFileName = Form("../Fits/InternalConstraintFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  
  RooAbsPdf* genModel = Nomws->pdf("model");
  RooWorkspace *wsgen = new RooWorkspace("workspace");
  wsgen->import(*genModel);
  
  double ups1smass = Nomws->var("m_{#Upsilon(1S)}")->getVal();
  double sigma1s_1_init = Nomws->var("sigma1s_1")->getVal();
  double x1s_init = Nomws->var("x1s")->getVal();
  double alpha1s_1_init = Nomws->var("alpha1s_1")->getVal();
  double n1s_1_init = Nomws->var("n1s_1")->getVal();
  double f1s_init = Nomws->var("f1s")->getVal();
  double nSig1s_init = Nomws->var("nSig1s")->getVal();
  double nSig2s_init = Nomws->var("nSig2s")->getVal();
  double nSig3s_init = Nomws->var("nSig3s")->getVal();
  double err_mu_init = 8;
  double err_sigma_init = 5;
  if (ptLow<5) {
    err_mu_init = Nomws->var("#mu")->getVal();
    err_sigma_init = Nomws->var("#sigma")->getVal();
  }
  double m_lambda_init = Nomws->var("#lambda")->getVal();
  double nBkg_init = Nomws->var("nBkg")->getVal();

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
    TString kineLabelpp = getKineLabel (kPPDATA, ptLow, ptHigh, yLowpp, yHighpp, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString ppFileName = Form("../Fits/InternalConstraintFits/nomfitresults_upsilon_%s.root",kineLabelpp.Data());
    cout << ppFileName << endl;
    TFile* ppFile = TFile::Open(ppFileName,"READ");
    RooWorkspace *ppws = (RooWorkspace*)ppFile->Get("workspace");
    ppFile->Close("R");
    float ffix = ppws->var("f1s")->getVal();
    float fdev = ppws->var("f1s")->getError();
    delete ppws;
    delete ppFile;
  }

  //create histograms
  TH1F* histo1s = new TH1F("1SDiff","1S %Diff in Yield",100,-100,100);
  TH1F* histo2s = new TH1F("2SDiff","2S %Diff in Yield",100,-100,100);
  TH1F* histo3s = new TH1F("3SDiff","3S %Diff in Yield",100,-100,100);
  TCanvas* cyields =  new TCanvas("canvas3","results",4,45,1100,400);
  cyields->Divide(3,1);

  cyields->cd(1);
  histo1s->SetXTitle("%Diff");
  histo1s->GetXaxis()->SetTitleSize(0.05);
  histo1s->GetYaxis()->SetLabelSize(0.05);
  histo1s->GetXaxis()->SetLabelSize(0.05);
  histo1s->SetStats(kTRUE);
  histo1s->Draw();

  cyields->cd(2);
  histo2s->SetXTitle("%Diff");
  histo2s->GetXaxis()->SetTitleSize(0.05);
  histo2s->GetYaxis()->SetLabelSize(0.05);
  histo2s->GetXaxis()->SetLabelSize(0.05);
  histo2s->Draw();

  cyields->cd(3);
  histo3s->SetXTitle("%Diff");
  histo3s->GetXaxis()->SetTitleSize(0.05);
  histo3s->GetYaxis()->SetLabelSize(0.05);
  histo3s->GetXaxis()->SetLabelSize(0.05);
  histo3s->Draw();

  TCanvas* cfit =  new TCanvas("canvas2","fitted data",600,500,550,520);

  TString outFileName = Form("PseudoExpResults_%s.root",kineLabel.Data());
  TFile outfile (outFileName, "RECREATE");
  TNtuple* ntuple = new TNtuple("ntuple","Data from fits","chisqNom:yield1sNom:yield2sNom:yield3sNom:chisqAlt:yield1sAlt:yield2sAlt:yield3sAlt:diff1s:diff2s:diff3s",numtrials);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.30);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad1->SetLeftMargin(0.18);
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  pad2->SetLeftMargin(0.18);
  pad2->SetTicks(1,1);

  int numrejected = 0;
  for (int itrial = 0; itrial<numtrials; itrial++) {

    cout << "Number of fits rejected = " << numrejected << endl;

    float chisqtest[2] = {0};
    float yield1s[2] = {0};
    float yield2s[2] = {0};
    float yield3s[2] = {0};

    cout << "*****************************************************" << endl;
    cout << "STARTING TRIAL " << itrial+1 << " OF " << numtrials << endl;
    cout << "*****************************************************" << endl;

    //Generate fake data from the model
    RooDataSet* reducedDS = genModel->generate(*(wsgen->var("mass")));
    reducedDS->SetName("reducedDS");

    //loop through the two fitting models
    for (int imodel=0; imodel<=1; imodel++){
  
      RooWorkspace *ws = new RooWorkspace("workspace");
      ws->import(*reducedDS);
  
      cfit->cd();
      pad1->Draw(); pad1->cd();
    
      RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",ups1smass, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
      RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
      RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
      RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
      RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
  
      //SIGNAL:
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
  
      if (imodel) {
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
  
      RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",nSig1s_init,0,1000000);
      RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",nSig2s_init,-20,360000);
      RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",nSig3s_init,-50,260000);
  
      //BACKGROUND
      RooRealVar err_mu("#mu","err_mu", err_mu_init, paramslower[5], paramsupper[5]) ;
      RooRealVar err_sigma("#sigma","err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
      RooRealVar m_lambda("#lambda","m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);
      
      RooGenericPdf *bkg;
      RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
    
      //THIS IS THE HIGH-PT BACKGROUND FUNCTION
      RooGenericPdf *bkgHighPt = new RooGenericPdf("bkgHighPt","Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda));
      
      if  (ptLow >= 5)        bkg = bkgHighPt ;
      else bkg = bkgLowPt;
    
      RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",nBkg_init,0,5000000);

      // Construct Gaussian constraints on parameters
  float nfix, xfix, alphafix, ffix;
  float ndev, xdev, alphadev, fdev;
  if (collId==kPADATA) {
    if (abs((yLow+yHigh)/2+0.01)<1) { //mid rapidity
      alphafix = 2.047006;
      ffix = 0.199500;
      nfix = 3.235437;
      xfix = 0.498710;
      alphadev = 0.696588;
      fdev = 0.095662;
      ndev = 1.197137;
      xdev = 0.052111;
    }
    else { //extreme rapidity
      alphafix = 1.725369;
      ffix = 0.473512;
      nfix = 1.353564;
      xfix = 0.588294;
      alphadev = 0.241747;
      fdev = 0.099628;
      ndev = 0.149906;
      xdev = 0.040675;
    }
  }
  else if (collId==kPPDATA) {
    if (abs((yLow+yHigh)/2+0.01)<1) {
      //normal
      nfix = 2.29348;
      xfix = 0.542762;
      alphafix = 2.48620;
      ndev = 0.451840;
      xdev = 0.0585783;
      alphadev = 0.395293;
    }
    else {
    //extreme rapidity bins
      nfix = 1.42638;
      xfix = 0.556545;
      alphafix = 1.97426;
      ndev = 0.0518141;
      xdev = 0.0096801;
      alphadev = 0.0205957;
    }
  }
  //n1s_1.setVal(nfix);
  //x1s->setVal(xfix);
  //alpha1s_1.setVal(alphafix);

  RooGaussian nconstraint("nconstraint","nconstraint", n1s_1,RooConst(nfix),RooConst(ndev));
  RooGaussian alphaconstraint("alphaconstraint","alphaconstraint", alpha1s_1,RooConst(alphafix),RooConst(alphadev));
  RooGaussian xconstraint("xconstraint","xconstraint", *x1s,RooConst(xfix),RooConst(xdev));
  RooGaussian fconstraint("fconstraint","fconstraint", *f1s,RooConst(ffix),RooConst(fdev));

  RooArgSet *allConstraints;
  RooArgSet *constParams;
  if (collId==kPADATA) {
    allConstraints = new RooArgSet(nconstraint,alphaconstraint,xconstraint,fconstraint);
    constParams = new RooArgSet(n1s_1,alpha1s_1,*x1s,*f1s);
    //constParams = new RooArgSet(n1s_1,*x1s);
  }
  else if (collId==kPPDATA) {
    allConstraints = new RooArgSet(nconstraint,alphaconstraint,xconstraint);
    constParams = new RooArgSet(n1s_1,alpha1s_1,*x1s);
  }
 
      //Build the model
      RooAddPdf* modelNoC = new RooAddPdf("modelNoC","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkg),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
      RooProdPdf model("model","model with constraint",RooArgSet(*modelNoC,*allConstraints));
      ws->import(model);
    
      RooPlot* myPlot2 = ws->var("mass")->frame(nMassBin);
      ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));
    
      //Fit the model to the data
      RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Constrain(*constParams),Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
      ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
      ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
      ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
      ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
      ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkgLowPt)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
    
      //make a pretty plot
      myPlot2->SetFillStyle(4000);
      myPlot2->SetAxisRange(massLow, massHigh,"X");
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
    
      //Write cuts on plot
      float pos_text_x = 0.43;
      float pos_text_y = 0.816;
      float pos_y_diff = 0.056;
      float text_size = 19;
      int text_color = 1;
      drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
      drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
      drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
      if (imodel>0){
        drawText(Form("Run %i Alternate", itrial+1 ), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
      }
      else {
        drawText(Form("Run %i Nominal", itrial+1 ), pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
      }
    
      //Make a legend
      TLegend* fitleg = new TLegend(0.76,0.4,0.91,0.7); fitleg->SetTextSize(19);
      fitleg->SetTextFont(43);
      fitleg->SetBorderSize(0);
      fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
      fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
      fitleg->AddEntry(myPlot2->findObject("Sig1S"),"Signal","l");
      fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
      fitleg->Draw("same");
    
      // PULL
      pad2->cd();
      
      RooHist* hpull = myPlot2->pullHist("dataOS_FIT","modelHist");
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
      chisqtest[imodel] = chisq/ndf;
      cout << "chisq/dof = " << chisqtest[imodel] << endl;
    
      //Draw the zero line in the pull
      TLine *l1 = new TLine(massLow,0,massHigh,0);
      l1->SetLineStyle(9);
      l1->Draw("same");
      pad1->Update();
    
      //Write "Preliminary" on the plot
      setTDRStyle();
      writeExtraText = true;
      extraText = "Preliminary";
    
      //Write "pPb (5.02 TeV)" on the plot
      if(collId == kPPDATA) CMS_lumi(pad1, 1 ,33);
      else if(collId == kAADATA && cLow < 60) CMS_lumi(pad1, 2 ,33);
      else if(collId == kPADATA) CMS_lumi(pad1, 3 ,33);
      else if(collId == kAADATA && cLow>=60) CMS_lumi(pad1, 21 ,33);
    
      //Update canvases
      pad1->Update();
      pad2->Update();
    
      cfit->cd();
      pad1->Draw();
      pad2->Draw();
    
      pad1->Update();
      pad2->Update();
    
      //Output results  
      float temp1 = ws->var("nSig1s")->getVal();  
      float temp1err = ws->var("nSig1s")->getError();  
      float temp2 = ws->var("nSig2s")->getVal();  
      float temp2err = ws->var("nSig2s")->getError();  
      float temp3 = ws->var("nSig3s")->getVal();  
      float temp3err = ws->var("nSig3s")->getError();  
      
      cout << "1S signal    =  " << temp1 << " +/- " << temp1err << endl;
      cout << "2S signal    =  " << temp2 << " +/- " << temp2err << endl;
      cout << "3S signal    =  " << temp3 << " +/- " << temp3err << endl;
    
      yield1s[imodel] = temp1;
      yield2s[imodel] = temp2;
      yield3s[imodel] = temp3;
    
    }//end of model fit loop
  
    //Reject this trial if the fit is bad.
    if (chisqtest[0]>10 || chisqtest[1]>10) {
      itrial--;
      cout << "MOST RECENT TRIAL REJECTED DUE TO BAD FIT." << endl;
      numrejected++;
      continue;
    }
  
    //Record results if the fits are good
    Float_t chisqtestNom = chisqtest[0];
    Float_t chisqtestAlt = chisqtest[1];
    Float_t yield1sNom = yield1s[0];
    Float_t yield1sAlt = yield1s[1];
    Float_t yield2sNom = yield2s[0];
    Float_t yield2sAlt = yield2s[1];
    Float_t yield3sNom = yield3s[0];
    Float_t yield3sAlt = yield3s[1];
    Float_t perdif1s = 100*(yield1sAlt-yield1sNom)/yield1sNom;
    Float_t perdif2s = 100*(yield2sAlt-yield2sNom)/yield2sNom;
    Float_t perdif3s = 100*(yield3sAlt-yield3sNom)/yield3sNom;
    histo1s->Fill(perdif1s);
    histo2s->Fill(perdif2s);
    histo3s->Fill(perdif3s);
    ntuple->Fill(chisqtestNom,yield1sNom,yield2sNom,yield3sNom,chisqtestAlt,yield1sAlt,yield2sAlt,yield3sAlt,perdif1s,perdif2s,perdif3s);
    cout << "nominal chi^2 = " << chisqtestNom << endl;
    cout << "yield1sNom = " << yield1sNom << endl;
    cout << "yield2sNom = " << yield2sNom << endl;
    cout << "yield3sNom = " << yield3sNom << endl;
    cout << "alternate chi^2 = " << chisqtestAlt << endl;
    cout << "yield1sAlt = " << yield1sAlt << endl;
    cout << "yield2sAlt = " << yield2sAlt << endl;
    cout << "yield3sAlt = " << yield3sAlt << endl;
    cout << "diff1s = " << perdif1s << endl;
    cout << "diff2s = " << perdif2s << endl;
    cout << "diff3s = " << perdif3s << endl;
    cout << "*****************************************************" << endl;
    cout << "TRIAL " << itrial+1 << " OF " << numtrials << " COMPLETED." << endl;
    cout << "*****************************************************" << endl;
  
    cyields->Update();
    cyields->cd(1);
    histo1s->Draw();
    cyields->cd(2);
    histo2s->Draw();
    cyields->cd(3);
    histo3s->Draw();
  
  }//end of main loop
  
  //save histograms
  outfile.cd();
  histo1s->Write();
  histo2s->Write();
  histo3s->Write();
  ntuple->Write();
  outfile.Close();

}
