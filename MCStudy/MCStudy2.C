//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the background shape.

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

#include "RooMCStudy.h"


using namespace std;
using namespace RooFit;
void MCStudy2( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
	   int numTrials = 100,
	   int mcbins = 40, //number of bins in param plots
	   bool binmc = false //use binning in MC generation?
			) 
{
  
	int cLow = 0;
	int cHigh = 200;
	float muPtCut=4.0;

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

	TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  
  
  
  //File for nominal fit
	TString NomFileName = Form("../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
	cout << NomFileName << endl;
	TFile* NomFile = new TFile(NomFileName,"READ");
	if (NomFile->IsZombie())
	{
		cout << "NOMINAL FIT FILE NOT FOUND" << endl;
		cout << "ABORTING" << endl;
		return;
	}
	RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
	RooAbsPdf* genModel = Nomws->pdf("model");
	RooWorkspace *wsgen = new RooWorkspace("workspace");
	wsgen->import(*genModel);
	NomFile->Close();
  
  
	/******************/
	/*  TOY MC STUDY  */
	/******************/
	
	RooMCStudy* mcstudy;
	if (binmc)
		mcstudy = new RooMCStudy(*genModel,*(Nomws->var("mass")),Binned(kTRUE),Silence(),Extended(),FitOptions(Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE)));
	else
		mcstudy = new RooMCStudy(*genModel,*(Nomws->var("mass")),Silence(),Extended(),FitOptions(Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE)));
	
	mcstudy->generateAndFit(numTrials);
	
	RooPlot* frameMC[9];
	frameMC[0] = mcstudy->plotParam(*(Nomws->var("nSig1s")),Bins(mcbins));
	frameMC[1] = mcstudy->plotError(*(Nomws->var("nSig1s")),Bins(mcbins));
	frameMC[2] = mcstudy->plotPull(*(Nomws->var("nSig1s")),Bins(mcbins),FitGauss(kTRUE));
	frameMC[3] = mcstudy->plotParam(*(Nomws->var("nSig2s")),Bins(mcbins));
	frameMC[4] = mcstudy->plotError(*(Nomws->var("nSig2s")),Bins(mcbins));
	frameMC[5] = mcstudy->plotPull(*(Nomws->var("nSig2s")),Bins(mcbins),FitGauss(kTRUE));
	frameMC[6] = mcstudy->plotParam(*(Nomws->var("nSig3s")),Bins(mcbins));
	frameMC[7] = mcstudy->plotError(*(Nomws->var("nSig3s")),Bins(mcbins));
	frameMC[8] = mcstudy->plotPull(*(Nomws->var("nSig3s")),Bins(mcbins),FitGauss(kTRUE));
	
	//gStyle->SetOptStat(0) ;
	TCanvas* c1 = new TCanvas("MCcanvas","MCcanvas",900,900) ;
	c1->Divide(3,3) ;
	for (int i=0; i<9; i++)
	{
		c1->cd(i+1);
		frameMC[i]->Draw();
	}
	
	gDirectory->Add(mcstudy);
	
	TString collIdLabel;
	if (collId == kPADATA) collIdLabel = "PA"; else if (collId == kPPDATA) collIdLabel = "PP";
	TString outName = Form("Results/mcstudy_") + collIdLabel + Form("_%.1f-%.1f_y%.2f-%.2f",ptLow,ptHigh,yLow,yHigh);
	TFile* outfile = new TFile(outName+".root","recreate");
	
	mcstudy->Write();
	c1->Write();
	wsgen->Write();
} 
 
