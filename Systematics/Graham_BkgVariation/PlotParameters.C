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

void PlotParameters()
{
	int cLow=0;
	int cHigh=200;
    float muPtCut=4.0;
	float dphiEp2Low = 0;
	float dphiEp2High = 100;
	
	float ptbins1sLow[4] = {0,2,4,6};
	float ptbins1sHigh[4] = {6,9,12,30};
	float ybins1s[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ptbins2sLow[2] = {0,4};
	float ptbins2sHigh[3] = {4,9,30};
	float ybins2s[5] = {-1.93,-0.8,0.0,0.8,1.93};
	float ptbins3sLow[2] = {0,6};
	float ptbins3sHigh[2] = {6,30};
	float ybins3s[3] = {-1.93,0.0,1.93};
	
	//////////Declare histograms
	
	/////Parameters of low pT bkg function
	//PA
	TH1D* hA1_pt1s_PA = new TH1D("hA1_pt1s_PA","A1 vs p_{T}",3,ptbins1sLow);
	TH1D* hA1_pt2s_PA = new TH1D("hA1_pt2s_PA","A1 vs p_{T}",1,ptbins2sLow);
	TH1D* hA1_pt3s_PA = new TH1D("hA1_pt3s_PA","A1 vs p_{T}",1,ptbins3sLow);
	TH1D* hA2_pt1s_PA = new TH1D("hA2_pt1s_PA","A2 vs p_{T}",3,ptbins1sLow);
	TH1D* hA2_pt2s_PA = new TH1D("hA2_pt2s_PA","A2 vs p_{T}",1,ptbins2sLow);
	TH1D* hA2_pt3s_PA = new TH1D("hA2_pt3s_PA","A2 vs p_{T}",1,ptbins3sLow);
	TH1D* hA3_pt1s_PA = new TH1D("hA3_pt1s_PA","A3 vs p_{T}",3,ptbins1sLow);
	TH1D* hA3_pt2s_PA = new TH1D("hA3_pt2s_PA","A3 vs p_{T}",1,ptbins2sLow);
	TH1D* hA3_pt3s_PA = new TH1D("hA3_pt3s_PA","A3 vs p_{T}",1,ptbins3sLow);
	TH1D* hA4_pt1s_PA = new TH1D("hA4_pt1s_PA","A4 vs p_{T}",3,ptbins1sLow);
	TH1D* hA4_pt2s_PA = new TH1D("hA4_pt2s_PA","A4 vs p_{T}",1,ptbins2sLow);
	TH1D* hA4_pt3s_PA = new TH1D("hA4_pt3s_PA","A4 vs p_{T}",1,ptbins3sLow);
	
	TH1D* hA1_y1s_PA = new TH1D("hA1_y1s_PA","A1 vs y_{cm}",8,ybins1s);
	TH1D* hA1_y2s_PA = new TH1D("hA1_y2s_PA","A1 vs y_{cm}",4,ybins2s);
	TH1D* hA1_y3s_PA = new TH1D("hA1_y3s_PA","A1 vs y_{cm}",2,ybins3s);
	TH1D* hA2_y1s_PA = new TH1D("hA2_y1s_PA","A2 vs y_{cm}",8,ybins1s);
	TH1D* hA2_y2s_PA = new TH1D("hA2_y2s_PA","A2 vs y_{cm}",4,ybins2s);
	TH1D* hA2_y3s_PA = new TH1D("hA2_y3s_PA","A2 vs y_{cm}",2,ybins3s);
	TH1D* hA3_y1s_PA = new TH1D("hA3_y1s_PA","A3 vs y_{cm}",8,ybins1s);
	TH1D* hA3_y2s_PA = new TH1D("hA3_y2s_PA","A3 vs y_{cm}",4,ybins2s);
	TH1D* hA3_y3s_PA = new TH1D("hA3_y3s_PA","A3 vs y_{cm}",2,ybins3s);
	TH1D* hA4_y1s_PA = new TH1D("hA4_y1s_PA","A4 vs y_{cm}",8,ybins1s);
	TH1D* hA4_y2s_PA = new TH1D("hA4_y2s_PA","A4 vs y_{cm}",4,ybins2s);
	TH1D* hA4_y3s_PA = new TH1D("hA4_y3s_PA","A4 vs y_{cm}",2,ybins3s);
	
	//PP
	TH1D* hA1_pt1s_PP = new TH1D("hA1_pt1s_PP","A1 vs p_{T}",3,ptbins1sLow);
	TH1D* hA1_pt2s_PP = new TH1D("hA1_pt2s_PP","A1 vs p_{T}",1,ptbins2sLow);
	TH1D* hA1_pt3s_PP = new TH1D("hA1_pt3s_PP","A1 vs p_{T}",1,ptbins3sLow);
	TH1D* hA2_pt1s_PP = new TH1D("hA2_pt1s_PP","A2 vs p_{T}",3,ptbins1sLow);
	TH1D* hA2_pt2s_PP = new TH1D("hA2_pt2s_PP","A2 vs p_{T}",1,ptbins2sLow);
	TH1D* hA2_pt3s_PP = new TH1D("hA2_pt3s_PP","A2 vs p_{T}",1,ptbins3sLow);
	TH1D* hA3_pt1s_PP = new TH1D("hA3_pt1s_PP","A3 vs p_{T}",3,ptbins1sLow);
	TH1D* hA3_pt2s_PP = new TH1D("hA3_pt2s_PP","A3 vs p_{T}",1,ptbins2sLow);
	TH1D* hA3_pt3s_PP = new TH1D("hA3_pt3s_PP","A3 vs p_{T}",1,ptbins3sLow);
	TH1D* hA4_pt1s_PP = new TH1D("hA4_pt1s_PP","A4 vs p_{T}",3,ptbins1sLow);
	TH1D* hA4_pt2s_PP = new TH1D("hA4_pt2s_PP","A4 vs p_{T}",1,ptbins2sLow);
	TH1D* hA4_pt3s_PP = new TH1D("hA4_pt3s_PP","A4 vs p_{T}",1,ptbins3sLow);
	
	TH1D* hA1_y1s_PP = new TH1D("hA1_y1s_PP","A1 vs y_{cm}",8,ybins1s);
	TH1D* hA1_y2s_PP = new TH1D("hA1_y2s_PP","A1 vs y_{cm}",4,ybins2s);
	TH1D* hA1_y3s_PP = new TH1D("hA1_y3s_PP","A1 vs y_{cm}",2,ybins3s);
	TH1D* hA2_y1s_PP = new TH1D("hA2_y1s_PP","A2 vs y_{cm}",8,ybins1s);
	TH1D* hA2_y2s_PP = new TH1D("hA2_y2s_PP","A2 vs y_{cm}",4,ybins2s);
	TH1D* hA2_y3s_PP = new TH1D("hA2_y3s_PP","A2 vs y_{cm}",2,ybins3s);
	TH1D* hA3_y1s_PP = new TH1D("hA3_y1s_PP","A3 vs y_{cm}",8,ybins1s);
	TH1D* hA3_y2s_PP = new TH1D("hA3_y2s_PP","A3 vs y_{cm}",4,ybins2s);
	TH1D* hA3_y3s_PP = new TH1D("hA3_y3s_PP","A3 vs y_{cm}",2,ybins3s);
	TH1D* hA4_y1s_PP = new TH1D("hA4_y1s_PP","A4 vs y_{cm}",8,ybins1s);
	TH1D* hA4_y2s_PP = new TH1D("hA4_y2s_PP","A4 vs y_{cm}",4,ybins2s);
	TH1D* hA4_y3s_PP = new TH1D("hA4_y3s_PP","A4 vs y_{cm}",2,ybins3s);
	
	/////Parameters of power law
	//PA
	TH1D* hAmp_pt1s_PA = new TH1D("hAmp_pt1s_PA","A vs p_{T}",3,ptbins1sHigh);
	TH1D* hAmp_pt2s_PA = new TH1D("hAmp_pt2s_PA","A vs p_{T}",2,ptbins2sHigh);
	TH1D* hAmp_pt3s_PA = new TH1D("hAmp_pt3s_PA","A vs p_{T}",1,ptbins3sHigh);
	TH1D* hM0_pt1s_PA = new TH1D("hM0_pt1s_PA","m_{0} vs p_{T}",3,ptbins1sHigh);
	TH1D* hM0_pt2s_PA = new TH1D("hM0_pt2s_PA","m_{0} vs p_{T}",2,ptbins2sHigh);
	TH1D* hM0_pt3s_PA = new TH1D("hM0_pt3s_PA","m_{0} vs p_{T}",1,ptbins3sHigh);
	TH1D* hPow_pt1s_PA = new TH1D("hPow_pt1s_PA","a vs p_{T}",3,ptbins1sHigh);
	TH1D* hPow_pt2s_PA = new TH1D("hPow_pt2s_PA","a vs p_{T}",2,ptbins2sHigh);
	TH1D* hPow_pt3s_PA = new TH1D("hPow_pt3s_PA","a vs p_{T}",1,ptbins3sHigh);
	TH1D* hMpow_pt1s_PA = new TH1D("hMpow_pt1s_PA","b vs p_{T}",3,ptbins1sHigh);
	TH1D* hMpow_pt2s_PA = new TH1D("hMpow_pt2s_PA","b vs p_{T}",2,ptbins2sHigh);
	TH1D* hMpow_pt3s_PA = new TH1D("hMpow_pt3s_PA","b vs p_{T}",1,ptbins3sHigh);
	
	TH1D* hAmp_y1s_PA = new TH1D("hAmp_y1s_PA","A vs y_{cm}",3,ybins1s);
	TH1D* hAmp_y2s_PA = new TH1D("hAmp_y2s_PA","A vs y_{cm}",2,ybins2s);
	TH1D* hAmp_y3s_PA = new TH1D("hAmp_y3s_PA","A vs y_{cm}",1,ybins3s);
	TH1D* hM0_y1s_PA = new TH1D("hM0_y1s_PA","m_{0} vs y_{cm}",3,ybins1s);
	TH1D* hM0_y2s_PA = new TH1D("hM0_y2s_PA","m_{0} vs y_{cm}",2,ybins2s);
	TH1D* hM0_y3s_PA = new TH1D("hM0_y3s_PA","m_{0} vs y_{cm}",1,ybins3s);
	TH1D* hPow_y1s_PA = new TH1D("hPow_y1s_PA","a vs y_{cm}",3,ybins1s);
	TH1D* hPow_y2s_PA = new TH1D("hPow_y2s_PA","a vs y_{cm}",2,ybins2s);
	TH1D* hPow_y3s_PA = new TH1D("hPow_y3s_PA","a vs y_{cm}",1,ybins3s);
	TH1D* hMpow_y1s_PA = new TH1D("hMpow_y1s_PA","b vs y_{cm}",3,ybins1s);
	TH1D* hMpow_y2s_PA = new TH1D("hMpow_y2s_PA","b vs y_{cm}",2,ybins2s);
	TH1D* hMpow_y3s_PA = new TH1D("hMpow_y3s_PA","b vs y_{cm}",1,ybins3s);
	
	//PP
	TH1D* hAmp_pt1s_PP = new TH1D("hAmp_pt1s_PP","A vs p_{T}",3,ptbins1sHigh);
	TH1D* hAmp_pt2s_PP = new TH1D("hAmp_pt2s_PP","A vs p_{T}",2,ptbins2sHigh);
	TH1D* hAmp_pt3s_PP = new TH1D("hAmp_pt3s_PP","A vs p_{T}",1,ptbins3sHigh);
	TH1D* hM0_pt1s_PP = new TH1D("hM0_pt1s_PP","m_{0} vs p_{T}",3,ptbins1sHigh);
	TH1D* hM0_pt2s_PP = new TH1D("hM0_pt2s_PP","m_{0} vs p_{T}",2,ptbins2sHigh);
	TH1D* hM0_pt3s_PP = new TH1D("hM0_pt3s_PP","m_{0} vs p_{T}",1,ptbins3sHigh);
	TH1D* hPow_pt1s_PP = new TH1D("hPow_pt1s_PP","a vs p_{T}",3,ptbins1sHigh);
	TH1D* hPow_pt2s_PP = new TH1D("hPow_pt2s_PP","a vs p_{T}",2,ptbins2sHigh);
	TH1D* hPow_pt3s_PP = new TH1D("hPow_pt3s_PP","a vs p_{T}",1,ptbins3sHigh);
	TH1D* hMpow_pt1s_PP = new TH1D("hMpow_pt1s_PP","b vs p_{T}",3,ptbins1sHigh);
	TH1D* hMpow_pt2s_PP = new TH1D("hMpow_pt2s_PP","b vs p_{T}",2,ptbins2sHigh);
	TH1D* hMpow_pt3s_PP = new TH1D("hMpow_pt3s_PP","b vs p_{T}",1,ptbins3sHigh);
	
	TH1D* hAmp_y1s_PP = new TH1D("hAmp_y1s_PP","A vs y_{cm}",3,ybins1s);
	TH1D* hAmp_y2s_PP = new TH1D("hAmp_y2s_PP","A vs y_{cm}",2,ybins2s);
	TH1D* hAmp_y3s_PP = new TH1D("hAmp_y3s_PP","A vs y_{cm}",1,ybins3s);
	TH1D* hM0_y1s_PP = new TH1D("hM0_y1s_PP","m_{0} vs y_{cm}",3,ybins1s);
	TH1D* hM0_y2s_PP = new TH1D("hM0_y2s_PP","m_{0} vs y_{cm}",2,ybins2s);
	TH1D* hM0_y3s_PP = new TH1D("hM0_y3s_PP","m_{0} vs y_{cm}",1,ybins3s);
	TH1D* hPow_y1s_PP = new TH1D("hPow_y1s_PP","a vs y_{cm}",3,ybins1s);
	TH1D* hPow_y2s_PP = new TH1D("hPow_y2s_PP","a vs y_{cm}",2,ybins2s);
	TH1D* hPow_y3s_PP = new TH1D("hPow_y3s_PP","a vs y_{cm}",1,ybins3s);
	TH1D* hMpow_y1s_PP = new TH1D("hMpow_y1s_PP","b vs y_{cm}",3,ybins1s);
	TH1D* hMpow_y2s_PP = new TH1D("hMpow_y2s_PP","b vs y_{cm}",2,ybins2s);
	TH1D* hMpow_y3s_PP = new TH1D("hMpow_y3s_PP","b vs y_{cm}",1,ybins3s);
	
	//////////Get parameters
	float ptLow;
	float ptHigh;
	float yLow;
	float yHigh;
	TString kineLabel;
	TFile* inFile;
	RooWorkspace* ws;
	
	//1s ptlow
	for (int i = 0; i < 3; i++)
	{
		ptLow = ptbins1sLow[i];
		ptHigh = ptbins1sLow[i+1];
		yLow = -1.93;
		yHigh = 1.93;
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_pt1s_PA->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_pt1s_PA->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_pt1s_PA->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_pt1s_PA->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_pt1s_PA->SetBinError(i+1,ws->var("A1")->getError());
		hA2_pt1s_PA->SetBinError(i+1,ws->var("A2")->getError());
		hA3_pt1s_PA->SetBinError(i+1,ws->var("A3")->getError());
		hA4_pt1s_PA->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
		yLow = 0;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_pt1s_PP->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_pt1s_PP->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_pt1s_PP->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_pt1s_PP->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_pt1s_PP->SetBinError(i+1,ws->var("A1")->getError());
		hA2_pt1s_PP->SetBinError(i+1,ws->var("A2")->getError());
		hA3_pt1s_PP->SetBinError(i+1,ws->var("A3")->getError());
		hA4_pt1s_PP->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
	}
	//2s ptlow
	ptLow = ptbins2sLow[0];
	ptHigh = ptbins2sLow[1];
	yLow = -1.93;
	yHigh = 1.93;
	kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hA1_pt2s_PA->SetBinContent(1,ws->var("A1")->getVal());
	hA2_pt2s_PA->SetBinContent(1,ws->var("A2")->getVal());
	hA3_pt2s_PA->SetBinContent(1,ws->var("A3")->getVal());
	hA4_pt2s_PA->SetBinContent(1,ws->var("A4")->getVal());
	hA1_pt2s_PA->SetBinError(1,ws->var("A1")->getError());
	hA2_pt2s_PA->SetBinError(1,ws->var("A2")->getError());
	hA3_pt2s_PA->SetBinError(1,ws->var("A3")->getError());
	hA4_pt2s_PA->SetBinError(1,ws->var("A4")->getError());
	delete ws;
	delete inFile;
	yLow = 0;
	kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hA1_pt2s_PP->SetBinContent(1,ws->var("A1")->getVal());
	hA2_pt2s_PP->SetBinContent(1,ws->var("A2")->getVal());
	hA3_pt2s_PP->SetBinContent(1,ws->var("A3")->getVal());
	hA4_pt2s_PP->SetBinContent(1,ws->var("A4")->getVal());
	hA1_pt2s_PP->SetBinError(1,ws->var("A1")->getError());
	hA2_pt2s_PP->SetBinError(1,ws->var("A2")->getError());
	hA3_pt2s_PP->SetBinError(1,ws->var("A3")->getError());
	hA4_pt2s_PP->SetBinError(1,ws->var("A4")->getError());
	delete ws;
	delete inFile;
	//3s ptlow
	ptLow = ptbins3sLow[0];
	ptHigh = ptbins3sLow[1];
	yLow = -1.93;
	yHigh = 1.93;
	kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hA1_pt3s_PA->SetBinContent(1,ws->var("A1")->getVal());
	hA2_pt3s_PA->SetBinContent(1,ws->var("A2")->getVal());
	hA3_pt3s_PA->SetBinContent(1,ws->var("A3")->getVal());
	hA4_pt3s_PA->SetBinContent(1,ws->var("A4")->getVal());
	hA1_pt3s_PA->SetBinError(1,ws->var("A1")->getError());
	hA2_pt3s_PA->SetBinError(1,ws->var("A2")->getError());
	hA3_pt3s_PA->SetBinError(1,ws->var("A3")->getError());
	hA4_pt3s_PA->SetBinError(1,ws->var("A4")->getError());
	delete ws;
	delete inFile;
	yLow = 0;
	kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hA1_pt3s_PP->SetBinContent(1,ws->var("A1")->getVal());
	hA2_pt3s_PP->SetBinContent(1,ws->var("A2")->getVal());
	hA3_pt3s_PP->SetBinContent(1,ws->var("A3")->getVal());
	hA4_pt3s_PP->SetBinContent(1,ws->var("A4")->getVal());
	hA1_pt3s_PP->SetBinError(1,ws->var("A1")->getError());
	hA2_pt3s_PP->SetBinError(1,ws->var("A2")->getError());
	hA3_pt3s_PP->SetBinError(1,ws->var("A3")->getError());
	hA4_pt3s_PP->SetBinError(1,ws->var("A4")->getError());
	delete ws;
	delete inFile;
	//1s y
	for (int i = 0; i < 8; i++)
	{
		ptLow = 0;
		ptHigh = 30;
		yLow = ybins1s[i];
		yHigh = ybins1s[i+1];
		float yLowPP = yLow;
		float yHighPP = yHigh;
		if (yLow < 0.0 && yHigh > 0.0)
			yLowPP = 0.0;
		else if (yLow < 0.0)
		{
			yLowPP = TMath::Abs(yHigh);
			yHighPP = TMath::Abs(yLow);
		}
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y1s_PA->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y1s_PA->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y1s_PA->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y1s_PA->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y1s_PA->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y1s_PA->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y1s_PA->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y1s_PA->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y1s_PP->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y1s_PP->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y1s_PP->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y1s_PP->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y1s_PP->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y1s_PP->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y1s_PP->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y1s_PP->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
	}
	//2s y
	for (int i = 0; i < 4; i++)
	{
		ptLow = 0;
		ptHigh = 30;
		yLow = ybins2s[i];
		yHigh = ybins2s[i+1];
		float yLowPP = yLow;
		float yHighPP = yHigh;
		if (yLow < 0.0 && yHigh > 0.0)
			yLowPP = 0.0;
		else if (yLow < 0.0)
		{
			yLowPP = TMath::Abs(yHigh);
			yHighPP = TMath::Abs(yLow);
		}
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y2s_PA->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y2s_PA->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y2s_PA->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y2s_PA->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y2s_PA->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y2s_PA->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y2s_PA->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y2s_PA->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y2s_PP->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y2s_PP->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y2s_PP->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y2s_PP->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y2s_PP->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y2s_PP->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y2s_PP->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y2s_PP->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
	}
	//3s y
	for (int i = 0; i < 2; i++)
	{
		ptLow = 0;
		ptHigh = 30;
		yLow = ybins3s[i];
		yHigh = ybins3s[i+1];
		float yLowPP = yLow;
		float yHighPP = yHigh;
		if (yLow < 0.0 && yHigh > 0.0)
			yLowPP = 0.0;
		else if (yLow < 0.0)
		{
			yLowPP = TMath::Abs(yHigh);
			yHighPP = TMath::Abs(yLow);
		}
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y2s_PA->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y2s_PA->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y2s_PA->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y2s_PA->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y2s_PA->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y2s_PA->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y2s_PA->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y2s_PA->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hA1_y2s_PP->SetBinContent(i+1,ws->var("A1")->getVal());
		hA2_y2s_PP->SetBinContent(i+1,ws->var("A2")->getVal());
		hA3_y2s_PP->SetBinContent(i+1,ws->var("A3")->getVal());
		hA4_y2s_PP->SetBinContent(i+1,ws->var("A4")->getVal());
		hA1_y2s_PP->SetBinError(i+1,ws->var("A1")->getError());
		hA2_y2s_PP->SetBinError(i+1,ws->var("A2")->getError());
		hA3_y2s_PP->SetBinError(i+1,ws->var("A3")->getError());
		hA4_y2s_PP->SetBinError(i+1,ws->var("A4")->getError());
		delete ws;
		delete inFile;
	}
	//1s pthigh
	for (int i = 0; i < 3; i++)
	{
		ptLow = ptbins1sHigh[i];
		ptHigh = ptbins1sHigh[i+1];
		yLow = -1.93;
		yHigh = 1.93;
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hAmp_pt1s_PA->SetBinContent(i+1,ws->var("amp")->getVal());
		hM0_pt1s_PA->SetBinContent(i+1,ws->var("m0")->getVal());
		hPow_pt1s_PA->SetBinContent(i+1,ws->var("pow")->getVal());
		hMpow_pt1s_PA->SetBinContent(i+1,ws->var("mpow")->getVal());
		hAmp_pt1s_PA->SetBinError(i+1,ws->var("amp")->getError());
		hM0_pt1s_PA->SetBinError(i+1,ws->var("m0")->getError());
		hPow_pt1s_PA->SetBinError(i+1,ws->var("pow")->getError());
		hMpow_pt1s_PA->SetBinError(i+1,ws->var("mpow")->getError());
		delete ws;
		delete inFile;
		yLow = 0;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hAmp_pt1s_PP->SetBinContent(i+1,ws->var("amp")->getVal());
		hM0_pt1s_PP->SetBinContent(i+1,ws->var("m0")->getVal());
		hPow_pt1s_PP->SetBinContent(i+1,ws->var("pow")->getVal());
		hMpow_pt1s_PP->SetBinContent(i+1,ws->var("mpow")->getVal());
		hAmp_pt1s_PP->SetBinError(i+1,ws->var("amp")->getError());
		hM0_pt1s_PP->SetBinError(i+1,ws->var("m0")->getError());
		hPow_pt1s_PP->SetBinError(i+1,ws->var("pow")->getError());
		hMpow_pt1s_PP->SetBinError(i+1,ws->var("mpow")->getError());
		delete ws;
		delete inFile;
	}
	//2s pthigh
	for (int i = 0; i < 2; i++)
	{
		ptLow = ptbins2sHigh[i];
		ptHigh = ptbins2sHigh[i+1];
		yLow = -1.93;
		yHigh = 1.93;
		
		kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hAmp_pt2s_PA->SetBinContent(i+1,ws->var("amp")->getVal());
		hM0_pt2s_PA->SetBinContent(i+1,ws->var("m0")->getVal());
		hPow_pt2s_PA->SetBinContent(i+1,ws->var("pow")->getVal());
		hMpow_pt2s_PA->SetBinContent(i+1,ws->var("mpow")->getVal());
		hAmp_pt2s_PA->SetBinError(i+1,ws->var("amp")->getError());
		hM0_pt2s_PA->SetBinError(i+1,ws->var("m0")->getError());
		hPow_pt2s_PA->SetBinError(i+1,ws->var("pow")->getError());
		hMpow_pt2s_PA->SetBinError(i+1,ws->var("mpow")->getError());
		delete ws;
		delete inFile;
		yLow = 0;
		kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
		inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
		ws = (RooWorkspace*)inFile->Get("workspace");
		hAmp_pt2s_PP->SetBinContent(i+1,ws->var("amp")->getVal());
		hM0_pt2s_PP->SetBinContent(i+1,ws->var("m0")->getVal());
		hPow_pt2s_PP->SetBinContent(i+1,ws->var("pow")->getVal());
		hMpow_pt2s_PP->SetBinContent(i+1,ws->var("mpow")->getVal());
		hAmp_pt2s_PP->SetBinError(i+1,ws->var("amp")->getError());
		hM0_pt2s_PP->SetBinError(i+1,ws->var("m0")->getError());
		hPow_pt2s_PP->SetBinError(i+1,ws->var("pow")->getError());
		hMpow_pt2s_PP->SetBinError(i+1,ws->var("mpow")->getError());
		delete ws;
		delete inFile;
	}
	//3s pthigh
	ptLow = ptbins3sHigh[0];
	ptHigh = ptbins3sHigh[1];
	yLow = -1.93;
	yHigh = 1.93;
	kineLabel = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hAmp_pt3s_PA->SetBinContent(1,ws->var("amp")->getVal());
	hM0_pt3s_PA->SetBinContent(1,ws->var("m0")->getVal());
	hPow_pt3s_PA->SetBinContent(1,ws->var("pow")->getVal());
	hMpow_pt3s_PA->SetBinContent(1,ws->var("mpow")->getVal());
	hAmp_pt3s_PA->SetBinError(1,ws->var("amp")->getError());
	hM0_pt3s_PA->SetBinError(1,ws->var("m0")->getError());
	hPow_pt3s_PA->SetBinError(1,ws->var("pow")->getError());
	hMpow_pt3s_PA->SetBinError(1,ws->var("mpow")->getError());
	delete ws;
	delete inFile;
	yLow = 0;
	kineLabel = getKineLabel (kPPDATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	inFile = new TFile(Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data()));
	ws = (RooWorkspace*)inFile->Get("workspace");
	hAmp_pt3s_PP->SetBinContent(1,ws->var("amp")->getVal());
	hM0_pt3s_PP->SetBinContent(1,ws->var("m0")->getVal());
	hPow_pt3s_PP->SetBinContent(1,ws->var("pow")->getVal());
	hMpow_pt3s_PP->SetBinContent(1,ws->var("mpow")->getVal());
	hAmp_pt3s_PP->SetBinError(1,ws->var("amp")->getError());
	hM0_pt3s_PP->SetBinError(1,ws->var("m0")->getError());
	hPow_pt3s_PP->SetBinError(1,ws->var("pow")->getError());
	hMpow_pt3s_PP->SetBinError(1,ws->var("mpow")->getError());
	delete ws;
	delete inFile;
	
	//////////Draw Plots
	TCanvas* c1 = new TCanvas("canv1","canv1",800,400);
	c1->Divide(2,1);
	TCanvas* c2 = new TCanvas("canv2","canv2",800,400);
	c2->Divide(2,1);
	TCanvas* c3 = new TCanvas("canv3","canv3",800,400);
	c3->Divide(2,1);
	TCanvas* c4 = new TCanvas("canv4","canv4",800,400);
	c4->Divide(2,1);
	TCanvas* c5 = new TCanvas("canv5","canv5",800,400);
	c5->Divide(2,1);
	TCanvas* c6 = new TCanvas("canv6","canv6",800,400);
	c6->Divide(2,1);
	
	gStyle->SetMarkerStyle(kFullCircle);
	//A1
	c1->cd(1);
	hA1_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hA1_pt1s_PA->SetYTitle("A1");
	hA1_pt1s_PA->SetStats(0);
	hA1_pt1s_PA->Draw("e1");
	hA1_pt2s_PA->Draw("SAME e1");
	hA1_pt3s_PA->Draw("SAME e1");
	hA1_pt1s_PP->SetLineColor(kGreen+1);
	hA1_pt2s_PP->SetLineColor(kGreen+1);
	hA1_pt3s_PP->SetLineColor(kGreen+1);
	hA1_pt1s_PP->Draw("SAME e1");
	hA1_pt2s_PP->Draw("SAME e1");
	hA1_pt3s_PP->Draw("SAME e1");
	c1->cd(2);
	hA1_y1s_PA->SetXTitle("y_{cm}");
	hA1_y1s_PA->SetYTitle("A1");
	hA1_y1s_PA->SetStats(0);
	hA1_y1s_PA->Draw("SAME e1");
	hA1_y2s_PA->Draw("SAME e1");
	hA1_y3s_PA->Draw("SAME e1");
	hA1_y1s_PP->SetLineColor(kGreen+1);
	hA1_y2s_PP->SetLineColor(kGreen+1);
	hA1_y3s_PP->SetLineColor(kGreen+1);
	hA1_y1s_PP->Draw("SAME e1");
	hA1_y2s_PP->Draw("SAME e1");
	hA1_y3s_PP->Draw("SAME e1");
	c1->SaveAs("ResultsBkg/ParamPlots_A1.pdf");
	
	//A2
	c2->cd(1);
	hA2_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hA2_pt1s_PA->SetYTitle("A2");
	hA2_pt1s_PA->SetStats(0);
	hA2_pt1s_PA->Draw("e1");
	hA2_pt2s_PA->Draw("SAME e1");
	hA2_pt3s_PA->Draw("SAME e1");
	hA2_pt1s_PP->SetLineColor(kGreen+1);
	hA2_pt2s_PP->SetLineColor(kGreen+1);
	hA2_pt3s_PP->SetLineColor(kGreen+1);
	hA2_pt1s_PP->Draw("SAME e1");
	hA2_pt2s_PP->Draw("SAME e1");
	hA2_pt3s_PP->Draw("SAME e1");
	c2->cd(2);
	hA2_y1s_PA->SetXTitle("y_{cm}");
	hA2_y1s_PA->SetYTitle("A2");
	hA2_y1s_PA->SetStats(0);
	hA2_y1s_PA->Draw("SAME e1");
	hA2_y2s_PA->Draw("SAME e1");
	hA2_y3s_PA->Draw("SAME e1");
	hA2_y1s_PP->SetLineColor(kGreen+1);
	hA2_y2s_PP->SetLineColor(kGreen+1);
	hA2_y3s_PP->SetLineColor(kGreen+1);
	hA2_y1s_PP->Draw("SAME e1");
	hA2_y2s_PP->Draw("SAME e1");
	hA2_y3s_PP->Draw("SAME e1");
	c2->SaveAs("ResultsBkg/ParamPlots_A2.pdf");
	
	//A3
	c3->cd(1);
	hA3_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hA3_pt1s_PA->SetYTitle("A3");
	hA3_pt1s_PA->SetStats(0);
	hA3_pt1s_PA->Draw("e1");
	hA3_pt2s_PA->Draw("SAME e1");
	hA3_pt3s_PA->Draw("SAME e1");
	hA3_pt1s_PP->SetLineColor(kGreen+1);
	hA3_pt2s_PP->SetLineColor(kGreen+1);
	hA3_pt3s_PP->SetLineColor(kGreen+1);
	hA3_pt1s_PP->Draw("SAME e1");
	hA3_pt2s_PP->Draw("SAME e1");
	hA3_pt3s_PP->Draw("SAME e1");
	c3->cd(2);
	hA3_y1s_PA->SetXTitle("y_{cm}");
	hA3_y1s_PA->SetYTitle("A3");
	hA3_y1s_PA->SetStats(0);
	hA3_y1s_PA->Draw("SAME e1");
	hA3_y2s_PA->Draw("SAME e1");
	hA3_y3s_PA->Draw("SAME e1");
	hA3_y1s_PP->SetLineColor(kGreen+1);
	hA3_y2s_PP->SetLineColor(kGreen+1);
	hA3_y3s_PP->SetLineColor(kGreen+1);
	hA3_y1s_PP->Draw("SAME e1");
	hA3_y2s_PP->Draw("SAME e1");
	hA3_y3s_PP->Draw("SAME e1");
	c3->SaveAs("ResultsBkg/ParamPlots_A3.pdf");
	
	//A4
	c4->cd(1);
	hA4_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hA4_pt1s_PA->SetYTitle("A4");
	hA4_pt1s_PA->SetStats(0);
	hA4_pt1s_PA->Draw("e1");
	hA4_pt2s_PA->Draw("SAME e1");
	hA4_pt3s_PA->Draw("SAME e1");
	hA4_pt1s_PP->SetLineColor(kGreen+1);
	hA4_pt2s_PP->SetLineColor(kGreen+1);
	hA4_pt3s_PP->SetLineColor(kGreen+1);
	hA4_pt1s_PP->Draw("SAME e1");
	hA4_pt2s_PP->Draw("SAME e1");
	hA4_pt3s_PP->Draw("SAME e1");
	c4->cd(2);
	hA4_y1s_PA->SetXTitle("y_{cm}");
	hA4_y1s_PA->SetYTitle("A4");
	hA4_y1s_PA->SetStats(0);
	hA4_y1s_PA->Draw("SAME e1");
	hA4_y2s_PA->Draw("SAME e1");
	hA4_y3s_PA->Draw("SAME e1");
	hA4_y1s_PP->SetLineColor(kGreen+1);
	hA4_y2s_PP->SetLineColor(kGreen+1);
	hA4_y3s_PP->SetLineColor(kGreen+1);
	hA4_y1s_PP->Draw("SAME e1");
	hA4_y2s_PP->Draw("SAME e1");
	hA4_y3s_PP->Draw("SAME e1");
	c4->SaveAs("ResultsBkg/ParamPlots_A4.pdf");
	
	//Amp
	c5->cd(1);
	hAmp_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hAmp_pt1s_PA->SetYTitle("Amp");
	hAmp_pt1s_PA->SetStats(0);
	hAmp_pt1s_PA->Draw("e1");
	hAmp_pt2s_PA->Draw("SAME e1");
	hAmp_pt3s_PA->Draw("SAME e1");
	hAmp_pt1s_PP->SetLineColor(kGreen+1);
	hAmp_pt2s_PP->SetLineColor(kGreen+1);
	hAmp_pt3s_PP->SetLineColor(kGreen+1);
	hAmp_pt1s_PP->Draw("SAME e1");
	hAmp_pt2s_PP->Draw("SAME e1");
	hAmp_pt3s_PP->Draw("SAME e1");
	c5->cd(2);
	hM0_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hM0_pt1s_PA->SetYTitle("m_{0}");
	hM0_pt1s_PA->SetStats(0);
	hM0_pt1s_PA->Draw();
	hM0_pt2s_PA->Draw("SAME e1");
	hM0_pt3s_PA->Draw("SAME e1");
	hM0_pt1s_PP->SetLineColor(kGreen+1);
	hM0_pt2s_PP->SetLineColor(kGreen+1);
	hM0_pt3s_PP->SetLineColor(kGreen+1);
	hM0_pt1s_PP->Draw("SAME e1");
	hM0_pt2s_PP->Draw("SAME e1");
	hM0_pt3s_PP->Draw("SAME e1");
	c5->SaveAs("ResultsBkg/ParamPlots_AmpM0.pdf");
	
	//Pow
	c6->cd(1);
	hPow_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hPow_pt1s_PA->SetYTitle("a");
	hPow_pt1s_PA->SetStats(0);
	hPow_pt1s_PA->Draw("e1");
	hPow_pt2s_PA->Draw("SAME e1");
	hPow_pt3s_PA->Draw("SAME e1");
	hPow_pt1s_PP->SetLineColor(kGreen+1);
	hPow_pt2s_PP->SetLineColor(kGreen+1);
	hPow_pt3s_PP->SetLineColor(kGreen+1);
	hPow_pt1s_PP->Draw("SAME e1");
	hPow_pt2s_PP->Draw("SAME e1");
	hPow_pt3s_PP->Draw("SAME e1");
	c6->cd(2);
	hMpow_pt1s_PA->SetXTitle("p_{T} (GeV/c)");
	hMpow_pt1s_PA->SetYTitle("b");
	hMpow_pt1s_PA->SetStats(0);
	hMpow_pt1s_PA->Draw();
	hMpow_pt2s_PA->Draw("SAME e1");
	hMpow_pt3s_PA->Draw("SAME e1");
	hMpow_pt1s_PP->SetLineColor(kGreen+1);
	hMpow_pt2s_PP->SetLineColor(kGreen+1);
	hMpow_pt3s_PP->SetLineColor(kGreen+1);
	hMpow_pt1s_PP->Draw("SAME e1");
	hMpow_pt2s_PP->Draw("SAME e1");
	hMpow_pt3s_PP->Draw("SAME e1");
	c6->SaveAs("ResultsBkg/ParamPlots_PowMpow.pdf");
}