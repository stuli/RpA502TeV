#include "../../HeaderFiles/cutsAndBin.h"

using namespace RooFit;

void altBkgFitHists(int collId = kPADATA)
{
	float dphiEp2Low = 0 ;
	float dphiEp2High = 100 ;
	float eta_low = -2.4;
	float eta_high = 2.4;
	float yLow = -1.93;
	float yHigh = 1.93;
	float yBoost = 0.47;
	float muPtCut=4.0;
	
	float massLow = 6;
	float massHigh = 20;
	float massLowForPlot = massLow;    
	float massHighForPlot = massHigh;
	float massLowMC = 6;
	float massHighMC = 20;
	int   nMassBin  = (massHigh-massLow)*100;
	int   nMassBinMC  = (massHighMC-massLowMC)*100;
	
	const int nPtBins = 4;
	float ptBinDelim[5] = {0.0,1.5,3.0,4.5,6.0};
	
	//File names
	TString fileName;
	if (collId == kPADATA)
		fileName = "altBkgData.root";
	else if (collId == kPPDATA)
		fileName = "altBkgDataPP.root";
	//TString fileNamePA = "altBkgData.root";
	//TString fileNamePP = "../../skimmedFiles_Graham/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMuOpen2016_SSign_20171118239_unIdentified.root";
	
	//Open data file
	TFile* file = new TFile(fileName,"READ");
	
	RooWorkspace* ws = (RooWorkspace*)file->Get("workspace");
	
	TH2F* hPtEtaData = (TH2F*)file->Get("hPtEtaData");
	TH2F* hPtEtaMC = (TH2F*)file->Get("hPtEtaMC");
	TH2F* hPtEtaRatio = (TH2F*)file->Get("hPtEtaRatio");
	TH2F* hPtEtaMCrw = (TH2F*)file->Get("hPtEtaMCrw");
	
	//////////
	//////////MAKE PLOTS
	//////////
	
	RooPlot* massPlotMC[nPtBins+1];
	for (int i = 1; i < nPtBins+1; i++)
	{
		massPlotMC[i] = ws->var("mass")->frame(nMassBinMC);
		//massPlotMC[i] = ws->var("UpsM")->frame(nMassBinMC);
		ws->data(Form("reducedMC_Bin%d_rw",i))->plotOn(massPlotMC[i],Name("MC Hist"));
		massPlotMC[i]->SetAxisRange(massLowMC,massHighMC,"X");
		massPlotMC[i]->SetTitle(Form("MC, %.1f < p_{T} < %.1f", ptBinDelim[i-1], ptBinDelim[i]));
	}
	RooPlot* massPlotData = ws->var("mass")->frame(nMassBin);
	ws->data("reducedDS_Int")->plotOn(massPlotData,Name("Data Hist"));
	massPlotData->SetAxisRange(massLow,massHigh,"X");
	massPlotData->SetTitle(Form("Data, %.1f < p_{T} < %.1f", ptBinDelim[0], ptBinDelim[nPtBins]));

	//////////
	//////////FITTING
	//////////
	
	RooWorkspace* outws = new RooWorkspace("altBkgWorkspace");
	
	TString strSumExpErf = "(@1*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1) + TMath::Max(0,1-TMath::Exp(-(@0-@4)/@5)))*TMath::Exp(-@0/@6)";
	
	RooRealVar* ratio[nPtBins+1];
	RooRealVar* erfMu[nPtBins+1];
	RooRealVar* erfSig[nPtBins+1];
	RooRealVar* expMu[nPtBins+1];
	RooRealVar* expLam[nPtBins+1];
	RooRealVar* decayLam[nPtBins+1];
	
	float expMuInit[nPtBins+1] = {8.0, 8.0, 8.0, 7.0, 6.0};
	float expMuMax[nPtBins+1] = {8.0, 8.0, 8.0, 7.0, 6.0};
	
	RooGenericPdf* bkgSumExpErf[nPtBins+1];
	for (int i = 1; i < nPtBins+1; i++)
	{
		ratio[i] = new RooRealVar(Form("R_%d",i),"Ratio", 0.8, 0.1, 10.0);
		erfMu[i] = new RooRealVar(Form("erfMu_%d",i),"Erf #mu", 8.0, 6.0, 10.0);
		erfSig[i] = new RooRealVar(Form("erfSig_%d",i),"Erf #sigma", 0.5, 0.0, 20.0);
		expMu[i] = new RooRealVar(Form("expMu_%d",i),"Exp #mu", expMuInit[i], 5.0, expMuMax[i]);
		expLam[i] = new RooRealVar(Form("expLam_%d",i),"Exp #lambda", 1.0, 0.0, 20.0);
		decayLam[i] = new RooRealVar(Form("decayLam_%d",i),"Decay #lambda", 4.4, 0.0, 100.0);
		bkgSumExpErf[i] = new RooGenericPdf(Form("bkgSumExpErf_%d",i),"Sum of Exp and Erf",strSumExpErf,RooArgList(*(ws->var("mass")),*ratio[i],*erfMu[i],*erfSig[i],*expMu[i],*expLam[i],*decayLam[i]));
		//bkgSumExpErf[i] = new RooGenericPdf(Form("bkgSumExpErf_%d",i),"Sum of Exp and Erf",strSumExpErf,RooArgList(*(ws->var("UpsM")),*ratio[i],*erfMu[i],*erfSig[i],*expMu[i],*expLam[i],*decayLam[i]));
		
		bkgSumExpErf[i]->fitTo(*(ws->data(Form("reducedMC_Bin%d_rw",i))),Save(), Hesse(kTRUE),Range(massLowMC, massHighMC),Timer(kTRUE));
		bkgSumExpErf[i]->plotOn(massPlotMC[i],Name(Form("model_%d",i)));
		bkgSumExpErf[i]->paramOn(massPlotMC[i]);
		
		ratio[i]->setConstant(kTRUE);
		erfMu[i]->setConstant(kTRUE);
		erfSig[i]->setConstant(kTRUE);
		expMu[i]->setConstant(kTRUE);
		expLam[i]->setConstant(kTRUE);
		decayLam[i]->setConstant(kTRUE);
		
		//outws->import(*bkgSumExpErf[i]);
	}
	
	RooRealVar* a1 = new RooRealVar("A1","A1",200,0,100000);
	RooRealVar* a2 = new RooRealVar("A2","A2",200,0,100000);
	RooRealVar* a3 = new RooRealVar("A3","A3",200,0,100000);
	RooRealVar* a4 = new RooRealVar("A4","A4",200,0,100000);
	
	RooAddPdf* bkgLinCom = new RooAddPdf("bkgLinCom","Linear Combination Background", RooArgList(*bkgSumExpErf[1],*bkgSumExpErf[2],*bkgSumExpErf[3],*bkgSumExpErf[4]), RooArgList(*a1,*a2,*a3,*a4));
	bkgLinCom->fitTo(*(ws->data("reducedDS_Int")),Save(), Hesse(kTRUE),Range(massLowMC, massHighMC),Timer(kTRUE));
	bkgLinCom->plotOn(massPlotData,Name("model"));
	bkgLinCom->paramOn(massPlotData);
	
	outws->import(*bkgLinCom);
	
	//////////
	//////////DRAWING
	//////////
	
	//Draw PtEta hists
	TCanvas* cPtEta = new TCanvas("cPtEta","Plots of single-muon pT vs Eta",800,800);
	cPtEta->Divide(2,2);
	
	cPtEta->cd(1);
	hPtEtaData->Draw("COL");
	cPtEta->cd(2);
	hPtEtaMC->Draw("COL");
	cPtEta->cd(3);
	hPtEtaRatio->Draw("COL");
	cPtEta->cd(4);
	hPtEtaMCrw->Draw("COL");

	//Draw datasets
	TCanvas* cMass = new TCanvas("cMass","Plots of same-sign dimuon mass",1200,800);
	cMass->Divide(2,3);
	for (int i = 1; i < nPtBins+1; i++)
	{
		cMass->cd(i+1);
		massPlotMC[i]->Draw();
	}
	cMass->cd(1);
	massPlotData->Draw();
	
	//////////OUTPUT
	TString outfileName;
	if (collId == kPADATA)
		outfileName = "altBkgModels.root";
	else if (collId == kPPDATA)
		outfileName = "altBkgModelsPP.root";
	
	TFile* outfile = new TFile(outfileName,"RECREATE");
	outws->Write();
}