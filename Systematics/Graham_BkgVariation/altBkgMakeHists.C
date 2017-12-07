#include "../../HeaderFiles/cutsAndBin.h"

using namespace RooFit;

void altBkgMakeHists(int collId = kPADATA)
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
	TString fileNamePA1 = "../../skimmedFiles_Graham/yskimPA1st_SSign_201711181848_unIdentified.root";
	TString fileNamePA2 = "../../skimmedFiles_Graham/yskimPA2nd_SSign_201711182048_unIdentified.root";
	TString fileNamePP = "../../skimmedFiles_Graham/yskimPP_L1DoubleMu0PD_Trig-L1DoubleMuOpen2016_SSign_20171118239_unIdentified.root";
	TString fileNameMC = "../../../upsilonFlatDimuonMass2BodyNtuple.root";
	
	//Open data files
	TFile* f1;
	TFile* f2;
	TString kineCut[nPtBins+1];
	if (collId == kPADATA)
	{
		f1 = new TFile(fileNamePA1,"READ");
		f2 = new TFile(fileNamePA2,"READ");
		kineCut[0] = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptBinDelim[0], ptBinDelim[nPtBins], yLow+yBoost, yHigh+yBoost, eta_high,eta_low, eta_high,eta_low );
		for (int i = 1; i<nPtBins+1; i++)
			kineCut[i] = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptBinDelim[i-1], ptBinDelim[i], yLow+yBoost, yHigh+yBoost, eta_high,eta_low, eta_high,eta_low );
	}
	else if (collId == kPPDATA)
	{
		f1 = new TFile(fileNamePP,"READ");
		kineCut[0] = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptBinDelim[0], ptBinDelim[nPtBins], 0.0, yHigh, eta_high,eta_low, eta_high,eta_low );
		for (int i = 1; i<nPtBins+1; i++)
			kineCut[i] = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptBinDelim[i-1], ptBinDelim[i], 0.0, yHigh, eta_high,eta_low, eta_high,eta_low );
	}
	if (muPtCut>0)
		for (int i = 0; i< nPtBins+1; i++)
			kineCut[i] = kineCut[i] + Form(" && (pt1>%.2f) && (pt2>%.2f) ", (float)muPtCut, (float)muPtCut);
	
	RooWorkspace* ws = new RooWorkspace("workspace");
	
	//Load dataset
	RooDataSet* dataset = (RooDataSet*)f1->Get("dataset");
	if (collId == kPADATA)
	{
		RooDataSet* dataset2 = (RooDataSet*)f2->Get("dataset");
		dataset->append(*dataset2);
	}
	ws->import(*dataset);
	cout << "Loaded dataset." << endl;
	//Reduced dataset with cuts
	RooDataSet* reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut[0].Data() );
	reducedDS->SetName("reducedDS_Int");
	ws->import(*reducedDS);
	
	//Load MC
	TFile* fMC = new TFile(fileNameMC,"READ");
	TNtuple* MCtuple = (TNtuple*)fMC->Get("upsilonNtuple");
	float UpsM, UpsRap, UpsPt, MuPPt, MuMPt, MuPy, MuMy;
	MCtuple->SetBranchAddress("UpsM",&UpsM);
	MCtuple->SetBranchAddress("UpsRap",&UpsRap);
	MCtuple->SetBranchAddress("UpsPt",&UpsPt);
	MCtuple->SetBranchAddress("MuPPt",&MuPPt);
	MCtuple->SetBranchAddress("MuMPt",&MuMPt);
	MCtuple->SetBranchAddress("MuPy",&MuPy);
	MCtuple->SetBranchAddress("MuMy",&MuMy);
	
	//Weighting functions for MC
	TF1* mfunc = new TF1("mfunc","[0]*exp(-x/[1])",6,20);
	mfunc->SetParameters(100,5.2);
	mfunc->SetParNames("A","#lambda");
	TF1* rfunc = new TF1("upsRapidityFunc","gaus(0)+gaus(3)",-5,5);
	rfunc->SetParameters(70,0.0,1.0,550,0.0,0.23);
	rfunc->SetParNames("A1","mean1","sigma1","A2","mean2","sigma2");
	TF1* ptfunc = new TF1("ptfunc","[0]*x/(exp(x/[1])+1)",0,5);
	ptfunc->SetParameters(100,2.71);
	ptfunc->SetParNames("A","T");
	
	//////////
	//////////PtEta hists
	//////////
	
	//Make PtEta hist from data
	float ptEta_ptLow = 4;
	float ptEta_ptHigh = 14;
	int ptEta_nEtaBins = 100;
	float ptEta_etaBinSize = (eta_high - eta_low) / ptEta_nEtaBins;
	int ptEta_nPtBins = 100;
	float ptEta_ptBinSize = (ptEta_ptHigh - ptEta_ptLow) / ptEta_nPtBins;
	
	TH2F* hPtEtaData = (TH2F*)dataset->createHistogram("hPtEtaData",*(ws->var("eta1")),RooFit::Binning(ptEta_nEtaBins,eta_low,eta_high),RooFit::YVar(*(ws->var("pt1")),RooFit::Binning(ptEta_nPtBins,ptEta_ptLow,ptEta_ptHigh)));
	hPtEtaData->SetName("hPtEtaData");
	hPtEtaData->SetTitle("Eta and p_{T} from Data");
	cout << "Made PtEta hist from data." << endl;
	
	//Make PtEta hist from MC
	TH2F* hPtEtaMC = new TH2F("hPtEtaMC","Rapidity and p_{T} from MC",ptEta_nEtaBins,eta_low,eta_high,ptEta_nPtBins,ptEta_ptLow,ptEta_ptHigh);
	
	int nEntriesMC = MCtuple->GetEntries();
	for (int i = 0; i < nEntriesMC; i++)
	{
		if (i % (nEntriesMC/10) == 0) cout << Form("Filling hPtEtaMC: %d/%d (%.0f%%)",i,nEntriesMC,(100.0*i)/nEntriesMC) << endl; //report progress at 10% intervals
		
		MCtuple->GetEntry(i);	
		
		//Apply weighting
		double mWeight = mfunc->Eval(UpsM);
		double rWeight = rfunc->Eval(UpsRap);
		double ptWeight = ptfunc->Eval(UpsPt);
		
		double weight = mWeight * rWeight * ptWeight;
		
		//Apply kinematic cuts
		if (MuPPt>muPtCut && MuMPt>muPtCut && MuPy<eta_high && MuPy>eta_low && MuMy<eta_high && MuMy>eta_low)
		{
			hPtEtaMC->Fill(MuPy,MuPPt,weight);
		}
	}
	cout << "Made PtEta hist from MC." << endl;
	
	//Take ratio of data to MC
	TH2F* hPtEtaRatio = (TH2F*)hPtEtaData->Clone("hPtEtaRatio");
	hPtEtaRatio->SetTitle("Ratio of data/MC");
	hPtEtaRatio->Divide(hPtEtaMC);
	cout << "Made PtEta ratio of data/MC." << endl;
	
	//Reweight
	TH2F* hPtEtaMCrw = new TH2F("hPtEtaMCrw","Reweighted Rapidity and p_{T} from MC",ptEta_nEtaBins,eta_low,eta_high,ptEta_nPtBins,ptEta_ptLow,ptEta_ptHigh);
	
	//Reminder of MC variables: UpsM, UpsRap, UpsPt, MuPPt, MuMPt, MuPy, MuMy;
	RooRealVar* UpsMVar = new RooRealVar("UpsM","MC mass variable",0,200,"GeV/c^{2}");
	RooRealVar* UpsRapVar = new RooRealVar("UpsRap","MC dimuon rapidity",-5,5,"");
	RooRealVar* UpsPtVar = new RooRealVar("UpsPt","MC dimuon pT",0,100,"GeV/c");
	RooRealVar* MuPPtVar = new RooRealVar("MuPPt","MC MuPlus pT",0,100,"GeV/c");
	RooRealVar* MuMPtVar = new RooRealVar("MuMPt","MC MuMinus pT",0,100,"GeV/c");
	RooRealVar* MuPyVar = new RooRealVar("MuPy","MC MuPlus rapidity",-5,5,"");
	RooRealVar* MuMyVar = new RooRealVar("MuMy","MC MuMinus rapidity",-5,5,"");
	RooRealVar* UpsWt = new RooRealVar("UpsWt","MC event weight",0,100,"");
	
	//RooArgSet* argSetMC = new RooArgSet(*UpsMVar, *UpsRapVar, *UpsPtVar, *MuPPtVar, *MuMPtVar, *MuPyVar, *MuMyVar, *UpsWt);
	RooArgSet* argSetMC = new RooArgSet(*(ws->var("mass")), *UpsRapVar, *UpsPtVar, *MuPPtVar, *MuMPtVar, *MuPyVar, *MuMyVar, *UpsWt);
	RooDataSet* reducedMC[nPtBins+1];
	reducedMC[0] = new RooDataSet("reducedMC_Int","reducedMC_Int",*argSetMC);
	ws->import(*reducedMC[0]);
	for (int i = 1; i < nPtBins+1; i++)
	{
		reducedMC[i] = new RooDataSet(Form("reducedMC_Bin%d",i),Form("reducedMC_Bin%d",i),*argSetMC);
		ws->import(*reducedMC[i]);
	}
	
	for (int i = 0; i < nEntriesMC; i++)
	{
		if (i % (nEntriesMC/10) == 0) cout << Form("Filling hPtEtaMCrw: %d/%d (%.0f%%)",i,nEntriesMC,(100.0*i)/nEntriesMC) << endl; //report progress at 10% intervals
		
		MCtuple->GetEntry(i);	
		
		//Apply weighting
		int rbin = TMath::Ceil((MuPy - eta_low)/ptEta_etaBinSize);
		int ptbin = TMath::Ceil((MuPPt - ptEta_ptLow)/ptEta_ptBinSize);
		double muplWeight = hPtEtaRatio->GetBinContent(rbin,ptbin);
		
		double mWeight = mfunc->Eval(UpsM);
		double rWeight = rfunc->Eval(UpsRap);
		double ptWeight = ptfunc->Eval(UpsPt);
		
		double weight = muplWeight * mWeight * rWeight * ptWeight;
		
		//Apply kinematic cuts
		if (MuPPt>muPtCut && MuMPt>muPtCut && MuPy<eta_high && MuPy>eta_low && MuMy<eta_high && MuMy>eta_low)
		{
			hPtEtaMCrw->Fill(MuPy,MuPPt,weight);
			
			ws->var("mass")->setVal(UpsM);
			//UpsMVar->setVal(UpsM);
			UpsRapVar->setVal(UpsRap);
			UpsPtVar->setVal(UpsPt);
			MuPPtVar->setVal(MuPPt);
			MuMPtVar->setVal(MuMPt);
			MuPyVar->setVal(MuPy);
			MuMyVar->setVal(MuMy);
			UpsWt->setVal(weight);
			
			if (UpsPt < ptBinDelim[nPtBins])
				reducedMC[0]->add(*argSetMC);
			
			for (int i = 1; i < nPtBins+1; i++)
				if (UpsPt > ptBinDelim[i-1] && UpsPt < ptBinDelim[i])
					reducedMC[i]->add(*argSetMC);
		}
	}
	cout << "Made reweighted PtEta hist." << endl;
	
	RooDataSet* reducedMCrw[nPtBins+1];
	for (int i = 0; i < nPtBins+1; i++)
	{
		reducedMCrw[i] = new RooDataSet(Form("%s_rw",reducedMC[i]->GetName()),Form("%s_rw",reducedMC[i]->GetTitle()),reducedMC[i],*reducedMC[i]->get(),0,"UpsWt");
		ws->import(*reducedMCrw[i]);
	}
	
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
	RooPlot* massPlotMC[nPtBins+1];
	for (int i = 1; i < nPtBins+1; i++)
	{
		cMass->cd(i+1);
		massPlotMC[i] = ws->var("mass")->frame(nMassBinMC);
		//massPlotMC[i] = ws->var("UpsM")->frame(nMassBinMC);
		reducedMCrw[i]->plotOn(massPlotMC[i],Name("MC Hist"));
		massPlotMC[i]->SetAxisRange(massLowMC,massHighMC,"X");
		massPlotMC[i]->Draw();
	}
	cMass->cd(1);
	RooPlot* massPlotData = ws->var("mass")->frame(nMassBin);
	reducedDS->plotOn(massPlotData,Name("Data Hist"));
	massPlotData->SetAxisRange(massLow,massHigh,"X");
	massPlotData->Draw();
	
	//Save
	TString outfileName;
	if (collId == kPADATA)
		outfileName = "altBkgData.root";
	else if (collId == kPPDATA)
		outfileName = "altBkgDataPP.root";
	TFile* outfile = new TFile(outfileName,"RECREATE");
	
	hPtEtaData->Write();
	hPtEtaMC->Write();
	hPtEtaRatio->Write();
	hPtEtaMCrw->Write();
	ws->Write();
}