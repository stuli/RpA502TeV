#include "FitData.C"
#include "RooMyPdf.h"
#include "RooMyPdfPP.h"

using std::vector;

void RunPseudoExpts(
		float ptLow=0, float ptHigh=30, 
		float yLow=0.0, float yHigh=0.4,//Run 1 has p going in -z direction
		int whichModel=1,   // Nominal = 0. Alternative = 1. Chebychev = 2. Power Law = 3. This is the model chosen to compare with nominal.
		int numTrials = 1,
		double maxchisq = 10
			)
{
	int cLow=0;
	int cHigh=200;
    float muPtCut=4.0;
	float dphiEp2Low = 0;
	float dphiEp2High = 100;
	
	gROOT->ProcessLine(".L RooMyPdf.cxx+");
	gROOT->ProcessLine(".L RooMyPdfPP.cxx+");
	
	TString modelLabel = "";
	if (whichModel == 1)
		modelLabel = "alt";
	else if (whichModel == 2)
		modelLabel = "cheb";
	else if (whichModel == 3)
		modelLabel = "pow";
	
	//Set up histograms
	TH1F* histo1sPA = new TH1F("PA_1SDiff","1S %Diff in Yield for PA",100,-20,20);
	TH1F* histo2sPA = new TH1F("PA_2SDiff","2S %Diff in Yield for PA",100,-40,40);
	TH1F* histo3sPA = new TH1F("PA_3SDiff","3S %Diff in Yield for PA",100,-100,100);
	TH1F* histo1sPP = new TH1F("PP_1SDiff","1S %Diff in Yield for PP",100,-20,20);
	TH1F* histo2sPP = new TH1F("PP_2SDiff","2S %Diff in Yield for PP",100,-40,40);
	TH1F* histo3sPP = new TH1F("PP_3SDiff","3S %Diff in Yield for PP",100,-100,100);
	TH1F* histo1sRpA = new TH1F("RpA_1SDiff","1S %Diff in RpA",100,-20,20);
	TH1F* histo2sRpA = new TH1F("RpA_2SDiff","2S %Diff in RpA",100,-40,40);
	TH1F* histo3sRpA = new TH1F("RpA_3SDiff","3S %Diff in RpA",100,-100,100);
	TCanvas* cPA =  new TCanvas("canvasPA","PA results",4,45,1100,400);
	TCanvas* cPP =  new TCanvas("canvasPP","PP results",4,45,1100,400);
	TCanvas* cRpA =  new TCanvas("canvasRpA","RpA results",4,45,1100,400);
	cPA->Divide(3,1);
	cPP->Divide(3,1);
	cRpA->Divide(3,1);

	cPA->cd(1);
	histo1sPA->SetXTitle("%Diff");
	histo1sPA->GetXaxis()->SetTitleSize(0.05);
	histo1sPA->GetYaxis()->SetLabelSize(0.05);
	histo1sPA->GetXaxis()->SetLabelSize(0.05);
	histo1sPA->SetStats(kTRUE);
	histo1sPA->Draw();
	cPA->cd(2);
	histo2sPA->SetXTitle("%Diff");
	histo2sPA->GetXaxis()->SetTitleSize(0.05);
	histo2sPA->GetYaxis()->SetLabelSize(0.05);
	histo2sPA->GetXaxis()->SetLabelSize(0.05);
	histo2sPA->Draw();
	cPA->cd(3);
	histo3sPA->SetXTitle("%Diff");
	histo3sPA->GetXaxis()->SetTitleSize(0.05);
	histo3sPA->GetYaxis()->SetLabelSize(0.05);
	histo3sPA->GetXaxis()->SetLabelSize(0.05);
	histo3sPA->Draw();
	
	cPP->cd(1);
	histo1sPP->SetXTitle("%Diff");
	histo1sPP->GetXaxis()->SetTitleSize(0.05);
	histo1sPP->GetYaxis()->SetLabelSize(0.05);
	histo1sPP->GetXaxis()->SetLabelSize(0.05);
	histo1sPP->SetStats(kTRUE);
	histo1sPP->Draw();
	cPP->cd(2);
	histo2sPP->SetXTitle("%Diff");
	histo2sPP->GetXaxis()->SetTitleSize(0.05);
	histo2sPP->GetYaxis()->SetLabelSize(0.05);
	histo2sPP->GetXaxis()->SetLabelSize(0.05);
	histo2sPP->Draw();
	cPP->cd(3);
	histo3sPP->SetXTitle("%Diff");
	histo3sPP->GetXaxis()->SetTitleSize(0.05);
	histo3sPP->GetYaxis()->SetLabelSize(0.05);
	histo3sPP->GetXaxis()->SetLabelSize(0.05);
	histo3sPP->Draw();
	
	cRpA->cd(1);
	histo1sRpA->SetXTitle("%Diff");
	histo1sRpA->GetXaxis()->SetTitleSize(0.05);
	histo1sRpA->GetYaxis()->SetLabelSize(0.05);
	histo1sRpA->GetXaxis()->SetLabelSize(0.05);
	histo1sRpA->SetStats(kTRUE);
	histo1sRpA->Draw();
	cRpA->cd(2);
	histo2sRpA->SetXTitle("%Diff");
	histo2sRpA->GetXaxis()->SetTitleSize(0.05);
	histo2sRpA->GetYaxis()->SetLabelSize(0.05);
	histo2sRpA->GetXaxis()->SetLabelSize(0.05);
	histo2sRpA->Draw();
	cRpA->cd(3);
	histo3sRpA->SetXTitle("%Diff");
	histo3sRpA->GetXaxis()->SetTitleSize(0.05);
	histo3sRpA->GetYaxis()->SetLabelSize(0.05);
	histo3sRpA->GetXaxis()->SetLabelSize(0.05);
	histo3sRpA->Draw();
	
	TNtuple* ntupleSig = new TNtuple("ntupleSig","Yields from fits","nSig1sPAalt:nSig2sPAalt:nSig3sPAalt:nSig1sPAnom:nSig2sPAnom:nSig3sPAnom:nSig1sPPalt:nSig2sPPalt:nSig3sPPalt:nSig1sPPnom:nSig2sPPnom:nSig3sPPnom",numTrials);
	TNtuple* ntupleDiff = new TNtuple("ntupleDiff","Differences from fits","diff1sPA:diff2sPA:diff3sPA:diff1sPP:diff2sPP:diff3sPP:diff1sRpA:diff2sRpA:diff3sRpA",numTrials);
	TNtuple* ntupleChisq = new TNtuple("ntupleChisq","Chisq/ndf from fits","chisqndfPAalt:chisqndfPAnom:chisqndfPPalt:chisqndfPPnom",numTrials);
	
	//Set up nominal files, workspaces, and generators
	TString kineLabelPA = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	TString NomFileNamePA = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabelPA.Data());
	cout << NomFileNamePA << endl;
	TFile* NomFilePA = TFile::Open(NomFileNamePA,"READ");
	RooWorkspace *NomwsPA = (RooWorkspace*)NomFilePA->Get("workspace");
	
	float yLowPP = yLow;
	float yHighPP = yHigh;
	if (yLow < 0.0 && yHigh > 0.0)
		yLowPP = 0.0;
	else if (yLow < 0.0)
	{
		yLowPP = TMath::Abs(yHigh);
		yHighPP = TMath::Abs(yLow);
	}
	TString kineLabelPP = getKineLabel (kPPDATA, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	TString NomFileNamePP = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabelPP.Data());
	cout << NomFileNamePP << endl;
	TFile* NomFilePP = TFile::Open(NomFileNamePP,"READ");
	RooWorkspace *NomwsPP = (RooWorkspace*)NomFilePP->Get("workspace");

	RooAbsPdf* genModelPA = NomwsPA->pdf("model");
	RooWorkspace *wsgenPA = new RooWorkspace("workspace");
	wsgenPA->import(*genModelPA);
	
	RooAbsPdf* genModelPP = NomwsPP->pdf("model");
	RooWorkspace *wsgenPP = new RooWorkspace("workspace");
	wsgenPP->import(*genModelPP);
	
	int nFailedTrials = 0;
	
	cout << "RUNNING PSEUDOEXPERIMENTS" << endl;
	for (int i = 0; i < numTrials; i++)
	{
		//RooWorkspace* mws;
		vector<double>* resultVector = new vector<double>(4);
		//Generate PA pseudodata
		RooDataSet* pseudoData = genModelPA->generate(*(wsgenPA->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//PA Alternate fit
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,resultVector,pseudoData);
		float nSig1sPAalt = (*resultVector)[1];
		float nSig2sPAalt = (*resultVector)[2];
		float nSig3sPAalt = (*resultVector)[3];
		float chisqndfPAalt = (*resultVector)[0];
		if (chisqndfPAalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//PA Nominal fit
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,0,resultVector,pseudoData);
		float nSig1sPAnom = (*resultVector)[1];
		float nSig2sPAnom = (*resultVector)[2];
		float nSig3sPAnom = (*resultVector)[3];
		float chisqndfPAnom = (*resultVector)[0];
		if (chisqndfPAnom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//Generate PP pseudodata
		delete pseudoData;
		pseudoData = genModelPP->generate(*(wsgenPP->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//PP Alternate fit
		FitData(kPPDATA,ptLow,ptHigh,yLowPP,yHighPP,cLow,cHigh,muPtCut,whichModel,resultVector,pseudoData);
		float nSig1sPPalt = (*resultVector)[1];
		float nSig2sPPalt = (*resultVector)[2];
		float nSig3sPPalt = (*resultVector)[3];
		float chisqndfPPalt = (*resultVector)[0];
		if (chisqndfPPalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//PP Nominal fit
		FitData(kPPDATA,ptLow,ptHigh,yLowPP,yHighPP,cLow,cHigh,muPtCut,0,resultVector,pseudoData);
		float nSig1sPPnom = (*resultVector)[1];
		float nSig2sPPnom = (*resultVector)[2];
		float nSig3sPPnom = (*resultVector)[3];
		float chisqndfPPnom = (*resultVector)[0];
		if (chisqndfPPnom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		/*
		//PA Alternate fit
		mws = new RooWorkspace("masterWorkspace");
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,mws,pseudoData);
		float nSig1sPAalt = mws->var("nSig1s")->getVal();
		float nSig2sPAalt = mws->var("nSig2s")->getVal();
		float nSig3sPAalt = mws->var("nSig3s")->getVal();
		float chisqndfPAalt = mws->var("chisqndf")->getVal();
		delete mws;
		if (chisqndfPAalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//PA Nominal fit
		mws = new RooWorkspace("masterWorkspace");
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,0,mws,pseudoData);
		float nSig1sPAnom = mws->var("nSig1s")->getVal();
		float nSig2sPAnom = mws->var("nSig2s")->getVal();
		float nSig3sPAnom = mws->var("nSig3s")->getVal();
		float chisqndfPAnom = mws->var("chisqndf")->getVal();
		delete mws;
		if (chisqndfPAnom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//Generate PP pseudodata
		delete pseudoData;
		pseudoData = genModelPP->generate(*(wsgenPP->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//PP Alternate fit
		mws = new RooWorkspace("masterWorkspace");
		FitData(kPPDATA,ptLow,ptHigh,yLowPP,yHighPP,cLow,cHigh,muPtCut,whichModel,mws,pseudoData);
		float nSig1sPPalt = mws->var("nSig1s")->getVal();
		float nSig2sPPalt = mws->var("nSig2s")->getVal();
		float nSig3sPPalt = mws->var("nSig3s")->getVal();
		float chisqndfPPalt = mws->var("chisqndf")->getVal();
		delete mws;
		if (chisqndfPPalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}
		
		//PP Nominal fit
		mws = new RooWorkspace("masterWorkspace");
		FitData(kPPDATA,ptLow,ptHigh,yLowPP,yHighPP,cLow,cHigh,muPtCut,0,mws,pseudoData);
		float nSig1sPPnom = mws->var("nSig1s")->getVal();
		float nSig2sPPnom = mws->var("nSig2s")->getVal();
		float nSig3sPPnom = mws->var("nSig3s")->getVal();
		float chisqndfPPnom = mws->var("chisqndf")->getVal();
		delete mws;
		if (chisqndfPPnom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			continue;
		}*/
		
		cout << "FINISHED FITS" << endl;
		
		float RpA1sAlt = nSig1sPAalt/nSig1sPPalt;
		float RpA2sAlt = nSig2sPAalt/nSig2sPPalt;
		float RpA3sAlt = nSig3sPAalt/nSig3sPPalt;
		
		float RpA1sNom = nSig1sPAnom/nSig1sPPnom;
		float RpA2sNom = nSig2sPAnom/nSig2sPPnom;
		float RpA3sNom = nSig3sPAnom/nSig3sPPnom;
		
		float diff1sPA = 100*(nSig1sPAalt-nSig1sPAnom)/nSig1sPAnom;
		float diff2sPA = 100*(nSig2sPAalt-nSig2sPAnom)/nSig2sPAnom;
		float diff3sPA = 100*(nSig3sPAalt-nSig3sPAnom)/nSig3sPAnom;
		
		float diff1sPP = 100*(nSig1sPPalt-nSig1sPPnom)/nSig1sPPnom;
		float diff2sPP = 100*(nSig2sPPalt-nSig2sPPnom)/nSig2sPPnom;
		float diff3sPP = 100*(nSig3sPPalt-nSig3sPPnom)/nSig3sPPnom;
		
		float diff1sRpA = 100*(RpA1sAlt-RpA1sNom)/RpA1sNom;
		float diff2sRpA = 100*(RpA2sAlt-RpA2sNom)/RpA2sNom;
		float diff3sRpA = 100*(RpA3sAlt-RpA3sNom)/RpA3sNom;
		
		cout << "CALCULATED VALUES" << endl;
		
		histo1sPA->Fill(diff1sPA);
		histo2sPA->Fill(diff2sPA);
		histo3sPA->Fill(diff3sPA);
		histo1sPP->Fill(diff1sPP);
		histo2sPP->Fill(diff2sPP);
		histo3sPP->Fill(diff3sPP);
		histo1sRpA->Fill(diff1sRpA);
		histo2sRpA->Fill(diff2sRpA);
		histo3sRpA->Fill(diff3sRpA);
		
		cout << "FILLED HISTOGRAMS" << endl;
		
		ntupleSig->Fill(nSig1sPAalt,nSig2sPAalt,nSig3sPAalt,nSig1sPAnom,nSig2sPAnom,nSig3sPAnom,nSig1sPPalt,nSig2sPPalt,nSig3sPPalt,nSig1sPPnom,nSig2sPPnom,nSig3sPPnom);
		ntupleDiff->Fill(diff1sPA,diff2sPA,diff3sPA,diff1sPP,diff2sPP,diff3sPP,diff1sRpA,diff2sRpA,diff3sRpA);
		ntupleChisq->Fill(chisqndfPAalt,chisqndfPAnom,chisqndfPPalt,chisqndfPPnom);
		
		cout << "FILLED TUPLES" << endl;
		
		cPA->Update();
		cPA->cd(1);
		histo1sPA->Draw();
		cPA->cd(2);
		histo2sPA->Draw();
		cPA->cd(3);
		histo3sPA->Draw();
		
		cPP->Update();
		cPP->cd(1);
		histo1sPP->Draw();
		cPP->cd(2);
		histo2sPP->Draw();
		cPP->cd(3);
		histo3sPP->Draw();
		
		cRpA->Update();
		cRpA->cd(1);
		histo1sRpA->Draw();
		cRpA->cd(2);
		histo2sRpA->Draw();
		cRpA->cd(3);
		histo3sRpA->Draw();
		
		cout << "DREW HISTOGRAMS" << endl << endl;
		
		cout << "RESULTS OF PSEUDOEXPERIMENT " << i << endl;
		cout << "PA ALT 1S = " << nSig1sPAalt << endl;
		cout << "PA ALT 2S = " << nSig2sPAalt << endl;
		cout << "PA ALT 3S = " << nSig3sPAalt << endl;
		cout << "PA Alt CHISQ/NDF = " << chisqndfPAalt << endl;
		cout << "PA NOM 1S = " << nSig1sPAnom << endl;
		cout << "PA NOM 2S = " << nSig2sPAnom << endl;
		cout << "PA NOM 3S = " << nSig3sPAnom << endl;
		cout << "PA NOM CHISQ/NDF = " << chisqndfPAnom << endl;
		cout << "PP ALT 1S = " << nSig1sPPalt << endl;
		cout << "PP ALT 2S = " << nSig2sPPalt << endl;
		cout << "PP ALT 3S = " << nSig3sPPalt << endl;
		cout << "PP Alt CHISQ/NDF = " << chisqndfPPalt << endl;
		cout << "PP NOM 1S = " << nSig1sPPnom << endl;
		cout << "PP NOM 2S = " << nSig2sPPnom << endl;
		cout << "PP NOM 3S = " << nSig3sPPnom << endl;
		cout << "PP NOM CHISQ/NDF = " << chisqndfPPnom << endl;
		
		delete pseudoData;
	}
	
	
	
	//Write results to file
	TFile* outfile = new TFile(Form("ResultsBkg/PseudoExpResults_") + modelLabel + Form("_pt%.1f-%.1f_y%.2f-%.2f.root",ptLow,ptHigh,yLow,yHigh),"RECREATE");
	histo1sPA->Write();
	histo2sPA->Write();
	histo3sPA->Write();
	histo1sPP->Write();
	histo2sPP->Write();
	histo3sPP->Write();
	histo1sRpA->Write();
	histo2sRpA->Write();
	histo3sRpA->Write();
	ntupleSig->Write();
	ntupleDiff->Write();
	ntupleChisq->Write();
	
	outfile->Close();
}