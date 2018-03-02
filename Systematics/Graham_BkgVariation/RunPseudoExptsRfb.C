#include "FitData.C"
#include "RooMyPdf.h"
#include "RooMyPdfPP.h"

using std::vector;

void RunPseudoExptsRfb(
		float ptLow=0, float ptHigh=30, 
		float yLow=0.0, float yHigh=0.4,//Run 1 has p going in -z direction
		int whichModel=1,   // Nominal = 0. Alternative = 1. Chebychev = 2. Power Law = 3. This is the model chosen to compare with nominal.
		int numTrials = 1,
		double maxchisq = 10,
		double widthMult = 1
			)
{
	int cLow=0;
	int cHigh=200;
    float muPtCut=4.0;
	float dphiEp2Low = 0;
	float dphiEp2High = 100;
	
	//Set up nominal files, workspaces, and generators
	TString kineLabelPL = getKineLabel (kPADATA, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	TString NomFileNamePL = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabelPL.Data());
	cout << NomFileNamePL << endl;
	TFile* NomFilePL = TFile::Open(NomFileNamePL,"READ");
	RooWorkspace *NomwsPL = (RooWorkspace*)NomFilePL->Get("workspace");
	
	float yLowMI = -yHigh;
	float yHighMI = -yLow;
	if (yLow == 0.00)
		yHighMI = 0.00;
	TString kineLabelMI = getKineLabel (kPADATA, ptLow, ptHigh, yLowMI, yHighMI, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	TString NomFileNameMI = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabelMI.Data());
	cout << NomFileNameMI << endl;
	TFile* NomFileMI = TFile::Open(NomFileNameMI,"READ");
	RooWorkspace *NomwsMI = (RooWorkspace*)NomFileMI->Get("workspace");

	RooAbsPdf* genModelPL = NomwsPL->pdf("model");
	RooWorkspace *wsgenPL = new RooWorkspace("workspace");
	wsgenPL->import(*genModelPL);
	
	RooAbsPdf* genModelMI = NomwsMI->pdf("model");
	RooWorkspace *wsgenMI = new RooWorkspace("workspace");
	wsgenMI->import(*genModelMI);
	
	NomFilePL->Close();
	delete NomFilePL;
	NomFileMI->Close();
	delete NomFileMI;
	
	TFile* outfile = new TFile(Form("ResultsBkg/PseudoExpResults_pt%.1f-%.1f_y%.2f-%.2f_rfb.root",ptLow,ptHigh,yLow,yHigh),"RECREATE");
	
	//Set up histograms
	TH1F* histo1sPL = new TH1F("PL_1SDiff","1S %Diff in Yield for PL",100,-20*widthMult,20*widthMult);
	TH1F* histo2sPL = new TH1F("PL_2SDiff","2S %Diff in Yield for PL",100,-40*widthMult,40*widthMult);
	TH1F* histo3sPL = new TH1F("PL_3SDiff","3S %Diff in Yield for PL",100,-100*widthMult,100*widthMult);
	TH1F* histo1sMI = new TH1F("MI_1SDiff","1S %Diff in Yield for MI",100,-20*widthMult,20*widthMult);
	TH1F* histo2sMI = new TH1F("MI_2SDiff","2S %Diff in Yield for MI",100,-40*widthMult,40*widthMult);
	TH1F* histo3sMI = new TH1F("MI_3SDiff","3S %Diff in Yield for MI",100,-100*widthMult,100*widthMult);
	TH1F* histo1sRfb = new TH1F("RPL_1SDiff","1S %Diff in Rfb",100,-20*widthMult,20*widthMult);
	TH1F* histo2sRfb = new TH1F("Rfb_2SDiff","2S %Diff in Rfb",100,-40*widthMult,40*widthMult);
	TH1F* histo3sRfb = new TH1F("Rfb_3SDiff","3S %Diff in Rfb",100,-100*widthMult,100*widthMult);
	TCanvas* cPL =  new TCanvas("canvasPL","PL results",4,45,1100,400);
	TCanvas* cMI =  new TCanvas("canvasMI","MI results",4,45,1100,400);
	TCanvas* cRfb =  new TCanvas("canvasRfb","Rfb results",4,45,1100,400);
	cPL->Divide(3,1);
	cMI->Divide(3,1);
	cRfb->Divide(3,1);

	cPL->cd(1);
	histo1sPL->SetXTitle("%Diff");
	histo1sPL->GetXaxis()->SetTitleSize(0.05);
	histo1sPL->GetYaxis()->SetLabelSize(0.05);
	histo1sPL->GetXaxis()->SetLabelSize(0.05);
	histo1sPL->SetStats(kTRUE);
	histo1sPL->Draw();
	cPL->cd(2);
	histo2sPL->SetXTitle("%Diff");
	histo2sPL->GetXaxis()->SetTitleSize(0.05);
	histo2sPL->GetYaxis()->SetLabelSize(0.05);
	histo2sPL->GetXaxis()->SetLabelSize(0.05);
	histo2sPL->Draw();
	cPL->cd(3);
	histo3sPL->SetXTitle("%Diff");
	histo3sPL->GetXaxis()->SetTitleSize(0.05);
	histo3sPL->GetYaxis()->SetLabelSize(0.05);
	histo3sPL->GetXaxis()->SetLabelSize(0.05);
	histo3sPL->Draw();
	
	cMI->cd(1);
	histo1sMI->SetXTitle("%Diff");
	histo1sMI->GetXaxis()->SetTitleSize(0.05);
	histo1sMI->GetYaxis()->SetLabelSize(0.05);
	histo1sMI->GetXaxis()->SetLabelSize(0.05);
	histo1sMI->SetStats(kTRUE);
	histo1sMI->Draw();
	cMI->cd(2);
	histo2sMI->SetXTitle("%Diff");
	histo2sMI->GetXaxis()->SetTitleSize(0.05);
	histo2sMI->GetYaxis()->SetLabelSize(0.05);
	histo2sMI->GetXaxis()->SetLabelSize(0.05);
	histo2sMI->Draw();
	cMI->cd(3);
	histo3sMI->SetXTitle("%Diff");
	histo3sMI->GetXaxis()->SetTitleSize(0.05);
	histo3sMI->GetYaxis()->SetLabelSize(0.05);
	histo3sMI->GetXaxis()->SetLabelSize(0.05);
	histo3sMI->Draw();
	
	cRfb->cd(1);
	histo1sRfb->SetXTitle("%Diff");
	histo1sRfb->GetXaxis()->SetTitleSize(0.05);
	histo1sRfb->GetYaxis()->SetLabelSize(0.05);
	histo1sRfb->GetXaxis()->SetLabelSize(0.05);
	histo1sRfb->SetStats(kTRUE);
	histo1sRfb->Draw();
	cRfb->cd(2);
	histo2sRfb->SetXTitle("%Diff");
	histo2sRfb->GetXaxis()->SetTitleSize(0.05);
	histo2sRfb->GetYaxis()->SetLabelSize(0.05);
	histo2sRfb->GetXaxis()->SetLabelSize(0.05);
	histo2sRfb->Draw();
	cRfb->cd(3);
	histo3sRfb->SetXTitle("%Diff");
	histo3sRfb->GetXaxis()->SetTitleSize(0.05);
	histo3sRfb->GetYaxis()->SetLabelSize(0.05);
	histo3sRfb->GetXaxis()->SetLabelSize(0.05);
	histo3sRfb->Draw();
	
	TNtuple* ntupleSig = new TNtuple("ntupleSig","Yields from fits","nSig1sPLalt:nSig2sPLalt:nSig3sPLalt:nSig1sPLnom:nSig2sPLnom:nSig3sPLnom:nSig1sMIalt:nSig2sMIalt:nSig3sMIalt:nSig1sMInom:nSig2sMInom:nSig3sMInom",numTrials);
	TNtuple* ntupleDiff = new TNtuple("ntupleDiff","Differences from fits","diff1sPL:diff2sPL:diff3sPL:diff1sMI:diff2sMI:diff3sMI:diff1sRfb:diff2sRfb:diff3sRfb",numTrials);
	TNtuple* ntupleChisq = new TNtuple("ntupleChisq","Chisq/ndf from fits","chisqndfPLalt:chisqndfPLnom:chisqndfMIalt:chisqndfMInom",numTrials);
	
	int nFailedTrials = 0;
	
	cout << "RUNNING PSEUDOEXPERIMENTS" << endl;
	for (int i = 0; i < numTrials; i++)
	{
		/*if (nFailedTrials > 10 && nFailedTrials > i*4)
		{
			cout << "TOO MANY FAILED TRIALS, ABORTING" << endl;
			break;
		}*/
		
		//RooWorkspace* mws;
		vector<double>* resultVector = new vector<double>(4);
		//Generate PL pseudodata
		RooDataSet* pseudoData = genModelPL->generate(*(wsgenPL->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//PL Alternate fit
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,resultVector,pseudoData);
		float nSig1sPLalt = (*resultVector)[1];
		float nSig2sPLalt = (*resultVector)[2];
		float nSig3sPLalt = (*resultVector)[3];
		float chisqndfPLalt = (*resultVector)[0];
		if (chisqndfPLalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		//PL Nominal fit
		FitData(kPADATA,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,0,resultVector,pseudoData);
		float nSig1sPLnom = (*resultVector)[1];
		float nSig2sPLnom = (*resultVector)[2];
		float nSig3sPLnom = (*resultVector)[3];
		float chisqndfPLnom = (*resultVector)[0];
		if (chisqndfPLnom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		//Generate MI pseudodata
		delete pseudoData;
		pseudoData = genModelMI->generate(*(wsgenMI->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//MI Alternate fit
		FitData(kPADATA,ptLow,ptHigh,yLowMI,yHighMI,cLow,cHigh,muPtCut,whichModel,resultVector,pseudoData);
		float nSig1sMIalt = (*resultVector)[1];
		float nSig2sMIalt = (*resultVector)[2];
		float nSig3sMIalt = (*resultVector)[3];
		float chisqndfMIalt = (*resultVector)[0];
		if (chisqndfMIalt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		//MI Nominal fit
		FitData(kPADATA,ptLow,ptHigh,yLowMI,yHighMI,cLow,cHigh,muPtCut,0,resultVector,pseudoData);
		float nSig1sMInom = (*resultVector)[1];
		float nSig2sMInom = (*resultVector)[2];
		float nSig3sMInom = (*resultVector)[3];
		float chisqndfMInom = (*resultVector)[0];
		if (chisqndfMInom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		cout << "FINISHED FITS" << endl;
		
		float Rfb1sAlt = nSig1sPLalt/nSig1sMIalt;
		float Rfb2sAlt = nSig2sPLalt/nSig2sMIalt;
		float Rfb3sAlt = nSig3sPLalt/nSig3sMIalt;
		
		float Rfb1sNom = nSig1sPLnom/nSig1sMInom;
		float Rfb2sNom = nSig2sPLnom/nSig2sMInom;
		float Rfb3sNom = nSig3sPLnom/nSig3sMInom;
		
		float diff1sPL = 100*(nSig1sPLalt-nSig1sPLnom)/nSig1sPLnom;
		float diff2sPL = 100*(nSig2sPLalt-nSig2sPLnom)/nSig2sPLnom;
		float diff3sPL = 100*(nSig3sPLalt-nSig3sPLnom)/nSig3sPLnom;
		
		float diff1sMI = 100*(nSig1sMIalt-nSig1sMInom)/nSig1sMInom;
		float diff2sMI = 100*(nSig2sMIalt-nSig2sMInom)/nSig2sMInom;
		float diff3sMI = 100*(nSig3sMIalt-nSig3sMInom)/nSig3sMInom;
		
		float diff1sRfb = 100*(Rfb1sAlt-Rfb1sNom)/Rfb1sNom;
		float diff2sRfb = 100*(Rfb2sAlt-Rfb2sNom)/Rfb2sNom;
		float diff3sRfb = 100*(Rfb3sAlt-Rfb3sNom)/Rfb3sNom;
		
		cout << "CALCULATED VALUES" << endl;
		
		//Store data in hists/tuples
		outfile->cd();
		
		histo1sPL->Fill(diff1sPL);
		histo2sPL->Fill(diff2sPL);
		histo3sPL->Fill(diff3sPL);
		histo1sMI->Fill(diff1sMI);
		histo2sMI->Fill(diff2sMI);
		histo3sMI->Fill(diff3sMI);
		histo1sRfb->Fill(diff1sRfb);
		histo2sRfb->Fill(diff2sRfb);
		histo3sRfb->Fill(diff3sRfb);
		
		cout << "FILLED HISTOGRAMS" << endl;
		
		ntupleSig->Fill(nSig1sPLalt,nSig2sPLalt,nSig3sPLalt,nSig1sPLnom,nSig2sPLnom,nSig3sPLnom,nSig1sMIalt,nSig2sMIalt,nSig3sMIalt,nSig1sMInom,nSig2sMInom,nSig3sMInom);
		ntupleDiff->Fill(diff1sPL,diff2sPL,diff3sPL,diff1sMI,diff2sMI,diff3sMI,diff1sRfb,diff2sRfb,diff3sRfb);
		ntupleChisq->Fill(chisqndfPLalt,chisqndfPLnom,chisqndfMIalt,chisqndfMInom);
		
		cout << "FILLED TUPLES" << endl;
		
		cPL->Update();
		cPL->cd(1);
		histo1sPL->Draw();
		cPL->cd(2);
		histo2sPL->Draw();
		cPL->cd(3);
		histo3sPL->Draw();
		
		cMI->Update();
		cMI->cd(1);
		histo1sMI->Draw();
		cMI->cd(2);
		histo2sMI->Draw();
		cMI->cd(3);
		histo3sMI->Draw();
		
		cRfb->Update();
		cRfb->cd(1);
		histo1sRfb->Draw();
		cRfb->cd(2);
		histo2sRfb->Draw();
		cRfb->cd(3);
		histo3sRfb->Draw();
		
		cout << "DREW HISTOGRAMS" << endl << endl;
		
		cout << "RESULTS OF PSEUDOEXPERIMENT " << i << endl;
		cout << "PL ALT 1S = " << nSig1sPLalt << endl;
		cout << "PL ALT 2S = " << nSig2sPLalt << endl;
		cout << "PL ALT 3S = " << nSig3sPLalt << endl;
		cout << "PL Alt CHISQ/NDF = " << chisqndfPLalt << endl;
		cout << "PL NOM 1S = " << nSig1sPLnom << endl;
		cout << "PL NOM 2S = " << nSig2sPLnom << endl;
		cout << "PL NOM 3S = " << nSig3sPLnom << endl;
		cout << "PL NOM CHISQ/NDF = " << chisqndfPLnom << endl;
		cout << "MI ALT 1S = " << nSig1sMIalt << endl;
		cout << "MI ALT 2S = " << nSig2sMIalt << endl;
		cout << "MI ALT 3S = " << nSig3sMIalt << endl;
		cout << "MI Alt CHISQ/NDF = " << chisqndfMIalt << endl;
		cout << "MI NOM 1S = " << nSig1sMInom << endl;
		cout << "MI NOM 2S = " << nSig2sMInom << endl;
		cout << "MI NOM 3S = " << nSig3sMInom << endl;
		cout << "MI NOM CHISQ/NDF = " << chisqndfMInom << endl;
		
		delete pseudoData;
	}
	
	delete NomwsMI;
	delete NomwsPL;
	
	//Write results to file
	outfile->cd();
	histo1sPL->Write();
	histo2sPL->Write();
	histo3sPL->Write();
	histo1sMI->Write();
	histo2sMI->Write();
	histo3sMI->Write();
	histo1sRfb->Write();
	histo2sRfb->Write();
	histo3sRfb->Write();
	ntupleSig->Write();
	ntupleDiff->Write();
	ntupleChisq->Write();
	
	cout << nFailedTrials << endl;
	
	outfile->Close();
}