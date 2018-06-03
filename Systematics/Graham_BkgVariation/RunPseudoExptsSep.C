/*
Run pseudoexperiments for a single bin, for either PA -or- PP.
*/

#include "FitData.C"

#include "../../HeaderFiles/cutsAndBin.h"

#include "TTimeStamp.h"
#include "RooRandom.h"

using std::vector;

void RunPseudoExptsSep(
		int collId = kPADATA,
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
	
	TTimeStamp startTime;
	
	//Set up nominal files, workspaces, and generators
	TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
	TString NomFileName = Form("../../../JaredNomFits/nomfitresults_upsilon_%s.root",kineLabel.Data());
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
	delete NomFile;
	
	
	//Data for output
	const char* collLabel;
	if (collId == kPADATA)
		collLabel = "PA";
	else if (collId == kPPDATA)
		collLabel = "PP";
	TFile* outfile = new TFile(Form("ResultsBkg/PseudoExpResults_%s_pt%.1f-%.1f_y%.2f-%.2f.root",collLabel,ptLow,ptHigh,yLow,yHigh),"RECREATE");
	
	TNtuple* ntupleSig = new TNtuple("ntupleSig","Yields from fits","nSig1sAlt:nSig2sAlt:nSig3sAlt:nSig1sNom:nSig2sNom:nSig3sNom",numTrials);
	TNtuple* ntupleChisq = new TNtuple("ntupleChisq","Chisq/ndf from fits","chisqndfAlt:chisqndfNom",numTrials);
	
	//Set Random Seed
	RooRandom::randomGenerator()->SetSeed(startTime.GetSec());
	
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
		//Generate pseudodata
		RooDataSet* pseudoData = genModel->generate(*(wsgen->var("mass")));
		pseudoData->SetName("reducedDS");
		
		//Alternate fit
		FitData(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,resultVector,pseudoData,NomFileName);
		float nSig1sAlt = (*resultVector)[1];
		float nSig2sAlt = (*resultVector)[2];
		float nSig3sAlt = (*resultVector)[3];
		float chisqndfAlt = (*resultVector)[0];
		if (chisqndfAlt > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		//Nominal fit
		FitData(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,0,resultVector,pseudoData,NomFileName);
		float nSig1sNom = (*resultVector)[1];
		float nSig2sNom = (*resultVector)[2];
		float nSig3sNom = (*resultVector)[3];
		float chisqndfNom = (*resultVector)[0];
		if (chisqndfNom > maxchisq)
		{
			i--;
			nFailedTrials++;
			cout << "CHISQ TOO HIGH, TRIAL FAILED" << endl;
			delete pseudoData;
			continue;
		}
		
		cout << "FINISHED FITS" << endl;
		
		//Store data in tuples
		outfile->cd();
		
		ntupleSig->Fill(nSig1sAlt,nSig2sAlt,nSig3sAlt,nSig1sNom,nSig2sNom,nSig3sNom);
		ntupleChisq->Fill(chisqndfAlt,chisqndfNom);
		
		cout << "FILLED TUPLES" << endl;
		
		cout << "RESULTS OF PSEUDOEXPERIMENT " << i << endl;
		cout << "ALT 1S = " << nSig1sAlt << endl;
		cout << "ALT 2S = " << nSig2sAlt << endl;
		cout << "ALT 3S = " << nSig3sAlt << endl;
		cout << "Alt CHISQ/NDF = " << chisqndfAlt << endl;
		cout << "NOM 1S = " << nSig1sNom << endl;
		cout << "NOM 2S = " << nSig2sNom << endl;
		cout << "NOM 3S = " << nSig3sNom << endl;
		cout << "NOM CHISQ/NDF = " << chisqndfNom << endl;
		
		delete pseudoData;
	}
	
	delete Nomws;
	
	//Write results to file
	outfile->cd();
	ntupleSig->Write();
	ntupleChisq->Write();
	
	cout << "Failed trials: " << nFailedTrials << endl;
	
	outfile->Close();
}