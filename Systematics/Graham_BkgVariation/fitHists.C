/******************************
Fits and draws several histograms loaded from a file.
Use makeHists.C to make the histogram file.
******************************/

#include "myFitFunctions.C"

void fitHists()
{
	//Program parameters
	const int nPtBins = 5; //For macros, Root wants array sizes to be constants
	const int mMin = 6;
	const int mMax = 20;
	
	//Open histogram file
	//TFile* histInFile = new TFile("../../../histograms.root","READ");
	TFile* histInFile = new TFile("../../../histogramsPP.root","READ");
	
	//Load histograms from file
	TH1D* mH[nPtBins];
	TH1D* mH_oppS[nPtBins];
	TH1D* MCmH[nPtBins];
	for (int h = 0; h < nPtBins; h++)
	{
		mH[h] = (TH1D*)histInFile->Get(Form("mH_%d",h));
		mH_oppS[h] = (TH1D*)histInFile->Get(Form("mH_oppS_%d",h));
		MCmH[h] = (TH1D*)histInFile->Get(Form("MCmH_%d",h));
	}
	TH2D* muplH = (TH2D*)histInFile->Get("muplH");
	TH2D* MCmuplH = (TH2D*)histInFile->Get("MCmuplH");
	TH2D* muplRatioH = (TH2D*)histInFile->Get("muplRatioH");
	TH2D* MCmuplWeightedH = (TH2D*)histInFile->Get("MCmuplWeightedH");
	
	////////////Fitting
	TF1* mFit[nPtBins];
	double mFitParam[7][nPtBins];
	
	//Fit MC in individual pT bins (index 0 is integrated pT bin)
	for (int h = 1; h < nPtBins; h++)
	{
		mFit[h] = new TF1(Form("mFit_%d",h),fSumErfExp,mMin,mMax,7);
		mFit[h]->SetParNames("Norm","Ratio","erf #mu","erf #sigma","exp #mu","exp #lambda","decay #lambda");
		mFit[h]->SetParameters(20,0.8,8,0.5,8,1,3.4);
		mFit[h]->SetParLimits(2,6,10);	//erf mu
		mFit[h]->SetParLimits(3,0,20);	//erf sigma
		mFit[h]->SetParLimits(4,6,8);	//exp mu
		mFit[h]->SetParLimits(5,0,20);	//exp lambda
		mFit[h]->SetParLimits(6,0,100); //decay lambda
		if (h == 1)
		{
			mFit[h]->SetParLimits(2,6,11);
			mFit[h]->SetParLimits(3,0.3,20);
		}
		if (h == 3)
		{
			mFit[h]->SetParLimits(2,6,11);
			mFit[h]->SetParLimits(5,0.5,20);
		}
		if (h == 5)
			mFit[h]->SetParLimits(4,6,7); //reduced upper limit on exp mu in 4-5 pt bin to prevent kink
		
		MCmH[h]->Fit(mFit[h]);
		
		//Store fit parameters to use later
		TF1* hf = MCmH[h]->GetFunction(Form("mFit_%d", h));
		for (int p = 0; p < 7; p++)
		{
			mFitParam[p][h] = hf->GetParameter(p);
		}
	}
	
	//Fit data in pT 0-5 bin
	mFit[0] = new TF1("mFit_0",fSumErfExpTotal3,mMin,mMax,35);
	mFit[0]->SetParName(0,"A1");
	mFit[0]->SetParName(7,"A2");
	mFit[0]->SetParName(14,"A3");
	mFit[0]->SetParName(21,"A4");
	mFit[0]->SetParName(28,"A5");
	for (int h = 1; h < nPtBins; h++)
	{
		for (int p = 1; p < 7; p++)
		{
			mFit[0]->FixParameter((h-1)*7+p,mFitParam[p][h]); //Fix all parameters from MC fits except amplitudes
		}
	}
	mFit[0]->SetParLimits(0,0,1000); //Limit amplitudes to be positive!
	mFit[0]->SetParLimits(7,0,1000);
	mFit[0]->SetParLimits(14,0,1000);
	mFit[0]->SetParLimits(21,0,1000);
	mFit[0]->SetParLimits(28,0,1000);
	for (int i = 28; i < 35; i++)
		mFit[0]->FixParameter(i,0);
	mH[0]->Fit(mFit[0]);
	
	//Print parameters in format that can be copied into code
	cout << "\nPARAMETERS\n";
	for (int h = 1; h < nPtBins; h++)
	{
		cout << "Double_t ratio" << h << " = " << mFitParam[1][h] << ";\n";
		cout << "Double_t muErf" << h << " = " << mFitParam[2][h] << ";\n";
		cout << "Double_t sigma" << h << " = " << mFitParam[3][h] << ";\n";
		cout << "Double_t muExp" << h << " = " << mFitParam[4][h] << ";\n";
		cout << "Double_t lambdaExp" << h << " = " << mFitParam[5][h] << ";\n";
		cout << "Double_t lambdaDecay" << h << " = " << mFitParam[6][h] << ";\n";
	}
	
	/*
	//Fit exp decay to tail of data
	TF1* mTailFit[nPtBins];
	for (int h = 0; h < nPtBins; h++)
	{
		mTailFit[h] = new TF1(Form("mTailFit_%d",h),"[0]*TMath::Exp(-x/[1])",11,20);
		mTailFit[h]->SetParameters(100,3.4);
		mTailFit[h]->SetParLimits(1,2,10);
		mH[h]->Fit(mTailFit[h],"R");
	}
	*/
	
	/////////////Draw histograms//////////
	
	//Draw PtRap histograms
	TCanvas* muc = new TCanvas("muc","MuPlus Pt and Eta",800,600);
	muc->Divide(2,2);
	muc->cd(1);
	muplH->Draw("COL");
	muplH->SetStats(0);
	muplH->SetXTitle("#eta");
	muplH->SetYTitle("pT");
	muc->cd(2);
	MCmuplH->Draw("COL");
	MCmuplH->SetStats(0);
	MCmuplH->SetXTitle("Rapidity");
	MCmuplH->SetYTitle("pT");
	muc->cd(3);
	muplRatioH->Draw("COL");
	muplRatioH->SetStats(0);
	muplRatioH->SetXTitle("Rapidity");
	muplRatioH->SetYTitle("pT");
	muc->cd(4);
	MCmuplWeightedH->Draw("COL");
	MCmuplWeightedH->SetStats(0);
	MCmuplWeightedH->SetXTitle("Rapidity");
	MCmuplWeightedH->SetYTitle("pT");
	
	//Draw mass histograms
	TCanvas* massC = new TCanvas("massC","Mass plots",800,600);
	massC->Divide(2,3);
	
	//Draw data for pT 0-5 and MC for individual bins
	massC->cd(1);
	gStyle->SetOptStat(1000000001);
	gStyle->SetOptFit(1);
	gStyle->SetTitleFontSize(.05);
	mH[0]->Draw();
	mH[0]->SetXTitle("Mass (GeV)");
	mH[0]->SetYTitle("Count");
	mH[0]->SetMinimum(0);
	for (int h = 1; h < nPtBins; h++)
	{
		massC->cd(h + 1);
		gStyle->SetOptStat(1000000001);
		gStyle->SetOptFit(1);
		gStyle->SetTitleFontSize(.06);
		MCmH[h]->Draw();
		
		MCmH[h]->SetLineColor(kGreen+1);
		MCmH[h]->SetXTitle("Mass (GeV)");
		MCmH[h]->SetYTitle("Count");
		MCmH[h]->SetMinimum(0);
	}
	
	/*
	//Draw data and MC together
	for (int h = 0; h < nPtBins; h++)
	{
		massC->cd(h + 1);
		gStyle->SetOptStat(0000000000);
		gStyle->SetOptFit(0);
		gStyle->SetTitleFontSize(.05);
		mH[h]->Draw();
		MCmH[h]->Draw("SAME");
		//mH_oppS[h]->Draw("SAME");
		
		mH[h]->SetXTitle("Mass (GeV)");
		mH[h]->SetYTitle("Count");
		mH[h]->SetMinimum(0);
		MCmH[h]->SetLineColor(kGreen+1);
		//mH_oppS[h]->SetLineColor(kViolet+1);
		//mH[h]->SetMaximum(mH_oppS[h]->GetMaximum());
	}
	*/
	
	double axTitleSize = 0.05;
	double axLabelSize = 0.06;
	for (int h = 0; h < nPtBins; h++)
	{
		mH[h]->GetXaxis()->SetTitleSize(axTitleSize);
		mH[h]->GetXaxis()->SetLabelSize(axLabelSize);
		mH[h]->GetYaxis()->SetTitleSize(axTitleSize);
		mH[h]->GetYaxis()->SetLabelSize(axLabelSize);
		MCmH[h]->GetXaxis()->SetTitleSize(axTitleSize);
		MCmH[h]->GetXaxis()->SetLabelSize(axLabelSize);
		MCmH[h]->GetYaxis()->SetTitleSize(axTitleSize);
		MCmH[h]->GetYaxis()->SetLabelSize(axLabelSize);
	}
	
}