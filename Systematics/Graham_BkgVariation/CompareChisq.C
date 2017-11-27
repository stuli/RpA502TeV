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


void CompareChisq(int collId = kPADATA, int whichState = 1, int firstBin = 0, int lastBin = -1)
{
	int cLow=0;
	int cHigh=200;
	float muPtCut=4.0;

	float dphiEp2Low = 0;
	float dphiEp2High = 100;
	
	float ptbins1s[7] = {0,2,4,6,9,12,30};
	float ybins1s[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ptbins2s[4] = {0,4,9,30};
	float ybins2s[5] = {-1.93,-0.8,0.0,0.8,1.93};
	float ptbins3s[3] = {0,6,30};
	float ybins3s[3] = {-1.93,0.0,1.93};
	
	//Choose which set of bins to use based on Upsilon state
	float* ptbins;
	float* ybins;
	int numptbins;
	int numybins;
	if (whichState == 1)
	{
		ptbins = ptbins1s;
		ybins = ybins1s;
		numptbins = 6;
		numybins = 8;
	}
	else if (whichState == 2)
	{
		ptbins = ptbins2s;
		ybins = ybins2s;
		numptbins = 3;
		numybins = 4;
	}
	else if (whichState == 3)
	{
		ptbins = ptbins3s;
		ybins = ybins3s;
		numptbins = 2;
		numybins = 2;
	}
	else if (whichState != 0)
	{
		cout << "Invalid state specified. Aborting" << endl;
		return;
	}
	int totalnumbins = numptbins + numybins;
	
	//If no lastBin is specified, do only the first bin
	if (lastBin < firstBin)
		lastBin = firstBin;
	if (lastBin > totalnumbins)
		lastBin = totalnumbins;
		
	TString lineValues[28];
	for (int i = 0; i < 28; i++)
		lineValues[i] = Form("& - & - & -");
	
	//If whichState passed in is 0, do everything
	int maxState = whichState;
	if (whichState == 0)
	{
		maxState = 3;
		whichState = 1;
		firstBin = 0;
		lastBin = 99;
	}
	for (int iState = whichState; iState <= maxState; iState++)
	{
		for (int i = firstBin; i <= lastBin; i++)
		{
			if (iState == 1)
			{
				ptbins = ptbins1s;
				ybins = ybins1s;
				numptbins = 6;
				numybins = 8;
			}
			else if (iState == 2)
			{
				ptbins = ptbins2s;
				ybins = ybins2s;
				numptbins = 3;
				numybins = 4;
			}
			else if (iState == 3)
			{
				ptbins = ptbins3s;
				ybins = ybins3s;
				numptbins = 2;
				numybins = 2;
			}
			else
			{
				cout << "Invalid state specified. Aborting" << endl;
				return;
			}
			int totalnumbins = numptbins + numybins;
			
			if (lastBin < firstBin)
				lastBin = firstBin;
			if (lastBin > totalnumbins)
				lastBin = totalnumbins;
			
			float ptLow, ptHigh, yLow, yHigh;
			
			//Set pt, y range for this bin
			if (i == 0)
			{
				ptLow = 0;
				ptHigh = 30;
				yLow = -1.93;
				yHigh = 1.93;
			}
			else if (i <= numptbins)
			{
				ptLow = ptbins[i-1];
				ptHigh = ptbins[i];
				yLow = -1.93;
				yHigh = 1.93;
			}
			else
			{
				ptLow = 0;
				ptHigh = 30;
				yLow = ybins[i-numptbins-1];
				yHigh = ybins[i-numptbins];
			}
			
			if (collId == kPPDATA)
			{
				float yLowPP = yLow;
				float yHighPP = yHigh;
				if (yLow < 0.0 && yHigh > 0.0)
					yLowPP = 0.0;
				else if (yLow < 0.0)
				{
					yLowPP = TMath::Abs(yHigh);
					yHighPP = TMath::Abs(yLow);
				}
				yLow = yLowPP;
				yHigh = yHighPP;
			}
			
			TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
			TString altFileName = Form("ResultsBkg/altfitresults_upsilon_%s.root",kineLabel.Data());
			TString powFileName = Form("ResultsBkg/powfitresults_upsilon_%s.root",kineLabel.Data());
			TString chebFileName = Form("ResultsBkg/chebfitresults_upsilon_%s.root",kineLabel.Data());
			
			TString altChisqndfStr = "-";
			TFile* altFile = new TFile(altFileName,"READ");
			if (!altFile->IsZombie())
			{
				RooWorkspace* altWorkspace = (RooWorkspace*)altFile->Get("workspace");
				double altChisqndf = altWorkspace->var("chisqndf")->getVal();
				altChisqndfStr = Form("%.2f",altChisqndf);
				delete altWorkspace;
				cout << "Found " << altFileName << endl;
			}
			TString powChisqndfStr = "-";
			TFile* powFile = new TFile(powFileName,"READ");
			if (!powFile->IsZombie())
			{
				RooWorkspace* powWorkspace = (RooWorkspace*)powFile->Get("workspace");
				double powChisqndf = powWorkspace->var("chisqndf")->getVal();
				powChisqndfStr = Form("%.2f",powChisqndf);
				delete powWorkspace;
				cout << "Found " << powFileName << endl;
			}
			TString chebChisqndfStr = "-";
			TFile* chebFile = new TFile(chebFileName,"READ");
			if (!chebFile->IsZombie())
			{
				RooWorkspace* chebWorkspace = (RooWorkspace*)chebFile->Get("workspace");
				double chebChisqndf = chebWorkspace->var("chisqndf")->getVal();
				chebChisqndfStr = Form("%.2f",chebChisqndf);
				delete chebWorkspace;
				cout << "Found " << chebFileName << endl;
			}
			
			int lineNumber = i;
			if (iState > 1) lineNumber += 15;
			if (iState > 2) lineNumber += 8;

			lineValues[lineNumber] = Form("& ") + altChisqndfStr + " & " + powChisqndfStr + " & " + chebChisqndfStr;
		}
	}
	cout << "PRINTING LATEX TABLES\n\n";
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:bkgAltVsChebChisqPA1S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & Combination & Power Law & Chebychev\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt$, y integrated " << lineValues[0] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 2 \\GeVc$ " << lineValues[1] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc$ " << lineValues[2] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc$ " << lineValues[3] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc$ " << lineValues[4] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc$ " << lineValues[5] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc$ " << lineValues[6] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-1.2$ " << lineValues[7] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8$ " << lineValues[8] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4$ " << lineValues[9] << "\\\\" << endl;
	cout << "$-1.4<y<0.0$ " << lineValues[10] << "\\\\" << endl;
	cout << "$0.0<y<0.4$ " << lineValues[11] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[12] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[13] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[14] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Comparison of chi^{2}/ndf between linear combination, power law, and 4th order Chebychev background PDFs in fits to pPb data for the Upsilon 1S bins.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:bkgAltVsChebChisqPA2S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & Combination & Power Law & Chebychev\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt$, y integrated " << lineValues[15] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 4 \\GeVc$ " << lineValues[16] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc$ " << lineValues[17] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc$ " << lineValues[18] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-0.8$ " << lineValues[19] << "\\\\" << endl;
	cout << "$-0.8<y<0.0$ " << lineValues[20] << "\\\\" << endl;
	cout << "$0.0<y<0.8$ " << lineValues[21] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[22] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Comparison of chi^{2}/ndf between linear combination, power law, and 4th order Chebychev background PDFs in fits to pPb data for the Upsilon 2S bins.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:bkgAltVsChebChisqPA3S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & Combination & Power Law & Chebychev\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt$, y integrated " << lineValues[23] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 6 \\GeVc$ " << lineValues[24] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc$ " << lineValues[25] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<0.0$ " << lineValues[26] << "\\\\" << endl;
	cout << "$0.0<y<1.93$ " << lineValues[27] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Comparison of chi^{2}/ndf between linear combination, power law, and 4th order Chebychev background PDFs in fits to pPb data for the Upsilon 3S bins.}" << endl;
	cout << "\\end{table}" << endl;
}