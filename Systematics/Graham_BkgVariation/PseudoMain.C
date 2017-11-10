#include "RunPseudoExpts.C"

void PseudoMain(int whichState = 1, int firstBin = 0, int lastBin = -1, int numTrials = 100, double chisqmax = 10.0)
{
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
	else
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
		
		
	for (int i = firstBin; i <= lastBin; i++)
	{
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
		
		//Set which alternative model to use.
		int whichModel = 1;
		if (ptLow > 5.0 || (ptLow == 4.0 && ptHigh == 9.0))
			whichModel = 3;
		
		cout << Form("RunPseudoExpts(%.1f,%.1f,%.2f,%.2f,%d,%d,%.1f)",ptLow,ptHigh,yLow,yHigh,whichModel,numTrials,chisqmax) << endl;
		RunPseudoExpts(ptLow,ptHigh,yLow,yHigh,whichModel,numTrials,chisqmax);
	}
}