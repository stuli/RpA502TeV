void calcError(TNtuple* tup1, double* results, int state, TH1F* h1);
void calcError(TNtuple* tup1, TNtuple* tup2, double* results, int state, TH1F* h1, TH1F* h2, TH1F* hR);

void GetErrorsFromPseudoSep(int whichState = 0, int firstBin = 0, int lastBin = -1)
{
	float ptbins1s[7] = {0,2,4,6,9,12,30};
	float ybins1s[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	//float ybins1s[10] = {-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ptbins2s[4] = {0,4,9,30};
	float ybins2s[5] = {-1.93,-0.8,0.0,0.8,1.93};
	//float ybins2s[6] = {-2.4,-1.93,-0.8,0.0,0.8,1.93};
	float ptbins3s[3] = {0,6,30};
	float ybins3s[3] = {-1.93,0.0,1.93};
	//float ybins3s[4] = {-2.4,-1.93,0.0,1.93};
	
	float ybins1s_db[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ybins2s_db[5] = {-1.93,-0.8,0.0,0.8,1.93};
	float ybins3s_db[3] = {-1.93,0.0,1.93};
	
	float ybins1s_cr[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ybins2s_cr[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
	float ybins3s_cr[4] = {-2.87,-1.93,0.0,1.93};
	
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
		//numybins = 9;
	}
	else if (whichState == 2)
	{
		ptbins = ptbins2s;
		ybins = ybins2s;
		numptbins = 3;
		numybins = 4;
		//numybins = 5;
	}
	else if (whichState == 3)
	{
		ptbins = ptbins3s;
		ybins = ybins3s;
		numptbins = 2;
		numybins = 2;
		//numybins = 3;
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
	
	const int numLines = 81+56+3+11+24+2; //standard and diff + nt and hf + extra cr (-2.87) + cr pt + nt and hf with yInt + int rfb
	TString lineValues[numLines];
	for (int i = 0; i < numLines; i++)
		lineValues[i] = Form("& - & - & -");
		
	TString outNames[numLines];
	TString figCaption[numLines];
	TString figureLabel[numLines];
	
	TNtuple* resultTuple = new TNtuple("bkgSystTuple","Systematic Errors from Background Variation","ptLow:ptHigh:yLow:yHigh:yieldErrPA:yieldErrPP:RpAErr");
	//resultTuple->Fill(0,6,-1.93,1.93,0.204,.0601,0.1358);
	//resultTuple->Fill(0,30,0,1.93,0.1682,0.0525,0.2330);
	
	//HISTOGRAMS
	//Standard bins
	TH1F* hInt1S_PADiff = new TH1F("hInt1S_PADiff","",1,-1.93,1.93);
	TH1F* hInt2S_PADiff = new TH1F("hInt2S_PADiff","",1,-1.93,1.93);
	TH1F* hInt3S_PADiff = new TH1F("hInt3S_PADiff","",1,-1.93,1.93);
	TH1F* hpt1S_PADiff = new TH1F("hpt1S_PADiff","",6,ptbins1s);
	TH1F* hpt2S_PADiff = new TH1F("hpt2S_PADiff","",3,ptbins2s);
	TH1F* hpt3S_PADiff = new TH1F("hpt3S_PADiff","",2,ptbins3s);
	TH1F* hy1S_PADiff = new TH1F("hy1S_PADiff","",8,ybins1s);
	TH1F* hy2S_PADiff = new TH1F("hy2S_PADiff","",4,ybins2s);
	TH1F* hy3S_PADiff = new TH1F("hy3S_PADiff","",2,ybins3s);
	
	TH1F* hInt1S_PPDiff = new TH1F("hInt1S_PPDiff","",1,-1.93,1.93);
	TH1F* hInt2S_PPDiff = new TH1F("hInt2S_PPDiff","",1,-1.93,1.93);
	TH1F* hInt3S_PPDiff = new TH1F("hInt3S_PPDiff","",1,-1.93,1.93);
	TH1F* hpt1S_PPDiff = new TH1F("hpt1S_PPDiff","",6,ptbins1s);
	TH1F* hpt2S_PPDiff = new TH1F("hpt2S_PPDiff","",3,ptbins2s);
	TH1F* hpt3S_PPDiff = new TH1F("hpt3S_PPDiff","",2,ptbins3s);
	TH1F* hy1S_PPDiff = new TH1F("hy1S_PPDiff","",8,ybins1s);
	TH1F* hy2S_PPDiff = new TH1F("hy2S_PPDiff","",4,ybins2s);
	TH1F* hy3S_PPDiff = new TH1F("hy3S_PPDiff","",2,ybins3s);
	
	TH1F* hInt1S_RpADiff = new TH1F("hInt1S_RpADiff","",1,-1.93,1.93);
	TH1F* hInt2S_RpADiff = new TH1F("hInt2S_RpADiff","",1,-1.93,1.93);
	TH1F* hInt3S_RpADiff = new TH1F("hInt3S_RpADiff","",1,-1.93,1.93);
	TH1F* hpt1S_RpADiff = new TH1F("hpt1S_RpADiff","",6,ptbins1s);
	TH1F* hpt2S_RpADiff = new TH1F("hpt2S_RpADiff","",3,ptbins2s);
	TH1F* hpt3S_RpADiff = new TH1F("hpt3S_RpADiff","",2,ptbins3s);
	TH1F* hy1S_RpADiff = new TH1F("hy1S_RpADiff","",8,ybins1s);
	TH1F* hy2S_RpADiff = new TH1F("hy2S_RpADiff","",4,ybins2s);
	TH1F* hy3S_RpADiff = new TH1F("hy3S_RpADiff","",2,ybins3s);
	
	//Cross section
	TH1F* hy1S_cr_PADiff = new TH1F("hy1S_cr_PADiff","",9,ybins1s_cr);
	TH1F* hy2S_cr_PADiff = new TH1F("hy2S_cr_PADiff","",5,ybins2s_cr);
	TH1F* hy3S_cr_PADiff = new TH1F("hy3S_cr_PADiff","",3,ybins3s_cr);
	
	TH1F* hpt1S_cr_PADiff = new TH1F("hpt1S_cr_PADiff","1S pt bins for y-2.87 to 1.93",6,ptbins1s);
	TH1F* hpt2S_cr_PADiff = new TH1F("hpt2S_cr_PADiff","2S pt bins for y-2.87 to 1.93",3,ptbins2s);
	TH1F* hpt3S_cr_PADiff = new TH1F("hpt3S_cr_PADiff","3S pt bins for y-2.87 to 1.93",2,ptbins3s);
	
	//Differential bins
	TH1F* hpt1S_ym_PADiff = new TH1F("hpt1S_ym_PADiff","",6,ptbins1s);
	TH1F* hpt2S_ym_PADiff = new TH1F("hpt2S_ym_PADiff","",3,ptbins2s);
	TH1F* hpt3S_ym_PADiff = new TH1F("hpt3S_ym_PADiff","",2,ptbins3s);
	TH1F* hpt1S_ym_PPDiff = new TH1F("hpt1S_ym_PPDiff","",6,ptbins1s);
	TH1F* hpt2S_ym_PPDiff = new TH1F("hpt2S_ym_PPDiff","",3,ptbins2s);
	TH1F* hpt3S_ym_PPDiff = new TH1F("hpt3S_ym_PPDiff","",2,ptbins3s);
	TH1F* hpt1S_ym_RpADiff = new TH1F("hpt1S_ym_RpADiff","",6,ptbins1s);
	TH1F* hpt2S_ym_RpADiff = new TH1F("hpt2S_ym_RpADiff","",3,ptbins2s);
	TH1F* hpt3S_ym_RpADiff = new TH1F("hpt3S_ym_RpADiff","",2,ptbins3s);
	TH1F* hpt1S_yp_PADiff = new TH1F("hpt1S_yp_PADiff","",6,ptbins1s);
	TH1F* hpt2S_yp_PADiff = new TH1F("hpt2S_yp_PADiff","",3,ptbins2s);
	TH1F* hpt3S_yp_PADiff = new TH1F("hpt3S_yp_PADiff","",2,ptbins3s);
	TH1F* hpt1S_yp_PPDiff = new TH1F("hpt1S_yp_PPDiff","",6,ptbins1s);
	TH1F* hpt2S_yp_PPDiff = new TH1F("hpt2S_yp_PPDiff","",3,ptbins2s);
	TH1F* hpt3S_yp_PPDiff = new TH1F("hpt3S_yp_PPDiff","",2,ptbins3s);
	TH1F* hpt1S_yp_RpADiff = new TH1F("hpt1S_yp_RpADiff","",6,ptbins1s);
	TH1F* hpt2S_yp_RpADiff = new TH1F("hpt2S_yp_RpADiff","",3,ptbins2s);
	TH1F* hpt3S_yp_RpADiff = new TH1F("hpt3S_yp_RpADiff","",2,ptbins3s);
	
	TH1F* hy1S_pt06_PADiff = new TH1F("hy1S_pt06_PADiff","",8,ybins1s_db);
	TH1F* hy2S_pt06_PADiff = new TH1F("hy2S_pt06_PADiff","",4,ybins2s_db);
	TH1F* hy3S_pt06_PADiff = new TH1F("hy3S_pt06_PADiff","",2,ybins3s_db);
	TH1F* hy1S_pt06_PPDiff = new TH1F("hy1S_pt06_PPDiff","",8,ybins1s_db);
	TH1F* hy2S_pt06_PPDiff = new TH1F("hy2S_pt06_PPDiff","",4,ybins2s_db);
	TH1F* hy3S_pt06_PPDiff = new TH1F("hy3S_pt06_PPDiff","",2,ybins3s_db);
	TH1F* hy1S_pt06_RpADiff = new TH1F("hy1S_pt06_RpADiff","",8,ybins1s_db);
	TH1F* hy2S_pt06_RpADiff = new TH1F("hy2S_pt06_RpADiff","",4,ybins2s_db);
	TH1F* hy3S_pt06_RpADiff = new TH1F("hy3S_pt06_RpADiff","",2,ybins3s_db);
	TH1F* hy1S_pt630_PADiff = new TH1F("hy1S_pt630_PADiff","",8,ybins1s_db);
	TH1F* hy2S_pt630_PADiff = new TH1F("hy2S_pt630_PADiff","",4,ybins2s_db);
	TH1F* hy3S_pt630_PADiff = new TH1F("hy3S_pt630_PADiff","",2,ybins3s_db);
	TH1F* hy1S_pt630_PPDiff = new TH1F("hy1S_pt630_PPDiff","",8,ybins1s_db);
	TH1F* hy2S_pt630_PPDiff = new TH1F("hy2S_pt630_PPDiff","",4,ybins2s_db);
	TH1F* hy3S_pt630_PPDiff = new TH1F("hy3S_pt630_PPDiff","",2,ybins3s_db);
	TH1F* hy1S_pt630_RpADiff = new TH1F("hy1S_pt630_RpADiff","",8,ybins1s_db);
	TH1F* hy2S_pt630_RpADiff = new TH1F("hy2S_pt630_RpADiff","",4,ybins2s_db);
	TH1F* hy3S_pt630_RpADiff = new TH1F("hy3S_pt630_RpADiff","",2,ybins3s_db);
	
	//Rfb bins
	float ybins1s_rfb[5] = {0,0.4,0.8,1.2,1.93};
	float ybins2s_rfb[3] = {0,0.8,1.93};
	float ybins3s_rfb[2] = {0,1.93};
	/*float ntbins[5] = {0,40,65,90,400};
	float hfbins[5] = {0,15,22,30,120};*/
	float ntbins[5] = {0,40,62,88,400};
	float hfbins[5] = {0,12,19,27,120};
	float ntbins3S[3] = {0,40,400};
	float hfbins3S[3] = {0,12,120};
	TH2F* hnt1S_yPLdiff = new TH2F("hnt1S_yPLdiff","",4,ybins1s_rfb,4,ntbins);
	TH2F* hnt2S_yPLdiff = new TH2F("hnt2S_yPLdiff","",2,ybins2s_rfb,4,ntbins);
	TH2F* hnt3S_yPLdiff = new TH2F("hnt3S_yPLdiff","",1,ybins3s_rfb,4,ntbins);
	TH2F* hnt1S_yMIdiff = new TH2F("hnt1S_yMIdiff","",4,ybins1s_rfb,4,ntbins);
	TH2F* hnt2S_yMIdiff = new TH2F("hnt2S_yMIdiff","",2,ybins2s_rfb,4,ntbins);
	TH2F* hnt3S_yMIdiff = new TH2F("hnt3S_yMIdiff","",1,ybins3s_rfb,4,ntbins);
	TH2F* hnt1S_yRfbdiff = new TH2F("hnt1S_yRfbdiff","",4,ybins1s_rfb,4,ntbins);
	TH2F* hnt2S_yRfbdiff = new TH2F("hnt2S_yRfbdiff","",2,ybins2s_rfb,4,ntbins);
	TH2F* hnt3S_yRfbdiff = new TH2F("hnt3S_yRfbdiff","",1,ybins3s_rfb,4,ntbins);
	
	TH2F* hhf1S_yPLdiff = new TH2F("hhf1S_yPLdiff","",4,ybins1s_rfb,4,hfbins);
	TH2F* hhf2S_yPLdiff = new TH2F("hhf2S_yPLdiff","",2,ybins2s_rfb,4,hfbins);
	TH2F* hhf3S_yPLdiff = new TH2F("hhf3S_yPLdiff","",1,ybins3s_rfb,4,hfbins);
	TH2F* hhf1S_yMIdiff = new TH2F("hhf1S_yMIdiff","",4,ybins1s_rfb,4,hfbins);
	TH2F* hhf2S_yMIdiff = new TH2F("hhf2S_yMIdiff","",2,ybins2s_rfb,4,hfbins);
	TH2F* hhf3S_yMIdiff = new TH2F("hhf3S_yMIdiff","",1,ybins3s_rfb,4,hfbins);
	TH2F* hhf1S_yRfbdiff = new TH2F("hhf1S_yRfbdiff","",4,ybins1s_rfb,4,hfbins);
	TH2F* hhf2S_yRfbdiff = new TH2F("hhf2S_yRfbdiff","",2,ybins2s_rfb,4,hfbins);
	TH2F* hhf3S_yRfbdiff = new TH2F("hhf3S_yRfbdiff","",1,ybins3s_rfb,4,hfbins);
	
	TH1F* hnt1S_yInt_PLdiff = new TH1F("hnt1S_yInt_PLdiff","Ntrack bins for 0<y<1.93",4,ntbins);
	TH1F* hnt2S_yInt_PLdiff = new TH1F("hnt2S_yInt_PLdiff","Ntrack bins for 0<y<1.93",4,ntbins);
	TH1F* hnt3S_yInt_PLdiff = new TH1F("hnt3S_yInt_PLdiff","Ntrack bins for 0<y<1.93",2,ntbins3S);
	TH1F* hnt1S_yInt_MIdiff = new TH1F("hnt1S_yInt_MIdiff","Ntrack bins for -1.93<y<0",4,ntbins);
	TH1F* hnt2S_yInt_MIdiff = new TH1F("hnt2S_yInt_MIdiff","Ntrack bins for -1.93<y<0",4,ntbins);
	TH1F* hnt3S_yInt_MIdiff = new TH1F("hnt3S_yInt_MIdiff","Ntrack bins for -1.93<y<0",2,ntbins3S);
	TH1F* hnt1S_yInt_Rfbdiff = new TH1F("hnt1S_yInt_Rfbdiff","Ntrack bins for Rfb, 0<|y|<1.93",4,ntbins);
	TH1F* hnt2S_yInt_Rfbdiff = new TH1F("hnt2S_yInt_Rfbdiff","Ntrack bins for Rfb, 0<|y|<1.93",4,ntbins);
	TH1F* hnt3S_yInt_Rfbdiff = new TH1F("hnt3S_yInt_Rfbdiff","Ntrack bins for Rfb, 0<|y|<1.93",2,ntbins3S);
	
	TH1F* hhf1S_yInt_PLdiff = new TH1F("hhf1S_yInt_PLdiff","HF bins for 0<y<1.93",4,hfbins);
	TH1F* hhf2S_yInt_PLdiff = new TH1F("hhf2S_yInt_PLdiff","HF bins for 0<y<1.93",4,hfbins);
	TH1F* hhf3S_yInt_PLdiff = new TH1F("hhf3S_yInt_PLdiff","HF bins for 0<y<1.93",2,hfbins3S);
	TH1F* hhf1S_yInt_MIdiff = new TH1F("hhf1S_yInt_MIdiff","HF bins for -1.93<y<0",4,hfbins);
	TH1F* hhf2S_yInt_MIdiff = new TH1F("hhf2S_yInt_MIdiff","HF bins for -1.93<y<0",4,hfbins);
	TH1F* hhf3S_yInt_MIdiff = new TH1F("hhf3S_yInt_MIdiff","HF bins for -1.93<y<0",2,hfbins3S);
	TH1F* hhf1S_yInt_Rfbdiff = new TH1F("hhf1S_yInt_Rfbdiff","HF bins for Rfb, 0<y<1.93",4,hfbins);
	TH1F* hhf2S_yInt_Rfbdiff = new TH1F("hhf2S_yInt_Rfbdiff","HF bins for Rfb, 0<y<1.93",4,hfbins);
	TH1F* hhf3S_yInt_Rfbdiff = new TH1F("hhf3S_yInt_Rfbdiff","HF bins for Rfb, 0<y<1.93",2,hfbins3S);
	/*
	TH1F* hInt1S_yPLdiff = new TH1F("hInt1S_yPLdiff","",1,-1.93,1.93);
	TH1F* hInt2S_yPLdiff = new TH1F("hInt2S_yPLdiff","",1,-1.93,1.93);
	TH1F* hInt3S_yPLdiff = new TH1F("hInt3S_yPLdiff","",1,-1.93,1.93);
	TH1F* hInt1S_yMIdiff = new TH1F("hInt1S_yMIdiff","",1,-1.93,1.93);
	TH1F* hInt2S_yMIdiff = new TH1F("hInt2S_yMIdiff","",1,-1.93,1.93);
	TH1F* hInt3S_yMIdiff = new TH1F("hInt3S_yMIdiff","",1,-1.93,1.93);
	TH1F* hInt1S_Rfbdiff = new TH1F("hInt1S_Rfbdiff","",1,-1.93,1.93);
	TH1F* hInt2S_Rfbdiff = new TH1F("hInt2S_Rfbdiff","",1,-1.93,1.93);
	TH1F* hInt3S_Rfbdiff = new TH1F("hInt3S_Rfbdiff","",1,-1.93,1.93);
	*/
	TH1F* hInt1S_yPLdiff = new TH1F("hInt1S_yPLdiff","",1,0,1.93);
	TH1F* hInt1S_yMIdiff = new TH1F("hInt1S_yMIdiff","",1,0,1.93);
	TH1F* hInt1S_Rfbdiff = new TH1F("hInt1S_Rfbdiff","",1,0,1.93);
	TH1F* hInt2S_yPLdiff = new TH1F("hInt2S_yPLdiff","",1,0,1.93);
	TH1F* hInt2S_yMIdiff = new TH1F("hInt2S_yMIdiff","",1,0,1.93);
	TH1F* hInt2S_Rfbdiff = new TH1F("hInt2S_Rfbdiff","",1,0,1.93);
	
	TH1F* hThisPA;
	TH1F* hThisPP;
	TH1F* hThisRpA;
	int outbin;
	
	//If whichState passed in is 0, do everything
	int maxState = whichState;
	if (whichState == 0)
	{
		maxState = 3;
		whichState = 1;
		firstBin = 0;
		lastBin = 99;
	}
	for (int diffBin = 0; diffBin < 3; diffBin++) //0=standard bins, 1=ptlow,ymi, 2=pthigh,ypl
	{
		lastBin = 99;
		for (int iState = whichState; iState <= maxState; iState++)
		{
			for (int i = firstBin; i <= lastBin+1; i++)
			{
				if (i==lastBin+1 && diffBin!=0)
					continue;
				
				if (iState == 1)
				{
					ptbins = ptbins1s;
					numptbins = 6;
					if (diffBin == 0)
						{ybins = ybins1s; numybins = 8;}
					else
						{ybins = ybins1s_db; numybins = 8;}
				}
				else if (iState == 2)
				{
					ptbins = ptbins2s;
					numptbins = 3;
					if (diffBin == 0)
						{ybins = ybins2s; numybins = 4;}
					else
						{ybins = ybins2s_db; numybins = 4;}
				}
				else if (iState == 3)
				{
					ptbins = ptbins3s;
					numptbins = 2;
					if (diffBin == 0)
						{ybins = ybins3s; numybins = 2;}
					else
						{ybins = ybins3s_db; numybins = 2;}
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
					if (diffBin != 0)
						continue; //no integrated differential bin
					ptLow = 0;
					ptHigh = 30;
					yLow = -1.93;
					yHigh = 1.93;
					if (iState == 1) {hThisPA = hInt1S_PADiff; hThisPP = hInt1S_PPDiff; hThisRpA = hInt1S_RpADiff;}
					else if (iState == 2) {hThisPA = hInt2S_PADiff; hThisPP = hInt2S_PPDiff; hThisRpA = hInt2S_RpADiff;}
					else if (iState == 3) {hThisPA = hInt3S_PADiff; hThisPP = hInt3S_PPDiff; hThisRpA = hInt3S_RpADiff;}
					outbin = 1;
				}
				else if (i <= numptbins)
				{
					ptLow = ptbins[i-1];
					ptHigh = ptbins[i];
					if (diffBin == 0)
					{
						yLow = -1.93;
						yHigh = 1.93;
						if (iState == 1) {hThisPA = hpt1S_PADiff; hThisPP = hpt1S_PPDiff; hThisRpA = hpt1S_RpADiff;}
						else if (iState == 2) {hThisPA = hpt2S_PADiff; hThisPP = hpt2S_PPDiff; hThisRpA = hpt2S_RpADiff;}
						else if (iState == 3) {hThisPA = hpt3S_PADiff; hThisPP = hpt3S_PPDiff; hThisRpA = hpt3S_RpADiff;}
					}
					else if (diffBin == 1)
					{
						yLow = -1.93;
						yHigh = 0;
						if (iState == 1) {hThisPA = hpt1S_ym_PADiff; hThisPP = hpt1S_ym_PPDiff; hThisRpA = hpt1S_ym_RpADiff;}
						else if (iState == 2) {hThisPA = hpt2S_ym_PADiff; hThisPP = hpt2S_ym_PPDiff; hThisRpA = hpt2S_ym_RpADiff;}
						else if (iState == 3) {hThisPA = hpt3S_ym_PADiff; hThisPP = hpt3S_ym_PPDiff; hThisRpA = hpt3S_ym_RpADiff;}
					}
					else if (diffBin == 2)
					{
						yLow = 0.00;
						yHigh = 1.93;
						if (iState == 1) {hThisPA = hpt1S_yp_PADiff; hThisPP = hpt1S_yp_PPDiff; hThisRpA = hpt1S_yp_RpADiff;}
						else if (iState == 2) {hThisPA = hpt2S_yp_PADiff; hThisPP = hpt2S_yp_PPDiff; hThisRpA = hpt2S_yp_RpADiff;}
						else if (iState == 3) {hThisPA = hpt3S_yp_PADiff; hThisPP = hpt3S_yp_PPDiff; hThisRpA = hpt3S_yp_RpADiff;}
					}
					outbin = i;
				}
				else if (i <= lastBin)
				{
					if (diffBin == 0)
					{
						ptLow = 0;
						ptHigh = 30;
						if (iState == 1) {hThisPA = hy1S_PADiff; hThisPP = hy1S_PPDiff; hThisRpA = hy1S_RpADiff;}
						else if (iState == 2) {hThisPA = hy2S_PADiff; hThisPP = hy2S_PPDiff; hThisRpA = hy2S_RpADiff;}
						else if (iState == 3) {hThisPA = hy3S_PADiff; hThisPP = hy3S_PPDiff; hThisRpA = hy3S_RpADiff;}
					}
					else if (diffBin == 1)
					{
						ptLow = 0;
						ptHigh = 6;
						if (iState == 1) {hThisPA = hy1S_pt06_PADiff; hThisPP = hy1S_pt06_PPDiff; hThisRpA = hy1S_pt06_RpADiff;}
						else if (iState == 2) {hThisPA = hy2S_pt06_PADiff; hThisPP = hy2S_pt06_PPDiff; hThisRpA = hy2S_pt06_RpADiff;}
						else if (iState == 3) {hThisPA = hy3S_pt06_PADiff; hThisPP = hy3S_pt06_PPDiff; hThisRpA = hy3S_pt06_RpADiff;}
					}
					else if (diffBin == 2)
					{
						ptLow = 6;
						ptHigh = 30;
						if (iState == 1) {hThisPA = hy1S_pt630_PADiff; hThisPP = hy1S_pt630_PPDiff; hThisRpA = hy1S_pt630_RpADiff;}
						else if (iState == 2) {hThisPA = hy2S_pt630_PADiff; hThisPP = hy2S_pt630_PPDiff; hThisRpA = hy2S_pt630_RpADiff;}
						else if (iState == 3) {hThisPA = hy3S_pt630_PADiff; hThisPP = hy3S_pt630_PPDiff; hThisRpA = hy3S_pt630_RpADiff;}
					}
					yLow = ybins[i-numptbins-1];
					yHigh = ybins[i-numptbins];
					outbin = i-numptbins;
				}
				else if (i == lastBin+1) //-2.87--2.4 cross
				{
					ptLow = 0;
					ptHigh = 30;
					yLow = -2.87;
					yHigh = -1.93;
					if (iState == 1) hThisPA = hy1S_cr_PADiff;
					else if (iState == 2) hThisPA = hy2S_cr_PADiff;
					else if (iState == 3) hThisPA = hy3S_cr_PADiff;
					outbin = 1;
				}
				
				int lineNumber = i;
				if (diffBin == 0)
				{
					if (i > numptbins && i <= lastBin) lineNumber++; //Due to omitted -2.4,-1.93 bin
					if (iState > 1) lineNumber += 16;
					if (iState > 2) lineNumber += 9;
				}
				else if (diffBin > 0)
				{
					lineNumber += 30;
					if (iState > 1) lineNumber += 14;
					if (iState > 2) lineNumber += 7;
					if (diffBin > 1)
						lineNumber += 25;
				}
				if (i == lastBin+1)
					lineNumber = 136 + iState; //cr bins: 1s=137, 2s=138, 3s=139
				
				cout << Form("diffbin = %d, state = %d, i = %d, Bin (%.1f<pt<%.1f, %.2f<y<%.2f)",diffBin,iState,i,ptLow,ptHigh,yLow,yHigh) << endl;
				
				bool hasPP = true;
				
				TString inFileNamePA = Form("ResultsBkg/PseudoExpResults_PA_pt%.1f-%.1f_y%.2f-%.2f.root",ptLow,ptHigh,yLow,yHigh);
				TFile* inFilePA = new TFile(inFileNamePA,"READ");
				if (inFilePA->IsZombie())
				{
					cout << "LINE " << lineNumber << " HAS NO PA FILE" << endl;
					lineValues[lineNumber] = "& * & * & *";
					continue;
				}
				float yLowPP = yLow;
				float yHighPP = yHigh;
				if (yLow < 0.0 && yHigh > 0.0)
					yLowPP = 0.0;
				else if (yLow < 0.0)
				{
					yLowPP = TMath::Abs(yHigh);
					yHighPP = TMath::Abs(yLow);
				}
				TString inFileNamePP = Form("ResultsBkg/PseudoExpResults_PP_pt%.1f-%.1f_y%.2f-%.2f.root",ptLow,ptHigh,yLowPP,yHighPP);
				TFile* inFilePP = new TFile(inFileNamePP,"READ");
				if (inFilePP->IsZombie())
				{
					cout << "LINE " << lineNumber << " HAS NO PP FILE" << endl;
					hasPP = false;
					//lineValues[lineNumber] = "& * & * & *";
					//continue; //Skip for now
				}
				
				/*TH1F* histoPA;
				TH1F* histoPP;
				TH1F* histoRpA;*/
				TH1F* histoPA = new TH1F("histoPA",Form("%dS %%Diff in Yields for PA",iState),100,0,20);
				TH1F* histoPP = new TH1F("histoPP",Form("%dS %%Diff in Yields for PP",iState),100,0,40);
				TH1F* histoRpA = new TH1F("histoRpA",Form("%dS %%Diff in Yields for R_{pA}",iState),100,0,100);
				histoPA->SetXTitle("%Diff");
				histoPP->SetXTitle("%Diff");
				histoRpA->SetXTitle("%Diff");
				
				TNtuple* tupPA = (TNtuple*)(inFilePA->Get("ntupleSig"));
				TNtuple* tupPP = (TNtuple*)(inFilePP->Get("ntupleSig"));
				
				double diffs[3] = {0,0,0};
				if (hasPP)
					calcError(tupPA,tupPP,diffs,iState,histoPA,histoPP,histoRpA);
				else
					calcError(tupPA,diffs,iState,histoPA);
				
				
				TCanvas* c1 = new TCanvas("cHists","%Diff histograms",4,45,1100,400);
				gStyle->SetOptStat("nemr");
				c1->Divide(3,1);
				c1->cd(1);
				histoPA->Draw();
				if (hasPP)
				{
					c1->cd(2);
					histoPP->Draw();
					c1->cd(3);
					histoRpA->Draw();
				}
				
				TPaveText* binBox[3];
				TText* binTextPt[3];
				TText* binTextRap[3];
				for (int cdi=0; cdi<3; cdi++)
				{
					c1->cd(cdi+1);
					if (hasPP || cdi==0)
					{
						binBox[cdi] = new TPaveText(0.1,0.8,0.35,0.9,"NDC");
						binBox[cdi]->SetFillColor(kWhite);
						binBox[cdi]->SetBorderSize(1);
						binTextPt[cdi] = binBox[cdi]->AddText(Form("%.0f<p_{T}<%.0f",ptLow,ptHigh));
						binTextPt[cdi]->SetTextAlign(12);
						binTextRap[cdi] = binBox[cdi]->AddText(Form("%.2f<y_{cm}<%.2f",yLow,yHigh));
						binTextRap[cdi]->SetTextAlign(12);
					}
					else
					{
						binBox[cdi] = new TPaveText(0.3,0.45,0.7,0.55,"NDC");
						if (cdi==1)
							binBox[cdi]->AddText("No PP data for this bin.");
						else if (cdi==2)
							binBox[cdi]->AddText("No PP data for this bin.");
					}
					binBox[cdi]->Draw();
				}
				
				TString outFileName = Form("PseudoExpResultHists_pt%.1f-%.1f_y%.2f-%.2f_%dS",ptLow,ptHigh,yLow,yHigh,iState);
				c1->SaveAs(Form("ResultsBkg/plots/") + outFileName + ".pdf");
				c1->SaveAs(Form("ResultsBkg/plots/pngs/") + outFileName + ".png");
				
				
				float PAdiff = diffs[0];
				float PPdiff = 0;
				float RpAdiff = 0;
				if (hasPP)
				{
					PPdiff = diffs[1];
					RpAdiff = diffs[2];
				}
				
				resultTuple->Fill(ptLow,ptHigh,yLow,yHigh,PPdiff,PAdiff,RpAdiff);
				
				hThisPA->SetBinContent(outbin,PAdiff);
				if (hasPP)
				{
				hThisPP->SetBinContent(outbin,PPdiff);
				hThisRpA->SetBinContent(outbin,RpAdiff);
				}
				
				//number of digits after decimal to display
				int sfad_pp = 2 - TMath::Ceil(TMath::Log10(PPdiff*100)); if (PPdiff==0) sfad_pp = 1;
				int sfad_pa = 2 - TMath::Ceil(TMath::Log10(PAdiff*100)); if (PAdiff==0) sfad_pa = 1;
				int sfad_rpa = 2 - TMath::Ceil(TMath::Log10(RpAdiff*100)); if (RpAdiff==0) sfad_rpa = 1;
				
				cout << "Setting line " << lineNumber << " to " << Form("& %.2f & %.2f & %.2f",PPdiff,PAdiff,RpAdiff) << endl;
				if (hasPP)
					lineValues[lineNumber] = Form(Form("& %%.%df & %%.%df & %%.%df",sfad_pp,sfad_pa,sfad_rpa),PPdiff*100,PAdiff*100,RpAdiff*100);
					//lineValues[lineNumber] = Form("& %.2f & %.2f & %.2f",PPdiff*100,PAdiff*100,RpAdiff*100);
				else
					lineValues[lineNumber] = Form(Form("& & %%.%df & ",sfad_pa),PAdiff*100);
					//lineValues[lineNumber] = Form("& & %.2f & ",PAdiff*100);
				/*
				outNames[lineNumber] = outFileName;
				figCaption[lineNumber] = Form("Percent differences in $\\Upsilon$ %dS yields and $R_{pA}$ between alternative and nominal background PDFs in fits to pseudodata for $\\pt$ $\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f].",iState,ptLow,ptHigh,yLow,yHigh);
				figureLabel[lineNumber] = Form("BkgPseudoExptResult_%dS_pt%.0fto%.0f_y%.0fto%.0f",iState,ptLow,ptHigh,yLow*100,yHigh*100);
				
				if (!hasPP)
					figCaption[lineNumber] = Form("Percent differences in $\\Upsilon$ %dS yields of alternative and nominal background PDFs in fits to pseudodata for $\\pt$ $\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f].",iState,ptLow,ptHigh,yLow,yHigh);
				
				for (int cdi=0; cdi<3; cdi++)
					delete binBox[cdi];
				*/
			}
		}
	}
	
	//Copy ybins into cross section, except 1st bin
	for (int i=1; i<9; i++)
		hy1S_cr_PADiff->SetBinContent(i+1,hy1S_PADiff->GetBinContent(i));
	for (int i=1; i<5; i++)
		hy2S_cr_PADiff->SetBinContent(i+1,hy2S_PADiff->GetBinContent(i));
	for (int i=1; i<3; i++)
		hy3S_cr_PADiff->SetBinContent(i+1,hy3S_PADiff->GetBinContent(i));
	
	//Cross section pt bins
	for (int iState=1; iState<=3; iState++)
	{
		if (iState == 1)
		{
			ptbins = ptbins1s;
			numptbins = 6;
			hThisPA = hpt1S_cr_PADiff;
		}
		else if (iState == 2)
		{
			ptbins = ptbins2s;
			numptbins = 3;
			hThisPA = hpt2S_cr_PADiff;
		}
		else if (iState == 3)
		{
			ptbins = ptbins3s;
			numptbins = 2;
			hThisPA = hpt3S_cr_PADiff;
		}
		else
		{
			cout << "Invalid state specified. Aborting" << endl;
			return;
		}
		for (int i=0; i<numptbins; i++)
		{
			float yLow = -2.87;
			float yHigh = 1.93;
			float ptLow = ptbins[i];
			float ptHigh = ptbins[i+1];
			
			int lineNumber = 140+i;
			if (iState > 1)
				lineNumber += 6;
			if (iState > 2)
				lineNumber += 3;
				
			outbin = i+1;
			
			cout << Form("state = %d, i = %d, Bin (%.1f<pt<%.1f, %.2f<y<%.2f)",iState,i,ptLow,ptHigh,yLow,yHigh) << endl;
			
			bool hasPP = true;
			
			TString inFileNamePA = Form("ResultsBkg/PseudoExpResults_PA_pt%.1f-%.1f_y%.2f-%.2f.root",ptLow,ptHigh,yLow,yHigh);
			TFile* inFilePA = new TFile(inFileNamePA,"READ");
			if (inFilePA->IsZombie())
			{
				cout << "LINE " << lineNumber << " HAS NO PA FILE" << endl;
				lineValues[lineNumber] = "& * & * & *";
				continue;
			}
			
			//TH1F* histoPA;
			TH1F* histoPA = new TH1F("histoPA","histoPA",100,0,100);
			
			TNtuple* tupPA = (TNtuple*)(inFilePA->Get("ntupleSig"));
			
			double diffs[3] = {0,0,0};
			calcError(tupPA,diffs,iState,histoPA);
			
			
			TCanvas* c1 = new TCanvas("cHists","%Diff histograms",4,45,1100,400);
			gStyle->SetOptStat("nemr");
			c1->Divide(3,1);
			c1->cd(1);
			histoPA->Draw();
			
			TPaveText* binBox[3];
			TText* binTextPt[3];
			TText* binTextRap[3];
			for (int cdi=0; cdi<3; cdi++)
			{
				c1->cd(cdi+1);
				if (cdi==0)
				{
					binBox[cdi] = new TPaveText(0.1,0.8,0.35,0.9,"NDC");
					binBox[cdi]->SetFillColor(kWhite);
					binBox[cdi]->SetBorderSize(1);
					binTextPt[cdi] = binBox[cdi]->AddText(Form("%.0f<p_{T}<%.0f",ptLow,ptHigh));
					binTextPt[cdi]->SetTextAlign(12);
					binTextRap[cdi] = binBox[cdi]->AddText(Form("%.2f<y_{cm}<%.2f",yLow,yHigh));
					binTextRap[cdi]->SetTextAlign(12);
				}
				else
				{
					binBox[cdi] = new TPaveText(0.3,0.45,0.7,0.55,"NDC");
					if (cdi==1)
						binBox[cdi]->AddText("No PP data for this bin.");
					else if (cdi==2)
						binBox[cdi]->AddText("No PP data for this bin.");
				}
				binBox[cdi]->Draw();
			}
			
			TString outFileName = Form("PseudoExpResultHists_pt%.1f-%.1f_y%.2f-%.2f_%dS",ptLow,ptHigh,yLow,yHigh,iState);
			c1->SaveAs(Form("ResultsBkg/plots/") + outFileName + ".pdf");
			c1->SaveAs(Form("ResultsBkg/plots/pngs/") + outFileName + ".png");
			
			
			float PAdiff = diffs[0];
			float PPdiff = 0;
			float RpAdiff = 0;
			
			resultTuple->Fill(ptLow,ptHigh,yLow,yHigh,PPdiff,PAdiff,RpAdiff);
			
			hThisPA->SetBinContent(outbin,PAdiff);
			
			//number of digits after decimal to display
			int sfad_pp = 2 - TMath::Ceil(TMath::Log10(PPdiff*100)); if (PPdiff==0) sfad_pp = 1;
			int sfad_pa = 2 - TMath::Ceil(TMath::Log10(PAdiff*100)); if (PAdiff==0) sfad_pa = 1;
			int sfad_rpa = 2 - TMath::Ceil(TMath::Log10(RpAdiff*100)); if (RpAdiff==0) sfad_rpa = 1;
			
			cout << "Setting line " << lineNumber << " to " << Form("& %.2f & %.2f & %.2f",PPdiff,PAdiff,RpAdiff) << endl;
			lineValues[lineNumber] = Form(Form("& & %%.%df & ",sfad_pa),PAdiff*100);
			//lineValues[lineNumber] = Form("& & %.2f & ",PAdiff*100);
			
			/*
			outNames[lineNumber] = outFileName;
			figCaption[lineNumber] = Form("Percent differences in $\\Upsilon$ %dS yields of alternative and nominal background PDFs in fits to pseudodata for $\\pt$ $\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f].",iState,ptLow,ptHigh,yLow,yHigh);
			figureLabel[lineNumber] = Form("BkgPseudoExptResult_%dS_pt%.0fto%.0f_y%.0fto%.0f",iState,ptLow,ptHigh,yLow*100,yHigh*100);
			
			for (int cdi=0; cdi<3; cdi++)
				delete binBox[cdi];
			*/
		}
	}
	
	//Rfb
	float yLow_rfbBins[7] = {0,0.4,0.8,1.2,0,0.8,0};
	float yHigh_rfbBins[7] = {0.4,0.8,1.2,1.93,0.8,1.93,1.93};
	/*
	int ntLow_rfbBins[8] = {0,40,65,90,0,0,0,0};
	int ntHigh_rfbBins[8] = {40,65,90,400,400,400,400,400};
	float hfLow_rfbBins[8] = {0,0,0,0,0,15,22,30};
	float hfHigh_rfbBins[8] = {400,400,400,400,15,22,30,120};
	*/
	int ntLow_rfbBins[8] = {0,40,62,88,0,0,0,0};
	int ntHigh_rfbBins[8] = {40,62,88,400,400,400,400,400};
	float hfLow_rfbBins[8] = {0,0,0,0,0,12,19,27};
	float hfHigh_rfbBins[8] = {120,120,120,120,12,19,27,120};
	TH2F* hThisPL;
	TH2F* hThisMI;
	TH2F* hThisRfb;
	TH1F* hThisPL1;
	TH1F* hThisMI1;
	TH1F* hThisRfb1;
	for (int i=0; i<8; i++)
	{
		int outbinnthf = i+1;
		if (i>3)
			outbinnthf = i-3;
		
		int ntLow = ntLow_rfbBins[i];
		int ntHigh = ntHigh_rfbBins[i];
		float hfLow = hfLow_rfbBins[i];
		float hfHigh = hfHigh_rfbBins[i];
		for (int j=0; j<10; j++)
		{
			float yLow;
			float yHigh;
			if (j>=7)
			{
				yLow = yLow_rfbBins[6];
				yHigh = yHigh_rfbBins[6];
			}
			else
			{
				yLow = yLow_rfbBins[j];
				yHigh = yHigh_rfbBins[j];
			}
			float ptLow = 0;
			float ptHigh = 30;
			
			int iState = 1;
			
			int outbiny = j+1;
			if (i>3) {hThisPL = hhf1S_yPLdiff; hThisMI = hhf1S_yMIdiff; hThisRfb = hhf1S_yRfbdiff;}
			else {hThisPL = hnt1S_yPLdiff; hThisMI = hnt1S_yMIdiff; hThisRfb = hnt1S_yRfbdiff;}
			if (j>3)
			{
				outbiny = j-3;
				iState = 2;
				if (i>3) {hThisPL = hhf2S_yPLdiff; hThisMI = hhf2S_yMIdiff; hThisRfb = hhf2S_yRfbdiff;}
				else {hThisPL = hnt2S_yPLdiff; hThisMI = hnt2S_yMIdiff; hThisRfb = hnt2S_yRfbdiff;}
			}
			if (j==6)
			{
				outbiny = 1;
				iState = 3;
				if (i>3) {hThisPL = hhf3S_yPLdiff; hThisMI = hhf3S_yMIdiff; hThisRfb = hhf3S_yRfbdiff;}
				else {hThisPL = hnt3S_yPLdiff; hThisMI = hnt3S_yMIdiff; hThisRfb = hnt3S_yRfbdiff;}
			}
			if (j>6)
			{
				iState = j-6;
				if (iState==1)
				{
					if (i>3) {hThisPL1 = hhf1S_yInt_PLdiff; hThisMI1 = hhf1S_yInt_MIdiff; hThisRfb1 = hhf1S_yInt_Rfbdiff;}
					else {hThisPL1 = hnt1S_yInt_PLdiff; hThisMI1 = hnt1S_yInt_MIdiff; hThisRfb1 = hnt1S_yInt_Rfbdiff;}
				}
				if (iState==2)
				{
					if (i>3) {hThisPL1 = hhf2S_yInt_PLdiff; hThisMI1 = hhf2S_yInt_MIdiff; hThisRfb1 = hhf2S_yInt_Rfbdiff;}
					else {hThisPL1 = hnt2S_yInt_PLdiff; hThisMI1 = hnt2S_yInt_MIdiff; hThisRfb1 = hnt2S_yInt_Rfbdiff;}
				}
				if (iState==3)
				{
					if (i>3) {hThisPL1 = hhf3S_yInt_PLdiff; hThisMI1 = hhf3S_yInt_MIdiff; hThisRfb1 = hhf3S_yInt_Rfbdiff;}
					else {hThisPL1 = hnt3S_yInt_PLdiff; hThisMI1 = hnt3S_yInt_MIdiff; hThisRfb1 = hnt3S_yInt_Rfbdiff;}
					
					if (i==1) ntHigh=400;
					else if (i==5) hfHigh=120;
					else if (i==2 || i==3 || i>5) continue;
				}
			}
			
			int lineNumber = 81+i*7+j;
			if (j>6)
			{
				lineNumber = 151+(iState-1)*4;
				if (i<=3)
					lineNumber += i;
				else
					lineNumber += i-4+12;
			}
			cout << Form("State = %d, i = %d, j=%d, Bin (%.2f<y<%.2f, %d<Ntracks<%d, %.0f<HF<%.0f)",iState,i,j,yLow,yHigh,ntLow,ntHigh,hfLow,hfHigh) << endl;
			
			TString inFileName = Form("ResultsBkg/PseudoExpResults_pt%.1f-%.1f_y%.2f-%.2f_hfsum%.2f-%.2f_ntracks%d-%d.root",ptLow,ptHigh,yLow,yHigh,hfLow,hfHigh,ntLow,ntHigh);
			TFile* inFile = new TFile(inFileName,"READ");
			if (inFile->IsZombie())
			{
				cout << "LINE " << lineNumber << " HAS NO FILE" << endl;
				lineValues[lineNumber] = "& * & * & *";
				continue;
			}
			
			TH1F* histoPL = (TH1F*)inFile->Get(Form("PL_%dSDiff",iState));
			TH1F* histoMI = (TH1F*)inFile->Get(Form("MI_%dSDiff",iState));
			TH1F* histoRfb = (TH1F*)inFile->Get(Form("Rfb_%dSDiff",iState));
			if (iState == 1) //Argh, typo in pseudoexperiments naming
				histoRfb = (TH1F*)inFile->Get(Form("RPL_%dSDiff",iState));
			
			TCanvas* c1 = new TCanvas("cHists","%Diff histograms",4,45,1100,400);
			//gStyle->SetOptFit(1);
			gStyle->SetOptStat("nemr");
			c1->Divide(3,1);
			c1->cd(1);
			histoPL->Draw();
			c1->cd(2);
			histoMI->Draw();
			c1->cd(3);
			histoRfb->Draw();
			
			TPaveText* binBox[3];
			TText* binTextPt[3];
			TText* binTextRap[3];
			TText* binTextNthf[3];
			for (int cdi=0; cdi<3; cdi++)
			{
				c1->cd(cdi+1);
				binBox[cdi] = new TPaveText(0.1,0.75,0.35,0.9,"NDC");
				binBox[cdi]->SetFillColor(kWhite);
				binBox[cdi]->SetBorderSize(1);
				binTextPt[cdi] = binBox[cdi]->AddText(Form("%.0f<p_{T}<%.0f",ptLow,ptHigh));
				binTextPt[cdi]->SetTextAlign(12);
				binTextRap[cdi] = binBox[cdi]->AddText(Form("%.2f<y_{cm}<%.2f",yLow,yHigh));
				binTextRap[cdi]->SetTextAlign(12);
				if (i<=3)
					binTextNthf[cdi] = binBox[cdi]->AddText(Form("%d<NTracks<%d",ntLow,ntHigh));
				else
					binTextNthf[cdi] = binBox[cdi]->AddText(Form("%.0f<HF<%.0f",hfLow,hfHigh));
				binTextNthf[cdi]->SetTextAlign(12);
				binBox[cdi]->Draw();
			}
			
			TString outFileName = Form("PseudoExpResultHists_pt%.1f-%.1f_y%.2f-%.2f_hfsum%.2f-%.2f_ntracks%d-%d_%dS",ptLow,ptHigh,yLow,yHigh,hfLow,hfHigh,ntLow,ntHigh,iState);
			c1->SaveAs(Form("ResultsBkg/plots/") + outFileName + ".pdf");
			c1->SaveAs(Form("ResultsBkg/plots/pngs/") + outFileName + ".png");
			
			TNtuple* tuple = (TNtuple*)(inFile->Get("ntupleDiff"));
			float diffPL, diffMI, diffRfb;
			tuple->SetBranchAddress(Form("diff%dsPL",iState),&diffPL);
			tuple->SetBranchAddress(Form("diff%dsMI",iState),&diffMI);
			tuple->SetBranchAddress(Form("diff%dsRfb",iState),&diffRfb);
			int ntEntries = tuple->GetEntries();
			int denom = ntEntries;
			
			float sumPL = 0;
			float sumMI = 0;
			float sumRfb = 0;
			float hiPL = 0; float hiMI = 0; float hiRfb = 0;
			float loPL = 0; float loMI = 0; float loRfb = 0;
			
			for (int ent=0; ent<ntEntries; ent++)
			{
				tuple->GetEntry(ent);
				if (diffPL>1000 || diffMI>1000 || diffRfb>1000) {denom--; continue;} //skip far outliers
				sumPL += TMath::Abs(diffPL);
				sumMI += TMath::Abs(diffMI);
				sumRfb += TMath::Abs(diffRfb);
				if (diffPL > hiPL) hiPL = diffPL; if (diffPL < loPL) loPL = diffPL;
				if (diffMI > hiMI) hiMI = diffMI; if (diffMI < loMI) loMI = diffMI;
				if (diffRfb > hiRfb) hiRfb = diffRfb; if (diffRfb < loRfb) loRfb = diffRfb;
			}
			float rmsPL = 0; float rmsMI = 0; float rmsRfb = 0;
			for (int ent=0; ent<ntEntries; ent++)
			{
				tuple->GetEntry(ent);
				if (diffPL>1000 || diffMI>1000 || diffRfb>1000) {continue;} //skip far outliers
				rmsPL += TMath::Power(TMath::Abs(diffPL)-sumPL/denom,2);
				rmsMI += TMath::Power(TMath::Abs(diffMI)-sumMI/denom,2);
				rmsRfb += TMath::Power(TMath::Abs(diffRfb)-sumRfb/denom,2);
			}
			rmsPL = TMath::Sqrt(rmsPL/denom); rmsMI = TMath::Sqrt(rmsMI/denom); rmsRfb = TMath::Sqrt(rmsRfb/denom);

			float PLdiff = TMath::Abs(0.01*sumPL/denom);
			float MIdiff = TMath::Abs(0.01*sumMI/denom);
			float Rfbdiff = TMath::Abs(0.01*sumRfb/denom);
			
			/*float PLdiff = TMath::Max(TMath::Abs(histoPL->GetMean()),TMath::Abs(histoPL->GetRMS()))*0.01;
			float MIdiff = TMath::Max(TMath::Abs(histoMI->GetMean()),TMath::Abs(histoMI->GetRMS()))*0.01;
			float Rfbdiff = TMath::Max(TMath::Abs(histoRfb->GetMean()),TMath::Abs(histoRfb->GetRMS()))*0.01;*/
			
			resultTuple->Fill(ptLow,ptHigh,yLow,yHigh,MIdiff,PLdiff,Rfbdiff);
			
			if (j>6)
			{
				hThisPL1->SetBinContent(outbinnthf,PLdiff);
				hThisMI1->SetBinContent(outbinnthf,MIdiff);
				hThisRfb1->SetBinContent(outbinnthf,Rfbdiff);
			}
			else
			{
				hThisPL->SetBinContent(outbiny,outbinnthf,PLdiff);
				hThisMI->SetBinContent(outbiny,outbinnthf,MIdiff);
				hThisRfb->SetBinContent(outbiny,outbinnthf,Rfbdiff);
			}
			
			//number of digits after decimal to display
			int sfad_mi = 2 - TMath::Ceil(TMath::Log10(MIdiff*100)); if (MIdiff==0) sfad_mi = 1;
			int sfad_pl = 2 - TMath::Ceil(TMath::Log10(PLdiff*100)); if (PLdiff==0) sfad_pl = 1;
			int sfad_rfb = 2 - TMath::Ceil(TMath::Log10(Rfbdiff*100)); if (Rfbdiff==0) sfad_rfb = 1;
			
			cout << "Setting line " << lineNumber << " to " << Form("& %.2f & %.2f & %.2f",MIdiff,PLdiff,Rfbdiff) << endl;
			lineValues[lineNumber] = Form(Form("& %%.%df & %%.%df & %%.%df",sfad_mi,sfad_pl,sfad_rfb),MIdiff*100,PLdiff*100,Rfbdiff*100);
			
			outNames[lineNumber] = outFileName;
			if (i<4)
			{
				figCaption[lineNumber] = Form("Percent differences in $\\Upsilon$ %dS yields and $R_{fb}$ between alternative and nominal background PDFs in fits to pseudodata for $\\pt$ $\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f], NTracks $\\in$ [%d,%d].",iState,ptLow,ptHigh,yLow,yHigh,ntLow,ntHigh);
				figureLabel[lineNumber] = Form("BkgPseudoExptResult_%dS_pt%.0fto%.0f_y%.0fto%.0f_nt%dto%d",iState,ptLow,ptHigh,yLow*100,yHigh*100,ntLow,ntHigh);
			}
			else
			{
				figCaption[lineNumber] = Form("Percent differences in $\\Upsilon$ %dS yields and $R_{fb}$ between alternative and nominal background PDFs in fits to pseudodata for $\\pt$ $\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f], HF $\\in$ [%.0f,%.0f].",iState,ptLow,ptHigh,yLow,yHigh,hfLow,hfHigh);
				figureLabel[lineNumber] = Form("BkgPseudoExptResult_%dS_pt%.0fto%.0f_y%.0fto%.0f_hf%.0fto%.0f",iState,ptLow,ptHigh,yLow*100,yHigh*100,hfLow,hfHigh);
			}
			
			//clean up
			for (int cdi=0; cdi<3; cdi++)
					delete binBox[cdi];
		}
	}
	
	//Integrated Rfb
	for (int iState = 1; iState <= 2; iState++)
	{
		int lineNumber = 174 + iState;
		
		TFile* inFilePL = new TFile("ResultsBkg/PseudoExpResults_PA_pt0.0-30.0_y0.00-1.93.root","READ");
		if (inFilePL->IsZombie())
		{
			cout << "INTEGRATED RFB HAS NO PL FILE" << endl;
			continue;
		}
		TFile* inFileMI = new TFile("ResultsBkg/PseudoExpResults_PA_pt0.0-30.0_y-1.93-0.00.root","READ");
		if (inFileMI->IsZombie())
		{
			cout << "INTEGRATED RFB HAS NO MI FILE" << endl;
			continue;
		}
		
		if (iState==1) {hThisPL1 = hInt1S_yPLdiff; hThisMI1 = hInt1S_yMIdiff; hThisRfb1 = hInt1S_Rfbdiff;}
		else if (iState==2) {hThisPL1 = hInt2S_yPLdiff; hThisMI1 = hInt2S_yMIdiff; hThisRfb1 = hInt2S_Rfbdiff;}
		
		TNtuple* tupPL = (TNtuple*)(inFilePL->Get("ntupleSig"));
		TNtuple* tupMI = (TNtuple*)(inFileMI->Get("ntupleSig"));
		
		TH1F* histoPL = new TH1F("histoPL","histoPL",100,0,100);
		TH1F* histoMI = new TH1F("histoMI","histoMI",100,0,100);
		TH1F* histoRfb = new TH1F("histoRfb","histoRfb",100,0,100);
		
		double diffs[3] = {0,0,0};
		calcError(tupPL,tupMI,diffs,iState,histoPL,histoMI,histoRfb);

		float PLdiff = diffs[0];
		float MIdiff = diffs[1];
		float Rfbdiff = diffs[2];
		
		/*float PLdiff = TMath::Max(TMath::Abs(histoPL->GetMean()),TMath::Abs(histoPL->GetRMS()))*0.01;
		float MIdiff = TMath::Max(TMath::Abs(histoMI->GetMean()),TMath::Abs(histoMI->GetRMS()))*0.01;
		float Rfbdiff = TMath::Max(TMath::Abs(histoRfb->GetMean()),TMath::Abs(histoRfb->GetRMS()))*0.01;*/
		
		hThisPL1->SetBinContent(1,PLdiff);
		hThisMI1->SetBinContent(1,MIdiff);
		hThisRfb1->SetBinContent(1,Rfbdiff);
		
		int sfad_mi = 2 - TMath::Ceil(TMath::Log10(MIdiff*100)); if (MIdiff==0) sfad_mi = 1;
		int sfad_pl = 2 - TMath::Ceil(TMath::Log10(PLdiff*100)); if (PLdiff==0) sfad_pl = 1;
		int sfad_rfb = 2 - TMath::Ceil(TMath::Log10(Rfbdiff*100)); if (Rfbdiff==0) sfad_rfb = 1;
		
		cout << "Setting line " << lineNumber << " to " << Form("& %.2f & %.2f & %.2f",MIdiff,PLdiff,Rfbdiff) << endl;
		lineValues[lineNumber] = Form(Form("& %%.%df & %%.%df & %%.%df",sfad_mi,sfad_pl,sfad_rfb),MIdiff*100,PLdiff*100,Rfbdiff*100);
	}
	
	//Make 1D hists for nt, hf
	TH1F* hnt1S_y1_yPLdiff = new TH1F("hnt1S_y1_yPLdiff","Ntrack bins, y0.00to0.40",4,ntbins);
	TH1F* hnt2S_y1_yPLdiff = new TH1F("hnt2S_y1_yPLdiff","Ntrack bins, y0.00to0.80",4,ntbins);
	TH1F* hnt3S_y1_yPLdiff = new TH1F("hnt3S_y1_yPLdiff","Ntrack bins, y0.00to1.93",4,ntbins);
	TH1F* hnt1S_y1_yMIdiff = new TH1F("hnt1S_y1_yMIdiff","Ntrack bins, y0.00to0.40",4,ntbins);
	TH1F* hnt2S_y1_yMIdiff = new TH1F("hnt2S_y1_yMIdiff","Ntrack bins, y0.00to0.80",4,ntbins);
	TH1F* hnt3S_y1_yMIdiff = new TH1F("hnt3S_y1_yMIdiff","Ntrack bins, y0.00to1.93",4,ntbins);
	TH1F* hnt1S_y1_yRfbdiff = new TH1F("hnt1S_y1_yRfbdiff","Ntrack bins, y0.00to0.40",4,ntbins);
	TH1F* hnt2S_y1_yRfbdiff = new TH1F("hnt2S_y1_yRfbdiff","Ntrack bins, y0.00to0.80",4,ntbins);
	TH1F* hnt3S_y1_yRfbdiff = new TH1F("hnt3S_y1_yRfbdiff","Ntrack bins, y0.00to1.93",4,ntbins);
	TH1F* hhf1S_y1_yPLdiff = new TH1F("hhf1S_y1_yPLdiff","HF bins, y0.00to0.40",4,hfbins);
	TH1F* hhf2S_y1_yPLdiff = new TH1F("hhf2S_y1_yPLdiff","HF bins, y0.00to0.80",4,hfbins);
	TH1F* hhf3S_y1_yPLdiff = new TH1F("hhf3S_y1_yPLdiff","HF bins, y0.00to1.93",4,hfbins);
	TH1F* hhf1S_y1_yMIdiff = new TH1F("hhf1S_y1_yMIdiff","HF bins, y0.00to0.40",4,hfbins);
	TH1F* hhf2S_y1_yMIdiff = new TH1F("hhf2S_y1_yMIdiff","HF bins, y0.00to0.80",4,hfbins);
	TH1F* hhf3S_y1_yMIdiff = new TH1F("hhf3S_y1_yMIdiff","HF bins, y0.00to1.93",4,hfbins);
	TH1F* hhf1S_y1_yRfbdiff = new TH1F("hhf1S_y1_yRfbdiff","HF bins, y0.00to0.40",4,hfbins);
	TH1F* hhf2S_y1_yRfbdiff = new TH1F("hhf2S_y1_yRfbdiff","HF bins, y0.00to0.80",4,hfbins);
	TH1F* hhf3S_y1_yRfbdiff = new TH1F("hhf3S_y1_yRfbdiff","HF bins, y0.00to1.93",4,hfbins);
	for (int i=1; i<=4; i++)
	{
		hnt1S_y1_yPLdiff->SetBinContent(i,hnt1S_yPLdiff->GetBinContent(1,i));
		hnt2S_y1_yPLdiff->SetBinContent(i,hnt2S_yPLdiff->GetBinContent(1,i));
		hnt3S_y1_yPLdiff->SetBinContent(i,hnt3S_yPLdiff->GetBinContent(1,i));
		hnt1S_y1_yMIdiff->SetBinContent(i,hnt1S_yMIdiff->GetBinContent(1,i));
		hnt2S_y1_yMIdiff->SetBinContent(i,hnt2S_yMIdiff->GetBinContent(1,i));
		hnt3S_y1_yMIdiff->SetBinContent(i,hnt3S_yMIdiff->GetBinContent(1,i));
		hnt1S_y1_yRfbdiff->SetBinContent(i,hnt1S_yRfbdiff->GetBinContent(1,i));
		hnt2S_y1_yRfbdiff->SetBinContent(i,hnt2S_yRfbdiff->GetBinContent(1,i));
		hnt3S_y1_yRfbdiff->SetBinContent(i,hnt3S_yRfbdiff->GetBinContent(1,i));
		hhf1S_y1_yPLdiff->SetBinContent(i,hhf1S_yPLdiff->GetBinContent(1,i));
		hhf2S_y1_yPLdiff->SetBinContent(i,hhf2S_yPLdiff->GetBinContent(1,i));
		hhf3S_y1_yPLdiff->SetBinContent(i,hhf3S_yPLdiff->GetBinContent(1,i));
		hhf1S_y1_yMIdiff->SetBinContent(i,hhf1S_yMIdiff->GetBinContent(1,i));
		hhf2S_y1_yMIdiff->SetBinContent(i,hhf2S_yMIdiff->GetBinContent(1,i));
		hhf3S_y1_yMIdiff->SetBinContent(i,hhf3S_yMIdiff->GetBinContent(1,i));
		hhf1S_y1_yRfbdiff->SetBinContent(i,hhf1S_yRfbdiff->GetBinContent(1,i));
		hhf2S_y1_yRfbdiff->SetBinContent(i,hhf2S_yRfbdiff->GetBinContent(1,i));
		hhf3S_y1_yRfbdiff->SetBinContent(i,hhf3S_yRfbdiff->GetBinContent(1,i));
	}
	TH1F* hnt1S_y2_yPLdiff = new TH1F("hnt1S_y2_yPLdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hnt2S_y2_yPLdiff = new TH1F("hnt2S_y2_yPLdiff","Ntrack bins, y0.80to1.93",4,ntbins);
	TH1F* hnt3S_y2_yPLdiff = new TH1F("hnt3S_y2_yPLdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hnt1S_y2_yMIdiff = new TH1F("hnt1S_y2_yMIdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hnt2S_y2_yMIdiff = new TH1F("hnt2S_y2_yMIdiff","Ntrack bins, y0.80to1.93",4,ntbins);
	TH1F* hnt3S_y2_yMIdiff = new TH1F("hnt3S_y2_yMIdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hnt1S_y2_yRfbdiff = new TH1F("hnt1S_y2_yRfbdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hnt2S_y2_yRfbdiff = new TH1F("hnt2S_y2_yRfbdiff","Ntrack bins, y0.80to1.93",4,ntbins);
	TH1F* hnt3S_y2_yRfbdiff = new TH1F("hnt3S_y2_yRfbdiff","Ntrack bins, y0.40to0.80",4,ntbins);
	TH1F* hhf1S_y2_yPLdiff = new TH1F("hhf1S_y2_yPLdiff","HF bins, y0.40to0.80",4,hfbins);
	TH1F* hhf2S_y2_yPLdiff = new TH1F("hhf2S_y2_yPLdiff","HF bins, y0.80to1.93",4,hfbins);
	TH1F* hhf3S_y2_yPLdiff = new TH1F("hhf3S_y2_yPLdiff","HF bins, y0.40to0.80",4,hfbins);
	TH1F* hhf1S_y2_yMIdiff = new TH1F("hhf1S_y2_yMIdiff","HF bins, y0.40to0.80",4,hfbins);
	TH1F* hhf2S_y2_yMIdiff = new TH1F("hhf2S_y2_yMIdiff","HF bins, y0.80to1.93",4,hfbins);
	TH1F* hhf3S_y2_yMIdiff = new TH1F("hhf3S_y2_yMIdiff","HF bins, y0.40to0.80",4,hfbins);
	TH1F* hhf1S_y2_yRfbdiff = new TH1F("hhf1S_y2_yRfbdiff","HF bins, y0.40to0.80",4,hfbins);
	TH1F* hhf2S_y2_yRfbdiff = new TH1F("hhf2S_y2_yRfbdiff","HF bins, y0.80to1.93",4,hfbins);
	TH1F* hhf3S_y2_yRfbdiff = new TH1F("hhf3S_y2_yRfbdiff","HF bins, y0.40to0.80",4,hfbins);
	for (int i=1; i<=4; i++)
	{
		hnt1S_y2_yPLdiff->SetBinContent(i,hnt1S_yPLdiff->GetBinContent(2,i));
		hnt2S_y2_yPLdiff->SetBinContent(i,hnt2S_yPLdiff->GetBinContent(2,i));
		hnt3S_y2_yPLdiff->SetBinContent(i,hnt3S_yPLdiff->GetBinContent(2,i));
		hnt1S_y2_yMIdiff->SetBinContent(i,hnt1S_yMIdiff->GetBinContent(2,i));
		hnt2S_y2_yMIdiff->SetBinContent(i,hnt2S_yMIdiff->GetBinContent(2,i));
		hnt3S_y2_yMIdiff->SetBinContent(i,hnt3S_yMIdiff->GetBinContent(2,i));
		hnt1S_y2_yRfbdiff->SetBinContent(i,hnt1S_yRfbdiff->GetBinContent(2,i));
		hnt2S_y2_yRfbdiff->SetBinContent(i,hnt2S_yRfbdiff->GetBinContent(2,i));
		hnt3S_y2_yRfbdiff->SetBinContent(i,hnt3S_yRfbdiff->GetBinContent(2,i));
		hhf1S_y2_yPLdiff->SetBinContent(i,hhf1S_yPLdiff->GetBinContent(2,i));
		hhf2S_y2_yPLdiff->SetBinContent(i,hhf2S_yPLdiff->GetBinContent(2,i));
		hhf3S_y2_yPLdiff->SetBinContent(i,hhf3S_yPLdiff->GetBinContent(2,i));
		hhf1S_y2_yMIdiff->SetBinContent(i,hhf1S_yMIdiff->GetBinContent(2,i));
		hhf2S_y2_yMIdiff->SetBinContent(i,hhf2S_yMIdiff->GetBinContent(2,i));
		hhf3S_y2_yMIdiff->SetBinContent(i,hhf3S_yMIdiff->GetBinContent(2,i));
		hhf1S_y2_yRfbdiff->SetBinContent(i,hhf1S_yRfbdiff->GetBinContent(2,i));
		hhf2S_y2_yRfbdiff->SetBinContent(i,hhf2S_yRfbdiff->GetBinContent(2,i));
		hhf3S_y2_yRfbdiff->SetBinContent(i,hhf3S_yRfbdiff->GetBinContent(2,i));
	}
	TH1F* hnt1S_y3_yPLdiff = new TH1F("hnt1S_y3_yPLdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt2S_y3_yPLdiff = new TH1F("hnt2S_y3_yPLdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt3S_y3_yPLdiff = new TH1F("hnt3S_y3_yPLdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt1S_y3_yMIdiff = new TH1F("hnt1S_y3_yMIdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt2S_y3_yMIdiff = new TH1F("hnt2S_y3_yMIdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt3S_y3_yMIdiff = new TH1F("hnt3S_y3_yMIdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt1S_y3_yRfbdiff = new TH1F("hnt1S_y3_yRfbdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt2S_y3_yRfbdiff = new TH1F("hnt2S_y3_yRfbdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hnt3S_y3_yRfbdiff = new TH1F("hnt3S_y3_yRfbdiff","Ntrack bins, y0.80to1.20",4,ntbins);
	TH1F* hhf1S_y3_yPLdiff = new TH1F("hhf1S_y3_yPLdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf2S_y3_yPLdiff = new TH1F("hhf2S_y3_yPLdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf3S_y3_yPLdiff = new TH1F("hhf3S_y3_yPLdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf1S_y3_yMIdiff = new TH1F("hhf1S_y3_yMIdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf2S_y3_yMIdiff = new TH1F("hhf2S_y3_yMIdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf3S_y3_yMIdiff = new TH1F("hhf3S_y3_yMIdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf1S_y3_yRfbdiff = new TH1F("hhf1S_y3_yRfbdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf2S_y3_yRfbdiff = new TH1F("hhf2S_y3_yRfbdiff","HF bins, y0.80to1.20",4,hfbins);
	TH1F* hhf3S_y3_yRfbdiff = new TH1F("hhf3S_y3_yRfbdiff","HF bins, y0.80to1.20",4,hfbins);
	for (int i=1; i<=4; i++)
	{
		hnt1S_y3_yPLdiff->SetBinContent(i,hnt1S_yPLdiff->GetBinContent(3,i));
		hnt2S_y3_yPLdiff->SetBinContent(i,hnt2S_yPLdiff->GetBinContent(3,i));
		hnt3S_y3_yPLdiff->SetBinContent(i,hnt3S_yPLdiff->GetBinContent(3,i));
		hnt1S_y3_yMIdiff->SetBinContent(i,hnt1S_yMIdiff->GetBinContent(3,i));
		hnt2S_y3_yMIdiff->SetBinContent(i,hnt2S_yMIdiff->GetBinContent(3,i));
		hnt3S_y3_yMIdiff->SetBinContent(i,hnt3S_yMIdiff->GetBinContent(3,i));
		hnt1S_y3_yRfbdiff->SetBinContent(i,hnt1S_yRfbdiff->GetBinContent(3,i));
		hnt2S_y3_yRfbdiff->SetBinContent(i,hnt2S_yRfbdiff->GetBinContent(3,i));
		hnt3S_y3_yRfbdiff->SetBinContent(i,hnt3S_yRfbdiff->GetBinContent(3,i));
		hhf1S_y3_yPLdiff->SetBinContent(i,hhf1S_yPLdiff->GetBinContent(3,i));
		hhf2S_y3_yPLdiff->SetBinContent(i,hhf2S_yPLdiff->GetBinContent(3,i));
		hhf3S_y3_yPLdiff->SetBinContent(i,hhf3S_yPLdiff->GetBinContent(3,i));
		hhf1S_y3_yMIdiff->SetBinContent(i,hhf1S_yMIdiff->GetBinContent(3,i));
		hhf2S_y3_yMIdiff->SetBinContent(i,hhf2S_yMIdiff->GetBinContent(3,i));
		hhf3S_y3_yMIdiff->SetBinContent(i,hhf3S_yMIdiff->GetBinContent(3,i));
		hhf1S_y3_yRfbdiff->SetBinContent(i,hhf1S_yRfbdiff->GetBinContent(3,i));
		hhf2S_y3_yRfbdiff->SetBinContent(i,hhf2S_yRfbdiff->GetBinContent(3,i));
		hhf3S_y3_yRfbdiff->SetBinContent(i,hhf3S_yRfbdiff->GetBinContent(3,i));
	}
	TH1F* hnt1S_y4_yPLdiff = new TH1F("hnt1S_y4_yPLdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt2S_y4_yPLdiff = new TH1F("hnt2S_y4_yPLdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt3S_y4_yPLdiff = new TH1F("hnt3S_y4_yPLdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt1S_y4_yMIdiff = new TH1F("hnt1S_y4_yMIdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt2S_y4_yMIdiff = new TH1F("hnt2S_y4_yMIdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt3S_y4_yMIdiff = new TH1F("hnt3S_y4_yMIdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt1S_y4_yRfbdiff = new TH1F("hnt1S_y4_yRfbdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt2S_y4_yRfbdiff = new TH1F("hnt2S_y4_yRfbdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hnt3S_y4_yRfbdiff = new TH1F("hnt3S_y4_yRfbdiff","Ntrack bins, y1.20to1.93",4,ntbins);
	TH1F* hhf1S_y4_yPLdiff = new TH1F("hhf1S_y4_yPLdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf2S_y4_yPLdiff = new TH1F("hhf2S_y4_yPLdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf3S_y4_yPLdiff = new TH1F("hhf3S_y4_yPLdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf1S_y4_yMIdiff = new TH1F("hhf1S_y4_yMIdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf2S_y4_yMIdiff = new TH1F("hhf2S_y4_yMIdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf3S_y4_yMIdiff = new TH1F("hhf3S_y4_yMIdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf1S_y4_yRfbdiff = new TH1F("hhf1S_y4_yRfbdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf2S_y4_yRfbdiff = new TH1F("hhf2S_y4_yRfbdiff","HF bins, y1.20to1.93",4,hfbins);
	TH1F* hhf3S_y4_yRfbdiff = new TH1F("hhf3S_y4_yRfbdiff","HF bins, y1.20to1.93",4,hfbins);
	for (int i=1; i<=4; i++)
	{
		hnt1S_y4_yPLdiff->SetBinContent(i,hnt1S_yPLdiff->GetBinContent(4,i));
		hnt2S_y4_yPLdiff->SetBinContent(i,hnt2S_yPLdiff->GetBinContent(4,i));
		hnt3S_y4_yPLdiff->SetBinContent(i,hnt3S_yPLdiff->GetBinContent(4,i));
		hnt1S_y4_yMIdiff->SetBinContent(i,hnt1S_yMIdiff->GetBinContent(4,i));
		hnt2S_y4_yMIdiff->SetBinContent(i,hnt2S_yMIdiff->GetBinContent(4,i));
		hnt3S_y4_yMIdiff->SetBinContent(i,hnt3S_yMIdiff->GetBinContent(4,i));
		hnt1S_y4_yRfbdiff->SetBinContent(i,hnt1S_yRfbdiff->GetBinContent(4,i));
		hnt2S_y4_yRfbdiff->SetBinContent(i,hnt2S_yRfbdiff->GetBinContent(4,i));
		hnt3S_y4_yRfbdiff->SetBinContent(i,hnt3S_yRfbdiff->GetBinContent(4,i));
		hhf1S_y4_yPLdiff->SetBinContent(i,hhf1S_yPLdiff->GetBinContent(4,i));
		hhf2S_y4_yPLdiff->SetBinContent(i,hhf2S_yPLdiff->GetBinContent(4,i));
		hhf3S_y4_yPLdiff->SetBinContent(i,hhf3S_yPLdiff->GetBinContent(4,i));
		hhf1S_y4_yMIdiff->SetBinContent(i,hhf1S_yMIdiff->GetBinContent(4,i));
		hhf2S_y4_yMIdiff->SetBinContent(i,hhf2S_yMIdiff->GetBinContent(4,i));
		hhf3S_y4_yMIdiff->SetBinContent(i,hhf3S_yMIdiff->GetBinContent(4,i));
		hhf1S_y4_yRfbdiff->SetBinContent(i,hhf1S_yRfbdiff->GetBinContent(4,i));
		hhf2S_y4_yRfbdiff->SetBinContent(i,hhf2S_yRfbdiff->GetBinContent(4,i));
		hhf3S_y4_yRfbdiff->SetBinContent(i,hhf3S_yRfbdiff->GetBinContent(4,i));
	}
	
	
	cout << "PRINTING LATEX TABLES\n\n";
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{1S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
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
	cout << "$-2.87<y<-1.93$ " << lineValues[137] << "\\\\" << endl;
	//cout << "$-2.4<y<-1.93$ " << lineValues[7] << "\\\\" << endl;
	cout << "$-1.93<y<-1.2$ " << lineValues[8] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8$ " << lineValues[9] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4$ " << lineValues[10] << "\\\\" << endl;
	cout << "$-0.4<y<0.0$ " << lineValues[11] << "\\\\" << endl;
	cout << "$0.0<y<0.4$ " << lineValues[12] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[13] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[14] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[15] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError1S}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 2S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{2S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt$, y integrated " << lineValues[16] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 4 \\GeVc$ " << lineValues[17] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc$ " << lineValues[18] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc$ " << lineValues[19] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-2.87<y<-1.93$ " << lineValues[138] << "\\\\" << endl;
	//cout << "$-2.4<y<-1.93$ " << lineValues[20] << "\\\\" << endl;
	cout << "$-1.93<y<-0.8$ " << lineValues[21] << "\\\\" << endl;
	cout << "$-0.8<y<0.0$ " << lineValues[22] << "\\\\" << endl;
	cout << "$0.0<y<0.8$ " << lineValues[23] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[24] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError2S}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{3S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt$, y integrated " << lineValues[25] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 6 \\GeVc$ " << lineValues[26] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc$ " << lineValues[27] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-2.87<y<-1.93$ " << lineValues[139] << "\\\\" << endl;
	//cout << "$-2.4<y<-1.93$ " << lineValues[28] << "\\\\" << endl;
	cout << "$-1.93<y<0.0$ " << lineValues[29] << "\\\\" << endl;
	cout << "$0.0<y<1.93$ " << lineValues[30] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError3S}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	//Cross section pT bins (y-2.87-1.93)
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields for $-2.87<y<1.93$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|ccc|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c}{Yield Dev.($\\%$) } & \\\\" << endl;
	cout << "& & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c}{} & \\\\" << endl;
	cout << "$\\pt < 2 \\GeVc$ " << lineValues[140] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc$ " << lineValues[141] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc$ " << lineValues[142] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc$ " << lineValues[143] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc$ " << lineValues[144] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc$ " << lineValues[145] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c}{} & \\\\" << endl;
	cout << "$\\pt < 4 \\GeVc$ " << lineValues[146] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc$ " << lineValues[147] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc$ " << lineValues[148] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c}{} & \\\\" << endl;
	cout << "$\\pt < 6 \\GeVc$ " << lineValues[149] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc$ " << lineValues[150] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_crosspt}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	//////////////////////////////////////////////////Differential bins////////////////////////////////////////////
	/*
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError1S_difflo}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{1S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 2 \\GeVc, -1.93<y<0$ " << lineValues[31] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc, -1.93<y<0$ " << lineValues[32] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc, -1.93<y<0$ " << lineValues[33] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc, -1.93<y<0$ " << lineValues[34] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc, -1.93<y<0$ " << lineValues[35] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[36] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-1.2, \\pt<6 \\GeVc$ " << lineValues[37] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8, \\pt<6 \\GeVc$ " << lineValues[38] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4, \\pt<6 \\GeVc$ " << lineValues[39] << "\\\\" << endl;
	cout << "$-0.4<y<0.0, \\pt<6 \\GeVc$ " << lineValues[40] << "\\\\" << endl;
	cout << "$0.0<y<0.4, \\pt<6 \\GeVc$ " << lineValues[41] << "\\\\" << endl;
	cout << "$0.4<y<0.8, \\pt<6 \\GeVc$ " << lineValues[42] << "\\\\" << endl;
	cout << "$0.8<y<1.2, \\pt<6 \\GeVc$ " << lineValues[43] << "\\\\" << endl;
	cout << "$1.2<y<1.93, \\pt<6 \\GeVc$ " << lineValues[44] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError2S_difflo}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{2S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 4 \\GeVc, -1.93<y<0$ " << lineValues[45] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc, -1.93<y<0$ " << lineValues[46] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[47] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-0.8, \\pt<6 \\GeVc$ " << lineValues[48] << "\\\\" << endl;
	cout << "$-0.8<y<0.8, \\pt<6 \\GeVc$ " << lineValues[49] << "\\\\" << endl;
	cout << "$0.8<y<0.8, \\pt<6 \\GeVc$ " << lineValues[50] << "\\\\" << endl;
	cout << "$0.8<y<1.93, \\pt<6 \\GeVc$ " << lineValues[51] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 2S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError3S_difflo}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{3S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 6 \\GeVc, -1.93<y<0$ " << lineValues[52] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[53] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<0.0, \\pt<6 \\GeVc$ " << lineValues[54] << "\\\\" << endl;
	cout << "$0.0<y<1.93, \\pt<6 \\GeVc$ " << lineValues[55] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError1S_diffhi}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{1S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 2 \\GeVc, 0<y<1.93$ " << lineValues[56] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc, 0<y<1.93$ " << lineValues[57] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc, 0<y<1.93$ " << lineValues[58] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc, 0<y<1.93$ " << lineValues[59] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc, 0<y<1.93$ " << lineValues[60] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[61] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-1.2, 6<\\pt<30 \\GeVc$ " << lineValues[62] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8, 6<\\pt<30 \\GeVc$ " << lineValues[63] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4, 6<\\pt<30 \\GeVc$ " << lineValues[64] << "\\\\" << endl;
	cout << "$-0.4<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[65] << "\\\\" << endl;
	cout << "$0.0<y<0.4, 6<\\pt<30 \\GeVc$ " << lineValues[66] << "\\\\" << endl;
	cout << "$0.4<y<0.8, 6<\\pt<30 \\GeVc$ " << lineValues[67] << "\\\\" << endl;
	cout << "$0.8<y<1.2, 6<\\pt<30 \\GeVc$ " << lineValues[68] << "\\\\" << endl;
	cout << "$1.2<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[69] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError2S_diffhi}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{2S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 4 \\GeVc, 0<y<1.93$ " << lineValues[70] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc, 0<y<1.93$ " << lineValues[71] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[72] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<-0.8, 6<\\pt<30 \\GeVc$ " << lineValues[73] << "\\\\" << endl;
	cout << "$-0.8<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[74] << "\\\\" << endl;
	cout << "$0.0<y<0.8, 6<\\pt<30 \\GeVc$ " << lineValues[75] << "\\\\" << endl;
	cout << "$0.8<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[76] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 2S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError3S_diffhi}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{3S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$\\pt < 6 \\GeVc, 0<y<1.93$ " << lineValues[77] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[78] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "$-1.93<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[79] << "\\\\" << endl;
	cout << "$0.0<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[80] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	*/
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 2 \\GeVc, -1.93<y<0$ " << lineValues[31] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc, -1.93<y<0$ " << lineValues[32] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc, -1.93<y<0$ " << lineValues[33] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc, -1.93<y<0$ " << lineValues[34] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc, -1.93<y<0$ " << lineValues[35] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[36] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 4 \\GeVc, -1.93<y<0$ " << lineValues[45] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc, -1.93<y<0$ " << lineValues[46] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[47] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 6 \\GeVc, -1.93<y<0$ " << lineValues[52] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc, -1.93<y<0$ " << lineValues[53] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_diffptlo}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 2 \\GeVc, 0<y<1.93$ " << lineValues[56] << "\\\\" << endl;
	cout << "$2<\\pt<4 \\GeVc, 0<y<1.93$ " << lineValues[57] << "\\\\" << endl;
	cout << "$4<\\pt<6 \\GeVc, 0<y<1.93$ " << lineValues[58] << "\\\\" << endl;
	cout << "$6<\\pt<9 \\GeVc, 0<y<1.93$ " << lineValues[59] << "\\\\" << endl;
	cout << "$9<\\pt<12 \\GeVc, 0<y<1.93$ " << lineValues[60] << "\\\\" << endl;
	cout << "$12<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[61] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 4 \\GeVc, 0<y<1.93$ " << lineValues[70] << "\\\\" << endl;
	cout << "$4<\\pt<9 \\GeVc, 0<y<1.93$ " << lineValues[71] << "\\\\" << endl;
	cout << "$9<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[72] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$\\pt < 6 \\GeVc, 0<y<1.93$ " << lineValues[77] << "\\\\" << endl;
	cout << "$6<\\pt<30 \\GeVc, 0<y<1.93$ " << lineValues[78] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_diffpthi}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<-1.2, \\pt<6 \\GeVc$ " << lineValues[37] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8, \\pt<6 \\GeVc$ " << lineValues[38] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4, \\pt<6 \\GeVc$ " << lineValues[39] << "\\\\" << endl;
	cout << "$-0.4<y<0.0, \\pt<6 \\GeVc$ " << lineValues[40] << "\\\\" << endl;
	cout << "$0.0<y<0.4, \\pt<6 \\GeVc$ " << lineValues[41] << "\\\\" << endl;
	cout << "$0.4<y<0.8, \\pt<6 \\GeVc$ " << lineValues[42] << "\\\\" << endl;
	cout << "$0.8<y<1.2, \\pt<6 \\GeVc$ " << lineValues[43] << "\\\\" << endl;
	cout << "$1.2<y<1.93, \\pt<6 \\GeVc$ " << lineValues[44] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<-0.8, \\pt<6 \\GeVc$ " << lineValues[48] << "\\\\" << endl;
	cout << "$-0.8<y<0.8, \\pt<6 \\GeVc$ " << lineValues[49] << "\\\\" << endl;
	cout << "$0.8<y<0.8, \\pt<6 \\GeVc$ " << lineValues[50] << "\\\\" << endl;
	cout << "$0.8<y<1.93, \\pt<6 \\GeVc$ " << lineValues[51] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<0.0, \\pt<6 \\GeVc$ " << lineValues[54] << "\\\\" << endl;
	cout << "$0.0<y<1.93, \\pt<6 \\GeVc$ " << lineValues[55] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_diffylo}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{pA}$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<-1.2, 6<\\pt<30 \\GeVc$ " << lineValues[62] << "\\\\" << endl;
	cout << "$-1.2<y<-0.8, 6<\\pt<30 \\GeVc$ " << lineValues[63] << "\\\\" << endl;
	cout << "$-0.8<y<-0.4, 6<\\pt<30 \\GeVc$ " << lineValues[64] << "\\\\" << endl;
	cout << "$-0.4<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[65] << "\\\\" << endl;
	cout << "$0.0<y<0.4, 6<\\pt<30 \\GeVc$ " << lineValues[66] << "\\\\" << endl;
	cout << "$0.4<y<0.8, 6<\\pt<30 \\GeVc$ " << lineValues[67] << "\\\\" << endl;
	cout << "$0.8<y<1.2, 6<\\pt<30 \\GeVc$ " << lineValues[68] << "\\\\" << endl;
	cout << "$1.2<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[69] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<-0.8, 6<\\pt<30 \\GeVc$ " << lineValues[73] << "\\\\" << endl;
	cout << "$-0.8<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[74] << "\\\\" << endl;
	cout << "$0.0<y<0.8, 6<\\pt<30 \\GeVc$ " << lineValues[75] << "\\\\" << endl;
	cout << "$0.8<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[76] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$-1.93<y<0.0, 6<\\pt<30 \\GeVc$ " << lineValues[79] << "\\\\" << endl;
	cout << "$0.0<y<1.93, 6<\\pt<30 \\GeVc$ " << lineValues[80] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_diffyhi}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	//////////////////////////////NTrack, HF bins///////////////////////////////
	/*
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[81] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[82] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[83] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[84] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[85] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[86] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[87] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $0<Ntracks<40$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_nt040}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[88] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[89] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[90] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[91] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[92] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[93] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[94] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $40<Ntracks<65$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_nt4065}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[95] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[96] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[97] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[98] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[99] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[100] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[101] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $65<Ntracks<90$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_nt6590}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[102] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[103] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[104] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[105] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[106] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[107] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[108] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $90<Ntracks<400$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_nt90400}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[109] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[110] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[111] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[112] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[113] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[114] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[115] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $0<HF<15$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_hf015}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[116] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[117] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[118] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[119] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[120] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[121] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[122] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $15<HF<22$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_hf1522}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[123] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[124] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[125] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[126] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[127] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[128] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[129] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $22<HF<30$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_hf2230}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.4$ " << lineValues[130] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[131] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[132] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[133] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<0.8$ " << lineValues[134] << "\\\\" << endl;
	cout << "$0.8<y<1.93$ " << lineValues[135] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<y<1.93$ " << lineValues[136] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $30<HF<120$ due to variation of background PDF.}" << endl;
	cout << "\\label{sys:backgroundPDFError_hf30120}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	*/
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $-1.93<y<0$ and $0<y<1.93$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<Ntracks<40$ " << lineValues[151] << "\\\\" << endl;
	cout << "$40<Ntracks<62$ " << lineValues[152] << "\\\\" << endl;
	cout << "$62<Ntracks<88$ " << lineValues[153] << "\\\\" << endl;
	cout << "$88<Ntracks<400$ " << lineValues[154] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<Ntracks<40$ " << lineValues[155] << "\\\\" << endl;
	cout << "$40<Ntracks<62$ " << lineValues[156] << "\\\\" << endl;
	cout << "$62<Ntracks<88$ " << lineValues[157] << "\\\\" << endl;
	cout << "$88<Ntracks<400$ " << lineValues[158] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<Ntracks<40$ " << lineValues[159] << "\\\\" << endl;
	cout << "$40<Ntracks<400$ " << lineValues[160] << "\\\\" << endl;
	//cout << "$40<Ntracks<62$ " << lineValues[160] << "\\\\" << endl;
	//cout << "$62<Ntracks<88$ " << lineValues[161] << "\\\\" << endl;
	//cout << "$88<Ntracks<400$ " << lineValues[162] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_nt_yInt}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $-1.93<y<0$ and $0<y<1.93$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<HF<12$ " << lineValues[163] << "\\\\" << endl;
	cout << "$12<HF<19$ " << lineValues[164] << "\\\\" << endl;
	cout << "$19<HF<27$ " << lineValues[165] << "\\\\" << endl;
	cout << "$27<HF<120$ " << lineValues[166] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<HF<12$ " << lineValues[167] << "\\\\" << endl;
	cout << "$12<HF<19$ " << lineValues[168] << "\\\\" << endl;
	cout << "$19<HF<27$ " << lineValues[169] << "\\\\" << endl;
	cout << "$27<HF<120$ " << lineValues[170] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "3S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$0<HF<12$ " << lineValues[171] << "\\\\" << endl;
	cout << "$12<HF<120$ " << lineValues[172] << "\\\\" << endl;
	//cout << "$12<HF<19$ " << lineValues[172] << "\\\\" << endl;
	//cout << "$19<HF<27$ " << lineValues[173] << "\\\\" << endl;
	//cout << "$27<HF<400$ " << lineValues[174] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_hf_yInt}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\centering" << endl;
	cout << "\\caption{Systematic uncertainties of 1S, 2S, and 3S yields and $R_{FB}$ for $-1.93<y<0$ and $0<y<1.93$ due to variation of background PDF.}" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{c|}{Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\" << endl;
	cout << "& Backward y & Forward y &  \\\\" << endl;
	cout << "\\hline" << endl;
	cout << "1S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$Ntracks, HF Integrated$ " << lineValues[175] << "\\\\" << endl;
	cout << "2S & \\multicolumn{2}{c|}{} & \\\\" << endl;
	cout << "$Ntracks, HF Integrated$ " << lineValues[176] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\label{sys:backgroundPDFError_rfb_int}" << endl;
	cout << "\\end{table}" << endl;
	
	
	/*
	cout << endl << endl << "LATEX FIGURE INCLUDES" << endl << endl;
	
	///////////////LATEX figures////////////
	
	
	cout << "\\subsection{Results of Pseudoexperiments in $p_T$ and y Bins}" << endl;
	cout << "\\captionsetup{width=0.8\\textwidth,font=footnotesize}" << endl;
	for (int i=0, j=0; i<=30; i++, j++)
	{
		if (j<3) cout << "\\begin{figure}[htpb]" << endl;
		else cout << "\\begin{figure}[H]" << endl;
		cout << "\\centering" << endl;
		cout << "  \\includegraphics[width=0.8\\textwidth]{{figures/systematics_bkg/" << outNames[i] << "}.pdf}" << endl;
		cout << "  \\caption{ " << figCaption[i] << " }" << endl;
		//cout << "  \\label{fig:" << figureLabel[i] << "}" << endl;
		cout << "\\end{figure}" << endl;
		cout << "%" << endl;
		if (i==6) i = 136; else if (i==137) i = 6;
		if (i==19) i = 137; else if (i==138) i = 19;
		if (i==27) i = 138; else if (i==139) i = 27;
	}
	cout << "\\clearpage" << endl << endl;
	
	cout << "\\subsection{Results of Pseudoexperiments in $p_T$ Bins for $-2.87<y<1.93$}" << endl;
	cout << "\\captionsetup{width=0.8\\textwidth,font=footnotesize}" << endl;
	for (int i=140, j=0; i<=150; i++, j++)
	{
		if (j<3) cout << "\\begin{figure}[htpb]" << endl;
		else cout << "\\begin{figure}[H]" << endl;
		cout << "\\centering" << endl;
		cout << "  \\includegraphics[width=0.8\\textwidth]{{figures/systematics_bkg/" << outNames[i] << "}.pdf}" << endl;
		cout << "  \\caption{ " << figCaption[i] << " }" << endl;
		//cout << "  \\label{fig:" << figureLabel[i] << "}" << endl;
		cout << "\\end{figure}" << endl;
		cout << "%" << endl;
	}
	cout << "\\clearpage" << endl << endl;
	
	cout << "\\subsection{Results of Pseudoexperiments in Differential Bins}" << endl;
	cout << "\\captionsetup{width=0.8\\textwidth,font=footnotesize}" << endl;
	for (int i=31, j=0; i<=80; i++, j++)
	{
		if (j<3) cout << "\\begin{figure}[htpb]" << endl;
		else cout << "\\begin{figure}[H]" << endl;
		cout << "\\centering" << endl;
		cout << "  \\includegraphics[width=0.8\\textwidth]{{figures/systematics_bkg/" << outNames[i] << "}.pdf}" << endl;
		cout << "  \\caption{ " << figCaption[i] << " }" << endl;
		//cout << "  \\label{fig:" << figureLabel[i] << "}" << endl;
		cout << "\\end{figure}" << endl;
		cout << "%" << endl;
	}
	cout << "\\clearpage" << endl << endl;
	
	cout << "\\subsection{Results of Pseudoexperiments in NTrack Bins}" << endl;
	cout << "\\captionsetup{width=0.8\\textwidth,font=footnotesize}" << endl;
	for (int i=81, j=0; i<=162; i++, j++)
	{
		if (j<3) cout << "\\begin{figure}[htpb]" << endl;
		else cout << "\\begin{figure}[H]" << endl;
		cout << "\\centering" << endl;
		cout << "  \\includegraphics[width=0.8\\textwidth]{{figures/systematics_bkg/" << outNames[i] << "}.pdf}" << endl;
		cout << "  \\caption{ " << figCaption[i] << " }" << endl;
		//cout << "  \\label{fig:" << figureLabel[i] << "}" << endl;
		cout << "\\end{figure}" << endl;
		cout << "%" << endl;
		if (i==108) i=150;
	}
	cout << "\\clearpage" << endl << endl;
	
	cout << "\\subsection{Results of Pseudoexperiments in HF Bins}" << endl;
	cout << "\\captionsetup{width=0.8\\textwidth,font=footnotesize}" << endl;
	for (int i=109, j=0; i<=174; i++, j++)
	{
		if (j<3) cout << "\\begin{figure}[htpb]" << endl;
		else cout << "\\begin{figure}[H]" << endl;
		cout << "\\centering" << endl;
		cout << "  \\includegraphics[width=0.8\\textwidth]{{figures/systematics_bkg/" << outNames[i] << "}.pdf}" << endl;
		cout << "  \\caption{ " << figCaption[i] << " }" << endl;
		//cout << "  \\label{fig:" << figureLabel[i] << "}" << endl;
		cout << "\\end{figure}" << endl;
		cout << "%" << endl;
		if (i==136) i=162;
	}
	cout << "\\clearpage" << endl << endl;
	*/
	
	
	/////////////////////////////////////////output numbers to root file/////////////////////////////////////
	TFile* outfile = new TFile("BkgPdfSystematics.root","RECREATE");
	//resultTuple->Write();
	hInt1S_PADiff->Write();
	hInt2S_PADiff->Write();
	hInt3S_PADiff->Write();
	hpt1S_PADiff->Write();
	hpt2S_PADiff->Write();
	hpt3S_PADiff->Write();
	hy1S_PADiff->Write();
	hy2S_PADiff->Write();
	hy3S_PADiff->Write();
	
	hInt1S_PPDiff->Write();
	hInt2S_PPDiff->Write();
	hInt3S_PPDiff->Write();
	hpt1S_PPDiff->Write();
	hpt2S_PPDiff->Write();
	hpt3S_PPDiff->Write();
	hy1S_PPDiff->Write();
	hy2S_PPDiff->Write();
	hy3S_PPDiff->Write();
	
	hInt1S_RpADiff->Write();
	hInt2S_RpADiff->Write();
	hInt3S_RpADiff->Write();
	hpt1S_RpADiff->Write();
	hpt2S_RpADiff->Write();
	hpt3S_RpADiff->Write();
	hy1S_RpADiff->Write();
	hy2S_RpADiff->Write();
	hy3S_RpADiff->Write();
	
	//Cross section bins
	hy1S_cr_PADiff->Write();
	hy2S_cr_PADiff->Write();
	hy3S_cr_PADiff->Write();
	
	hpt1S_cr_PADiff->Write();
	hpt2S_cr_PADiff->Write();
	hpt3S_cr_PADiff->Write();
	
	//Differential bins
	hpt1S_ym_PADiff->Write();
	hpt2S_ym_PADiff->Write();
	hpt3S_ym_PADiff->Write();
	hpt1S_ym_PPDiff->Write();
	hpt2S_ym_PPDiff->Write();
	hpt3S_ym_PPDiff->Write();
	hpt1S_ym_RpADiff->Write();
	hpt2S_ym_RpADiff->Write();
	hpt3S_ym_RpADiff->Write();
	hpt1S_yp_PADiff->Write();
	hpt2S_yp_PADiff->Write();
	hpt3S_yp_PADiff->Write();
	hpt1S_yp_PPDiff->Write();
	hpt2S_yp_PPDiff->Write();
	hpt3S_yp_PPDiff->Write();
	hpt1S_yp_RpADiff->Write();
	hpt2S_yp_RpADiff->Write();
	hpt3S_yp_RpADiff->Write();
	
	hy1S_pt06_PADiff->Write();
	hy2S_pt06_PADiff->Write();
	hy3S_pt06_PADiff->Write();
	hy1S_pt06_PPDiff->Write();
	hy2S_pt06_PPDiff->Write();
	hy3S_pt06_PPDiff->Write();
	hy1S_pt06_RpADiff->Write();
	hy2S_pt06_RpADiff->Write();
	hy3S_pt06_RpADiff->Write();
	hy1S_pt630_PADiff->Write();
	hy2S_pt630_PADiff->Write();
	hy3S_pt630_PADiff->Write();
	hy1S_pt630_PPDiff->Write();
	hy2S_pt630_PPDiff->Write();
	hy3S_pt630_PPDiff->Write();
	hy1S_pt630_RpADiff->Write();
	hy2S_pt630_RpADiff->Write();
	hy3S_pt630_RpADiff->Write();
	
	//Rfb bins
	//2d hists
	/*hnt1S_yPLdiff->Write();
	hnt2S_yPLdiff->Write();
	hnt3S_yPLdiff->Write();
	hnt1S_yMIdiff->Write();
	hnt2S_yMIdiff->Write();
	hnt3S_yMIdiff->Write();
	hnt1S_yRfbdiff->Write();
	hnt2S_yRfbdiff->Write();
	hnt3S_yRfbdiff->Write();
	hhf1S_yPLdiff->Write();
	hhf2S_yPLdiff->Write();
	hhf3S_yPLdiff->Write();
	hhf1S_yMIdiff->Write();
	hhf2S_yMIdiff->Write();
	hhf3S_yMIdiff->Write();
	hhf1S_yRfbdiff->Write();
	hhf2S_yRfbdiff->Write();
	hhf3S_yRfbdiff->Write();*/
	
	//1D nt, hf bins
	/*hnt1S_y1_yPLdiff->Write();
	hnt2S_y1_yPLdiff->Write();
	hnt3S_y1_yPLdiff->Write();
	hnt1S_y1_yMIdiff->Write();
	hnt2S_y1_yMIdiff->Write();
	hnt3S_y1_yMIdiff->Write();
	hnt1S_y1_yRfbdiff->Write();
	hnt2S_y1_yRfbdiff->Write();
	hnt3S_y1_yRfbdiff->Write();
	hhf1S_y1_yPLdiff->Write();
	hhf2S_y1_yPLdiff->Write();
	hhf3S_y1_yPLdiff->Write();
	hhf1S_y1_yMIdiff->Write();
	hhf2S_y1_yMIdiff->Write();
	hhf3S_y1_yMIdiff->Write();
	hhf1S_y1_yRfbdiff->Write();
	hhf2S_y1_yRfbdiff->Write();
	hhf3S_y1_yRfbdiff->Write();
	
	hnt1S_y2_yPLdiff->Write();
	hnt2S_y2_yPLdiff->Write();
	//hnt3S_y2_yPLdiff->Write();
	hnt1S_y2_yMIdiff->Write();
	hnt2S_y2_yMIdiff->Write();
	//hnt3S_y2_yMIdiff->Write();
	hnt1S_y2_yRfbdiff->Write();
	hnt2S_y2_yRfbdiff->Write();
	//hnt3S_y2_yRfbdiff->Write();
	hhf1S_y2_yPLdiff->Write();
	hhf2S_y2_yPLdiff->Write();
	//hhf3S_y2_yPLdiff->Write();
	hhf1S_y2_yMIdiff->Write();
	hhf2S_y2_yMIdiff->Write();
	//hhf3S_y2_yMIdiff->Write();
	hhf1S_y2_yRfbdiff->Write();
	hhf2S_y2_yRfbdiff->Write();
	//hhf3S_y2_yRfbdiff->Write();
	
	hnt1S_y3_yPLdiff->Write();
	//hnt2S_y3_yPLdiff->Write();
	//hnt3S_y3_yPLdiff->Write();
	hnt1S_y3_yMIdiff->Write();
	//hnt2S_y3_yMIdiff->Write();
	//hnt3S_y3_yMIdiff->Write();
	hnt1S_y3_yRfbdiff->Write();
	//hnt2S_y3_yRfbdiff->Write();
	//hnt3S_y3_yRfbdiff->Write();
	hhf1S_y3_yPLdiff->Write();
	//hhf2S_y3_yPLdiff->Write();
	//hhf3S_y3_yPLdiff->Write();
	hhf1S_y3_yMIdiff->Write();
	//hhf2S_y3_yMIdiff->Write();
	//hhf3S_y3_yMIdiff->Write();
	hhf1S_y3_yRfbdiff->Write();
	//hhf2S_y3_yRfbdiff->Write();
	//hhf3S_y3_yRfbdiff->Write();
	
	hnt1S_y4_yPLdiff->Write();
	//hnt2S_y4_yPLdiff->Write();
	//hnt3S_y4_yPLdiff->Write();
	hnt1S_y4_yMIdiff->Write();
	//hnt2S_y4_yMIdiff->Write();
	//hnt3S_y4_yMIdiff->Write();
	hnt1S_y4_yRfbdiff->Write();
	//hnt2S_y4_yRfbdiff->Write();
	//hnt3S_y4_yRfbdiff->Write();
	hhf1S_y4_yPLdiff->Write();
	//hhf2S_y4_yPLdiff->Write();
	//hhf3S_y4_yPLdiff->Write();
	hhf1S_y4_yMIdiff->Write();
	//hhf2S_y4_yMIdiff->Write();
	//hhf3S_y4_yMIdiff->Write();
	hhf1S_y4_yRfbdiff->Write();
	//hhf2S_y4_yRfbdiff->Write();
	//hhf3S_y4_yRfbdiff->Write();*/
	
	hnt1S_yInt_PLdiff->Write();
	hnt2S_yInt_PLdiff->Write();
	hnt3S_yInt_PLdiff->Write();
	hnt1S_yInt_MIdiff->Write();
	hnt2S_yInt_MIdiff->Write();
	hnt3S_yInt_MIdiff->Write();
	hnt1S_yInt_Rfbdiff->Write();
	hnt2S_yInt_Rfbdiff->Write();
	hnt3S_yInt_Rfbdiff->Write();
	
	hhf1S_yInt_PLdiff->Write();
	hhf2S_yInt_PLdiff->Write();
	hhf3S_yInt_PLdiff->Write();
	hhf1S_yInt_MIdiff->Write();
	hhf2S_yInt_MIdiff->Write();
	hhf3S_yInt_MIdiff->Write();
	hhf1S_yInt_Rfbdiff->Write();
	hhf2S_yInt_Rfbdiff->Write();
	hhf3S_yInt_Rfbdiff->Write();
	
	hInt1S_yPLdiff->Write();
	hInt1S_yMIdiff->Write();
	hInt1S_Rfbdiff->Write();
	hInt2S_yPLdiff->Write();
	hInt2S_yMIdiff->Write();
	hInt2S_Rfbdiff->Write();
	
	outfile->Close();
}

void calcError(TNtuple* tup1, double* results, int state, TH1F* h1)
{
	float sigAlt1, sigNom1;
	tup1->SetBranchAddress(Form("nSig%dsAlt",state),&sigAlt1);
	tup1->SetBranchAddress(Form("nSig%dsNom",state),&sigNom1);
	
	int ntEntries = tup1->GetEntries();
	int denom = ntEntries;
	
	float sum1 = 0;
	float hi1 = 20;
	float lo1 = 0;
	
	for (int ent=0; ent<ntEntries; ent++)
	{
		tup1->GetEntry(ent);
		
		//Calculate differences from yields
		float diff1 = TMath::Abs((sigAlt1-sigNom1)/sigNom1);
		
		if (diff1>10) {denom--; continue;} //skip far outliers
		
		sum1 += diff1;
		if (diff1 > hi1) hi1 += 20; if (diff1 < lo1) lo1 -= 20;
	}
	float mean1 = sum1/denom;
	
	//h1 = new TH1F("hist1","hist1",100,lo1,hi1);
	
	float rms1 = 0;
	for (int ent=0; ent<ntEntries; ent++)
	{
		tup1->GetEntry(ent);
		
		//Calculate differences from yields
		float diff1 = TMath::Abs((sigAlt1-sigNom1)/sigNom1);
		
		if (diff1>10) {continue;} //skip far outliers
		h1->Fill(diff1*100.0);
		rms1 += TMath::Power(diff1-mean1,2);
	}
	rms1 = TMath::Sqrt(rms1/denom);;
	
	results[0] = mean1;
}

void calcError(TNtuple* tup1, TNtuple* tup2, double* results, int state, TH1F* h1, TH1F* h2, TH1F* hR)
{
	float sigAlt1, sigNom1;
	tup1->SetBranchAddress(Form("nSig%dsAlt",state),&sigAlt1);
	tup1->SetBranchAddress(Form("nSig%dsNom",state),&sigNom1);
	float sigAlt2, sigNom2;
	tup2->SetBranchAddress(Form("nSig%dsAlt",state),&sigAlt2);
	tup2->SetBranchAddress(Form("nSig%dsNom",state),&sigNom2);
	
	int ntEntries = tup1->GetEntries();
	int denom = ntEntries;
	
	float sum1 = 0;
	float sum2 = 0;
	float sumR = 0;
	float hi1 = 20; float hi2 = 20; float hiR = 20;
	
	for (int ent=0; ent<ntEntries; ent++)
	{
		tup1->GetEntry(ent);
		tup2->GetEntry(ent);
		
		//Calculate differences from yields
		float diff1 = TMath::Abs((sigAlt1-sigNom1)/sigNom1);
		float diff2 = TMath::Abs((sigAlt2-sigNom2)/sigNom2);
		float rAlt = sigAlt1/sigAlt2;
		float rNom = sigNom1/sigNom2;
		float diffR = TMath::Abs((rAlt-rNom)/rNom);
		
		if (diff1>10 || diff2>10 || diffR>10) {denom--; continue;} //skip far outliers
		
		sum1 += diff1;
		sum2 += diff2;
		sumR += diffR;
		if (diff1 > hi1) hi1 += 20;
		if (diff2 > hi2) hi2 += 20;
		if (diffR > hiR) hiR += 20;
	}
	float mean1 = sum1/denom;
	float mean2 = sum2/denom;
	float meanR = sumR/denom;
	
	/*h1 = new TH1F("hist1","hist1",100,lo1,hi1);
	h2 = new TH1F("hist2","hist2",100,lo2,hi2);
	hR = new TH1F("histR","histR",100,loR,hiR);*/
	/*int bindiv1 = 100/hi1; int bindiv2 = 100/hi2; int bindivR = 100/hiR;
	h1->Rebin(bindiv1);
	h2->Rebin(bindiv2);
	hR->Rebin(bindivR);*/
	
	float rms1 = 0; float rms2 = 0; float rmsR = 0;
	for (int ent=0; ent<ntEntries; ent++)
	{
		tup1->GetEntry(ent);
		tup2->GetEntry(ent);
		
		//Calculate differences from yields
		float diff1 = TMath::Abs((sigAlt1-sigNom1)/sigNom1);
		float diff2 = TMath::Abs((sigAlt2-sigNom2)/sigNom2);
		float rAlt = sigAlt1/sigAlt2;
		float rNom = sigNom1/sigNom2;
		float diffR = TMath::Abs((rAlt-rNom)/rNom);
		
		
		if (diff1>10 || diff2>10 || diffR>10) {continue;} //skip far outliers
		
		h1->Fill(diff1*100.0);
		h2->Fill(diff2*100.0);
		hR->Fill(diffR*100.0);
		rms1 += TMath::Power(diff1-mean1,2);
		rms2 += TMath::Power(diff2-mean2,2);
		rmsR += TMath::Power(diffR-meanR,2);
	}
	rms1 = TMath::Sqrt(rms1/denom); rms2 = TMath::Sqrt(rms2/denom); rmsR = TMath::Sqrt(rmsR/denom);
	
	results[0] = mean1;
	results[1] = mean2;
	results[2] = meanR;
}