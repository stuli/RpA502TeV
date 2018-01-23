void FitPseudoHists(int whichState = 1, int firstBin = 0, int lastBin = -1)
{
	float ptbins1s[7] = {0,2,4,6,9,12,30};
	//float ybins1s[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ybins1s[10] = {-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
	float ptbins2s[4] = {0,4,9,30};
	//float ybins2s[5] = {-1.93,-0.8,0.0,0.8,1.93};
	float ybins2s[6] = {-2.4,-1.93,-0.8,0.0,0.8,1.93};
	float ptbins3s[3] = {0,6,30};
	//float ybins3s[3] = {-1.93,0.0,1.93};
	float ybins3s[4] = {-2.4,-1.93,0.0,1.93};
	
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
		//numybins = 8;
		numybins = 9;
	}
	else if (whichState == 2)
	{
		ptbins = ptbins2s;
		ybins = ybins2s;
		numptbins = 3;
		//numybins = 4;
		numybins = 5;
	}
	else if (whichState == 3)
	{
		ptbins = ptbins3s;
		ybins = ybins3s;
		numptbins = 2;
		//numybins = 2;
		numybins = 3;
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
	
	TNtuple* resultTuple = new TNtuple("bkgSystTuple","Systematic Errors from Background Variation","ptLow:ptHigh:yLow:yHigh:yieldErrPA:yieldErrPP:RpAErr");
	//resultTuple->Fill(0,6,-1.93,1.93,0.204,.0601,0.1358);
	//resultTuple->Fill(0,30,0,1.93,0.1682,0.0525,0.2330);
	
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
				if (iState == 1) {hThisPA = hInt1S_PADiff; hThisPP = hInt1S_PPDiff; hThisRpA = hInt1S_RpADiff;}
				else if (iState == 2) {hThisPA = hInt2S_PADiff; hThisPP = hInt2S_PPDiff; hThisRpA = hInt2S_RpADiff;}
				else if (iState == 3) {hThisPA = hInt3S_PADiff; hThisPP = hInt3S_PPDiff; hThisRpA = hInt3S_RpADiff;}
				outbin = 1;
			}
			else if (i <= numptbins)
			{
				ptLow = ptbins[i-1];
				ptHigh = ptbins[i];
				yLow = -1.93;
				yHigh = 1.93;
				if (iState == 1) {hThisPA = hpt1S_PADiff; hThisPP = hpt1S_PPDiff; hThisRpA = hpt1S_RpADiff;}
				else if (iState == 2) {hThisPA = hpt2S_PADiff; hThisPP = hpt2S_PPDiff; hThisRpA = hpt2S_RpADiff;}
				else if (iState == 3) {hThisPA = hpt3S_PADiff; hThisPP = hpt3S_PPDiff; hThisRpA = hpt3S_RpADiff;}
				outbin = i;
			}
			else
			{
				ptLow = 0;
				ptHigh = 30;
				yLow = ybins[i-numptbins-1];
				yHigh = ybins[i-numptbins];
				if (iState == 1) {hThisPA = hy1S_PADiff; hThisPP = hy1S_PPDiff; hThisRpA = hy1S_RpADiff;}
				else if (iState == 2) {hThisPA = hy2S_PADiff; hThisPP = hy2S_PPDiff; hThisRpA = hy2S_RpADiff;}
				else if (iState == 3) {hThisPA = hy3S_PADiff; hThisPP = hy3S_PPDiff; hThisRpA = hy3S_RpADiff;}
				outbin = i-numptbins;
			}
			
			TString inFileName = Form("ResultsBkg/PseudoExpResults_pt%.1f-%.1f_y%.2f-%.2f.root",ptLow,ptHigh,yLow,yHigh);
			TFile* inFile = new TFile(inFileName,"READ");
			if (inFile->IsZombie())
				continue;
			
			TH1F* histoPA = (TH1F*)inFile->Get(Form("PA_%dSDiff",iState));
			TH1F* histoPP = (TH1F*)inFile->Get(Form("PP_%dSDiff",iState));
			TH1F* histoRpA = (TH1F*)inFile->Get(Form("RpA_%dSDiff",iState));
			
			/*double hPAxmin = histoPA->GetXaxis()->GetXmin();
			double hPAxmax = histoPA->GetXaxis()->GetXmin();
			double hPArange = hPAxmax - hPAxmin;
			double hPAmax = histoPA->GetMaximum();
			int hPAmaxbin = histoPA->GetMaximumBin();
			double hPAmuGuess = (hPAmaxbin - histoPA->GetNbinsX()/2)*hPArange/histoPA->GetNbinsX();*/
			
			/*TF1* fitPA = new TF1("fitPA","gaus(0)",hPAmuGuess-hPArange/5,hPAmuGuess+hPArange/5);
			fitPA->SetParNames("Amplitude","#mu","sigma");
			//fitPA->SetParameters(hPAmax,hPAmuGuess,1);
			fitPA->SetParameters(14,11,1);
			//fitPA->SetParLimits(1,hPAmuGuess-hPArange/10,hPAmuGuess+hPArange/10);*/
			TF1* fitPA = new TF1("fitPA","gaus(0)",histoPA->GetXaxis()->GetXmin(),histoPA->GetXaxis()->GetXmax());
			fitPA->SetParName(0,"Amplitude");
			fitPA->SetParName(1,"#mu");
			fitPA->SetParName(2,"#sigma");
			fitPA->SetParameters(histoPA->GetMaximum(),(histoPA->GetMaximumBin()-histoPA->GetNbinsX()/2)*(histoPA->GetXaxis()->GetXmax() - histoPA->GetXaxis()->GetXmin())/histoPA->GetNbinsX(),1);
			TF1* fitPP = new TF1("fitPP","gaus(0)",histoPP->GetXaxis()->GetXmin(),histoPP->GetXaxis()->GetXmax());
			fitPP->SetParName(0,"Amplitude");
			fitPP->SetParName(1,"#mu");
			fitPP->SetParName(2,"#sigma");
			fitPP->SetParameters(histoPP->GetMaximum(),(histoPP->GetMaximumBin()-histoPP->GetNbinsX()/2)*(histoPP->GetXaxis()->GetXmax() - histoPP->GetXaxis()->GetXmin())/histoPP->GetNbinsX(),1);
			TF1* fitRpA = new TF1("fitRpA","gaus(0)",histoRpA->GetXaxis()->GetXmin(),histoRpA->GetXaxis()->GetXmax());
			fitRpA->SetParName(0,"Amplitude");
			fitRpA->SetParName(1,"#mu");
			fitRpA->SetParName(2,"#sigma");
			fitRpA->SetParameters(histoRpA->GetMaximum(),(histoRpA->GetMaximumBin()-histoRpA->GetNbinsX()/2)*(histoRpA->GetXaxis()->GetXmax() - histoRpA->GetXaxis()->GetXmin())/histoRpA->GetNbinsX(),1);
			fitRpA->SetParLimits(2,0,50);
			
			histoPA->Fit(fitPA,"B");
			histoPP->Fit(fitPP,"B");
			histoRpA->Fit(fitRpA,"B");
			
			TCanvas* c1 = new TCanvas("cHists","%Diff histograms",4,45,1100,400);
			gStyle->SetOptFit(1);
			c1->Divide(3,1);
			c1->cd(1);
			histoPA->Draw();
			c1->cd(2);
			histoPP->Draw();
			c1->cd(3);
			histoRpA->Draw();
			
			TString outFileName = Form("ResultsBkg/PseudoExpResultsFitted_pt%.1f-%.1f_y%.2f-%.2f",ptLow,ptHigh,yLow,yHigh);
			c1->SaveAs(outFileName + ".pdf");
			
			TF1* fittedPA = histoPA->GetFunction("fitPA");
			TF1* fittedPP = histoPP->GetFunction("fitPP");
			TF1* fittedRpA = histoRpA->GetFunction("fitRpA");
			
			float PAdiff = TMath::Max(TMath::Max(TMath::Abs(histoPA->GetMean()),TMath::Abs(fittedPA->GetParameter(1))),TMath::Abs(histoPA->GetRMS()))*0.01;
			float PPdiff = TMath::Max(TMath::Max(TMath::Abs(histoPP->GetMean()),TMath::Abs(fittedPP->GetParameter(1))),TMath::Abs(histoPP->GetRMS()))*0.01;
			float RpAdiff = TMath::Max(TMath::Max(TMath::Abs(histoRpA->GetMean()),TMath::Abs(fittedRpA->GetParameter(1))),TMath::Abs(histoRpA->GetRMS()))*0.01;
			
			resultTuple->Fill(ptLow,ptHigh,yLow,yHigh,PPdiff,PAdiff,RpAdiff);
			
			hThisPA->SetBinContent(outbin,PAdiff);
			hThisPP->SetBinContent(outbin,PPdiff);
			hThisRpA->SetBinContent(outbin,RpAdiff);
			
			int lineNumber = i;
			if (iState > 1) lineNumber += 15;
			if (iState > 2) lineNumber += 8;
			
			//cout << "Setting line " << lineNumber << " to " << Form("& %.2f & %.2f & %.2f",finalPAdiff,finalPPdiff,finalRpAdiff) << endl;
			lineValues[lineNumber] = Form("& %.2f & %.2f & %.2f",PPdiff*100,PAdiff*100,RpAdiff*100);
		}
	}
	
	cout << "PRINTING LATEX TABLES\n\n";
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError1S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{1S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
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
	cout << "$-0.4<y<0.0$ " << lineValues[10] << "\\\\" << endl;
	cout << "$0.0<y<0.4$ " << lineValues[11] << "\\\\" << endl;
	cout << "$0.4<y<0.8$ " << lineValues[12] << "\\\\" << endl;
	cout << "$0.8<y<1.2$ " << lineValues[13] << "\\\\" << endl;
	cout << "$1.2<y<1.93$ " << lineValues[14] << "\\\\" << endl;
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\caption{Systematic uncertainties of 1S yields and RpA due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError2S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{2S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
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
	cout << "\\caption{Systematic uncertainties of 2S yields and RpA due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	cout << endl;
	
	cout << "\\begin{table}[hbtp]" << endl;
	cout << "\\label{sys:backgroundPDFError3S}" << endl;
	cout << "\\centering" << endl;
	cout << "\\begin{tabular}{|c|cc|c|}" << endl;
	cout << "\\hline" << endl;
	cout << "Bin  & \\multicolumn{2}{l}{3S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
	cout << "& pp & pPb &  \\\\" << endl;
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
	cout << "\\caption{Systematic uncertainties of 3S yields and RpA due to variation of background PDF.}" << endl;
	cout << "\\end{table}" << endl;
	
	//output numbers to root file
	TFile* outfile = new TFile("BkgPdfSystematics_2.4.root","RECREATE");
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
	
	outfile->Close();
}