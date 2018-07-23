using namespace RooFit;

vector<int> getMultInt(bool os = 1, int run=1, float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	float eta_low = -2.4;
	float eta_high = 2.4;
	
	int nbins = 100;
	float massLow = 9;
	float massHigh = 11;
	float binsize = (massHigh - massLow) / nbins;
	
	float massLow1s = 9.4;
	float massHigh1s = 9.8;
	float massLow2s = 9.8;
	float massHigh2s = 10.2;
	float massLow3s = 10.2;
	float massHigh3s = 10.4;
	
	TFile* f;
	if (run==1)
		if (os)
			f = new TFile("../../../yskimPA1st_OpSign_20177262037_unIdentified.root","READ");
		else
			f = new TFile("../../../yskimPA1st_SSign_201711181848_unIdentified.root","READ");
	else if (run==2)
		if (os)
			f = new TFile("../../../yskimPA2nd_OpSign_20177262044_unIdentified.root","READ");
		else
			f = new TFile("../../../yskimPA2nd_SSign_201711182048_unIdentified.root","READ");
	else
	{
		cout << "INVALID RUN: " << run << endl;
	}
		
	float yLowLab = yLow+0.47;
    float yHighLab = yHigh+0.47;
	
	TString kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && mass>%.2f && mass<%.2f", 
							ptLow, ptHigh, yLowLab, yHighLab, eta_high, eta_low, eta_high, eta_low, massLow, massHigh);
							
	kineCut += (" && pt1>4.0 && pt2>4.0");

	cout << "Getting tree" << endl;
	
	TTree* mmtree = (TTree*)f->Get("mm");
	
	TH1D* hist = new TH1D("hist","hist",nbins,massLow,massHigh);
	
	TCanvas* c1 = new TCanvas("c1","c1",600,400);
	mmtree->Draw("mass>>hist",kineCut);
	
	cerr << "Calculating numbers" << endl;
	
	vector<int> results = {-1,-2,-3};
	int binlow, binhigh;
	
	binlow = 2+(massLow1s-massLow)/binsize;
	binhigh = (massHigh1s-massLow)/binsize;
	cout << "1S: Integrating from " << binlow << " to " << binhigh << endl;
	results[0] = hist->Integral(binlow,binhigh);
	binlow = 2+(massLow2s-massLow)/binsize;
	binhigh = (massHigh2s-massLow)/binsize;
	cout << "2S: Integrating from " << binlow << " to " << binhigh << endl;
	results[1] = hist->Integral(binlow,binhigh);
	binlow = 2+(massLow3s-massLow)/binsize;
	binhigh = (massHigh3s-massLow)/binsize;
	cout << "3S: Integrating from " << binlow << " to " << binhigh << endl;
	results[2] = hist->Integral(binlow,binhigh);
	
	cout << "Numbers are " << results[0] << ", " << results[1] << ", " << results[2] << ", " << endl;
	
	//delete f;
	//delete mmtree;
	//delete hist;
	//delete c1;
	
	return results;
}

TString multLine(TString bin="blank", float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	TString errline = "ERROR";
	
	float lum1 = 20.644;
	float lum2 = 13.958;
	bool os = 0;
	
	if (bin.EqualTo("pt"))
		bin = Form("%.1f<pt<%.1f",ptLow,ptHigh);
	else if (bin.EqualTo("y"))
		bin = Form("%.2f<y<%.2f",yLow,yHigh);
	
	cerr << "*****  Doing bin: " << bin << "  *****" << endl;
	
	cerr << "Getting multiplicities..." << endl;
	
	vector<int> multsRun1SS = getMultInt(os,1,ptLow,ptHigh,yLow,yHigh);
	cerr << "Returned" << endl;
	double mult1sRun1SS = multsRun1SS[0] / lum1;
	double mult2sRun1SS = multsRun1SS[1] / lum1;
	double mult3sRun1SS = multsRun1SS[2] / lum1;
	vector<int> multsRun2SS = getMultInt(os,2,ptLow,ptHigh,yLow,yHigh);
	double mult1sRun2SS = multsRun2SS[0] / lum2;
	double mult2sRun2SS = multsRun2SS[1] / lum2;
	double mult3sRun2SS = multsRun2SS[2] / lum2;
	
	double avgErr1s = (TMath::Sqrt(mult1sRun1SS)+TMath::Sqrt(mult1sRun2SS))/2;
	double avgErr2s = (TMath::Sqrt(mult2sRun1SS)+TMath::Sqrt(mult2sRun2SS))/2;
	double avgErr3s = (TMath::Sqrt(mult3sRun1SS)+TMath::Sqrt(mult3sRun2SS))/2;
	
	double mult1sDiffSS = (mult1sRun1SS-mult1sRun2SS)/avgErr1s;
	double mult2sDiffSS = (mult2sRun1SS-mult2sRun2SS)/avgErr2s;
	double mult3sDiffSS = (mult3sRun1SS-mult3sRun2SS)/avgErr3s;
	
	/*double mult1sDiffSS = 100.0*mult1sRun2SS/mult1sRun1SS;
	double mult2sDiffSS = 100.0*mult2sRun2SS/mult2sRun1SS;
	double mult3sDiffSS = 100.0*mult3sRun2SS/mult3sRun1SS;*/
		
	TString line = bin + Form("\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",mult1sRun1SS,mult1sRun2SS,mult1sDiffSS,mult2sRun1SS,mult2sRun2SS,mult2sDiffSS,mult3sRun1SS,mult3sRun2SS,mult3sDiffSS);
	return line;
}

void CompareRunsIntegrals()
{
	vector<TString> outlines;
	outlines.push_back("Bin\t\t1S\t \t \t2S\t \t \t3S\t \t ");
	outlines.push_back("Bin\t\tRun1\tRun2\tDiff\tRun1\tRun2\tDiff\tRun1\tRun2\tDiff");
	outlines.push_back(multLine("Int.\t",0,30,-1.93,1.93));
	outlines.push_back("1S bins");
	outlines.push_back(multLine("pt",0,2,-1.93,1.93));
	outlines.push_back(multLine("pt",2,4,-1.93,1.93));
	outlines.push_back(multLine("pt",4,6,-1.93,1.93));
	outlines.push_back(multLine("pt",6,9,-1.93,1.93));
	outlines.push_back(multLine("pt",9,12,-1.93,1.93));
	outlines.push_back(multLine("pt",12,30,-1.93,1.93));
	outlines.push_back(multLine("y",0,30,-1.93,-1.2));
	outlines.push_back(multLine("y",0,30,-1.2,-0.8));
	outlines.push_back(multLine("y",0,30,-0.8,-0.4));
	outlines.push_back(multLine("y",0,30,-0.4,0));
	outlines.push_back(multLine("y",0,30,0,0.4));
	outlines.push_back(multLine("y",0,30,0.4,0.8));
	outlines.push_back(multLine("y",0,30,0.8,1.2));
	outlines.push_back(multLine("y",0,30,1.2,1.93));
	outlines.push_back("2S bins");
	outlines.push_back(multLine("pt",0,4,-1.93,1.93));
	outlines.push_back(multLine("pt",4,9,-1.93,1.93));
	outlines.push_back(multLine("pt",9,30,-1.93,1.93));
	outlines.push_back(multLine("y",0,30,-1.93,-0.8));
	outlines.push_back(multLine("y",0,30,-0.8,0));
	outlines.push_back(multLine("y",0,30,0,0.8));
	outlines.push_back(multLine("y",0,30,0.8,1.93));
	outlines.push_back("3S bins");
	outlines.push_back(multLine("pt",0,6,-1.93,1.93));
	outlines.push_back(multLine("pt",6,30,-1.93,1.93));
	outlines.push_back(multLine("y",0,30,-1.93,0));
	outlines.push_back(multLine("y",0,30,0,1.93));
	
	int numLines = outlines.size();
	for (int i = 0; i < numLines; i++)
		cout << outlines[i] << endl;
}