#include<fstream>
//#include<iostream>

using namespace RooFit;

int getMult(bool os = 1, int run=1, float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	float eta_low = -2.4;
	float eta_high = 2.4;
	
	float massLow = 8;
	float massHigh = 14;
	
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
		return 0;
	}
		
	float yLowLab = yLow+0.47;
    float yHighLab = yHigh+0.47;
	
	TString kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && mass>%.2f && mass<%.2f", 
							ptLow, ptHigh, yLowLab, yHighLab, eta_high, eta_low, eta_high, eta_low, massLow, massHigh);
							
	kineCut += (" && pt1>4.0 && pt2>4.0");
	
	RooWorkspace *ws = new RooWorkspace("workspace");
	RooDataSet* dataset = (RooDataSet*)f->Get("dataset");
	ws->import(*dataset);
	RooDataSet* reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
	reducedDS->SetName("reducedDS");
	ws->import(*reducedDS);
	
	int num = reducedDS->sumEntries();
	
	delete f;
	delete ws;
	delete dataset;
	delete reducedDS;
	
	return num;
}

TString multLine(TString bin="blank", float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	TString errline = "ERROR";
	
	float lum1 = 20.644;
	float lum2 = 13.958;
	
	if (bin.EqualTo("pt"))
		bin = Form("%.1f<pt<%.1f",ptLow,ptHigh);
	else if (bin.EqualTo("y"))
		bin = Form("%.2f<y<%.2f",yLow,yHigh);
	
	cerr << endl << "*****  Doing bin: " << bin << "  *****" << endl << endl;
	
	cerr << "Getting multiplicities..." << endl;
	
	int mult1SS = getMult(0,1,ptLow,ptHigh,yLow,yHigh) / lum1;
	int mult2SS = getMult(0,2,ptLow,ptHigh,yLow,yHigh) / lum2;
	double multDiffSS = 100.0*mult2SS/mult1SS;
	
	int mult1OS = getMult(1,1,ptLow,ptHigh,yLow,yHigh) / lum1;
	int mult2OS = getMult(1,2,ptLow,ptHigh,yLow,yHigh) / lum2;
	double multDiffOS = 100.0*mult2OS/mult1OS;
	
	cerr << "Loading constrained fit 1..." << endl;
	
	TFile* file1 = new TFile(Form("Run1FitsConstrained/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0.root",ptLow,ptHigh,yLow,yHigh),"READ");
	if (file1->IsZombie())
		{cerr << "FIT FILE 1 NOT FOUND" << endl; return errline;}
	else
		cerr << "Found fit file 1" << endl;
	RooWorkspace* ws1 = (RooWorkspace*)file1->Get("workspace");
	cerr << "Getting numbers" << endl;
	double sigs1 = ( ws1->var("nSig1s")->getVal() + ws1->var("nSig2s")->getVal() + ws1->var("nSig3s")->getVal() ) / lum1;
	double nBkg1 = ws1->var("nBkg")->getVal() / lum1;
	double nBkgErr1 = ws1->var("nBkg")->getError() / lum1;
	delete ws1;
	delete file1;
	
	cerr << "Loading constrained fit 2..." << endl;
	
	TFile* file2 = new TFile(Form("Run2FitsConstrained/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0.root",ptLow,ptHigh,yLow,yHigh),"READ");
	if (file2->IsZombie())
		{cerr << "FIT FILE 2 NOT FOUND" << endl; return errline;}
	else
		cerr << "Found fit file 2" << endl;
	RooWorkspace* ws2 = (RooWorkspace*)file2->Get("workspace");
	cerr << "Getting numbers" << endl;
	double sigs2 = ( ws2->var("nSig1s")->getVal() + ws2->var("nSig2s")->getVal() + ws2->var("nSig3s")->getVal() ) / lum2;
	double nBkg2 = ws2->var("nBkg")->getVal() / lum2;
	double nBkgErr2 = ws2->var("nBkg")->getError() / lum1;
	delete ws2;
	delete file2;
	
	double nBkgDiff = 100.0*nBkg2/nBkg1;
		
	TString line = bin + Form("\t%d\t%d\t%.1f%%\t%d\t%d\t%.1f%%\t%.0f\t%.0f\t%.1f%%\t%.1f\t%.1f\t%.0f\t%.0f",mult1SS,mult2SS,multDiffSS,mult1OS,mult2OS,multDiffOS,nBkg1,nBkg2,nBkgDiff,nBkgErr1,nBkgErr2,sigs1,sigs2);
	//TString line = bin + Form("\t%d\t%d\t%.1f%%\t%d\t%d\t%.1f%%",mult1SS,mult2SS,multDiffSS,mult1OS,mult2OS,multDiffOS);
	return line;
}

void CompareRuns()
{
	vector<TString> outlines;
	outlines.push_back("Bin\t\tRun1ss\tRun2ss\t%Diff\tRun1os\tRun2os\t%Diff\tnBkg1\tnBkg2\t%Diff\tBErr1\tBErr2\tSigs1\tSigs2");
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