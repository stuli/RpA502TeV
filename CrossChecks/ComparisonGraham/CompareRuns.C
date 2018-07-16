#include<fstream>
//#include<iostream>

using namespace RooFit;

int getMult(int run=1, float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	float eta_low = -2.4;
	float eta_high = 2.4;
	
	float massLow = 8;
	float massHigh = 14;
	
	TFile* f;
	if (run==1)
		//f = new TFile("../../../yskimPA1st_SSign_201711181848_unIdentified.root","READ");
		f = new TFile("../../../yskimPA1st_OpSign_20177262037_unIdentified.root","READ");
	else if (run==2)
		//f = new TFile("../../../yskimPA2nd_SSign_201711182048_unIdentified.root","READ");
		f = new TFile("../../../yskimPA2nd_OpSign_20177262044_unIdentified.root","READ");
	else
	{
		cout << "INVALID RUN: " << run << endl;
		return 0;
	}
		
	float yLowLab = yLow+0.47;
    float yHighLab = yHigh+0.47;
	
	TString kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f && mass>%.2f && mass<%.2f", 
							ptLow, ptHigh, yLowLab, yHighLab, eta_high, eta_low, eta_high, eta_low, massLow, massHigh);
	
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

void getMults(int* mult1, int* mult2, int index, float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	mult1[index] = getMult(1,ptLow,ptHigh,yLow,yHigh);
	mult2[index] = getMult(2,ptLow,ptHigh,yLow,yHigh);
}

TString multLine(TString bin="blank", float ptLow=0, float ptHigh=30, float yLow=-1.93, float yHigh=1.93)
{
	if (bin.EqualTo("pt"))
		bin = Form("%.1f<pt<%.1f",ptLow,ptHigh);
	else if (bin.EqualTo("y"))
		bin = Form("%.2f<y<%.2f",yLow,yHigh);
	
	cout << endl << "*****  Doing bin: " << bin << "  *****" << endl << endl;
	
	int mult1 = getMult(1,ptLow,ptHigh,yLow,yHigh);
	int mult2 = getMult(2,ptLow,ptHigh,yLow,yHigh);
	double diff = 100.0*mult2/mult1;
	
	TFile* file1 = new TFile(Form("Run1FitsConstrained/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0.root",ptLow,ptHigh,yLow,yHigh),"READ");
	RooWorkspace* ws1 = (RooWorkspace*)file1->Get("workspace");
	double nBkg1 = ws1->var("nBkg")->getVal();
	delete ws1;
	delete file1;
	
	TFile* file2 = new TFile(Form("Run2FitsConstrained/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0.root",ptLow,ptHigh,yLow,yHigh),"READ");
	RooWorkspace* ws2 = (RooWorkspace*)file2->Get("workspace");
	double nBkg2 = ws2->var("nBkg")->getVal();
	delete ws2;
	delete file2;
	
	double diffbkg = 100.0*nBkg2/nBkg1;
		
	TString line = bin + Form("\t%d\t%d\t%.1f%%\t%.0f\t%.0f\t%.1f%%",mult1,mult2,diff,nBkg1,nBkg2,diffbkg);
	return line;
}

void CompareRuns()
{
	/*
	int mult1[7];
	int mult2[7];
	*/
	
	/*
	//int bin
	getMults(mult1,mult2,0,0,30,-1.93,1.93);
	//pt bins
	getMults(mult1,mult2,1,0,2,-1.93,1.93);
	getMults(mult1,mult2,2,2,4,-1.93,1.93);
	getMults(mult1,mult2,3,4,6,-1.93,1.93);
	getMults(mult1,mult2,4,6,9,-1.93,1.93);
	getMults(mult1,mult2,5,9,12,-1.93,1.93);
	getMults(mult1,mult2,6,12,30,-1.93,1.93);
	*/
	vector<TString> outlines;
	outlines.push_back("Bin\t\tRun1\tRun2\t%Diff\tnBkg1\tnBkg2\t%Diff");
	outlines.push_back(multLine("Int.\t",0,30,-1.93,1.93));
	outlines.push_back("1S bins");
	outlines.push_back(multLine("pt",0,2,-1.93,1.93));
	outlines.push_back(multLine("pt",2,4,-1.93,1.93));
	outlines.push_back(multLine("pt",4,6,-1.93,1.93));
	outlines.push_back(multLine("pt",6,9,-1.93,1.93));
	outlines.push_back(multLine("pt",9,12,-1.93,1.93));
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
	
	//cout << "Bin\t\tRun1\tRun2" << endl;
	//int bin
	//cout << "Int." << "\t\t" << mult1[0] << "\t" << mult2[0] << endl;
	//pt bins
	/*cout << "0<pt<2" << "\t\t" << mult1[1] << "\t" << mult2[1] << endl;
	cout << "2<pt<4" << "\t\t" << mult1[2] << "\t" << mult2[2] << endl;
	cout << "4<pt<6" << "\t\t" << mult1[3] << "\t" << mult2[3] << endl;
	cout << "6<pt<9" << "\t\t" << mult1[4] << "\t" << mult2[4] << endl;
	cout << "9<pt<12" << "\t\t" << mult1[5] << "\t" << mult2[5] << endl;
	cout << "12<pt<30" << "\t" << mult1[6] << "\t" << mult2[6] << endl;*/
	
	int numLines = outlines.size();
	for (int i = 0; i < numLines; i++)
		cout << outlines[i] << endl;
}