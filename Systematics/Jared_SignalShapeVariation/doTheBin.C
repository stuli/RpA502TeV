using namespace std;

#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
//#include "RooGaussian.h"
//#include "RooCBShape.h"
#include "RooWorkspace.h"
//#include "RooPlot.h"
//#include "TText.h"
//#include "TArrow.h"
#include "TFile.h"
//#include "../../HeaderFiles/cutsAndBin.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"

TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

TString kineLabel, NomFileName;

TString collIdPA = "PA_DATA";
float muPtCut=4.0;

void doTheBin(
	int whichUpsilon=1,
	float ptLow=0.0, float ptHigh=30.0,
	float yLow=-1.93, float yHigh=1.93,
       	float hfLow=0, float hfHigh=120,
       	int ntracksLow=0, int ntracksHigh=400,
       	bool whichModel=0) {

  float PANomYield = 0;

    //import fitted yields
    kineLabel = Form("%s_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f",collIdPA.Data(), ptLow,ptHigh, yLow, yHigh, muPtCut);
    NomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabel.Data());
    cout << Form("nomfitresults_upsilon_%s.root",kineLabel.Data()) << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace* Nomws = (RooWorkspace*)NomFile->Get("workspace");
    PANomYield = Nomws->var(Form("nSig%is",whichUpsilon))->getVal();
    //RooAbsData* reducedDS = Nomws->data("reducedDS");
    //delete reducedDS;
    //Nomws->SetOwner(kTRUE);
    //Nomws->Delete();
    delete Nomws;
    NomFile->Close("R");
    delete NomFile;

  cout << "This is what's still in memory:" << endl;
  gDirectory->ls();

}

