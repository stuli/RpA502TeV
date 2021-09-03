using namespace std;

#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"

using namespace RooFit;

TString getMyKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut) {
  TString kineLabeltemp = Form("%s_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, muPtCut);
  return kineLabeltemp.Data();
}

bool isAbout(float a, float b) {
  if (abs(a-b)<0.01) return kTRUE;
  else return kFALSE;
}

int CheckFitQuality( 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-1.93, float yHigh=1.93,
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int whichModel=0,   // Nominal = 0. Alternative = 1.
       int changeParam=0
			) 
{

  TString checkKineLabel, checkNomFileName;
  TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

  bool PPtoo, RFBBin;
  float muPtCut=4.0;

  bool changen = kFALSE;
  bool changex = kFALSE;
  bool changealpha = kFALSE;
  bool changef = kFALSE;
  if (changeParam==1) changealpha = kTRUE;
  else if (changeParam==2) changef = kTRUE;
  else if (changeParam==3) changen = kTRUE;
  else if (changeParam==4) changex = kTRUE;
  TString directory = nominalDir;
  TString logFileName = "log.txt";
  if (changen) {
    directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changen/";
    logFileName = "log_changen.txt";
  }
  if (changex) {
    directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changex/";
    logFileName = "log_changex.txt";
  }
  if (changealpha) {
    directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changealpha/";
    logFileName = "log_changealpha.txt";
  }
  if (changef) {
    directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changef/";
    logFileName = "log_changef.txt";
  }

  TString Params[5] = {"sigma1s_1","x1s","alpha1s_1","n1s_1","f1s"};

  //Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  int binmode = 0;//The original skim file
  if (hfHigh-hfLow<120) binmode = 1;
  else if (ntracksHigh-ntracksLow<400) binmode = 2;

  //import the model
  cout << "Importing workspace" << endl;
  checkKineLabel = getMyKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut);
  if (binmode>0) checkKineLabel = checkKineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
  checkNomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),checkKineLabel.Data());
  cout << "Checking fit: " << checkNomFileName << endl;
  if (gSystem->AccessPathName(checkNomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(checkNomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");
  RooAbsData* reducedDS = ws->data("reducedDS");
  NomFile->Close("R");
  delete NomFile;

  RooFitResult* fitRes2 = (RooFitResult*)ws->obj("fitresult_model_reducedDS");

  //Plot it
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));

  ws->pdf("model")->plotOn(myPlot,Name("modelHist"));

  // PULL
  RooHist* hpull = myPlot->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);

  //calculate chi-squared
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin; 
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  cout << "numFitPar = " << numFitPar << endl;
  int ndf = nFullBinsPull - numFitPar;
  cout << "chisq/dof = " << chisq/ndf << endl;

  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  float temp2 = ws->var("nSig2s")->getVal();  
  float temp2err = ws->var("nSig2s")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();

  cout << "1S signal    =  " << temp1 << " +/- " << temp1err << endl;
  cout << "2S signal    =  " << temp2 << " +/- " << temp2err << endl;
  cout << "3S signal    =  " << temp3 << " +/- " << temp3err << endl;
  cout << "Total signal =  " << temp1+temp2+temp3 << endl;

  //QUALITY CHECK

  //chi-squared requirements:
  bool goodChi2 = kTRUE;
  double chisqUpperCut = 2.0 + temp1/10000;
  double chisqLowerCut = 0.5 + temp1/30000;

  //Signal error requirements:
  bool goodSigErr = kTRUE;
  double errUpperLimit = 0.1;
  double errLowerLimit = 0.005;
  if (collId==kPADATA) errUpperLimit = 0.15;

  //parameter requirements:
  bool goodParams = kTRUE;
  //TString badParamsList = "";
  std::ostringstream sstream;
  double buffer = 0.03;
  for (int i=0;i<5;i++) {
    double fittedval = ws->var(Params[i])->getVal();
    double width = paramsupper[i]-paramslower[i];
    double upper = paramsupper[i]-(buffer*width);
    double lower = paramslower[i]+(buffer*width);
    if (fittedval>upper || fittedval<lower) {
      goodParams = kFALSE;
      //badParamsList = badParamsList + " " + Params[i] + "=" + fittedval;
      sstream << " " << Params[i] << "=" << fittedval;
    }
  }
  std::string badParamsList = sstream.str();

  //Special bins
  if (collId==kPPDATA && isAbout(yLow,0.0) && isAbout(yHigh,0.4) && isAbout(ptHigh-ptLow,30)) chisqUpperCut = 3.5;
  if (collId==kPPDATA && isAbout(ptLow,0) && isAbout(ptHigh,6) && isAbout(yLow,0.0) && isAbout(yHigh,0.4)) chisqUpperCut = 2.8;
  if (collId==kPPDATA && isAbout(ptLow,0) && isAbout(ptHigh,6) && isAbout(yLow,0.0) && isAbout(yHigh,0.4)) chisqUpperCut = 3.5;
  if (collId==kPPDATA && isAbout(ptLow,0) && isAbout(ptHigh,6) && isAbout(yLow,0.0) && isAbout(yHigh,0.8)) chisqUpperCut = 3.5;
  if (collId==kPPDATA && isAbout(yLow,0.0) && isAbout(yHigh,0.8) && isAbout(ptHigh-ptLow,30)) chisqUpperCut = 4.6;
  //if (collId==kAADATA && isAbout(ptLow,4.0) && isAbout(ptHigh,9.0) && isAbout(yHigh-yLow,2.4)) chisqUpperCut = 2.3;

  delete myPlot;
  delete reducedDS;
  delete fitRes2;

  //Check
  ofstream logFile;
  logFile.open(logFileName, fstream::in | fstream::out | fstream::app);
  if (chisq/ndf>chisqUpperCut || chisq/ndf<chisqLowerCut) goodChi2 = kFALSE;
  if (temp1err/temp1>errUpperLimit || temp1err/temp1<errLowerLimit) goodSigErr = kFALSE;
  cout << "This is what's still in memory:" << endl;
  gDirectory->ls("-m");
  if (goodChi2 && goodSigErr && goodParams){
    //logFile << "THE FIT PASSED THE QUALITY CHECK! :)" << endl;
    logFile.close();
    return 1;
  }
  else{
    logFile << endl << checkNomFileName << endl;
    logFile << "THE FIT FAILED THE QUALITY CHECK! :(" << endl;
    if (!goodChi2) logFile << "  -->bad chi^2 value = " << chisq/ndf << " is not within [" << chisqLowerCut << "-" << chisqUpperCut << "]" << endl;
    if (!goodSigErr) logFile << "  -->bad signal error = " << temp1err/temp1*100 << "\%" << endl;
    if (!goodParams) logFile << "  -->bad parameter (" << badParamsList << ") hits limit" << endl;

    logFile.close();
    return 0;
  }
} 
