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
#include "RoundsHeader.h"

using namespace RooFit;

bool isAbout(float a, float b) {
  if (abs(a-b)<0.01) return kTRUE;
  else return kFALSE;
}

int CheckFitQuality( 
       int collId = kPPDATA,
       float ptLow=6, float ptHigh=9,
       float yLow=0.0, float yHigh=1.93,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int whichModel=0,   // Nominal = 0. Alternative = 1.
       int whichRound=R1a
			) 
{

  TString directory = Form("RoundFits_%s/",roundLabel[whichRound].Data());
  TString logFileName = Form("log_%s.txt",roundLabel[whichRound].Data());

  TString Params[5] = {"sigma1s_1","x1s","alpha1s_1","n1s_1","f1s"};
  //TString Params[2] = {"alpha1s_1","n1s_1"};
  const int numParams = 5;

  //Limits: {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  float massLow = 8; 
  float massHigh = 14;
  int   nMassBin  = (massHigh-massLow)*10;

  //import the model
  int binmode = 0;//The original skim file
  if (hfHigh-hfLow<120) binmode = 1;
  else if (ntracksHigh-ntracksLow<400) binmode = 2;
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  if (binmode>0) kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
  TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");

  RooAbsData* reducedDS = ws->data("reducedDS");

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
  double errUpperLimit = 0.15;
  double errLowerLimit = 0.005;
  if (collId==kPADATA) errUpperLimit = 0.15;

  //parameter requirements:
  bool goodParams = kTRUE;
  //TString badParamsList = "";
  std::ostringstream sstream;
  double buffer = 0.03;
  //for (int i=0;i<numParams;i++) {
  for (int i=1;i<5;i++) {
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

  //Check
  int good = 0;
  if (chisq/ndf>chisqUpperCut || chisq/ndf<chisqLowerCut) goodChi2 = kFALSE;
  else goodChi2 = kTRUE;
  if (temp1err/temp1>errUpperLimit || temp1err/temp1<errLowerLimit) goodSigErr = kFALSE;
  else goodSigErr = kTRUE;

  if (goodChi2 && goodSigErr && goodParams){
    cout << "THE FIT PASSED THE QUALITY CHECK! :)" << endl;
    good = 1;
  }
  else{
    cout << "THE FIT FAILED THE QUALITY CHECK! :(" << endl;
    if (!goodChi2) cout << "  -->bad chi^2 value = " << chisq/ndf << " is not within [" << chisqLowerCut << "-" << chisqUpperCut << "]" << endl;
    if (!goodSigErr) cout << "  -->bad signal error = " << temp1err/temp1*100 << "\%" << endl;
    if (!goodParams) cout << "  -->bad parameter (" << badParamsList << ") hits limit" << endl;

    good = 0;
  }
  cout << "good = " << good << endl;

  delete myPlot;
  //cout << "This is what's still in memory:" << endl;
  //gDirectory->ls("-m");

  return good;
} 
