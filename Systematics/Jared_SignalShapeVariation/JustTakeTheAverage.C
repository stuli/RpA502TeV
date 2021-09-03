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

TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

//4 bins for (PA, PP, RpA, RFB).
float alphasum[4] = {0};
float fsum[4] = {0};
float nsum[4] = {0};
float xsum[4] = {0};

int alphabins[4] = {0};
int fbins[4] = {0};
int nbins[4] = {0};
int xbins[4] = {0};

bool PPtoo, RFBBin;

float muPtCut=4.0;

TString getMyKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut) {
  TString kineLabeltemp = Form("%s_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, muPtCut);
  return kineLabeltemp.Data();
}

bool isAbout(float a, float b) {
  if (abs(a-b)<0.01) return kTRUE;
  else return kFALSE;
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

using namespace RooFit;

TString checkKineLabel, checkNomFileName;

int CheckFitQuality(RooWorkspace* ws, 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-1.93, float yHigh=1.93,
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int whichModel=0,   // Nominal = 0. Alternative = 1.
       int changeParam=0
			) 
{

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
  ws = (RooWorkspace*)NomFile->Get("workspace");
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
 
TString kineLabel, NomFileName;
TString PPkineLabel, PPFileName;
TString AltalphakineLabel, AltalphaFileName;
TString AltfkineLabel, AltfFileName;
TString AltnkineLabel, AltnFileName;
TString AltxkineLabel, AltxFileName;

int collIdPA = kPADATA;
int collIdPP = kPPDATA;

void doBin(RooWorkspace* ws,
	int whichUpsilon=1,
	float ptLow=0.0, float ptHigh=30.0,
	float yLow=-1.93, float yHigh=1.93,
       	float hfLow=0, float hfHigh=120,
       	int ntracksLow=0, int ntracksHigh=400,
       	bool whichModel=0) {

  PPtoo = kTRUE;

  RFBBin = kFALSE;
  if ((hfHigh-hfLow<120) || (ntracksHigh-ntracksLow<400)) RFBBin = kTRUE;

  float yLowPP, yHighPP;

  float PANomYield = 0;
  float PAAltalphaYield = 0;
  float PAAltfYield = 0;
  float PAAltnYield = 0;
  float PAAltxYield = 0;
  float PPNomYield = 0;
  float PPAltalphaYield = 0;
  float PPAltfYield = 0;
  float PPAltnYield = 0;
  float PPAltxYield = 0;

  if (yLow>=0) {
    yLowPP = yLow;
    yHighPP = yHigh;
  }
  else if (yHigh-yLow>3.8) {
    yLowPP = 0.0;
    yHighPP = 1.93;
  }
  else {
    yLowPP = abs(yHigh);
    yHighPP = abs(yLow);
  }
    if (yLow<-2.5) {
      PPtoo = kFALSE;
      //cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;

    if (RFBBin) {
      collIdPP = kPADATA;
      yLowPP = -1.93;
      yHighPP = 0.0;
      PPtoo = kTRUE;
    }
    else collIdPP = kPPDATA;

    //import fitted yields
    /*kineLabel = getMyKineLabel (collIdPA, ptLow, ptHigh, yLow, yHigh, muPtCut);
    if (RFBBin) kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    NomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabel.Data());
    cout << Form("nomfitresults_upsilon_%s.root",kineLabel.Data()) << endl;
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    NomFile->Close("R");
    PANomYield = Nomws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Nomws;
    delete NomFile;

    AltalphakineLabel = getMyKineLabel (collIdPA, ptLow, ptHigh, yLow, yHigh, muPtCut);
    if (RFBBin) AltalphakineLabel = AltalphakineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltalphaFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changealpha/nomfitresults_upsilon_%s.root",AltalphakineLabel.Data());
    TFile* AltalphaFile = TFile::Open(AltalphaFileName,"READ");
    RooWorkspace *Altalphaws = (RooWorkspace*)AltalphaFile->Get("workspace");
    AltalphaFile->Close("R");
    PAAltalphaYield = Altalphaws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Altalphaws;
    delete AltalphaFile;

    AltfkineLabel = getMyKineLabel (collIdPA, ptLow, ptHigh, yLow, yHigh, muPtCut);
    if (RFBBin) AltfkineLabel = AltfkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltfFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changef/nomfitresults_upsilon_%s.root",AltfkineLabel.Data());
    TFile* AltfFile = TFile::Open(AltfFileName,"READ");
    RooWorkspace *Altfws = (RooWorkspace*)AltfFile->Get("workspace");
    AltfFile->Close("R");
    PAAltfYield = Altfws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Altfws;
    delete AltfFile;

    AltnkineLabel = getMyKineLabel (collIdPA, ptLow, ptHigh, yLow, yHigh, muPtCut);
    if (RFBBin) AltnkineLabel = AltnkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltnFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changen/nomfitresults_upsilon_%s.root",AltnkineLabel.Data());
    TFile* AltnFile = TFile::Open(AltnFileName,"READ");
    RooWorkspace *Altnws = (RooWorkspace*)AltnFile->Get("workspace");
    AltnFile->Close("R");
    PAAltnYield = Altnws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Altnws;
    delete AltnFile;

    AltxkineLabel = getMyKineLabel (collIdPA, ptLow, ptHigh, yLow, yHigh, muPtCut);
    if (RFBBin) AltxkineLabel = AltxkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltxFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changex/nomfitresults_upsilon_%s.root",AltxkineLabel.Data());
    TFile* AltxFile = TFile::Open(AltxFileName,"READ");
    RooWorkspace *Altxws = (RooWorkspace*)AltxFile->Get("workspace");
    AltxFile->Close("R");
    PAAltxYield = Altxws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Altxws;
    delete AltxFile;

    if (PPtoo) {
      PPkineLabel = getMyKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut);
      if (RFBBin) PPkineLabel = PPkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
      TString PPNomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),PPkineLabel.Data());
      cout << Form("nomfitresults_upsilon_%s.root",PPkineLabel.Data()) << endl;
      TFile* PPNomFile = TFile::Open(PPNomFileName,"READ");
      RooWorkspace *PPNomws = (RooWorkspace*)PPNomFile->Get("workspace");
      PPNomFile->Close("R");
      PPNomYield = PPNomws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPNomws;
      delete PPNomFile;

      AltalphakineLabel = getMyKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut);
      if (RFBBin) AltalphakineLabel = AltalphakineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
      AltalphaFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changealpha/nomfitresults_upsilon_%s.root",AltalphakineLabel.Data());
      TFile* PPAltalphaFile = TFile::Open(AltalphaFileName,"READ");
      RooWorkspace *PPAltalphaws = (RooWorkspace*)PPAltalphaFile->Get("workspace");
      PPAltalphaFile->Close("R");
      PPAltalphaYield = PPAltalphaws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPAltalphaws;
      delete PPAltalphaFile;

      AltfkineLabel = getMyKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut);
      if (RFBBin) AltfkineLabel = AltfkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
      AltfFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changef/nomfitresults_upsilon_%s.root",AltfkineLabel.Data());
      TFile* PPAltfFile = TFile::Open(AltfFileName,"READ");
      RooWorkspace *PPAltfws = (RooWorkspace*)PPAltfFile->Get("workspace");
      PPAltfFile->Close("R");
      PPAltfYield = PPAltfws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPAltfws;
      delete PPAltfFile;

      AltnkineLabel = getMyKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut);
      if (RFBBin) AltnkineLabel = AltnkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
      AltnFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changen/nomfitresults_upsilon_%s.root",AltnkineLabel.Data());
      TFile* PPAltnFile = TFile::Open(AltnFileName,"READ");
      RooWorkspace *PPAltnws = (RooWorkspace*)PPAltnFile->Get("workspace");
      PPAltnFile->Close("R");
      PPAltnYield = PPAltnws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPAltnws;
      delete PPAltnFile;

      AltxkineLabel = getMyKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut);
      if (RFBBin) AltxkineLabel = AltxkineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
      AltxFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_changex/nomfitresults_upsilon_%s.root",AltxkineLabel.Data());
      TFile* PPAltxFile = TFile::Open(AltxFileName,"READ");
      RooWorkspace *PPAltxws = (RooWorkspace*)PPAltxFile->Get("workspace");
      PPAltxFile->Close("R");
      PPAltxYield = PPAltxws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPAltxws;
      delete PPAltxFile;
    }*/

    //Calculate errors
    float AltalphaPAerr = TMath::Abs((PAAltalphaYield-PANomYield)/PANomYield*100);
    float AltalphaPPerr = 0;
    float AltalphaRpAerr = 0;
    float AltfPAerr = TMath::Abs((PAAltfYield-PANomYield)/PANomYield*100);
    float AltfPPerr = 0;
    float AltfRpAerr = 0;
    float AltnPAerr = TMath::Abs((PAAltnYield-PANomYield)/PANomYield*100);
    float AltnPPerr = 0;
    float AltnRpAerr = 0;
    float AltxPAerr = TMath::Abs((PAAltxYield-PANomYield)/PANomYield*100);
    float AltxPPerr = 0;
    float AltxRpAerr = 0;
    if (PPtoo) {
      float NomRpA = PANomYield/PPNomYield;

      AltalphaPPerr = TMath::Abs((PPAltalphaYield-PPNomYield)/PPNomYield*100);
      float AltalphaRpA = PAAltalphaYield/PPAltalphaYield;
      AltalphaRpAerr = TMath::Abs((AltalphaRpA-NomRpA)/NomRpA*100);

      float AltfRpA = PAAltfYield/PPNomYield;
      if (RFBBin) {
        AltfPPerr = TMath::Abs((PPAltfYield-PPNomYield)/PPNomYield*100);
        AltfRpA = PAAltfYield/PPAltfYield;
      }
      AltfRpAerr = TMath::Abs((AltfRpA-NomRpA)/NomRpA*100);

      AltnPPerr = TMath::Abs((PPAltnYield-PPNomYield)/PPNomYield*100);
      float AltnRpA = PAAltnYield/PPAltnYield;
      AltnRpAerr = TMath::Abs((AltnRpA-NomRpA)/NomRpA*100);

      AltxPPerr = TMath::Abs((PPAltxYield-PPNomYield)/PPNomYield*100);
      float AltxRpA = PAAltxYield/PPAltxYield;
      AltxRpAerr = TMath::Abs((AltxRpA-NomRpA)/NomRpA*100);
    }

  int PAgood[4] = {0};
  int PPgood[4] = {0};

  for (int i=0;i<4;i++) {
    PAgood[i] = CheckFitQuality(ws,collIdPA,ptLow,ptHigh,yLow,yHigh, hfLow,hfHigh,ntracksLow,ntracksHigh,whichModel,i+1);
    if (PPtoo) PPgood[i] = CheckFitQuality(ws,collIdPP,ptLow,ptHigh,yLowPP,yHighPP, hfLow,hfHigh,ntracksLow,ntracksHigh,whichModel,i+1);
  }
  if (PAgood[0]) {
    alphabins[0] = alphabins[0] + 1;
    alphasum[0] = alphasum[0] + AltalphaPAerr;
  }
  if (PAgood[1]) {
    fbins[0] = fbins[0] + 1;
    fsum[0] = fsum[0] + AltfPAerr;
  }
  if (PAgood[2]) {
    nbins[0] = nbins[0] + 1;
    nsum[0] = nsum[0] + AltnPAerr;
  }
  if (PAgood[3]) {
    xbins[0] = xbins[0] + 1;
    xsum[0] = xsum[0] + AltxPAerr;
  }

  if (PPtoo && !RFBBin) {

    if (PPgood[0]) {
      alphabins[1] = alphabins[1] + 1;
      alphasum[1] = alphasum[1] + AltalphaPPerr;
    }
    //if (PPgood[1]) { //PP doesn't have an f constraint.
      //fbins[1] = fbins[1] + 1;
      //fsum[1] = fsum[1] + AltfPPerr;
    //}
    if (PPgood[2]) {
      nbins[1] = nbins[1] + 1;
      nsum[1] = nsum[1] + AltnPPerr;
    }
    if (PPgood[3]) {
      xbins[1] = xbins[1] + 1;
      xsum[1] = xsum[1] + AltxPPerr;
    }

    if (PAgood[0] && PPgood[0]) {
      alphabins[2] = alphabins[2] + 1;
      alphasum[2] = alphasum[2] + AltalphaRpAerr;
    }
    if (PAgood[1]) {
      fbins[2] = fbins[2] + 1;
      fsum[2] = fsum[2] + AltfRpAerr;
    }
    if (PAgood[2] && PPgood[2]) {
      nbins[2] = nbins[2] + 1;
      nsum[2] = nsum[2] + AltnRpAerr;
    }
    if (PAgood[3] && PPgood[3]) {
      xbins[2] = xbins[2] + 1;
      xsum[2] = xsum[2] + AltxRpAerr;
    }

  }

  if (RFBBin) {

    if (PAgood[0]) {
      alphabins[0] = alphabins[0] + 1;
      alphasum[0] = alphasum[0] + AltalphaPAerr;
    }
    if (PAgood[1]) {
      fbins[0] = fbins[0] + 1;
      fsum[0] = fsum[0] + AltfPAerr;
    }
    if (PAgood[2]) {
      nbins[0] = nbins[0] + 1;
      nsum[0] = nsum[0] + AltnPAerr;
    }
    if (PAgood[3]) {
      xbins[0] = xbins[0] + 1;
      xsum[0] = xsum[0] + AltxPAerr;
    }
    if (PPgood[0]) {
      alphabins[0] = alphabins[0] + 1;
      alphasum[0] = alphasum[0] + AltalphaPPerr;
    }
    if (PPgood[1]) {
      fbins[0] = fbins[0] + 1;
      fsum[0] = fsum[0] + AltfPPerr;
    }
    if (PPgood[2]) {
      nbins[0] = nbins[0] + 1;
      nsum[0] = nsum[0] + AltnPPerr;
    }
    if (PPgood[3]) {
      xbins[0] = xbins[0] + 1;
      xsum[0] = xsum[0] + AltxPPerr;
    }

    if (PAgood[0] && PPgood[0]) {
      alphabins[3] = alphabins[3] + 1;
      alphasum[3] = alphasum[3] + AltalphaRpAerr;
    }
    if (PAgood[1] && PPgood[1]) {
      fbins[3] = fbins[3] + 1;
      fsum[3] = fsum[3] + AltfRpAerr;
    }
    if (PAgood[2] && PPgood[2]) {
      nbins[3] = nbins[3] + 1;
      nsum[3] = nsum[3] + AltnRpAerr;
    }
    if (PAgood[3] && PPgood[3]) {
      xbins[3] = xbins[3] + 1;
      xsum[3] = xsum[3] + AltxRpAerr;
    }
  }

  ofstream logFile;
  logFile.open(Form("JustTakeTheAverageLog_%is.txt",whichUpsilon), fstream::in | fstream::out | fstream::app);
  logFile << endl;
  logFile << Form("nomfitresults_upsilon_%s.root",kineLabel.Data()) << endl;
  logFile << Form("nomfitresults_upsilon_%s.root",PPkineLabel.Data()) << endl;

  logFile << "Format = " << "(PAerr, PPerr, RpAerr)" << endl;

  logFile << "alphaerr = (" << AltalphaPAerr << ", " << AltalphaPPerr << ", " << AltalphaRpAerr << ")" << endl;
  logFile << "ferr = (" << AltfPAerr << ", " << AltfPPerr << ", " << AltfRpAerr << ")" << endl;
  logFile << "nerr = (" << AltnPAerr << ", " << AltnPPerr << ", " << AltnRpAerr << ")" << endl;
  logFile << "xerr = (" << AltxPAerr << ", " << AltxPPerr << ", " << AltxRpAerr << ")" << endl;
  logFile << endl;

  logFile.close();

  cout << "This is what's still in memory:" << endl;
  gDirectory->ls("-m");

}


void JustTakeTheAverage(int whichUpsilon=1) {

  ofstream logFile;
  //logFile.open("log.txt", fstream::in | fstream::out | fstream::app);
  logFile.open(Form("JustTakeTheAverageLog_%is.txt",whichUpsilon));
  logFile << "This is the log file" << endl;
  logFile.close();

  RooWorkspace* ws = new RooWorkspace();

if (whichUpsilon==1) {
//integrated bin
  /*doBin(whichUpsilon,0.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular pt bins
  doBin(whichUpsilon,0.0,2.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,2.0,4.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,4.0,6.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,9.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,9.0,12.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,12.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular y bins
  doBin(whichUpsilon,0.0,30.0,-1.93,-1.2,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,-1.2,-0.8,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,-0.8,-0.4,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,-0.4,0.0,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,0.4,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.4,0.8,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.8,1.2,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,1.2,1.93,0,120,0,400,0);

//pt bins in backward and forward rapitidy
  doBin(whichUpsilon,0.0,6.0,-1.93,0.0,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,-1.93,0.0,0,120,0,400,0);

  doBin(whichUpsilon,0.0,6.0,0.0,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,0.0,1.93,0,120,0,400,0);

//integrated bin in -2.87<y<1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,1.93,0,120,0,400,0);

//pt bins in rapidity range [-2.87,1.93]
  doBin(whichUpsilon,0.0,2.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,2.0,4.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,4.0,6.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,9.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,9.0,12.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,12.0,30.0,-2.87,1.93,0,120,0,400,0);

//-2.87<y<-1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,-1.93,0,120,0,400,0);

//HF bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,12,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,12,19,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,19,27,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,27,120,0,400,0);

//Ntracks bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,0,40,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,40,62,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,62,88,0);*/
  doBin(ws,whichUpsilon,0.0,30.0,0.0,1.93,0,120,88,400,0);
}
/*else if (whichUpsilon==2) {
//integrated bin
  doBin(whichUpsilon,0.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular pt bins
  doBin(whichUpsilon,0.0,4.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,4.0,9.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,9.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular y bins
  doBin(whichUpsilon,0.0,30.0,-1.93,-0.8,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,-0.8,0.0,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,0.8,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.8,1.93,0,120,0,400,0);

//pt bins in backward and forward rapitidy
  doBin(whichUpsilon,0.0,6.0,-1.93,0.0,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,-1.93,0.0,0,120,0,400,0);

  doBin(whichUpsilon,0.0,6.0,0.0,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,0.0,1.93,0,120,0,400,0);

//integrated bin in -2.87<y<1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,1.93,0,120,0,400,0);

//pt bins in rapidity range [-2.87,1.93]
  doBin(whichUpsilon,0.0,4.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,4.0,9.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,9.0,30.0,-2.87,1.93,0,120,0,400,0);

//-2.87<y<-1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,-1.93,0,120,0,400,0);

//HF bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,12,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,12,19,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,19,27,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,27,120,0,400,0);

//Ntracks bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,0,40,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,40,62,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,62,88,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,88,400,0);
}
else if (whichUpsilon==3) {
//integrated bin
  doBin(whichUpsilon,0.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular pt bins
  doBin(whichUpsilon,0.0,6.0,-1.93,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,-1.93,1.93,0,120,0,400,0);

//regular y bins
  doBin(whichUpsilon,0.0,30.0,-1.93,0.0,0,120,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,0,400,0);

//pt bins in backward and forward rapitidy
  doBin(whichUpsilon,0.0,6.0,-1.93,0.0,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,-1.93,0.0,0,120,0,400,0);

  doBin(whichUpsilon,0.0,6.0,0.0,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,0.0,1.93,0,120,0,400,0);

//integrated bin in -2.87<y<1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,1.93,0,120,0,400,0);

//pt bins in rapidity range [-2.87,1.93]
  doBin(whichUpsilon,0.0,6.0,-2.87,1.93,0,120,0,400,0);
  doBin(whichUpsilon,6.0,30.0,-2.87,1.93,0,120,0,400,0);

//-2.87<y<-1.93 in pPb.
  doBin(whichUpsilon,0.0,30.0,-2.87,-1.93,0,120,0,400,0);

//HF bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,12,0,400,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,12,120,0,400,0);

//Ntracks bins in backward and forward y
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,0,40,0);
  doBin(whichUpsilon,0.0,30.0,0.0,1.93,0,120,40,400,0);
}*/

  //calculate average errors
  for (int i=0;i<4;i++) {
    if (alphabins[i]>0) alphasum[i] = alphasum[i]/alphabins[i];
    else alphasum[i] = 0;
    if (fbins[i]>0) fsum[i] = fsum[i]/fbins[i];
    else fsum[i] = 0;
    if (nbins[i]>0) nsum[i] = nsum[i]/nbins[i];
    else nsum[i] = 0;
    if (xbins[i]>0) xsum[i] = xsum[i]/xbins[i];
    else xsum[i] = 0;
  }

  //ofstream logFile;
  logFile.open(Form("JustTakeTheAverageLog_%is.txt",whichUpsilon), fstream::in | fstream::out | fstream::app);

  logFile << endl << "Format = " << "(PAerr, PPerr, RpAerr, RFBerr)" << endl;

  logFile << "alphabins = (" << alphabins[0] << ", " << alphabins[1] << ", " << alphabins[2] << ", " << alphabins[3] << ")" << endl;
  logFile << "alphaerr = (" << alphasum[0] << ", " << alphasum[1] << ", " << alphasum[2] << ", " << alphasum[3] << ")" << endl;

  logFile << "fbins = (" << fbins[0] << ", " << fbins[1] << ", " << fbins[2] << ", " << fbins[3] << ")" << endl;
  logFile << "ferr = (" << fsum[0] << ", " << fsum[1] << ", " << fsum[2] << ", " << fsum[3] << ")" << endl;

  logFile << "nbins = (" << nbins[0] << ", " << nbins[1] << ", " << nbins[2] << ", " << nbins[3] << ")" << endl;
  logFile << "nerr = (" << nsum[0] << ", " << nsum[1] << ", " << nsum[2] << ", " << nsum[3] << ")" << endl;

  logFile << "xbins = (" << xbins[0] << ", " << xbins[1] << ", " << xbins[2] << ", " << xbins[3] << ")" << endl;
  logFile << "xerr = (" << xsum[0] << ", " << xsum[1] << ", " << xsum[2] << ", " << xsum[3] << ")" << endl;

  logFile.close();

  cout << endl << "Format = " << "(PAerr, PPerr, RpAerr, RFBerr)" << endl;

  cout << "alphabins = (" << alphabins[0] << ", " << alphabins[1] << ", " << alphabins[2] << ", " << alphabins[3] << ")" << endl;
  cout << "alphaerr = (" << alphasum[0] << ", " << alphasum[1] << ", " << alphasum[2] << ", " << alphasum[3] << ")" << endl;

  cout << "fbins = (" << fbins[0] << ", " << fbins[1] << ", " << fbins[2] << ", " << fbins[3] << ")" << endl;
  cout << "ferr = (" << fsum[0] << ", " << fsum[1] << ", " << fsum[2] << ", " << fsum[3] << ")" << endl;

  cout << "nbins = (" << nbins[0] << ", " << nbins[1] << ", " << nbins[2] << ", " << nbins[3] << ")" << endl;
  cout << "nerr = (" << nsum[0] << ", " << nsum[1] << ", " << nsum[2] << ", " << nsum[3] << ")" << endl;

  cout << "xbins = (" << xbins[0] << ", " << xbins[1] << ", " << xbins[2] << ", " << xbins[3] << ")" << endl;
  cout << "xerr = (" << xsum[0] << ", " << xsum[1] << ", " << xsum[2] << ", " << xsum[3] << ")" << endl;

  cout << "This is what's still in memory:" << endl;
  gDirectory->ls("-m");

  TH1F* halpha = new TH1F("halpha","halpha",4,0,4);
  TH1F* hf = new TH1F("hf","hf",4,0,4);
  TH1F* hn = new TH1F("hn","hn",4,0,4);
  TH1F* hx = new TH1F("hx","hx",4,0,4);
  for (int i=0;i<4;i++) {
    halpha->SetBinContent(i+1, alphasum[i]);
    hf->SetBinContent(i+1, fsum[i]);
    hn->SetBinContent(i+1, nsum[i]);
    hx->SetBinContent(i+1, xsum[i]);
  }

  //save histograms
  TFile outFile(Form("FixParamDeviationsAveraged_%is.root",whichUpsilon), "RECREATE");

  halpha->Write();
  hf->Write();
  hn->Write();
  hx->Write();

  delete ws;
  delete halpha;
  delete hf;
  delete hn;
  delete hx;

  outFile.Close();


}
