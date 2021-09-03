#include <iostream>
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "RoundsHeader.h"

using namespace std;
int CheckFitExists( 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-1.93, float yHigh=1.93,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       float hfLow=0, float hfHigh=120,
       int ntracksLow=0, int ntracksHigh=400,
       int whichModel=0,   // Nominal = 0. Alternative = 1.
       int whichRound=R1a
			) 
{

  TString directory = Form("RoundFits_%s/",roundLabel[whichRound].Data());

  int binmode = 0;//The original skim file
  if (hfHigh-hfLow<120) binmode = 1;
  else if (ntracksHigh-ntracksLow<400) binmode = 2;

  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  if (binmode>0) kineLabel = kineLabel + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
  TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
  cout << NomFileName << endl;
  if (gSystem->AccessPathName(NomFileName)) {
    cout << "THE FIT DOES NOT EXIST! :O";
    return 0;
  }
  else {
    return 1;
  }
} 
 
