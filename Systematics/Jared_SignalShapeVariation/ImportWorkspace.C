//Author: Jared Jay
#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "../../HeaderFiles/PsetCollection.h"
#include "../../HeaderFiles/CMS_lumi.C"
#include "../../HeaderFiles/tdrstyle.C"
#include "../../HeaderFiles/StyleSetting.h"



using namespace std;
using namespace RooFit;
void ImportWorkspace( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-2.4, float yHigh=1.46,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
			) 
{

  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  //import generating model
  cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("nomfitresults_upsilon_%s.root",kineLabel.Data());
  cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *ws = (RooWorkspace*)NomFile->Get("workspace");

  float muval = ws->var("#mu")->getVal();
  cout << "muval = " << muval << endl;
}
