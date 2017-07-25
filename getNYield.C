#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include <sstream>
#include "TH1.h"
#include "cutsAndBin.h"

using namespace std;

valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, 	float dphiEp2Low=0,  float dphiEp2High=0) ;
void getNYield(int state=1)
{
  
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double *yBin;

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
    nCentBins = nCentBins1s;  centBin = centBin1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
    nYBins = nYBins2S;  yBin = yBin2S;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
    nYBins = nYBins3S;  yBin = yBin3S;
  }
  
  //for pT 
  for(int i=1;i<=nPtBins;i++)
  {
    valErr yieldPP;
    valErr yieldAA;
    yieldPP = getYield(state, 0, ptBin[i-1],ptBin[i],0,2.4,0,200,0,100);
    yieldAA = getYield(state, 2, ptBin[i-1],ptBin[i],0,2.4,0,200,0,100);
    cout << Form(" PP %dS yield, pt  ",state) << ptBin[i-1] << " - " << ptBin[i] << " : " << yieldPP.val << endl;
    cout << Form(" AA %dS yield, pt  ",state) << ptBin[i-1] << " - " << ptBin[i] << " : " << yieldAA.val << endl;
  }

  //for rap
  for(int i=1;i<=nYBins;i++)
  {
    valErr yieldPP_y;
    valErr yieldAA_y;
    yieldPP_y = getYield(state, 0, 0,30,yBin[i-1],yBin[i],0,200,0,100);
    yieldAA_y = getYield(state, 2, 0,30,yBin[i-1],yBin[i],0,200,0,100);
    cout << Form(" PP %dS yield, y  ",state) << yBin[i-1] << " - " << yBin[i] << " : " << yieldPP_y.val << endl;
    cout << Form(" AA %dS yield, y  ",state) << yBin[i-1] << " - " << yBin[i] << " : " << yieldAA_y.val << endl;
  }

  //for pT 
  for(int i=1;i<=nCentBins;i++)
  {
    valErr yieldPP_c;
    valErr yieldAA_c;
    yieldPP_c = getYield(state, 0, 0,30,0,2.4,0,200,0,100);
    if(centBin[i-1] <60) yieldAA_c = getYield(state, 2, 0,30,0,2.4,centBin[i-1],centBin[i],0,100);
    else if(centBin[i-1] <60) yieldAA_c = getYield(state, 6, 0,30,0,2.4,centBin[i-1],centBin[i],0,100);
    cout << Form(" PP %dS yield, cent  ",state) << centBin[i-1]/2 << " - " << centBin[i]/2 << " : " << yieldPP_c.val << endl;
    cout << Form(" AA %dS yield, cent  ",state) << centBin[i-1]/2 << " - " << centBin[i]/2 << " : " << yieldAA_c.val << endl;
  }
}


valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
    float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString SignalCB = "Double";
  TFile* inf = new TFile(Form("../TESTPLOT_FORNOMINAL/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
//TFile* inf = new TFile(Form("../fitResults/dataFit_fixParam1MuPt4_2016_08_30/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}
