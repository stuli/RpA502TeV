#include <iostream>
#include "FitMCData.C"


using namespace std;
using namespace RooFit;
void FitAllMC() {

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;
  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  int whichModel = 0;

  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[3] = {-1.93,0.0,1.93};

  int numptbins, numybins;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = -1.93;
  float yHigh = 1.93;

  //start at y=0 if using PP data
  if (collId==kPPDATA) yLow = 0.0;

  //Fit the pt, y integrated bin
  //if (collId==kPADATA) FitMCData(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,hfLow,hfHigh,ntracksLow,ntracksHigh,whichModel,1);

  //Define pointers to sets of bins
  float *ptbinsptr;
  float *ybinsptr;
  int ystart;
  
  double chisqpt[6] = {0};
  double chisqy[8] = {0};

  //loop through all the bins and do the fit in each one.
  for (int whichUpsilon =1; whichUpsilon<=1; whichUpsilon++) {

    if (whichUpsilon==1) {
      numptbins = 6;
      numybins = 8;
      ptbinsptr = &ptbins1[0];
      ybinsptr = &ybins1[0];
    }
    else if (whichUpsilon==2) {
      numptbins = 3;
      numybins = 4;
      ptbinsptr = &ptbins2[0];
      ybinsptr = &ybins2[0];
    }
    else if (whichUpsilon==3) {
      numptbins = 2;
      numybins = 2;
      ptbinsptr = &ptbins3[0];
      ybinsptr = &ybins3[0];
    }

    if (collId==kPADATA) ystart = 0;
    else if (collId==kPPDATA) ystart = numybins/2;

    for (int ipt = 0; ipt<numptbins; ipt++) {

      ptLow = *(ptbinsptr+ipt);
      ptHigh = *(ptbinsptr+ipt+1);
      //yLow = *(ybinsptr+ystart);//-1.93 for pPb, 0.0 for pp
      yLow = 0.0;
      yHigh = 1.93;

      cout << "(" << ptLow << "," << ptHigh << "," << yLow << "," << yHigh << ")" << endl;
      chisqpt[ipt] = FitMCData(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,hfLow,hfHigh,ntracksLow,ntracksHigh,whichModel,1);
    }

    for (int iy = ystart; iy<numybins; iy++) {

      ptLow = 6;
      ptHigh = 30;
      yLow = *(ybinsptr+iy);
      yHigh = *(ybinsptr+iy+1);

      cout << "(" << ptLow << "," << ptHigh << "," << yLow << "," << yHigh << ")" << endl;
      chisqy[iy] = FitMCData(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,hfLow,hfHigh,ntracksLow,ntracksHigh,whichModel,1);
    }

    for (int ipt = 0; ipt<numptbins; ipt++) cout << chisqpt[ipt] << endl;
    for (int iy = ystart; iy<numybins; iy++) cout << chisqy[iy] << endl;
  }
  gROOT->ProcessLine(".q");
} 
 
