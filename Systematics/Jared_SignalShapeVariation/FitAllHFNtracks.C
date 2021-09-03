#include <iostream>
#include "FitDataWithNominalSeedsHFNtracks.C"


using namespace std;
using namespace RooFit;
void FitAllHFNtracks() {

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;
  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  int whichModel = 0;

  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ybins3[3] = {-1.93,0.0,1.93};
  //float ybins3[2] = {-1.93,1.93};
  float hfbins1[5] =  {0,12,19,27,120};
  int ntracksbins1[5] = {0,40,62,88,400};
  float hfbins2[5] =  {0,120};
  int ntracksbins2[5] = {0,400};
  float hfbins3[5] =  {0,12,120};
  int ntracksbins3[5] = {0,40,400};

  int numybins, numhfbins, numntracksbins;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = -1.93;
  float yHigh = 0.0;
  float hfLow = 0;
  float hfHigh = 120;
  int ntracksLow = 0;
  int ntracksHigh = 400;

  //Define pointers to sets of bins
  float *ybinsptr;
  float *hfbinsptr;
  int *ntracksbinsptr;
  int ystart;

  numhfbins = 4;
  numntracksbins = 4;
  //loop through all the bins and do the fit in each one.
  for (int whichUpsilon = 1; whichUpsilon<=3; whichUpsilon++) {

    numybins = 2;
    ybinsptr = &ybins3[0];
    if (whichUpsilon==1) {
      numhfbins = 4;
      hfbinsptr = &hfbins1[0];
      numntracksbins = 4;
      ntracksbinsptr = &ntracksbins1[0];
    }
    else if (whichUpsilon==2) {
      numhfbins = 1;
      hfbinsptr = &hfbins2[0];
      numntracksbins = 1;
      ntracksbinsptr = &ntracksbins2[0];
    }
    else if (whichUpsilon==3) {
      numhfbins = 2;
      hfbinsptr = &hfbins3[0];
      numntracksbins = 2;
      ntracksbinsptr = &ntracksbins3[0];
    }

    if (collId==kPADATA) ystart = 0;
    else if (collId==kPPDATA) ystart = numybins/2;

    for (int ihf = 0; ihf<numhfbins; ihf++) {
      for (int iy = ystart; iy<numybins; iy++) {

        ptLow = 0;
        ptHigh = 30;
        yLow = *(ybinsptr+iy);
        yHigh = *(ybinsptr+iy+1);
        hfLow = *(hfbinsptr+ihf);
        hfHigh = *(hfbinsptr+ihf+1);
        ntracksLow = 0;
        ntracksHigh = 400;

        cout << "(" << ptLow << "," << ptHigh << "," << yLow << "," << yHigh << "," << hfLow << "," << hfHigh << "," << ntracksLow << "," << ntracksHigh << ")" << endl;
        FitDataWithNominalSeedsHFNtracks( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, hfLow, hfHigh, ntracksLow, ntracksHigh, whichModel);
      }
    }

    for (int intracks = 0; intracks<numntracksbins; intracks++) {
      for (int iy = ystart; iy<numybins; iy++) {

        ptLow = 0;
        ptHigh = 30;
        yLow = *(ybinsptr+iy);
        yHigh = *(ybinsptr+iy+1);
        hfLow = 0;
        hfHigh = 120;
        ntracksLow = *(ntracksbinsptr+intracks);
        ntracksHigh = *(ntracksbinsptr+intracks+1);

        cout << "(" << ptLow << "," << ptHigh << "," << yLow << "," << yHigh << "," << hfLow << "," << hfHigh << "," << ntracksLow << "," << ntracksHigh << ")" << endl;
        FitDataWithNominalSeedsHFNtracks( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut, hfLow, hfHigh, ntracksLow, ntracksHigh, whichModel);
      }
    }

  }
  gROOT->ProcessLine(".q");
} 
 
