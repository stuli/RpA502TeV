//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

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

//TString directory = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithSimICs_AllParamFree_March15/";
TString directory = "RoundFits_R1b/";

  float cLow = 0; float cHigh = 200;
  float dphiEp2Low = 0; float dphiEp2High = 100;
  float muPtCut=4.0;

void GetFitResult(int collId=kPADATA, float ptLow=0.0, float ptHigh=30.0, float yLow=-1.93, float yHigh=1.93, double* alphasumptr=0, double* fsumptr=0, double* nsumptr=0, double* xsumptr=0, int* totalbinsptr=0) {
      TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString NomFileName = Form("%snomfitresults_upsilon_%s.root",directory.Data(),kineLabel.Data());
      TFile* NomFile = TFile::Open(NomFileName,"READ");
      cout << Form("pt%.1f-%.1f_y%.2f-%.2f:\n",ptLow,ptHigh,yLow,yHigh);
      RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");

      double alphaval = Nomws->var("alpha1s_1")->getVal();
      double alphaerr = Nomws->var("alpha1s_1")->getError();
      cout << "alpha1s_1 = " << alphaval << " +/- " << alphaerr << endl;
      double fval = Nomws->var("f1s")->getVal();
      double ferr = Nomws->var("f1s")->getError();
      cout << "f1s = " << fval << " +/- " << ferr << endl;
      double nval = Nomws->var("n1s_1")->getVal();
      double nerr = Nomws->var("n1s_1")->getError();
      cout << "n1s_1 = " << nval << " +/- " << nerr << endl;
      double xval = Nomws->var("x1s")->getVal();
      double xerr = Nomws->var("x1s")->getError();
      cout << "x1s = " << xval << " +/- " << xerr << endl;

      //old method
      if (alphaval<1.2 || alphaval>4.8) {
        cout << "bad alpha value -> skipping this bin." << endl;
        return;
      }
      if (nval<1.2 || nval>4.8) {
        cout << "bad n value -> skipping this bin." << endl;
        return;
      }
      //new method
      /*if (alphaval+alphaerr>=5.0 || alphaval-alphaerr<=1.0) {
        cout << "bad alpha value -> skipping this bin." << endl;
        return;
      }
      if (nval+nerr>=5.0 || nval-nerr<=1.0) {
        cout << "bad n value -> skipping this bin." << endl;
        return;
      }*/
      //For running on round 3: Actually it's too tight a constraint.
      /*if (fval+ferr/2>=1.0 || fval-ferr/2<=0.0) {
        cout << "bad f value -> skipping this bin." << endl;
        return;
      }
      if (xval+xerr>=1.0 || xval-xerr<=0.0) {
        cout << "bad x value -> skipping this bin." << endl;
        return;
      }*/
      /*if (fval>0.95 || fval<0.05) {
        cout << "bad f value -> skipping this bin." << endl;
        return;
      }
      if (xval>0.95 || xval<0.05) {
        cout << "bad x value -> skipping this bin." << endl;
        return;
      }*/

      //Sum to get mean
      *alphasumptr += alphaval;
      *fsumptr += fval;
      *nsumptr += nval;
      *xsumptr += xval;

      //Sum of squares to get RMS
      *(alphasumptr+1) += alphaval*alphaval;
      *(fsumptr+1) += fval*fval;
      *(nsumptr+1) += nval*nval;
      *(xsumptr+1) += xval*xval;

      //Sum of statistical errors to get mean statistical error
      *(alphasumptr+2) += alphaerr;
      *(fsumptr+2) += ferr;
      *(nsumptr+2) += nerr;
      *(xsumptr+2) += xerr;

      *totalbinsptr += 1;

      delete Nomws;
      NomFile->Close("R");
      delete NomFile;
}


TString CalculateAverage(double sum=0.0, double sumsqr=0.0, double sumerr=0.0, int totalbins=1) {

  double avg = sum/totalbins;
  double ms = sumsqr/totalbins;
  double sd = sqrt(ms-avg*avg);
  double meanerr = sumerr/totalbins;
  TString outString = Form("%f/%i = %f +/- %f; meanerr=%f",sum,totalbins,avg,sd,meanerr);
  return outString;
}


void FindAverageNSplit(int collId=kPADATA) {

  float ptbinsmid[6] = {0,2,4,6,9,12};
  float ybinsmid[4] = {-0.8,-0.4,0.0,0.4};

  float ptbinsext[2] = {12,30};
  float ybinsext[5] = {-1.93,-1.2,0.8,1.2,1.93};

  int numptbinsmid, numybinsmid, numptbinsext, numybinsext;
  float ptLow = 0;
  float ptHigh = 30;
  float yLow = -1.93;
  float yHigh = 1.93;

  double alphasum[3] = {0,0}; //sum, sum of squares, sum of errors.
  double fsum[3] = {0,0};
  double nsum[3] = {0,0};
  double xsum[3] = {0,0};

  //double nsumsqr = 0;
  //double xsumsqr = 0;
  //double alphasumsqr = 0;
  //double fsumsqr = 0;

  int totalbins = 0;
  double totalptbinsize = 30*3;
  double totalybinsize = 1.93*2*3;

  double* alphasumptr;
  double* fsumptr;
  double* nsumptr;
  double* xsumptr;
  int* totalbinsptr;

  alphasumptr = &alphasum[0];
  fsumptr = &fsum[0];
  nsumptr = &nsum[0];
  xsumptr = &xsum[0];
  totalbinsptr = &totalbins;

//Midrapidity:
  //GetFitResult(collId, 0, 2, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  //GetFitResult(collId, 2, 4, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  //GetFitResult(collId, 4, 6, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  //GetFitResult(collId, 6, 9, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  //GetFitResult(collId, 9, 12, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, -1.2, -0.8, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, -0.8, -0.4, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, -0.4, 0.0, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 0.0, 0.4, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 0.4, 0.8, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);

  GetFitResult(collId, 0, 30, -0.8, 0.0, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 0.0, 0.8, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);

  cout << "\nResults for midrapidity:\n";
  cout << "alphaavg = " << CalculateAverage(alphasum[0], alphasum[1], alphasum[2], totalbins) << endl;
  cout << "favg = " << CalculateAverage(fsum[0], fsum[1], fsum[2], totalbins) << endl;
  cout << "navg = " << CalculateAverage(nsum[0], nsum[1], nsum[2], totalbins) << endl;
  cout << "xavg = " << CalculateAverage(xsum[0], xsum[1], xsum[2], totalbins) << endl << endl;

//Extreme rapidity:
  alphasum[0]=0; alphasum[1]=0; alphasum[2]=0;
  fsum[0]=0; fsum[1]=0; fsum[2]=0;
  nsum[0]=0; nsum[1]=0; nsum[2]=0;
  xsum[0]=0; xsum[1]=0; xsum[2]=0;
  totalbins=0;
  //GetFitResult(collId, 12, 30, -1.93, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, -2.87, -1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, -1.93, -1.2, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 0.8, 1.2, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 1.2, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);

  GetFitResult(collId, 0, 30, -1.93, -0.8, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);
  GetFitResult(collId, 0, 30, 0.8, 1.93, alphasumptr, fsumptr, nsumptr, xsumptr, totalbinsptr);

  cout << "\nResults for extreme rapidity:\n";
  cout << "alphaavg = " << CalculateAverage(alphasum[0], alphasum[1], alphasum[2], totalbins) << endl;
  cout << "favg = " << CalculateAverage(fsum[0], fsum[1], fsum[2], totalbins) << endl;
  cout << "navg = " << CalculateAverage(nsum[0], nsum[1], nsum[2], totalbins) << endl;
  cout << "xavg = " << CalculateAverage(xsum[0], xsum[1], xsum[2], totalbins) << endl << endl;

}
