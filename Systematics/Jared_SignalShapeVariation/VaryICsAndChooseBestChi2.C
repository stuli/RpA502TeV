#include <iostream>
#include "TRandom3.h"
#include "VaryICsFitter.C"


using namespace std;
using namespace RooFit;
void VaryICsAndChooseBestChi2( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=30, 
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{
  double chi2 = 2.0;
  double chi2min = 2.0;
  int goodreq = 7; //There are 8 parameters.
  //if (ptLow<=5) goodreq = 5;

  //array of initial parameter values.
  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_sigma,err_mu,m_lambda}
  double paramsICs[8] = {0};
  double* paramsICsptr;
  paramsICsptr = &paramsICs[0];

  double bestICs[8] = {0};
  double* bestICsptr;
  bestICsptr = &bestICs[0];

  //array that will hold the final parameter values.
  double params[8] = {0};
  double* paramsptr;
  paramsptr = &params[0];

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.2, 7.0, 3.321, 5.0, 1.0, 25.0, 25.0, 25.0};
  double paramslower[8] = {0.02, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  for (int j = 0; j<10; j++) { //Fit 20 times and choose the best

    //choose random initial values for all parameters within their ranges.
    TRandom3 rnd3(0);
    cout << "Starting fit attempt #" << j+1 << " with the following ICs:" << endl;
    for (int i = 0; i<8; i++) {
      paramsICs[i] = rnd3.Rndm()*(paramsupper[i]-paramslower[i])+paramslower[i];
      cout << "paramsICs[" << i << "] = " << paramsICs[i] << endl;
    }

  cout << "paramsptr = " << paramsptr << endl;
  cout << "*paramsptr = " << *paramsptr << endl;

    //Run the fit
    chi2 = VaryICsFitter(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,paramsICsptr,paramsptr);
    cout << "paramsptr = " << paramsptr << endl;
    cout << "*paramsptr = " << *paramsptr << endl;

    //Check that it's okay. If not, try again.
    if (chi2<chi2min) {
      for (int i = 0; i<8; i++) {
        cout << "*(paramsptr+" << i << ") = " << *(paramsptr+i) << endl;
        bestICs[i] = *(paramsICsptr+i);
      }
    }
  }//end of loop

  cout << "THE BEST ICS WERE: " << endl;
  for (int i = 0; i<8; i++) {
    cout << "bestICs[" << i << "] = " << bestICs[i] << endl;
  }

  //Run the fit again with the best ICs
  chi2 = VaryICsFitter(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,bestICsptr,paramsptr);

} 
 
