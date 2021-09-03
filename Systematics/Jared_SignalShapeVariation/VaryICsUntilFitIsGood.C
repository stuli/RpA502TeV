#include <iostream>
#include "TRandom3.h"
#include "VaryICsFitter.C"


using namespace std;
using namespace RooFit;
void VaryICsUntilFitIsGood( 
       int collId = kPADATA,  
       float ptLow=6, float ptHigh=30, 
       float yLow=-1.93, float yHigh=-0.8,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			) 
{

  int goodreq = 6; //There are 8 parameters.
  //if (ptLow<=5) goodreq = 5;

  //array of initial parameter values.
  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_sigma,err_mu,m_lambda}
  double paramsICs[8] = {0};
  double* paramsICsptr;
  paramsICsptr = &paramsICs[0];

  //array that will hold the final parameter values.
  double params[8] = {0};
  double* paramsptr;
  paramsptr = &params[0];

  //arrays of upper and lower limits.
  double paramsupper[8] = {0.2, 3.0, 3.321, 5.0, 1.0, 25.0, 25.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0};

  bool GoodFit = kFALSE;
  int trycount = 0;

  while(!GoodFit) { //Fit until parameters are within bounds

    trycount += 1;

    //choose random initial values for all parameters within their ranges.
    TRandom3 rnd3(0);
    cout << "Starting fit attempt #" << trycount << " with the following ICs:" << endl;
    for (int i = 0; i<8; i++) {
      paramsICs[i] = rnd3.Rndm()*(paramsupper[i]-paramslower[i])+paramslower[i];
      cout << "paramsICs[" << i << "] = " << paramsICs[i] << endl;
    }

  cout << "paramsptr = " << paramsptr << endl;
  cout << "*paramsptr = " << *paramsptr << endl;

    //Run the fit
    VaryICsFitter(collId,ptLow,ptHigh,yLow,yHigh,cLow,cHigh,muPtCut,whichModel,paramsICsptr,paramsptr);
    cout << "paramsptr = " << paramsptr << endl;
    cout << "*paramsptr = " << *paramsptr << endl;

    //Check that it's okay. If not, try again.
    int goodcount = 0;
    for (int i = 0; i<8; i++) {
      cout << "*(paramsptr+" << i << ") = " << *(paramsptr+i) << endl;
      params[i] = *(paramsptr+i);
      double reldistupper = (paramsupper[i]-params[i])/(paramsupper[i]-paramslower[i]);
      double reldistlower = (params[i]-paramslower[i])/(paramsupper[i]-paramslower[i]);
      if ((reldistupper>0.05) && (reldistlower>0.05)) goodcount++;
    }
    cout << "goodcount = " << goodcount << endl;
    if (goodcount >= goodreq) {
      cout << "THE FIT IS A SUCCESS!!!!!!!!!!!!" << endl;
      break;
    }
    else cout << "THE FIT IS UNSATISFACTORY" << endl;
  }//end of while loop
  cout << "THE SUCCESSFUL ICS WERE: " << endl;
  for (int i = 0; i<8; i++) {
    cout << "paramsICs[" << i << "] = " << paramsICs[i] << endl;
  }
} 
 
