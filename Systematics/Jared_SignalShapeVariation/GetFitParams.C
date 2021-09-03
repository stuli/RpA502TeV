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
void GetFitParams( 
       int collId = kPADATA,
       float ptLow=0, float ptHigh=30,
       float yLow=-1.93, float yHigh=1.93,
       float* paramsptr, float* errorsptr
			) 
{
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  //import generating models
  //cout << "Importing workspace" << endl;
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString NomFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2018_03_20/nomfitresults_upsilon_%s.root",kineLabel.Data());
  //cout << NomFileName << endl;
  TFile* NomFile = TFile::Open(NomFileName,"READ");
  RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
  NomFile->Close("R");

  float ups1smass = Nomws->var("m_{#Upsilon(1S)}")->getVal();
  //cout << "m_{#Upsilon(1S)} = " << ups1smass << endl;

  //SIGNAL:
  double sigma1s_1_init = Nomws->var("sigma1s_1")->getVal();
  double x1s_init = Nomws->var("x1s")->getVal();
  double alpha1s_1_init = Nomws->var("alpha1s_1")->getVal();
  double n1s_1_init = Nomws->var("n1s_1")->getVal();
  double f1s_init = Nomws->var("f1s")->getVal();
  double nSig1s_init = Nomws->var("nSig1s")->getVal();
  double nSig2s_init = Nomws->var("nSig2s")->getVal();
  double nSig3s_init = Nomws->var("nSig3s")->getVal();
  /*cout << "sigma1s_1 = " << sigma1s_1_init << endl;
  cout << "x1s = " << x1s_init << endl;
  cout << "alpha1s_1 = " << alpha1s_1_init << endl;
  cout << "n1s_1 = " << n1s_1_init << endl;
  cout << "f1s = " << f1s_init << endl;
  cout << "nSig1s = " << nSig1s_init << endl;
  cout << "nSig2s = " << nSig2s_init << endl;
  cout << "nSig3s = " << nSig3s_init << endl;*/

  //BACKGROUND
  double err_mu_init = 0;
  double err_sigma_init = 0;
  if (ptLow<5) {
    err_mu_init = Nomws->var("#mu")->getVal();
    err_sigma_init = Nomws->var("#sigma")->getVal();
  }
  double m_lambda_init = Nomws->var("#lambda")->getVal();
  double nBkg_init = Nomws->var("nBkg")->getVal(); 
  /*cout << "#mu = " << err_mu_init << endl;
  cout << "#sigma = " << err_sigma_init << endl;
  cout << "#lambda = " << m_lambda_init << endl;
  cout << "nBkg = " << nBkg_init << endl;
*/
  //"m0:n:alpha:sigma0:f:x:mu:sigma:lambda"

  //free signal parameters
  *paramsptr = ups1smass;
  *(paramsptr+1) = sigma1s_1_init;
  *(paramsptr+2) = f1s_init;
  //fixed signal parameters
  *(paramsptr+3) = n1s_1_init;
  *(paramsptr+4) = alpha1s_1_init;
  *(paramsptr+5) = x1s_init;
  //yields
  *(paramsptr+6) = nSig1s_init;
  *(paramsptr+7) = nSig2s_init;
  *(paramsptr+8) = nSig3s_init;
  //background parameters
  *(paramsptr+9) = err_mu_init;
  *(paramsptr+10) = err_sigma_init;
  *(paramsptr+11) = m_lambda_init;
  *(paramsptr+12) = nBkg_init;

  //Get Errors
  float ups1smass_err = Nomws->var("m_{#Upsilon(1S)}")->getError();
  //cout << "m_{#Upsilon(1S)} = " << ups1smass << endl;

  //SIGNAL:
  double sigma1s_1_err = Nomws->var("sigma1s_1")->getError();
  double x1s_err = Nomws->var("x1s")->getError();
  double alpha1s_1_err = Nomws->var("alpha1s_1")->getError();
  double n1s_1_err = Nomws->var("n1s_1")->getError();
  double f1s_err = Nomws->var("f1s")->getError();
  double nSig1s_err = Nomws->var("nSig1s")->getError();
  double nSig2s_err = Nomws->var("nSig2s")->getError();
  double nSig3s_err = Nomws->var("nSig3s")->getError();
  /*cout << "sigma1s_1 = " << sigma1s_1_err << endl;
  cout << "x1s = " << x1s_err << endl;
  cout << "alpha1s_1 = " << alpha1s_1_err << endl;
  cout << "n1s_1 = " << n1s_1_err << endl;
  cout << "f1s = " << f1s_err << endl;
  cout << "nSig1s = " << nSig1s_err << endl;
  cout << "nSig2s = " << nSig2s_err << endl;
  cout << "nSig3s = " << nSig3s_err << endl;*/

  //BACKGROUND
  double err_mu_err = 0;
  double err_sigma_err = 0;
  if (ptLow<5) {
    err_mu_err = Nomws->var("#mu")->getError();
    err_sigma_err = Nomws->var("#sigma")->getError();
  }
  double m_lambda_err = Nomws->var("#lambda")->getError();
  double nBkg_err = Nomws->var("nBkg")->getError(); 
  /*cout << "#mu = " << err_mu_err << endl;
  cout << "#sigma = " << err_sigma_err << endl;
  cout << "#lambda = " << m_lambda_err << endl;
  cout << "nBkg = " << nBkg_err << endl;
*/
  //"m0:n:alpha:sigma0:f:x:mu:sigma:lambda"

  //free signal parameters
  *errorsptr = ups1smass_err;
  *(errorsptr+1) = sigma1s_1_err;
  *(errorsptr+2) = f1s_err;
  //fixed signal parameters
  *(errorsptr+3) = n1s_1_err;
  *(errorsptr+4) = alpha1s_1_err;
  *(errorsptr+5) = x1s_err;
  //yields
  *(errorsptr+6) = nSig1s_err;
  *(errorsptr+7) = nSig2s_err;
  *(errorsptr+8) = nSig3s_err;
  //background parameters
  *(errorsptr+9) = err_mu_err;
  *(errorsptr+10) = err_sigma_err;
  *(errorsptr+11) = m_lambda_err;
  *(errorsptr+12) = nBkg_err;
} 
 
