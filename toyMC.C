#include "rootFitHeaders.h"
#include "commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include "PsetCollection.h"
#include "RooRandom.h"

int kChPol3 = 1 ;
int kErrExp = 2 ;
int kChPol4 = 3 ;
int kErrExpExp = 4 ;

using namespace std;
using namespace RooFit;
void toyMC(
	   int collId = kAADATA,
	   float ptLow=0, float ptHigh=5,
	   float yLow=0, float yHigh=2.4,
	   int cLow=0, int cHigh=200,
	   float muPtCut=4.0,
	   int inputOption=kChPol4, //kChPol3,
	   int nGen = 10000,
	   int useCentIntBkgShape = 1,
     int nToys = 1000,
     int rdmseed = 111
	    ) 
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(rdmseed);
  gStyle->SetEndErrorSize(0);
  float Val_1S_nom = 0;
  float Val_1S_alt = 0;
  float Dev_1S = 0;
  float Val_2S_nom = 0;
  float Val_2S_alt = 0;
  float Dev_2S = 0;
  float Val_3S_nom = 0;
  float Val_3S_alt = 0;
  float Dev_3S = 0;
 
  int nGenCent = 0;

  if(collId == kAADATAPeri) collId =2; 
  //nGen = nGenCent;
  TString fcoll;
  TString finput;
  if(collId == kAADATA) fcoll = "AA";
  else if(collId == kPPDATA) fcoll = "PP";
  if(inputOption == 3) finput = "4th poly";
  else if(inputOption == 4) finput = "Nominal+Exp";
  
  TFile *wf = new TFile(Form("%s_fit_pt%.1f-%.1f_rap%.1f-%.1f_cent%d-%d_Gen1000000_input%d_useCentBkg%d_nToys%d_%d.root",fcoll.Data(),ptLow,ptHigh,yLow,yHigh,cLow,cHigh,inputOption,useCentIntBkgShape,nToys, rdmseed),"recreate");
  
  float massLow = 8. ;
  float massHigh = 14.;
  int   nMassBin  = (massHigh-massLow)*10;
  
  RooWorkspace *ws = new RooWorkspace("ws");
  RooWorkspace *wsinp = new RooWorkspace("wsinp");

  RooRealVar mass("mass","mass", massLow, massHigh);
  
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
          
  PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  initPset.SetMCSgl();

  RooRealVar sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma2s_1("sigma2s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma3s_1("sigma3s_1","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma1s_2("sigma1s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma2s_2("sigma2s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);
  RooRealVar sigma3s_2("sigma3s_2","width/sigma of the signal gaussian mass PDF",0.05, 0.05, 0.3);

  RooRealVar alpha1s_1("alpha1s_1","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha2s_1("alpha2s_1","tail shift", 5. , 1.5, 9.2);  
  RooRealVar alpha3s_1("alpha3s_1","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha1s_2("alpha1s_2","tail shift", 5. , 1.5, 9.8);
  RooRealVar alpha2s_2("alpha2s_2","tail shift", 5. , 1.5, 9.2);  
  RooRealVar alpha3s_2("alpha3s_2","tail shift", 5. , 1.5, 9.8);

  RooRealVar n1s_1("n1s_1","power order", 5. , 1.5, 9.8);
  RooRealVar n2s_1("n2s_1","power order", 5. , 1.5, 10.);
  RooRealVar n3s_1("n3s_1","power order", 4. , 1.5, 9.8);
  RooRealVar n1s_2("n1s_2","power order", 5. , 1.5, 9.8);
  RooRealVar n2s_2("n2s_2","power order", 5. , 1.5, 10.);
  RooRealVar n3s_2("n3s_2","power order", 4. , 1.5, 9.8);
  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);
  RooRealVar *f2s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);
  RooRealVar *f3s = new RooRealVar("f1s","1S CB fraction", 0.5, 0, 1);

  if ( initPset.n1s_1 == -1 )
  {
    cout << endl << endl << endl << "#########################  ERROR!!!! ##################" << endl;
    cout << "Fitting macro is stopped!" << endl << endl << endl;
    return;
  }
  else { 
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << endl << "Fixing the parameters..." << endl << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "initPset.n1s_1 = " << initPset.n1s_1 << endl;
    n1s_1.setVal(initPset.n1s_1);  n1s_1.setConstant(); 
    n2s_1.setVal(initPset.n2s_1);  n2s_1.setConstant();  
    n3s_1.setVal(initPset.n3s_1);  n3s_1.setConstant();  
    n1s_2.setVal(initPset.n1s_2);  n1s_2.setConstant(); 
    n2s_2.setVal(initPset.n2s_2);  n2s_2.setConstant();  
    n3s_2.setVal(initPset.n3s_2);  n3s_2.setConstant();  
    alpha1s_1.setVal(initPset.alpha1s_1);  alpha1s_1.setConstant();  
    alpha2s_1.setVal(initPset.alpha2s_1);  alpha2s_1.setConstant();  
    alpha3s_1.setVal(initPset.alpha3s_1);  alpha3s_1.setConstant();  
    alpha1s_2.setVal(initPset.alpha1s_2);  alpha1s_2.setConstant();  
    alpha2s_2.setVal(initPset.alpha2s_2);  alpha2s_2.setConstant();  
    alpha3s_2.setVal(initPset.alpha3s_2);  alpha3s_2.setConstant();  
    sigma1s_1.setVal(initPset.sigma1s_1);  sigma1s_1.setConstant();  
    sigma2s_1.setVal(initPset.sigma2s_1);  sigma2s_1.setConstant();  
    sigma3s_1.setVal(initPset.sigma3s_1);  sigma3s_1.setConstant();  
    sigma1s_2.setVal(initPset.sigma1s_2);  sigma1s_2.setConstant();  
    sigma2s_2.setVal(initPset.sigma2s_2);  sigma2s_2.setConstant();  
    sigma3s_2.setVal(initPset.sigma3s_2);  sigma3s_2.setConstant();  
    f1s->setVal(initPset.f1s);  f1s->setConstant();  
    f2s->setVal(initPset.f2s);  f2s->setConstant();  
    f3s->setVal(initPset.f3s);  f3s->setConstant();  
  } 

  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", mass, mean1s, sigma1s_1, alpha1s_1, n1s_1);
  cout << " n1s_1.getVal() = " << n1s_1.getVal() << endl;
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", mass, mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", mass, mean3s, sigma3s_1, alpha3s_1, n3s_1);
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", mass, mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", mass,mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", mass,mean3s, sigma3s_2, alpha3s_2, n3s_2);

  RooAddPdf*  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  RooAddPdf*  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1s) );
  RooAddPdf*  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*cb3s_2), RooArgList(*f1s) );


  //----------------------------------------------------------------------------------------
  //Generating function from nominal fit
  PSetUpsAndBkg bkgParm = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut) ; 
  bkgParm.SetToyMCParm();

  RooRealVar err_mu_gen("err_mu_gen","err_mu_gen",  bkgParm.bkg_mu_res) ;
  RooRealVar err_sigma_gen("err_sigma_gen","err_sigma_gen", bkgParm.bkg_sigma_res);
  RooRealVar m_decay_gen("err_lambda_gen","m_decay_gen", bkgParm.bkg_lambda_res);

  err_mu_gen.setVal(bkgParm.bkg_mu_res); err_mu_gen.setConstant();
  err_sigma_gen.setVal(bkgParm.bkg_sigma_res); err_sigma_gen.setConstant();
  m_decay_gen.setVal(bkgParm.bkg_lambda_res); m_decay_gen.setConstant();
  
  RooGenericPdf* bkgInp_gen;
  RooGenericPdf *bkgInp_in;
  if ( ptLow < 5)  { 
    bkgInp_in = new RooGenericPdf("bkgInp_gen","Background Gen","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu_gen,err_sigma_gen,m_decay_gen));
  }
  else {
    bkgInp_in = new RooGenericPdf("bkgInp_gen","Background Gen","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay_gen));
  }
  
  bkgInp_gen = bkgInp_in;
  
  float r1S_overTot = bkgParm.nSignal1s / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nSignal3s + bkgParm.nBkg ) ; // Numbers obtained from the real data
  float r2S_overTot = bkgParm.nSignal2s / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nSignal3s + bkgParm.nBkg ) ; 
  float r3S_overTot = bkgParm.nSignal3s / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nSignal3s + bkgParm.nBkg ) ; 
  float rBkg_overTot = bkgParm.nBkg / ( bkgParm.nSignal1s + bkgParm.nSignal2s + bkgParm.nSignal3s + bkgParm.nBkg ) ; 
  
  RooRealVar *nSig1sInp  = new RooRealVar("nSig1sInp","nSig1sInp", nGen * r1S_overTot,  0,   nGen);
  RooRealVar *nSig2sInp  = new RooRealVar("nSig2sInp","nSig2sInp", nGen * r2S_overTot, 0,   nGen);
  RooRealVar *nSig3sInp  = new RooRealVar("nSig3sInp","nSig3sInp", nGen * r3S_overTot, 0,   nGen);
  RooRealVar *nBkgInp  = new RooRealVar("nBkgInp","n_bkgInp",      nGen * rBkg_overTot,  0,   nGen);
  
  RooAddPdf* modelInput_gen; 
  modelInput_gen = new RooAddPdf("modelInput_gen","1S+2S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkgInp_gen),RooArgList(*nSig1sInp,*nSig2sInp,*nSig3sInp,*nBkgInp));
  
  TH1D *hNom_1S = new TH1D("hNom_1S",Form("SR Nominal, %d toys, %d events, cent %d-%d;N(1S) nom;Counts",nToys,nGen,cLow,cHigh),500,nGen * r1S_overTot*0.5 ,nGen * r1S_overTot*2);
  TH1D *hAlt_1S = new TH1D("hAlt_1S",Form("SR %s, %d toys, %d events, cent %d-%d;N(1S) alt;Counts",finput.Data(),nToys,nGen,cLow,cHigh),500,nGen * r1S_overTot*0.5, nGen * r1S_overTot*2);
  TH1D *hDev_1S = new TH1D("hDev_1S","Deviation;1S dev;Counts",2000,-100,100);

  TH1D *hNom_2S = new TH1D("hNom_2S",Form("SR Nominal, %d toys, %d events, cent %d-%d;N(2S) nom;Counts",nToys,nGen,cLow,cHigh),500,nGen * r2S_overTot*0.5,nGen * r2S_overTot*2);
  TH1D *hAlt_2S = new TH1D("hAlt_2S",Form("SR %s, %d toys, %d events, cent %d-%d;N(2S) alt;Counts",finput.Data(),nToys,nGen,cLow,cHigh),500,nGen * r2S_overTot*0.5,nGen * r2S_overTot*2);
  TH1D *hDev_2S = new TH1D("hDev_2S","Deviation;2S dev;Counts",2000,-100,100);

  TH1D *hNom_3S = new TH1D("hNom_3S",Form("SR Nominal, %d toys, %d events, cent %d-%d;N(3S) nom;Counts",nToys,nGen,cLow,cHigh),500,nGen * r3S_overTot*0.5,nGen * r3S_overTot * 2);
  TH1D *hAlt_3S = new TH1D("hAlt_3S",Form("SR %s, %d toys, %d events, cent %d-%d;N(3S) alt;Counts",finput.Data(),nToys,nGen,cLow,cHigh),500,nGen * r3S_overTot*0.5,nGen * r3S_overTot *2);
  TH1D *hDev_3S = new TH1D("hDev_3S","Deviation;3S dev;Counts",200000,-10000,10000);


  cout << "nSig1s : " << nGen * r1S_overTot << endl;
  cout << "nSig2s : " << nGen * r2S_overTot << endl;
  cout << "nSig3s : " << nGen * r3S_overTot << endl;
  cout << "nBkg : " << nGen * rBkg_overTot << endl;

  //----------------------------------------------------------------------------------------
 

  //Alternative Fit PDF 
  float the_ch4_k1 = bkgParm.ch4_k1 ;   float the_ch4_k2 = bkgParm.ch4_k2 ;  float the_ch4_k3 = bkgParm.ch4_k3 ; float the_ch4_k4 = bkgParm.ch4_k4 ;

  RooRealVar ch4_k1("pol4_k1","pol4_k1", 0.1 , -1.1, 1.1) ;
  RooRealVar ch4_k2("pol4_k2","pol4_k2", -0.1 , -1.1, 1.1) ;
  RooRealVar ch4_k3("pol4_k3","pol4_k3", -0.1 , -1.1, 1.1) ;
  RooRealVar ch4_k4("pol4_k4","pol4_k4", 0.1 , -1.1, 1.1) ;

  RooChebychev * bkgChPol4 = new RooChebychev("cPol4Bkg","Background4",mass,RooArgSet(ch4_k1,ch4_k2,ch4_k3,ch4_k4));  // if ( inputOption == kChPol3 )

  RooAddPdf* modelAltFit; 
  RooGenericPdf* bkgAltFit;
  bkgAltFit = (RooGenericPdf*) bkgChPol4;
  
  modelAltFit = new RooAddPdf("modelAltFit","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkgAltFit),RooArgList(*nSig1sInp,*nSig2sInp,*nSig3sInp,*nBkgInp));
  wsinp->import(*modelAltFit);
  
  initPset.SetMCBkg();
  double init_mu = initPset.bkg_mu ;
  double init_sigma = initPset.bkg_sigma ;
  double init_lambda = initPset.bkg_lambda ;

  RooRealVar err_mu("err_mu","err_mu", init_mu, 0, 30);
  RooRealVar err_sigma("err_sigma","err_sigma", init_sigma,0,30);
  RooRealVar m_decay("m_decay","m_decay", init_lambda, 0., 30.);
 /* 
  if( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4) && collId==kPPDATA) 
  {
    err_sigma.setVal(1.055); 
    err_sigma.setConstant();
  }
  if( ( ptLow == (float)0 ) && (ptHigh == (float)30 ) && (yLow == (float)0 ) && (yHigh == (float)2.4) && collId==kAADATA) 
  {
    err_sigma.setVal(1.103); 
    err_sigma.setConstant();
  }
*/
  RooGenericPdf *bkgFitOut;
  if ( ptLow < 5)  { 
    bkgFitOut = new RooGenericPdf("bkgFitOut","BackgroundOut","(TMath::Erf((@0-@1)/(TMath::Sqrt(2)*@2))+1)*0.5*TMath::Exp(-@0/@3)",RooArgList(mass,err_mu,err_sigma,m_decay));
  }
  else {
    bkgFitOut = new RooGenericPdf("bkgFitOut","BackgroundOut","TMath::Exp(-@0/@1)",RooArgList(mass,m_decay));
  }
  

  RooRealVar *nSig1sOut  = new RooRealVar("nSig1sOut","nSig1sOut", r1S_overTot*nGen, 0,  nGen);
  RooRealVar *nSig2sOut  = new RooRealVar("nSig2sOut","nSig2sOut", r2S_overTot*nGen, 0, nGen);
  RooRealVar *nSig3sOut  = new RooRealVar("nSig3sOut","nSig3sOut", r3S_overTot*nGen, 0, nGen);
  RooRealVar *nBkgOut  = new RooRealVar("nBkgOut","n_bkgOut",nGen * rBkg_overTot, 0, nGen);

  RooAddPdf*  cb1sOut = (RooAddPdf*)cb1s->Clone("cb1sOutput");
  RooAddPdf*  cb2sOut = (RooAddPdf*)cb2s->Clone("cb2sOutput");
  RooAddPdf*  cb3sOut = (RooAddPdf*)cb3s->Clone("cb3sOutput");
  RooAddPdf* modelOutput = new RooAddPdf("modelOutput","1S+2S+3S + Bkg",RooArgList(*cb1sOut, *cb2sOut, *cb3sOut,*bkgFitOut),RooArgList(*nSig1sOut,*nSig2sOut,*nSig3sOut,*nBkgOut));
  ws->import(*modelOutput);
  
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  //----****************--------for loop -----*******************-----------------
  
  for(int i=0;i<nToys;i++){
  
  ws->import(mass);
  wsinp->import(mass);
  mass.Print();
  
  Val_1S_nom=0;
  Val_1S_alt=0;
  Dev_1S=0;

  Val_2S_nom=0;
  Val_2S_alt=0;
  Dev_2S=0;

  Val_3S_nom=0;
  Val_3S_alt=0;
  Dev_3S=0;


  RooDataSet *data = modelInput_gen->generate(mass,nGen) ;

  RooPlot* xframe  = ws->var("mass")->frame(nMassBin); // bins
  xframe->SetXTitle("mass (Gev/c^{2})");
  xframe->GetXaxis()->CenterTitle();
  xframe->GetYaxis()->CenterTitle();
  RooPlot* xframe2 = (RooPlot*)xframe->Clone("xframe2");
  
  RooFitResult* fitResInput = wsinp->pdf("modelAltFit")->fitTo(*data,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  data->plotOn(xframe,Name("dataHist"),MarkerSize(0.7)) ;
  wsinp->pdf("modelAltFit")->plotOn(xframe, Name("inputModelHist"));
  wsinp->pdf("modelAltFit")->plotOn(xframe, Components(RooArgSet(*bkgAltFit)),LineColor(kBlack),LineStyle(kDashed));

  
  // New fit 
  
  RooFitResult* fitRes = ws->pdf("modelOutput")->fitTo(*data,Save(), Hesse(kTRUE),Range(massLow, massHigh),Minos(0), SumW2Error(kTRUE));
  data->plotOn(xframe2,Name("dataHist2"),MarkerSize(0.7)) ;
  ws->pdf("modelOutput")->plotOn(xframe2, Name("outputModelHist"));
  ws->pdf("modelOutput")->plotOn(xframe2, Components(RooArgSet(*bkgFitOut)),LineColor(kBlack),LineStyle(kDashed));
  
  Val_1S_nom = (float)(ws->var("nSig1sOut"))->getVal();
  Val_1S_alt = (float)(wsinp->var("nSig1sInp"))->getVal();
  Dev_1S = (1-Val_1S_alt/Val_1S_nom) * 100;

  Val_2S_nom = (float)(ws->var("nSig2sOut"))->getVal();
  Val_2S_alt = (float)(wsinp->var("nSig2sInp"))->getVal();
  Dev_2S = (1-Val_2S_alt/Val_2S_nom) * 100;

  Val_3S_nom = (float)(ws->var("nSig3sOut"))->getVal();
  Val_3S_alt = (float)(wsinp->var("nSig3sInp"))->getVal();
  Dev_3S = (1-Val_3S_alt/Val_3S_nom) * 100;

  // DRAW! 
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,800,400);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.49, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  pad1->SetBottomMargin(0); // Upper and lower plot are joined

  xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  drawText(Form("#Upsilon(1S) = %.5f",(float)(wsinp->var("nSig1sInp")->getVal())),0.2,0.75,1,16) ;
  
  if (inputOption==kChPol4 ) 
    drawText("4th order poly. Bkg.",0.2,0.82,2,15) ;
  
  if(collId == kAADATA)
    drawText("PbPb",0.3,0.45,1,15);
  if(collId == kPPDATA)
    drawText("pp", 0.3,0.45,1,15);

  drawText(Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV",ptLow,ptHigh ),0.31,0.60,1,12);
  drawText(Form("%.1f < y^{#mu#mu} < %.1f",yLow,yHigh ), 0.31,0.55,1,12);
  TString perc = "%";
  if(collId == kAADATA)
    drawText(Form("Cent %d-%d%s",cLow/2,cHigh/2,perc.Data()),0.31,0.5,4,12);
  
  TLatex *tex = new TLatex(0.4,0.88,"Toy MC generated");
  tex->SetTextFont(43);
  tex->SetTextSize(15);
  tex->SetNDC();
  //  tex->SetTextAngle(180);
  tex->Draw();

  RooArgList paramListinp = fitResInput->floatParsFinal();
  paramListinp.Print("v");
  RooPlot* legFrameinp = wsinp->var("mass")->frame(Name("Fit Results"), Title("Fit Results"));
  wsinp->pdf("modelAltFit")->paramOn(legFrameinp,Layout(.6,.9, .5),Parameters(paramListinp));
  legFrameinp->getAttText()->SetTextAlign(11);
  legFrameinp->getAttText()->SetTextSize(0.028);
  TPaveText* hhinp = (TPaveText*)legFrameinp->findObject(Form("%s_paramBox",wsinp->pdf("modelAltFit")->GetName()));
  hhinp->SetY1(0.35); hhinp->SetY2(0.83);
  hhinp->Draw();
  // PULL

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.49, 0.25);
  c1->cd();
  pad2->Draw();
  pad2->cd();
  RooHist* hpull = xframe->pullHist("dataHist","inputModelHist");
  RooPlot* pullFrame = wsinp->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(2.57);
  pullFrame->GetYaxis()->SetTitleOffset(1.8) ;
  pullFrame->GetYaxis()->SetLabelSize(0.16) ;
  pullFrame->GetYaxis()->SetRange(-10,10) ;
  pullFrame->GetXaxis()->SetTitleOffset(0.7) ;
  pullFrame->GetXaxis()->SetLabelSize(0.1) ;
  pullFrame->GetXaxis()->SetTitleSize(0.13) ;
  pullFrame->Draw() ;

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

  int numFitPar = fitRes->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;


  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  drawText(Form("chi^{2}/ndf : %.3f / %d ",chisq,ndf ),0.15,0.95,1,12);

  TPad *pad3 = new TPad("pad3", "pad3", 0.51, 0.25, 0.99, 1);
  pad3->SetTicks(1,1);
  pad3->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();
  pad3->Draw(); pad3->cd();

  xframe2->GetYaxis()->SetTitleOffset(1.4) ; xframe2->Draw() ;
  TLatex *tex2 = new TLatex(0.4,0.9,"Fitted by Nominal function");
  tex2->SetTextFont(43);
  tex2->SetTextSize(15);
  tex2->SetTextColor(2);
  tex2->SetNDC();
  tex2->Draw();
  drawText(Form("#Upsilon(2S)/#Upsilon(1S) = %.5f",(float)(ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal())), 0.4,0.85,1,16 );

  // *~*~*~*~*~*~*~* Draw the parameters in the plot  *~*~*~*~*~*~*~* //
  RooArgList paramList = fitRes->floatParsFinal();
  paramList.Print("v");
  RooPlot* legFrame = ws->var("mass")->frame(Name("Fit Results"), Title("Fit Results"));
  ws->pdf("modelOutput")->paramOn(legFrame,Layout(.6,.9, .5),Parameters(paramList));
  legFrame->getAttText()->SetTextAlign(11);
  legFrame->getAttText()->SetTextSize(0.028);
  TPaveText* hh = (TPaveText*)legFrame->findObject(Form("%s_paramBox",ws->pdf("modelOutput")->GetName()));
  hh->SetY1(0.35); hh->SetY2(0.83);
  hh->Draw();

  TPad *pad4 = new TPad("pad4", "pad4", 0.51, 0.05, 0.99, 0.25);
  // pad4->SetBottomMargin(0); // Upper and lower plot are joined
  c1->cd();
  pad4->Draw();
  pad4->cd();
  RooHist* hpullOut = xframe2->pullHist("dataHist2","outputModelHist");
  RooPlot* pullOutFrm = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullOutFrm->addPlotable(hpullOut,"P") ;
  pullOutFrm->SetTitleSize(2.57);
  pullOutFrm->GetYaxis()->SetTitleOffset(1.8) ;
  pullOutFrm->GetYaxis()->SetLabelSize(0.16) ;
  pullOutFrm->GetYaxis()->SetRange(-10,10) ;
  pullOutFrm->GetXaxis()->SetTitleOffset(0.7) ;
  pullOutFrm->GetXaxis()->SetLabelSize(0.1) ;
  pullOutFrm->GetXaxis()->SetTitleSize(0.13) ;
  pullOutFrm->Draw() ;

  hNom_1S->Fill(Val_1S_nom);
  hAlt_1S->Fill(Val_1S_alt);
  hDev_1S->Fill(Dev_1S);
  
  hNom_2S->Fill(Val_2S_nom);
  hAlt_2S->Fill(Val_2S_alt);
  hDev_2S->Fill(Dev_2S);
  
  hNom_3S->Fill(Val_3S_nom);
  hAlt_3S->Fill(Val_3S_alt);
  hDev_3S->Fill(Dev_3S);
  
  //c1->SaveAs(Form( "toyMCFit_collId%d_pt%.0f-%.0fGeV_y%.0f-%.0f_cBin%d-%d_muPtCut%.0fGeV_BkgPDFOpt%d_nGen%d_useCentIntBkgShape%d_%d.png", collId, ptLow, ptHigh, yLow*10, yHigh*10, cLow, cHigh, muPtCut, inputOption, nGen,useCentIntBkgShape, rdmseed) );
  if((double)chisq/ndf>5) 
  {
     c1->SaveAs(Form( "toyMCFit_collId%d_pt%.0f-%.0fGeV_y%.0f-%.0f_cBin%d-%d_muPtCut%.0fGeV_BkgPDFOpt%d_nGen%d_useCentIntBkgShape%d_%d.png", collId, ptLow, ptHigh, yLow*10, yHigh*10, cLow, cHigh, muPtCut, inputOption, nGen,useCentIntBkgShape, rdmseed) );
    continue;
  }
  


  // *~*~*~*~*~*~*~* Print the results *~*~*~*~*~*~*~* //

  //cout << "nSig2sInp/nSig1sInp = " << nSig2sInp->getVal() / nSig1sInp->getVal() << endl;
  cout << "input fit ratio = " << wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal() << endl;
  cout << "output fit ratio = " << ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal() << endl;
  
  float r1 =  wsinp->var("nSig2sInp")->getVal() / wsinp->var("nSig1sInp")->getVal() ; 
  float r2 =  ws->var("nSig2sOut")->getVal() / ws->var("nSig1sOut")->getVal() ; 
  cout << Form( "collId: %d,    pt: %.0f - %.0fGeV,   y: %.1f - %.1f,  cBin: %d - %d", collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh ) << endl;
  cout << "Uncertainty = "  << (r2 - r1 ) / r1 << endl;
//  }  
  
  }

  wf->cd();
  hNom_1S->Write();
  hAlt_1S->Write();
  hDev_1S->Write();

  hNom_2S->Write();
  hAlt_2S->Write();
  hDev_2S->Write();

  hNom_3S->Write();
  hAlt_3S->Write();
  hDev_3S->Write();

   
} 

