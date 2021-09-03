//This code fits the upsilon data with either the nominal fit or an alternative fit. The difference between the two fits is the signal shape. The nominal fit fits the signals with double CB functions, while the alternative fit fits them with just a gaussian.

#include <iostream>
#include "../HeaderFiles/rootFitHeaders.h"
#include "../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../HeaderFiles/cutsAndBin.h"
#include "../HeaderFiles/PsetCollection.h"
#include "../HeaderFiles/CMS_lumi.C"
#include "../HeaderFiles/tdrstyle.C"
#include "../HeaderFiles/StyleSetting.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"


using namespace std;
using namespace RooFit;
void SimultaneousFitEverything( 
       int collId = kPADATA,  
       float yLow=-1.93, float yHigh=1.93,//Run 1 has p going in -z direction
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       bool whichModel=0,   // Nominal = 0. Alternative = 1.
       int ICset = 1
			) 
{


  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  float eta_low = -2.4;
  float eta_high = 2.4;
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8;
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  TFile* f1;
  TFile* f2;
  float yLowLab;
  float yHighLab;

  //The order is {sigma1s_1,x1s,alpha1s_1,n1s_1,f1s,err_mu,err_sigma,m_lambda}
  double paramsupper[8] = {0.2, 1.0, 5.0, 5.0, 1.0, 15.0, 15.0, 25.0};
  double paramslower[8] = {0.02, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double sigma1s_1_init = 0.1;
  double x1s_init = 0.6;
  double alpha1s_1_init = 1.5;
  double n1s_1_init = 2.0;
  double f1s_init = 0.5;
  if (ICset>1 && ICset<4) {
    sigma1s_1_init = 0.3;
    x1s_init = 0.3;
    alpha1s_1_init = 2.6;
    n1s_1_init = 3.0;
    f1s_init = 0.1;
  }
  double err_sigma_init = 1;
  double err_mu_init = 8;
  double m_lambda_init = 8;
  if (ICset>2) {
    err_mu_init = 5;
    m_lambda_init = 5;
  }

  float ptLow[3], ptHigh[3];
  float ptLow[0] = 0;
  float ptHigh[0] = 6;
  float ptLow[1] = 6;
  float ptHigh[1] = 30;
  float ptLow[2] = 0;
  float ptHigh[2] = 30;
  //Select Data Set
  if (collId==kPADATA) {
    f1 = new TFile("../yskimPA1st_OpSign_20177262037_unIdentified.root");
    f2 = new TFile("../yskimPA2nd_OpSign_20177262044_unIdentified.root");
    yLowLab = yLow+0.47;
    yHighLab = yHigh+0.47;
  }
  else if (collId==kPPDATA) {
    f1 = new TFile("../yskimPP_L1DoubleMu0PD_Trig-L1DoubleMu0_OpSign_20177262158_.root");
    yLowLab = yLow;
    yHighLab = yHigh;
  }

  //import and merge datasets
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  if (collId==kPADATA) {
    RooDataSet *dataset2 = (RooDataSet*)f2->Get("dataset");
    dataset->append(*dataset2);
  }
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS[3];
  for (int idataset=0; idataset<3; idataset++) {
    TString kineCut;
    if (collId==kPADATA) {
      kineCut = Form("pt>%.2f && pt<%.2f && y>%.2f && y<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow[idataset], ptHigh[idataset], yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    }
    else if (collId==kPPDATA) {
      kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f && eta1<%.2f && eta1>%.2f && eta2<%.2f && eta2>%.2f",ptLow[idataset], ptHigh[idataset], yLowLab, yHighLab, eta_high,eta_low, eta_high,eta_low );
    }
    if (muPtCut>0){
      kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) ", muPtCut, muPtCut);
    }
    cout << "Creating dataset " << idataset << endl;
    reducedDS[idataset] = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
    cout << "Setting dataset name" << endl;
    reducedDS[idataset]->SetName(Form("reducedDS_%i",idataset));
    cout << "Importing dataset " << idataset << endl;
    ws->import(*reducedDS[idataset]);
    cout << "Imported dataset " << idataset << endl;
  }
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  cout << "Loaded all the datasets!!!!!!!!!!!!!!!!!!!!" << endl;
  RooCategory tp("tp","tp");
  tp.defineType("1");
  tp.defineType("2");
  tp.defineType("3");

  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );

  // Create a dataset that imports contents of all the above datasets mapped by index category tp
  RooDataSet* dsABC = new RooDataSet("dsABC","dsABC",RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))),Index(tp),Import("1",*reducedDS[0]),Import("2",*reducedDS[1]),Import("3",*reducedDS[2]));
  cout << "******** New Combined Dataset ***********" << endl;
  dsABC->Print();
  ws->import(*dsABC);

  TCanvas* c[3];
  RooPlot* myPlot[3];
  RooRealVar mean1s[3];
  RooFormulaVar mean2s[3];
  RooFormulaVar mean3s[3];
  RooRealVar sigma1s_1[3];
  RooFormulaVar sigma2s_1[3];
  RooFormulaVar sigma3s_1[3];
  RooRealVar* x1s[3];
  RooFormulaVar sigma1s_2[3];
  RooFormulaVar sigma2s_2[3];
  RooFormulaVar sigma3s_2[3];
  RooRealVar alpha1s_1[3];
  RooFormulaVar alpha2s_1[3];
  RooFormulaVar alpha3s_1[3];
  RooFormulaVar alpha1s_2[3];
  RooFormulaVar alpha2s_2[3];
  RooFormulaVar alpha3s_2[3];
  RooRealVar n1s_1[3];
  RooFormulaVar n2s_1[3];
  RooFormulaVar n3s_1[3];
  RooFormulaVar n1s_2[3];
  RooFormulaVar n2s_2[3];
  RooFormulaVar n3s_2[3];
  RooRealVar *f1s[3];
  RooFormulaVar f2s[3];
  RooFormulaVar f3s[3];
  RooCBShape* cb1s_1[3];
  RooCBShape* cb2s_1[3];
  RooCBShape* cb3s_1[3];
  RooCBShape* cb1s_2[3];
  RooCBShape* cb2s_2[3];
  RooCBShape* cb3s_2[3];
  RooAddPdf* cb1s[3];
  RooAddPdf* cb2s[3];
  RooAddPdf* cb3s[3];
  RooRealVar *nSig1sFree[3];
  RooRealVar *nSig2sFree[3];
  RooRealVar *nSig3sFree[3];
  RooRealVar err_mu[3];
  RooRealVar err_sigma[3];
  RooRealVar m_lambda[3];
  RooGenericPdf *bkg[3];
  RooGenericPdf *bkgLowPt[3];
  RooGenericPdf *bkgHighPt[3];
  RooRealVar *nBkgFree[3];
  RooFormulaVar *nSig1s[3];
  RooFormulaVar *nSig2s[3];
  RooFormulaVar *nSig3s[3];
  RooFormulaVar *nBkg[3];
  RooAddPdf* model[3];
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",tp);
  for (int i=0; i<3; i++) {
    cout << "making canvas " << i << endl;
    c[i] = new TCanvas(Form("canvas[%i]",i),"My plots",4,45,550,520);
    c[i]->cd();
    cout << "making plot " << i << endl;
    myPlot[i] = ws->var("mass")->frame(nMassBin); // bins
    cout << "Plotting data " << i << endl;
    ws->data(Form("reducedDS_%i",i))->plotOn(myPlot[i],Name(Form("dataHist[%i]",i)));

    //SIGNAL:
    mean1s[i] = RooRealVar(Form("m_{#Upsilon(1S)}[%i]",i),"mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
    mean2s[i] = RooFormulaVar(Form("mean2s[%i]",i),Form("m_{#Upsilon(1S)}[%i]*mRatio21",i), RooArgSet(mean1s[i],mRatio21) );
    mean3s[i] = RooFormulaVar(Form("mean3s[%i]",i),Form("m_{#Upsilon(1S)}[%i]*mRatio31",i), RooArgSet(mean1s[i],mRatio31) );

    sigma1s_1[i] = RooRealVar(Form("sigma1s_1[%i]",i),"width/sigma of the signal gaussian mass PDF",sigma1s_1_init, paramslower[0], paramsupper[0]);
    sigma2s_1[i] = RooFormulaVar(Form("sigma2s_1[%i]",i),"@0*@1",RooArgList(sigma1s_1[i],mRatio21) );
    sigma3s_1[i] = RooFormulaVar(Form("sigma3s_1[%i]",i),"@0*@1",RooArgList(sigma1s_1[i],mRatio31) );

    x1s[i] = new RooRealVar(Form("x1s[%i]",i),"sigma ratio ", x1s_init, paramslower[1], paramsupper[1]);

    sigma1s_2[i] = RooFormulaVar(Form("sigma1s_2[%i]",i),"@0*@1",RooArgList(sigma1s_1[i], *x1s[i]) );
    sigma2s_2[i] = RooFormulaVar(Form("sigma2s_2[%i]",i),"@0*@1",RooArgList(sigma1s_2[i],mRatio21) );
    sigma3s_2[i] = RooFormulaVar(Form("sigma3s_2[%i]",i),"@0*@1",RooArgList(sigma1s_2[i],mRatio31) );

    alpha1s_1[i] = RooRealVar(Form("alpha1s_1[%i]",i),"tail shift", alpha1s_1_init, paramslower[2], paramsupper[2]);
    alpha2s_1[i] = RooFormulaVar(Form("alpha2s_1[%i]",i),"1.0*@0",RooArgList(alpha1s_1[i]) );
    alpha3s_1[i] = RooFormulaVar(Form("alpha3s_1[%i]",i),"1.0*@0",RooArgList(alpha1s_1[i]) );
    alpha1s_2[i] = RooFormulaVar(Form("alpha1s_2[%i]",i),"1.0*@0",RooArgList(alpha1s_1[i]) );
    alpha2s_2[i] = RooFormulaVar(Form("alpha2s_2[%i]",i),"1.0*@0",RooArgList(alpha1s_1[i]) );
    alpha3s_2[i] = RooFormulaVar(Form("alpha3s_2[%i]",i),"1.0*@0",RooArgList(alpha1s_1[i]) );

    n1s_1[i] = RooRealVar(Form("n1s_1[%i]",i),"power order", n1s_1_init , paramslower[3], paramsupper[3]);
    n2s_1[i] = RooFormulaVar(Form("n2s_1[%i]",i),"1.0*@0",RooArgList(n1s_1[i]) );
    n3s_1[i] = RooFormulaVar(Form("n3s_1[%i]",i),"1.0*@0",RooArgList(n1s_1[i]) );
    n1s_2[i] = RooFormulaVar(Form("n1s_2[%i]",i),"1.0*@0",RooArgList(n1s_1[i]) );
    n2s_2[i] = RooFormulaVar(Form("n2s_2[%i]",i),"1.0*@0",RooArgList(n1s_1[i]) );
    n3s_2[i] = RooFormulaVar(Form("n3s_2[%i]",i),"1.0*@0",RooArgList(n1s_1[i]) );

    f1s[i] = new RooRealVar(Form("f1s[%i]",i),"1S CB fraction", f1s_init, paramslower[4], paramsupper[4]);
    f2s[i] = RooFormulaVar(Form("f2s[%i]",i),"1.0*@0",RooArgList(*f1s[i]) );
    f3s[i] = RooFormulaVar(Form("f3s[%i]",i),"1.0*@0",RooArgList(*f1s[i]) );

    // Set up crystal ball shapes
    cb1s_1[i] = new RooCBShape(Form("cball1s_1[%i]",i), "cystal Ball", *(ws->var("mass")), mean1s[i], sigma1s_1[i], alpha1s_1[i], n1s_1[i]);
    cb2s_1[i] = new RooCBShape(Form("cball2s_1[%i]",i), "cystal Ball", *(ws->var("mass")), mean2s[i], sigma2s_1[i], alpha2s_1[i], n2s_1[i]);
    cb3s_1[i] = new RooCBShape(Form("cball3s_1[%i]",i), "cystal Ball", *(ws->var("mass")), mean3s[i], sigma3s_1[i], alpha3s_1[i], n3s_1[i]);

  //DOUBLE CRYSTAL BALL
    cb1s_2[i] = new RooCBShape(Form("cball1s_2[%i]",i), "cystal Ball", *(ws->var("mass")), mean1s[i], sigma1s_2[i], alpha1s_2[i], n1s_2[i]);
    cb2s_2[i] = new RooCBShape(Form("cball2s_2[%i]",i), "cystal Ball", *(ws->var("mass")), mean2s[i], sigma2s_2[i], alpha2s_2[i], n2s_2[i]);
    cb3s_2[i] = new RooCBShape(Form("cball3s_2[%i]",i), "cystal Ball", *(ws->var("mass")), mean3s[i], sigma3s_2[i], alpha3s_2[i], n3s_2[i]);
    cb1s[i] = new RooAddPdf(Form("cb1s[%i]",i),"Signal 1S",RooArgList(*cb1s_1[i],*cb1s_2[i]), RooArgList(*f1s[i]) );
    cb2s[i] = new RooAddPdf(Form("cb2s[%i]",i),"Signal 2S",RooArgList(*cb2s_1[i],*cb2s_2[i]), RooArgList(*f1s[i]) );
    cb3s[i] = new RooAddPdf(Form("cb3s[%i]",i),"Signal 3S",RooArgList(*cb3s_1[i],*cb3s_2[i]), RooArgList(*f1s[i]) );

    nSig1s[i]= new RooRealVar(Form("nSig1sFree[%i]",i)," 1S signals",0,1000000);
    nSig2s[i]= new RooRealVar(Form("nSig2sFree[%i]",i)," 2S signals",-20,360000);
    nSig3s[i]= new RooRealVar(Form("nSig3sFree[%i]",i)," 3S signals",-50,260000);

    //BACKGROUND
    err_mu[i] = RooRealVar(Form("#mu[%i]",i),"err_mu", err_mu_init,  paramslower[5], paramsupper[5]) ;
    err_sigma[i] = RooRealVar(Form("#sigma[%i]",i),"err_sigma", err_sigma_init, paramslower[6], paramsupper[6]);
    m_lambda[i] = RooRealVar(Form("#lambda[%i]",i),"m_lambda",  m_lambda_init, paramslower[7], paramsupper[7]);

    bkgLowPt[i] = new RooGenericPdf(Form("bkgLowPt[%i]",i),"Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda[i], err_mu[i], err_sigma[i]) );

    //THIS IS THE HIGH-PT BACKGROUND FUNCTION
    bkgHighPt[i] = new RooGenericPdf(Form("bkgHighPt[%i]",i),"Background","TMath::Exp(-@0/@1)",RooArgList(*(ws->var("mass")),m_lambda[i]));
    if  (ptLow[i] >= 5)        bkg[i] = bkgHighPt[i] ;
    else bkg[i] = bkgLowPt[i];

    nBkgFree[i] = new RooRealVar(Form("nBkgFree[%i]",i),"fraction of component 1 in bkg",10000,0,5000000); 

    //Constraints
    nSig1s[i] = new RooFormulaVar(Form("nSig1s[%i]",i)," 1S signals","@0",RooArgList(*nSig1sFree[i]) );
    nSig2s[i] = new RooFormulaVar(Form("nSig2s[%i]",i)," 2S signals","@0",RooArgList(*nSig2sFree[i]) );
    nSig3s[i] = new RooFormulaVar(Form("nSig3s[%i]",i)," 3S signals","@0",RooArgList(*nSig3sFree[i]) );
    nBkg[i] = new RooFormulaVar(Form("nBkg[%i]",i)," 3S signals","@0",RooArgList(*nBkgFree[i]) );

  //Build the model
    model[i] = new RooAddPdf(Form("model[%i]",i),"1S+2S+3S + Bkg",RooArgList(*cb1s[i], *cb2s[i], *cb3s[i], *bkg[i]),RooArgList(*nSig1s[i],*nSig2s[i],*nSig3s[i],*nBkg[i]));
    ws->import(*model[i]);
    c1[i]->cd();
    RooPlot* myPlot2[i] = (RooPlot*)myPlot[i]->Clone();
    ws->data(Form("reducedDS[%i]",i))->plotOn(myPlot2[i],Name(Form("dataOS_FIT[%i]",i)),MarkerSize(.8));

    // Construct simultaneous PDF
    simPdf->addPdf(*(ws->pdf(Form("model[%i]",i))),Form("%i",i)) ;
  }

  ws->import(*simPdf);

  cout << endl << "********* Starting Simutaneous Fit **************" << endl << endl;
  RooFitResult* fitResSim = ws->pdf("simPdf")->fitTo(*dsABC,Save(), Hesse(kTRUE),Range(massLow,massHigh),Timer(kTRUE),Extended(kTRUE));
  cout << endl << "********* Finished Simutaneous Fit **************" << endl << endl;
  ws->import(*fitResSim);

  c1[1]->cd();
  ws->pdf("model[1]")->plotOn(myPlot2[1],Name("modelHist[1]"));
  ws->pdf("model[1]")->plotOn(myPlot2[1],Name("Sig1S[1]"),Components(RooArgSet(*cb1s[1])),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model[1]")->plotOn(myPlot2[1],Components(RooArgSet(*cb2s[1])),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model[1]")->plotOn(myPlot2[1],Components(RooArgSet(*cb3s[1])),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model[1]")->plotOn(myPlot2[1],Name("bkgPDF[1]"),Components(RooArgSet(*bkg[1])),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  //make a pretty plot
  myPlot2[1]->SetFillStyle(4000);
  myPlot2[1]->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2[1]->GetYaxis()->SetTitleOffset(1.43);
  myPlot2[1]->GetYaxis()->CenterTitle();
  //myPlot2[1]->GetYaxis()->SetTitleSize(0.058);
  myPlot2[1]->GetYaxis()->SetLabelSize(0.054);
  //myPlot2[1]->GetXaxis()->SetLabelSize(0);
  myPlot2[1]->GetXaxis()->SetRangeUser(8,14);
  //myPlot2[1]->GetXaxis()->SetTitleSize(0);
  myPlot2[1]->Draw();

  c[1]->Update();

  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString outFileName;
  if (whichModel){
    outFileName = Form("SimultaneousFits/sim_ptSplit[1]ltfitresults_upsilon_%s.root",kineLabel.Data());
  }
  else {
    outFileName = Form("SimultaneousFits/sim_ptSplit_nomfitresults_upsilon_%s.root",kineLabel.Data());
  }
  TFile* outf = new TFile(outFileName,"recreate");
  ws->Write();
  outf->Close();

  delete c[1];
  delete x1s[1];
  delete f1s[1];
  delete cb1s_1[1];
  delete cb1s_2[1];
  delete cb2s_1[1];
  delete cb2s_2[1];
  delete cb3s_1[1];
  delete cb3s_2[1];
  delete cb1s[1];
  delete cb2s[1];
  delete cb3s[1];
  delete nSig1s[1];
  delete nSig2s[1];
  delete nSig3s[1];
  delete nBkg[1];
  delete bkgLowPt[1];
  delete bkgHighPt[1];
  delete model[1];

  delete c_B;
  delete x1s_B;
  delete f1s_B;
  delete cb1s_1_B;
  delete cb1s_2_B;
  delete cb2s_1_B;
  delete cb2s_2_B;
  delete cb3s_1_B;
  delete cb3s_2_B;
  delete cb1s_B;
  delete cb2s_B;
  delete cb3s_B;
  delete nSig1s_B;
  delete nSig2s_B;
  delete nSig3s_B;
  delete nBkg_B;
  delete bkgLowPt_B;
  delete bkgHighPt_B;
  delete model_B;

  delete c_C;
  delete x1s_C;
  delete f1s_C;
  delete cb1s_1_C;
  delete cb1s_2_C;
  delete cb2s_1_C;
  delete cb2s_2_C;
  delete cb3s_1_C;
  delete cb3s_2_C;
  delete cb1s_C;
  delete cb2s_C;
  delete cb3s_C;
  delete nSig1s_C;
  delete nSig2s_C;
  delete nSig3s_C;
  delete nBkg_C;
  delete bkgLowPt_C;
  delete bkgHighPt_C;
  delete model_C;

  delete f1;
  if (collId==kPADATA) delete f2;
  delete ws;
  delete simPdf;
  delete dsABC;
  delete outf;

} 
 
