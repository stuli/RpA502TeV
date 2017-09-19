//Author: Jared Jay
//This code is based on the fitting code in simple_Jared.C, which is a modified version of simple_Santona.C.
//This code generates pseudo-data from the nominal fit, and then fits it with the nominal fit and the new alternative fit. It then compares the estimated upsilon yields from each fit and plots the percent difference for each of the upsilons 1S, 2S, and 3S. It repeats this process as many times as you want. It also plots the most recent fit, and all the plots are updated after every fit.

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



using namespace std;
using namespace RooFit;
void Fit_Fake_Data( 
       int collId = kPADATA,  
       float ptLow=0, float ptHigh=5, 
       float yLow=-1.93, float yHigh=1.93,
       int cLow=0, int cHigh=200,
       float muPtCut=4.0,
       const int numtrials = 100
			) 
{
  gROOT->ProcessLine(".L RooMyPdf.cxx+");

  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  

  gStyle->SetEndErrorSize(0);
 
  TString SignalCB = "Double";

  float massLow = 8; 
  float massHigh = 14;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;

  int   nMassBin  = (massHigh-massLow)*10;

  //import generating model
  const char *inputFile="oldBkgModel.root";
  TFile *thefile = new TFile(inputFile);
  RooAddPdf* genModel = (RooAddPdf*)thefile->Get("model;1");
  RooWorkspace *wsgen = new RooWorkspace("workspace");
  wsgen->import(*genModel);

  //create histograms
  TH1D* histo1s = new TH1D("1S %Diff in Yield","1S %Diff in Yield",100,-100,100);
  TH1D* histo2s = new TH1D("2S %Diff in Yield","2S %Diff in Yield",100,-100,100);
  TH1D* histo3s = new TH1D("3S %Diff in Yield","3S %Diff in Yield",100,-100,100);
  TCanvas* cyields =  new TCanvas("canvas3","results",4,45,1100,400);
  cyields->Divide(3,1);

  gStyle->SetStatW(0.4);// Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.2);// Set height of stat-box (fraction of pad size)

  cyields->cd(1);
  histo1s->SetXTitle("%Diff");
  histo1s->GetXaxis()->SetTitleSize(0.05);
  histo1s->GetYaxis()->SetLabelSize(0.05);
  histo1s->GetXaxis()->SetLabelSize(0.05);
  histo1s->SetStats(kTRUE);
  histo1s->Draw();

  cyields->cd(2);
  histo2s->SetXTitle("%Diff");
  histo2s->GetXaxis()->SetTitleSize(0.05);
  histo2s->GetYaxis()->SetLabelSize(0.05);
  histo2s->GetXaxis()->SetLabelSize(0.05);
  histo2s->Draw();

  cyields->cd(3);
  histo3s->SetXTitle("%Diff");
  histo3s->GetXaxis()->SetTitleSize(0.05);
  histo3s->GetYaxis()->SetLabelSize(0.05);
  histo3s->GetXaxis()->SetLabelSize(0.05);
  histo3s->Draw();

  TCanvas* cfit =  new TCanvas("canvas2","fitted data",600,500,550,520);

  double chisqtest[2] = {0};
  TNtuple* ntuple = new TNtuple("ntuple","Data from fits","chisqnom:chisqalt:diff1s:diff2s:diff3s",numtrials);

for (int itrial = 0; itrial<numtrials; itrial++) {

  double yield1s[2] = {0};
  double yield2s[2] = {0};
  double yield3s[2] = {0};

  cout << "Starting trial " << itrial+1 << " of " << numtrials << endl;

  //Generate fake data from the model
  //The real data set had 12954 events.
  RooDataSet* reducedDS = genModel->generate(*(wsgen->var("mass")));
  reducedDS->SetName("reducedDS");

//loop through the two fitting models
for (int imodel = 0; imodel<=1; imodel++){

  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*reducedDS);

  cfit->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  
  RooPlot* myPlot = ws->var("mass")->frame(nMassBin);
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
          
  PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, muPtCut); //, muyCut) ; 
  initPset.SetMCSgl();


  //Set initial parameters
  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass PDF",5.6329e-02, 0.02, 0.3);
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", 1.7898, 0, 2.4);

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );

  RooRealVar alpha1s_1("alpha1s_1","tail shift", 2.1849 , 1.429, 3.321);
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", 1.4176 , 1.416, 3.357);
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 1.9832e-01, 0, 1);
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(*f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(*f1s) );

  //fix upsilon signal shape parameters
  sigma1s_1.setConstant(kTRUE);
  x1s->setConstant(kTRUE);
  alpha1s_1.setConstant(kTRUE);
  n1s_1.setConstant(kTRUE);
  f1s->setConstant(kTRUE);
  
  // Set initial parameters
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb2s_1 = new RooCBShape("cball2s_1", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_1, alpha2s_1, n2s_1);
  RooCBShape* cb3s_1 = new RooCBShape("cball3s_1", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_1, alpha3s_1, n3s_1);
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);
  RooCBShape* cb2s_2 = new RooCBShape("cball2s_2", "cystal Ball", *(ws->var("mass")), mean2s, sigma2s_2, alpha2s_2, n2s_2);
  RooCBShape* cb3s_2 = new RooCBShape("cball3s_2", "cystal Ball", *(ws->var("mass")), mean3s, sigma3s_2, alpha3s_2, n3s_2);

  RooAddPdf*  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  RooAddPdf*  cb2s = new RooAddPdf("cb2s","Signal 2S",RooArgList(*cb2s_1,*cb2s_2), RooArgList(*f1s) );
  RooAddPdf*  cb3s = new RooAddPdf("cb3s","Signal 3S",RooArgList(*cb3s_1,*cb3s_2), RooArgList(*f1s) );

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",0,1000000);
  RooRealVar *nSig2s= new RooRealVar("nSig2s"," 2S signals",-20,360000);
  RooRealVar *nSig3s= new RooRealVar("nSig3s"," 3S signals",-50,260000);


  // background : 
  initPset.SetMCBkg();
  double init_mu = initPset.bkg_mu ;
  double init_sigma = initPset.bkg_sigma ;
  double init_lambda = initPset.bkg_lambda ;

  double init_mu_min = init_mu - 10; double init_mu_max = init_mu + 10;
  double init_sigma_min = init_sigma - 10.; double init_sigma_max = init_sigma + 10;
  double init_lambda_min = init_lambda - 10; double init_lambda_max = init_lambda + 10;
  if(init_mu_min <0) init_mu_min = 0;
  if(init_sigma_min <0) init_sigma_min = 0;
  if(init_lambda_min <0) init_lambda_min = 0;
 
  RooRealVar err_mu("#mu","err_mu", 8,  0, 25) ;
  RooRealVar err_sigma("#sigma","err_sigma", 8, 0,25);
  RooRealVar m_lambda("#lambda","m_lambda",  8, 0,25);
  
  //THIS IS THE NEW LOW-PT BACKGROUND FUNCTION
  if (imodel>0){
    RooRealVar a1("A1","A1",100,0,1000);
    RooRealVar a2("A2","A2",100,0,1000);
    RooRealVar a3("A3","A3",100,0,1000);
    RooRealVar a4("A4","A4",100,0,1000);
    RooRealVar a5("A5","A5",100,0,1000);
    RooMyPdf *bkgLowPtnew = new RooMyPdf("bkgLowPt","Background",*(ws->var("mass")),a1,a2,a3,a4,a5);
  }
  else {
  //THIS IS THE OLD LOW-PT BACKGROUND FUNCTION
    RooGenericPdf *bkgLowPt = new RooGenericPdf("bkgLowPt","Background","TMath::Exp(-@0/@1)*(TMath::Erf((@0-@2)/(TMath::Sqrt(2)*@3))+1)*0.5",RooArgList( *(ws->var("mass")), m_lambda, err_mu, err_sigma) );
  }

  RooRealVar *nBkg = new RooRealVar("nBkg","fraction of component 1 in bkg",10000,0,5000000);  

  //Build the model
  RooAddPdf* model = new RooAddPdf();
  if (imodel>0) {
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkgLowPtnew),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  }
  else {
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *cb2s, *cb3s, *bkgLowPt),RooArgList(*nSig1s,*nSig2s,*nSig3s,*nBkg));
  }

  ws->import(*model);

  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));

  //Fit the model to the data
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb2s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Components(RooArgSet(*cb3s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  if (imodel>0) {
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkgLowPtnew)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  }
  else {
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkgLowPt)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));
  }

  //make a pretty plot
  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  myPlot2->GetYaxis()->SetLabelSize(0.054);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  float pos_text_x = 0.43;
  float pos_text_y = 0.816;
  float pos_y_diff = 0.056;
  float text_size = 19;
  int text_color = 1;
  if(ptLow==0 && ptHigh!=2.5) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 2.5 && ptHigh==5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 0 && ptHigh==2.5) drawText(Form("p_{T}^{#mu#mu} < %.1f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  drawText(Form("%.2f < y^{#mu#mu} < %.2f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);    // for pPb
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S)
  {
      drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.f GeV/c", muPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
  }


  TLegend* fitleg = new TLegend(0.76,0.4,0.91,0.7); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg->AddEntry(myPlot2->findObject("Sig1S"),"Signal","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"Background","l");
  fitleg->Draw("same");


  // PULL
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.30);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.67);
  pad1->SetLeftMargin(0.18);
  pad1->SetRightMargin(0.02);
  pad2->SetRightMargin(0.02);
  pad2->SetLeftMargin(0.18);
  pad2->SetTicks(1,1);
  pad2->cd();
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.43) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.18) ; //19
  pullFrame->GetYaxis()->SetLabelSize(0.113) ; // 113
  pullFrame->GetYaxis()->SetRangeUser(-3.8,3.8) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.05) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.20) ; //23
  pullFrame->GetXaxis()->SetTitleSize(0.25) ;  //28
  pullFrame->GetXaxis()->CenterTitle();

  pullFrame->GetYaxis()->SetTickSize(0.04);
  pullFrame->GetYaxis()->SetNdivisions(404);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  //calculate chi-squared
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
  cout << "chisq = " << chisq << endl;

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;
  chisqtest[imodel] = chisq/ndf;
  cout << "chisq/dof = " << chisqtest[imodel] << endl;

  //continue beautifying the plot and print out results
  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  pad1->Update();

              
  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  TString label;
  label="";
  if(collId == kPPDATA) CMS_lumi(pad1, 1 ,33);
  else if(collId == kAADATA && cLow < 60) CMS_lumi(pad1, 2 ,33);
  else if(collId == kPADATA) CMS_lumi(pad1, 3 ,33);
  else if(collId == kAADATA && cLow>=60) CMS_lumi(pad1, 21 ,33);


  pad1->Update();
  pad2->Update();

  cfit->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();

  
  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  outh->GetXaxis()->SetBinLabel(2,"Upsilon2S");
  outh->GetXaxis()->SetBinLabel(3,"Upsilon3S");
  
  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  float temp2 = ws->var("nSig2s")->getVal();  
  float temp2err = ws->var("nSig2s")->getError();  
  float temp3 = ws->var("nSig3s")->getVal();  
  float temp3err = ws->var("nSig3s")->getError();  
  
  outh->SetBinContent(1,  temp1 ) ;
  outh->SetBinError  (1,  temp1err ) ;
  outh->SetBinContent(2,  temp2 ) ;
  outh->SetBinError  (2,  temp2err ) ;
  outh->SetBinContent(3,  temp3 ) ;
  outh->SetBinError  (3,  temp3err ) ;

  cout << "1S signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;
  cout << "2S signal    =  " << outh->GetBinContent(2) << " +/- " << outh->GetBinError(2) << endl;
  cout << "3S signal    =  " << outh->GetBinContent(3) << " +/- " << outh->GetBinError(3) << endl;

	cout << "if ( binMatched( "<<muPtCut<<",  " << ptLow <<", "<< ptHigh<<", "<< yLow<<", "<< yHigh << ", " << cLow << ", " << cHigh << ") ) " ; 
  cout << "  { setSignalParMC( " ;
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("alpha1s_1")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " ;
  cout <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;

  yield1s[imodel] = outh->GetBinContent(1);
  yield2s[imodel] = outh->GetBinContent(2);
  yield3s[imodel] = outh->GetBinContent(3);

}//end of model fit loop

  //record % differences if the fits are good
  if (chisqtest[0]<10 && chisqtest[1]<10) {
    double perdif1s = 100*(yield1s[1]-yield1s[0])/yield1s[0];
    double perdif2s = 100*(yield2s[1]-yield2s[0])/yield2s[0];
    double perdif3s = 100*(yield3s[1]-yield3s[0])/yield3s[0];
    histo1s->Fill(perdif1s);
    histo2s->Fill(perdif2s);
    histo3s->Fill(perdif3s);
    ntuple->Fill(chisqtest[0],chisqtest[1],perdif1s,perdif2s,perdif3s);
    cout << "nominal chi^2 = " << chisqtest[0] << endl;
    cout << "alternate chi^2 = " << chisqtest[1] << endl;
    cout << "diff1s = " << perdif1s << endl;
    cout << "diff2s = " << perdif2s << endl;
    cout << "diff3s = " << perdif3s << endl;
    cout << "Trial " << itrial+1 << " of " << numtrials << " completed." << endl;

  }
  else {
    itrial--;
    cout << "Most recent trial rejected due to bad fit." << endl;
  }

  histo1s->SetStats(kTRUE);
  cyields->Update();
  cyields->cd(1);
  histo1s->SetStats(kTRUE);
  histo1s->Draw();
  cyields->cd(2);
  histo2s->Draw();
  cyields->cd(3);
  histo3s->Draw();

}//end of main loop

  //save histograms
  TFile outfile ("systematic_histos.root", "RECREATE");
  histo1s->Write();
  histo2s->Write();
  histo3s->Write();
  ntuple->Write();
  outfile.Close();

} 
 
