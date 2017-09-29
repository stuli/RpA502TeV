#include <iostream>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TFile.h"


using namespace RooFit;
void Simple_Fake_Data() {
  
  gStyle->SetEndErrorSize(0);

  float massLow = 8; 
  float massHigh = 14;
  int   nMassBin  = (massHigh-massLow)*10;

  //import model
  const char *inputFile="oldBkgModel.root";
  TFile *thefile = new TFile(inputFile);
  RooAddPdf* genModel = (RooAddPdf*)thefile->Get("model;1");
  RooWorkspace *wsgen = new RooWorkspace("workspace");
  wsgen->import(*genModel);

  //Generate fake data from the model
  //The real data set had 12954 events.
  RooDataSet* fakedata = genModel->generate(*(wsgen->var("mass")));
  fakedata->SetName("fakedata");
  wsgen->import(*fakedata);
  TFile fakefile ("FakeData.root", "RECREATE");
  fakedata->Write();
  fakefile.Close();

  //Plot it
  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  RooPlot* myPlot2 = wsgen->var("mass")->frame(nMassBin); // 60 bins
  wsgen->data("fakedata")->plotOn(myPlot2,Name("fakeDataHist"),MarkerSize(.8));
  myPlot2->SetFillStyle(4000);
  myPlot2->GetYaxis()->SetTitleOffset(1.43);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.058);
  //myPlot2->GetYaxis()->SetLabelSize(0.054);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->Draw();


} 
 
