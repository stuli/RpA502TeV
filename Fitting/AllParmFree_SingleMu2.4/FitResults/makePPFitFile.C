#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMath.h"
#include "../../../cutsAndBin.h"

using namespace std;

void makePPFitFile(int state=2)
{

  int nPtBins=0;
  int nYBins=0;
  double* ptBin;
  double* yBin;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;


  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; 
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S; 
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S; 
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  int BinOne, BinTwo;
  double SetBinOne,SetBinTwo,SetBinThree;
  if(state ==1 ){BinOne = 7; BinTwo = 8; SetBinOne = 1.2; SetBinTwo = 1.6; SetBinThree = 2.0;}
  else if(state ==2 ){BinOne = 3; BinTwo = 4; SetBinOne = 0.8; SetBinTwo = 1.6; SetBinThree = 2.4;}
  else if(state ==3 ){BinOne = 1; BinTwo = 2; SetBinOne = 0.0; SetBinTwo = 2.4;}
  int numBinOne[2] = {BinOne, BinTwo};

  TFile *f1_1 = new TFile(Form("PAS_fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",SetBinOne,SetBinTwo),"read");
  TFile *f1_2 = new TFile(Form("PAS_fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0.root",SetBinOne,SetBinTwo),"read");

  TFile *f2 = new TFile(Form("New_AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PP_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",yBin[BinOne],yBin[BinTwo]),"read");

  TFile *f_new = new TFile(Form("AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PP_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",yBin[BinOne],yBin[BinTwo]),"recreate");
  TH1D *h1 = (TH1D*) f1_1 -> Get("fitResults");
  TH1D *h2 = (TH1D*) f1_2 -> Get("fitResults");
  TH1D *h3 = (TH1D*) f2 -> Get("fitResults");

  RooWorkspace* ws = (RooWorkspace*) f2 -> Get("workspace");

  double newval,newerr,ratio,newerr_ratio;
  if(state!=3){
    newval = h1->GetBinContent(state)+h2->GetBinContent(state);
    newerr = h1->GetBinError(state)*h1->GetBinError(state)+h2->GetBinError(state)*h2->GetBinError(state);
    ratio = (double)h3->GetBinContent(state)/newval;
    newerr_ratio = TMath::Sqrt(newerr*ratio);

    cout << Form("Yield %.1f-%.1f : ",SetBinOne,SetBinTwo) <<  h1->GetBinContent(state) << endl;
    cout << Form("Yield %.1f-%.1f : ",SetBinTwo,SetBinThree) <<  h2->GetBinContent(state) << endl;
    cout << Form("Yield merge %.1f-%.1f : ",SetBinOne,SetBinThree) << newval << endl;
    cout << Form("Yield %.2f-%.2f : ",yBin[BinOne],yBin[BinTwo]) << h3->GetBinContent(state) << endl;
    cout << Form("Err %.1f-%.1f : ",SetBinOne,SetBinTwo) << h1->GetBinError(state) << endl;
    cout << Form("Err %.1f-%.1f : ",SetBinTwo,SetBinThree) << h2->GetBinError(state) << endl;
    cout << "Err merge : " << TMath::Sqrt(newerr) << endl;
    cout << "newerr_ratio : " << newerr_ratio << endl;
  }
  else if(state==3){
    newval = h1->GetBinContent(state);
    newerr = h1->GetBinError(state)*h1->GetBinError(state);
    ratio = (double)h3->GetBinContent(state)/newval;
    newerr_ratio = TMath::Sqrt(newerr*ratio);

    cout << "Yield 0.0-2.4 : " <<  h1->GetBinContent(state) << endl;
    cout << "Yield 0.00-1.93 : " << h3->GetBinContent(state) << endl;
    cout << "Err 0.0-2.4 : " << h1->GetBinError(state) << endl;
    cout << "newerr_ratio : " << newerr_ratio << endl;
  }


  TH1D *h_ = (TH1D*) h1->Clone("fitResults");
  h_->SetBinContent(state,h3->GetBinContent(state));
  h_->SetBinError(state,newerr_ratio);
  f_new->cd();
  h_->Write();
  ws->Write();
}
