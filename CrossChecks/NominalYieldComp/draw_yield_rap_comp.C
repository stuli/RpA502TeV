#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>

#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>
#include "../../SONGKYO.h"
#include "../../commonUtility.h"

using namespace std;
using namespace RooFit;

int draw_yield_rap_comp(TString szAA = "PA", int states =2)
{
  /////////////////////////////////////////////////////////
  //// set style
  /////////////////////////////////////////////////////////
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");
  gStyle->SetTitleOffset(1.6,"y"); // KYO for yield

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12) ; 
  gStyle->SetPadLeftMargin(0.16) ; // KYO for yield
  
  /////////////////////////////////////////////////////////
  //// binning setting
  /////////////////////////////////////////////////////////
  double tmpArr1s[9] = {-1.93,-1.2,-0.8,-0.4,0, 0.4, 0.8, 1.2, 1.93};
  double tmpArr2s[5] = {-1.93,-0.8,0.0, 0.8, 1.93};
  double tmpArr3s[3] = {-1.93, 0.0, 1.93};
  
  int tmpBin;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin = 8;
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin = 4;
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin = 2;
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl; return 0;
  }
 
  const int nStates = 3; 
  const int nFit = 2;
  const int nBin = tmpBin; // number of bin 
  const int nArrNum = nBin+1; // number of array
  double binArr[nArrNum]; // array
  cout << "nBin = " << nBin << endl;
  for (int ib =0; ib < nArrNum; ib ++ ) {
    if (states ==1) { binArr[ib] = tmpArr1s[ib]; }
    else if (states ==2) { binArr[ib] = tmpArr2s[ib]; }
    else if (states ==3) { binArr[ib] = tmpArr3s[ib]; }
    cout << ib <<"th bin = " << binArr[ib] << endl;
  }

  /////////////////////////////////////////////////////////
  //// Open RooDataFile
  /////////////////////////////////////////////////////////
 
  //file and ws
  TFile *fileIn[nFit][nBin];
  RooWorkspace* ws[nFit][nBin];
  // parameters 
  double nSig1s[nFit][nBin];
  double nSig1sErr[nFit][nBin];
  double nSig2s[nFit][nBin];
  double nSig2sErr[nFit][nBin];
  double nSig3s[nFit][nBin];
  double nSig3sErr[nFit][nBin];
  double nBkg[nFit][nBin];
  double nBkgErr[nFit][nBin];
 
  TString fileLoc[nFit] = {" ../../NominalFitResult/jaebeomFit/", "../../NominalFitResult/jaredFit/"}; 
  TString fitName[nFit] = {"JaeBeom", "Jared"}; 
  
  Int_t fitColorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };

  for (int ib =0; ib < nBin; ib ++ ) {
    for(int ifit = 0; ifit < nFit; ifit++){
      //// read files
      if (szAA == "PP" ) { 
        if(ifit==0 && ib < nBin/2) fileIn[ifit][ib]= new TFile(Form("%sAllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[nBin-1-ib],binArr[nBin-ib]));
        else if(ifit==0 && ib>=nBin/2) fileIn[ifit][ib]= new TFile(Form("%sAllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
        else {fileIn[ifit][ib]= new TFile(Form("%sAllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));}
      }
      else if (szAA == "PA" ) { fileIn[ifit][ib]= new TFile(Form("%sAllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s_DATA_pt0.0-30.0_y%.2f-%.2f_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1])); }
      else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
      cout << ib << "th file = " << fileIn[ifit][ib]->GetName() << endl;
      if (fileIn[ifit][ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
      fileIn[ifit][ib]->cd();
      ws[ifit][ib]= (RooWorkspace*)fileIn[ifit][ib]->Get("workspace");
      //ws[ifit][ib]->Print();

      //// get parameters
      if(ifit==0 && szAA == "PP"){
        nSig1s[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getVal()/2;
        nSig1sErr[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getError()/2;
        nSig2s[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getVal()/2;
        nSig2sErr[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getError()/2;
        nSig3s[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getVal()/2;
        nSig3sErr[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getError()/2;
        nBkg[ifit][ib]=ws[ifit][ib]->var("nBkg")->getVal()/2;
        nBkgErr[ifit][ib]=ws[ifit][ib]->var("nBkg")->getError()/2;
      }
      else {
        nSig1s[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getVal();
        nSig1sErr[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getError();
        nSig2s[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getVal();
        nSig2sErr[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getError();
        nSig3s[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getVal();
        nSig3sErr[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getError();
        nBkg[ifit][ib]=ws[ifit][ib]->var("nBkg")->getVal();
        nBkgErr[ifit][ib]=ws[ifit][ib]->var("nBkg")->getError();
      }
      //cout << ib << "th nSig1s = " << nSig1s[ifit][ib] << endl;
      //cout << ib << "th nSig2s = " << nSig2s[ifit][ib] << endl;
      //cout << ib << "th nSig3s = " << nSig3s[ifit][ib] << endl;
      //cout << ib << "th nBkg = " << nBkg[ifit][ib] << endl;
    }
  }

  //// histogram
  TH1D* h1_nSig[nFit]; 
  TH1D* h1_nBkg[nFit]; 
  
  for(int ifit=0; ifit<nFit; ifit++){
    h1_nSig[ifit] = new TH1D(Form("h1_nSig%ds_%d",states,ifit+1),Form("h1_nSig%ds;y_{CM};dN_{(#Upsilon%dS)}/dy_{CM}",states,states),nBin,binArr); 
    h1_nBkg[ifit] = new TH1D(Form("h1_nBkg_%d",ifit+1),"h1_nBkg;y_{CM};dN_{Bkg}/dy_{CM}",nBin,binArr); 
    for (int ib =0; ib < nBin; ib ++ ) {
      if(states ==1) { h1_nSig[ifit]->SetBinContent(ib+1,nSig1s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig1sErr[ifit][ib]);}   
      else if(states ==2) { h1_nSig[ifit]->SetBinContent(ib+1,nSig2s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig2sErr[ifit][ib]);}   
      else if(states ==3) { h1_nSig[ifit]->SetBinContent(ib+1,nSig3s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig3sErr[ifit][ib]);}   
      h1_nBkg[ifit]->SetBinContent(ib+1,nBkg[ifit][ib]);   
      h1_nBkg[ifit]->SetBinError(ib+1,nBkgErr[ifit][ib]);   
    }
  }

  //// normalization
  for(int ifit=0; ifit<nFit; ifit++){
    TH1ScaleByWidth(h1_nSig[ifit]);
    TH1ScaleByWidth(h1_nBkg[ifit]);
    SetHistStyle(h1_nSig[ifit],ifit, ifit);
    SetHistStyle(h1_nBkg[ifit],ifit, ifit);
  }

  int binmax = h1_nSig[0]->GetMaximumBin();
  double valmax = h1_nSig[0]->GetBinContent(binmax);
  int binmax_bkg = h1_nBkg[0]->GetMaximumBin();
  double valmax_bkg = h1_nBkg[0]->GetBinContent(binmax_bkg);
  cout << " binmax : " << binmax << endl;
  cout << " valmax : " << valmax << endl;
  h1_nSig[0]->GetYaxis()->SetRangeUser(10,valmax*100);
  h1_nBkg[0]->GetYaxis()->SetRangeUser(10,valmax_bkg*100);
  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.035);

  double pos_y_diff = 0.05; 

  TLegend* fitleg = new TLegend(0.56,0.65,0.71,0.75); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);

  TCanvas* c_nSigs = new TCanvas("c_nSigs","c_nSigs",600,600);
  c_nSigs->cd();
  gPad->SetLogy();  // KTO for yield
  for(int ifit=0; ifit<nFit; ifit++){ 
    h1_nSig[0]->Draw("pe"); if(ifit>0) h1_nSig[ifit]->Draw("pe same");
  }
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.86,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  for(int ifit=0;ifit<nFit; ifit++){
    fitleg->AddEntry(h1_nSig[ifit],Form("%s fit",fitName[ifit].Data()),"pe");
  }
  fitleg->Draw("same");
  c_nSigs->SaveAs(Form("yield/rap_nSig%ds_%s.pdf",states,szAA.Data()));
  
    //latex->SetTextColor(fitColorArr[ifit]);
    //latex->DrawLatex(0.55,0.81-pos_y_diff*(ifit+1),Form("%s fit",fitName[ifit].Data()));

  TCanvas* c_nBkg = new TCanvas("c_nBkg","c_nBkg",600,600);
  c_nBkg->cd();
  gPad->SetLogy();  // KTO for yield
  for(int ifit=0; ifit<nFit; ifit++){ 
    h1_nBkg[0]->Draw("pe"); if(ifit>0) h1_nBkg[ifit]->Draw("pe same");   
  }
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.86,Form("%s Background",szAA.Data()));
  for(int ifit=0;ifit<nFit; ifit++){
  }
  fitleg->Draw("same");
  c_nBkg->SaveAs(Form("yield/rap_nBkg%ds_%s.pdf",states,szAA.Data()));
  
  return 0;
}

