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

using namespace std;
using namespace RooFit;

int draw_yield_pt(TString szAA = "PP", int states =1)
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
  
  double tmpArr1s[7] = {0.0, 2.0, 4.0, 6.0, 9.0, 12.0, 30.0};
  double tmpArr2s[4] = {0.0, 4.0, 9.0, 30.0};
  double tmpArr3s[3] = {0.0, 6.0, 30.0};

  int tmpBin;
  if ( states ==1) {
    cout << " ***** 1S *****" << endl; tmpBin = 6;
  }else if (states ==2){
    cout << " ***** 2S *****" << endl; tmpBin = 3;
  }else if (states ==3){
    cout << " ***** 3S *****" << endl; tmpBin = 2;
  }else {
    cout << " Error ::: Select among 1S, 2S, and 3S" << endl; return 0;
  }
  
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
  TFile *fileIn[nBin];
  RooWorkspace* ws[nBin];
  // parameters 
  double nSig1s[nBin];
  double nSig1sErr[nBin];
  double nSig2s[nBin];
  double nSig2sErr[nBin];
  double nSig3s[nBin];
  double nSig3sErr[nBin];
  double nBkg[nBin];
  double nBkgErr[nBin];
  
  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    if (szAA == "PP" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_%s_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else if (szAA == "AA" ) { fileIn[ib]= new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_%s_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",szAA.Data(),binArr[ib],binArr[ib+1])); }
    else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
    //cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    //if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    //ws[ib]->Print();
    //// get parameters
    //cout << ws[ib]->var("nSig1s")->getVal() << endl;
    nSig1s[ib]=ws[ib]->var("nSig1s")->getVal();
    nSig1sErr[ib]=ws[ib]->var("nSig1s")->getError();
    nSig2s[ib]=ws[ib]->var("nSig2s")->getVal();
    nSig2sErr[ib]=ws[ib]->var("nSig2s")->getError();
    nSig3s[ib]=ws[ib]->var("nSig3s")->getVal();
    nSig3sErr[ib]=ws[ib]->var("nSig3s")->getError();
    nBkg[ib]=ws[ib]->var("nBkg")->getVal();
    nBkgErr[ib]=ws[ib]->var("nBkg")->getError();
    //cout << ib << "th nSig1s = " << nSig1s[ib] << endl;
    //cout << ib << "th nSig2s = " << nSig2s[ib] << endl;
    //cout << ib << "th nSig3s = " << nSig3s[ib] << endl;
    //cout << ib << "th nBkg = " << nBkg[ib] << endl;
  }
 
  //// histogram
  TH1D* h1_nSig1s = new TH1D("h1_nSig1s","h1_nSig1s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h1_nSig2s = new TH1D("h1_nSig2s","h1_nSig2s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h1_nSig3s = new TH1D("h1_nSig3s","h1_nSig3s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h1_nBkg = new TH1D("h1_nBkg","h1_nBkg;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  
  for (int ib =0; ib < nBin; ib ++ ) {
    h1_nSig1s->SetBinContent(ib+1,nSig1s[ib]);   
    h1_nSig1s->SetBinError(ib+1,nSig1sErr[ib]);   
    h1_nSig2s->SetBinContent(ib+1,nSig2s[ib]);   
    h1_nSig2s->SetBinError(ib+1,nSig2sErr[ib]);   
    h1_nSig3s->SetBinContent(ib+1,nSig3s[ib]);   
    h1_nSig3s->SetBinError(ib+1,nSig3sErr[ib]);   
    h1_nBkg->SetBinContent(ib+1,nBkg[ib]);   
    h1_nBkg->SetBinError(ib+1,nBkgErr[ib]);   
  }

  //// normalization
  h1_nSig1s->Scale(1,"width");
  h1_nSig2s->Scale(1,"width");
  h1_nSig3s->Scale(1,"width");
  h1_nBkg->Scale(1,"width");

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
 
  TCanvas* c_nSig1s = new TCanvas("c_nSig1s","c_nSig1s",600,600);
  c_nSig1s->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nSig1s->GetXaxis()->CenterTitle(1);
  h1_nSig1s->GetYaxis()->CenterTitle(1);
  h1_nSig1s->GetYaxis()->SetRangeUser(1,100000);
  h1_nSig1s->SetMarkerColor(kRed);
  h1_nSig1s->SetLineColor(kRed);
  h1_nSig1s->SetMarkerStyle(kFullCircle);
  h1_nSig1s->Draw("pe");
  h1_nBkg->SetMarkerColor(kBlue);
  h1_nBkg->SetLineColor(kBlue);
  h1_nBkg->SetMarkerStyle(kOpenCircle);
  h1_nBkg->Draw("pe same");
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",szAA.Data()));
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.55,0.79,Form("#Upsilon(%dS) yields",states));
  latex->SetTextColor(kBlue);
  latex->DrawLatex(0.55,0.73,"Backgrounds");
  
  c_nSig1s->SaveAs(Form("yield/pt_nSig1s_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_nSig2s = new TCanvas("c_nSig2s","c_nSig2s",600,600);
  c_nSig2s->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nSig2s->GetXaxis()->CenterTitle(1);
  h1_nSig2s->GetYaxis()->CenterTitle(1);
  h1_nSig2s->GetYaxis()->SetRangeUser(1,100000);
  h1_nSig2s->SetMarkerColor(kRed);
  h1_nSig2s->SetLineColor(kRed);
  h1_nSig2s->SetMarkerStyle(kFullCircle);
  h1_nSig2s->Draw("pe");
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",szAA.Data()));
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.55,0.79,Form("#Upsilon(%dS) yields",states));
  c_nSig2s->SaveAs(Form("yield/pt_nSig2s_%s_%dS.pdf",szAA.Data(),states));
  
  TCanvas* c_nSig3s = new TCanvas("c_nSig3s","c_nSig3s",600,600);
  c_nSig3s->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nSig3s->GetXaxis()->CenterTitle(1);
  h1_nSig3s->GetYaxis()->CenterTitle(1);
  h1_nSig3s->GetYaxis()->SetRangeUser(1,100000);
  h1_nSig3s->SetMarkerColor(kRed);
  h1_nSig3s->SetLineColor(kRed);
  h1_nSig3s->SetMarkerStyle(kFullCircle);
  h1_nSig3s->Draw("pe");
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",szAA.Data()));
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.55,0.79,Form("#Upsilon(%dS) yields",states));
  c_nSig3s->SaveAs(Form("yield/pt_nSig3s_%s_%dS.pdf",szAA.Data(),states));
 
 /* 
  TCanvas* c_nBkg = new TCanvas("c_nBkg","c_nBkg",600,600);
  c_nBkg->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nBkg->GetXaxis()->CenterTitle(1);
  h1_nBkg->GetYaxis()->CenterTitle(1);
  h1_nBkg->GetYaxis()->SetRangeUser(1,100000);
  h1_nBkg->SetMarkerColor(kRed);
  h1_nBkg->SetLineColor(kRed);
  h1_nBkg->SetMarkerStyle(kFullCircle);
  h1_nBkg->Draw("pe");
  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_nBkg->SaveAs(Form("yield/pt_nBkg_%s_%dS.pdf",szAA.Data(),states));
*/

  return 0;

}

