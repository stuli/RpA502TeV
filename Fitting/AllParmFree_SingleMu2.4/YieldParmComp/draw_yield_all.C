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
#include "../../../cutsAndBin.h"

using namespace std;
using namespace RooFit;

void draw_yield_all()
{
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


  TFile *frap_1 = new TFile("Yield_Rap_PA_1s.root","read");
  TFile *frap_2 = new TFile("Yield_Rap_PA_2s.root","read");
  TFile *frap_3 = new TFile("Yield_Rap_PA_3s.root","read");
  TFile *fpt_1 = new TFile("Yield_Pt_PA_1s.root","read");
  TFile *fpt_2 = new TFile("Yield_Pt_PA_2s.root","read");
  TFile *fpt_3 = new TFile("Yield_Pt_PA_3s.root","read");

  TH1D* h_rap_1 = (TH1D*) frap_1 -> Get("h1_nSig1s");
  TH1D* h_rap_2 = (TH1D*) frap_2 -> Get("h1_nSig2s");
  TH1D* h_rap_3 = (TH1D*) frap_3 -> Get("h1_nSig3s");

  TH1D* h_pt_1 = (TH1D*) fpt_1 -> Get("h1_nSig1s");
  TH1D* h_pt_2 = (TH1D*) fpt_2 -> Get("h1_nSig2s");
  TH1D* h_pt_3 = (TH1D*) fpt_3 -> Get("h1_nSig3s");

  TH1D* h_rap_all = (TH1D*) h_rap_1->Clone("h_rap_all");
  TH1D* h_pt_all = (TH1D*) h_pt_1->Clone("h_pt_all");

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
  TCanvas* c_rap = new TCanvas("c_rap","c_rap",600,600);
  c_rap->cd();
  c_rap->SetLogy(1);
  h_rap_all->GetXaxis()->CenterTitle(1);
  h_rap_all->GetYaxis()->CenterTitle(1);
  h_rap_all->GetYaxis()->SetRangeUser(1,100000);
  h_rap_all->SetMarkerColor(kRed);
  h_rap_all->SetLineColor(kRed);
  h_rap_all->SetMarkerStyle(kFullCircle);
  h_rap_all->Draw("pe");
  h_rap_2->SetMarkerColor(kBlue);
  h_rap_2->SetLineColor(kBlue);
  h_rap_2->SetMarkerStyle(kFullCircle);
  h_rap_3->SetMarkerColor(kGreen);
  h_rap_3->SetLineColor(kGreen);
  h_rap_3->SetMarkerStyle(kFullCircle);
  h_rap_2->Draw("same pe");
  h_rap_3->Draw("same pe");
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,"pPb");
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.55,0.79,"#Upsilon(1S)");
  latex->SetTextColor(kBlue);
  latex->DrawLatex(0.55,0.73,"#Upsilon(2S)");
  latex->SetTextColor(kGreen);
  latex->DrawLatex(0.55,0.67,"#Upsilon(3S)");

  TCanvas* c_pt = new TCanvas("c_pt","c_pt",600,600);
  c_pt->cd();
  c_pt->SetLogy(1);
  h_pt_all->GetXaxis()->CenterTitle(1);
  h_pt_all->GetYaxis()->CenterTitle(1);
  h_pt_all->GetYaxis()->SetRangeUser(1,10000);
  h_pt_all->SetMarkerColor(kRed);
  h_pt_all->SetLineColor(kRed);
  h_pt_all->SetMarkerStyle(kFullCircle);
  h_pt_all->Draw("pe");
  h_pt_2->SetMarkerColor(kBlue);
  h_pt_2->SetLineColor(kBlue);
  h_pt_2->SetMarkerStyle(kFullCircle);
  h_pt_3->SetMarkerColor(kGreen);
  h_pt_3->SetLineColor(kGreen);
  h_pt_3->SetMarkerStyle(kFullCircle);
  h_pt_2->Draw("same pe");
  h_pt_3->Draw("same pe");
  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,"pPb");
  latex->SetTextColor(kRed);
  latex->DrawLatex(0.55,0.79,"#Upsilon(1S)");
  latex->SetTextColor(kBlue);
  latex->DrawLatex(0.55,0.73,"#Upsilon(2S)");
  latex->SetTextColor(kGreen);
  latex->DrawLatex(0.55,0.67,"#Upsilon(3S)");

}




