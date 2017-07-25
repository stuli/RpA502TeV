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

int draw_yield_pt_OldvsNew_ratio(int states =1, TString ColId = "PP")
{
  /////////////////////////////////////////////////////////
  //// set style
  /////////////////////////////////////////////////////////
  
  TH1::SetDefaultSumw2();

  TString szAA = "old";

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
  TFile *fileIn1[nBin];
  RooWorkspace* ws1[nBin];
  TFile *fileIn2[nBin];
  RooWorkspace* ws2[nBin];
  // parameters 
  double nSig1s_old[nBin];
  double nSig1sErr_old[nBin];
  double nSig2s_old[nBin];
  double nSig2sErr_old[nBin];
  double nSig3s_old[nBin];
  double nSig3sErr_old[nBin];
  double nBkg_old[nBin];
  double nBkgErr_old[nBin];

  double nSig1s_new[nBin];
  double nSig1sErr_new[nBin];
  double nSig2s_new[nBin];
  double nSig2sErr_new[nBin];
  double nSig3s_new[nBin];
  double nSig3sErr_new[nBin];
  double nBkg_new[nBin];
  double nBkgErr_new[nBin];

  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    //if (szAA == "old" ) { fileIn[ib]= new TFile(Form("/cms/home/goni/work/Analysis/Upsilon_Raa_8_0_24/src/usercode/upsilonRAA5TeV/fitResults/NominalFits/fitresults_upsilon_DoubleCB_AA_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",binArr[ib],binArr[ib+1])); }
    //else if (szAA == "AA" ) { fileIn[ib]= new TFile(Form("/cms/home/goni/work/Analysis/Upsilon_Raa_8_0_24/src/usercode/upsilonRAA5TeV/fitResults/NominalFits/PAS_fitresults_upsilon_DoubleCB_%s_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",szAA.Data(),binArr[ib],binArr[ib+1])); }

    if(ColId=="AA"){
    fileIn1[ib]= new TFile(Form("NomPlot/fitresults_upsilon_DoubleCB_AA_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",binArr[ib],binArr[ib+1]));
    fileIn2[ib]= new TFile(Form("fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_AA_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI.root",binArr[ib],binArr[ib+1])); }
    if(ColId=="PP"){
    fileIn1[ib]= new TFile(Form("NomPlot/fitresults_upsilon_DoubleCB_PP_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0.root",binArr[ib],binArr[ib+1]));
    fileIn2[ib]= new TFile(Form("fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_DoubleCB_PP_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0.root",binArr[ib],binArr[ib+1])); }
    //else { cout << " Error ::: Select among old and AA" << endl; return 0; }

    fileIn1[ib]->cd();
    ws1[ib]= (RooWorkspace*)fileIn1[ib]->Get("workspace");
    nSig1s_old[ib]=ws1[ib]->var("nSig1s")->getVal();
    nSig1sErr_old[ib]=ws1[ib]->var("nSig1s")->getError();
    nSig2s_old[ib]=ws1[ib]->var("nSig2s")->getVal();
    nSig2sErr_old[ib]=ws1[ib]->var("nSig2s")->getError();
    nSig3s_old[ib]=ws1[ib]->var("nSig3s")->getVal();
    nSig3sErr_old[ib]=ws1[ib]->var("nSig3s")->getError();
    nBkg_old[ib]=ws1[ib]->var("nBkg")->getVal();
    nBkgErr_old[ib]=ws1[ib]->var("nBkg")->getError();

    fileIn2[ib]->cd();
    ws2[ib]= (RooWorkspace*)fileIn2[ib]->Get("workspace");
    nSig1s_new[ib]=ws2[ib]->var("nSig1s")->getVal();
    nSig1sErr_new[ib]=ws2[ib]->var("nSig1s")->getError();
    nSig2s_new[ib]=ws2[ib]->var("nSig2s")->getVal();
    nSig2sErr_new[ib]=ws2[ib]->var("nSig2s")->getError();
    nSig3s_new[ib]=ws2[ib]->var("nSig3s")->getVal();
    nSig3sErr_new[ib]=ws2[ib]->var("nSig3s")->getError();
    nBkg_new[ib]=ws2[ib]->var("nBkg")->getVal();
    nBkgErr_new[ib]=ws2[ib]->var("nBkg")->getError();
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

  TH1D* h2_nSig1s = new TH1D("h2_nSig1s","h2_nSig1s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h2_nSig2s = new TH1D("h2_nSig2s","h2_nSig2s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h2_nSig3s = new TH1D("h2_nSig3s","h2_nSig3s;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 
  TH1D* h2_nBkg = new TH1D("h2_nBkg","h2_nBkg;p_{T} (GeV/c);events / ( GeV/c )",nBin,binArr); 

  TFile *fout=new TFile(Form("yield/Yields_%ds.root", states),"recreate");
  
  for (int ib =0; ib < nBin; ib ++ ) {
    h1_nSig1s->SetBinContent(ib+1,nSig1s_old[ib]);   
    h1_nSig1s->SetBinError(ib+1,nSig1sErr_old[ib]);   
    h1_nSig2s->SetBinContent(ib+1,nSig2s_old[ib]);   
    h1_nSig2s->SetBinError(ib+1,nSig2sErr_old[ib]);   
    h1_nSig3s->SetBinContent(ib+1,nSig3s_old[ib]);   
    h1_nSig3s->SetBinError(ib+1,nSig3sErr_old[ib]);   
    h1_nBkg->SetBinContent(ib+1,nBkg_old[ib]);   
    h1_nBkg->SetBinError(ib+1,nBkgErr_old[ib]);   

    h2_nSig1s->SetBinContent(ib+1,nSig1s_new[ib]);   
    h2_nSig1s->SetBinError(ib+1,nSig1sErr_new[ib]);   
    h2_nSig2s->SetBinContent(ib+1,nSig2s_new[ib]);   
    h2_nSig2s->SetBinError(ib+1,nSig2sErr_new[ib]);   
    h2_nSig3s->SetBinContent(ib+1,nSig3s_new[ib]);   
    h2_nSig3s->SetBinError(ib+1,nSig3sErr_new[ib]);   
    h2_nBkg->SetBinContent(ib+1,nBkg_new[ib]);   
    h2_nBkg->SetBinError(ib+1,nBkgErr_new[ib]);   

    fout->cd();
    h1_nSig1s->Write();
    h1_nSig2s->Write();
    h1_nSig3s->Write();
    //h1_nbkg->Write();
    h2_nSig1s->Write();
    h2_nSig2s->Write();
    h2_nSig3s->Write();
    //h2_nbkg->Write();
  }

  cout<< "nSig1S(1): " << h1_nSig1s->GetBinContent(1) <<endl;
  cout<< "nSig1S(2): " << h1_nSig1s->GetBinContent(2) <<endl;
  cout<< "nSig1S(3): " << h1_nSig1s->GetBinContent(3) <<endl;
  cout<< "nSig1S(4): " << h1_nSig1s->GetBinContent(4) <<endl;
  cout<< "nSig1S(5): " << h1_nSig1s->GetBinContent(5) <<endl;
  cout<< "nSig1S(6): " << h1_nSig1s->GetBinContent(6) <<endl;
  
  cout<< "nSig2S(1): " << h1_nSig2s->GetBinContent(1) <<endl;
  cout<< "nSig2S(2): " << h1_nSig2s->GetBinContent(2) <<endl;
  cout<< "nSig2S(3): " << h1_nSig2s->GetBinContent(3) <<endl;
  
  cout<< "nSig3S(1): " << h1_nSig3s->GetBinContent(1) <<endl;
  cout<< "nSig3S(2): " << h1_nSig3s->GetBinContent(2) <<endl;

  //// normalization
  //h1_nSig1s->Scale(1,"width");
  //h1_nSig2s->Scale(1,"width");
  //h1_nSig3s->Scale(1,"width");
  //h1_nBkg->Scale(1,"width");
  //h2_nSig1s->Scale(1,"width");
  //h2_nSig2s->Scale(1,"width");
  //h2_nSig3s->Scale(1,"width");
  //h2_nBkg->Scale(1,"width");

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
 
  TCanvas* c_nSig1s = new TCanvas("c_nSig1s","c_nSig1s",700,700);
  c_nSig1s->cd();
  TPad *p1 = new TPad("p1","p1", 0, 0.16, 0.98, 1.0);
  p1->SetTicks(1,1);
  p1->Draw(); p1->cd();
  p1->SetLogy(1);
  
  
//  gPad->SetLogy(1);  // KTO for yield
  h1_nSig1s->GetXaxis()->SetLabelSize(0);
  h1_nSig1s->GetYaxis()->CenterTitle(1);
  h1_nSig1s->GetYaxis()->SetRangeUser(1,1000000);
  h1_nSig1s->SetMarkerColor(kBlue+2);
  h1_nSig1s->SetLineColor(kBlue+2);
  h1_nSig1s->SetMarkerStyle(kFullCircle);
  h1_nSig1s->Draw("pe");

  h2_nSig1s->GetXaxis()->CenterTitle(1);
  h2_nSig1s->GetYaxis()->CenterTitle(1);
  h2_nSig1s->GetYaxis()->SetRangeUser(1,1000000);
  h2_nSig1s->SetMarkerColor(kRed+2);
  h2_nSig1s->SetLineColor(kRed+2);
  h2_nSig1s->SetMarkerStyle(kOpenCircle);
  h2_nSig1s->Draw("pe same");

  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",ColId.Data()));
  latex->SetTextColor(kRed+2);
  latex->DrawLatex(0.55,0.79,Form("(New) #Upsilon(%dS) yields",states));
  latex->SetTextColor(kBlue+2);
  latex->DrawLatex(0.55,0.73,Form("(Old) #Upsilon(%dS) yields",states));
  
  TPad *p2 = new TPad("p2","p2", 0, 0.03, 0.98, 0.24);
  p2->SetTopMargin(0); // Upper and lower plot are joined
  p2->SetBottomMargin(0.43);
  p1->SetLeftMargin(0.15);
  p2->SetLeftMargin(0.15);
  p2->SetTicks(1,1);
  p2->cd();

  TH1D* h_nSig1s_Ratio = (TH1D*) h2_nSig1s->Clone("h_nSig1s_Ratio");
  h_nSig1s_Ratio->Divide(h1_nSig1s);
  h_nSig1s_Ratio->SetAxisRange(0.8,1.2,"Y");
  h_nSig1s_Ratio->SetTitleSize(0);
  h_nSig1s_Ratio->GetYaxis()->SetTitleOffset(0.36) ;
  h_nSig1s_Ratio->GetYaxis()->SetTitle("Ratio") ;
  h_nSig1s_Ratio->GetYaxis()->SetTitleSize(0.14) ;
  h_nSig1s_Ratio->GetYaxis()->SetLabelSize(0.13) ;
  h_nSig1s_Ratio->SetMarkerColor(kBlack);
  h_nSig1s_Ratio->SetMarkerStyle(20);
  h_nSig1s_Ratio->SetLineColor(kBlack);

  h_nSig1s_Ratio->GetYaxis()->SetRangeUser(0.8,1.2) ;
//  h_nSig1s_Ratio->GetYaxis()->SetLimits(-6,6) ;
  h_nSig1s_Ratio->GetYaxis()->CenterTitle();

  h_nSig1s_Ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_nSig1s_Ratio->GetXaxis()->SetTitleOffset(1.04) ;
  h_nSig1s_Ratio->GetXaxis()->SetLabelOffset(0.04) ;
  h_nSig1s_Ratio->GetXaxis()->SetLabelSize(0.180) ;
  h_nSig1s_Ratio->GetXaxis()->SetTitleSize(0.196) ;
  
  h_nSig1s_Ratio->GetXaxis()->CenterTitle();
 // h_nSig1s_Ratio->GetXaxis()->SetTitleFont(43);
 // h_nSig1s_Ratio->GetYaxis()->SetTitleFont(43);
  
  h_nSig1s_Ratio->GetYaxis()->SetTickSize(0.02);
  h_nSig1s_Ratio->GetYaxis()->SetNdivisions(505);
  h_nSig1s_Ratio->GetXaxis()->SetTickSize(0.03);
  h_nSig1s_Ratio->Draw() ;
  
  TLine *l1 = new TLine(0,1,30,1);
  l1->SetLineStyle(9);
  l1->Draw("same");

  p1->Update();
  p2->Update();
  c_nSig1s->cd();
  p1->Draw();
  p2->Draw();
  p1->Update();
  p2->Update();



  c_nSig1s->SaveAs(Form("yield/pt_nSig1s_%s_%dS.pdf",ColId.Data(),states));
  
  TCanvas* c_nSig2s = new TCanvas("c_nSig2s","c_nSig2s",600,600);
  c_nSig2s->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nSig2s->GetXaxis()->CenterTitle(1);
  h1_nSig2s->GetYaxis()->CenterTitle(1);
  h1_nSig2s->GetYaxis()->SetRangeUser(1,8000);
  h1_nSig2s->SetMarkerColor(kBlue+2);
  h1_nSig2s->SetLineColor(kBlue+2);
  h1_nSig2s->SetMarkerStyle(kFullCircle);
  h1_nSig2s->Draw("pe");

  h2_nSig2s->GetXaxis()->CenterTitle(1);
  h2_nSig2s->GetYaxis()->CenterTitle(1);
  h2_nSig2s->GetYaxis()->SetRangeUser(1,8000);
  h2_nSig2s->SetMarkerColor(kRed+2);
  h2_nSig2s->SetLineColor(kRed+2);
  h2_nSig2s->SetMarkerStyle(kOpenCircle);
  h2_nSig2s->Draw("pe same");

  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",ColId.Data()));
  latex->SetTextColor(kRed+2);
  latex->DrawLatex(0.55,0.79,Form("(New) #Upsilon(%dS) yields",states));
  latex->SetTextColor(kBlue+2);
  latex->DrawLatex(0.55,0.73,Form("(Old) #Upsilon(%dS) yields",states));
  c_nSig2s->SaveAs(Form("yield/pt_nSig2s_%s_%dS.pdf",ColId.Data(),states));
  
  TCanvas* c_nSig3s = new TCanvas("c_nSig3s","c_nSig3s",600,600);
  c_nSig3s->cd();
  gPad->SetLogy(1);  // KTO for yield
  h1_nSig3s->GetXaxis()->CenterTitle(1);
  h1_nSig3s->GetYaxis()->CenterTitle(1);
  h1_nSig3s->GetYaxis()->SetRangeUser(1,3000);
  h1_nSig3s->SetMarkerColor(kBlue+2);
  h1_nSig3s->SetLineColor(kBlue+2);
  h1_nSig3s->SetMarkerStyle(kFullCircle);
  h1_nSig3s->Draw("pe");

  h2_nSig3s->SetMarkerColor(kRed+2);
  h2_nSig3s->SetLineColor(kRed+2);
  h2_nSig3s->SetMarkerStyle(kOpenCircle);
  h2_nSig3s->Draw("pe same");

  latex->SetTextColor(kBlack);
  latex->DrawLatex(0.55,0.85,Form("%s",ColId.Data()));
  latex->SetTextColor(kRed+2);
  latex->DrawLatex(0.55,0.79,Form("(New) #Upsilon(%dS) yields",states));
  latex->SetTextColor(kBlue+2);
  latex->DrawLatex(0.55,0.73,Form("(Old) #Upsilon(%dS) yields",states));
  c_nSig3s->SaveAs(Form("yield/pt_nSig3s_%s_%dS.pdf",ColId.Data(),states));
 
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

