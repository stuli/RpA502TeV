#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <cstring>
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
#include "../SONGKYO.h"
#include "../commonUtility.h"
#include "../cutsAndBin.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_internal.C"

using namespace std;
using namespace RooFit;

int draw_yield_pt_comp_constrain(TString szAA = "PA", int states =1) 
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod;
  if(szAA=="PP") iPeriod=1; 
  else if(szAA=="PA") iPeriod=3; 
  int iPos = 33;
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

  /*gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
*/
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
 
  const int nStates = 3; 
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

  ///////////////////////////////////////////////////////////// Open RooDataFile
  /////////////////////////////////////////////////////////
 
  int nFit = 2;

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
 
  char *Fit_loc[2] = {"../NominalFitResult/jaredFit/NominalFits/Pbp", "../NominalFitResult/jaredFit/NominalFits/pPb"};
  char *Name_Fit[2] = {"Run1", "Run2"};
  TString fileLoc[nFit];
  TString fitName[nFit]; 
  for(int i=0; i<nFit; i++)
  {
      fileLoc[i] = Fit_loc[i];
      fitName[i] = Name_Fit[i];
  }

  Int_t fitColorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };

  for (int ib =0; ib < nBin; ib ++ ) {
    for(int ifit = 0; ifit < nFit; ifit++){
      //// read files
      if (szAA == "PP" ) { 
        fileIn[ifit][ib]= new TFile(Form("%s/nomfitresults_upsilon_%s_DATA_pt%.1f-%.1f_y0.00-1.93_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
      }
      else if (szAA == "PA" ) { 
        fileIn[ifit][ib]= new TFile(Form("%s/nomfitresults_upsilon_%s_DATA_pt%.1f-%.1f_y-1.93-1.93_muPt4.0.root",fileLoc[ifit].Data(),szAA.Data(),binArr[ib],binArr[ib+1]));
      }
      else { cout << " Error ::: Select among PP and AA" << endl; return 0; }
      cout << ib << "th file = " << fileIn[ifit][ib]->GetName() << endl;
      if (fileIn[ifit][ib]->IsZombie()) { cout << "CANNOT open data root file\n"; return 1; }
      fileIn[ifit][ib]->cd();
      ws[ifit][ib]= (RooWorkspace*)fileIn[ifit][ib]->Get("workspace");
      //ws[ifit][ib]->Print();

      //// get parameters
      nSig1s[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getVal();
      nSig1sErr[ifit][ib]=ws[ifit][ib]->var("nSig1s")->getError();
      nSig2s[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getVal();
      nSig2sErr[ifit][ib]=ws[ifit][ib]->var("nSig2s")->getError();
      nSig3s[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getVal();
      nSig3sErr[ifit][ib]=ws[ifit][ib]->var("nSig3s")->getError();
      nBkg[ifit][ib]=ws[ifit][ib]->var("nBkg")->getVal();
      nBkgErr[ifit][ib]=ws[ifit][ib]->var("nBkg")->getError();
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
    h1_nSig[ifit] = new TH1D(Form("h1_nSig%ds_%d",states,ifit+1),Form("h1_nSig%ds;p_{T} (GeV/c);dN_{(#Upsilon%dS)}/dp_{T}",states,states),nBin,binArr); 
    h1_nBkg[ifit] = new TH1D(Form("h1_nBkg_%d",ifit+1),"h1_nBkg;p_{T} (GeV/c);dN_{Bkg}/dp_{T}",nBin,binArr); 
    for (int ib =0; ib < nBin; ib ++ ) {
      if(states ==1) { h1_nSig[ifit]->SetBinContent(ib+1,nSig1s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig1sErr[ifit][ib]);}   
      else if(states ==2) { h1_nSig[ifit]->SetBinContent(ib+1,nSig2s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig2sErr[ifit][ib]);}   
      else if(states ==3) { h1_nSig[ifit]->SetBinContent(ib+1,nSig3s[ifit][ib]); h1_nSig[ifit]->SetBinError(ib+1,nSig3sErr[ifit][ib]);}   
      h1_nBkg[ifit]->SetBinContent(ib+1,nBkg[ifit][ib]);   
      h1_nBkg[ifit]->SetBinError(ib+1,nBkgErr[ifit][ib]);   
    }
  }

  h1_nSig[0]->Scale(1./(lumi_pa_run1*1000));
  h1_nSig[1]->Scale(1./(lumi_pa_run2*1000));

  //// normalization
  for(int ifit=0; ifit<nFit; ifit++){
    TH1ScaleByWidth(h1_nSig[ifit]);
    TH1ScaleByWidth(h1_nBkg[ifit]);
    SetHistStyle(h1_nSig[ifit],ifit, ifit);
    SetHistStyle(h1_nBkg[ifit],ifit, ifit);
  }
  TGraphErrors *gSig[nFit];
  TGraphErrors *gBkg[nFit];
  double pxtmp, pytmp, extmp, eytmp;
  double shift_x=0.32;
  double shift_x_diff=shift_x*2/3;
  
  for(int i=0;i<nFit;i++){
    gSig[i] = new TGraphErrors(h1_nSig[i]);
    gBkg[i] = new TGraphErrors(h1_nBkg[i]);
    for(int j=0;j<gSig[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      gSig[i]->GetPoint(j,pxtmp,pytmp);
      extmp=gSig[i]->GetErrorX(j);
      eytmp=gSig[i]->GetErrorY(j);
      gSig[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      gSig[i]->SetPointError(j,0,eytmp);
    }
    for(int j=0;j<gBkg[i]->GetN();j++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      gBkg[i]->GetPoint(j,pxtmp,pytmp);
      extmp=gBkg[i]->GetErrorX(j);
      eytmp=gBkg[i]->GetErrorY(j);
      gBkg[i]->SetPoint(j,pxtmp-shift_x+shift_x_diff*i,pytmp);
      gBkg[i]->SetPointError(j,0,eytmp);
    }
  }
  int binmax = h1_nSig[0]->GetMaximumBin();
  double valmax = h1_nSig[0]->GetBinContent(binmax);
  int binmax_bkg = h1_nBkg[0]->GetMaximumBin();
  double valmax_bkg = h1_nBkg[0]->GetBinContent(binmax_bkg);
  cout << " binmax : " << binmax << endl;
  cout << " valmax : " << valmax << endl;
  gSig[0]->GetYaxis()->SetRangeUser(0,valmax*1.7);
  gBkg[0]->GetYaxis()->SetRangeUser(0,valmax_bkg*1.7);

  gSig[0]->GetYaxis()->SetTitle("B #frac{d#sigma}{ dp_{T}} (nb/ GeV/c)");
  gSig[0]->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  gSig[0]->GetXaxis()->SetRangeUser(0,30);
  gSig[0]->GetXaxis()->SetLimits(0,30);

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.035);

  double pos_y_diff = 0.05; 
  double pos_y = 0.837; 
  double pos_x = 0.51; 
  double line_x_diff = 0.50;
  
  double leg_posx1 = 0.47;
  double leg_posy1 = 0.57;
  double leg_posx2 = 0.83;
  double leg_posy2 = 0.67;

  TLegend* fitleg = new TLegend(leg_posx1,leg_posy1,leg_posx2,leg_posy2); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);

  TCanvas* c_nSigs = new TCanvas("c_nSigs","c_nSigs",600,600);
  c_nSigs->cd();
//  gPad->SetLogy();  // KTO for yield
  for(int ifit=0; ifit<nFit; ifit++){ gSig[0]->Draw("AP"); if(ifit>0) gSig[ifit]->Draw("P");}
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s #Upsilon(%dS)",szAA.Data(),states));
  for(int ifit=0;ifit<nFit; ifit++){fitleg->AddEntry(h1_nSig[ifit],Form("%s",fitName[ifit].Data()),"pe");}
  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_nSigs->SetTicks(1,1);
  c_nSigs->Modified();
  c_nSigs->Update();
  CMS_lumi_internal( c_nSigs, iPeriod, iPos );
  c_nSigs->SaveAs(Form("pt_nSig%ds_%s_constrain.pdf",states,szAA.Data()));
  
    //latex->SetTextColor(fitColorArr[ifit]);
    //latex->DrawLatex(0.55,0.81-pos_y_diff*(ifit+1),Form("%s fit",fitName[ifit].Data()));

  TCanvas* c_nBkg = new TCanvas("c_nBkg","c_nBkg",600,600);
  c_nBkg->cd();
//  gPad->SetLogy();  // KTO for yield
  for(int ifit=0; ifit<nFit; ifit++){ gBkg[0]->Draw("AP"); if(ifit>0) gBkg[ifit]->Draw("P");   }
  latex->SetTextColor(kBlack);
  latex->DrawLatex(pos_x,pos_y,Form("%s Background",szAA.Data()));
  fitleg->Draw("same");
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.1);
  c_nBkg->SetTicks(1,1);
  c_nBkg->Modified();
  c_nBkg->Update();
  CMS_lumi_internal( c_nBkg, iPeriod, iPos ); 
  c_nBkg->SaveAs(Form("pt_nBkg%ds_%s_constrain.pdf",states,szAA.Data()));
  
  return 0;
}

