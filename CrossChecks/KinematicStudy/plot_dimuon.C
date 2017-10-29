#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../cutsAndBin.h"
#include "../../PsetCollection.h"
#include "../../CMS_lumi.C"
#include "../../tdrstyle.C"
#include "../../SONGKYO.h"

using namespace std;
using namespace RooFit;

void plot_dimuon()
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;

  gStyle->SetEndErrorSize(0);

  float rapLow = -2.5; 
  float rapHigh = 2.5; 
  const int nRapBin = 50;

  float y;
  
  TFile *f1 = new TFile("/afs/cern.ch/work/j/jaebeom/public/RpPb5TeV/SkimmedFiles/yskimPA1st_OpSign_20177262037_unIdentified.root","read");
  TFile *f2 = new TFile("/afs/cern.ch/work/j/jaebeom/public/RpPb5TeV/SkimmedFiles/yskimPA2nd_OpSign_20177262044_unIdentified.root","read");
 
  float SimuPtCut = 4;
  float MassCutLow = 8;
  float MassCutHigh = 14;
  
  TString kineCut = Form("pt1>%.2f && pt2>%.2f && mass>%.f && mass<%.f", SimuPtCut, SimuPtCut, MassCutLow, MassCutHigh);
  
  TTree *tree1 = (TTree*) f1->Get("mm");
  TTree *tree2 = (TTree*) f2->Get("mm");
  
  RooDataSet *dataset1 = (RooDataSet*) f1->Get("dataset");
  RooDataSet *dataset2 = (RooDataSet*) f2->Get("dataset");

  RooWorkspace *ws1 = new RooWorkspace("workspace1");
  RooWorkspace *ws2 = new RooWorkspace("workspace2");

  ws1->import(*dataset1);
  ws2->import(*dataset2);

  cout << "<<<<<<<<<<<<<< 1st run workspace >>>>>>>>>>>>> " << endl;
  ws1->data("dataset")->Print();
  cout << endl;
  cout << "<<<<<<<<<<<<<< 2nd run workspace >>>>>>>>>>>>> " << endl;
  ws2->data("dataset")->Print();
  
  RooDataSet *reducedDS1 = (RooDataSet*)dataset1->reduce(RooArgSet(*(ws1->var("mass")), *(ws1->var("pt")), *(ws1->var("y"))), kineCut.Data() );
  RooDataSet *reducedDS2 = (RooDataSet*)dataset2->reduce(RooArgSet(*(ws2->var("mass")), *(ws2->var("pt")), *(ws2->var("y"))), kineCut.Data() );

  reducedDS1->SetName("reducedDS1");
  ws1->import(*reducedDS1);
  ws1->var("y")->setRange(rapLow, rapHigh);
  reducedDS2->SetName("reducedDS2");
  ws2->import(*reducedDS2);
  ws2->var("y")->setRange(rapLow, rapHigh);
  
  TCanvas* c1 = new TCanvas("canvas1","Plots",650,710);
  c1->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
  pad1->Range(-0.5526315,-0.2278481,2.078947,1.291139);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetLeftMargin(0.21);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.15);
  pad1->SetFrameLineColor(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameLineColor(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();

  RooPlot* myPlot1 = ws1->var("y")->frame(nRapBin);
  ws1->data("reducedDS1")->plotOn(myPlot1,Name("hist1"),MarkerSize(0.7),MarkerColor(kRed));
  RooPlot* myPlot2 = ws2->var("y")->frame(nRapBin);
  ws2->data("reducedDS2")->plotOn(myPlot2,Name("hist2"),MarkerSize(0.7),MarkerColor(kBlue));

  myPlot1->SetMinimum(0);
  myPlot1->SetMaximum(2000);
  myPlot1->SetStats(0);
  myPlot1->SetLineWidth(0);
  myPlot1->SetTitleSize(0);
  myPlot1->SetTitle(" ");
  myPlot1->GetXaxis()->SetTitle("y_{lab}");
  myPlot1->GetXaxis()->CenterTitle(true);
  myPlot1->GetXaxis()->SetLabelFont(42);
  myPlot1->GetXaxis()->SetLabelSize(0.035);
  myPlot1->GetXaxis()->SetTitleSize(0.0525);
  myPlot1->GetXaxis()->SetTitleOffset(1.25);
  myPlot1->GetXaxis()->SetTitleFont(42);
  myPlot1->GetYaxis()->SetTitle("Counts");
  myPlot1->GetYaxis()->CenterTitle(true);
  myPlot1->GetYaxis()->SetLabelFont(42);
  myPlot1->GetYaxis()->SetLabelSize(0.035);
  myPlot1->GetYaxis()->SetTitleSize(0.0525);
  myPlot1->GetYaxis()->SetTitleOffset(1.5);
  myPlot1->GetYaxis()->SetTitleFont(42);

  myPlot2->SetMinimum(0);
  myPlot2->SetMaximum(1000);
  myPlot2->SetStats(0);
  myPlot2->SetLineWidth(0);
  myPlot2->SetTitleSize(0);
  myPlot2->SetTitle(" ");
  myPlot2->GetXaxis()->SetTitle("y_{lab}");
  myPlot2->GetXaxis()->CenterTitle(true);
  myPlot2->GetXaxis()->SetLabelFont(42);
  myPlot2->GetXaxis()->SetLabelSize(0.035);
  myPlot2->GetXaxis()->SetTitleSize(0.0525);
  myPlot2->GetXaxis()->SetTitleOffset(1.25);
  myPlot2->GetXaxis()->SetTitleFont(42);
  myPlot2->GetYaxis()->SetTitle("Counts");
  myPlot2->GetYaxis()->CenterTitle(true);
  myPlot2->GetYaxis()->SetLabelFont(42);
  myPlot2->GetYaxis()->SetLabelSize(0.035);
  myPlot2->GetYaxis()->SetTitleSize(0.0525);
  myPlot2->GetYaxis()->SetTitleOffset(1.5);
  myPlot2->GetYaxis()->SetTitleFont(42);
  myPlot1->Draw();
  myPlot2->Draw("same");

  TLatex* tex = 0;
  tex = new TLatex(0.22,0.96,"#sqrt{s_{NN}} = 5.02 TeV");
  tex->SetNDC();
  tex->SetTextFont(43);
  tex->SetTextSize(30);
  tex->SetLineWidth(2);
  tex->Draw("same");
  tex = new TLatex(0.96,0.96,"pPb 34.6 nb^{-1}");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(43);
  tex->SetTextSize(30);
  tex->SetLineWidth(2);
  tex->Draw("same");
  tex = new TLatex(0.25,0.86,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.06);
  tex->SetLineWidth(2);
  tex->Draw("same");
  tex = new TLatex(0.38,0.86,"Internal");
  tex->SetNDC();
  tex->SetTextFont(52);
  tex->SetLineWidth(2);
  tex->Draw("same");

  float pos_text_x = 0.63;
  float pos_text_y = 0.796;
  float pos_y_diff = 0.056;
  float text_size = 19;
  int text_color = 1;
  
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", SimuPtCut ), pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("%.f < m_{#mu^{+}#mu^{-}} < %.f GeV/c^{2}", MassCutLow,MassCutHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  
  TLegend* fitleg = new TLegend(0.66,0.61,0.81,0.69); fitleg->SetTextSize(19);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot1->findObject("hist1"),"1st skim","pe");
  fitleg->AddEntry(myPlot2->findObject("hist2"),"2nd skim","pe");
  fitleg->Draw("same");

  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  c1->Update();
  c1->SaveAs("KinePlot.pdf");

  TCanvas *c2 = new TCanvas("c2", "",0,0,800,800);
  c2->Range(-0.5526315,-0.2278481,2.078947,1.291139);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetBorderSize(0);
  c2->SetTickx(1);
  c2->SetTicky(1);
  c2->SetLeftMargin(0.21);
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.06);
  c2->SetBottomMargin(0.15);
  c2->SetFrameLineColor(0);
  c2->SetFrameBorderMode(0);
  c2->SetFrameLineColor(0);
  c2->SetFrameBorderMode(0);

  TH1D *h1D = new TH1D("h1D","",nRapBin,rapLow,rapHigh);
  TH1D *h2D = new TH1D("h2D","",nRapBin,rapLow,rapHigh);
  h1D->Reset();
  h1D->SetMinimum(0);
  h1D->SetMaximum(2000);
  h1D->SetStats(0);
  //h1D->SetLineWidth(0);
  h1D->GetXaxis()->SetTitle("");
  h1D->GetXaxis()->CenterTitle(true);
  h1D->GetXaxis()->SetLabelFont(42);
  h1D->GetXaxis()->SetLabelSize(0.035);
  h1D->GetXaxis()->SetTitleSize(0.0525);
  h1D->GetXaxis()->SetTitleOffset(1.25);
  h1D->GetXaxis()->SetTitleFont(42);
  h1D->GetXaxis()->SetTitle("y_{lab}");
  h1D->GetYaxis()->CenterTitle(true);
  h1D->GetYaxis()->SetLabelFont(42);
  h1D->GetYaxis()->SetLabelSize(0.035);
  h1D->GetYaxis()->SetTitleSize(0.0525);
  h1D->GetYaxis()->SetTitleOffset(1.5);
  h1D->GetYaxis()->SetTitleFont(42);
  h2D->Reset();
  h2D->SetMinimum(0);
  h2D->SetMaximum(2000);
  //h2D->SetLineWidth(0);
  h2D->SetStats(0);
  h2D->GetXaxis()->SetTitle("");
  h2D->GetXaxis()->CenterTitle(true);
  h2D->GetXaxis()->SetLabelFont(42);
  h2D->GetXaxis()->SetLabelSize(0.035);
  h2D->GetXaxis()->SetTitleSize(0.0525);
  h2D->GetXaxis()->SetTitleOffset(1.25);
  h2D->GetXaxis()->SetTitleFont(42);
  h2D->GetXaxis()->SetTitle("y_{lab}");
  h2D->GetYaxis()->CenterTitle(true);
  h2D->GetYaxis()->SetLabelFont(42);
  h2D->GetYaxis()->SetLabelSize(0.035);
  h2D->GetYaxis()->SetTitleSize(0.0525);
  h2D->GetYaxis()->SetTitleOffset(1.5);
  h2D->GetYaxis()->SetTitleFont(42);
 
  tree1->Draw("y>>h1D","pt1>4&&pt2>4 && mass<14 && mass>8",""); 
  tree2->Draw("y>>h2D","pt1>4&&pt2>4 && mass<14 && mass>8",""); 
  c2->cd();
  h1D->SetMarkerColor(kRed);
  h1D->SetMarkerStyle(20);
  h1D->SetLineColor(kBlack);
  h2D->SetMarkerColor(kBlue);
  h2D->SetMarkerStyle(20);
  h2D->SetLineColor(kBlack);
  h1D->Scale(1/h1D->Integral());
  h2D->Scale(1/h2D->Integral());
  h1D->GetYaxis()->SetRangeUser(0,0.06);
  h2D->GetYaxis()->SetRangeUser(0,0.06);
  h1D->Draw("e");
  h2D->Draw("e same");
  drawText(Form("p_{T}^{#mu} > %.f GeV/c", SimuPtCut ), pos_text_x,pos_text_y,text_color,text_size);
  drawText(Form("%.f < m_{#mu^{+}#mu^{-}} < %.f GeV/c^{2}", MassCutLow,MassCutHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);

  TLatex* tex2 = 0;
  tex2 = new TLatex(0.22,0.96,"#sqrt{s_{NN}} = 5.02 TeV");
  tex2->SetNDC();
  tex2->SetTextFont(43);
  tex2->SetTextSize(30);
  tex2->SetLineWidth(2);
  tex2->Draw("same");
  tex2 = new TLatex(0.96,0.96,"pPb 34.6 nb^{-1}");
  tex2->SetNDC();
  tex2->SetTextAlign(31);
  tex2->SetTextFont(43);
  tex2->SetTextSize(30);
  tex2->SetLineWidth(2);
  tex2->Draw("same");
  tex2 = new TLatex(0.25,0.86,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.06);
  tex2->SetLineWidth(2);
  tex2->Draw("same");
  tex2 = new TLatex(0.38,0.86,"Internal");
  tex2->SetNDC();
  tex2->SetTextFont(52);
  tex2->SetLineWidth(2);
  tex2->Draw("same");

  fitleg->Draw("same");
  
  c2->Modified();
  c2->SetSelected(c2);
  c2->Update();
  c2->SaveAs("KinePlot_tree.pdf");
}

