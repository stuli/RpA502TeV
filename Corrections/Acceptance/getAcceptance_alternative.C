#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
using namespace std;

int kNoWeight=0;
int kPP = 1;
int kAA = 2;

int k1S = 1; 
int k2S = 2; 
int k3S = 3; 

int kNoVar = 0;
int kPtPlus = 1;
int kPtMinus = 2;
int kYPlus = 3;
int kYMinus = 4;
//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
double pTAcc1[10][10] = {0.0}; // Upsilon 1S acceptance for n-th pT and n-th y bins 
double pTAcc2[10][10] = {0.0}; // Upsilon 2S acceptance for n-th pT and n-th y bins
double pTAcc3[10][10] = {0.0}; // Upsilon 3S acceptance for n-th pT and n-th y bins
double pTErr1[10][10] = {0.0}; // Upsilon 1S acceptance error for n-th pT and n-th y bins 
double pTErr2[10][10] = {0.0}; // Upsilon 2S acceptance error for n-th pT and n-th y bins
double pTErr3[10][10] = {0.0}; // Upsilon 3S acceptance error for n-th pT and n-th y bins


TH1F *hAccPt1;
TH1F *hAccPt2;
TH1F *hAccPt3;

TH1F *hAccY1;
TH1F *hAccY2;
TH1F *hAccY3;

TH1F *hPadPt = new TH1F("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1F *hPadY = new TH1F("hPadY",";|y|;Acceptance",10,0.0,2.4);

TLatex *lt1 = new TLatex();

float getAcceptanceSingleBin(int state=k1S, int icoll=kNoWeight, TString kineCut="", int varianceID=0);
void getAcceptance(int bin = 1, int doDraw = 1) { // doDraw == 1 : draw, == 0 : no draw
  gStyle->SetOptStat(0);

  hPadPt->GetYaxis()->CenterTitle();
  hPadPt->GetYaxis()->SetTitleOffset(1.3);
  hPadPt->GetXaxis()->CenterTitle();
  hPadY->GetYaxis()->CenterTitle();
  hPadY->GetYaxis()->SetTitleOffset(1.3);
  hPadY->GetXaxis()->CenterTitle();
  lt1->SetNDC();

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample
  int nPtBins = 1;  double* ptBin = NULL;  int nYBins = 1;  double *yBin = NULL;
  const int nPtBins1 = 1;   double ptBin1[nPtBins1+1] = {0,30};  const int nYBins1 = 1;   double yBin1[nYBins1+1] = {0,2.4};
  const int nPtBins2 = 6;   double ptBin2[nPtBins2+1] = {0,2,4,6,9,12,30};  const int nYBins2 = 1;   double yBin2[nYBins2+1] = {0,2.4};
  //const int nPtBins2 = 3;   double ptBin2[nPtBins2+1] = {0,5,12,30};  const int nYBins2 = 1;   double yBin2[nYBins2+1] = {0,2.4};
  const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,30};  const int nYBins3 = 6;   double yBin3[nYBins3+1] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  //const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,30};  const int nYBins3 = 2;   double yBin3[nYBins3+1] = {0, 1.2, 2.4};
  
  if ( bin==1 ) {
    nPtBins = nPtBins1;  ptBin=ptBin1;     nYBins = nYBins1; yBin=yBin1;
  }
  else if ( bin==2 ) {
    nPtBins = nPtBins2;  ptBin=ptBin2;     nYBins = nYBins2; yBin=yBin2;
  }
  else if ( bin==3 ) {
    nPtBins = nPtBins3;  ptBin=ptBin3;     nYBins = nYBins3; yBin=yBin3;
  }

  hAccPt1 = new TH1F("hAccPt1",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,0.0,30.0);
  hAccPt2 = new TH1F("hAccPt2",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,0.0,30.0);
  hAccPt3 = new TH1F("hAccPt3",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,0.0,30.0);

  hAccY1 = new TH1F("hAccY1",";|y|;Acceptance of #Upsilon(1S)",nYBins,0.0,2.4);
  hAccY2 = new TH1F("hAccY2",";|y|;Acceptance of #Upsilon(2S)",nYBins,0.0,2.4);
  hAccY3 = new TH1F("hAccY3",";|y|;Acceptance of #Upsilon(3S)",nYBins,0.0,2.4);

  hAccPt1->Sumw2();
  hAccPt2->Sumw2();
  hAccPt3->Sumw2();
  hAccY1->Sumw2();
  hAccY2->Sumw2();
  hAccY3->Sumw2();

  hAccPt1->SetMarkerStyle(20);hAccPt1->SetMarkerColor(kBlue+2);
  hAccPt2->SetMarkerStyle(20);hAccPt2->SetMarkerColor(kBlue+2);
  hAccPt3->SetMarkerStyle(20);hAccPt3->SetMarkerColor(kBlue+2);

  hAccY1->SetMarkerStyle(20);hAccY1->SetMarkerColor(kBlue+2);
  hAccY2->SetMarkerStyle(20);hAccY2->SetMarkerColor(kBlue+2);
  hAccY3->SetMarkerStyle(20);hAccY3->SetMarkerColor(kBlue+2);

  for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) { 
    for ( int iy =1 ; iy<=nYBins ; iy++) {
      TString kRange = Form("pt>%0.1f && pt<%0.1f && abs(y)>%0.1f && abs(y)<%0.1f",(float)ptBin[ipt-1],(float)ptBin[ipt], (float)yBin[iy-1], (float)yBin[iy] ) ;
      cout << "*===*===*===*===*  " << kRange << "  *===*===*===*===*" << endl;
      float y1s_noRwt = getAcceptanceSingleBin( 1, kNoWeight, kRange,kNoVar);
      float y1s_ppRwt = getAcceptanceSingleBin( 1, kPP, kRange,kNoVar);
      float y1s_ppRwt_ptPlus  = getAcceptanceSingleBin( 1, kPP, kRange, kPtPlus);
      float y1s_ppRwt_ptMinus = getAcceptanceSingleBin( 1, kPP, kRange, kPtMinus);
      float y1s_ppRwt_yPlus  = getAcceptanceSingleBin( 1, kPP, kRange, kYPlus);
      float y1s_ppRwt_yMinus = getAcceptanceSingleBin( 1, kPP, kRange, kYMinus);

      float y2s_noRwt = getAcceptanceSingleBin( 2, kNoWeight, kRange,kNoVar);
      float y2s_ppRwt = getAcceptanceSingleBin( 2, kPP, kRange,kNoVar);
      float y2s_ppRwt_ptPlus  = getAcceptanceSingleBin( 2, kPP, kRange, kPtPlus);
      float y2s_ppRwt_ptMinus = getAcceptanceSingleBin( 2, kPP, kRange, kPtMinus);
      float y2s_ppRwt_yPlus  = getAcceptanceSingleBin( 2, kPP, kRange, kYPlus);
      float y2s_ppRwt_yMinus = getAcceptanceSingleBin( 2, kPP, kRange, kYMinus);

      float y3s_noRwt = getAcceptanceSingleBin( 3, kNoWeight, kRange,kNoVar);
      float y3s_ppRwt = getAcceptanceSingleBin( 3, kPP, kRange,kNoVar);
      float y3s_ppRwt_ptPlus  = getAcceptanceSingleBin( 3, kPP, kRange, kPtPlus);
      float y3s_ppRwt_ptMinus = getAcceptanceSingleBin( 3, kPP, kRange, kPtMinus);
      float y3s_ppRwt_yPlus  = getAcceptanceSingleBin( 3, kPP, kRange, kYPlus);
      float y3s_ppRwt_yMinus = getAcceptanceSingleBin( 3, kPP, kRange, kYMinus);

      float err1s_ppPt = max ( abs(y1s_ppRwt_ptPlus/y1s_ppRwt-1), abs(y1s_ppRwt_ptMinus/y1s_ppRwt-1) ) ;
      float err2s_ppPt = max ( abs(y2s_ppRwt_ptPlus/y2s_ppRwt-1), abs(y2s_ppRwt_ptMinus/y2s_ppRwt-1) ) ;
      float err3s_ppPt = max ( abs(y3s_ppRwt_ptPlus/y3s_ppRwt-1), abs(y3s_ppRwt_ptMinus/y3s_ppRwt-1) ) ;
      float err1s_ppY = max ( abs(y1s_ppRwt_yPlus/y1s_ppRwt-1), abs(y1s_ppRwt_yMinus/y1s_ppRwt-1) ) ;
      float err2s_ppY = max ( abs(y2s_ppRwt_yPlus/y2s_ppRwt-1), abs(y2s_ppRwt_yMinus/y2s_ppRwt-1) ) ;
      float err3s_ppY = max ( abs(y3s_ppRwt_yPlus/y3s_ppRwt-1), abs(y3s_ppRwt_yMinus/y3s_ppRwt-1) ) ;

      float y1s_aaRwt = getAcceptanceSingleBin( 1, kAA, kRange,kNoVar);
      float y1s_aaRwt_ptPlus  = getAcceptanceSingleBin( 1, kAA, kRange, kPtPlus);
      float y1s_aaRwt_ptMinus = getAcceptanceSingleBin( 1, kAA, kRange, kPtMinus);
      float y1s_aaRwt_yPlus  = getAcceptanceSingleBin( 1, kAA, kRange, kYPlus);
      float y1s_aaRwt_yMinus = getAcceptanceSingleBin( 1, kAA, kRange, kYMinus);

      float y2s_aaRwt = getAcceptanceSingleBin( 2, kAA, kRange,kNoVar);
      float y2s_aaRwt_ptPlus  = getAcceptanceSingleBin( 2, kAA, kRange, kPtPlus);
      float y2s_aaRwt_ptMinus = getAcceptanceSingleBin( 2, kAA, kRange, kPtMinus);
      float y2s_aaRwt_yPlus  = getAcceptanceSingleBin( 2, kAA, kRange, kYPlus);
      float y2s_aaRwt_yMinus = getAcceptanceSingleBin( 2, kAA, kRange, kYMinus);

      float y3s_aaRwt = getAcceptanceSingleBin( 3, kAA, kRange,kNoVar);
      float y3s_aaRwt_ptPlus  = getAcceptanceSingleBin( 3, kAA, kRange, kPtPlus);
      float y3s_aaRwt_ptMinus = getAcceptanceSingleBin( 3, kAA, kRange, kPtMinus);
      float y3s_aaRwt_yPlus  = getAcceptanceSingleBin( 3, kAA, kRange, kYPlus);
      float y3s_aaRwt_yMinus = getAcceptanceSingleBin( 3, kAA, kRange, kYMinus);

      float err1s_aaPt = max ( abs(y1s_aaRwt_ptPlus/y1s_aaRwt-1), abs(y1s_aaRwt_ptMinus/y1s_aaRwt-1) ) ;
      float err2s_aaPt = max ( abs(y2s_aaRwt_ptPlus/y2s_aaRwt-1), abs(y2s_aaRwt_ptMinus/y2s_aaRwt-1) ) ;
      float err3s_aaPt = max ( abs(y3s_aaRwt_ptPlus/y3s_aaRwt-1), abs(y3s_aaRwt_ptMinus/y3s_aaRwt-1) ) ;
      float err1s_aaY = max ( abs(y1s_aaRwt_yPlus/y1s_aaRwt-1), abs(y1s_aaRwt_yMinus/y1s_aaRwt-1) ) ;
      float err2s_aaY = max ( abs(y2s_aaRwt_yPlus/y2s_aaRwt-1), abs(y2s_aaRwt_yMinus/y2s_aaRwt-1) ) ;
      float err3s_aaY = max ( abs(y3s_aaRwt_yPlus/y3s_aaRwt-1), abs(y3s_aaRwt_yMinus/y3s_aaRwt-1) ) ;

      pTAcc1[ipt][iy] = y1s_noRwt;
      pTAcc2[ipt][iy] = y2s_noRwt;
      pTAcc3[ipt][iy] = y3s_noRwt;
      pTErr1[ipt][iy] = err1s_ppPt;
      pTErr2[ipt][iy] = err2s_ppPt;
      pTErr3[ipt][iy] = err3s_ppPt;
      cout<<"dmoon chk 1 : "<<pTAcc1[ipt][iy]<<endl;
      
      cout << "Un-weighed acceptance : (1S) = " << y1s_noRwt << ",   (2S) = " << y2s_noRwt << ",   (3S) = "<< y3s_noRwt<<endl;
      cout << "Re-weighed Acceptance : " << "pp(1S) = " << y1s_ppRwt << ",  PbPb(1S) = " << y1s_aaRwt ;
      cout << " pp(2S) = " << y2s_ppRwt << ",  PbPb(2S) = " << y2s_aaRwt ;
      cout << " pp(3S) = " << y3s_ppRwt << ",  PbPb(3S) = " << y3s_aaRwt ;
      //cout << ",     (2S/1S,AA/pp) acc. ratio = " <<  (y2s_aaRwt/y2s_ppRwt) /(y1s_aaRwt/y1s_ppRwt)  << endl;
      cout<<""<<endl;
      
      cout << int(1000*y1s_noRwt)*0.001 << " & " << int(1000*y1s_ppRwt)*0.001<< " & " << int(1000*y1s_aaRwt)*0.001<< " & " ;
      cout << int(1000*y2s_noRwt)*0.001 << " & " << int(1000*y2s_ppRwt)*0.001<< " & " << int(1000*y2s_aaRwt)*0.001 << " & ";
      cout << int(1000*y3s_noRwt)*0.001 << " & " << int(1000*y3s_ppRwt)*0.001<< " & " << int(1000*y3s_aaRwt)*0.001 << " & ";
      cout <<""<<endl;
      //cout << int(1000* (y2s_aaRwt/y2s_ppRwt) / (y1s_aaRwt/y1s_ppRwt) )*0.001 << endl; // Final correction
      cout << int(1000*err1s_ppPt)*0.1 << "\% & " << int(1000*err1s_aaPt)*0.1 << "\% & " ;
      cout << int(1000*err2s_ppPt)*0.1 << "\% & " << int(1000*err2s_aaPt)*0.1 << "\% & " ;
      cout << int(1000*err3s_ppPt)*0.1 << "\% & " << int(1000*err3s_aaPt)*0.1 << "\% & " ;
      cout << int(1000*err1s_ppY)*0.1 << "\% & " << int(1000*err1s_aaY)*0.1 << "\% & " ;
      cout << int(1000*err2s_ppY)*0.1 << "\% & " << int(1000*err2s_aaY)*0.1 << "\% & " ;
      cout << int(1000*err3s_ppY)*0.1 << "\% & " << int(1000*err3s_aaY)*0.1 << "\% & " ;
      //cout << int(1000*sqrt( err1s_ppPt*err1s_ppPt +  err2s_ppPt* err2s_ppPt + err1s_ppY*err1s_ppY + err2s_ppY*err2s_ppY + err1s_aaPt*err1s_aaPt +  err2s_aaPt* err2s_aaPt + err1s_aaY*err1s_aaY + err2s_aaY*err2s_aaY ) ) *0.1 << "\%" << endl;
      cout<<""<<endl;

      if(bin == 1 || bin == 2){
        hAccPt1->SetBinContent(ipt,pTAcc1[ipt][iy]);
        hAccPt2->SetBinContent(ipt,pTAcc2[ipt][iy]);
        hAccPt3->SetBinContent(ipt,pTAcc3[ipt][iy]);
        hAccPt1->SetBinError(ipt,pTErr1[ipt][iy]);
        hAccPt2->SetBinError(ipt,pTErr2[ipt][iy]);
        hAccPt3->SetBinError(ipt,pTErr3[ipt][iy]);
        cout<<"dmoon chk 2 : "<<pTAcc1[ipt][iy]<<endl;
      }

      if(bin == 3){
        hAccY1->SetBinContent(iy,pTAcc1[ipt][iy]);
        hAccY2->SetBinContent(iy,pTAcc2[ipt][iy]);
        hAccY3->SetBinContent(iy,pTAcc3[ipt][iy]);
        hAccY1->SetBinError(iy,pTErr1[ipt][iy]);
        hAccY2->SetBinError(iy,pTErr2[ipt][iy]);
        hAccY3->SetBinError(iy,pTErr3[ipt][iy]);
      }
    } // ipt
  } // iy
  if(doDraw == 1){
    TCanvas *c1 = new TCanvas("c1","",1800,600);
    c1->Divide(3,1);
    if(bin == 1 || bin == 2){
      c1->cd(1); hPadPt->Draw(); hAccPt1->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (1S)");
      c1->cd(2); hPadPt->Draw(); hAccPt2->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (2S)");
      c1->cd(3); hPadPt->Draw(); hAccPt3->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (3S)");
      if(bin == 1) c1->SaveAs("plot_acc_int.png");
      if(bin == 2) c1->SaveAs("plot_acc_pT.png");
    }
    if(bin == 3){
      c1->cd(1); hPadY->Draw(); hAccY1->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (1S)");
      c1->cd(2); hPadY->Draw(); hAccY2->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (2S)");
      c1->cd(3); hPadY->Draw(); hAccY3->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (3S)");
      c1->SaveAs("plot_acc_y.png");
    }
  }
}

float getAcceptanceSingleBin(int state, int icoll, TString kineCut, int varianceID) {
  
  TH1::SetDefaultSumw2();
  
  float muPtCut = 4; // for acceptance 
  float muEtaCut = 2.4;
  
  TH1D* hGen; // in rapidity bins, centrality bins
  TH1D* hGenAcc; // in rapidity bins, centrality bins
  TH1D* hAcc; // Acceptance rate : in rapidity bins, centrality bins 
  
  TFile* f = NULL;
  if ( state == k1S){ 
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016106221_unIdentified.root");
  } else if ( state == k2S) {
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016106227_unIdentified.root");
  } else if ( state == k3S) {
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161062212_unIdentified.root");
  }

  TTree* mmGen = (TTree*)f->Get("mmGen");
  TString accCut = "pt1>4 && pt2>4 && abs(eta1)<2.4 && abs(eta2)<2.4";
  
  TString ptWeight;
  if ( icoll == kNoWeight ) ptWeight = "1";
  else if ( icoll == kPP ) { 
    if (state == k1S)  ptWeight = "1.09 - 0.022*pt";
    if (state == k2S)  ptWeight = "0.69 + 0.073*pt";
    if (state == k3S)  ptWeight = "0.69 + 0.073*pt"; // k3S check
  }    
  else if ( icoll == kAA ) { 
    if (state == k1S)  ptWeight = "1.08 - 0.015*pt";
    if (state == k2S)  ptWeight = "0.68 + 0.1*pt";
    if (state == k3S)  ptWeight = "0.69 + 0.073*pt"; // k3S check
  }    

  if ( varianceID == kPtPlus ) 
    ptWeight = "("+ptWeight+")* (0.8 + 0.0133*pt)" ;
  else if ( varianceID == kPtMinus ) 
    ptWeight = "("+ptWeight+")* (1.2 - 0.0133*pt)" ;
  else if ( varianceID == kYPlus ) 
    ptWeight = "("+ptWeight+")* (0.8 + 0.167*abs(y))" ;
  else if ( varianceID == kYMinus ) 
    ptWeight = "("+ptWeight+")* (1.2 - 0.167*abs(y))" ;
  //  cout << " ptWeight = " << ptWeight << endl;
  
  hGen = new TH1D(Form("hGen_state%d_icoll%d",state,icoll),";Rapidity;",1,-2.4,2.4);
  hGenAcc = (TH1D*)hGen->Clone(Form("hGenAcc_state%d_icoll%d",state,icoll));
  //// all GEN
  mmGen->Draw(Form("y>>%s",hGen->GetName()), Form("(%s)*(%s)", kineCut.Data(), ptWeight.Data()));
  //// GEN with sglMuAccCut 
  mmGen->Draw(Form("y>>%s",hGenAcc->GetName()), Form("((%s)&&(%s))*(%s)", kineCut.Data(), accCut.Data(), ptWeight.Data()));
  handsomeTH1(hGenAcc,1);
  //// --- Acceptance rate : 
  hAcc = (TH1D*)hGenAcc->Clone(Form("hAcc_state%d_icoll%d",state,icoll));
  hAcc->Divide(hGen);
  
  
  /*
  TCanvas* c1 = new TCanvas("c1","",400,400);
  hGen->Draw("hist");
  hGenAcc->Draw("same e");
  
  TCanvas* c2 = new TCanvas("c2","",400,400);
  hAcc->SetAxisRange(0,1,"Y");
  hAcc->Draw();
  */
    
  return hAcc->GetBinContent(1);
}


