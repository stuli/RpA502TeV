#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "cutsAndBin.h"
#include "multiTreeUtil.h"
using namespace std;

TString ResultDir  = "nominalFits";


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0) ;

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 

void getSpectra_RpPb_highPt(int state = 1) {  

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  int nPtBins=0;
  int nYBins=0;
  int nYBins_cr=0;
  double* ptBin;
  double* yBin;
  double* yBin_cr;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;
  //  double nColl1s[nCentBins] = {1819,1432,1005,606,349,186,90.7,40.1,7.67};
  //  double nPart1s[nCentBins] = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};



  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; nYBins_cr = nYBins1S_cr;   yBin_cr = yBin1S_cr;
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S;  nYBins_cr = nYBins2S_cr;   yBin_cr = yBin2S_cr;
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S;   nYBins_cr = nYBins3S_cr;   yBin_cr = yBin3S_cr;
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  

  // signals :
  TH1D* hRPAraw_rap;   // w/o efficiency correction
  TH1D* hRPAraw_int;   // w/o efficiency correction
  TH1D* hRPA_rap;   // w efficiency correction
  TH1D* hRPA_int;   // w efficiency correction
  TH1D* hRPAraw_pt;   // w/o efficiency correction
  TH1D* hRPA_pt;   // w efficiency correction

  TH1D* hrapAccPA;
  TH1D* hrapAccPA_cross;
  TH1D* hrapAccPP;
  TH1D* hptAccPA;
  TH1D* hptAccPP;
  TH1D* hintAccPA;
  TH1D* hintAccPP;

  TH1D* hrapEffPA;
  TH1D* hrapEffPA_cross;
  TH1D* hrapEffPP;
  TH1D* hptEffPA;
  TH1D* hptEffPP;
  TH1D* hintEffPA;
  TH1D* hintEffPP;

  TH1D* hrapEffPA_gen;
  TH1D* hrapEffPA_gen_cross;
  TH1D* hrapEffPP_gen;
  TH1D* hptEffPA_gen;
  TH1D* hptEffPP_gen;
  TH1D* hintEffPA_gen;
  TH1D* hintEffPP_gen;

  TFile* infacc = new TFile(Form("Acceptance/acceptance_wgt_%dS_20171121.root",state),"read");
  hrapAccPA  = (TH1D*)infacc->Get(Form("hrapAccPA%dS",state));
  hrapAccPP  = (TH1D*)infacc->Get(Form("hrapAccPP%dS",state));
  hptAccPA  = (TH1D*) infacc->Get(Form("hptAccPA%dS",state));
  hptAccPP  = (TH1D*) infacc->Get(Form("hptAccPP%dS",state));
  hintAccPA  = (TH1D*) infacc->Get(Form("hIntAccPA%dS",state));
  hintAccPP  = (TH1D*) infacc->Get(Form("hIntAccPP%dS",state));
  hrapAccPA_cross = (TH1D*) infacc->Get(Form("hrapAccXsPA%dS",state)); 

  //TFile* infeff_pPb = new TFile(Form("Efficiency_rootfiles/pPb/Eff_pPb_%dS_8_22_NewPtReweights.root",state),"read");
  TFile* infeff_pPb = new TFile(Form("Efficiency_rootfiles/pPb/Eff_pPb_%dS_11_20_NewRpABin.root",state),"read");
  TFile* infeff_pp = new TFile(Form("Efficiency_rootfiles/pp/Eff_pp_%dS_11_20_NewRpABin.root",state),"read");
  hrapEffPA  = (TH1D*)infeff_pPb->Get("RecoEventsRap");
  hrapEffPA_gen  = (TH1D*)infeff_pPb->Get("GenEventsRap");
  hptEffPA  = (TH1D*) infeff_pPb->Get("RecoEventsPtRpA");
  hptEffPA_gen  = (TH1D*) infeff_pPb->Get("GenEventsPtRpA");
  hintEffPA  = (TH1D*) infeff_pPb->Get("RecoEventsIntRpA");
  hintEffPA_gen  = (TH1D*) infeff_pPb->Get("GenEventsIntRpA");
  hrapEffPP  = (TH1D*)infeff_pp->Get("RecoEventsRap");
  hrapEffPP_gen  = (TH1D*)infeff_pp->Get("GenEventsRap");
  hptEffPP  = (TH1D*) infeff_pp->Get("RecoEventsPt");
  hptEffPP_gen  = (TH1D*) infeff_pp->Get("GenEventsPt");
  hintEffPP  = (TH1D*) infeff_pp->Get("RecoEventsInt");
  hintEffPP_gen  = (TH1D*) infeff_pp->Get("GenEventsInt");

  hrapEffPA->Divide(hrapEffPA_gen);
  hrapEffPP->Divide(hrapEffPP_gen);
  hptEffPA->Divide(hptEffPA_gen);
  hptEffPP->Divide(hptEffPP_gen);
  hintEffPA->Divide(hintEffPA_gen);
  hintEffPP->Divide(hintEffPP_gen);


  stripErrorBars(hrapAccPA_cross);
  stripErrorBars(hrapAccPA);
  stripErrorBars(hrapAccPP);
  stripErrorBars(hptAccPA);
  stripErrorBars(hptAccPP);
  stripErrorBars(hintAccPA);
  stripErrorBars(hintAccPP);

  stripErrorBars(hrapEffPA);
  stripErrorBars(hrapEffPP);
  stripErrorBars(hptEffPA);
  stripErrorBars(hptEffPP);
  stripErrorBars(hintEffPA);
  stripErrorBars(hintEffPP);

  hrapEffPA_cross = (TH1D*) hrapEffPA->Clone("hrapEffPA_Cross");
  TH1D* hrapEffPA_pa = (TH1D*) hrapAccPA->Clone("hrapEff_pa");
  TH1D* hrapEffPA_pp = (TH1D*) hrapAccPP->Clone("hrapEff_pp");
  hrapEffPA_pp->Reset();
  cout << " hrapEffPA_pp->GetNbinsX() : " << hrapEffPA_pp->GetNbinsX() << endl;
  cout << "hrapEffPP->GetNbinsX() : " << hrapEffPP->GetNbinsX() << endl;
  for(int i=1;i<=hrapEffPA_pp->GetNbinsX();i++){
    cout << "i : " << i << endl;
    if(i<((double)hrapEffPA_pp->GetNbinsX()/2+1)) hrapEffPA_pp->SetBinContent(i,hrapEffPP->GetBinContent(hrapEffPP->GetNbinsX()-i+1));
    else if(i>((double)hrapEffPA_pp->GetNbinsX()/2+1)) hrapEffPA_pp->SetBinContent(i,hrapEffPP->GetBinContent(i-hrapEffPP->GetNbinsX()));
    cout << " Bin : " << hrapEffPA_pp->GetBinContent(i) << endl;
    if(i==5) cout << " hrapEffPP->GetBinContent(hrapEffPP->GetNbinsX()-i+1) : " << hrapEffPP->GetBinContent(hrapEffPP->GetNbinsX()-i+1) << endl;
    cout << "hrapEffPA_pp->GetNbinsX()/2+1 " << hrapEffPA_pp->GetNbinsX()/2+1 << endl;
  }
  for(int i=1; i<=nYBins; i++){
    hrapEffPA_pa->SetBinContent(i,hrapEffPA->GetBinContent(i+1));
    hrapEffPA_pa->SetBinError(i,0);
  }

  TH1D* hrapSigPP = (TH1D*) hrapAccPP -> Clone("hrapPP");
  TH1D* hrapSigPA = (TH1D*) hrapAccPA -> Clone("hrapPA");
  TH1D* hrapSigPA_cross = (TH1D*) hrapAccPA_cross->Clone("hrapPA_Cross");
  TH1D* hintSigPP = (TH1D*) hintAccPA -> Clone("hintPP");
  TH1D* hintSigPA = (TH1D*) hintAccPA -> Clone("hintPA");
  TH1D* hptSigPP = (TH1D*) hptAccPP -> Clone("hptPP");
  TH1D* hptSigPA = (TH1D*) hptAccPA -> Clone("hptPA");
  hptSigPP->Reset();
  hptSigPA->Reset();
  hrapSigPA_cross->Reset();
  hrapSigPP->Reset();
  hrapSigPA->Reset();
  hintSigPP->Reset();
  hintSigPA->Reset();



  //***Int***
  TCanvas* c_int = new TCanvas("c_int","",400,400);
  valErr yield_intPP = getYield(state,kPPDATA,0,30,0.00,1.93,0,200,0,100);
  valErr yield_intPA = getYield(state,kPADATA,0,30,-1.93,1.93,0,200,0,100);
  hintSigPP->SetBinContent(1,yield_intPP.val/2);
  hintSigPP->SetBinError(1,yield_intPP.err/2);
  hintSigPA->SetBinContent(1,yield_intPA.val);
  hintSigPA->SetBinError(1,yield_intPA.err);
  hintSigPP->SetAxisRange(10,1e5,"Y");
  hintSigPA->SetAxisRange(10,1e5,"Y");
  handsomeTH1(hintSigPP,2);
  handsomeTH1(hintSigPA,2);
  hintSigPP->SetMarkerStyle(24);

  valErr yieldPP;

  //***yCM***
  TCanvas* c_rap =  new TCanvas("c_rap","",400,400);
  for ( int irap = 1 ; irap<= nYBins ; irap++) {
    valErr yieldPP;
    if(irap <= (nYBins/2+1)) {
      yieldPP = getYield(state, kPPDATA, 6,30, TMath::Abs(yBin[irap]), TMath::Abs(yBin[irap-1]), 0,200,0,100);
    }
    else if(irap > (nYBins/2+1)) {
      cout << "irap : " << irap << endl;
      yieldPP = getYield(state, kPPDATA, 6,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    }
    valErr yieldPA = getYield(state, kPADATA, 6,30, yBin[irap-1], yBin[irap], 0,200,0,100);
    hrapSigPP->SetBinContent( irap, yieldPP.val/2 ) ;
    hrapSigPP->SetBinError( irap, yieldPP.err/2 ) ;
    hrapSigPA->SetBinContent( irap, yieldPA.val ) ;
    hrapSigPA->SetBinError( irap, yieldPA.err ) ;
  }

  for(int irap = 1; irap<=nYBins_cr; irap++){
    cout << "nYBins_cr : " << nYBins_cr<< endl;
    cout << "yBin_cr[irap-1] to yBin_cr[irap] : " << yBin_cr[irap-1] << " - " << yBin_cr[irap] << endl;
    valErr yieldPA = getYield(state, kPADATA, 6,30, yBin_cr[irap-1],yBin_cr[irap],0,200,0,100);
    cout << "yieldPA.val : " << yieldPA.val << endl;
    cout << "yieldPA.err : " << yieldPA.err << endl;
    hrapSigPA_cross->SetBinContent(irap, yieldPA.val);
    hrapSigPA_cross->SetBinError(irap, yieldPA.err);
  }


  //yCM Yield 
  hrapSigPP->SetAxisRange(10,1e5,"Y");
  hrapSigPA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hrapSigPP,2);
  handsomeTH1(hrapSigPA,2);
  hrapSigPP->SetMarkerStyle(24);
 
  //***pt***
  TCanvas* c_pt =  new TCanvas("c_pt","",400,400);
  for ( int ipt = 1 ; ipt<= nPtBins ; ipt++) {
    valErr yieldPP = getYield(state, kPPDATA, ptBin[ipt-1],ptBin[ipt], 0.00,1.93, 0,200,0,100);
    valErr yieldPA = getYield(state, kPADATA, ptBin[ipt-1],ptBin[ipt], 0.00,1.93, 0,200,0,100);
    hptSigPA->SetBinContent( ipt, yieldPA.val ) ;
    hptSigPA->SetBinError( ipt, yieldPA.err ) ;
    hptSigPP->SetBinContent( ipt, yieldPP.val/2 ) ;
    hptSigPP->SetBinError( ipt, yieldPP.err/2 ) ;
  }

  //pt Yield
  hptSigPP->SetAxisRange(10,1e5,"Y");
  hptSigPA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hptSigPP,2);
  handsomeTH1(hptSigPA,2);
  hptSigPP->SetMarkerStyle(24);


  double Factor_Scale = lumi_pp*1000/(lumi_pa*208);

 //*****Draw int*****
  c_int->cd();
  hintSigPP->Draw();
  gPad->SetLogy();
  hintSigPA->Draw("same");

  hRPAraw_int = (TH1D*)hintSigPA->Clone("rpa_vs_int_raw");
  hRPAraw_int -> Divide(hintSigPP);
  TH1D* hrel_Acc_int = (TH1D*) hintAccPA -> Clone("hrel_Acc_int");
  TH1D* hrel_Eff_int;
  if(state!=3) hrel_Eff_int = (TH1D*) hintEffPA -> Clone("hrel_Eff_int");
  hrel_Acc_int->Scale(1./hintAccPP->GetBinContent(1));
  if(state!=3) hrel_Eff_int->Divide(hintEffPP);
  hRPAraw_int->Divide(hrel_Acc_int);
  if(state!=3) hRPAraw_int->Scale(hrel_Eff_int->GetBinContent(1));
  hRPAraw_int->Scale(Factor_Scale);

 //*****Draw rap***** 
  c_rap->cd(); 
  hrapSigPP->Draw();
//  gPad->SetLogy();
  hrapSigPA->Draw("same");

  hRPAraw_rap = (TH1D*)hrapSigPA->Clone("rpa_vs_rap_raw");
  hRPAraw_rap->Divide(hrapSigPP);

  TH1D* hrel_Acc_rap = (TH1D*) hrapAccPA -> Clone("hrel_Acc_rap");
  TH1D* hrel_Eff_rap = (TH1D*) hrapEffPA_pa -> Clone("hrel_Eff_rap");
  
  //Cros sec
  TH1D* hrap_cross_pA = (TH1D*)hrapSigPA_cross->Clone("rpa_vs_rap_cross");
  hrap_cross_pA->Divide(hrapAccPA_cross);
  cout << " :: " << hrapAccPA_cross->GetNbinsX() << endl;
  cout << " :: " << hrapEffPA_cross->GetNbinsX() << endl;
  hrap_cross_pA->Divide(hrapEffPA_cross);
  hrap_cross_pA->Scale(1./(1000.*lumi_pp));
  TH1ScaleByWidth(hrap_cross_pA);
  //  

  hrel_Acc_rap->Divide(hrapAccPP);
  hrel_Eff_rap->Divide(hrapEffPA_pp);
/*  hRPAraw_rap->Divide(hrel_Acc_rap);
  hRPAraw_rap->Divide(hrel_Eff_rap);
*/  hRPAraw_rap->Scale(Factor_Scale);


 //*****Draw pt***** 
  c_pt->cd(); 
  //hptSigPP->Draw();
  //gPad->SetLogy();
  //hptSigPA->Draw("same");

  hRPAraw_pt = (TH1D*)hptSigPA->Clone("raa_vs_pt_raw");
  hRPAraw_pt->Divide(hptSigPP);

  TH1D* hrel_Acc_pt = (TH1D*) hptAccPA -> Clone("hrel_Acc_pt");
  TH1D* hrel_Eff_pt = (TH1D*) hptAccPA -> Clone("hrel_Eff_pt");
  //Cros sec
  TH1D* hpt_cross_pA = (TH1D*)hptSigPA->Clone("rpa_vs_pt_cross");
  hpt_cross_pA->Divide(hrel_Acc_pt);
  hpt_cross_pA->Divide(hrel_Eff_pt);
  hpt_cross_pA->Scale(1./(1000.*lumi_pp*1.93*2));
  TH1ScaleByWidth(hpt_cross_pA);
  //
  hrel_Acc_pt->Divide(hptAccPP);
  hrel_Eff_pt->Divide(hptEffPP);
/*
  hRPAraw_pt->Divide(hrel_Acc_pt);
  hRPAraw_pt->Divide(hrel_Eff_pt);
*/  hRPAraw_pt->Scale(Factor_Scale);

  hrel_Eff_pt->Draw();
  
  //***** Save ******
  //***** int *****

  TCanvas* cIntRPA = new TCanvas("cRpA_int","",400,400);
  hRPA_int = (TH1D*)hRPAraw_int->Clone("rpa_vs_int");

  hRPA_int->SetMarkerStyle(24);
  hRPA_int->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_int->GetYaxis()->SetTitle("R_{pPb} corrected");
  cIntRPA->cd();
  hRPA_int->Draw("pe");
  jumSun(0,1,1,1,1,1);

  TGraphErrors *gRPA_int = new TGraphErrors(hRPA_int);
  gRPA_int->SetName("gRPA_int");
  gRPA_int->SetTitle("rpa_vs_int");


  //***** rap *****
  TCanvas* cRapRPA =  new TCanvas("cRpA_rap","",400,400);
  hRPA_rap = (TH1D*)hRPAraw_rap->Clone("rpa_vs_rap");

  hRPA_rap->SetMarkerStyle(24);
  hRPA_rap->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_rap->GetYaxis()->SetTitle("R_{pPb} corrected");
  cRapRPA->cd();
  hRPA_rap->Draw("pe");
  jumSun(yMin,1,yMax,1,1,1);

  TGraphErrors *gRPA_rap = new TGraphErrors(hRPAraw_rap);
  gRPA_rap->SetName("gRPA_rap");
  gRPA_rap->SetTitle("rpa_vs_rap");

  TGraphErrors *gCross_rap = new TGraphErrors(hrap_cross_pA);
  gCross_rap->SetName("gCross_rap");


  //******pt*******
  TCanvas* cPtRPA =  new TCanvas("cRpA_pt","",400,400);
  hRPA_pt = (TH1D*)hRPAraw_pt->Clone("rpa_vs_pt");

  hRPA_pt->SetMarkerStyle(24);
  hRPA_pt->GetYaxis()->SetRangeUser(0,1.5);
  hRPA_pt->GetYaxis()->SetTitle("R_{pPb} uncorrected");
  cPtRPA->cd();
  hRPA_pt->Draw("pe");
  jumSun(ptMin,1,ptMax,1,1,1);

  TGraphErrors *gRPA_pt = new TGraphErrors(hRPAraw_pt);
  gRPA_pt->SetName("gRPA_pt");
  gRPA_pt->SetTitle("rpa_vs_pt");
  
  TGraphErrors *gCross_pt = new TGraphErrors(hpt_cross_pA);
  gCross_pt->SetName("gCross_pt");

  TFile *wf = new TFile(Form("finalResults/Ups_%d_RPA_highPt.root",state),"recreate");
  gRPA_rap->Write();
  gRPA_pt->Write();
  gRPA_int->Write();
  gCross_rap->Write();
  gCross_pt->Write();
  wf->Close();
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh,   float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TFile* inf = new TFile(Form("NominalFitResult/jaredFit/NominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/afs/cern.ch/work/j/jaebeom/private/Analysis/RpA502TeV/NominalFitResult/jaredFit/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //Santona
  //TFile* inf = new TFile(Form("/home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/WorkingPlots/FitResults/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //Jaebeom
  //TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/samba/UpsilonAnalysis/fitResultFiles/mcFit_MuPt4_2016_11_04/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

void stripErrorBars( TH1* h, double defaultErr  ) {
  
  for ( int i=1;  i<= h->GetNbinsX() ; i++) {
    h->SetBinError( i, defaultErr);
  }
}

double getScale(int fTAA, double* TAA, double* centBin, int nCentBins)
{
  double flumi_;
  if(fTAA == nCentBins+1) flumi_ = 368;
  else if(centBin[fTAA-1]>=60 && centBin[fTAA-1]<120 && fTAA !=nCentBins+1) flumi_ = 464;
  else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 464;
  //else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 334.82249848;
  else if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1) flumi_ = 368;
  double nMBColl = NumberOfMBColl;
  if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1)  nMBColl = NumberOfMBColl;
  else if(centBin[fTAA-1]>=60 && fTAA !=nCentBins+1) nMBColl = NumberOfMBColl1;
  double scaleFactor;
  if(fTAA == nCentBins+1) scaleFactor = 28000000000000./(nMBColl*TAA[fTAA-1]*1000);
  else scaleFactor = 28000000000000./(nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200.*TAA[fTAA-1]*1000);
  
  if(fTAA!=nCentBins+1){
  cout << endl;
  cout << "icent : " << centBin[fTAA-1] << " - " << centBin[fTAA] << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "nMBColl*(centBin[fTAA]-centBin[fTAA-1])/100. : " << nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200. << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }

  else if(fTAA==nCentBins+1){
  cout << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }
  return scaleFactor;
}

