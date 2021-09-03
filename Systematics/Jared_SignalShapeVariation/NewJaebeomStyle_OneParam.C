#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void NewJaebeomStyle_OneParam(int whichUpsilon=1, int param=1, bool takeAvg=kFALSE) {

  TString avg = "";
  if (takeAvg) avg = "avg";

  TString whichParam;
  if (param==1) whichParam = "alpha";
  else if (param==2) whichParam = "f";
  else if (param==3) whichParam = "n";
  else if (param==4) whichParam = "x";
  
  TString filename = Form("FixParamErrors/Param%sDev%is.root",whichParam.Data(),whichUpsilon);
  TString filename2 = Form("FixParamErrors/Param%sDev%is_pt0to6to30.root",whichParam.Data(),whichUpsilon);
  TString filename3 = Form("FixParamErrors/Param%sDev%is_y193to000to193.root",whichParam.Data(),whichUpsilon);
  TString filename4 = Form("FixParamErrors/Param%sDev%is_hfbins.root",whichParam.Data(),whichUpsilon);
  TString filename5 = Form("FixParamErrors/Param%sDev%is_ntracksbins.root",whichParam.Data(),whichUpsilon);
  TString filename6 = Form("FixParamErrors/Param%sDev%isXSBins.root",whichParam.Data(),whichUpsilon);
  TString filename7 = Form("FixParamErrors/Param%sDev%is_pt0to6to30_in3Sbins.root",whichParam.Data(),whichUpsilon);
  //TString filename7 = Form("ErrorEstimates/SystematicErrorSignal%is_hfbins.root",whichUpsilon);
  //TString filename8 = Form("ErrorEstimates/SystematicErrorSignal%is_ntracksbins.root",whichUpsilon);
  //TString filename9 = Form("ErrorEstimates/SystematicErrorSignal%is_y287to193.root",whichUpsilon);
  //TString filename10 = Form("ErrorEstimates/SystematicErrorSignal%is_hfbinsIntBin.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/SysSig%is_%s.root",whichUpsilon,whichParam.Data());
  cout << outfilename << endl;

  TString ntupleptname = Form("ntuple%ispt%s;1",1,avg.Data());
  TString ntupleyname = Form("ntuple%isy%s;1",1,avg.Data());
  TString ntuplename = Form("ntuple%is%s;1",1,avg.Data());

  TFile *inFile = new TFile(filename);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);

  TFile *inFile2 = new TFile(filename2);
  TNtuple* ntuple1syptLowHigh = (TNtuple*)inFile2->Get(ntupleyname);

  TFile *inFile3 = new TFile(filename3);
  TNtuple* ntuple1sptyBackForw = (TNtuple*)inFile3->Get(ntupleptname);

  TFile *inFile4 = new TFile(filename4);
  TNtuple* ntuple1shf = (TNtuple*)inFile4->Get(ntuplename);

  TFile *inFile5 = new TFile(filename5);
  TNtuple* ntuple1sntracks = (TNtuple*)inFile5->Get(ntuplename);

  TFile *inFile6 = new TFile(filename6);
  TNtuple* ntuple1sptXSBins = (TNtuple*)inFile6->Get(ntupleptname);

  TFile *inFile7 = new TFile(filename7);
  TNtuple* ntuple1syptLowHigh_in3Sbins = (TNtuple*)inFile7->Get(ntupleyname);

  //TFile *inFile7 = new TFile(filename7);
  //TNtuple* ntuplehf0to193 = (TNtuple*)inFile7->Get("ntuple1s;1");

  //TFile *inFile8 = new TFile(filename8);
  //TNtuple* ntuplent0to193 = (TNtuple*)inFile8->Get("ntuple1s;1");

  //TFile *inFile9 = new TFile(filename9);
  //TNtuple* ntuple1sy287to193 = (TNtuple*)inFile9->Get(ntupleyname);

  //TFile *inFile10 = new TFile(filename10);
  //TNtuple* ntupleFullHF = (TNtuple*)inFile10->Get("ntuple;1");

  float binlowlim, binuplim;
  float binlowlimy, binuplimy;
  double pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float y287bins1[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float hfbins1[5] = {0,12,19,27,120};
  float ntbins1[5] = {0,40,62,88,400};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float y287bins2[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  float hfbins2[5] = {0,12,19,27,120};
  float ntbins2[5] = {0,40,62,88,400};
  float ptbins3[3] = {0,6,30};
  float ybins3[3] = {-1.93,0.0,1.93};
  float y287bins3[4] = {-2.87,-1.93,0.0,1.93};
  float hfbins3[3] = {0,12,120};
  float ntbins3[3] = {0,40,400};

  float ybins3s[3] = {-1.93,0.0,1.93};

  const int numptbins1 = sizeof(ptbins1)/sizeof(float)-1;
  const int numybins1 = sizeof(ybins1)/sizeof(float)-1;
  const int numy287bins1 = sizeof(y287bins1)/sizeof(float)-1;
  const int numhfbins1 = sizeof(hfbins1)/sizeof(float)-1;
  const int numntbins1 = sizeof(ntbins1)/sizeof(float)-1;
  const int numbins1 = (numy287bins1-2)/2;

  const int numptbins2 = sizeof(ptbins2)/sizeof(float)-1;
  const int numybins2 = sizeof(ybins2)/sizeof(float)-1;
  const int numy287bins2 = sizeof(y287bins2)/sizeof(float)-1;
  const int numhfbins2 = sizeof(hfbins2)/sizeof(float)-1;
  const int numntbins2 = sizeof(ntbins2)/sizeof(float)-1;
  const int numbins2 = (numy287bins2-2)/2;

  const int numptbins3 = sizeof(ptbins3)/sizeof(float)-1;
  const int numybins3 = sizeof(ybins3)/sizeof(float)-1;
  const int numy287bins3 = sizeof(y287bins3)/sizeof(float)-1;
  const int numhfbins3 = sizeof(hfbins3)/sizeof(float)-1;
  const int numntbins3 = sizeof(ntbins3)/sizeof(float)-1;
  const int numbins3 = (numy287bins3-2)/2;

  float* ptbinsptr;
  float* ybinsptr;
  float* y287binsptr;
  float* hfbinsptr;
  float* ntbinsptr;

  int numptbinstmp, numybinstmp, numy287binstmp, numhfbinstmp, numntbinstmp, numbinstmp;

  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    ybinsptr = &ybins1[0];
    y287binsptr = &y287bins1[0];
    hfbinsptr = &hfbins1[0];
    ntbinsptr = &ntbins1[0];
    numptbinstmp = numptbins1;
    numybinstmp = numybins1;
    numy287binstmp = numy287bins1;
    numhfbinstmp = numhfbins1;
    numntbinstmp = numntbins1;
    numbinstmp = numbins1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    ybinsptr = &ybins2[0];
    y287binsptr = &y287bins2[0];
    hfbinsptr = &hfbins2[0];
    ntbinsptr = &ntbins2[0];
    numptbinstmp = numptbins2;
    numybinstmp = numybins2;
    numy287binstmp = numy287bins2;
    numhfbinstmp = numhfbins2;
    numntbinstmp = numntbins2;
    numbinstmp = numbins2;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    ybinsptr = &ybins3[0];
    y287binsptr = &y287bins3[0];
    hfbinsptr = &hfbins3[0];
    ntbinsptr = &ntbins3[0];
    numptbinstmp = numptbins3;
    numybinstmp = numybins3;
    numy287binstmp = numy287bins3;
    numhfbinstmp = numhfbins3;
    numntbinstmp = numntbins3;
    numbinstmp = numbins3;
  }

  const int numptbins = numptbinstmp;
  const int numybins = numybinstmp;
  const int numy287bins = numy287binstmp;
  const int numhfbins = numhfbinstmp;
  const int numntbins = numntbinstmp;
  const int numbins = numbinstmp;

  //Declare histograms
  //Integrated bin
  TH1D* hSignalErrIntPP = new TH1D("hintSysSigPP","Systematic for PP in integrated bin",1,0,30);
  TH1D* hSignalErrIntPA = new TH1D("hintSysSigPA","Systematic for PA in integrated bin",1,0,30);
  TH1D* hSignalErrIntRpA = new TH1D("hintSysSigRpA","Systematic for RpA in integrated bin",1,0,30);

  //Normal pt and y bins
  TH1D* hSignalErrptPP = new TH1D("hptSysSigPP","Systematic for PP in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPP = new TH1D("hySysSigPP","Systematic for PP in y",numybins,ybinsptr);
  TH1D* hSignalErrptPA = new TH1D("hptSysSigPA","Systematic for PA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPA = new TH1D("hySysSigPA","Systematic for PA in y",numybins,ybinsptr);
  TH1D* hSignalErrptRpA = new TH1D("hptSysSigRpA","Systematic for RpA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryRpA = new TH1D("hySysSigRpA","Systematic for RpA in y",numybins,ybinsptr);

  //Differential bins
  TH1D* hSignalErrptPPBackwardY = new TH1D("hptSysSigPPBackwardY","Systematic for PP in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPPLowPt = new TH1D("hySysSigPPLowPt","Systematic for PP in y",numybins,ybinsptr);
  TH1D* hSignalErrptPABackwardY = new TH1D("hptSysSigPABackwardY","Systematic for PA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPALowPt = new TH1D("hySysSigPALowPt","Systematic for PA in y",numybins,ybinsptr);
  TH1D* hSignalErrptRpABackwardY = new TH1D("hptSysSigRpABackwardY","Systematic for RpA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryRpALowPt = new TH1D("hySysSigRpALowPt","Systematic for RpA in y",numybins,ybinsptr);

  TH1D* hSignalErrptPPForwardY = new TH1D("hptSysSigPPForwardY","Systematic for PP in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPPHighPt = new TH1D("hySysSigPPHighPt","Systematic for PP in y",numybins,ybinsptr);
  TH1D* hSignalErrptPAForwardY = new TH1D("hptSysSigPAForwardY","Systematic for PA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryPAHighPt = new TH1D("hySysSigPAHighPt","Systematic for PA in y",numybins,ybinsptr);
  TH1D* hSignalErrptRpAForwardY = new TH1D("hptSysSigRpAForwardY","Systematic for RpA in pt",numptbins,ptbinsptr);
  TH1D* hSignalErryRpAHighPt = new TH1D("hySysSigRpAHighPt","Systematic for RpA in y",numybins,ybinsptr);

  TH1D* hSignalErryPPLowPt_in3Sbins = new TH1D("hySysSigPPLowPt_in3Sbins","Systematic for PP in y",2,ybins3s);
  TH1D* hSignalErryPALowPt_in3Sbins = new TH1D("hySysSigPALowPt_in3Sbins","Systematic for PA in y",2,ybins3s);
  TH1D* hSignalErryRpALowPt_in3Sbins = new TH1D("hySysSigRpALowPt_in3Sbins","Systematic for RpA in y",2,ybins3s);
  TH1D* hSignalErryPPHighPt_in3Sbins = new TH1D("hySysSigPPHighPt_in3Sbins","Systematic for PP in y",2,ybins3s);
  TH1D* hSignalErryPAHighPt_in3Sbins = new TH1D("hySysSigPAHighPt_in3Sbins","Systematic for PA in y",2,ybins3s);
  TH1D* hSignalErryRpAHighPt_in3Sbins = new TH1D("hySysSigRpAHighPt_in3Sbins","Systematic for RpA in y",2,ybins3s);

  //Cross section histograms
  //TH1D* hrapSysSigCrossLowPt = new TH1D("hrapSysSigCrossLowPt","Systematic for PA in y",numy287bins,y287binsptr);
  //TH1D* hrapSysSigCrossHighPt = new TH1D("hrapSysSigCrossHighPt","Systematic for PA in pt for y in [-2.87,1.93]",numy287bins,y287binsptr);
  TH1D* hSignalErrptPAy287to193 = new TH1D("hptSysSigPA_y287to193","Systematic for PA in pt",numptbins,ptbinsptr);
  TH1D* hrapSysSigCross = new TH1D("hrapSysSigCross","Systematic for PA in y including backward bin [-2.87,-1.93]",numy287bins,y287binsptr);

  //RFB histograms
  //TH1D* hHFSysSigRFB000to040 = new TH1D("hHFSysSigRFB000to040","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hHFSysSigRFB040to080 = new TH1D("hHFSysSigRFB040to080","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hHFSysSigRFB080to120 = new TH1D("hHFSysSigRFB080to120","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hHFSysSigRFB120to193 = new TH1D("hHFSysSigRFB120to193","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hHFSysSigRFB000to080 = new TH1D("hHFSysSigRFB000to080","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hHFSysSigRFB080to193 = new TH1D("hHFSysSigRFB080to193","Systematic for PA in hf",numhfbins,hfbinsptr);
  TH1D* hHFSysSigRFB000to193 = new TH1D("hHFSysSigRFB000to193","Systematic for PA in hf",numhfbins,hfbinsptr);
  //TH1D* hNtracksSysSigRFB = new TH1D("hNtracksSysSigRFB","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB000to040 = new TH1D("hNtracksSysSigRFB000to040","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB040to080 = new TH1D("hNtracksSysSigRFB040to080","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB080to120 = new TH1D("hNtracksSysSigRFB080to120","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB120to193 = new TH1D("hNtracksSysSigRFB120to193","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB000to080 = new TH1D("hNtracksSysSigRFB000to080","Systematic for PA in ntracks",numntbins,ntbinsptr);
  //TH1D* hNtracksSysSigRFB080to193 = new TH1D("hNtracksSysSigRFB080to193","Systematic for PA in ntracks",numntbins,ntbinsptr);
  TH1D* hNtracksSysSigRFB000to193 = new TH1D("hNtracksSysSigRFB000to193","Systematic for PA in ntracks",numntbins,ntbinsptr);
  TH1D* hSysSigRFBIntActivity = new TH1D("hSysSigRFBIntActivity","Systematic for PA in integrated activity",1,0,120);

  //Get errors in integrated bin
  cout << "--------------------------------" << endl;
  cout << "Integrated bin" << endl;
  TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf(strpp1sErr);
  TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf(strpPb1sErr);
  TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf(strRpPbErr);
  TLeaf *binlowleaf = ntuple1spt->GetLeaf("binlowlim");
  TLeaf *binupleaf = ntuple1spt->GetLeaf("binuplim");
    ntuple1spt->GetEntry(0);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("pt%.1f-%.1f:\t",binlowlim,binuplim);
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << endl;
    hSignalErrIntPP->SetBinContent(1, TMath::Abs(pp1sErr)/100);
    hSignalErrIntPA->SetBinContent(1, TMath::Abs(pPb1sErr)/100);
    hSignalErrIntRpA->SetBinContent(1, TMath::Abs(RpPbErr)/100);

  //Get errors in pt bins
  cout << "--------------------------------" << endl;
  cout << "pt bins" << endl;
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    ntuple1spt->GetEntry(ipt);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("pt%.1f-%.1f:\t",binlowlim,binuplim);
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << endl;
    hSignalErrptPP->SetBinContent(ipt, TMath::Abs(pp1sErr)/100);
    hSignalErrptPA->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
    hSignalErrptRpA->SetBinContent(ipt, TMath::Abs(RpPbErr)/100);
  }

  //Get errors in pt bins with y in [-2.87,1.93]
  pPb1sErrLeaf = ntuple1sptXSBins->GetLeaf(strpPb1sErr);
  binlowleaf = ntuple1sptXSBins->GetLeaf("binlowlim");
  binupleaf = ntuple1sptXSBins->GetLeaf("binuplim");
  cout << "--------------------------------" << endl;
  cout << "pt bins with y in [-2.87,1.93]" << endl;
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    ntuple1sptXSBins->GetEntry(ipt);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("pt%.1f-%.1f:\t",binlowlim,binuplim);
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    cout << Form("%.4f",pPb1sErr) << endl;
    hSignalErrptPAy287to193->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
  }

  //Get errors in y bins
  cout << "--------------------------------" << endl;
  cout << "y bins" << endl;
  pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
  pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
  RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
  binlowleaf = ntuple1sy->GetLeaf("binlowlim");
  binupleaf = ntuple1sy->GetLeaf("binuplim");
  for (int iy = 1; iy<numy287bins; iy++) {
    ntuple1sy->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f:\t",binlowlim,binuplim);
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << endl;
    hSignalErryPA->SetBinContent(iy, TMath::Abs(pPb1sErr)/100);
    if (iy>0) {
    hSignalErryPP->SetBinContent(iy, TMath::Abs(pp1sErr)/100);
    hSignalErryRpA->SetBinContent(iy, TMath::Abs(RpPbErr)/100);
    }
  }

  //Get errors in y bins for cross-section
  cout << "--------------------------------" << endl;
  cout << "y bins for cross-section" << endl;
  pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
  pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
  RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
  binlowleaf = ntuple1sy->GetLeaf("binlowlim");
  binupleaf = ntuple1sy->GetLeaf("binuplim");
  for (int iy = 0; iy<numy287bins; iy++) {
    ntuple1sy->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f:\t",binlowlim,binuplim);
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    cout << Form("%.4f",pPb1sErr) << endl;
    hrapSysSigCross->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
  }
/*
  //Get errors in pt bins in forward and backward y
  cout << "--------------------------------" << endl;
  cout << "pt bins in backward and forward y" << endl;
  pp1sErrLeaf = ntuple1sptyBackForw->GetLeaf(strpp1sErr);
  pPb1sErrLeaf = ntuple1sptyBackForw->GetLeaf(strpPb1sErr);
  RpPbErrLeaf = ntuple1sptyBackForw->GetLeaf(strRpPbErr);
  binlowleaf = ntuple1sptyBackForw->GetLeaf("binlowlim");
  binupleaf = ntuple1sptyBackForw->GetLeaf("binuplim");
  for (int ipt = 1; ipt<numptbins*2+1; ipt++) {
    ntuple1sptyBackForw->GetEntry(ipt-1);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("pt%.1f-%.1f:\t",binlowlim,binuplim);
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << "\t";
    if (ipt<numptbins+1) {
      cout << "Backward" << endl;
      hSignalErrptPPBackwardY->SetBinContent(ipt, TMath::Abs(pp1sErr)/100);
      hSignalErrptPABackwardY->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
      hSignalErrptRpABackwardY->SetBinContent(ipt, TMath::Abs(RpPbErr)/100);
    }
    else {
      cout << "Forward" << endl;
      hSignalErrptPPForwardY->SetBinContent(ipt-numptbins, TMath::Abs(pp1sErr)/100);
      hSignalErrptPAForwardY->SetBinContent(ipt-numptbins, TMath::Abs(pPb1sErr)/100);
      hSignalErrptRpAForwardY->SetBinContent(ipt-numptbins, TMath::Abs(RpPbErr)/100);
    }
  }

  //Get errors in y bins in low and high pt
  cout << "--------------------------------" << endl;
  cout << "y bins in low and high pt" << endl;
  pp1sErrLeaf = ntuple1syptLowHigh->GetLeaf(strpp1sErr);
  pPb1sErrLeaf = ntuple1syptLowHigh->GetLeaf(strpPb1sErr);
  RpPbErrLeaf = ntuple1syptLowHigh->GetLeaf(strRpPbErr);
  binlowleaf = ntuple1syptLowHigh->GetLeaf("binlowlim");
  binupleaf = ntuple1syptLowHigh->GetLeaf("binuplim");
  for (int iy = 0; iy<numybins*2; iy++) {
    //if (iy==numy287bins) iy++;
    ntuple1syptLowHigh->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f:\t",binlowlim,binuplim);
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << "\t";
    if (iy<numybins) {
      cout << " Low Pt " << iy+1 << endl;
      hSignalErryPPLowPt->SetBinContent(iy+1, TMath::Abs(pp1sErr)/100);
      hSignalErryPALowPt->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpALowPt->SetBinContent(iy+1, TMath::Abs(RpPbErr)/100);
    }
    else {
      cout << " High Pt " << iy+1-numybins << endl;
      hSignalErryPPHighPt->SetBinContent(iy+1-numybins, TMath::Abs(pp1sErr)/100);
      hSignalErryPAHighPt->SetBinContent(iy+1-numybins, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpAHighPt->SetBinContent(iy+1-numybins, TMath::Abs(RpPbErr)/100);
    }
  }*/

  //Get errors in y bins in low and high pt in 3S y bins
  cout << "--------------------------------" << endl;
  cout << "y bins in low and high pt in 3S bins" << endl;
  pp1sErrLeaf = ntuple1syptLowHigh_in3Sbins->GetLeaf(strpp1sErr);
  pPb1sErrLeaf = ntuple1syptLowHigh_in3Sbins->GetLeaf(strpPb1sErr);
  RpPbErrLeaf = ntuple1syptLowHigh_in3Sbins->GetLeaf(strRpPbErr);
  binlowleaf = ntuple1syptLowHigh_in3Sbins->GetLeaf("binlowlim");
  binupleaf = ntuple1syptLowHigh_in3Sbins->GetLeaf("binuplim");
  for (int iy = 0; iy<2*2; iy++) {
    //if (iy==numy287bins) iy++;
    ntuple1syptLowHigh_in3Sbins->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f:\t",binlowlim,binuplim);
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",pp1sErr) << "\t";
    cout << Form("%.4f",pPb1sErr) << "\t";
    cout << Form("%.4f",RpPbErr) << "\t";
    if (iy<2) {
      cout << " Low Pt " << iy+1 << endl;
      hSignalErryPPLowPt_in3Sbins->SetBinContent(iy+1, TMath::Abs(pp1sErr)/100);
      hSignalErryPALowPt_in3Sbins->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpALowPt_in3Sbins->SetBinContent(iy+1, TMath::Abs(RpPbErr)/100);
    }
    else {
      cout << " High Pt " << iy+1-2 << endl;
      hSignalErryPPHighPt_in3Sbins->SetBinContent(iy+1-2, TMath::Abs(pp1sErr)/100);
      hSignalErryPAHighPt_in3Sbins->SetBinContent(iy+1-2, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpAHighPt_in3Sbins->SetBinContent(iy+1-2, TMath::Abs(RpPbErr)/100);
    }
  }

  //Get errors in hf bins
  cout << "--------------------------------" << endl;
  cout << "HF bins" << endl;
  RpPbErrLeaf = ntuple1shf->GetLeaf("RFBErr");
  binlowleaf = ntuple1shf->GetLeaf("binlowlim");
  binupleaf = ntuple1shf->GetLeaf("binuplim");
  TLeaf* binlowleafy = ntuple1shf->GetLeaf("binlowlimy");
  TLeaf* binupleafy = ntuple1shf->GetLeaf("binuplimy");
  for (int ihf = 1; ihf<numhfbins+1; ihf++) {
    ntuple1shf->GetEntry(ihf);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("hf%.0f-%.0f_y%.2f-%.2f:\t",binlowlim,binuplim,binlowlimy,binuplimy);
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",RpPbErr) << endl;
    hHFSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
  }


  //Get errors in ntracks bins
  cout << "--------------------------------" << endl;
  cout << "ntracks bins" << endl;
  RpPbErrLeaf = ntuple1sntracks->GetLeaf("RFBErr");
  binlowleaf = ntuple1sntracks->GetLeaf("binlowlim");
  binupleaf = ntuple1sntracks->GetLeaf("binuplim");
  binlowleafy = ntuple1sntracks->GetLeaf("binlowlimy");
  binupleafy = ntuple1sntracks->GetLeaf("binuplimy");
  for (int ihf = 1; ihf<numntbins+1; ihf++) {
    ntuple1sntracks->GetEntry(ihf);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("ntracks%.0f-%.0f_y%.2f-%.2f: \t",binlowlim,binuplim,binlowlimy,binuplimy);
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",RpPbErr) << endl;
    hNtracksSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
  }

  //Get errors in integrated activity bin
  cout << "--------------------------------" << endl;
  cout << "Integrated activity bin" << endl;
    ntuple1sntracks->GetEntry(0);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("hfsum%.2f-%.2f_y%.2f-%.2f:\t",binlowlim,binuplim,binlowlimy,binuplimy);
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << Form("%.4f",RpPbErr) << endl;
    hSysSigRFBIntActivity->SetBinContent(1, TMath::Abs(RpPbErr)/100);


  //save histograms
  TFile outFile(outfilename, "RECREATE");

  hSignalErrIntPP->Write();
  hSignalErrIntPA->Write();
  hSignalErrIntRpA->Write();

  hSignalErrptPP->Write();
  hSignalErrptPA->Write();
  hSignalErrptRpA->Write();
  hSignalErryPP->Write();
  hSignalErryPA->Write();
  hSignalErryRpA->Write();
/*
  hSignalErrptPPBackwardY->Write();
  hSignalErrptPABackwardY->Write();
  hSignalErrptRpABackwardY->Write();
  hSignalErrptPPForwardY->Write();
  hSignalErrptPAForwardY->Write();
  hSignalErrptRpAForwardY->Write();

  hSignalErryPPLowPt->Write();
  hSignalErryPALowPt->Write();
  hSignalErryRpALowPt->Write();
  hSignalErryPPHighPt->Write();
  hSignalErryPAHighPt->Write();
  hSignalErryRpAHighPt->Write();
*/
  hSignalErryPPLowPt_in3Sbins->Write();
  hSignalErryPALowPt_in3Sbins->Write();
  hSignalErryRpALowPt_in3Sbins->Write();
  hSignalErryPPHighPt_in3Sbins->Write();
  hSignalErryPAHighPt_in3Sbins->Write();
  hSignalErryRpAHighPt_in3Sbins->Write();

  //hrapSysSigCrossLowPt->Write();
  //hrapSysSigCrossHighPt->Write();
  hSignalErrptPAy287to193->Write();
  hrapSysSigCross->Write();
/*
  if (whichUpsilon==1) {
    hHFSysSigRFB000to040->Write();
    hHFSysSigRFB040to080->Write();
    hHFSysSigRFB080to120->Write();
    hHFSysSigRFB120to193->Write();
    hNtracksSysSigRFB000to040->Write();
    hNtracksSysSigRFB040to080->Write();
    hNtracksSysSigRFB080to120->Write();
    hNtracksSysSigRFB120to193->Write();
  }
  else if (whichUpsilon==2) {
    hHFSysSigRFB000to080->Write();
    hHFSysSigRFB080to193->Write();
    hNtracksSysSigRFB000to080->Write();
    hNtracksSysSigRFB080to193->Write();
  }*/
  hHFSysSigRFB000to193->Write();
  hNtracksSysSigRFB000to193->Write();
  hSysSigRFBIntActivity->Write();

  outFile.Close();
}
