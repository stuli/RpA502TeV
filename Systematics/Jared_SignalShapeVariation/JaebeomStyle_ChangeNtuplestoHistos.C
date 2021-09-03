#include <iostream>
#include <iomanip>
#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TH1.h"

using namespace std;
void JaebeomStyle_ChangeNtuplestoHistos(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is.root",whichUpsilon);
  TString filename2 = Form("ErrorEstimates/SystematicErrorSignal%is_pt0to6to30.root",whichUpsilon);
  TString filename3 = Form("ErrorEstimates/SystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  TString filename4 = Form("ErrorEstimates/SystematicErrorSignal%is_hfbins.root",whichUpsilon);
  TString filename5 = Form("ErrorEstimates/SystematicErrorSignal%is_ntracksbins.root",whichUpsilon);
  TString filename6 = Form("ErrorEstimates/SystematicErrorSignal%isXSBins.root",whichUpsilon);
  TString filename7 = Form("ErrorEstimates/SystematicErrorSignal12s_hfbins.root");
  TString filename8 = Form("ErrorEstimates/SystematicErrorSignal12s_ntracksbins.root");
  TString filename9 = Form("ErrorEstimates/SystematicErrorSignal%is_y287to193.root",whichUpsilon);
  TString filename10 = Form("ErrorEstimates/SystematicErrorSignal%is_hfbinsIntBin.root",whichUpsilon);
  TString outfilename = Form("ErrorEstimates/SysSig%is.root",whichUpsilon);
  cout << outfilename << endl;

  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TString ntuplename = Form("ntuple%is;1",whichUpsilon);

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
  TNtuple* ntuplehf0to193 = (TNtuple*)inFile7->Get("ntuple12;1");

  TFile *inFile8 = new TFile(filename8);
  TNtuple* ntuplent0to193 = (TNtuple*)inFile8->Get("ntuple12;1");

  TFile *inFile9 = new TFile(filename9);
  TNtuple* ntuple1sy287to193 = (TNtuple*)inFile9->Get(ntupleyname);

  TFile *inFile10 = new TFile(filename10);
  TNtuple* ntupleFullHF = (TNtuple*)inFile10->Get("ntuple;1");

  float binlowlim, binuplim;
  float binlowlimy, binuplimy;
  double pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[10] = {-2.4,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float y287bins[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
    float hfbins[5] = {0,12,19,27,120};
    float ntbins[5] = {0,40,62,88,400};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[6] = {-2.4,-1.93,-0.8,0.0,0.8,1.93};
    float y287bins[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
    float hfbins[5] = {0,12,19,27,120};
    float ntbins[5] = {0,40,62,88,400};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[4] = {-2.4,-1.93,0.0,1.93};
    float y287bins[4] = {-2.87,-1.93,0.0,1.93};
    float hfbins[3] = {0,12,120};
    float ntbins[3] = {0,40,400};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numy287bins = sizeof(y287bins)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntbins = sizeof(ntbins)/sizeof(float)-1;
  const int numbins = (numy287bins-2)/2;

  //Declare histograms
  //Integrated bin
  TH1D* hSignalErrIntPP = new TH1D("hintSysSigPP","Systematic for PP in integrated bin",1,0,30);
  TH1D* hSignalErrIntPA = new TH1D("hintSysSigPA","Systematic for PA in integrated bin",1,0,30);
  TH1D* hSignalErrIntRpA = new TH1D("hintSysSigRpA","Systematic for RpA in integrated bin",1,0,30);

  //Normal pt and y bins
  TH1D* hSignalErrptPP = new TH1D("hptSysSigPP","Systematic for PP in pt",numptbins,ptbins);
  TH1D* hSignalErryPP = new TH1D("hySysSigPP","Systematic for PP in y",numybins,ybins);
  TH1D* hSignalErrptPA = new TH1D("hptSysSigPA","Systematic for PA in pt",numptbins,ptbins);
  TH1D* hSignalErryPA = new TH1D("hySysSigPA","Systematic for PA in y",numybins,ybins);
  TH1D* hSignalErrptRpA = new TH1D("hptSysSigRpA","Systematic for RpA in pt",numptbins,ptbins);
  TH1D* hSignalErryRpA = new TH1D("hySysSigRpA","Systematic for RpA in y",numybins,ybins);

  //Differential bins
  TH1D* hSignalErrptPPBackwardY = new TH1D("hptSysSigPPBackwardY","Systematic for PP in pt",numptbins,ptbins);
  TH1D* hSignalErryPPLowPt = new TH1D("hySysSigPPLowPt","Systematic for PP in y",numybins,ybins);
  TH1D* hSignalErrptPABackwardY = new TH1D("hptSysSigPABackwardY","Systematic for PA in pt",numptbins,ptbins);
  TH1D* hSignalErryPALowPt = new TH1D("hySysSigPALowPt","Systematic for PA in y",numybins,ybins);
  TH1D* hSignalErrptRpABackwardY = new TH1D("hptSysSigRpABackwardY","Systematic for RpA in pt",numptbins,ptbins);
  TH1D* hSignalErryRpALowPt = new TH1D("hySysSigRpALowPt","Systematic for RpA in y",numybins,ybins);

  TH1D* hSignalErrptPPForwardY = new TH1D("hptSysSigPPForwardY","Systematic for PP in pt",numptbins,ptbins);
  TH1D* hSignalErryPPHighPt = new TH1D("hySysSigPPHighPt","Systematic for PP in y",numybins,ybins);
  TH1D* hSignalErrptPAForwardY = new TH1D("hptSysSigPAForwardY","Systematic for PA in pt",numptbins,ptbins);
  TH1D* hSignalErryPAHighPt = new TH1D("hySysSigPAHighPt","Systematic for PA in y",numybins,ybins);
  TH1D* hSignalErrptRpAForwardY = new TH1D("hptSysSigRpAForwardY","Systematic for RpA in pt",numptbins,ptbins);
  TH1D* hSignalErryRpAHighPt = new TH1D("hySysSigRpAHighPt","Systematic for RpA in y",numybins,ybins);

  //Cross section histograms
  TH1D* hrapSysSigCrossLowPt = new TH1D("hrapSysSigCrossLowPt","Systematic for PA in y",numy287bins,y287bins);
  TH1D* hrapSysSigCrossHighPt = new TH1D("hrapSysSigCrossHighPt","Systematic for PA in pt for y in [-2.87,1.93]",numy287bins,y287bins);
  TH1D* hSignalErrptPAy287to193 = new TH1D("hptSysSigPA_y287to193","Systematic for PA in pt",numptbins,ptbins);
  TH1D* hrapSysSigCross = new TH1D("hrapSysSigCross","Systematic for PA in y including backward bin [-2.87,-1.93]",numy287bins,y287bins);

  //RFB histograms
  TH1D* hHFSysSigRFB000to040 = new TH1D("hHFSysSigRFB000to040","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB040to080 = new TH1D("hHFSysSigRFB040to080","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB080to120 = new TH1D("hHFSysSigRFB080to120","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB120to193 = new TH1D("hHFSysSigRFB120to193","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB000to080 = new TH1D("hHFSysSigRFB000to080","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB080to193 = new TH1D("hHFSysSigRFB080to193","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hHFSysSigRFB000to193 = new TH1D("hHFSysSigRFB000to193","Systematic for PA in hf",numhfbins,hfbins);
  TH1D* hNtracksSysSigRFB = new TH1D("hNtracksSysSigRFB","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB000to040 = new TH1D("hNtracksSysSigRFB000to040","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB040to080 = new TH1D("hNtracksSysSigRFB040to080","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB080to120 = new TH1D("hNtracksSysSigRFB080to120","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB120to193 = new TH1D("hNtracksSysSigRFB120to193","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB000to080 = new TH1D("hNtracksSysSigRFB000to080","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB080to193 = new TH1D("hNtracksSysSigRFB080to193","Systematic for PA in ntracks",numntbins,ntbins);
  TH1D* hNtracksSysSigRFB000to193 = new TH1D("hNtracksSysSigRFB000to193","Systematic for PA in ntracks",numntbins,ntbins);
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
    cout << Form("pt%.1f-%.1f",binlowlim,binuplim) << endl;
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
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
    cout << Form("pt%.1f-%.1f",binlowlim,binuplim) << endl;
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
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
    cout << Form("pt%.1f-%.1f",binlowlim,binuplim) << endl;
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
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
  for (int iy = 0; iy<numy287bins+1; iy++) {
    ntuple1sy->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f",binlowlim,binuplim) << endl;
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    cout << "pp1sErr = " << pp1sErr << endl;
    hSignalErryPA->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
    if (iy>0) {
    hSignalErryPP->SetBinContent(iy, TMath::Abs(pp1sErr)/100);
    hSignalErryRpA->SetBinContent(iy, TMath::Abs(RpPbErr)/100);
    }
  }

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
    cout << Form("pt%.1f-%.1f",binlowlim,binuplim) << endl;
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    if (ipt<numptbins+1) {
      hSignalErrptPPBackwardY->SetBinContent(ipt, TMath::Abs(pp1sErr)/100);
      hSignalErrptPABackwardY->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
      hSignalErrptRpABackwardY->SetBinContent(ipt, TMath::Abs(RpPbErr)/100);
    }
    else {
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
  for (int iy = 1; iy<numybins*2+2; iy++) {
    if (iy==numy287bins) iy++;
    ntuple1syptLowHigh->GetEntry(iy);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f",binlowlim,binuplim) << endl;
    pp1sErr = (double)pp1sErrLeaf->GetValue();
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    if (iy<numybins+1) {
      hSignalErryPPLowPt->SetBinContent(iy, TMath::Abs(pp1sErr)/100);
      hSignalErryPALowPt->SetBinContent(iy, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpALowPt->SetBinContent(iy, TMath::Abs(RpPbErr)/100);
    }
    else {
      hSignalErryPPHighPt->SetBinContent(iy-numy287bins, TMath::Abs(pp1sErr)/100);
      hSignalErryPAHighPt->SetBinContent(iy-numy287bins, TMath::Abs(pPb1sErr)/100);
      hSignalErryRpAHighPt->SetBinContent(iy-numy287bins, TMath::Abs(RpPbErr)/100);
    }
  }

  //Get errors in y bins including -2.87--1.93 for cross section
  pPb1sErrLeaf = ntuple1sy287to193->GetLeaf(strpPb1sErr);
  binlowleaf = ntuple1sy287to193->GetLeaf("binlowlim");
  binupleaf = ntuple1sy287to193->GetLeaf("binuplim");
  cout << "--------------------------------" << endl;
  cout << "y bins in including [-2.87,-1.93]" << endl;
  for (int iy = 1; iy<numy287bins+1; iy++) {
    ntuple1sy287to193->GetEntry(iy-1);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f",binlowlim,binuplim) << endl;
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    cout << "Error = " << pPb1sErr << endl;
    hrapSysSigCross->SetBinContent(iy, TMath::Abs(pPb1sErr)/100);
  }

  //Get errors in y bins in low and high pt for cross section
  cout << "--------------------------------" << endl;
  cout << "y bins in low and high pt for cross-section" << endl;
  for (int iy = 1; iy<numy287bins*2+1; iy++) {
    ntuple1syptLowHigh->GetEntry(iy-1);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    cout << Form("y%.2f-%.2f",binlowlim,binuplim);
    pPb1sErr = (double)pPb1sErrLeaf->GetValue();
    if (iy<numy287bins+1) {
      cout << " Low Pt" << endl;
      hrapSysSigCrossLowPt->SetBinContent(iy, TMath::Abs(pPb1sErr)/100);
    }
    else {
      cout << " High Pt" << endl;
      hrapSysSigCrossHighPt->SetBinContent(iy-numy287bins, TMath::Abs(pPb1sErr)/100);
    }
  }

  //Get errors in hf bins
  cout << "--------------------------------" << endl;
  cout << "HF bins" << endl;
/*  RpPbErrLeaf = ntuple1shf->GetLeaf("RFBErr");
  binlowleaf = ntuple1shf->GetLeaf("binlowlim");
  binupleaf = ntuple1shf->GetLeaf("binuplim");
  TLeaf* binlowleafy = ntuple1shf->GetLeaf("binlowlimy");
  TLeaf* binupleafy = ntuple1shf->GetLeaf("binuplimy");
  for (int whichBin=1; whichBin<=numbins; whichBin++) {
  for (int ihf = 1; ihf<numhfbins+1; ihf++) {
    if (whichUpsilon==1) {
      ntuple1shf->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-5+whichBin);
    }
    if (whichUpsilon==2) {
      ntuple1shf->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-3+whichBin);
    }
    if (whichUpsilon==3) {
      ntuple1shf->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-1);
    }
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("hf%.0f-%.0f_y%.2f-%.2f",binlowlim,binuplim,binlowlimy,binuplimy) << endl;
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    if (whichBin==1) {
      if (whichUpsilon==1) hHFSysSigRFB000to040->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==2) hHFSysSigRFB000to080->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==3) hHFSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==2) {
      if (whichUpsilon==1) hHFSysSigRFB040to080->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==2) hHFSysSigRFB080to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==3) {
      hHFSysSigRFB080to120->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==4) {
      hHFSysSigRFB120to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
  }
  }*/

  if (whichUpsilon<3) {
  //Get errors in hf bins with |y|<1.93 for 1s and 2s
  RpPbErrLeaf = ntuplehf0to193->GetLeaf("RFBErr");
  binlowleaf = ntuplehf0to193->GetLeaf("binlowlim");
  binupleaf = ntuplehf0to193->GetLeaf("binuplim");
  binlowleafy = ntuplehf0to193->GetLeaf("binlowlimy");
  binupleafy = ntuplehf0to193->GetLeaf("binuplimy");
  }
  else {
  RpPbErrLeaf = ntuple1shf->GetLeaf("RFBErr");
  binlowleaf = ntuple1shf->GetLeaf("binlowlim");
  binupleaf = ntuple1shf->GetLeaf("binuplim");
  binlowleafy = ntuple1shf->GetLeaf("binlowlimy");
  binupleafy = ntuple1shf->GetLeaf("binuplimy");
  }
  //ntuple={0-15_1S, 0-15_2S, 15-22_1S, 15-22_2S, 22-30_1S, 22-30_2S, 30-120_1S, 30-120_2S}
  for (int ihf = 1; ihf<numhfbins+1; ihf++) {
    if (whichUpsilon==1) ntuplehf0to193->GetEntry(2*ihf-1);
    else if (whichUpsilon==2) ntuplehf0to193->GetEntry(2*ihf);
    else if (whichUpsilon==3) ntuple1shf->GetEntry(ihf);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("hf%.0f-%.0f_y%.2f-%.2f",binlowlim,binuplim,binlowlimy,binuplimy) << endl;
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    hHFSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
  }


  //Get errors in ntracks bins
  cout << "--------------------------------" << endl;
  cout << "ntracks bins" << endl;
/*  RpPbErrLeaf = ntuple1sntracks->GetLeaf("RFBErr");
  binlowleaf = ntuple1sntracks->GetLeaf("binlowlim");
  binupleaf = ntuple1sntracks->GetLeaf("binuplim");
  binlowleafy = ntuple1sntracks->GetLeaf("binlowlimy");
  binupleafy = ntuple1sntracks->GetLeaf("binuplimy");
  for (int whichBin=1; whichBin<=numbins; whichBin++) {
  for (int ihf = 1; ihf<numntbins+1; ihf++) {
    if (whichUpsilon==1) {
      ntuple1sntracks->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-5+whichBin);
    }
    if (whichUpsilon==2) {
      ntuple1sntracks->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-3+whichBin);
    }
    if (whichUpsilon==3) {
      ntuple1sntracks->GetEntry(ihf*4/pow(2,(whichUpsilon-1))-1);
    }
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("ntracks%.0f-%.0f_y%.2f-%.2f",binlowlim,binuplim,binlowlimy,binuplimy) << endl;
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    if (whichBin==1) {
      if (whichUpsilon==1) hNtracksSysSigRFB000to040->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==2) hNtracksSysSigRFB000to080->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==3) hNtracksSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==2) {
      if (whichUpsilon==1) hNtracksSysSigRFB040to080->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
      else if (whichUpsilon==2) hNtracksSysSigRFB080to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==3) {
      hNtracksSysSigRFB080to120->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
    else if (whichBin==4) {
      hNtracksSysSigRFB120to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
    }
  }
  }*/

  if (whichUpsilon<3){
  //Get errors in ntracks bins with |y|<1.93 for 1s and 2s
  RpPbErrLeaf = ntuplent0to193->GetLeaf("RFBErr");
  binlowleaf = ntuplent0to193->GetLeaf("binlowlim");
  binupleaf = ntuplent0to193->GetLeaf("binuplim");
  binlowleafy = ntuplent0to193->GetLeaf("binlowlimy");
  binupleafy = ntuplent0to193->GetLeaf("binuplimy");
  }
  else {
  RpPbErrLeaf = ntuple1sntracks->GetLeaf("RFBErr");
  binlowleaf = ntuple1sntracks->GetLeaf("binlowlim");
  binupleaf = ntuple1sntracks->GetLeaf("binuplim");
  binlowleafy = ntuple1sntracks->GetLeaf("binlowlimy");
  binupleafy = ntuple1sntracks->GetLeaf("binuplimy");
  }
  for (int ihf = 1; ihf<numntbins+1; ihf++) {
    if (whichUpsilon==1) ntuplent0to193->GetEntry(2*ihf-1);
    if (whichUpsilon==2) ntuplent0to193->GetEntry(2*ihf);
    if (whichUpsilon==3) ntuple1sntracks->GetEntry(ihf);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("ntracks%.0f-%.0f_y%.2f-%.2f",binlowlim,binuplim,binlowlimy,binuplimy) << endl;
    RpPbErr = (double)RpPbErrLeaf->GetValue();
    hNtracksSysSigRFB000to193->SetBinContent(ihf, TMath::Abs(RpPbErr)/100);
  }

  //Get errors in integrated activity bin
  cout << "--------------------------------" << endl;
  cout << "Integrated activity bin" << endl;
  RpPbErrLeaf = ntupleFullHF->GetLeaf("RFBErr");
  binlowleaf = ntupleFullHF->GetLeaf("binlowlim");
  binupleaf = ntupleFullHF->GetLeaf("binuplim");
  binlowleafy = ntupleFullHF->GetLeaf("binlowlimy");
  binupleafy = ntupleFullHF->GetLeaf("binuplimy");
    ntupleFullHF->GetEntry(0);
    binlowlim = (float)binlowleaf->GetValue();
    binuplim = (float)binupleaf->GetValue();
    binlowlimy = (float)binlowleafy->GetValue();
    binuplimy = (float)binupleafy->GetValue();
    cout << Form("hfsum%.2f-%.2f_y%.2f-%.2f",binlowlim,binuplim,binlowlimy,binuplimy) << endl;
    RpPbErr = (double)RpPbErrLeaf->GetValue();
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

  hrapSysSigCrossLowPt->Write();
  hrapSysSigCrossHighPt->Write();
  hSignalErrptPAy287to193->Write();
  hrapSysSigCross->Write();

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
  }
  hHFSysSigRFB000to193->Write();
  hNtracksSysSigRFB000to193->Write();
  hSysSigRFBIntActivity->Write();

  outFile.Close();
}
