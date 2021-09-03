#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviations_Combine(int whichUpsilon=1, bool takeAvg=kFALSE) {

  TString avg = "";
  if (takeAvg) avg = "avg";

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  bool PPtoo = kTRUE;

  int collId = kPADATA;
  int collIdPP = kPPDATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  TString ntupleptname = Form("ntuple%ispt%s;1",1,avg.Data());
  TString ntupleyname = Form("ntuple%isy%s;1",1,avg.Data());
  TString ntuplename = Form("ntuple%is%s;1",1,avg.Data());

  TString outfilename = Form("FixParamErrors/ParamDev%isCombined.root",whichUpsilon);
  TFile outFile(outfilename, "RECREATE");
  cout << outfilename << endl;

  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //choose a set of bins
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ybins1[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ptbins2[4] = {0,4,9,30};
  float ybins2[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  float ptbins3[3] = {0,6,30};
  float ybins3[4] = {-2.87,-1.93,0.0,1.93};

  float* ptbinsptr;
  float* ybinsptr;
  int numptbinstemp, numybinstemp;
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    ybinsptr = &ybins1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
    numybinstemp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    ybinsptr = &ybins2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
    numybinstemp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    ybinsptr = &ybins3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
    numybinstemp = sizeof(ybins3)/sizeof(float)-1;
  }

  const int numptbins = numptbinstemp;
  const int numybins = numybinstemp;
  const int numtot = numptbins + numybins;

  //Define ntuples
  TNtuple* ntuple1spt = new TNtuple(Form("ntuple1spt%s",avg.Data()),"Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple1sy = new TNtuple(Form("ntuple1sy%s",avg.Data()),"Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

  //Load the individual ntuples
  TString fileNamealpha = Form("FixParamErrors/ParamalphaDev%is.root",whichUpsilon);
  TFile *inFilealpha = TFile::Open(fileNamealpha,"READ");
  TNtuple* ntuple1sptalpha = (TNtuple*)inFilealpha->Get(ntupleptname);
  TNtuple* ntuple1syalpha = (TNtuple*)inFilealpha->Get(ntupleyname);

  TString fileNamen = Form("FixParamErrors/ParamnDev%is.root",whichUpsilon);
  TFile *inFilen = TFile::Open(fileNamen,"READ");
  TNtuple* ntuple1sptn = (TNtuple*)inFilen->Get(ntupleptname);
  TNtuple* ntuple1syn = (TNtuple*)inFilen->Get(ntupleyname);

  TString fileNamex = Form("FixParamErrors/ParamxDev%is.root",whichUpsilon);
  TFile *inFilex = TFile::Open(fileNamex,"READ");
  TNtuple* ntuple1sptx = (TNtuple*)inFilex->Get(ntupleptname);
  TNtuple* ntuple1syx = (TNtuple*)inFilex->Get(ntupleyname);

  TString fileNamef = Form("FixParamErrors/ParamfDev%is.root",whichUpsilon);
  TFile *inFilef = TFile::Open(fileNamef,"READ");
  TNtuple* ntuple1sptf = (TNtuple*)inFilef->Get(ntupleptname);
  TNtuple* ntuple1syf = (TNtuple*)inFilef->Get(ntupleyname);

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "PT AND RAPIDITY BINS" << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel;

  //BIN LOOP********************************************************
  for (int ipt = -1; ipt<numtot; ipt++) {

    if (ipt<numptbins){
      yLowCM = -1.93;
      yHighCM = 1.93;
      yLowPP = 0.00;
      yHighPP = 1.93;
      if (ipt<0) {
        ptLow = 0;
        ptHigh = 30;
      }
      else {
        ptLow = *(ptbinsptr+ipt);
        ptHigh = *(ptbinsptr+ipt+1);
      }
      binLow = ptLow;
      binHigh = ptHigh;
      binvar = "pt";
    }
    else {
      ptLow = 0;
      ptHigh = 30;
      yLowCM = *(ybinsptr+ipt-numptbins);
      yHighCM = *(ybinsptr+ipt-numptbins+1);
      binLow = yLowCM;
      binHigh = yHighCM;
      binvar = "y";
      if (yLowCM<0) {
        yLowPP = TMath::Abs(yHighCM);
        yHighPP = TMath::Abs(yLowCM);
      }
      else {
        yLowPP = yLowCM;
        yHighPP = yHighCM;
      }
    }

    if (yLowCM<-2.5) {
      PPtoo = kFALSE;
      //cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    if (ipt<numptbins) binLabel = strbinLow+"<pt<"+strbinHigh;
    else binLabel = strbinLow+"<y<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //Take largest of the errors
    float temperr = 0;
    float PPerr = 0;
    float RpAerr = 0;

    TLeaf *pp1sErrLeafalpha;
    TLeaf *pPb1sErrLeafalpha;
    TLeaf *RpPbErrLeafalpha;
    TLeaf *pp1sErrLeafn;
    TLeaf *pPb1sErrLeafn;
    TLeaf *RpPbErrLeafn;
    TLeaf *pp1sErrLeafx;
    TLeaf *pPb1sErrLeafx;
    TLeaf *RpPbErrLeafx;
    TLeaf *pp1sErrLeaff;
    TLeaf *pPb1sErrLeaff;
    TLeaf *RpPbErrLeaff;
    if (ipt<numptbins){
      pp1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpp1sErr);
      pPb1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpPb1sErr);
      RpPbErrLeafalpha = ntuple1sptalpha->GetLeaf(strRpPbErr);
      ntuple1sptalpha->GetEntry(ipt+1);
      pp1sErrLeafn = ntuple1sptn->GetLeaf(strpp1sErr);
      pPb1sErrLeafn = ntuple1sptn->GetLeaf(strpPb1sErr);
      RpPbErrLeafn = ntuple1sptn->GetLeaf(strRpPbErr);
      ntuple1sptn->GetEntry(ipt+1);
      pp1sErrLeafx = ntuple1sptx->GetLeaf(strpp1sErr);
      pPb1sErrLeafx = ntuple1sptx->GetLeaf(strpPb1sErr);
      RpPbErrLeafx = ntuple1sptx->GetLeaf(strRpPbErr);
      ntuple1sptx->GetEntry(ipt+1);
      pp1sErrLeaff = ntuple1sptf->GetLeaf(strpp1sErr);
      pPb1sErrLeaff = ntuple1sptf->GetLeaf(strpPb1sErr);
      RpPbErrLeaff = ntuple1sptf->GetLeaf(strRpPbErr);
      ntuple1sptf->GetEntry(ipt+1);
    }
    else {
      pp1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpp1sErr);
      pPb1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpPb1sErr);
      RpPbErrLeafalpha = ntuple1syalpha->GetLeaf(strRpPbErr);
      ntuple1syalpha->GetEntry(ipt-numptbins);
      pp1sErrLeafn = ntuple1syn->GetLeaf(strpp1sErr);
      pPb1sErrLeafn = ntuple1syn->GetLeaf(strpPb1sErr);
      RpPbErrLeafn = ntuple1syn->GetLeaf(strRpPbErr);
      ntuple1syn->GetEntry(ipt-numptbins);
      pp1sErrLeafx = ntuple1syx->GetLeaf(strpp1sErr);
      pPb1sErrLeafx = ntuple1syx->GetLeaf(strpPb1sErr);
      RpPbErrLeafx = ntuple1syx->GetLeaf(strRpPbErr);
      ntuple1syx->GetEntry(ipt-numptbins);
      pPb1sErrLeaff = ntuple1syf->GetLeaf(strpPb1sErr);
      ntuple1syf->GetEntry(ipt-numptbins);
    }

    double pp1sErralpha = (double)pp1sErrLeafalpha->GetValue();
    double pPb1sErralpha = (double)pPb1sErrLeafalpha->GetValue();
    double RpPbErralpha = (double)RpPbErrLeafalpha->GetValue();
    if (pp1sErralpha>PPerr) PPerr = pp1sErralpha;
    if (pPb1sErralpha>temperr) temperr = pPb1sErralpha;
    if (RpPbErralpha>RpAerr) RpAerr = RpPbErralpha;

    double pp1sErrn = (double)pp1sErrLeafn->GetValue();
    double pPb1sErrn = (double)pPb1sErrLeafn->GetValue();
    double RpPbErrn = (double)RpPbErrLeafn->GetValue();
    if (pp1sErrn>PPerr) PPerr = pp1sErrn;
    if (pPb1sErrn>temperr) temperr = pPb1sErrn;
    if (RpPbErrn>RpAerr) RpAerr = RpPbErrn;

    double pp1sErrx = (double)pp1sErrLeafx->GetValue();
    double pPb1sErrx = (double)pPb1sErrLeafx->GetValue();
    double RpPbErrx = (double)RpPbErrLeafx->GetValue();
    if (pp1sErrx>PPerr) PPerr = pp1sErrx;
    if (pPb1sErrx>temperr) temperr = pPb1sErrx;
    if (RpPbErrx>RpAerr) RpAerr = RpPbErrx;

    double pPb1sErrf = (double)pPb1sErrLeaff->GetValue();
    if (pPb1sErrf>temperr) temperr = pPb1sErrf;

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    cout << endl;

    //put errors in ntuple
      if (ipt<numptbins)
        ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
  }

  //Write errors in ntuple file
  outFile.cd();
  ntuple1spt->Write();
  ntuple1sy->Write();
  delete ntuple1spt;
  delete ntuple1sy;

  outFile.Close();
}
