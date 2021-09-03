#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviationsXSBins_Combine(int whichUpsilon=1, bool takeAvg=kFALSE) {

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

  TString outfilename = Form("FixParamErrors/ParamDev%isCombinedXSBins.root",whichUpsilon);
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
  TString fileNamealpha = Form("FixParamErrors/ParamalphaDev%isXSBins.root",whichUpsilon);
  TFile *inFilealpha = TFile::Open(fileNamealpha,"READ");
  TNtuple* ntuple1sptalpha = (TNtuple*)inFilealpha->Get(ntupleptname);

  TString fileNamen = Form("FixParamErrors/ParamnDev%isXSBins.root",whichUpsilon);
  TFile *inFilen = TFile::Open(fileNamen,"READ");
  TNtuple* ntuple1sptn = (TNtuple*)inFilen->Get(ntupleptname);

  TString fileNamex = Form("FixParamErrors/ParamxDev%isXSBins.root",whichUpsilon);
  TFile *inFilex = TFile::Open(fileNamex,"READ");
  TNtuple* ntuple1sptx = (TNtuple*)inFilex->Get(ntupleptname);

  TString fileNamef = Form("FixParamErrors/ParamfDev%isXSBins.root",whichUpsilon);
  TFile *inFilef = TFile::Open(fileNamef,"READ");
  TNtuple* ntuple1sptf = (TNtuple*)inFilef->Get(ntupleptname);

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "PT BINS IN FULL Y RANGE" << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << endl;

  float yLowCM, yHighCM, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel;

  //BIN LOOP********************************************************
  for (int ipt = -1; ipt<numptbins; ipt++) {

    yLowCM = -2.87;
    yHighCM = 1.93;
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

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    binLabel = strbinLow+"<pt<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //Take largest of the errors
    float temperr = 0;
    float PPerr = 0;
    float RpAerr = 0;

    TLeaf *pPb1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpPb1sErr);
    ntuple1sptalpha->GetEntry(ipt+1);
    TLeaf *pPb1sErrLeafn = ntuple1sptn->GetLeaf(strpPb1sErr);
    ntuple1sptn->GetEntry(ipt+1);
    TLeaf *pPb1sErrLeafx = ntuple1sptx->GetLeaf(strpPb1sErr);
    ntuple1sptx->GetEntry(ipt+1);
    TLeaf *pPb1sErrLeaff = ntuple1sptf->GetLeaf(strpPb1sErr);
    ntuple1sptf->GetEntry(ipt+1);

    //print bin label again
    cout << endl;
    TLeaf *binlowLeaf = ntuple1sptalpha->GetLeaf("binlowlim");
    TLeaf *binupLeaf = ntuple1sptalpha->GetLeaf("binuplim");
    ptLow = (double)binlowLeaf->GetValue();
    ptHigh = (double)binupLeaf->GetValue();
    stringstream stream3;
    stream3 << fixed << setprecision(2) << binLow;
    strbinLow = stream3.str();
    stringstream stream4;
    stream4 << fixed << setprecision(2) << binHigh;
    strbinHigh = stream4.str();
    binLabel = strbinLow+"<pt<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    double pPb1sErralpha = (double)pPb1sErrLeafalpha->GetValue();
    if (pPb1sErralpha>temperr) temperr = pPb1sErralpha;

    double pPb1sErrn = (double)pPb1sErrLeafn->GetValue();
    if (pPb1sErrn>temperr) temperr = pPb1sErrn;

    double pPb1sErrx = (double)pPb1sErrLeafx->GetValue();
    if (pPb1sErrx>temperr) temperr = pPb1sErrx;

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
