#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviationsy193to000to193_Combine(int whichUpsilon=1) {

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

  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TString ntuplename = Form("ntuple%is;1",whichUpsilon);

  TString outfilename = Form("FixParamErrors/ParamDev%isCombined_y193to000to193.root",whichUpsilon);
  TFile outFile(outfilename, "RECREATE");
  cout << outfilename << endl;

  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //choose a set of bins
  float ptbins1[7] = {0,2,4,6,9,12,30};
  float ptbins2[4] = {0,4,9,30};
  float ptbins3[3] = {0,6,30};

  float* ptbinsptr;
  int numptbinstemp;
  if (whichUpsilon==1) {
    ptbinsptr = &ptbins1[0];
    numptbinstemp = sizeof(ptbins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ptbinsptr = &ptbins2[0];
    numptbinstemp = sizeof(ptbins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ptbinsptr = &ptbins3[0];
    numptbinstemp = sizeof(ptbins3)/sizeof(float)-1;
  }

  const int numptbins = numptbinstemp;

  //Define ntuples
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);

  //Load the individual ntuples
  TString fileNamealpha = Form("FixParamErrors/ParamalphaDev%is_y193to000to193.root",whichUpsilon);
  TFile *inFilealpha = TFile::Open(fileNamealpha,"READ");
  TNtuple* ntuple1sptalpha = (TNtuple*)inFilealpha->Get(ntupleptname);

  TString fileNamen = Form("FixParamErrors/ParamnDev%is_y193to000to193.root",whichUpsilon);
  TFile *inFilen = TFile::Open(fileNamen,"READ");
  TNtuple* ntuple1sptn = (TNtuple*)inFilen->Get(ntupleptname);

  TString fileNamex = Form("FixParamErrors/ParamxDev%is_y193to000to193.root",whichUpsilon);
  TFile *inFilex = TFile::Open(fileNamex,"READ");
  TNtuple* ntuple1sptx = (TNtuple*)inFilex->Get(ntupleptname);

  TString fileNamef = Form("FixParamErrors/ParamfDev%is_y193to000to193.root",whichUpsilon);
  TFile *inFilef = TFile::Open(fileNamef,"READ");
  TNtuple* ntuple1sptf = (TNtuple*)inFilef->Get(ntupleptname);

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "PT BINS IN FORWARD AND BACKWARD Y" << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel;

  //BIN LOOP********************************************************
  yLowPP = 0.00;
  yHighPP = 1.93;
  for (int iy = 0; iy<2; iy++) {
    if (iy==0) {
      yLowCM = -1.93;
      yHighCM = 0.0;
    }
    else {
      yLowCM = 0.0;
      yHighCM = 1.93;
    }
  for (int ipt = 0; ipt<numptbins; ipt++) {

    ptLow = *(ptbinsptr+ipt);
    ptHigh = *(ptbinsptr+ipt+1);
    binLow = ptLow;
    binHigh = ptHigh;
    binvar = "pt";

    if (yLowCM<-2.5) {
      PPtoo = kFALSE;
      //cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;

    //print bin label
    binLabel = Form("pt%.2f-%.2f_y%.2f-%.2f",ptLow,ptHigh,yLowCM,yHighCM);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //Take largest of the errors
    float temperr = 0;
    float PPerr = 0;
    float RpAerr = 0;

    TLeaf *pp1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafalpha = ntuple1sptalpha->GetLeaf(strRpPbErr);
    ntuple1sptalpha->GetEntry(ipt+iy*numptbins);
    TLeaf *pp1sErrLeafn = ntuple1sptn->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafn = ntuple1sptn->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafn = ntuple1sptn->GetLeaf(strRpPbErr);
    ntuple1sptn->GetEntry(ipt+iy*numptbins);
    TLeaf *pp1sErrLeafx = ntuple1sptx->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafx = ntuple1sptx->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafx = ntuple1sptx->GetLeaf(strRpPbErr);
    ntuple1sptx->GetEntry(ipt+iy*numptbins);
    TLeaf *pPb1sErrLeaff = ntuple1sptf->GetLeaf(strpPb1sErr);
    ntuple1sptf->GetEntry(ipt+iy*numptbins);

    TLeaf *binlowLeaf = ntuple1sptalpha->GetLeaf("binlowlim");
    TLeaf *binupLeaf = ntuple1sptalpha->GetLeaf("binuplim");
    ptLow = (double)binlowLeaf->GetValue();
    ptHigh = (double)binupLeaf->GetValue();
    //print bin label again
    binLabel = Form("pt%.2f-%.2f_y%.2f-%.2f",ptLow,ptHigh,yLowCM,yHighCM);
    cout << endl << setw(binColWidth) << setfill(separator) << binLabel;

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
    if (whichUpsilon==1) ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2) ntuple2spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==3) ntuple3spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);

  }// end of pt loop
  }// end of y loop

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1spt->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2spt->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3spt->Write();
  }
  delete ntuple1spt;
  delete ntuple2spt;
  delete ntuple3spt;

  outFile.Close();
}
