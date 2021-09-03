#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviations_Average(int whichUpsilon=1) {

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

  TString outfilename = Form("FixParamErrors/ParamDev%isAverage.root",whichUpsilon);
  TFile outFile(outfilename, "RECREATE");
  cout << outfilename << endl;

  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybinsCM[10] = {-2.87,-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybinsCM[6] = {-2.87,-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybinsCM[4] = {-2.87,-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;
  const int numtot = numptbins + numybins;

  //Define ntuples
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple1sy = new TNtuple("ntuple1sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2sy = new TNtuple("ntuple2sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3sy = new TNtuple("ntuple3sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

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
        ptLow = ptbins[ipt];
        ptHigh = ptbins[ipt+1];
      }
      binLow = ptLow;
      binHigh = ptHigh;
      binvar = "pt";
    }
    else {
      ptLow = 0;
      ptHigh = 30;
      yLowCM = ybinsCM[ipt-numptbins];
      yHighCM = ybinsCM[ipt-numptbins+1];
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

    if (ipt<numptbins){
      TLeaf *pp1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafalpha = ntuple1sptalpha->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafalpha = ntuple1sptalpha->GetLeaf(strRpPbErr);
      ntuple1sptalpha->GetEntry(ipt+1);
      TLeaf *pp1sErrLeafn = ntuple1sptn->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafn = ntuple1sptn->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafn = ntuple1sptn->GetLeaf(strRpPbErr);
      ntuple1sptn->GetEntry(ipt+1);
      TLeaf *pp1sErrLeafx = ntuple1sptx->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafx = ntuple1sptx->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafx = ntuple1sptx->GetLeaf(strRpPbErr);
      ntuple1sptx->GetEntry(ipt+1);
      TLeaf *pp1sErrLeaff = ntuple1sptf->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeaff = ntuple1sptf->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeaff = ntuple1sptf->GetLeaf(strRpPbErr);
      ntuple1sptf->GetEntry(ipt+1);
    }
    else {
      TLeaf *pp1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafalpha = ntuple1syalpha->GetLeaf(strRpPbErr);
      ntuple1syalpha->GetEntry(ipt-numptbins);
      TLeaf *pp1sErrLeafn = ntuple1syn->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafn = ntuple1syn->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafn = ntuple1syn->GetLeaf(strRpPbErr);
      ntuple1syn->GetEntry(ipt-numptbins);
      TLeaf *pp1sErrLeafx = ntuple1syx->GetLeaf(strpp1sErr);
      TLeaf *pPb1sErrLeafx = ntuple1syx->GetLeaf(strpPb1sErr);
      TLeaf *RpPbErrLeafx = ntuple1syx->GetLeaf(strRpPbErr);
      ntuple1syx->GetEntry(ipt-numptbins);
      TLeaf *pPb1sErrLeaff = ntuple1syf->GetLeaf(strpPb1sErr);
      ntuple1syf->GetEntry(ipt-numptbins);
    }

    double pp1sErralpha = (double)pp1sErrLeafalpha->GetValue();
    double pPb1sErralpha = (double)pPb1sErrLeafalpha->GetValue();
    double RpPbErralpha = (double)RpPbErrLeafalpha->GetValue();
    PPerr += pp1sErralpha;
    temperr += pPb1sErralpha;
    RpAerr += RpPbErralpha;

    double pp1sErrn = (double)pp1sErrLeafn->GetValue();
    double pPb1sErrn = (double)pPb1sErrLeafn->GetValue();
    double RpPbErrn = (double)RpPbErrLeafn->GetValue();
    PPerr += pp1sErrn;
    temperr += pPb1sErrn;
    RpAerr += RpPbErrn;

    double pp1sErrx = (double)pp1sErrLeafx->GetValue();
    double pPb1sErrx = (double)pPb1sErrLeafx->GetValue();
    double RpPbErrx = (double)RpPbErrLeafx->GetValue();
    PPerr += pp1sErrx;
    temperr += pPb1sErrx;
    RpAerr += RpPbErrx;

    double pPb1sErrf = (double)pPb1sErrLeaff->GetValue();
    temperr += pPb1sErrf;

    PPerr = PPerr/3;
    temperr = temperr/4;
    RpAerr = RpAerr/3;

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    cout << endl;

    //put errors in ntuple
    if (whichUpsilon==1)
      if (ipt<numptbins)
        ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2)
      if (ipt<numptbins)
        ntuple2spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple2sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==3)
      if (ipt<numptbins)
        ntuple3spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple3sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
  }

  //Write errors in ntuple file
  outFile.cd();
  if (whichUpsilon==1) {
    ntuple1spt->Write();
    ntuple1sy->Write();
  }
  else if (whichUpsilon==2) {
    ntuple2spt->Write();
    ntuple2sy->Write();
  }
  else if (whichUpsilon==3) {
    ntuple3spt->Write();
    ntuple3sy->Write();
  }
  delete ntuple1spt;
  delete ntuple1sy;
  delete ntuple2spt;
  delete ntuple2sy;
  delete ntuple3spt;
  delete ntuple3sy;

  outFile.Close();
}
