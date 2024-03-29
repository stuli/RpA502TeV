#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviationspt0to6to30_in3Sbins_Combine(int whichUpsilon=1, bool takeAvg=kFALSE) {

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

  TString outfilename = Form("FixParamErrors/ParamDev%isCombined_pt0to6to30_in3Sbins.root",whichUpsilon);
  TFile outFile(outfilename, "RECREATE");
  cout << outfilename << endl;

  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //choose a set of bins
  float ybinsCM[3] = {-1.93,0.0,1.93};

  const int numybins = sizeof(ybinsCM)/sizeof(float)-1;

  //Define ntuples
  TNtuple* ntuple1sy = new TNtuple(Form("ntuple1sy%s",avg.Data()),"Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

  //Load the individual ntuples
  TString fileNamealpha = Form("FixParamErrors/ParamalphaDev%is_pt0to6to30_in3Sbins.root",whichUpsilon);
  TFile *inFilealpha = TFile::Open(fileNamealpha,"READ");
  TNtuple* ntuple1syalpha = (TNtuple*)inFilealpha->Get(ntupleyname);

  TString fileNamen = Form("FixParamErrors/ParamnDev%is_pt0to6to30_in3Sbins.root",whichUpsilon);
  TFile *inFilen = TFile::Open(fileNamen,"READ");
  TNtuple* ntuple1syn = (TNtuple*)inFilen->Get(ntupleyname);

  TString fileNamex = Form("FixParamErrors/ParamxDev%is_pt0to6to30_in3Sbins.root",whichUpsilon);
  TFile *inFilex = TFile::Open(fileNamex,"READ");
  TNtuple* ntuple1syx = (TNtuple*)inFilex->Get(ntupleyname);

  TString fileNamef = Form("FixParamErrors/ParamfDev%is_pt0to6to30_in3Sbins.root",whichUpsilon);
  TFile *inFilef = TFile::Open(fileNamef,"READ");
  TNtuple* ntuple1syf = (TNtuple*)inFilef->Get(ntupleyname);

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "RAPIDITY BINS IN LOW AND HIGH PT" << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel;

  //BIN LOOP********************************************************
  for (int ipt = 0; ipt<2; ipt++) {
    if (ipt==0){
      ptLow = 0;
      ptHigh = 6;
    }
    else {
      ptLow = 6;
      ptHigh = 30;
    }
  for (int iy = 0; iy<numybins; iy++) {

    yLowCM = ybinsCM[iy];
    yHighCM = ybinsCM[iy+1];
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
    binLabel = strbinLow+"<y<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //Take largest of the errors
    float temperr = 0;
    float PPerr = 0;
    float RpAerr = 0;

    TLeaf *pp1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafalpha = ntuple1syalpha->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafalpha = ntuple1syalpha->GetLeaf(strRpPbErr);
    ntuple1syalpha->GetEntry(iy+ipt*numybins);
    TLeaf *pp1sErrLeafn = ntuple1syn->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafn = ntuple1syn->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafn = ntuple1syn->GetLeaf(strRpPbErr);
    ntuple1syn->GetEntry(iy+ipt*numybins);
    TLeaf *pp1sErrLeafx = ntuple1syx->GetLeaf(strpp1sErr);
    TLeaf *pPb1sErrLeafx = ntuple1syx->GetLeaf(strpPb1sErr);
    TLeaf *RpPbErrLeafx = ntuple1syx->GetLeaf(strRpPbErr);
    ntuple1syx->GetEntry(iy+ipt*numybins);
    TLeaf *pPb1sErrLeaff = ntuple1syf->GetLeaf(strpPb1sErr);
    ntuple1syf->GetEntry(iy+ipt*numybins);

    //print bin label again
    cout << endl;
    TLeaf *binlowLeaf = ntuple1syalpha->GetLeaf("binlowlim");
    TLeaf *binupLeaf = ntuple1syalpha->GetLeaf("binuplim");
    binLow = (double)binlowLeaf->GetValue();
    binHigh = (double)binupLeaf->GetValue();
    stringstream stream3;
    stream3 << fixed << setprecision(2) << binLow;
    stringstream stream4;
    stream4 << fixed << setprecision(2) << binHigh;
    binLabel = strbinLow+"<y<"+strbinHigh;
    cout << setw(binColWidth) << setfill(separator) << binLabel;

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
    ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);

  }// end of y loop
  }// end of pt loop

  //Write errors in ntuple file
  outFile.cd();
  ntuple1sy->Write();
  delete ntuple1sy;

  outFile.Close();
}
