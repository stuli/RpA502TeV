#include "TNtuple.h"
#include "TFile.h"
#include "TString.h"
#include "../../HeaderFiles/cutsAndBin.h"

void GetFixParamDeviations_HFNtracks_Combine(int whichUpsilon=1,int binmode=0, bool takeAvg=kFALSE) {

  TString avg = "";
  if (takeAvg) avg = "avg";

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;

  TString ntupleptname = Form("ntuple%ispt%s;1",1,avg.Data());
  TString ntupleyname = Form("ntuple%isy%s;1",1,avg.Data());
  TString ntuplename = Form("ntuple%is%s;1",1,avg.Data());

  //choose a set of bins
  float ybinsCM[3] = {-1.93,0.0,1.93};
  float hfbins3[3] = {0,12,120};
  int ntbins3[3] = {0,40,400};
  float hfbins1[5] = {0,12,19,27,120};
  int ntbins1[5] = {0,40,62,88,400};

  float* hfbinsptr;
  int* ntbinsptr;
  int numhfbinstemp, numntbinstemp;
  if (whichUpsilon==3) {
    hfbinsptr = &hfbins3[0];
    ntbinsptr = &ntbins3[0];
    numhfbinstemp = sizeof(hfbins3)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins3)/sizeof(float)-1;
  }
  else {
    hfbinsptr = &hfbins1[0];
    ntbinsptr = &ntbins1[0];
    numhfbinstemp = sizeof(hfbins1)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins1)/sizeof(float)-1;
  }

  const int numybins = 2;
  const int numhfbins = numhfbinstemp;
  const int numntbins = numntbinstemp;

  TString hfntracksbins = "_hfbins.root";
  int numbins = numhfbins;
  if (binmode==1) {
    hfntracksbins = "_ntracksbins.root";
    numbins = numntbins;
  }

  TString outfilename = Form("FixParamErrors/ParamDev%isCombined",whichUpsilon);
  outfilename = outfilename + hfntracksbins;
  TFile outFile(outfilename, "RECREATE");
  cout << outfilename << endl;

  TString strpPb1sBErr = Form("pPb%isBerr",1);
  TString strpPb1sFErr = Form("pPb%isFerr",1);
  TString strRFBErr = "RFBErr";

  //Define ntuples
  TNtuple* ntuple1s = new TNtuple(Form("ntuple1s%s",avg.Data()),"Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb1sFerr:pPb1sBerr:RFBErr",numybins);

  //Load the individual ntuples
  TString fileNamealpha = Form("FixParamErrors/ParamalphaDev%is",whichUpsilon);
  fileNamealpha = fileNamealpha + hfntracksbins;
  TFile *inFilealpha = TFile::Open(fileNamealpha,"READ");
  TNtuple* ntuple1salpha = (TNtuple*)inFilealpha->Get(ntuplename);

  TString fileNamen = Form("FixParamErrors/ParamnDev%is",whichUpsilon);
  fileNamen = fileNamen + hfntracksbins;
  TFile *inFilen = TFile::Open(fileNamen,"READ");
  TNtuple* ntuple1sn = (TNtuple*)inFilen->Get(ntuplename);

  TString fileNamex = Form("FixParamErrors/ParamxDev%is",whichUpsilon);
  fileNamex = fileNamex + hfntracksbins;
  TFile *inFilex = TFile::Open(fileNamex,"READ");
  TNtuple* ntuple1sx = (TNtuple*)inFilex->Get(ntuplename);

  TString fileNamef = Form("FixParamErrors/ParamfDev%is",whichUpsilon);
  fileNamef = fileNamef + hfntracksbins;
  TFile *inFilef = TFile::Open(fileNamef,"READ");
  TNtuple* ntuple1sf = (TNtuple*)inFilef->Get(ntuplename);

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -30;
  int max = 30;

  TString whichBins = "HF";
  if (binmode==1) whichBins = "NTRACKS";
  cout << whichUpsilon << "S " << whichBins << " BINS IN FORWARD AND BACKWARD Y" << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA F %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA B %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RFB ERR";
  cout << endl;

  float yLowB, yHighB, yLowF, yHighF, hfLow, hfHigh;
  int ntracksLow, ntracksHigh;
  float ptLow = 0;
  float ptHigh = 30;
  TString binLabel, kineLabelF, kineLabelB, histFileNameF, histFileNameB;
  TFile* theFileF;
  TFile* theFileB;

  //BIN LOOP********************************************************
  for (int ihf = -1; ihf<numbins; ihf++) {
    //Choose the hf bin
    if (ihf<0) {
      hfLow = 0;
      hfHigh = 120;
      ntracksLow = 0;
      ntracksHigh = 400;
    }
    else {
      if (binmode==1) {
        hfLow = 0;
        hfHigh = 120;
        ntracksLow = *(ntbinsptr+ihf);
        ntracksHigh = *(ntbinsptr+ihf+1);
      }
      else {
        hfLow = *(hfbinsptr+ihf);
        hfHigh = *(hfbinsptr+ihf+1);
        ntracksLow = 0;
        ntracksHigh = 400;
      }
    }
  for (int iy = 0; iy<numybins/2; iy++) {

    //Choose the rapidity bin
    yLowB = ybinsCM[numybins/2-iy-1];
    yHighB = ybinsCM[numybins/2-iy];
    yLowF = ybinsCM[iy+numybins/2];
    yHighF = ybinsCM[iy+numybins/2+1];

    //print bin label
    if (binmode==0) binLabel = Form("hf%.1f-%.1f_y%.2f-%.2f",hfLow,hfHigh,yLowF,yHighF);
    else if (binmode==1) binLabel = Form("ntracks%i-%i_y%.2f-%.2f",ntracksLow,ntracksHigh,yLowF,yHighF);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //Take largest of the errors
    float temperrF = 0;
    float temperrB = 0;
    float RFBErr = 0;

    TLeaf *pPb1sBErrLeafalpha = ntuple1salpha->GetLeaf(strpPb1sBErr);
    TLeaf *pPb1sFErrLeafalpha = ntuple1salpha->GetLeaf(strpPb1sFErr);
    TLeaf *RFBErrLeafalpha = ntuple1salpha->GetLeaf(strRFBErr);
    ntuple1salpha->GetEntry(ihf+1+iy*numhfbins);
    TLeaf *pPb1sBErrLeafn = ntuple1sn->GetLeaf(strpPb1sBErr);
    TLeaf *pPb1sFErrLeafn = ntuple1sn->GetLeaf(strpPb1sFErr);
    TLeaf *RFBErrLeafn = ntuple1sn->GetLeaf(strRFBErr);
    ntuple1sn->GetEntry(ihf+1+iy*numhfbins);
    TLeaf *pPb1sBErrLeafx = ntuple1sx->GetLeaf(strpPb1sBErr);
    TLeaf *pPb1sFErrLeafx = ntuple1sx->GetLeaf(strpPb1sFErr);
    TLeaf *RFBErrLeafx = ntuple1sx->GetLeaf(strRFBErr);
    ntuple1sx->GetEntry(ihf+1+iy*numhfbins);
    TLeaf *pPb1sBErrLeaff = ntuple1sf->GetLeaf(strpPb1sBErr);
    TLeaf *pPb1sFErrLeaff = ntuple1sf->GetLeaf(strpPb1sFErr);
    TLeaf *RFBErrLeaff = ntuple1sf->GetLeaf(strRFBErr);
    ntuple1sf->GetEntry(ihf+1+iy*numhfbins);

    //print bin label again
    cout << endl;
    TLeaf *binlowLeaf = ntuple1salpha->GetLeaf("binlowlim");
    TLeaf *binupLeaf = ntuple1salpha->GetLeaf("binuplim");
    ptLow = (double)binlowLeaf->GetValue();
    ptHigh = (double)binupLeaf->GetValue();
    if (binmode==0) binLabel = Form("hf%.1f-%.1f_y%.2f-%.2f",hfLow,hfHigh,yLowF,yHighF);
    else if (binmode==1) binLabel = Form("ntracks%i-%i_y%.2f-%.2f",ntracksLow,ntracksHigh,yLowF,yHighF);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    double pPb1sBErralpha = (double)pPb1sBErrLeafalpha->GetValue();
    double pPb1sFErralpha = (double)pPb1sFErrLeafalpha->GetValue();
    double RFBErralpha = (double)RFBErrLeafalpha->GetValue();
    if (pPb1sBErralpha>temperrB) temperrB = pPb1sBErralpha;
    if (pPb1sFErralpha>temperrF) temperrF = pPb1sFErralpha;
    if (RFBErralpha>RFBErr) RFBErr = RFBErralpha;

    double pPb1sBErrn = (double)pPb1sBErrLeafn->GetValue();
    double pPb1sFErrn = (double)pPb1sFErrLeafn->GetValue();
    double RFBErrn = (double)RFBErrLeafn->GetValue();
    if (pPb1sBErrn>temperrB) temperrB = pPb1sBErrn;
    if (pPb1sFErrn>temperrF) temperrF = pPb1sFErrn;
    if (RFBErrn>RFBErr) RFBErr = RFBErrn;

    double pPb1sBErrx = (double)pPb1sBErrLeafx->GetValue();
    double pPb1sFErrx = (double)pPb1sFErrLeafx->GetValue();
    double RFBErrx = (double)RFBErrLeafx->GetValue();
    if (pPb1sBErrx>temperrB) temperrB = pPb1sBErrx;
    if (pPb1sFErrx>temperrF) temperrF = pPb1sFErrx;
    if (RFBErrx>RFBErr) RFBErr = RFBErrx;

    double pPb1sBErrf = (double)pPb1sBErrLeaff->GetValue();
    double pPb1sFErrf = (double)pPb1sFErrLeaff->GetValue();
    double RFBErrf = (double)RFBErrLeaff->GetValue();
    if (pPb1sBErrf>temperrB) temperrB = pPb1sBErrf;
    if (pPb1sFErrf>temperrF) temperrF = pPb1sFErrf;
    if (RFBErrf>RFBErr) RFBErr = RFBErrf;

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << temperrF;
    cout << setw(Width1S) << setfill(separator) << temperrB;
    cout << setw(Width1S) << setfill(separator) << RFBErr;
    cout << endl;

    //put errors in ntuple
    int binlow = (int)hfLow;
    int binhigh = (int)hfHigh;
    if (binmode==1) {
      binlow = ntracksLow;
      binhigh = ntracksHigh;
    }
    ntuple1s->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBErr);

  }//end of rapidity bin loop
  }//end of hf bin loop

  //Write errors in ntuple file
  outFile.cd();
  ntuple1s->Write();
  delete ntuple1s;

  outFile.Close();

}
