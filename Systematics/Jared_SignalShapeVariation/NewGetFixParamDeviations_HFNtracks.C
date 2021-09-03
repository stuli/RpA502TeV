#include <iomanip>
#include <sstream>
#include "CheckFitQualitySlim.C"

using namespace std;
void NewGetFixParamDeviations_HFNtracks(int whichUpsilon=1,int binmode=0, int param=1) {

  TString whichParam;
  if (param==1) whichParam = "alpha";
  else if (param==2) whichParam = "f";
  else if (param==3) whichParam = "n";
  else if (param==4) whichParam = "x";

//0=hfmode, 1=ntmode

  TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_05_08/";

  gStyle->SetOptFit();
  gStyle->SetStatW(0.4);

  int collId = kPADATA;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

  float dphiEp2Low = 0;
  float dphiEp2High = 100;


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

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -30;
  int max = 30;

  TString whichBins = "HF";
  if (binmode==1) whichBins = "NTRACKS";
  cout << whichUpsilon << "S " << whichBins << " BINS IN FORWARD AND BACKWARD Y FOR " << whichParam.Data() << endl;
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

  TString hfntracksbins = "_hfbins.root";
  int numbins = numhfbins;
  if (binmode==1) {
    hfntracksbins = "_ntracksbins.root";
    numbins = numntbins;
  }

  TString outFileName = Form("FixParamErrors/Param%sDev%is",whichParam.Data(),whichUpsilon);
  outFileName = outFileName + hfntracksbins;
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1s = new TNtuple("ntuple1s","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb1sFerr:pPb1sBerr:RFBErr",numybins);
  TNtuple* ntuple1savg = new TNtuple("ntuple1savg","Error estimates in y bins","binlowlim:binuplim:binlowlimy:binuplimy:pPb1sFerr:pPb1sBerr:RFBErr",numybins);

  float temperrFMean = 0;
  float temperrBMean = 0;
  float RFBerrMean = 0;
  int FAltGood = 0;
  int BAltGood = 0;
  int Fbins = 0;
  int Bbins = 0;
  int RFBbins = 0;

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

    //import results of pseudo-experiments
    kineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelF = kineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString NomFileNameF = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabelF.Data());
    TFile* NomFileF = TFile::Open(NomFileNameF,"READ");
    RooWorkspace *NomwsF = (RooWorkspace*)NomFileF->Get("workspace");
    NomFileF->Close("R");
    float NomYieldF = NomwsF->var(Form("nSig%is",whichUpsilon))->getVal();
    delete NomwsF;
    delete NomFileF;

    kineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) kineLabelB = kineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString NomFileNameB = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabelB.Data());
    TFile* NomFileB = TFile::Open(NomFileNameB,"READ");
    RooWorkspace *NomwsB = (RooWorkspace*)NomFileB->Get("workspace");
    NomFileB->Close("R");
    float NomYieldB = NomwsB->var(Form("nSig%is",whichUpsilon))->getVal();
    delete NomwsB;
    delete NomFileB;

    TString AltkineLabelF = getKineLabel (collId, ptLow, ptHigh, yLowF, yHighF, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) AltkineLabelF = AltkineLabelF + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltFileNameF = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabelF.Data());
    TFile* AltFileF = TFile::Open(AltFileNameF,"READ");
    RooWorkspace *AltwsF = (RooWorkspace*)AltFileF->Get("workspace");
    AltFileF->Close("R");
    float AltYieldF = AltwsF->var(Form("nSig%is",whichUpsilon))->getVal();
    FAltGood = CheckFitQualitySlim(collId,ptLow,ptHigh,yLowF,yHighF,AltwsF);
    delete AltwsF;
    delete AltFileF;

    TString AltkineLabelB = getKineLabel (collId, ptLow, ptHigh, yLowB, yHighB, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    if (ihf>=0) AltkineLabelB = AltkineLabelB + Form("_hfsum%.2f-%.2f_ntracks%i-%i",hfLow, hfHigh, ntracksLow, ntracksHigh );
    TString AltFileNameB = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabelB.Data());
    TFile* AltFileB = TFile::Open(AltFileNameB,"READ");
    RooWorkspace *AltwsB = (RooWorkspace*)AltFileB->Get("workspace");
    AltFileB->Close("R");
    float AltYieldB = AltwsB->var(Form("nSig%is",whichUpsilon))->getVal();
    BAltGood = CheckFitQualitySlim(collId,ptLow,ptHigh,yLowF,yHighF,AltwsB);
    delete AltwsB;
    delete AltFileB;

    //Calculate errors
    float temperrF = TMath::Abs((AltYieldF-NomYieldF)/NomYieldF*100);
    float temperrB = TMath::Abs((AltYieldB-NomYieldB)/NomYieldB*100);
    float NomRFB = NomYieldF/NomYieldB;
    float AltRFB = AltYieldF/AltYieldB;
    float RFBerr = TMath::Abs((AltRFB-NomRFB)/NomRFB*100);

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << temperrF;
    cout << setw(Width1S) << setfill(separator) << temperrB;
    cout << setw(Width1S) << setfill(separator) << RFBerr;
    cout << endl;

    //Sum up errors to get mean values
    if (FAltGood) {Fbins++; temperrFMean = temperrFMean + temperrF;}
    if (BAltGood) {Bbins++; temperrBMean = temperrBMean + temperrB;}
    if (FAltGood && BAltGood) {RFBbins++; RFBerrMean = RFBerrMean + RFBerr;}
    //put errors in ntuples
    int binlow = (int)hfLow;
    int binhigh = (int)hfHigh;
    if (binmode==1) {
      binlow = ntracksLow;
      binhigh = ntracksHigh;
    }
    ntuple1s->Fill(binlow,binhigh,yLowF,yHighF,temperrF,temperrB,RFBerr);
  }//end of rapidity bin loop
  }//end of hf bin loop

  cout << "Bbins = " << Bbins << endl;
  cout << "Fbins = " << Fbins << endl;
  cout << "RFBbins = " << RFBbins << endl;

  if (Fbins>0) temperrFMean = temperrFMean/Fbins;
  else temperrFMean = 0;
  if (Bbins>0) temperrBMean = temperrBMean/Bbins;
  else temperrBMean = 0;
  if (RFBbins>0) RFBerrMean = RFBerrMean/RFBbins;
  else RFBerrMean = 0;

  //print mean values
  cout << setw(binColWidth) << setfill(separator) << "errMean";
  cout << setw(Width1S) << setfill(separator) << temperrBMean;
  cout << setw(Width1S) << setfill(separator) << temperrFMean;
  cout << setw(Width1S) << setfill(separator) << RFBerrMean;
  cout << endl;

  //put average errors in ntuples
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

    //put errors in ntuple
    int binlow = (int)hfLow;
    int binhigh = (int)hfHigh;
    if (binmode==1) {
      binlow = ntracksLow;
      binhigh = ntracksHigh;
    }
    ntuple1savg->Fill(binlow,binhigh,yLowF,yHighF,temperrFMean,temperrBMean,RFBerrMean);
  }//end of rapidity bin loop
  }//end of hf bin loop

  //Write errors in ntuple file
  outFile.cd();
  ntuple1s->Write();
  ntuple1savg->Write();
  delete ntuple1s;
  delete ntuple1savg;

  outFile.Close();

}
