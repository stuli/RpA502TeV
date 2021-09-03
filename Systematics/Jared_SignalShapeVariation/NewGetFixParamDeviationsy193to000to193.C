#include <iomanip>
#include <sstream>
#include "CheckFitQualitySlim.C"

using namespace std;
void NewGetFixParamDeviationsy193to000to193(int whichUpsilon = 1, int param=1) {

  TString whichParam;
  if (param==1) whichParam = "alpha";
  else if (param==2) whichParam = "f";
  else if (param==3) whichParam = "n";
  else if (param==4) whichParam = "x";

  TString nominalDir = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/NominalFits_2019_03_14/";

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

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -10;
  int max = 10;

  cout << whichUpsilon << "S " << "PT BINS IN FORWARD AND BACKWARD Y FOR " << whichParam.Data() << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPFileName;
  TFile* theFile;

  TString outFileName = Form("FixParamErrors/Param%sDev%is_y193to000to193.root",whichParam.Data(),whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);

  float temperrMeanpt = 0;
  float PPerrMeanpt = 0;
  float RpAerrMeanpt = 0;
  int PAAltGood = 0;
  int PPAltGood = 0;
  int PAptbins = 0;
  int PPptbins = 0;
  int RpAptbins = 0;

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
    if (whichParam=="f") PPtoo = kFALSE;
    if (whichParam=="MCf") PPtoo = kFALSE;

    //print bin label
    binLabel = Form("pt%.2f-%.2f_y%.2f-%.2f",ptLow,ptHigh,yLowCM,yHighCM);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import fitted yields
    //cout << "Importing fit results" << endl;
    TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString NomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),kineLabel.Data());
    TFile* NomFile = TFile::Open(NomFileName,"READ");
    RooWorkspace *Nomws = (RooWorkspace*)NomFile->Get("workspace");
    NomFile->Close("R");
    float PANomYield = Nomws->var(Form("nSig%is",whichUpsilon))->getVal();
    delete Nomws;
    delete NomFile;

    TString AltkineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    TString AltFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),AltkineLabel.Data());
    TFile* AltFile = TFile::Open(AltFileName,"READ");
    RooWorkspace *Altws = (RooWorkspace*)AltFile->Get("workspace");
    AltFile->Close("R");
    float PAAltYield = Altws->var(Form("nSig%is",whichUpsilon))->getVal();
    PAAltGood = CheckFitQualitySlim(collId,ptLow,ptHigh,yLowCM,yHighCM,Altws);
    delete Altws;
    delete AltFile;

    float PPNomYield, PPAltYield;
    if (PPtoo) {
      TString PPkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString PPNomFileName = Form("%snomfitresults_upsilon_%s.root",nominalDir.Data(),PPkineLabel.Data());
      TFile* PPNomFile = TFile::Open(PPNomFileName,"READ");
      RooWorkspace *PPNomws = (RooWorkspace*)PPNomFile->Get("workspace");
      PPNomFile->Close("R");
      PPNomYield = PPNomws->var(Form("nSig%is",whichUpsilon))->getVal();
      delete PPNomws;
      delete PPNomFile;

      TString PPAltkineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      TString PPAltFileName = Form("/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/Fits/FitsWithConstraints_change%s/nomfitresults_upsilon_%s.root",whichParam.Data(),PPAltkineLabel.Data());
      TFile* PPAltFile = TFile::Open(PPAltFileName,"READ");
      RooWorkspace *PPAltws = (RooWorkspace*)PPAltFile->Get("workspace");
      PPAltFile->Close("R");
      PPAltYield = PPAltws->var(Form("nSig%is",whichUpsilon))->getVal();
      PPAltGood = CheckFitQualitySlim(collIdPP,ptLow,ptHigh,yLowPP,yHighPP,PPAltws);
      delete PPAltws;
      delete PPAltFile;
    }

    //Calculate errors
    float temperr = TMath::Abs((PAAltYield-PANomYield)/PANomYield*100);
    float PPerr = 0;
    float RpAerr = 0;
    if (PPtoo) {
      PPerr = TMath::Abs((PPAltYield-PPNomYield)/PPNomYield*100);
      float NomRpA = PANomYield/PPNomYield;
      float AltRpA = PAAltYield/PPAltYield;
      RpAerr = TMath::Abs((AltRpA-NomRpA)/NomRpA*100);
    }

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    cout << endl;

    //Sum up errors to get mean values
    if (PAAltGood) {PAptbins++; temperrMeanpt = temperrMeanpt + temperr;}
    if (PPAltGood) {PPptbins++; PPerrMeanpt = PPerrMeanpt + PPerr;}
    if (PAAltGood && PPAltGood) {RpAptbins++; RpAerrMeanpt = RpAerrMeanpt + RpAerr;}
  }// end of pt loop
  }// end of y loop

  cout << "PPptbins = " << PPptbins << endl;
  cout << "PAptbins = " << PAptbins << endl;
  cout << "RpAptbins = " << RpAptbins << endl;

  if (PAptbins>0) temperrMeanpt = temperrMeanpt/PAptbins;
  else temperrMeanpt = 0;
  if (PPptbins>0) PPerrMeanpt = PPerrMeanpt/PPptbins;
  else PPerrMeanpt = 0;
  if (RpAptbins>0) RpAerrMeanpt = RpAerrMeanpt/RpAptbins;
  else RpAerrMeanpt = 0;

  //print mean values
  cout << setw(binColWidth) << setfill(separator) << "errMeanpt";
  cout << setw(Width1S) << setfill(separator) << PPerrMeanpt;
  cout << setw(Width1S) << setfill(separator) << temperrMeanpt;
  cout << setw(Width1S) << setfill(separator) << RpAerrMeanpt;
  cout << endl;

  //put errors in ntuples
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

    //put errors in ntuple
    if (whichUpsilon==1) ntuple1spt->Fill(binLow,binHigh,PPerrMeanpt,temperrMeanpt,RpAerrMeanpt);
    else if (whichUpsilon==2) ntuple2spt->Fill(binLow,binHigh,PPerrMeanpt,temperrMeanpt,RpAerrMeanpt);
    else if (whichUpsilon==3) ntuple3spt->Fill(binLow,binHigh,PPerrMeanpt,temperrMeanpt,RpAerrMeanpt);
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
