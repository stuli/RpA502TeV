#include <iomanip>
#include <sstream>
#include "CheckFitQualitySlim.C"

using namespace std;
void NewGetFixParamDeviationspt0to6to30(int whichUpsilon = 1, int param=1) {

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
  float ybins1[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  float ybins2[5] = {-1.93,-0.8,0.0,0.8,1.93};
  float ybins3[3] = {-1.93,0.0,1.93};

  float* ybinsptr;
  int numybinstemp;
  if (whichUpsilon==1) {
    ybinsptr = &ybins1[0];
    numybinstemp = sizeof(ybins1)/sizeof(float)-1;
  }
  else if (whichUpsilon==2) {
    ybinsptr = &ybins2[0];
    numybinstemp = sizeof(ybins2)/sizeof(float)-1;
  }
  else if (whichUpsilon==3) {
    ybinsptr = &ybins3[0];
    numybinstemp = sizeof(ybins3)/sizeof(float)-1;
  }

  const int numybins = numybinstemp;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -10;
  int max = 10;

  cout << whichUpsilon << "S " << "RAPIDITY BINS IN LOW AND HIGH PT FOR "  << whichParam.Data() << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPFileName;
  TFile* theFile;

  TString outFileName = Form("FixParamErrors/Param%sDev%is_pt0to6to30.root",whichParam.Data(),whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1sy = new TNtuple("ntuple1sy","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);
  TNtuple* ntuple1syavg = new TNtuple("ntuple1syavg","Error estimates in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

  float temperrMeany = 0;
  float PPerrMeany = 0;
  float RpAerrMeany = 0;
  int PAAltGood = 0;
  int PPAltGood = 0;
  int PAybins = 0;
  int PPybins = 0;
  int RpAybins = 0;

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

    yLowCM = *(ybinsptr+iy);
    yHighCM = *(ybinsptr+iy+1);
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
    if (whichParam=="f") PPtoo = kFALSE;
    if (whichParam=="MCf") PPtoo = kFALSE;

    //print bin label
    stringstream stream1;
    stream1 << fixed << setprecision(2) << binLow;
    string strbinLow = stream1.str();
    stringstream stream2;
    stream2 << fixed << setprecision(2) << binHigh;
    string strbinHigh = stream2.str();
    binLabel = strbinLow+"<y<"+strbinHigh;
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
    if (PAAltGood) {PAybins++; temperrMeany = temperrMeany + temperr;}
    if (PPAltGood) {PPybins++; PPerrMeany = PPerrMeany + PPerr;}
    if (PAAltGood && PPAltGood) {RpAybins++; RpAerrMeany = RpAerrMeany + RpAerr;}
    //put errors in ntuple
    ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
  }// end of y loop
  }// end of pt loop

  cout << "PPybins = " << PPybins << endl;
  cout << "PAybins = " << PAybins << endl;
  cout << "RpAybins = " << RpAybins << endl;

  if (PAybins>0) temperrMeany = temperrMeany/PAybins;
  else temperrMeany = 0;
  if (PPybins>0) PPerrMeany = PPerrMeany/PPybins;
  else PPerrMeany = 0;
  if (RpAybins>0) RpAerrMeany = RpAerrMeany/RpAybins;
  else RpAerrMeany = 0;

  //print mean values
  cout << setw(binColWidth) << setfill(separator) << "errMeany";
  cout << setw(Width1S) << setfill(separator) << PPerrMeany;
  cout << setw(Width1S) << setfill(separator) << temperrMeany;
  cout << setw(Width1S) << setfill(separator) << RpAerrMeany;
  cout << endl;

  //put errors in ntuples
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

    yLowCM = *(ybinsptr+iy);
    yHighCM = *(ybinsptr+iy+1);
    binLow = yLowCM;
    binHigh = yHighCM;
    binvar = "y";
    //put errors in ntuple
    ntuple1syavg->Fill(binLow,binHigh,PPerrMeany,temperrMeany,RpAerrMeany);
  }// end of y loop
  }// end of pt loop

  //Write errors in ntuple file
  outFile.cd();
  ntuple1sy->Write();
  ntuple1syavg->Write();
  delete ntuple1sy;
  delete ntuple1syavg;

  outFile.Close();


}
