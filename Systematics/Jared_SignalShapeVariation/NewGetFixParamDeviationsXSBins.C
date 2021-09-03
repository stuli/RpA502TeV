#include <iomanip>
#include <sstream>
#include "CheckFitQualitySlim.C"

using namespace std;
void NewGetFixParamDeviationsXSBins(int whichUpsilon = 1, int param=1) {

  TString whichParam;
  if (param==1) whichParam = "alpha";
  else if (param==2) whichParam = "f";
  else if (param==3) whichParam = "n";
  else if (param==4) whichParam = "x";

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

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  float min = -10;
  float max = 10;

  cout << whichUpsilon << "S " << "PT BINS IN FULL Y RANGE FOR " << whichParam.Data() << endl;
  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << endl;

  float yLowCM, yHighCM, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName;
  TFile* theFile;

  TString outFileName = Form("FixParamErrors/Param%sDev%isXSBins.root",whichParam.Data(),whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pPb1sErr",numptbins);
  TNtuple* ntuple1sptavg = new TNtuple("ntuple1sptavg","Error estimates in pt bins","binlowlim:binuplim:pPb1sErr",numptbins);

  float temperrInt = 0;
  float temperrMeanpt = 0;
  int PAAltGood = 0;
  int PAptbins = 0;

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

    //Calculate errors
    float temperr = TMath::Abs((PAAltYield-PANomYield)/PANomYield*100);
    float PPerr = 0;
    float RpAerr = 0;

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << endl;

    //Sum up errors to get mean values
    if (ipt<0) {
      temperrInt = temperr;
    }
    else {
      if (PAAltGood) {PAptbins++; temperrMeanpt = temperrMeanpt + temperr;}
    }
    //put errors in ntuple
    ntuple1spt->Fill(binLow,binHigh,temperr);
  }//end of main loop

  cout << "PAptbins = " << PAptbins << endl;

  if (PAptbins>0) temperrMeanpt = temperrMeanpt/PAptbins;
  else temperrMeanpt = 0;

  //print mean values
  cout << setw(binColWidth) << setfill(separator) << "errInt";
  cout << setw(Width1S) << setfill(separator) << temperrInt;
  cout << endl;
  cout << setw(binColWidth) << setfill(separator) << "errMeanpt";
  cout << setw(Width1S) << setfill(separator) << temperrMeanpt;
  cout << endl;

  //put errors in ntuples
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

    //put errors in ntuple
    ntuple1sptavg->Fill(binLow,binHigh,temperrMeanpt);
  }

  //Write errors in ntuple file
  outFile.cd();
  ntuple1spt->Write();
  ntuple1sptavg->Write();
  delete ntuple1spt;
  delete ntuple1sptavg;

  outFile.Close();


}
