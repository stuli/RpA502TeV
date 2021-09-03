#include <iostream>
#include "../../HeaderFiles/rootFitHeaders.h"
#include "../../HeaderFiles/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../HeaderFiles/cutsAndBin.h"
#include "../../HeaderFiles/PsetCollection.h"
#include "../../HeaderFiles/CMS_lumi.C"
#include "../../HeaderFiles/tdrstyle.C"
#include "../../HeaderFiles/StyleSetting.h"

void CombineChisqHistos(int collId=kPPDATA, int whichUpsilon=1, int whichn=2) {


  TString ptHistoName = "hchisqpt";
  TString yHistoName = "hchisqy";

  //Extract histograms
  TString HistoFileNamept = Form("ChisqHisto_PP_%is_fixedn%ipt.root",whichUpsilon,whichn);
  cout << HistoFileNamept << endl;
  HistoFilept = TFile::Open(HistoFileNamept,"READ");
  TH1F* hchisqpt = (TH1F*)HistoFilept->Get(ptHistoName);

  TString HistoFileNamey = Form("ChisqHisto_PP_%is_fixedn%iy.root",whichUpsilon,whichn);
  cout << HistoFileNamey << endl;
  HistoFiley = TFile::Open(HistoFileNamey,"READ");
  TH1F* hchisqy = (TH1F*)HistoFiley->Get(yHistoName);

  TString outfilename = Form("ChisqHisto_PP_%is_fixedn%i.root",whichUpsilon,whichn);
  TFile outFile(outfilename, "RECREATE");

  //save histograms
  outFile.cd();
  hchisqpt->Write();
  hchisqy->Write();
  outFile.Close();

}
