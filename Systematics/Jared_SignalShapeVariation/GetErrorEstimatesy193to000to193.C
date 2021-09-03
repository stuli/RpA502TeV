#include <iostream>
#include <iomanip>
#include <sstream>
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

using namespace std;
void GetErrorEstimatesy193to000to193(int whichUpsilon = 1) {

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
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;

  const char separator    = ' ';
  const int binColWidth     = 14;
  const int Width1S      = 12;
  int min = -10;
  int max = 10;

  TCanvas* cntuple =  new TCanvas("cntuple","ntupleresults",4,545,1100,400);
  cntuple->Divide(3,1);
  TCanvas *c1 = new TCanvas("c1","c1",4,45,400,400);

  cout << setw(binColWidth) << setfill(separator) << "     BIN     ";
  cout << setw(Width1S) << setfill(separator) << Form("PP %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << Form("PA %iS ERR",whichUpsilon);
  cout << setw(Width1S) << setfill(separator) << "RpA ERR";
  if (whichUpsilon>1) {
    cout << setw(Width1S) << setfill(separator) << Form("PP R%i1 ERR",whichUpsilon);
    cout << setw(Width1S) << setfill(separator) << Form("PA R%i1 ERR",whichUpsilon);
    cout << setw(Width1S) << setfill(separator) << "DR ERR";
  }
  cout << endl;

  float yLowCM, yHighCM, yLowPP, yHighPP, ptLow, ptHigh, binLow, binHigh;
  string binvar;
  TString binLabel, kineLabel, histFileName, PPFileName, canvasFileName;
  TFile* theFile;

  TString outFileName = Form("ErrorEstimates/SystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  TFile outFile(outFileName, "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);

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

    ptLow = ptbins[ipt];
    ptHigh = ptbins[ipt+1];
    binLow = ptLow;
    binHigh = ptHigh;
    binvar = "pt";

    if (yLowCM<-2.5) {
      PPtoo = kFALSE;
      cout << "No PP here." << endl;
    }
    else PPtoo = kTRUE;

    //print bin label
    binLabel = Form("pt%.2f-%.2f_y%.2f-%.2f",ptLow,ptHigh,yLowCM,yHighCM);
    cout << setw(binColWidth) << setfill(separator) << binLabel;

    //import results of pseudo-experiments
    kineLabel = getKineLabel (collId, ptLow, ptHigh, yLowCM, yHighCM, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
    histFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpResults_%s.root",kineLabel.Data());
    theFile = new TFile(histFileName);
    TNtuple* ntupleResults = (TNtuple*)theFile->Get("ntuple;1");

    if (PPtoo) {
      kineLabel = getKineLabel (collIdPP, ptLow, ptHigh, yLowPP, yHighPP, muPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
      PPFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpResults_%s.root",kineLabel.Data());
      theFile = new TFile(PPFileName);
      TNtuple* PPntupleResults = (TNtuple*)theFile->Get("ntuple;1");

      //extract error values
      if (PPntupleResults->GetMinimum("diff1s")<min) min = -20;
      else min = -10;
      if (PPntupleResults->GetMaximum("diff1s")>max) max = 20;
      else max = 10;
      cntuple->cd(1);
      TH1F * ntuplehistoPP = new TH1F("ntuplehistoPP", "pp %Diff in Yield", 100,min,max);
      ntuplehistoPP->SetXTitle("%Diff");
      ntuplehistoPP->GetXaxis()->SetTitleSize(0.05);
      ntuplehistoPP->GetYaxis()->SetLabelSize(0.05);
      ntuplehistoPP->GetXaxis()->SetLabelSize(0.05);
      ntuplehistoPP->GetXaxis()->SetRangeUser(min,max);
      if (whichUpsilon==1) PPntupleResults->Draw("diff1s>>ntuplehistoPP");
      else if (whichUpsilon==2) PPntupleResults->Draw("diff2s>>ntuplehistoPP");
      else if (whichUpsilon==3) PPntupleResults->Draw("diff3s>>ntuplehistoPP");
      float PPerr = ntuplehistoPP->GetMean();
      float PPrms = ntuplehistoPP->GetRMS();
      if (TMath::Abs(PPrms)>TMath::Abs(PPerr)) PPerr = PPrms;
      //ntuplehistoPP->Fit("gaus","Q");
      //float PPmu = ntuplehistoPP->GetFunction("gaus")->GetParameter(1);//fitted mean
      //if (TMath::Abs(PPmu)>TMath::Abs(PPerr)) PPerr = PPmu;
      //float PPsigma = ntuplehistoPP->GetFunction("gaus")->GetParameter(2);//fitted sigma
      //if (TMath::Abs(PPsigma)>TMath::Abs(PPerr)) PPerr = PPsigma;
    }

    if (ntupleResults->GetMinimum("diff1s")<min) min = -20;
    if (ntupleResults->GetMaximum("diff1s")>max) max = 20;
    cntuple->cd(2);
    TH1F * ntuplehisto = new TH1F("ntuplehisto", "pPb %Diff in Yield", 100,min,max);
    ntuplehisto->SetXTitle("%Diff");
    ntuplehisto->GetXaxis()->SetTitleSize(0.05);
    ntuplehisto->GetYaxis()->SetLabelSize(0.05);
    ntuplehisto->GetXaxis()->SetLabelSize(0.05);
    ntuplehisto->GetXaxis()->SetRangeUser(min,max);
    if (whichUpsilon==1) ntupleResults->Draw("diff1s>>ntuplehisto");
    else if (whichUpsilon==2) ntupleResults->Draw("diff2s>>ntuplehisto");
    else if (whichUpsilon==3) ntupleResults->Draw("diff3s>>ntuplehisto");
    float temperr = ntuplehisto->GetMean();  
    float temprms = ntuplehisto->GetRMS();
    if (TMath::Abs(temprms)>TMath::Abs(temperr)) temperr = temprms;
    //ntuplehisto->Fit("gaus","Q");
    //float tempmu = ntuplehisto->GetFunction("gaus")->GetParameter(1);//fitted mean
    //if (TMath::Abs(tempmu)>TMath::Abs(temperr)) temperr = tempmu;
    //float tempsigma = ntuplehisto->GetFunction("gaus")->GetParameter(2);//fitted sigma
    //if (TMath::Abs(tempsigma)>TMath::Abs(temperr)) temperr = tempsigma;

    //Calculate single and double ratios
    TLeaf *yield1sNomLeaf = ntupleResults->GetLeaf("yield1sNom");
    TLeaf *yield1sAltLeaf = ntupleResults->GetLeaf("yield1sAlt");
    if (whichUpsilon==2) {
      TLeaf *yield2sNomLeaf = ntupleResults->GetLeaf("yield2sNom");
      TLeaf *yield2sAltLeaf = ntupleResults->GetLeaf("yield2sAlt");
    }
    else if (whichUpsilon==3) {
      TLeaf *yield2sNomLeaf = ntupleResults->GetLeaf("yield3sNom");
      TLeaf *yield2sAltLeaf = ntupleResults->GetLeaf("yield3sAlt");
    }
  if (PPtoo) {
    TLeaf *PPyield1sNomLeaf = PPntupleResults->GetLeaf("yield1sNom");
    TLeaf *PPyield1sAltLeaf = PPntupleResults->GetLeaf("yield1sAlt");
    if (whichUpsilon==2) {
      TLeaf *PPyield2sNomLeaf = PPntupleResults->GetLeaf("yield2sNom");
      TLeaf *PPyield2sAltLeaf = PPntupleResults->GetLeaf("yield2sAlt");
    }
    else if (whichUpsilon==3) {
      TLeaf *PPyield2sNomLeaf = PPntupleResults->GetLeaf("yield3sNom");
      TLeaf *PPyield2sAltLeaf = PPntupleResults->GetLeaf("yield3sAlt");
    }
    TH1F* hratioPP = new TH1F("hratioPP","hratioPP",100,min,max);
    TH1F* hRpA = new TH1F("hRpA","hRpA",100,min,max);
    TH1F* hDoubleRatio = new TH1F("hDoubleRatio","hDoubleRatio",100,min,max);
  }
    const int NEvents = 100;
    TH1F* hratio = new TH1F("hratio","hratio",100,min,max);
    for (int i = 0; i<NEvents; i++) {
      ntupleResults->GetEntry(i);
      float yield1sNom = (float)yield1sNomLeaf->GetValue();
      float yield1sAlt = (float)yield1sAltLeaf->GetValue();
      if (whichUpsilon>1) {
        float yield2sNom = (float)yield2sNomLeaf->GetValue();
        float yield2sAlt = (float)yield2sAltLeaf->GetValue();
        //PA R21
        float r21Nom = yield2sNom/yield1sNom;
        float r21Alt = yield2sAlt/yield1sAlt;
        float r21Diff = (r21Alt-r21Nom)/r21Nom*100;
        hratio->Fill(r21Diff);
      }
      if (PPtoo) {
        PPntupleResults->GetEntry(i);
        float PPyield1sNom = (float)PPyield1sNomLeaf->GetValue();
        float PPyield1sAlt = (float)PPyield1sAltLeaf->GetValue();
        //RpA
        float RpANom = yield1sNom/PPyield1sNom;
        float RpAAlt = yield1sAlt/PPyield1sAlt;
        if (whichUpsilon>1) {
          float PPyield2sNom = (float)PPyield2sNomLeaf->GetValue();
          float PPyield2sAlt = (float)PPyield2sAltLeaf->GetValue();
          RpANom = yield2sNom/PPyield2sNom;
          RpAAlt = yield2sAlt/PPyield2sAlt;
          //PP R21
          float r21NomPP = PPyield2sNom/PPyield1sNom;
          float r21AltPP = PPyield2sAlt/PPyield1sAlt;
          float r21DiffPP = (r21AltPP-r21NomPP)/r21NomPP*100;
          hratioPP->Fill(r21DiffPP);
          //double ratio
          float DR21Nom = r21Nom/r21NomPP;
          float DR21Alt = r21Alt/r21AltPP;
          float DR21Diff = (DR21Alt-DR21Nom)/DR21Nom*100;
          hDoubleRatio->Fill(DR21Diff);
        }
        float RpADiff = (RpAAlt-RpANom)/RpANom*100;
        hRpA->Fill(RpADiff);
      }
    }

    if (whichUpsilon>1) {
      c1->cd();
      hratio->Draw();
      float ratioerr = hratio->GetMean();
      float ratiorms = hratio->GetRMS();
      if (TMath::Abs(ratiorms)>TMath::Abs(ratioerr)) ratioerr = ratiorms;
      //hratio->Fit("gaus","Q");
      //float ratiomu = hratio->GetFunction("gaus")->GetParameter(1);//fitted mean
      //if (TMath::Abs(ratiomu)>TMath::Abs(ratioerr)) ratioerr = ratiomu;
      //float ratiosigma = hratio->GetFunction("gaus")->GetParameter(2);//fitted sigma
      //if (TMath::Abs(ratiosigma)>TMath::Abs(ratioerr)) ratioerr = ratiosigma;

      if (PPtoo) {
        hratioPP->Draw();
        float ratioerrPP = hratioPP->GetMean();
        float ratiormsPP = hratioPP->GetRMS();
        if (TMath::Abs(ratiormsPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiormsPP;
        //hratioPP->Fit("gaus","Q");
        //float ratiomuPP = hratioPP->GetFunction("gaus")->GetParameter(1);//fitted mean
        //if (TMath::Abs(ratiomuPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiomuPP;
        //float ratiosigmaPP = hratioPP->GetFunction("gaus")->GetParameter(2);//fitted sigma
        //if (TMath::Abs(ratiosigmaPP)>TMath::Abs(ratioerrPP)) ratioerrPP = ratiosigmaPP;

        hDoubleRatio->Draw();
        float DRerr = hDoubleRatio->GetMean();
        float DRrms = hDoubleRatio->GetRMS();
        if (TMath::Abs(DRrms)>TMath::Abs(DRerr)) DRerr = DRrms;
        //hDoubleRatio->Fit("gaus","Q");
        //float DRmu = hDoubleRatio->GetFunction("gaus")->GetParameter(1);//fitted mean
        //if (TMath::Abs(DRmu)>TMath::Abs(DRerr)) DRerr = DRmu;
        //float DRsigma = hDoubleRatio->GetFunction("gaus")->GetParameter(2);//fitted sigma
        //if (TMath::Abs(DRsigma)>TMath::Abs(DRerr)) DRerr = DRsigma;
      }
    }

  if (PPtoo) {
    cntuple->cd(3);
    hRpA->SetXTitle("%Diff in RpA");
    hRpA->GetXaxis()->SetTitleSize(0.05);
    hRpA->GetYaxis()->SetLabelSize(0.05);
    hRpA->GetXaxis()->SetLabelSize(0.05);
    hRpA->GetXaxis()->SetRangeUser(min,max);
    hRpA->Draw();
    float RpAerr = hRpA->GetMean();
    float RpArms = hRpA->GetRMS();
    if (TMath::Abs(RpArms)>TMath::Abs(RpAerr)) RpAerr = RpArms;
    //hRpA->Fit("gaus","Q");
    //float RpAmu = hRpA->GetFunction("gaus")->GetParameter(1);//fitted mean
    //if (TMath::Abs(RpAmu)>TMath::Abs(RpAerr)) RpAerr = RpAmu;
    //float RpAsigma = hRpA->GetFunction("gaus")->GetParameter(2);//fitted sigma
    //if (TMath::Abs(RpAsigma)>TMath::Abs(RpAerr)) RpAerr = RpAsigma;
  }
  else {
    PPerr = 0;
    RpAerr = 0;
    ratioerrPP = 0;
    DRerr = 0;
  }

    //Save the plots of yields and RpA errors
    canvasFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f.png",whichUpsilon,ptLow,ptHigh,yLowCM,yHighCM,muPtCut);
    cntuple->SaveAs(canvasFileName);
    canvasFileName = Form("PseudoExperimentResultsFeb2018/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt%.1f.pdf",whichUpsilon,ptLow,ptHigh,yLowCM,yHighCM,muPtCut);
    cntuple->SaveAs(canvasFileName);

    //Take absolute values
    PPerr = TMath::Abs(PPerr);
    temperr = TMath::Abs(temperr);
    RpAerr = TMath::Abs(RpAerr);
    ratioerrPP = TMath::Abs(ratioerrPP);
    ratioerr = TMath::Abs(ratioerr);
    DRerr = TMath::Abs(DRerr);

    //print error estimate
    cout << setw(Width1S) << setfill(separator) << PPerr;
    cout << setw(Width1S) << setfill(separator) << temperr;
    cout << setw(Width1S) << setfill(separator) << RpAerr;
    if (whichUpsilon>1) {
      cout << setw(Width1S) << setfill(separator) << ratioerrPP;
      cout << setw(Width1S) << setfill(separator) << ratioerr;
      cout << setw(Width1S) << setfill(separator) << DRerr;
    }
    cout << endl;

    //put errors in ntuple
    if (whichUpsilon==1)
      if (ipt<numptbins)
        ntuple1spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
      else
        ntuple1sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr);
    else if (whichUpsilon==2)
      if (ipt<numptbins)
        ntuple2spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr,ratioerrPP,ratioerr,DRerr);
      else
        ntuple2sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr,ratioerrPP,ratioerr,DRerr);
    else if (whichUpsilon==3)
      if (ipt<numptbins)
        ntuple3spt->Fill(binLow,binHigh,PPerr,temperr,RpAerr,ratioerrPP,ratioerr,DRerr);
      else
        ntuple3sy->Fill(binLow,binHigh,PPerr,temperr,RpAerr,ratioerrPP,ratioerr,DRerr);

    theFile->Close();
  }// end of pt loop
  }// end of y loop

  c1->Close();
  cntuple->Close();

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
  outFile.Close();


}
