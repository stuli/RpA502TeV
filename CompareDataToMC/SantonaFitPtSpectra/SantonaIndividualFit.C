#include "../../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../cutsAndBin.h"
using namespace std;


// Ratio fitting functions

/*
Double_t fTsallis1SR(Double_t *x, Double_t *fpar)
{

  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y1S));
  Double_t mT = TMath::Sqrt(pdgMass.Y1S*pdgMass.Y1S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y1S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = c*xx*pow;


  Double_t c1 = (fpar[2]-1)*(fpar[2]-2)/(fpar[2]*fpar[3]*(fpar[2]*fpar[3]+(fpar[2]-2)*pdgMass.Y1S));
  Double_t mT1 = TMath::Sqrt(pdgMass.Y1S*pdgMass.Y1S+xx*xx);
  Double_t pow1 = TMath::Power((1+(mT1-pdgMass.Y1S)/(fpar[2]*fpar[3])),-fpar[2]);

  Double_t f1 = c1*xx*pow1;

  Double_t fr = f/f1;

                   
  return fr;
}

Double_t fTsallis2SR(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y2S));
  Double_t mT = TMath::Sqrt(pdgMass.Y2S*pdgMass.Y2S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y2S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = c*xx*pow;


  Double_t c1 = (fpar[2]-1)*(fpar[2]-2)/(fpar[2]*fpar[3]*(fpar[2]*fpar[3]+(fpar[2]-2)*pdgMass.Y2S));
  Double_t mT1 = TMath::Sqrt(pdgMass.Y2S*pdgMass.Y2S+xx*xx);
  Double_t pow1 = TMath::Power((1+(mT1-pdgMass.Y2S)/(fpar[2]*fpar[3])),-fpar[2]);

  Double_t f1 = c1*xx*pow1;

  Double_t fr = f/f1;

                   
  return fr;
}
// */
//

// Individual fitting functions: 
// Tsallis
/*
Double_t fTsallis1S(Double_t *x, Double_t *fpar)
{

  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y1S));
  Double_t mT = TMath::Sqrt(pdgMass.Y1S*pdgMass.Y1S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y1S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = fpar[2]*c*xx*pow;
//  Double_t f = c*xx*pow;

  return f;
}
// */

Double_t fTsallis2S(Double_t *x, Double_t *fpar)
{

  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y2S));
  Double_t mT = TMath::Sqrt(pdgMass.Y2S*pdgMass.Y2S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y2S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = fpar[2]*c*xx*pow;
//  Double_t f = c*xx*pow;

  return f;
}

// Exponential * pT:

Double_t fTsallis1S(Double_t *x, Double_t *fpar)
{

  Float_t xx = x[0];
  Double_t exptail = TMath::Exp(-xx/fpar[1]);
  Double_t f = fpar[0]*xx*exptail;

  return f;
  }






void SantonaIndividualFit(int state = 1, int collId= kAADATA) {
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1); 
  TH1::SetDefaultSumw2();
  
  TH1D* hptMc = NULL;
  TH1D* hptSig = NULL;
//  TH1D* hRatio;   // final Ratio w/ efficiency correctdion

  TH1::SetDefaultSumw2();
  
  TString fcollId;
  if(collId == kPPDATA) fcollId = "PP";
  else if(collId == kAADATA) fcollId = "AA";


  //  TFile* inf = new TFile(Form("../efficiency/efficiency_ups%ds_MC_noWeight.root",state));
  TFile* inf = new TFile(Form("../ratioDataMC_%s_DATA_%dsState.root",fcollId.Data(),state));
  //hRatio = (TH1D*)inf->Get("hRatio");
  hptMc  = (TH1D*)inf->Get("hMC");
  hptSig = (TH1D*)inf->Get("hData"); 

//  TH1D* haccBinwise = new TH1D( "haccBinwise", "", hptMc->GetNbinsX(), 0, 30);
//  TH1D* heffBinwise = new TH1D( "heffBinwise", "", hptMc->GetNbinsX(), 0, 30);
//  TH1::SetDefaultSumw2();

  TH1D* haccBinwise;
  TH1D* heffBinwise;
  haccBinwise = (TH1D*)  hptMc->Clone("haccBinwise");
  haccBinwise->Reset();
  heffBinwise = (TH1D*)  hptMc->Clone("heffBinwise");
  heffBinwise->Reset();

  // Acceptance correction : 
  TFile* accf = new TFile(Form("../../acceptance/acceptance_wgt_final_20170106.root",state));
  TH1D* hacc =  (TH1D*)accf->Get(Form("hptAcc%s%dS",fcollId.Data(),state));
  for ( int ii =1 ; ii<= hptMc->GetNbinsX() ; ii++ ) {
  	haccBinwise->SetBinContent( ii, hacc->Interpolate(hptMc->GetBinCenter(ii)) ); 
  }

  // Efficiency Correction:
  TFile* efff = new TFile(Form("../../efficiency/efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state));
  TH1D* heff =  (TH1D*)efff->Get(Form("hptEff%s",fcollId.Data()));
  for ( int ii =1 ; ii<= hptMc->GetNbinsX() ; ii++ ) {
        heffBinwise->SetBinContent( ii, heff->Interpolate(hptMc->GetBinCenter(ii)) );
  }
// */


  TCanvas* c1 =  new TCanvas("c1","",800,800);
  c1->Divide(1,2);
  c1->cd(1);
  hptMc->Divide(haccBinwise); // acceptance correction
  hptMc->Divide(heffBinwise); // efficiency correction
  TH1ScaleByWidth(hptMc);
  scaleInt(hptMc);

  handsomeTH1(hptMc,1);
  hptMc->SetAxisRange(0,0.95,"Y");
  if ( state == 2 ) hptMc ->SetAxisRange(0,1.50,"Y");
  hptMc->Draw("hist");


  // Fit :

  TF1* functMC;
  TF1* functData;

  int n = 2; // change number of parameters
  if(state==1) functMC = new TF1("MC",fTsallis1S, 0.0,30.0,n);
  else if(state==2) functMC = new TF1("MC",fTsallis2S, 0,15.0,n);
  if(state==1) functData = new TF1("data",fTsallis1S, 0.0,30.0,n);
  else if(state==2) functData = new TF1("data",fTsallis2S, 0,15.0,n);

/*  funct->SetParameters(0.06123,1.023,2.123,1);
  funct->SetParLimits(0,-1.,1.3);
  funct->SetParLimits(1, 0.001,3.1);
  funct->SetParLimits(2,0.12,6);
  funct->SetParLimits(3,-1.,20);
// */
//  functMC->SetParameters(1.8,0.2,0.02);
  functMC->SetParameters(0.4,2.3);
  functMC->SetParLimits(0,0,5);
  functMC->SetParLimits(1,0,5);
//  functMC->SetParLimits(2,0,100);
 
//  functData->SetParameters(1.8,0.5,0.02);
  functData->SetParameters(0.4,2.3);
  functData->SetParLimits(0,0,5);
  functData->SetParLimits(1,0,5);
//  functData->SetParLimits(2,0,2);
// 0-2 is a good range for data

  //hptMc->Fit ( functMC, "REM");
  hptMc->Fit ( functMC,"IR");
  hptMc->Draw();
  jumSun(0,1,30,1);

  c1->cd(2);
  hptSig->Divide(haccBinwise); // acceptance correction
  hptSig->Divide(heffBinwise); // efficiency correction
  TH1ScaleByWidth(hptSig);
  scaleInt(hptSig);

  handsomeTH1(hptSig,1);
  hptSig->SetAxisRange(0,0.95,"Y");
  if ( state == 2 ) hptSig ->SetAxisRange(0,1.50,"Y");
  hptSig->Draw("hist");

  hptSig->Fit ( functData,"IR");
  hptSig->Draw();
  jumSun(0,1,30,1);



// Writing all the parameters to a text file only for 2 parameters:
/*
  double onesigmaparametersData[2] = {0,0};
  double parametersData[2] = {0,0};
  double onesigmaparametersMC[2] = {0,0};
  double parametersMC[2] = {0,0};


  ofstream parameterfile;
  parameterfile.open (Form("Individualparameterfile%s_%dsState.txt",fcollId.Data(),state));

  for ( int i=0 ; i<2 ; i++)
  {
	parametersData[i] = functData->GetParameter(i);
  	parameterfile << parametersData[i] << "\n" ;
  }

  for ( int i=0 ; i<2 ; i++)
  {
        onesigmaparametersData[i] = functData->GetParameter(i)+functData->GetParError(i);
        parameterfile << onesigmaparametersData[i] << "\n" ;
  }
  
  parameterfile << "MC: \n" ;

  for ( int i=0 ; i<2 ; i++)
  {
        parametersMC[i] = functMC->GetParameter(i);
        parameterfile << parametersMC[i] << "\n" ;
  }

  for ( int i=0 ; i<2 ; i++)
  {
        onesigmaparametersMC[i] = functMC->GetParameter(i)+functMC->GetParError(i);
        parameterfile << onesigmaparametersMC[i] << "\n" ;
  }


  parameterfile.close();
// */

//  TF1 *func = (TF1*)hRatio->GetFunction("REM");

  c1->SaveAs(Form("IndividualdNdpT_dataMC_%s_DATA_%ds_Exp.png",fcollId.Data(),state));
  TFile* outf = new TFile(Form("IndividualFitResultsratioDataMC_%s_%ds_Exp.root",getCollID(collId).Data(), state),"recreate");
//  hRatio->Write();
  hptSig->Write();
  hptMc->Write();

  functData->Write();
  functMC->Write();

  outf->Write();

}

