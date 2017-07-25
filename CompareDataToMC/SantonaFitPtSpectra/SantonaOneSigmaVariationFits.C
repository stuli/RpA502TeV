#include "../../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../../cutsAndBin.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
//valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, 	float dphiEp2Low=0,  float dphiEp2High=0) ;


Double_t alt3(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  //  return fpar[0] +  exp( -(xx-fpar[1])/fpar[2] ) * (xx *fpar[3]);
  return  (  ( fpar[0] + fpar[1]*xx ) * ( 1 + fpar[2] / (xx-fpar[3]) ) );
}
Double_t alt4(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  //  return fpar[0] +  exp( -(xx-fpar[1])/fpar[2] ) * (xx *fpar[3]);
  return  (  fpar[0]*exp( -xx / fpar[1] ) - 1 / (xx-fpar[2]) + fpar[3] ) ;
}
Double_t alt5(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  //  return fpar[0] +  exp( -(xx-fpar[1])/fpar[2] ) * (xx *fpar[3]);
  return  (   fpar[0]*exp( -xx / fpar[1] )  +   fpar[4] * ( TMath::Erf( (xx -fpar[2]) / fpar[3] ) + 1 ) );
}
Double_t alt6(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  //  return fpar[0] +  exp( -(xx-fpar[1])/fpar[2] ) * (xx *fpar[3]);
  return  (   fpar[0]*exp( -xx / fpar[1] )  +   fpar[4] * ( TMath::Erf( (xx -fpar[2]) / fpar[3] ) + 1 ) + fpar[5]);
}

Double_t polyy(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  return fpar[0] + xx*fpar[1] + xx*xx*fpar[2] + xx*xx*xx*fpar[3] + xx*xx*xx*xx*fpar[4];// + xx*xx*xx*xx*xx*fpar[5];
}
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

void SantonaOneSigmaVariationFits(int state = 1, int collId= kAADATA) {
//  cout << "Run getMcSepctra.C before running this macro" << endl;
  
  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  int nPtBins=0;
  double* ptBin;

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
  }
  
//  TH1D* hptMc;
//  TH1D* hptSig;
  TH1D* hRatio;   // final Ratio w/ efficiency correctdion
  TF1* funct_write;
//  TH1D* hRatioTrue;
  TF1* old_funct;

  TString fcollId;
  if(collId == kPPDATA) fcollId = "PP";
  else if(collId == kAADATA) fcollId = "AA";


  //  TFile* inf = new TFile(Form("../efficiency/efficiency_ups%ds_MC_noWeight.root",state));
  TFile* inf = new TFile(Form("FitResultsratioDataMC_%s_DATA_%dsState.root",fcollId.Data(),state), "Open");
  hRatio = (TH1D*)inf->Get("hRatio");
  funct_write = (TF1*)inf->Get("funct_write");
  old_funct = (TF1*)inf->Get("funct");

  TF1 *func = (TF1*)hRatio->GetFunction("REM");

  cout << "BLOOP" << func->GetParameter(1) << endl;

  TCanvas* c1 =  new TCanvas("c1","",800,800);
  c1->cd();
  TH1ScaleByWidth(hRatio);
  hRatio ->SetAxisRange(0,2,"Y");
  if ( state == 2 ) hRatio ->SetAxisRange(0,5,"Y");
  handsomeTH1(hRatio,1);

  // Fit :
  //TF1 *fVar1pp = new TF1("f1srpp",fTsallis1SR,0,30,4);

  TF1* funct;
  float NewValues[4]={0,0,0,0};

/*  for ( int i=0 ; i<4 ; i++)
  {
  NewValues[i] = funct_write->GetParameter(i) + funct_write->GetParError(i);
  cout << NewValues[i] << " , " << endl;
  }
// */
 

  if(state==1) funct = new TF1("dataMcRatio",fTsallis1SR, 0,30,4);
  else if(state==2) funct = new TF1("dataMcRatio",fTsallis2SR, 0,30,4);
  funct->SetParameters(NewValues[0],NewValues[1],NewValues[2],NewValues[3]);
  funct->SetParLimits(0,-1.,1.3);
  funct->SetParLimits(1, 0.001,3.1);
  funct->SetParLimits(2,0.12,6);
  funct->SetParLimits(3,-1.,20);



  hRatio->Fit ( funct, "REM");
  hRatio->Draw();
  jumSun(0,1,30,1);

  int nRatioBin = 100;
  TH1D* hRatioUp = new TH1D("hRatioUp","",nRatioBin,0,30);
  TH1D* hRatioDown = new TH1D("hRatioDown","",nRatioBin,0,30);
  float variationCenter = 0;
  float upRatio = 0.3;
  if ( (state == 2) && (collId==kAADATA)) 
    { variationCenter = 15;  upRatio = 0.5;}
  if ( (state == 2) && (collId==kPPDATA))  
    { variationCenter = 15; upRatio = 0.1;}
  float downRatio =  - upRatio;
  for ( int ii=1; ii<= nRatioBin ;  ii++) {
    float xx = hRatioUp->GetBinCenter(ii);
    hRatioUp->SetBinContent(ii,   funct->Eval(xx) * (  (upRatio/15.) * (xx-variationCenter) +1 ) );
  }
  handsomeTH1(hRatioUp, 4);  
for ( int ii=1; ii<= nRatioBin ;  ii++) {
    float xx = hRatioDown->GetBinCenter(ii);
    hRatioDown->SetBinContent(ii,   funct->Eval(xx) * (  (downRatio/15.) * (xx-variationCenter) +1 ) );
  }
  handsomeTH1(hRatioUp, 4);
  handsomeTH1(hRatioDown, 4);

  hRatioUp->Draw("same hist");
  hRatioDown->Draw("same hist");


  if(state==1)
  { 
    cout << "TF1* weightFunction = new TF1(@weightCurve@,@(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2]))));" << endl;   cout << " weightFunction->SetParameters( " ;
    for ( int i=0 ; i<4 ; i++) 
    { 
      if ( i!=3 ) 
        cout << funct->GetParameter(i) << ", " ;
      else   cout << funct->GetParameter(i) << "); "<< endl ;
    }
  }
  else if(state==2)
  {
    //cout << "TF1* weightFunction = new TF1(@weightCurve@,@[0]*x+[1]);" << endl;   cout << " weightFunction->SetParameters( " ;
    cout << "TF1* weightFunction = new TF1(@weightCurve@,@(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2]))));" << endl;   cout << " weightFunction->SetParameters( " ;
    for ( int i=0 ; i<4 ; i++) 
    { 
      if ( i!=3 ) 
        cout << funct->GetParameter(i) << ", " ;
      else   cout << funct->GetParameter(i) << "); "<< endl ;
    }
  }

  //  cout << " This macro MUST be ran on ROOT5!!!! " << endl;



  c1->SaveAs(Form("VarieddNdpT_dataMC_%s_DATA_%dsState.png",fcollId.Data(),state));
  TFile* outf = new TFile(Form("VariedFitResultsratioDataMC_%s_%dsState.root",getCollID(collId).Data(), state),"recreate");
  hRatio->Write();
//  hptSig->SetName("hData");
//  hptMc->SetName("hMC");
//  hptSig->Write();
//  hptMc->Write();
  funct->Write();

  outf->Write();

  inf->Close();

}

/*
valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
    float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString SignalCB = "Double";
  //  TFile* inf = new TFile(Form("../fitResults/ptDependence/PAS_fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TFile* inf = new TFile(Form("../TEST_newNom/PAS_fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

// */
