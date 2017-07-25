#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
using namespace std;

int kNoWeight=0;
int kPP = 1;
int kAA = 2;

int k1S = 1; 
int k2S = 2; 
int k3S = 3; 

int kNoVar = 0;
int kPtPlus = 1;
int kPtMinus = 2;
int kYPlus = 3;
int kYMinus = 4;
//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
double Acc1S_NoW[10][10] = {0.0};
double Acc2S_NoW[10][10] = {0.0};
double Acc3S_NoW[10][10] = {0.0};
double Acc1S_NoW_Err[10][10] = {0.0};
double Acc2S_NoW_Err[10][10] = {0.0};
double Acc3S_NoW_Err[10][10] = {0.0};

double Acc1S_pp[10][10] = {0.0};
double Acc2S_pp[10][10] = {0.0};
double Acc3S_pp[10][10] = {0.0};
double Acc1S_pp_Err[10][10] = {0.0};
double Acc2S_pp_Err[10][10] = {0.0};
double Acc3S_pp_Err[10][10] = {0.0};

double Acc1S_AA[10][10] = {0.0};
double Acc2S_AA[10][10] = {0.0};
double Acc3S_AA[10][10] = {0.0};
double Acc1S_AA_Err[10][10] = {0.0};
double Acc2S_AA_Err[10][10] = {0.0};
double Acc3S_AA_Err[10][10] = {0.0};

TH1F *hAccInt1SNoW;
TH1F *hAccInt1Spp;
TH1F *hAccInt1SAA;

TH1F *hAccInt2SNoW;
TH1F *hAccInt2Spp;
TH1F *hAccInt2SAA;

TH1F *hAccInt3SNoW;
TH1F *hAccInt3Spp;
TH1F *hAccInt3SAA;

TH1F *hAccPt1SNoW;
TH1F *hAccPt1Spp;
TH1F *hAccPt1SAA;

TH1F *hAccPt2SNoW;
TH1F *hAccPt2Spp;
TH1F *hAccPt2SAA;

TH1F *hAccPt3SNoW;
TH1F *hAccPt3Spp;
TH1F *hAccPt3SAA;

TH1F *hAccY1SNoW;
TH1F *hAccY1Spp;
TH1F *hAccY1SAA;

TH1F *hAccY2SNoW;
TH1F *hAccY2Spp;
TH1F *hAccY2SAA;

TH1F *hAccY3SNoW;
TH1F *hAccY3Spp;
TH1F *hAccY3SAA;

// New weighting function 20161207
TF1* fWgtPP1 = new TF1("fWgtPP1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtPP2 = new TF1("fWgtPP2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");
TF1* fWgtAA1 = new TF1("fWgtAA1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtAA2 = new TF1("fWgtAA2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");

TH1F *hPadPt = new TH1F("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1F *hPadY = new TH1F("hPadY",";|y|;Acceptance",10,0.0,2.4);

TLatex *lt1 = new TLatex();

float getAcceptanceSingleBin(int state=k1S, int icoll=kNoWeight, TString kineCut="", int varianceID=0);
void getAcceptance_6Bins2S(int doDraw = 0) { // doDraw == 1 : draw, == 0 : no draw
  gStyle->SetOptStat(0);

  hPadPt->GetYaxis()->CenterTitle();
  hPadPt->GetYaxis()->SetTitleOffset(1.3);
  hPadPt->GetXaxis()->CenterTitle();
  hPadY->GetYaxis()->CenterTitle();
  hPadY->GetYaxis()->SetTitleOffset(1.3);
  hPadY->GetXaxis()->CenterTitle();
  lt1->SetNDC();


  fWgtPP1->SetParameters( 0.988141, 3.0971, 1.81891, 10.0239);
  fWgtPP2->SetParameters( 11.518, 7.53196, 2.38444, 2.68481);
  fWgtAA1->SetParameters( 1.0001, 5.1, 2.0024, 12.4243);
  fWgtAA2->SetParameters( 3.46994, 11.8612, 2.10006, 3.25859);

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample
  int nPtBins = 1;  double* ptBin = NULL;  int nYBins = 1;  double *yBin = NULL;
  const int nPtBins1 = 1;   double ptBin1[nPtBins1+1] = {0,30};  const int nYBins1 = 1;   double yBin1[nYBins1+1] = {0,2.4}; // integrated
  const int nPtBins2 = 6;   double ptBin2[nPtBins2+1] = {0,2,4,6,9,12,30};  const int nYBins2 = 1;   double yBin2[nYBins2+1] = {0,2.4}; // 1S
  const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,30};  const int nYBins3 = 6;   double yBin3[nYBins3+1] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4}; // 1S
  const int nPtBins4 = 6;   double ptBin4[nPtBins4+1] = {0,2,4,6,9,12,30};  const int nYBins4 = 1;   double yBin4[nYBins4+1] = {0,2.4}; // 2S
  const int nPtBins5 = 1;   double ptBin5[nPtBins5+1] = {0,30};  const int nYBins5 = 3;   double yBin5[nYBins5+1] = {0,0.8,1.6,2.4}; // 2S
  const int nPtBins6 = 2;   double ptBin6[nPtBins6+1] = {0,6,30};  const int nYBins6 = 1;   double yBin6[nYBins6+1] = {0,2.4}; // 3S
  const int nPtBins7 = 1;   double ptBin7[nPtBins7+1] = {0,30};  const int nYBins7 = 2;   double yBin7[nYBins7+1] = {0,1.2,2.4}; // 3S

  // bin 1 : integrated
  // bin 2, 4 : pT 
  // bin 3, 5 : y
  
  TFile *out = new TFile("acceptance_wgt_6BinsFor2S.root","RECREATE");
  //TFile *out = new TFile(Form("acceptance_bin%d_wgt.root",bin),"RECREATE");

  int bin = 0;

  for(int iCat = 0; iCat < 7; iCat++){
    bin = iCat+1;
    cout<<"Which bin : "<<bin<<endl;

    if ( bin==1 ) {
      nPtBins = nPtBins1;  ptBin=ptBin1;     nYBins = nYBins1; yBin=yBin1;
    }
    else if ( bin==2 ) {
      nPtBins = nPtBins2;  ptBin=ptBin2;     nYBins = nYBins2; yBin=yBin2;
    }
    else if ( bin==3 ) {
      nPtBins = nPtBins3;  ptBin=ptBin3;     nYBins = nYBins3; yBin=yBin3;
    }
    else if ( bin==4 ) {
      nPtBins = nPtBins4;  ptBin=ptBin4;     nYBins = nYBins4; yBin=yBin4;
    }
    else if ( bin==5 ) {
      nPtBins = nPtBins5;  ptBin=ptBin5;     nYBins = nYBins5; yBin=yBin5;
    }
    else if ( bin==6 ) {
      nPtBins = nPtBins6;  ptBin=ptBin6;     nYBins = nYBins6; yBin=yBin6;
    }
    else if ( bin==7 ) {
      nPtBins = nPtBins7;  ptBin=ptBin7;     nYBins = nYBins7; yBin=yBin7;
    }

    if( bin == 1 ){
      hAccInt1SNoW = new TH1F("hAccInt1SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin1);
      hAccInt1Spp = new TH1F("hAccInt1Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin1);
      hAccInt1SAA = new TH1F("hAccInt1SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin1);
      hAccInt2SNoW = new TH1F("hAccInt2SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin1);
      hAccInt3SNoW = new TH1F("hAccInt3SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin1);
      hAccInt2Spp = new TH1F("hAccInt2Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin1);
      hAccInt3Spp = new TH1F("hAccInt3Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin1);
      hAccInt2SAA = new TH1F("hAccInt2SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin1);
      hAccInt3SAA = new TH1F("hAccInt3SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin1);
      hAccInt1SNoW->Sumw2();
      hAccInt1Spp->Sumw2();
      hAccInt1SAA->Sumw2();
      hAccInt2SNoW->Sumw2();
      hAccInt2Spp->Sumw2();
      hAccInt2SAA->Sumw2();
      hAccInt3SNoW->Sumw2();
      hAccInt3Spp->Sumw2();
      hAccInt3SAA->Sumw2();
      hAccInt1SNoW->SetMarkerStyle(20);hAccInt1SNoW->SetMarkerColor(kBlue+2);
      hAccInt1Spp->SetMarkerStyle(20);hAccInt1Spp->SetMarkerColor(kBlue+2);
      hAccInt1SAA->SetMarkerStyle(20);hAccInt1SAA->SetMarkerColor(kBlue+2);
      hAccInt2SNoW->SetMarkerStyle(20);hAccInt2SNoW->SetMarkerColor(kBlue+2);
      hAccInt3SNoW->SetMarkerStyle(20);hAccInt3SNoW->SetMarkerColor(kBlue+2);
      hAccInt2Spp->SetMarkerStyle(20);hAccInt2Spp->SetMarkerColor(kBlue+2);
      hAccInt3Spp->SetMarkerStyle(20);hAccInt3Spp->SetMarkerColor(kBlue+2);
      hAccInt2SAA->SetMarkerStyle(20);hAccInt2SAA->SetMarkerColor(kBlue+2);
      hAccInt3SAA->SetMarkerStyle(20);hAccInt3SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 2 ){
      hAccPt1SNoW = new TH1F("hAccPt1SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin2);
      hAccPt1Spp = new TH1F("hAccPt1Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin2);
      hAccPt1SAA = new TH1F("hAccPt1SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(1S)",nPtBins,ptBin2);
      hAccPt1SNoW->Sumw2();
      hAccPt1Spp->Sumw2();
      hAccPt1SAA->Sumw2();
      hAccPt1SNoW->SetMarkerStyle(20);hAccPt1SNoW->SetMarkerColor(kBlue+2);
      hAccPt1Spp->SetMarkerStyle(20);hAccPt1Spp->SetMarkerColor(kBlue+2);
      hAccPt1SAA->SetMarkerStyle(20);hAccPt1SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 4 ){
      hAccPt2SNoW = new TH1F("hAccPt2SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin4);
      hAccPt2Spp = new TH1F("hAccPt2Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin4);
      hAccPt2SAA = new TH1F("hAccPt2SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(2S)",nPtBins,ptBin4);
      hAccPt2SNoW->Sumw2();
      hAccPt2Spp->Sumw2();
      hAccPt2SAA->Sumw2();
      hAccPt2SNoW->SetMarkerStyle(20);hAccPt2SNoW->SetMarkerColor(kBlue+2);
      hAccPt2Spp->SetMarkerStyle(20);hAccPt2Spp->SetMarkerColor(kBlue+2);
      hAccPt2SAA->SetMarkerStyle(20);hAccPt2SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 3 ) {
      hAccY1SNoW = new TH1F("hAccY1SNoW",";|y|;Acceptance of #Upsilon(1S)",nYBins,yBin3);
      hAccY1Spp = new TH1F("hAccY1Spp",";|y|;Acceptance of #Upsilon(1S)",nYBins,yBin3);
      hAccY1SAA = new TH1F("hAccY1SAA",";|y|;Acceptance of #Upsilon(1S)",nYBins,yBin3);
      hAccY1SNoW->Sumw2();
      hAccY1Spp->Sumw2();
      hAccY1SAA->Sumw2();
      hAccY1SNoW->SetMarkerStyle(20);hAccY1SNoW->SetMarkerColor(kBlue+2);
      hAccY1Spp->SetMarkerStyle(20);hAccY1Spp->SetMarkerColor(kBlue+2);
      hAccY1SAA->SetMarkerStyle(20);hAccY1SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 5 ) {
      hAccY2SNoW = new TH1F("hAccY2SNoW",";|y|;Acceptance of #Upsilon(2S)",nYBins,yBin5);
      hAccY2Spp = new TH1F("hAccY2Spp",";|y|;Acceptance of #Upsilon(2S)",nYBins,yBin5);
      hAccY2SAA = new TH1F("hAccY2SAA",";|y|;Acceptance of #Upsilon(2S)",nYBins,yBin5);
      hAccY2SNoW->Sumw2();
      hAccY2Spp->Sumw2();
      hAccY2SAA->Sumw2();
      hAccY2SNoW->SetMarkerStyle(20);hAccY2SNoW->SetMarkerColor(kBlue+2);
      hAccY2Spp->SetMarkerStyle(20);hAccY2Spp->SetMarkerColor(kBlue+2);
      hAccY2SAA->SetMarkerStyle(20);hAccY2SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 6 ) {
      hAccPt3SNoW = new TH1F("hAccPt3SNoW",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin4);
      hAccPt3Spp = new TH1F("hAccPt3Spp",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin4);
      hAccPt3SAA = new TH1F("hAccPt3SAA",";p_{T} (GeV/c);Acceptance of #Upsilon(3S)",nPtBins,ptBin4);
      hAccPt3SNoW->Sumw2();
      hAccPt3Spp->Sumw2();
      hAccPt3SAA->Sumw2();
      hAccPt3SNoW->SetMarkerStyle(20);hAccPt3SNoW->SetMarkerColor(kBlue+2);
      hAccPt3Spp->SetMarkerStyle(20);hAccPt3Spp->SetMarkerColor(kBlue+2);
      hAccPt3SAA->SetMarkerStyle(20);hAccPt3SAA->SetMarkerColor(kBlue+2);
    }else if( bin == 7 ) {
      hAccY3SNoW = new TH1F("hAccY3SNoW",";|y|;Acceptance of #Upsilon(3S)",nYBins,yBin5);
      hAccY3Spp = new TH1F("hAccY3Spp",";|y|;Acceptance of #Upsilon(3S)",nYBins,yBin5);
      hAccY3SAA = new TH1F("hAccY3SAA",";|y|;Acceptance of #Upsilon(3S)",nYBins,yBin5);
      hAccY3SNoW->Sumw2();
      hAccY3Spp->Sumw2();
      hAccY3SAA->Sumw2();
      hAccY3SNoW->SetMarkerStyle(20);hAccY3SNoW->SetMarkerColor(kBlue+2);
      hAccY3Spp->SetMarkerStyle(20);hAccY3Spp->SetMarkerColor(kBlue+2);
      hAccY3SAA->SetMarkerStyle(20);hAccY3SAA->SetMarkerColor(kBlue+2);
    }


    for ( int ipt = 1 ; ipt<=nPtBins ; ipt++ ) { 
      for ( int iy =1 ; iy<=nYBins ; iy++) {
        TString kRange = Form("pt>%0.1f && pt<%0.1f && abs(y)>%0.1f && abs(y)<%0.1f",(float)ptBin[ipt-1],(float)ptBin[ipt], (float)yBin[iy-1], (float)yBin[iy] ) ;
        cout << "*===*===*===*===*  " << kRange << "  *===*===*===*===*" << endl;
        float y1s_noRwt = getAcceptanceSingleBin( 1, kNoWeight, kRange,kNoVar);
        float y1s_ppRwt = getAcceptanceSingleBin( 1, kPP, kRange,kNoVar);
        float y1s_ppRwt_ptPlus  = getAcceptanceSingleBin( 1, kPP, kRange, kPtPlus);
        float y1s_ppRwt_ptMinus = getAcceptanceSingleBin( 1, kPP, kRange, kPtMinus);
        float y1s_ppRwt_yPlus  = getAcceptanceSingleBin( 1, kPP, kRange, kYPlus);
        float y1s_ppRwt_yMinus = getAcceptanceSingleBin( 1, kPP, kRange, kYMinus);

        float y2s_noRwt = getAcceptanceSingleBin( 2, kNoWeight, kRange,kNoVar);
        float y2s_ppRwt = getAcceptanceSingleBin( 2, kPP, kRange,kNoVar);
        float y2s_ppRwt_ptPlus  = getAcceptanceSingleBin( 2, kPP, kRange, kPtPlus);
        float y2s_ppRwt_ptMinus = getAcceptanceSingleBin( 2, kPP, kRange, kPtMinus);
        float y2s_ppRwt_yPlus  = getAcceptanceSingleBin( 2, kPP, kRange, kYPlus);
        float y2s_ppRwt_yMinus = getAcceptanceSingleBin( 2, kPP, kRange, kYMinus);

        float y3s_noRwt = getAcceptanceSingleBin( 3, kNoWeight, kRange,kNoVar);
        float y3s_ppRwt = getAcceptanceSingleBin( 3, kPP, kRange,kNoVar);
        float y3s_ppRwt_ptPlus  = getAcceptanceSingleBin( 3, kPP, kRange, kPtPlus);
        float y3s_ppRwt_ptMinus = getAcceptanceSingleBin( 3, kPP, kRange, kPtMinus);
        float y3s_ppRwt_yPlus  = getAcceptanceSingleBin( 3, kPP, kRange, kYPlus);
        float y3s_ppRwt_yMinus = getAcceptanceSingleBin( 3, kPP, kRange, kYMinus);

        float err1s_ppPt = max ( abs(y1s_ppRwt_ptPlus/y1s_ppRwt-1), abs(y1s_ppRwt_ptMinus/y1s_ppRwt-1) ) ;
        float err2s_ppPt = max ( abs(y2s_ppRwt_ptPlus/y2s_ppRwt-1), abs(y2s_ppRwt_ptMinus/y2s_ppRwt-1) ) ;
        float err3s_ppPt = max ( abs(y3s_ppRwt_ptPlus/y3s_ppRwt-1), abs(y3s_ppRwt_ptMinus/y3s_ppRwt-1) ) ;
        float err1s_ppY = max ( abs(y1s_ppRwt_yPlus/y1s_ppRwt-1), abs(y1s_ppRwt_yMinus/y1s_ppRwt-1) ) ;
        float err2s_ppY = max ( abs(y2s_ppRwt_yPlus/y2s_ppRwt-1), abs(y2s_ppRwt_yMinus/y2s_ppRwt-1) ) ;
        float err3s_ppY = max ( abs(y3s_ppRwt_yPlus/y3s_ppRwt-1), abs(y3s_ppRwt_yMinus/y3s_ppRwt-1) ) ;

        float y1s_aaRwt = getAcceptanceSingleBin( 1, kAA, kRange,kNoVar);
        float y1s_aaRwt_ptPlus  = getAcceptanceSingleBin( 1, kAA, kRange, kPtPlus);
        float y1s_aaRwt_ptMinus = getAcceptanceSingleBin( 1, kAA, kRange, kPtMinus);
        float y1s_aaRwt_yPlus  = getAcceptanceSingleBin( 1, kAA, kRange, kYPlus);
        float y1s_aaRwt_yMinus = getAcceptanceSingleBin( 1, kAA, kRange, kYMinus);

        float y2s_aaRwt = getAcceptanceSingleBin( 2, kAA, kRange,kNoVar);
        float y2s_aaRwt_ptPlus  = getAcceptanceSingleBin( 2, kAA, kRange, kPtPlus);
        float y2s_aaRwt_ptMinus = getAcceptanceSingleBin( 2, kAA, kRange, kPtMinus);
        float y2s_aaRwt_yPlus  = getAcceptanceSingleBin( 2, kAA, kRange, kYPlus);
        float y2s_aaRwt_yMinus = getAcceptanceSingleBin( 2, kAA, kRange, kYMinus);

        float y3s_aaRwt = getAcceptanceSingleBin( 3, kAA, kRange,kNoVar);
        float y3s_aaRwt_ptPlus  = getAcceptanceSingleBin( 3, kAA, kRange, kPtPlus);
        float y3s_aaRwt_ptMinus = getAcceptanceSingleBin( 3, kAA, kRange, kPtMinus);
        float y3s_aaRwt_yPlus  = getAcceptanceSingleBin( 3, kAA, kRange, kYPlus);
        float y3s_aaRwt_yMinus = getAcceptanceSingleBin( 3, kAA, kRange, kYMinus);

        float err1s_aaPt = max ( abs(y1s_aaRwt_ptPlus/y1s_aaRwt-1), abs(y1s_aaRwt_ptMinus/y1s_aaRwt-1) ) ;
        float err2s_aaPt = max ( abs(y2s_aaRwt_ptPlus/y2s_aaRwt-1), abs(y2s_aaRwt_ptMinus/y2s_aaRwt-1) ) ;
        float err3s_aaPt = max ( abs(y3s_aaRwt_ptPlus/y3s_aaRwt-1), abs(y3s_aaRwt_ptMinus/y3s_aaRwt-1) ) ;
        float err1s_aaY = max ( abs(y1s_aaRwt_yPlus/y1s_aaRwt-1), abs(y1s_aaRwt_yMinus/y1s_aaRwt-1) ) ;
        float err2s_aaY = max ( abs(y2s_aaRwt_yPlus/y2s_aaRwt-1), abs(y2s_aaRwt_yMinus/y2s_aaRwt-1) ) ;
        float err3s_aaY = max ( abs(y3s_aaRwt_yPlus/y3s_aaRwt-1), abs(y3s_aaRwt_yMinus/y3s_aaRwt-1) ) ;

        // dmoon added systematcis : unweighted - weighted
        float err1s_pp = abs(y1s_noRwt - y1s_ppRwt)/y1s_ppRwt;
        float err1s_aa = abs(y1s_noRwt - y1s_aaRwt)/y1s_aaRwt;
        float err2s_pp = abs(y2s_noRwt - y2s_ppRwt)/y2s_ppRwt;
        float err2s_aa = abs(y2s_noRwt - y2s_aaRwt)/y2s_aaRwt;
        float err3s_pp = abs(y3s_noRwt - y3s_ppRwt)/y3s_ppRwt;
        float err3s_aa = abs(y3s_noRwt - y3s_aaRwt)/y3s_aaRwt;

        Acc1S_NoW[ipt][iy] = y1s_noRwt;
        Acc2S_NoW[ipt][iy] = y2s_noRwt;
        Acc3S_NoW[ipt][iy] = y3s_noRwt;
        Acc1S_NoW_Err[ipt][iy] = err1s_ppPt;
        Acc2S_NoW_Err[ipt][iy] = err2s_ppPt;
        Acc3S_NoW_Err[ipt][iy] = err3s_ppPt;

        Acc1S_pp[ipt][iy] = y1s_ppRwt;
        Acc2S_pp[ipt][iy] = y2s_ppRwt;
        Acc3S_pp[ipt][iy] = y3s_ppRwt;
        Acc1S_pp_Err[ipt][iy] = err1s_pp;
        Acc2S_pp_Err[ipt][iy] = err2s_pp;
        Acc3S_pp_Err[ipt][iy] = err3s_pp;

        Acc1S_AA[ipt][iy] = y1s_aaRwt;
        Acc2S_AA[ipt][iy] = y2s_aaRwt;
        Acc3S_AA[ipt][iy] = y3s_aaRwt;
        Acc1S_AA_Err[ipt][iy] = err1s_aa;
        Acc2S_AA_Err[ipt][iy] = err2s_aa;
        Acc3S_AA_Err[ipt][iy] = err3s_aa;


        cout << "Un-weighed acceptance : (1S) = " << y1s_noRwt << ",   (2S) = " << y2s_noRwt << ",   (3S) = "<< y3s_noRwt<<endl;
        cout << "Re-weighed Acceptance : " << "pp(1S) = " << y1s_ppRwt << ",  PbPb(1S) = " << y1s_aaRwt ;
        cout << " pp(2S) = " << y2s_ppRwt << ",  PbPb(2S) = " << y2s_aaRwt ;
        cout << " pp(3S) = " << y3s_ppRwt << ",  PbPb(3S) = " << y3s_aaRwt ;
        //cout << ",     (2S/1S,AA/pp) acc. ratio = " <<  (y2s_aaRwt/y2s_ppRwt) /(y1s_aaRwt/y1s_ppRwt)  << endl;
        cout<<""<<endl;

        cout << int(1000*y1s_noRwt)*0.001 << " & " << int(1000*y1s_ppRwt)*0.001<< " & " << int(1000*y1s_aaRwt)*0.001<< " & " ;
        cout << int(1000*y2s_noRwt)*0.001 << " & " << int(1000*y2s_ppRwt)*0.001<< " & " << int(1000*y2s_aaRwt)*0.001 << " & ";
        cout << int(1000*y3s_noRwt)*0.001 << " & " << int(1000*y3s_ppRwt)*0.001<< " & " << int(1000*y3s_aaRwt)*0.001 << " & ";
        cout <<""<<endl;
        //cout << int(1000* (y2s_aaRwt/y2s_ppRwt) / (y1s_aaRwt/y1s_ppRwt) )*0.001 << endl; // Final correction
        cout << int(1000*err1s_ppPt)*0.1 << "\% & " << int(1000*err1s_aaPt)*0.1 << "\% & " ;
        cout << int(1000*err2s_ppPt)*0.1 << "\% & " << int(1000*err2s_aaPt)*0.1 << "\% & " ;
        cout << int(1000*err3s_ppPt)*0.1 << "\% & " << int(1000*err3s_aaPt)*0.1 << "\% & " ;
        cout << int(1000*err1s_ppY)*0.1 << "\% & " << int(1000*err1s_aaY)*0.1 << "\% & " ;
        cout << int(1000*err2s_ppY)*0.1 << "\% & " << int(1000*err2s_aaY)*0.1 << "\% & " ;
        cout << int(1000*err3s_ppY)*0.1 << "\% & " << int(1000*err3s_aaY)*0.1 << "\% & " ;
        cout <<""<<endl;

        cout<<"##### Systematics #####"<<endl;
        cout<<"pp(1S) "<<"AA(1S) "<<"pp(2S) "<<"AA(2S) "<<"pp(3S) "<<"AA(3S) "<<endl;
        cout<< int(1000*err1s_pp)*0.1 << "\% & " << int(1000*err1s_aa)*0.1 << "\% & ";
        cout<< int(1000*err2s_pp)*0.1 << "\% & " << int(1000*err2s_aa)*0.1 << "\% & ";
        cout<< int(1000*err3s_pp)*0.1 << "\% & " << int(1000*err3s_aa)*0.1 << "\% & ";
        //cout << int(1000*sqrt( err1s_ppPt*err1s_ppPt +  err2s_ppPt* err2s_ppPt + err1s_ppY*err1s_ppY + err2s_ppY*err2s_ppY + err1s_aaPt*err1s_aaPt +  err2s_aaPt* err2s_aaPt + err1s_aaY*err1s_aaY + err2s_aaY*err2s_aaY ) ) *0.1 << "\%" << endl;
        cout<<""<<endl;

        if(bin == 1){
          hAccInt1SNoW->SetBinContent(ipt,Acc1S_NoW[ipt][iy]);
          hAccInt1Spp->SetBinContent(ipt,Acc1S_pp[ipt][iy]);
          hAccInt1SAA->SetBinContent(ipt,Acc1S_AA[ipt][iy]);
          hAccInt2SNoW->SetBinContent(ipt,Acc2S_NoW[ipt][iy]);
          hAccInt3SNoW->SetBinContent(ipt,Acc3S_NoW[ipt][iy]);
          hAccInt2Spp->SetBinContent(ipt,Acc2S_pp[ipt][iy]);
          hAccInt3Spp->SetBinContent(ipt,Acc3S_pp[ipt][iy]);
          hAccInt2SAA->SetBinContent(ipt,Acc2S_AA[ipt][iy]);
          hAccInt3SAA->SetBinContent(ipt,Acc3S_AA[ipt][iy]);
          hAccInt1SNoW->SetBinError(ipt,Acc1S_NoW_Err[ipt][iy]);
          hAccInt1Spp->SetBinError(ipt,Acc1S_pp_Err[ipt][iy]);
          hAccInt1SAA->SetBinError(ipt,Acc1S_AA_Err[ipt][iy]);
          hAccInt2SNoW->SetBinError(ipt,Acc2S_NoW_Err[ipt][iy]);
          hAccInt3SNoW->SetBinError(ipt,Acc3S_NoW_Err[ipt][iy]);
          hAccInt2Spp->SetBinError(ipt,Acc2S_pp_Err[ipt][iy]);
          hAccInt3Spp->SetBinError(ipt,Acc3S_pp_Err[ipt][iy]);
          hAccInt2SAA->SetBinError(ipt,Acc2S_AA_Err[ipt][iy]);
          hAccInt3SAA->SetBinError(ipt,Acc3S_AA_Err[ipt][iy]);

        }else if(bin == 2){
          hAccPt1SNoW->SetBinContent(ipt,Acc1S_NoW[ipt][iy]);
          hAccPt1Spp->SetBinContent(ipt,Acc1S_pp[ipt][iy]);
          hAccPt1SAA->SetBinContent(ipt,Acc1S_AA[ipt][iy]);
          hAccPt1SNoW->SetBinError(ipt,Acc1S_NoW_Err[ipt][iy]);
          hAccPt1Spp->SetBinError(ipt,Acc1S_pp_Err[ipt][iy]);
          hAccPt1SAA->SetBinError(ipt,Acc1S_AA_Err[ipt][iy]);
        }else if(bin == 3){
          hAccY1SNoW->SetBinContent(iy,Acc1S_NoW[ipt][iy]);
          hAccY1Spp->SetBinContent(iy,Acc1S_pp[ipt][iy]);
          hAccY1SAA->SetBinContent(iy,Acc1S_AA[ipt][iy]);
          hAccY1SNoW->SetBinError(iy,Acc1S_NoW_Err[ipt][iy]);
          hAccY1Spp->SetBinError(iy,Acc1S_pp_Err[ipt][iy]);
          hAccY1SAA->SetBinError(iy,Acc1S_AA_Err[ipt][iy]);
        }else if(bin == 4){
          hAccPt2SNoW->SetBinContent(ipt,Acc2S_NoW[ipt][iy]);
          hAccPt2Spp->SetBinContent(ipt,Acc2S_pp[ipt][iy]);
          hAccPt2SAA->SetBinContent(ipt,Acc2S_AA[ipt][iy]);
          hAccPt2SNoW->SetBinError(ipt,Acc2S_NoW_Err[ipt][iy]);
          hAccPt2Spp->SetBinError(ipt,Acc2S_pp_Err[ipt][iy]);
          hAccPt2SAA->SetBinError(ipt,Acc2S_AA_Err[ipt][iy]);
        }else if(bin == 5){
          hAccY2SNoW->SetBinContent(iy,Acc2S_NoW[ipt][iy]);
          hAccY2Spp->SetBinContent(iy,Acc2S_pp[ipt][iy]);
          hAccY2SAA->SetBinContent(iy,Acc2S_AA[ipt][iy]);
          hAccY2SNoW->SetBinError(iy,Acc2S_NoW_Err[ipt][iy]);
          hAccY2Spp->SetBinError(iy,Acc2S_pp_Err[ipt][iy]);
          hAccY2SAA->SetBinError(iy,Acc2S_AA_Err[ipt][iy]);
        }else if(bin == 6){
          hAccPt3SNoW->SetBinContent(ipt,Acc3S_NoW[ipt][iy]);
          hAccPt3Spp->SetBinContent(ipt,Acc3S_pp[ipt][iy]);
          hAccPt3SAA->SetBinContent(ipt,Acc3S_AA[ipt][iy]);
          hAccPt3SNoW->SetBinError(ipt,Acc3S_NoW_Err[ipt][iy]);
          hAccPt3Spp->SetBinError(ipt,Acc3S_pp_Err[ipt][iy]);
          hAccPt3SAA->SetBinError(ipt,Acc3S_AA_Err[ipt][iy]);
        }else if(bin == 7){
          hAccY3SNoW->SetBinContent(iy,Acc3S_NoW[ipt][iy]);
          hAccY3Spp->SetBinContent(iy,Acc3S_pp[ipt][iy]);
          hAccY3SAA->SetBinContent(iy,Acc3S_AA[ipt][iy]);
          hAccY3SNoW->SetBinError(iy,Acc3S_NoW_Err[ipt][iy]);
          hAccY3Spp->SetBinError(iy,Acc3S_pp_Err[ipt][iy]);
          hAccY3SAA->SetBinError(iy,Acc3S_AA_Err[ipt][iy]);
        }
      } // ipt
    } // iy
  } // Cat_
  if(doDraw == 1){
    /*
    TCanvas *c1 = new TCanvas("c1","",1800,600);
    c1->Divide(3,1);
    if(bin == 1 || bin == 2){
      c1->cd(1); hPadPt->Draw(); hAccPt1->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (1S)");
      c1->cd(2); hPadPt->Draw(); hAccPt2->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (2S)");
      c1->cd(3); hPadPt->Draw(); hAccPt3->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (3S)");
      if(bin == 1) c1->SaveAs("plot_acc_int_wgt.png");
      if(bin == 2) c1->SaveAs("plot_acc_pT_wgt.png");
    }
    if(bin == 3){
      c1->cd(1); hPadY->Draw(); hAccY1->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (1S)");
      c1->cd(2); hPadY->Draw(); hAccY2->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (2S)");
      c1->cd(3); hPadY->Draw(); hAccY3->Draw("same"); lt1->DrawLatex(0.7,0.8,"#Upsilon (3S)");
      c1->SaveAs("plot_acc_y_wgt.png");
    }
    */
  }
  out->cd();
  /*
  if( bin == 2 ) {
    hAccPt1SNoW->SetName("hptAccNoW1S");
    hAccPt1Spp->SetName("hptAccPP1S");
    hAccPt1SAA->SetName("hptAccAA1S");
    hAccPt1SNoW->Write();
    hAccPt1Spp->Write();
    hAccPt1SAA->Write();
  }else if( bin == 3 ) {
    hAccY1SNoW->SetName("hrapAccNoW1S");
    hAccY1Spp->SetName("hrapAccPP1S");
    hAccY1SAA->SetName("hrapAccAA1S");
    hAccY1SNoW->Write();
    hAccY1Spp->Write();
    hAccY1SAA->Write();
  }else if( bin == 4 ) {
    hAccPt2SNoW->SetName("hptAccNoW2S");
    hAccPt3SNoW->SetName("hptAccNoW3S");
    hAccPt2Spp->SetName("hptAccPP2S");
    hAccPt3Spp->SetName("hptAccPP3S");
    hAccPt2SAA->SetName("hptAccAA2S");
    hAccPt3SAA->SetName("hptAccAA3S");
    hAccPt2SNoW->Write();
    hAccPt3SNoW->Write();
    hAccPt2Spp->Write();
    hAccPt3Spp->Write();
    hAccPt2SAA->Write();
    hAccPt3SAA->Write();
  }else if( bin == 5 ) {
    hAccY2SNoW->SetName("hrapAccNoW2S");
    hAccY3SNoW->SetName("hrapAccNoW3S");
    hAccY2Spp->SetName("hrapAccPP2S");
    hAccY3Spp->SetName("hrapAccPP3S");
    hAccY2SAA->SetName("hrapAccAA2S");
    hAccY3SAA->SetName("hrapAccAA3S");
    hAccY2SNoW->Write();
    hAccY3SNoW->Write();
    hAccY2Spp->Write();
    hAccY3Spp->Write();
    hAccY2SAA->Write();
    hAccY3SAA->Write();
  }
  */
  hAccInt1SNoW->SetName("hIntAccNoW1S");
  hAccInt1Spp->SetName("hIntAccPP1S");
  hAccInt1SAA->SetName("hIntAccAA1S");
  hAccInt2SNoW->SetName("hIntAccNoW2S");
  hAccInt2Spp->SetName("hIntAccPP2S");
  hAccInt2SAA->SetName("hIntAccAA2S");
  hAccInt3SNoW->SetName("hIntAccNoW3S");
  hAccInt3Spp->SetName("hIntAccPP3S");
  hAccInt3SAA->SetName("hIntAccAA3S");
  hAccInt1SNoW->Write();
  hAccInt1Spp->Write();
  hAccInt1SAA->Write();
  hAccInt2SNoW->Write();
  hAccInt2Spp->Write();
  hAccInt2SAA->Write();
  hAccInt3SNoW->Write();
  hAccInt3Spp->Write();
  hAccInt3SAA->Write();
  hAccPt1SNoW->SetName("hptAccNoW1S");
  hAccPt1Spp->SetName("hptAccPP1S");
  hAccPt1SAA->SetName("hptAccAA1S");
  hAccPt1SNoW->Write();
  hAccPt1Spp->Write();
  hAccPt1SAA->Write();
  hAccY1SNoW->SetName("hrapAccNoW1S");
  hAccY1Spp->SetName("hrapAccPP1S");
  hAccY1SAA->SetName("hrapAccAA1S");
  hAccY1SNoW->Write();
  hAccY1Spp->Write();
  hAccY1SAA->Write();
  hAccPt2SNoW->SetName("hptAccNoW2S");
  hAccPt3SNoW->SetName("hptAccNoW3S");
  hAccPt2Spp->SetName("hptAccPP2S");
  hAccPt3Spp->SetName("hptAccPP3S");
  hAccPt2SAA->SetName("hptAccAA2S");
  hAccPt3SAA->SetName("hptAccAA3S");
  hAccPt2SNoW->Write();
  hAccPt3SNoW->Write();
  hAccPt2Spp->Write();
  hAccPt3Spp->Write();
  hAccPt2SAA->Write();
  hAccPt3SAA->Write();
  hAccY2SNoW->SetName("hrapAccNoW2S");
  hAccY3SNoW->SetName("hrapAccNoW3S");
  hAccY2Spp->SetName("hrapAccPP2S");
  hAccY3Spp->SetName("hrapAccPP3S");
  hAccY2SAA->SetName("hrapAccAA2S");
  hAccY3SAA->SetName("hrapAccAA3S");
  hAccY2SNoW->Write();
  hAccY3SNoW->Write();
  hAccY2Spp->Write();
  hAccY3Spp->Write();
  hAccY2SAA->Write();
  hAccY3SAA->Write();
  out->Write();
}

float getAcceptanceSingleBin(int state, int icoll, TString kineCut, int varianceID) {
  
  TH1::SetDefaultSumw2();
  
  float muPtCut = 4; // for acceptance 
  float muEtaCut = 2.4;
  
  TH1D* hGen; // in rapidity bins, centrality bins
  TH1D* hGenAcc; // in rapidity bins, centrality bins
  TH1D* hAcc; // Acceptance rate : in rapidity bins, centrality bins 
  
  TFile* f = NULL;
  if ( state == k1S){ 
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161015_unIdentified_wgt.root");
  } else if ( state == k2S) {
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161015_unIdentified_wgt.root");
  } else if ( state == k3S) {
    f = new TFile("../skimmedFiles/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161015_unIdentified_wgt.root");
  }

  TTree* mmGen = (TTree*)f->Get("mmGen");
  TString accCut = "pt1>4 && pt2>4 && abs(eta1)<2.4 && abs(eta2)<2.4";
  
  TString ptWeight;
  if ( icoll == kNoWeight ) ptWeight = "1";
  else if ( icoll == kPP ) { 
    if (state == k1S)  ptWeight = "weightPP";
    if (state == k2S)  ptWeight = "weightPP";
    if (state == k3S)  ptWeight = "weightPP"; // k3S check
    //if (state == k1S)  pTweight = fWgtPP1->Eval(ptval);
    //if (state == k2S)  pTweight = fWgtPP2->Eval(ptval);
    //if (state == k3S)  pTweight = fWgtPP2->Eval(ptval); 
    //if (state == k1S)  ptWeight = "1";
    //if (state == k2S)  ptWeight = "1";
    //if (state == k3S)  ptWeight = "1"; // k3S check
    //if (state == k1S)  ptWeight = "1.09 - 0.022*pt";
    //if (state == k2S)  ptWeight = "0.69 + 0.073*pt";
    //if (state == k3S)  ptWeight = "0.69 + 0.073*pt"; // k3S check
  }    
  else if ( icoll == kAA ) { 
    if (state == k1S)  ptWeight = "weightAA";
    if (state == k2S)  ptWeight = "weightAA";
    if (state == k3S)  ptWeight = "weightAA"; // k3S check
    //if (state == k1S)  pTweight = fWgtAA1->Eval(ptval);
    //if (state == k2S)  pTweight = fWgtAA2->Eval(ptval);
    //if (state == k3S)  pTweight = fWgtAA2->Eval(ptval); 
    //if (state == k1S)  ptWeight = "1";
    //if (state == k2S)  ptWeight = "1";
    //if (state == k3S)  ptWeight = "1"; // k3S check
    //if (state == k1S)  ptWeight = "1.08 - 0.015*pt";
    //if (state == k2S)  ptWeight = "0.68 + 0.1*pt";
    //if (state == k3S)  ptWeight = "0.69 + 0.073*pt"; // k3S check
  }    

  //cout<<"dmoon chk :: state = "<<state<<", pTweight = "<<pTweight<<endl;

  if ( varianceID == kPtPlus ) 
    ptWeight = "("+ptWeight+") * (1.0)" ;
  else if ( varianceID == kPtMinus ) 
    ptWeight = "("+ptWeight+") * (1.0)" ;
  else if ( varianceID == kYPlus ) 
    ptWeight = "("+ptWeight+") * (1.0)" ;
  else if ( varianceID == kYMinus ) 
    ptWeight = "("+ptWeight+") * (1.0)" ;

  /* // plus & minus
  if ( varianceID == kPtPlus ) 
    ptWeight = "("+ptWeight+")* (0.8 + 0.0133*pt)" ;
  else if ( varianceID == kPtMinus ) 
    ptWeight = "("+ptWeight+")* (1.2 - 0.0133*pt)" ;
  else if ( varianceID == kYPlus ) 
    ptWeight = "("+ptWeight+")* (0.8 + 0.167*abs(y))" ;
  else if ( varianceID == kYMinus ) 
    ptWeight = "("+ptWeight+")* (1.2 - 0.167*abs(y))" ;
  */

  //  cout << " ptWeight = " << ptWeight << endl;
  
  hGen = new TH1D(Form("hGen_state%d_icoll%d",state,icoll),";Rapidity;",1,-2.4,2.4);
  hGenAcc = (TH1D*)hGen->Clone(Form("hGenAcc_state%d_icoll%d",state,icoll));
  //// all GEN
  mmGen->Draw(Form("y>>%s",hGen->GetName()), Form("(%s)*(%s)", kineCut.Data(), ptWeight.Data()));
  //// GEN with sglMuAccCut 
  mmGen->Draw(Form("y>>%s",hGenAcc->GetName()), Form("((%s)&&(%s))*(%s)", kineCut.Data(), accCut.Data(), ptWeight.Data()));
  handsomeTH1(hGenAcc,1);
  //// --- Acceptance rate : 
  hAcc = (TH1D*)hGenAcc->Clone(Form("hAcc_state%d_icoll%d",state,icoll));
  hAcc->Divide(hGen);
  
  
  /*
  TCanvas* c1 = new TCanvas("c1","",400,400);
  hGen->Draw("hist");
  hGenAcc->Draw("same e");
  
  TCanvas* c2 = new TCanvas("c2","",400,400);
  hAcc->SetAxisRange(0,1,"Y");
  hAcc->Draw();
  */
    
  return hAcc->GetBinContent(1);
}


