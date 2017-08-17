#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
//#include "../../cutsAndBin.h"
using namespace std;
// Editor : Dongho Moon (dmoon@cern.ch), Geonhee Oh (geonhee.oh@cern.ch)

// Define weighting functions for systematics
TF1 *fWgtHigh = new TF1("fWgtHigh","[0]",0.0,30.0);
TF1 *fWgtLow = new TF1("fWgtLow","[0]",0.0,30.0);

// Define up and down 20 % with pt
TF1 *fWgtPtp = new TF1("fWgtPtp","0.8 + 0.0133*x",0.0,30.0); // 20 % up
TF1 *fWgtPtm = new TF1("fWgtPtm","1.2 - 0.0133*x",0.0,30.0); // 20 % down
TF1 *fWgtYp = new TF1("fWgtYp","0.8 + 0.167*x",0.0,2.4); // 20 % up
TF1 *fWgtYm = new TF1("fWgtYm","1.2 - 0.167*x",0.0,2.4); // 20 % down

TH1F *hPadPt = new TH1F("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1F *hPadPt1 = new TH1F("hPadPt1",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1F *hPadY = new TH1F("hPadY",";y_{CM};Acceptance",10, -1.93, 1.93);
TH1F *hPadY1 = new TH1F("hPadY1",";y_{CM};Acceptance",10, -1.93, 1.93);

TLatex *lt1 = new TLatex();

void studyAccAna_KEEP_16_023(int Cat_ = 1) { // Cat_ == 0 (1S),  Cat_ == 1 (2S), Cat_ == 2 (3S)
  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptStat(0);

  cout<<" Calculating acceptance starts !!!"<<endl;
  if(Cat_ == 0) cout<<" %%%%% Upsilon 1S %%%%% "<<endl;
  if(Cat_ == 1) cout<<" %%%%% Upsilon 2S %%%%% "<<endl;
  if(Cat_ == 2) cout<<" %%%%% Upsilon 3S %%%%% "<<endl;

  hPadPt->GetYaxis()->CenterTitle();
  hPadPt->GetYaxis()->SetTitleOffset(1.3);
  hPadPt->GetXaxis()->CenterTitle();
  hPadY->GetYaxis()->CenterTitle();
  hPadY->GetYaxis()->SetTitleOffset(1.3);
  hPadY->GetXaxis()->CenterTitle();

  hPadPt1->GetYaxis()->CenterTitle();
  hPadPt1->GetYaxis()->SetTitleOffset(1.3);
  hPadPt1->GetXaxis()->CenterTitle();
  hPadY1->GetYaxis()->CenterTitle();
  hPadY1->GetYaxis()->SetTitleOffset(1.3);
  hPadY1->GetXaxis()->CenterTitle();

  lt1->SetNDC();

  TF1* fWgtPP;   //Nominal
  TF1* fWgtPP_Ap;//Variation A positive error for PP 
  TF1* fWgtPP_Bp;//Variation B positive error for PP
  TF1* fWgtPP_Cp;//Variation C positive error for PP
  TF1* fWgtPP_Dp;//Variation D positive error for PP
  TF1* fWgtPP_Am;//Variatino A negative error for PP
  TF1* fWgtPP_Bm;//Variation B negative error for PP
  TF1* fWgtPP_Cm;//Variation C negative error for PP
  TF1* fWgtPP_Dm;//Variation D negative error for PP
  TF1* fWgtPA;   //Nominal
  TF1* fWgtPA_Ap;//Variation A positive error for PA 
  TF1* fWgtPA_Bp;//Variation B positive error for PA
  TF1* fWgtPA_Cp;//Variation C positive error for PA
  TF1* fWgtPA_Dp;//Variation D positive error for PA
  TF1* fWgtPA_Am;//Variatino A negative error for PA
  TF1* fWgtPA_Bm;//Variation B negative error for PA
  TF1* fWgtPA_Cm;//Variation C negative error for PA
  TF1* fWgtPA_Dm;//Variation D negative error for PA

  if(Cat_==0){
    fWgtPP = new TF1("fWgtPP","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Cp = new TF1("fWgtPP_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Dp = new TF1("fWgtPP_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Cm = new TF1("fWgtPP_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP_Dm = new TF1("fWgtPP_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPP->SetParameters(    213.702, 5.33983, 25.9812, -4.57806  );
    fWgtPP_Ap->SetParameters( 278.673, 5.33983, 25.9812, -4.57806  );
    fWgtPP_Bp->SetParameters( 213.702, 14.9208, 25.9812, -4.57806  );
    fWgtPP_Cp->SetParameters( 213.702, 5.33983, 28.0564, -4.57806  );
    fWgtPP_Dp->SetParameters( 213.702, 5.33983, 25.9812, -4.03882  );
    fWgtPP_Am->SetParameters( 148.73,  5.33983, 25.9812, -4.57806  );
    fWgtPP_Bm->SetParameters( 213.702,-4.24115, 25.9812, -4.57806  );
    fWgtPP_Cm->SetParameters( 213.702, 5.33983, 23.9059, -4.57806  );
    fWgtPP_Dm->SetParameters( 213.702, 5.33983, 25.9812, -5.11731  );

    fWgtPA = new TF1("fWgtPA","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Cp = new TF1("fWgtPA_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Dp = new TF1("fWgtPA_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Cm = new TF1("fWgtPA_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Dm = new TF1("fWgtPA_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA->SetParameters(    241.832, -59.6455, 47.8526, -5.36701  );
    fWgtPA_Ap->SetParameters( 305.423, -59.6455, 47.8526, -5.36701  );
    fWgtPA_Bp->SetParameters( 241.832, -29.9962, 47.8526, -5.36701  );
    fWgtPA_Cp->SetParameters( 241.832, -59.6455, 52.7035, -5.36701  );
    fWgtPA_Dp->SetParameters( 241.832, -59.6455, 47.8526, -4.50973  );
    fWgtPA_Am->SetParameters( 178.241, -59.6455, 47.8526, -5.36701  );
    fWgtPA_Bm->SetParameters( 241.832, -89.2948, 47.8526, -5.36701  );
    fWgtPA_Cm->SetParameters( 241.832, -59.6455, 43.0018, -5.36701  );
    fWgtPA_Dm->SetParameters( 241.832, -59.6455, 47.8526, -6.2243   );
  }
  else if(Cat_==1){
    fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
    fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
    fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
    fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
    fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
    fWgtPP->SetParameters(    0.572979, 0.0694361  );
    fWgtPP_Ap->SetParameters( 0.596206, 0.0694361  );
    fWgtPP_Bp->SetParameters( 0.572979, 0.0727404  );
    fWgtPP_Am->SetParameters( 0.549751, 0.0694361  );
    fWgtPP_Bm->SetParameters( 0.572979, 0.0661317  );
    fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
    //fWgtPA->SetParameters(    0.713383, 0.0389179  );
    //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
    //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
    //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
    //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
    fWgtPA->SetParameters(    0.424411, 0.0966854  );
    fWgtPA_Ap->SetParameters( 0.487855, 0.0966854  );
    fWgtPA_Bp->SetParameters( 0.424411, 0.10523    );
    fWgtPA_Am->SetParameters( 0.360966, 0.0966854  );
    fWgtPA_Bm->SetParameters( 0.424411, 0.0881409  );
  }
  else if(Cat_==2){
    fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
    fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
    fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
    fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
    fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
    fWgtPP->SetParameters(    0.572979, 0.0694361  );
    fWgtPP_Ap->SetParameters( 0.596206, 0.0694361  );
    fWgtPP_Bp->SetParameters( 0.572979, 0.0727404  );
    fWgtPP_Am->SetParameters( 0.549751, 0.0694361  );
    fWgtPP_Bm->SetParameters( 0.572979, 0.0661317  );
    fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
    //fWgtPA->SetParameters(    0.713383, 0.0389179  );
    //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
    //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
    //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
    //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
    fWgtPA->SetParameters(    0.698849, 0.0322766  );
    fWgtPA_Ap->SetParameters( 0.854923, 0.0322766  );
    fWgtPA_Bp->SetParameters( 0.698849, 0.0443489  );
    fWgtPA_Am->SetParameters( 0.542774, 0.0322766  );
    fWgtPA_Bm->SetParameters( 0.698849, 0.0202042  );
  }

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample
  int nPtBins = 1;  double* ptBin = NULL;  int nYBins = 1;  double *yBin = NULL; int nYBinsPP = 1; double *yBinPP = NULL;
  const int nPtBins1 = 1;   double ptBin1[nPtBins1+1] = {0,30};             const int nYBins1 = 1;   double yBin1[nYBins1+1] = {-1.93,1.93}; // integrated
  const int nPtBins2 = 6;   double ptBin2[nPtBins2+1] = {0,2,4,6,9,12,30};  const int nYBins2 = 1;   double yBin2[nYBins2+1] = {-1.93,1.93}; // 1S
  //const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,30};             const int nYBins3 = 10;   double yBin3[nYBins3+1] = {-1.93, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 1.93}; // 1S
  const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,30};             const int nYBins3 = 8;   double yBin3[nYBins3+1] = {-1.93, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.93}; // 1S
  const int nPtBins4 = 3;   double ptBin4[nPtBins4+1] = {0,4,9,30};         const int nYBins4 = 1;   double yBin4[nYBins4+1] = {-1.93,1.93}; // 2S
  const int nPtBins5 = 1;   double ptBin5[nPtBins5+1] = {0,30};             const int nYBins5 = 4;   double yBin5[nYBins5+1] = {-1.93, -0.8, 0, 0.8, 1.93}; // 2S
  const int nPtBins6 = 2;   double ptBin6[nPtBins6+1] = {0,6,30};           const int nYBins6 = 1;   double yBin6[nYBins6+1] = {-1.93, 1.93}; // 3S
  const int nPtBins7 = 1;   double ptBin7[nPtBins7+1] = {0,30};             const int nYBins7 = 2;   double yBin7[nYBins7+1] = {-1.93, 0, 1.93}; // 3S
  //for PP
  const int nYBins8 = 1;   double yBin8[nYBins8+1] = {-1.93, 1.93}; // integrated
  const int nYBins9 = 8;   double yBin9[nYBins9+1] = {-1.93, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.93}; // integrated
  const int nYBins10 = 4;   double yBin10[nYBins10+1] = {-1.93, -0.8, 0, 0.8, 1.93}; // integrated
  const int nYBins11 = 2;   double yBin11[nYBins11+1] = {-1.93, 0, 1.93}; // integrated

  // bin 1 : integrated
  // bin 2, 4 : pT 
  // bin 3, 5 : y

  TFile *in; // input skimed files
  if(Cat_ == 0) in = new TFile("skimedForAcc_MC_Ups1S_20170808.root");
  if(Cat_ == 1) in = new TFile("skimedForAcc_MC_Ups2S_20170808.root");
  if(Cat_ == 2) in = new TFile("skimedForAcc_MC_Ups3S_20170808.root");

  TFile *out; // define output file
  if(Cat_ == 0) out = new TFile("acceptance_wgt_1S_20170817.root","RECREATE");
  if(Cat_ == 1) out = new TFile("acceptance_wgt_2S_20170817.root","RECREATE");
  if(Cat_ == 2) out = new TFile("acceptance_wgt_3S_20170817.root","RECREATE");

  // 1 : denominator, 2 : numerator 
  TH1F *hIntAccPPNoW = new TH1F("hIntAccPPNoW",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPPNoW1 = new TH1F("hIntAccPPNoW1",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPPNoW2 = new TH1F("hIntAccPPNoW2",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPANoW = new TH1F("hIntAccPANoW",";;",nYBins1,yBin1);
  TH1F *hIntAccPANoW1 = new TH1F("hIntAccPANoW1",";;",nYBins1,yBin1);
  TH1F *hIntAccPANoW2 = new TH1F("hIntAccPANoW2",";;",nYBins1,yBin1);
  TH1F *hIntAccPP = new TH1F("hIntAccPP",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPP1 = new TH1F("hIntAccPP1",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPP2 = new TH1F("hIntAccPP2",";|y_{CM}|;",nYBins8,yBin8);
  TH1F *hIntAccPA = new TH1F("hIntAccPA",";;",nYBins1,yBin1);
  TH1F *hIntAccPA1 = new TH1F("hIntAccPA1",";;",nYBins1,yBin1);
  TH1F *hIntAccPA2 = new TH1F("hIntAccPA2",";;",nYBins1,yBin1);
  TH1D* hIntSysPP=new TH1D("hIntSysPP",";|y_{CM}|;",nYBins8,yBin8);
  TH1D* hIntSysPA=new TH1D("hIntSysPA",";y_{CM};",nYBins1,yBin1);

  hIntAccPPNoW->Sumw2();
  hIntAccPPNoW1->Sumw2();
  hIntAccPPNoW2->Sumw2();
  hIntAccPANoW->Sumw2();
  hIntAccPANoW1->Sumw2();
  hIntAccPANoW2->Sumw2();
  hIntAccPP->Sumw2();
  hIntAccPP1->Sumw2();
  hIntAccPP2->Sumw2();
  hIntAccPA->Sumw2();
  hIntAccPA1->Sumw2();
  hIntAccPA2->Sumw2();

  if(Cat_ == 0) {
    nPtBins = nPtBins2; ptBin = ptBin2;
    nYBins = nYBins3; yBin = yBin3;
    nYBinsPP = nYBins9; yBinPP = yBin9;
  }else if(Cat_ == 1){
    nPtBins = nPtBins4; ptBin = ptBin4;
    nYBins = nYBins5; yBin = yBin5;
    nYBinsPP = nYBins10; yBinPP = yBin10;
  }else if(Cat_ == 2){
    nPtBins = nPtBins6; ptBin = ptBin6;
    nYBins = nYBins7; yBin = yBin7;
    nYBinsPP = nYBins11; yBinPP = yBin11;
  }

  TH1D* hptSysPP=new TH1D("hptSysPP",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hrapSysPP=new TH1D("hrapSysPP",";|y_{CM}|;",nYBinsPP,yBinPP);
  TH1D* hptSysPA=new TH1D("hptSysPA",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hrapSysPA=new TH1D("hrapSysPA",";y_{CM};",nYBins,yBin);

  TH1F *hptAccNoW = new TH1F("hptAccNoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPPNoW = new TH1F("hptAccPPNoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPPNoW1 = new TH1F("hptAccPPNoW1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPPNoW2 = new TH1F("hptAccPPNoW2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPANoW = new TH1F("hptAccPANoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPANoW1 = new TH1F("hptAccPANoW1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPANoW2 = new TH1F("hptAccPANoW2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPP = new TH1F("hptAccPP",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPP1 = new TH1F("hptAccPP1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPP2 = new TH1F("hptAccPP2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPA = new TH1F("hptAccPA",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPA1 = new TH1F("hptAccPA1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1F *hptAccPA2 = new TH1F("hptAccPA2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);

  hptAccPPNoW->Sumw2();
  hptAccPPNoW1->Sumw2();
  hptAccPPNoW2->Sumw2();
  hptAccPANoW->Sumw2();
  hptAccPANoW1->Sumw2();
  hptAccPANoW2->Sumw2();
  hptAccPP->Sumw2();
  hptAccPP1->Sumw2();
  hptAccPP2->Sumw2();
  hptAccPA->Sumw2();
  hptAccPA1->Sumw2();
  hptAccPA2->Sumw2();

  TH1F *hptAccPP_Ap1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Bp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Cp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Dp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Am1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Bm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Cm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Dm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Ap2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Bp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Cp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Dp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Am2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Bm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Cm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPP_Dm2 = (TH1F*)hptAccNoW->Clone();

  TH1F *hptAccPA_Ap1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Bp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Cp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Dp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Am1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Bm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Cm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Dm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Ap2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Bp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Cp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Dp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Am2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Bm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Cm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPA_Dm2 = (TH1F*)hptAccNoW->Clone();

  TH1F *hrapAccNoW = new TH1F("hrapAccNoW",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPPNoW = new TH1F("hrapAccPPNoW",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPPNoW1 = new TH1F("hrapAccPPNoW1",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPPNoW2 = new TH1F("hrapAccPPNoW2",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPANoW = new TH1F("hrapAccPANoW",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPANoW1 = new TH1F("hrapAccPANoW1",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPANoW2 = new TH1F("hrapAccPANoW2",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPP = new TH1F("hrapAccPP",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPP1 = new TH1F("hrapAccPP1",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPP2 = new TH1F("hrapAccPP2",";|y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1F *hrapAccPA = new TH1F("hrapAccPA",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPA1 = new TH1F("hrapAccPA1",";y_{CM};Acceptance",nYBins,yBin);
  TH1F *hrapAccPA2 = new TH1F("hrapAccPA2",";y_{CM};Acceptance",nYBins,yBin);

  hrapAccPPNoW->Sumw2();
  hrapAccPPNoW1->Sumw2();
  hrapAccPPNoW2->Sumw2();
  hrapAccPANoW->Sumw2();
  hrapAccPANoW1->Sumw2();
  hrapAccPANoW2->Sumw2();
  hrapAccPP->Sumw2();
  hrapAccPP1->Sumw2();
  hrapAccPP2->Sumw2();
  hrapAccPA->Sumw2();
  hrapAccPA1->Sumw2();
  hrapAccPA2->Sumw2();

  TH1F *hIntAccPP_Ap1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Bp1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Cp1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Dp1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Am1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Bm1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Cm1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Dm1 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Ap2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Bp2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Cp2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Dp2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Am2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Bm2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Cm2 = (TH1F*)hIntAccPPNoW->Clone();
  TH1F *hIntAccPP_Dm2 = (TH1F*)hIntAccPPNoW->Clone();


  TH1F *hrapAccPP_Ap1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Bp1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Cp1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Dp1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Am1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Bm1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Cm1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Dm1 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Ap2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Bp2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Cp2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Dp2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Am2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Bm2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Cm2 = (TH1F*)hrapAccPPNoW->Clone();
  TH1F *hrapAccPP_Dm2 = (TH1F*)hrapAccPPNoW->Clone();

  TH1F *hIntAccPA_Ap1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Bp1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Cp1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Dp1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Am1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Bm1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Cm1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Dm1 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Ap2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Bp2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Cp2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Dp2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Am2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Bm2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Cm2 = (TH1F*)hIntAccPANoW->Clone();
  TH1F *hIntAccPA_Dm2 = (TH1F*)hIntAccPANoW->Clone();

  TH1F *hrapAccPA_Ap1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Bp1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Cp1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Dp1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Am1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Bm1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Cm1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Dm1 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Ap2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Bp2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Cp2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Dp2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Am2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Bm2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Cm2 = (TH1F*)hrapAccPANoW->Clone();
  TH1F *hrapAccPA_Dm2 = (TH1F*)hrapAccPANoW->Clone();

  //getting trees from skimed files
  float mass = 0.0, pt = 0.0, phi = 0.0, y = 0.0, eta = 0.0;
  float pt1 = 0.0, phi1 = 0.0, eta1 = 0.0;
  float pt2 = 0.0, phi2 = 0.0, eta2 = 0.0;

  TBranch *b_mass;
  TBranch *b_y;
  TBranch *b_pt;
  TBranch *b_phi;
  TBranch *b_eta;
  TBranch *b_pt1;
  TBranch *b_phi1;
  TBranch *b_eta1;
  TBranch *b_pt2;
  TBranch *b_phi2;
  TBranch *b_eta2;

  TTree *mm = (TTree*)in->Get("mmGen");
  mm->SetBranchAddress("mass", &mass, &b_mass);
  mm->SetBranchAddress("pt", &pt, &b_pt);
  mm->SetBranchAddress("phi", &phi, &b_phi);
  mm->SetBranchAddress("y", &y, &b_y);
  mm->SetBranchAddress("eta", &eta, &b_eta);
  mm->SetBranchAddress("pt1", &pt1, &b_pt1);
  mm->SetBranchAddress("eta1", &eta1, &b_eta1);
  mm->SetBranchAddress("phi1", &phi1, &b_phi1);
  mm->SetBranchAddress("pt2", &pt2, &b_pt2);
  mm->SetBranchAddress("eta2", &eta2, &b_eta2);
  mm->SetBranchAddress("phi2", &phi2, &b_phi2);

  float boost = 0.47;

  int nEntries = mm->GetEntries();
  //nEntries = 100;

  int cnt1 = 0; // just for counting for denominator for pp
  int cnt2 = 0; // just for counting for numerator for pp
  int cnt3 = 0; // just for counting for denominator for PA
  int cnt4 = 0; // just for counting for numerator for PA

  for(int i = 0; i < nEntries; i++){
    mm->GetEntry(i);
    if (i%300000==0) cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    //cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    if( fabs(y) < 1.93 && pt<30.0 ){
      double ppwgt = fWgtPP->Eval(pt); // apply weighting factor from functions above-defined for pp
      double ppwgtAp = fWgtPP_Ap->Eval(pt);
      double ppwgtBp = fWgtPP_Bp->Eval(pt);
      double ppwgtAm = fWgtPP_Am->Eval(pt);
      double ppwgtBm = fWgtPP_Bm->Eval(pt);
      double ppwgtCp;
      double ppwgtDp;
      double ppwgtCm;
      double ppwgtDm;
      if(Cat_==0){
        ppwgtCp = fWgtPP_Cp->Eval(pt);
        ppwgtDp = fWgtPP_Dp->Eval(pt);
        ppwgtCm = fWgtPP_Cm->Eval(pt);
        ppwgtDm = fWgtPP_Dm->Eval(pt);
      }
      hIntAccPPNoW1->Fill( y );
      hptAccPPNoW1->Fill(pt);
      hrapAccPPNoW1->Fill( y );
      hIntAccPP1->Fill(y,ppwgt);
      hptAccPP1->Fill(pt,ppwgt);
      hrapAccPP1->Fill(y,ppwgt);
      //denominator for systematics with variations
      hIntAccPP_Ap1->Fill( y, ppwgtAp );
      hIntAccPP_Bp1->Fill( y, ppwgtBp );
      hIntAccPP_Am1->Fill( y, ppwgtAm );
      hIntAccPP_Bm1->Fill( y, ppwgtBm );
      hptAccPP_Ap1->Fill( pt, ppwgtAp );
      hptAccPP_Bp1->Fill( pt, ppwgtBp );
      hptAccPP_Am1->Fill( pt, ppwgtAm );
      hptAccPP_Bm1->Fill( pt, ppwgtBm );
      hrapAccPP_Ap1->Fill( y, ppwgtAp );
      hrapAccPP_Bp1->Fill( y, ppwgtBp );
      hrapAccPP_Am1->Fill( y, ppwgtAm );
      hrapAccPP_Bm1->Fill( y, ppwgtBm );
      if(Cat_==0){
        hIntAccPP_Cp1->Fill( y, ppwgtCp );
        hIntAccPP_Dp1->Fill( y, ppwgtDp );
        hIntAccPP_Cm1->Fill( y, ppwgtCm );
        hIntAccPP_Dm1->Fill( y, ppwgtDm );
        hptAccPP_Cp1->Fill( pt, ppwgtCp );
        hptAccPP_Dp1->Fill( pt, ppwgtDp );
        hptAccPP_Cm1->Fill( pt, ppwgtCm );
        hptAccPP_Dm1->Fill( pt, ppwgtDm );
        hrapAccPP_Cp1->Fill( y, ppwgtCp );
        hrapAccPP_Dp1->Fill( y, ppwgtDp );
        hrapAccPP_Cm1->Fill( y, ppwgtCm );
        hrapAccPP_Dm1->Fill( y, ppwgtDm );
      }
      cnt1++;
      if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
        hIntAccPPNoW2->Fill( y );
        hptAccPPNoW2->Fill(pt);
        hrapAccPPNoW2->Fill( y );
        hIntAccPP2->Fill(y,ppwgt);
        hptAccPP2->Fill(pt,ppwgt);
        hrapAccPP2->Fill(y,ppwgt);
        //numerator for systematics with variations 
        //
        hIntAccPP_Ap2->Fill( y, ppwgtAp );
        hIntAccPP_Bp2->Fill( y, ppwgtBp );
        hIntAccPP_Am2->Fill( y, ppwgtAm );
        hIntAccPP_Bm2->Fill( y, ppwgtBm );
        hptAccPP_Ap2->Fill( pt, ppwgtAp );
        hptAccPP_Bp2->Fill( pt, ppwgtBp );
        hptAccPP_Am2->Fill( pt, ppwgtAm );
        hptAccPP_Bm2->Fill( pt, ppwgtBm );
        hrapAccPP_Ap2->Fill( y, ppwgtAp );
        hrapAccPP_Bp2->Fill( y, ppwgtBp );
        hrapAccPP_Am2->Fill( y, ppwgtAm );
        hrapAccPP_Bm2->Fill( y, ppwgtBm );
        if(Cat_==0){
          hIntAccPP_Cp2->Fill( y, ppwgtCp );
          hIntAccPP_Dp2->Fill( y, ppwgtDp );
          hIntAccPP_Cm2->Fill( y, ppwgtCm );
          hIntAccPP_Dm2->Fill( y, ppwgtDm );
          hptAccPP_Cp2->Fill( pt, ppwgtCp );
          hptAccPP_Dp2->Fill( pt, ppwgtDp );
          hptAccPP_Cm2->Fill( pt, ppwgtCm );
          hptAccPP_Dm2->Fill( pt, ppwgtDm );
          hrapAccPP_Cp2->Fill( y, ppwgtCp );
          hrapAccPP_Dp2->Fill( y, ppwgtDp );
          hrapAccPP_Cm2->Fill( y, ppwgtCm );
          hrapAccPP_Dm2->Fill( y, ppwgtDm );
        }
        cnt2++;
      }
    }
    if( y > -1.46 && y < 2.4 && pt<30.0 ){
      double pawgt = fWgtPA->Eval(pt); // apply weighting factor from functions above-defined for PA
      double pawgtAp = fWgtPA_Ap->Eval(pt);
      double pawgtBp = fWgtPA_Bp->Eval(pt);
      double pawgtAm = fWgtPA_Am->Eval(pt);
      double pawgtBm = fWgtPA_Bm->Eval(pt);
      double pawgtCp;
      double pawgtDp;
      double pawgtCm;
      double pawgtDm;
      if(Cat_==0){
        pawgtCp = fWgtPA_Cp->Eval(pt);
        pawgtDp = fWgtPA_Dp->Eval(pt);
        pawgtCm = fWgtPA_Cm->Eval(pt);
        pawgtDm = fWgtPA_Dm->Eval(pt);
      }
      hIntAccPANoW1->Fill(y-boost);
      hptAccPANoW1->Fill(pt);
      hrapAccPANoW1->Fill(y-boost);
      hIntAccPA1->Fill(y-boost,pawgt);
      hptAccPA1->Fill(pt,pawgt);
      hrapAccPA1->Fill(y-boost,pawgt);
      //denominator for systematics with variations
      hIntAccPA_Ap1->Fill( y-boost, pawgtAp );
      hIntAccPA_Bp1->Fill( y-boost, pawgtBp );
      hIntAccPA_Am1->Fill( y-boost, pawgtAm );
      hIntAccPA_Bm1->Fill( y-boost, pawgtBm );
      hptAccPA_Ap1->Fill( pt, pawgtAp );
      hptAccPA_Bp1->Fill( pt, pawgtBp );
      hptAccPA_Am1->Fill( pt, pawgtAm );
      hptAccPA_Bm1->Fill( pt, pawgtBm );
      hrapAccPA_Ap1->Fill( y-boost, pawgtAp );
      hrapAccPA_Bp1->Fill( y-boost, pawgtBp );
      hrapAccPA_Am1->Fill( y-boost, pawgtAm );
      hrapAccPA_Bm1->Fill( y-boost, pawgtBm );
      if(Cat_==0){
        hIntAccPA_Cp1->Fill( y-boost, pawgtCp );
        hIntAccPA_Dp1->Fill( y-boost, pawgtDp );
        hIntAccPA_Cm1->Fill( y-boost, pawgtCm );
        hIntAccPA_Dm1->Fill( y-boost, pawgtDm );
        hptAccPA_Cp1->Fill( pt, pawgtCp );
        hptAccPA_Dp1->Fill( pt, pawgtDp );
        hptAccPA_Cm1->Fill( pt, pawgtCm );
        hptAccPA_Dm1->Fill( pt, pawgtDm );
        hrapAccPA_Cp1->Fill( y-boost, pawgtCp );
        hrapAccPA_Dp1->Fill( y-boost, pawgtDp );
        hrapAccPA_Cm1->Fill( y-boost, pawgtCm );
        hrapAccPA_Dm1->Fill( y-boost, pawgtDm );
      }
      cnt3++;
      if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
        hIntAccPANoW2->Fill(y-boost);
        hptAccPANoW2->Fill(pt);
        hrapAccPANoW2->Fill(y-boost);
        hIntAccPA2->Fill(y-boost,pawgt);
        hptAccPA2->Fill(pt,pawgt);
        hrapAccPA2->Fill(y-boost,pawgt);
        //numerator for systematics with variations 
        hIntAccPA_Ap2->Fill( y-boost, pawgtAp );
        hIntAccPA_Bp2->Fill( y-boost, pawgtBp );
        hIntAccPA_Am2->Fill( y-boost, pawgtAm );
        hIntAccPA_Bm2->Fill( y-boost, pawgtBm );
        hptAccPA_Ap2->Fill( pt, pawgtAp );
        hptAccPA_Bp2->Fill( pt, pawgtBp );
        hptAccPA_Am2->Fill( pt, pawgtAm );
        hptAccPA_Bm2->Fill( pt, pawgtBm );
        hrapAccPA_Ap2->Fill( y-boost, pawgtAp );
        hrapAccPA_Bp2->Fill( y-boost, pawgtBp );
        hrapAccPA_Am2->Fill( y-boost, pawgtAm );
        hrapAccPA_Bm2->Fill( y-boost, pawgtBm );
        if(Cat_==0){
          hIntAccPA_Cp2->Fill( y-boost, pawgtCp );
          hIntAccPA_Dp2->Fill( y-boost, pawgtDp );
          hIntAccPA_Cm2->Fill( y-boost, pawgtCm );
          hIntAccPA_Dm2->Fill( y-boost, pawgtDm );
          hptAccPA_Cp2->Fill( pt, pawgtCp );
          hptAccPA_Dp2->Fill( pt, pawgtDp );
          hptAccPA_Cm2->Fill( pt, pawgtCm );
          hptAccPA_Dm2->Fill( pt, pawgtDm );
          hrapAccPA_Cp2->Fill( y-boost, pawgtCp );
          hrapAccPA_Dp2->Fill( y-boost, pawgtDp );
          hrapAccPA_Cm2->Fill( y-boost, pawgtCm );
          hrapAccPA_Dm2->Fill( y-boost, pawgtDm );
        }
        cnt4++;
      }
    }
  }
  cout<<"cnt1 : "<<cnt1<<", cnt2 : "<<cnt2<<", cnt3 : "<<cnt3<<", cnt4 : "<<cnt4<<endl;
  //cout<<"\n Acceptance No W.F   : "<<hptAccPANoW2->GetBinContent(1)/hptAccPANoW1->GetBinContent(1)<<", "<<hptAccPANoW2->GetBinContent(2)/hptAccPANoW1->GetBinContent(2)<<
  //      "\n Acceptance Nominal  : "<<hptAccPA2->GetBinContent(1)/hptAccPA1->GetBinContent(1)<<", "<<hptAccPA2->GetBinContent(2)/hptAccPA1->GetBinContent(2)<<
  //      "\n Acceptance varied A : "<<hptAccPA_Ap2->GetBinContent(1)/hptAccPA_Ap1->GetBinContent(1)<<", "<<hptAccPA_Ap2->GetBinContent(2)/hptAccPA_Ap1->GetBinContent(2)<<
  //      "\n Acceptance varied B : "<<hptAccPA_Bp2->GetBinContent(1)/hptAccPA_Bp1->GetBinContent(1)<<", "<<hptAccPA_Bp2->GetBinContent(2)/hptAccPA_Bp1->GetBinContent(2)<<
  //      "\n Acceptance varied C : "<<hptAccPA_Cp2->GetBinContent(1)/hptAccPA_Cp1->GetBinContent(1)<<", "<<hptAccPA_Cp2->GetBinContent(2)/hptAccPA_Cp1->GetBinContent(2)<<
  //      "\n Acceptance varied D : "<<hptAccPA_Dp2->GetBinContent(1)/hptAccPA_Dp1->GetBinContent(1)<<", "<<hptAccPA_Dp2->GetBinContent(2)/hptAccPA_Dp1->GetBinContent(2)<<endl;
  //cout<<"\n Accerapance Nominal  : "<<hrapAccPA2->GetBinContent(1)/hrapAccPA1->GetBinContent(1)<<
  //      "\n Accerapance varied A : "<<hrapAccPA_Ap2->GetBinContent(1)/hrapAccPA_Ap1->GetBinContent(1)<<
  //      "\n Accerapance varied B : "<<hrapAccPA_Bp2->GetBinContent(1)/hrapAccPA_Bp1->GetBinContent(1)<<
  //      "\n Accerapance varied C : "<<hrapAccPA_Cp2->GetBinContent(1)/hrapAccPA_Cp1->GetBinContent(1)<<
  //      "\n Accerapance varied D : "<<hrapAccPA_Dp2->GetBinContent(1)/hrapAccPA_Dp1->GetBinContent(1)<<endl;

  //hptAccPA2->Divide(hptAccPA1);
  //hptAccPA_Ap2->Divide(hptAccPA_Ap1);
  //hptAccPA_Bp2->Divide(hptAccPA_Bp1);
  //hptAccPA_Cp2->Divide(hptAccPA_Cp1);
  //hptAccPA_Dp2->Divide(hptAccPA_Dp1);

  //hrapAccPA2->Divide(hrapAccPA1);
  //hrapAccPA_Ap2->Divide(hrapAccPA_Ap1);
  //hrapAccPA_Bp2->Divide(hrapAccPA_Bp1);
  //hrapAccPA_Cp2->Divide(hrapAccPA_Cp1);
  //hrapAccPA_Dp2->Divide(hrapAccPA_Dp1);

  //hrapAccPA2->Draw();
  //hrapAccPA_Ap2->Draw("same");
  //hrapAccPA_Bp2->Draw("same");
  //hrapAccPA_Cp2->Draw("same");
  //hrapAccPA_Dp2->Draw("same");

  //hptAccPA2->Draw();
  //hptAccPA_Ap2->Draw("same");
  //hptAccPA_Bp2->Draw("same");
  //hptAccPA_Cp2->Draw("same");
  //hptAccPA_Dp2->Draw("same");

  // Cloning for dividing
  TH1F * hIntAccPPNoW2Cp = (TH1F*)hIntAccPPNoW2->Clone();
  TH1F * hptAccPPNoW2Cp = (TH1F*)hptAccPPNoW2->Clone();
  TH1F * hrapAccPPNoW2Cp = (TH1F*)hrapAccPPNoW2->Clone();
  TH1F * hIntAccPP2Cp = (TH1F*)hIntAccPP2->Clone();
  TH1F * hptAccPP2Cp = (TH1F*)hptAccPP2->Clone();
  TH1F * hrapAccPP2Cp = (TH1F*)hrapAccPP2->Clone();
  TH1F * hIntAccPANoW2Cp = (TH1F*)hIntAccPANoW2->Clone();
  TH1F * hptAccPANoW2Cp = (TH1F*)hptAccPANoW2->Clone();
  TH1F * hrapAccPANoW2Cp = (TH1F*)hrapAccPANoW2->Clone();
  TH1F * hIntAccPA2Cp = (TH1F*)hIntAccPA2->Clone();
  TH1F * hptAccPA2Cp = (TH1F*)hptAccPA2->Clone();
  TH1F * hrapAccPA2Cp = (TH1F*)hrapAccPA2->Clone();
  hIntAccPPNoW2Cp->Divide(hIntAccPPNoW1);
  hptAccPPNoW2Cp->Divide(hptAccPPNoW1);
  hrapAccPPNoW2Cp->Divide(hrapAccPPNoW1);
  hIntAccPP2Cp->Divide(hIntAccPP1);
  hptAccPP2Cp->Divide(hptAccPP1);
  hrapAccPP2Cp->Divide(hrapAccPP1);
  hIntAccPANoW2Cp->Divide(hIntAccPANoW1);
  hptAccPANoW2Cp->Divide(hptAccPANoW1);
  hrapAccPANoW2Cp->Divide(hrapAccPANoW1);
  hIntAccPA2Cp->Divide(hIntAccPA1);
  hptAccPA2Cp->Divide(hptAccPA1);
  hrapAccPA2Cp->Divide(hrapAccPA1);
  //// for systematics with 20 % up and down
  hIntAccPP_Ap2->Divide(hIntAccPP_Ap1);
  hIntAccPP_Bp2->Divide(hIntAccPP_Bp1);
  hptAccPP_Ap2->Divide(hptAccPP_Ap1);
  hptAccPP_Bp2->Divide(hptAccPP_Bp1);
  hrapAccPP_Ap2->Divide(hrapAccPP_Ap1);
  hrapAccPP_Bp2->Divide(hrapAccPP_Bp1);
  hIntAccPP_Am2->Divide(hIntAccPP_Am1);
  hIntAccPP_Bm2->Divide(hIntAccPP_Bm1);
  hptAccPP_Am2->Divide(hptAccPP_Am1);
  hptAccPP_Bm2->Divide(hptAccPP_Bm1);
  hrapAccPP_Am2->Divide(hrapAccPP_Am1);
  hrapAccPP_Bm2->Divide(hrapAccPP_Bm1);
  if(Cat_==0){
    hIntAccPP_Cp2->Divide(hIntAccPP_Cp1);
    hIntAccPP_Dp2->Divide(hIntAccPP_Dp1);
    hptAccPP_Cp2->Divide(hptAccPP_Cp1);
    hptAccPP_Dp2->Divide(hptAccPP_Dp1);
    hrapAccPP_Cp2->Divide(hrapAccPP_Cp1);
    hrapAccPP_Dp2->Divide(hrapAccPP_Dp1);
    hIntAccPP_Cm2->Divide(hIntAccPP_Cm1);
    hIntAccPP_Dm2->Divide(hIntAccPP_Dm1);
    hptAccPP_Cm2->Divide(hptAccPP_Cm1);
    hptAccPP_Dm2->Divide(hptAccPP_Dm1);
    hrapAccPP_Cm2->Divide(hrapAccPP_Cm1);
    hrapAccPP_Dm2->Divide(hrapAccPP_Dm1);
  }
  hIntAccPA_Ap2->Divide(hIntAccPA_Ap1);
  hIntAccPA_Bp2->Divide(hIntAccPA_Bp1);
  hptAccPA_Ap2->Divide(hptAccPA_Ap1);
  hptAccPA_Bp2->Divide(hptAccPA_Bp1);
  hrapAccPA_Ap2->Divide(hrapAccPA_Ap1);
  hrapAccPA_Bp2->Divide(hrapAccPA_Bp1);
  hIntAccPA_Am2->Divide(hIntAccPA_Am1);
  hIntAccPA_Bm2->Divide(hIntAccPA_Bm1);
  hptAccPA_Am2->Divide(hptAccPA_Am1);
  hptAccPA_Bm2->Divide(hptAccPA_Bm1);
  hrapAccPA_Am2->Divide(hrapAccPA_Am1);
  hrapAccPA_Bm2->Divide(hrapAccPA_Bm1);
  if(Cat_==0){
    hIntAccPA_Cp2->Divide(hIntAccPA_Cp1);
    hIntAccPA_Dp2->Divide(hIntAccPA_Dp1);
    hptAccPA_Cp2->Divide(hptAccPA_Cp1);
    hptAccPA_Dp2->Divide(hptAccPA_Dp1);
    hrapAccPA_Cp2->Divide(hrapAccPA_Cp1);
    hrapAccPA_Dp2->Divide(hrapAccPA_Dp1);
    hIntAccPA_Cm2->Divide(hIntAccPA_Cm1);
    hIntAccPA_Dm2->Divide(hIntAccPA_Dm1);
    hptAccPA_Cm2->Divide(hptAccPA_Cm1);
    hptAccPA_Dm2->Divide(hptAccPA_Dm1);
    hrapAccPA_Cm2->Divide(hrapAccPA_Cm1);
    hrapAccPA_Dm2->Divide(hrapAccPA_Dm1);
  }
  hIntAccPP_Ap2->SetName(  "hIntAccPP_Sys_Ap"  );
  hIntAccPP_Bp2->SetName(  "hIntAccPP_Sys_Bp"  );
  hptAccPP_Ap2->SetName(    "hptAccPP_Sys_Ap"  );
  hptAccPP_Bp2->SetName(    "hptAccPP_Sys_Bp"  );
  hrapAccPP_Ap2->SetName(  "hrapAccPP_Sys_Ap"  );
  hrapAccPP_Bp2->SetName(  "hrapAccPP_Sys_Bp"  );
  hIntAccPP_Am2->SetName(  "hIntAccPP_Sys_Am"  );
  hIntAccPP_Bm2->SetName(  "hIntAccPP_Sys_Bm"  );
  hptAccPP_Am2->SetName(    "hptAccPP_Sys_Am"  );
  hptAccPP_Bm2->SetName(    "hptAccPP_Sys_Bm"  );
  hrapAccPP_Am2->SetName(  "hrapAccPP_Sys_Am"  );
  hrapAccPP_Bm2->SetName(  "hrapAccPP_Sys_Bm"  );
  if(Cat_==0){
    hIntAccPP_Cp2->SetName(  "hIntAccPP_Sys_Cp"  );
    hIntAccPP_Dp2->SetName(  "hIntAccPP_Sys_Dp"  );
    hptAccPP_Cp2->SetName(    "hptAccPP_Sys_Cp"  );
    hptAccPP_Dp2->SetName(    "hptAccPP_Sys_Dp"  );
    hrapAccPP_Cp2->SetName(  "hrapAccPP_Sys_Cp"  );
    hrapAccPP_Dp2->SetName(  "hrapAccPP_Sys_Dp"  );
    hIntAccPP_Cm2->SetName(  "hIntAccPP_Sys_Cm"  );
    hIntAccPP_Dm2->SetName(  "hIntAccPP_Sys_Dm"  );
    hptAccPP_Cm2->SetName(    "hptAccPP_Sys_Cm"  );
    hptAccPP_Dm2->SetName(    "hptAccPP_Sys_Dm"  );
    hrapAccPP_Cm2->SetName(  "hrapAccPP_Sys_Cm"  );
    hrapAccPP_Dm2->SetName(  "hrapAccPP_Sys_Dm"  );
  }
  hIntAccPA_Ap2->SetName(  "hIntAccPA_Sys_Ap"  );
  hIntAccPA_Bp2->SetName(  "hIntAccPA_Sys_Bp"  );
  hptAccPA_Ap2->SetName(    "hptAccPA_Sys_Ap"  );
  hptAccPA_Bp2->SetName(    "hptAccPA_Sys_Bp"  );
  hrapAccPA_Ap2->SetName(  "hrapAccPA_Sys_Ap"  );
  hrapAccPA_Bp2->SetName(  "hrapAccPA_Sys_Bp"  );
  hIntAccPA_Am2->SetName(  "hIntAccPA_Sys_Am"  );
  hIntAccPA_Bm2->SetName(  "hIntAccPA_Sys_Bm"  );
  hptAccPA_Am2->SetName(    "hptAccPA_Sys_Am"  );
  hptAccPA_Bm2->SetName(    "hptAccPA_Sys_Bm"  );
  hrapAccPA_Am2->SetName(  "hrapAccPA_Sys_Am"  );
  hrapAccPA_Bm2->SetName(  "hrapAccPA_Sys_Bm"  );
  if(Cat_==0){
    hIntAccPA_Cp2->SetName(  "hIntAccPA_Sys_Cp"  );
    hIntAccPA_Dp2->SetName(  "hIntAccPA_Sys_Dp"  );
    hptAccPA_Cp2->SetName(    "hptAccPA_Sys_Cp"  );
    hptAccPA_Dp2->SetName(    "hptAccPA_Sys_Dp"  );
    hrapAccPA_Cp2->SetName(  "hrapAccPA_Sys_Cp"  );
    hrapAccPA_Dp2->SetName(  "hrapAccPA_Sys_Dp"  );
    hIntAccPA_Cm2->SetName(  "hIntAccPA_Sys_Cm"  );
    hIntAccPA_Dm2->SetName(  "hIntAccPA_Sys_Dm"  );
    hptAccPA_Cm2->SetName(    "hptAccPA_Sys_Cm"  );
    hptAccPA_Dm2->SetName(    "hptAccPA_Sys_Dm"  );
    hrapAccPA_Cm2->SetName(  "hrapAccPA_Sys_Cm"  );
    hrapAccPA_Dm2->SetName(  "hrapAccPA_Sys_Dm"  );
  }

  TCanvas *c1 = new TCanvas("c1","",1200,400);
  c1->Divide(4,1);
  c1->cd(1);
  hPadPt->Draw();
  hptAccPPNoW2Cp->Draw("same E");
  c1->cd(2);
  hPadPt->Draw();
  hptAccPP2Cp->Draw("same E");
  c1->cd(3);
  hPadPt->Draw();
  hptAccPANoW2Cp->Draw("same E");
  c1->cd(4);
  hPadPt->Draw();
  hptAccPA2Cp->Draw("same E");

  if(Cat_ == 0) {c1->SaveAs("plot_acc_1S_pt.png");c1->SaveAs("plot_acc_1S_pt.pdf");}
  if(Cat_ == 1) {c1->SaveAs("plot_acc_2S_pt.png");c1->SaveAs("plot_acc_2S_pt.pdf");}
  if(Cat_ == 2) {c1->SaveAs("plot_acc_3S_pt.png");c1->SaveAs("plot_acc_3S_pt.pdf");}

  c1->cd(1);
  hPadY->Draw();
  hrapAccPPNoW2Cp->Draw("same E");
  c1->cd(2);
  hPadY->Draw();
  hrapAccPP2Cp->Draw("same E");
  c1->cd(3);
  hPadY->Draw();
  hrapAccPANoW2Cp->Draw("same E");
  c1->cd(4);
  hPadY->Draw();
  hrapAccPA2Cp->Draw("same E");

  if(Cat_ == 0) {c1->SaveAs("plot_acc_1S_rap.png");c1->SaveAs("plot_acc_1S_rap.pdf");}
  if(Cat_ == 1) {c1->SaveAs("plot_acc_2S_rap.png");c1->SaveAs("plot_acc_2S_rap.pdf");}
  if(Cat_ == 2) {c1->SaveAs("plot_acc_3S_rap.png");c1->SaveAs("plot_acc_3S_rap.pdf");}

  TF1 *fline1 = new TF1("fline1","1",0.0,30);
  fline1->SetLineColor(kRed+1);
  fline1->SetLineStyle(2);
  fline1->SetLineWidth(2);
  TF1 *fline2 = new TF1("fline2","1",-1.93,1.93);
  fline2->SetLineColor(kRed+1);
  fline2->SetLineStyle(2);
  fline2->SetLineWidth(2);
  TF1 *fline3 = new TF1("fline3","1",0,1.93);
  fline3->SetLineColor(kRed+1);
  fline3->SetLineStyle(2);
  fline3->SetLineWidth(2);

  // Ratio
  TH1F *hptAccPPNoWRat = (TH1F*)hptAccPPNoW2Cp->Clone();
  TH1F *hptAccPANoWRat = (TH1F*)hptAccPANoW2Cp->Clone();
  TH1F *hptAccPPRat = (TH1F*)hptAccPP2Cp->Clone();
  TH1F *hptAccPARat = (TH1F*)hptAccPA2Cp->Clone();
  hptAccPPRat->Divide(hptAccPPNoWRat);
  hptAccPARat->Divide(hptAccPANoWRat);

  TH1F *hrapAccPPNoWRat = (TH1F*)hrapAccPPNoW2Cp->Clone();
  TH1F *hrapAccPANoWRat = (TH1F*)hrapAccPANoW2Cp->Clone();
  TH1F *hrapAccPPRat = (TH1F*)hrapAccPP2Cp->Clone();
  TH1F *hrapAccPARat = (TH1F*)hrapAccPA2Cp->Clone();
  hrapAccPPRat->Divide(hrapAccPPNoWRat);
  hrapAccPARat->Divide(hrapAccPANoWRat);

  TCanvas *c2 = new TCanvas("c2","",800,400);
  c2->Divide(2,1);

  hPadPt1->GetYaxis()->SetRangeUser(0.8,1.2);
  hPadY1->GetYaxis()->SetRangeUser(0.8,1.2);

  c2->cd(1);
  hPadPt1->Draw();
  fline1->Draw("same");
  hptAccPPRat->Draw("same");
  c2->cd(2);
  hPadPt1->Draw();
  fline1->Draw("same");
  hptAccPARat->Draw("same");

  if(Cat_ == 0) {c2->SaveAs("plot_acc_ratio_1S_pt.png");c2->SaveAs("plot_acc_ratio_1S_pt.pdf");}
  if(Cat_ == 1) {c2->SaveAs("plot_acc_ratio_2S_pt.png");c2->SaveAs("plot_acc_ratio_2S_pt.pdf");}
  if(Cat_ == 2) {c2->SaveAs("plot_acc_ratio_3S_pt.png");c2->SaveAs("plot_acc_ratio_3S_pt.pdf");}

  c2->cd(1);
  hPadY1->Draw();
  fline2->Draw("same");
  hrapAccPPRat->Draw("same");
  c2->cd(2);
  hPadY1->Draw();
  fline2->Draw("same");
  hrapAccPARat->Draw("same");

  if(Cat_ == 0) {c2->SaveAs("plot_acc_ratio_1S_rap.png");c2->SaveAs("plot_acc_ratio_1S_rap.pdf");}
  if(Cat_ == 1) {c2->SaveAs("plot_acc_ratio_2S_rap.png");c2->SaveAs("plot_acc_ratio_2S_rap.pdf");}
  if(Cat_ == 2) {c2->SaveAs("plot_acc_ratio_3S_rap.png");c2->SaveAs("plot_acc_ratio_3S_rap.pdf");}

  if(Cat_ == 0) cout<<" %%%%% Upsilon 1S %%%%% "<<endl;
  if(Cat_ == 1) cout<<" %%%%% Upsilon 2S %%%%% "<<endl;
  if(Cat_ == 2) cout<<" %%%%% Upsilon 3S %%%%% "<<endl;

  // Acceptance values and systematics
  for(int i = 0; i < hIntAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Acceptance :: Integrated %%%"<<endl;
    cout<<"No Weighted pp :: "<<Form("%0.3f",hIntAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccPPNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp "<<Form("%0.3f",hIntAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccPP2Cp->GetBinError(i+1))<<endl;
  }
  for(int i = 0; i < hIntAccPANoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted PA :: "<<Form("%0.3f",hIntAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccPANoW2Cp->GetBinError(i+1))<<endl;
    cout<<"PA "<<Form("%0.3f",hIntAccPA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccPA2Cp->GetBinError(i+1))<<endl;
  }
  //cout<<""<<endl;
  //for(int i = 0; i < hIntAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
  //  cout<<"%%% Systematics :: Integrated %%%"<<endl;
  //  cout<<"pp "<<Form("%0.3f",fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPPNoW2Cp->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100)<<endl;
  //}
  //for(int i = 0; i < hIntAccPANoW2Cp->GetXaxis()->GetNbins(); i++){
  //  cout<<"PA "<<Form("%0.3f",fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPANoW2Cp->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100)<<endl;
  //}
  cout<<""<<endl;
  cout<<"%%% Pt  %%%"<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted pp :: "<<i<<" :: "<<Form("%0.3f",hptAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPPNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"No Weighted PA :: "<<i<<" :: "<<Form("%0.3f",hptAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPANoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp :: "<<i<<" :: "<<Form("%0.3f",hptAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPP2Cp->GetBinError(i+1))<<endl;
    cout<<"PA :: "<<i<<" :: "<<Form("%0.3f",hptAccPA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPA2Cp->GetBinError(i+1))<<endl;
  }
  //cout<<""<<endl;
  //for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
  //  cout<<"%%% Systematics :: Pt %%%"<<endl;
  //  cout<<"pp "<<Form("%0.3f",fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPPNoW2Cp->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100)<<endl;
  //  cout<<"PA "<<Form("%0.3f",fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPANoW2Cp->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100)<<endl;
  //}
  cout<<""<<endl;
  cout<<"%%% Rapidity  %%%"<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted pp :: "<<i<<" "<<Form("%0.3f",hrapAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPPNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp :: "<<i<<" :: "<<Form("%0.3f",hrapAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPP2Cp->GetBinError(i+1))<<endl;
  }
  for(int i = 0; i < hrapAccPANoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted PA :: "<<i<<" "<<Form("%0.3f",hrapAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPANoW2Cp->GetBinError(i+1))<<endl;
    cout<<"PA :: "<<i<<" :: "<<Form("%0.3f",hrapAccPA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPA2Cp->GetBinError(i+1))<<endl;
  }
  cout<<""<<endl;
  //for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
  //  cout<<"%%% Systematics :: Rapidity %%%"<<endl;
  //  cout<<"pp "<<Form("%0.3f",fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPNoW2Cp->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100)<<endl;
  //  cout<<"PA "<<Form("%0.3f",fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPANoW2Cp->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100)<<endl;
  //}
  // Calculation of Systematics for variation up and down
  cout<<"\n %%% Calculation of Systematics for 20 % up and down %%% "<<endl;
  cout<<"\n %%% Systematics :: Pt Bins %%%"<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Ap2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Am2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys2 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Bp2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Bm2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys5 = max(fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Ap2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100,fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Am2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100);
    double sys6 = max(fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Bp2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100,fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Bm2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100);
    if(Cat_==0){
    double sys3 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Cp2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Cm2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys4 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Dp2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPP_Dm2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys7 = max(fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Cp2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100,fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Cm2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100);
    double sys8 = max(fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Dp2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100,fabs(hptAccPA2Cp->GetBinContent(i+1)-hptAccPA_Dm2->GetBinContent(i+1))/hptAccPA2Cp->GetBinContent(i+1)*100);
    double Sys_pt_PP = max(sys1,max(sys2,max(sys3,sys4)));
    double Sys_pt_PA = max(sys5,max(sys6,max(sys7,sys8)));
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP 3 : "<<Form("%0.3f",sys3)<<", PP 4 : "<<Form("%0.3f",sys4)<<", PP total : "<<Form("%0.3f",abs(Sys_pt_PP))<<endl;
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA 3 : "<<Form("%0.3f",sys7)<<", PA 4 : "<<Form("%0.3f",sys8)<<", PA total : "<<Form("%0.3f",abs(Sys_pt_PA))<<endl;
    hptSysPP->SetBinContent(i+1,Sys_pt_PP);
    hptSysPA->SetBinContent(i+1,Sys_pt_PA);
    }
    else if(Cat_!=0){
    double Sys_pt_PP = max(sys1,sys2);
    double Sys_pt_PA = max(sys5,sys6);
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP total : "<<Form("%0.3f",abs(Sys_pt_PP))<<endl;
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA total : "<<Form("%0.3f",abs(Sys_pt_PA))<<endl;
    hptSysPP->SetBinContent(i+1,Sys_pt_PP);
    hptSysPA->SetBinContent(i+1,Sys_pt_PA);
    }
  }
  /////////////////////////////////////////////////////////////////
  cout<<"\n %%% Systematics :: Rapidity Bins %%%"<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Ap2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Am2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    double sys2 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Bp2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Bm2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    if(Cat_==0){
    double sys3 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Cp2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Cm2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    double sys4 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Dp2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPP_Dm2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    double Sys_rap_PP = max(sys1,max(sys2,max(sys3,sys4)));
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP 3 : "<<Form("%0.3f",sys3)<<", PP 4 : "<<Form("%0.3f",sys4)<<", PP total : "<<Form("%0.3f",abs(Sys_rap_PP))<<endl;
    hrapSysPP->SetBinContent(i+1,Sys_rap_PP);
    }
    else if(Cat_!=0){
    double Sys_rap_PP = max(sys1,sys2);
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP total : "<<Form("%0.3f",abs(Sys_rap_PP))<<endl;
    hrapSysPP->SetBinContent(i+1,Sys_rap_PP);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hrapAccPANoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys5 = max(fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Ap2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100,fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Am2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100);
    double sys6 = max(fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Bp2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100,fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Bm2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100);
    if(Cat_==0){
    double sys7 = max(fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Cp2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100,fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Cm2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100);
    double sys8 = max(fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Dp2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100,fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPA_Dm2->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100);
    double Sys_rap_PA = max(sys5,max(sys6,max(sys7,sys8)));
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA 3 : "<<Form("%0.3f",sys7)<<", PA 4 : "<<Form("%0.3f",sys8)<<", PA total : "<<Form("%0.3f",abs(Sys_rap_PA))<<endl;
    hrapSysPA->SetBinContent(i+1,Sys_rap_PA);
    }
    else if(Cat_!=0){
    double Sys_rap_PA = max(sys5,sys6);
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA total : "<<Form("%0.3f",abs(Sys_rap_PA))<<endl;
    hrapSysPA->SetBinContent(i+1,Sys_rap_PA);
    }
  }
  /////////////////////////////////////////////////////////////////
  cout<<"\n %%% Systematics :: Integrated Bins %%%"<<endl;
  for(int i = 0; i < hIntAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Ap2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100,fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Am2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100);
    double sys2 = max(fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Bp2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100,fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Bm2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100);
    if(Cat_==0){
    double sys3 = max(fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Cp2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100,fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Cm2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100);
    double sys4 = max(fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Dp2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100,fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccPP_Dm2->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100);
    double Sys_Int_PP = max(sys1,max(sys2,max(sys3,sys4)));
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP 3 : "<<Form("%0.3f",sys3)<<", PP 4 : "<<Form("%0.3f",sys4)<<", PP total : "<<Form("%0.3f",abs(Sys_Int_PP))<<endl;
    hIntSysPP->SetBinContent(i+1,Sys_Int_PP);
    }
    else if(Cat_!=0){
    double Sys_Int_PP = max(sys1,sys2);
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP total : "<<Form("%0.3f",abs(Sys_Int_PP))<<endl;
    hIntSysPP->SetBinContent(i+1,Sys_Int_PP);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hIntAccPANoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys5 = max(fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Ap2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100,fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Am2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100);
    double sys6 = max(fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Bp2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100,fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Bm2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100);
    if(Cat_==0){
    double sys7 = max(fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Cp2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100,fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Cm2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100);
    double sys8 = max(fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Dp2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100,fabs(hIntAccPA2Cp->GetBinContent(i+1)-hIntAccPA_Dm2->GetBinContent(i+1))/hIntAccPA2Cp->GetBinContent(i+1)*100);
    double Sys_Int_PA = max(sys5,max(sys6,max(sys7,sys8)));
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA 3 : "<<Form("%0.3f",sys7)<<", PA 4 : "<<Form("%0.3f",sys8)<<", PA total : "<<Form("%0.3f",abs(Sys_Int_PA))<<endl;
    hIntSysPA->SetBinContent(i+1,Sys_Int_PA);
    }
    else if(Cat_!=0){
    double Sys_Int_PA = max(sys5,sys6);
    cout<<"PA 1 : "<<Form("%0.3f",sys5)<<", PA 2 : "<<Form("%0.3f",sys6)<<", PA total : "<<Form("%0.3f",abs(Sys_Int_PA))<<endl;
    hIntSysPA->SetBinContent(i+1,Sys_Int_PA);
    }
  }
  //////////////////////////////////Make Systematics/////////////////////////////////
  TFile *fout = new TFile(Form("sys_acceptance_ups%dS_test.root",Cat_+1),"recreate");

  fout->cd();

  TH1F * hrapSysPP_Clone = (TH1F*)hrapSysPP->Clone();
  hptSysPP->Write();
  if(Cat_==0){
    hrapSysPP_Clone->SetName("hrapSysPP");
    hrapSysPP_Clone->SetBinContent( 1, hrapSysPP->GetBinContent(1) );
    hrapSysPP_Clone->SetBinContent( 2, 4.256719e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 3, 4.627279e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 4, 4.524264e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 5, 4.524264e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 6, 4.627279e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 7, 4.256719e-01 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 8, hrapSysPP->GetBinContent(8) );
    hrapSysPP_Clone->Write();
  }
  if(Cat_==1){
    hrapSysPP_Clone->SetName("hrapSysPP");
    hrapSysPP_Clone->SetBinContent( 1, hrapSysPP->GetBinContent(1) );
    hrapSysPP_Clone->SetBinContent( 2, 4.843599e-02 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 3, 4.843599e-02 ); //Fixed from HIN-16-023
    hrapSysPP_Clone->SetBinContent( 4, hrapSysPP->GetBinContent(4) );
    hrapSysPP_Clone->Write();
  }
  if(Cat_==2){
    hrapSysPP->Write();
  }
  hIntSysPP->Write();
  hptSysPA->Write();
  hrapSysPA->Write();
  hIntSysPA->Write();

  fout->Close();

  out->cd();
  hIntAccPPNoW1->Write();
  hptAccPPNoW1->Write();
  hrapAccPPNoW1->Write();
  hIntAccPP1->Write();
  hptAccPP1->Write();
  hrapAccPP1->Write();
  hIntAccPANoW1->Write();
  hptAccPANoW1->Write();
  hrapAccPANoW1->Write();
  hIntAccPA1->Write();
  hptAccPA1->Write();
  hrapAccPA1->Write();

  hIntAccPPNoW2->Write();
  hptAccPPNoW2->Write();
  hrapAccPPNoW2->Write();
  hIntAccPP2->Write();
  hptAccPP2->Write();
  hrapAccPP2->Write();
  hIntAccPANoW2->Write();
  hptAccPANoW2->Write();
  hrapAccPANoW2->Write();
  hIntAccPA2->Write();
  hptAccPA2->Write();
  hrapAccPA2->Write();

  if(Cat_ == 0) {
    hIntAccPPNoW2Cp->SetName("hIntAccNoW1S_PP");
    hptAccPPNoW2Cp->SetName("hptAccNoW1S_PP");
    hrapAccPPNoW2Cp->SetName("hrapAccNoW1S_PP");
    hIntAccPP2Cp->SetName("hIntAccPP1S");
    hptAccPP2Cp->SetName("hptAccPP1S");
    hrapAccPP2Cp->SetName("hrapAccPP1S");
    hIntAccPANoW2Cp->SetName("hIntAccNoW1S_PA");
    hptAccPANoW2Cp->SetName("hptAccNoW1S_PA");
    hrapAccPANoW2Cp->SetName("hrapAccNoW1S_PA");
    hIntAccPA2Cp->SetName("hIntAccPA1S");
    hptAccPA2Cp->SetName("hptAccPA1S");
    hrapAccPA2Cp->SetName("hrapAccPA1S");
  }

  if(Cat_ == 1) {
    hIntAccPPNoW2Cp->SetName("hIntAccNoW2S_PP");
    hptAccPPNoW2Cp->SetName("hptAccNoW2S_PP");
    hrapAccPPNoW2Cp->SetName("hrapAccNoW2S_PP");
    hIntAccPP2Cp->SetName("hIntAccPP2S");
    hptAccPP2Cp->SetName("hptAccPP2S");
    hrapAccPP2Cp->SetName("hrapAccPP2S");
    hIntAccPANoW2Cp->SetName("hIntAccNoW2S_PA");
    hptAccPANoW2Cp->SetName("hptAccNoW2S_PA");
    hrapAccPANoW2Cp->SetName("hrapAccNoW2S_PA");
    hIntAccPA2Cp->SetName("hIntAccPA2S");
    hptAccPA2Cp->SetName("hptAccPA2S");
    hrapAccPA2Cp->SetName("hrapAccPA2S");
  }

  if(Cat_ == 2) {
    hIntAccPPNoW2Cp->SetName("hIntAccNoW3S_PP");
    hptAccPPNoW2Cp->SetName("hptAccNoW3S_PP");
    hrapAccPPNoW2Cp->SetName("hrapAccNoW3S_PP");
    hIntAccPP2Cp->SetName("hIntAccPP3S");
    hptAccPP2Cp->SetName("hptAccPP3S");
    hrapAccPP2Cp->SetName("hrapAccPP3S");
    hIntAccPANoW2Cp->SetName("hIntAccNoW3S_PA");
    hptAccPANoW2Cp->SetName("hptAccNoW3S_PA");
    hrapAccPANoW2Cp->SetName("hrapAccNoW3S_PA");
    hIntAccPA2Cp->SetName("hIntAccPA3S");
    hptAccPA2Cp->SetName("hptAccPA3S");
    hrapAccPA2Cp->SetName("hrapAccPA3S");
  }

  hIntAccPPNoW2Cp->Write();
  hptAccPPNoW2Cp->Write();
  hrapAccPPNoW2Cp->Write();
  hIntAccPP2Cp->Write();
  hptAccPP2Cp->Write();
  hIntAccPANoW2Cp->Write();
  hptAccPANoW2Cp->Write();
  hrapAccPANoW2Cp->Write();
  hIntAccPA2Cp->Write();
  hptAccPA2Cp->Write();
  hrapAccPA2Cp->Write();
  TH1F * hrapAccPP2Cp_Clone = (TH1F*)hrapAccPP2Cp->Clone();
  if(Cat_==0){
  hrapAccPP2Cp_Clone->SetName("hrapAccPP1S");
  hrapAccPP2Cp_Clone->SetBinContent( 1, hrapAccPP2Cp->GetBinContent(1) );
  hrapAccPP2Cp_Clone->SetBinContent( 2, 2.475569e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 3, 2.480554e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 4, 2.494152e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 5, 2.494152e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 6, 2.480554e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 7, 2.475569e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 8, hrapAccPP2Cp->GetBinContent(8) );
  hrapAccPP2Cp_Clone->Write();
  }
  if(Cat_==1){
  hrapAccPP2Cp_Clone->SetName("hrapAccPP2S");
  hrapAccPP2Cp_Clone->SetBinContent( 1, hrapAccPP2Cp->GetBinContent(1) );
  hrapAccPP2Cp_Clone->SetBinContent( 2, 3.095539e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 3, 3.095539e-01 );
  hrapAccPP2Cp_Clone->SetBinContent( 4, hrapAccPP2Cp->GetBinContent(4) );
  hrapAccPP2Cp_Clone->Write();
  }

  hptAccPPNoWRat->Write();
  hptAccPANoWRat->Write();
  hptAccPPRat->Write();
  hptAccPARat->Write();
  hrapAccPPNoWRat->Write();
  hrapAccPANoWRat->Write();
  hrapAccPPRat->Write();
  hrapAccPARat->Write();

  hIntAccPP_Ap2->Write();
  hIntAccPP_Bp2->Write();
  hptAccPP_Ap2->Write( );
  hptAccPP_Bp2->Write( );
  hrapAccPP_Ap2->Write();
  hrapAccPP_Bp2->Write();
  hIntAccPP_Am2->Write();
  hIntAccPP_Bm2->Write();
  hptAccPP_Am2->Write( );
  hptAccPP_Bm2->Write( );
  hrapAccPP_Am2->Write();
  hrapAccPP_Bm2->Write();
  if(Cat_==0){
    hIntAccPP_Cp2->Write();
    hIntAccPP_Dp2->Write();
    hptAccPP_Cp2->Write( );
    hptAccPP_Dp2->Write( );
    hrapAccPP_Cp2->Write();
    hrapAccPP_Dp2->Write();
    hIntAccPP_Cm2->Write();
    hIntAccPP_Dm2->Write();
    hptAccPP_Cm2->Write( );
    hptAccPP_Dm2->Write( );
    hrapAccPP_Cm2->Write();
    hrapAccPP_Dm2->Write();
  }
  hIntAccPA_Ap2->Write();
  hIntAccPA_Bp2->Write();
  hptAccPA_Ap2->Write( );
  hptAccPA_Bp2->Write( );
  hrapAccPA_Ap2->Write();
  hrapAccPA_Bp2->Write();
  hIntAccPA_Am2->Write();
  hIntAccPA_Bm2->Write();
  hptAccPA_Am2->Write( );
  hptAccPA_Bm2->Write( );
  hrapAccPA_Am2->Write();
  hrapAccPA_Bm2->Write();
  if(Cat_==0){
    hIntAccPA_Cp2->Write();
    hIntAccPA_Dp2->Write();
    hptAccPA_Cp2->Write( );
    hptAccPA_Dp2->Write( );
    hrapAccPA_Cp2->Write();
    hrapAccPA_Dp2->Write();
    hIntAccPA_Cm2->Write();
    hIntAccPA_Dm2->Write();
    hptAccPA_Cm2->Write( );
    hptAccPA_Dm2->Write( );
    hrapAccPA_Cm2->Write();
    hrapAccPA_Dm2->Write();
  }

  out->Write();

}
