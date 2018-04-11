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

TH1D *hPadPt = new TH1D("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1D *hPadPt1 = new TH1D("hPadPt1",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1D *hPadY = new TH1D("hPadY",";y_{CM};Acceptance",10, -1.93, 1.93);
TH1D *hPadY1 = new TH1D("hPadY1",";y_{CM};Acceptance",10, -1.93, 1.93);

TLatex *lt1 = new TLatex();

double getMaxDev(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4);
void ColMar(TH1D* h1, Color_t kBlack, int a=0);

void AccAna(int Cat_ = 0, int date = 20180328) { // Cat_ == 0 (1S),  Cat_ == 1 (2S), Cat_ == 2 (3S)
  gROOT->Macro("~/.rootlogon.C");
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

  //Define variation by parameters error 
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
  TF1* fWgtPAcross;   //Nominal
  TF1* fWgtPAcross_Ap;//Variation A positive error for PAcross 
  TF1* fWgtPAcross_Bp;//Variation B positive error for PAcross
  TF1* fWgtPAcross_Cp;//Variation C positive error for PAcross
  TF1* fWgtPAcross_Dp;//Variation D positive error for PAcross
  TF1* fWgtPAcross_Am;//Variatino A negative error for PAcross
  TF1* fWgtPAcross_Bm;//Variation B negative error for PAcross
  TF1* fWgtPAcross_Cm;//Variation C negative error for PAcross
  TF1* fWgtPAcross_Dm;//Variation D negative error for PAcross
  //OLD Nominal
  //if(Cat_==0){
  //  fWgtPP = new TF1("fWgtPP","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Cp = new TF1("fWgtPP_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Dp = new TF1("fWgtPP_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Cm = new TF1("fWgtPP_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Dm = new TF1("fWgtPP_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP->SetParameters(    213.702, 5.33983, 25.9812, -4.57806  );
  //  fWgtPP_Ap->SetParameters( 278.673, 5.33983, 25.9812, -4.57806  );
  //  fWgtPP_Bp->SetParameters( 213.702, 14.9208, 25.9812, -4.57806  );
  //  fWgtPP_Cp->SetParameters( 213.702, 5.33983, 28.0564, -4.57806  );
  //  fWgtPP_Dp->SetParameters( 213.702, 5.33983, 25.9812, -4.03882  );
  //  fWgtPP_Am->SetParameters( 148.73,  5.33983, 25.9812, -4.57806  );
  //  fWgtPP_Bm->SetParameters( 213.702,-4.24115, 25.9812, -4.57806  );
  //  fWgtPP_Cm->SetParameters( 213.702, 5.33983, 23.9059, -4.57806  );
  //  fWgtPP_Dm->SetParameters( 213.702, 5.33983, 25.9812, -5.11731  );

  //  fWgtPA = new TF1("fWgtPA","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Cp = new TF1("fWgtPA_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Dp = new TF1("fWgtPA_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Cm = new TF1("fWgtPA_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Dm = new TF1("fWgtPA_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA->SetParameters(    241.832, -59.6455, 47.8526, -5.36701  );
  //  fWgtPA_Ap->SetParameters( 305.423, -59.6455, 47.8526, -5.36701  );
  //  fWgtPA_Bp->SetParameters( 241.832, -29.9962, 47.8526, -5.36701  );
  //  fWgtPA_Cp->SetParameters( 241.832, -59.6455, 52.7035, -5.36701  );
  //  fWgtPA_Dp->SetParameters( 241.832, -59.6455, 47.8526, -4.50973  );
  //  fWgtPA_Am->SetParameters( 178.241, -59.6455, 47.8526, -5.36701  );
  //  fWgtPA_Bm->SetParameters( 241.832, -89.2948, 47.8526, -5.36701  );
  //  fWgtPA_Cm->SetParameters( 241.832, -59.6455, 43.0018, -5.36701  );
  //  fWgtPA_Dm->SetParameters( 241.832, -59.6455, 47.8526, -6.2243   );
  //}
  //else if(Cat_==1){
  //  fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
  //  fWgtPP->SetParameters(    0.572979, 0.0694361  );
  //  fWgtPP_Ap->SetParameters( 0.596206, 0.0694361  );
  //  fWgtPP_Bp->SetParameters( 0.572979, 0.0727404  );
  //  fWgtPP_Am->SetParameters( 0.549751, 0.0694361  );
  //  fWgtPP_Bm->SetParameters( 0.572979, 0.0661317  );
  //  fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
  //  //fWgtPA->SetParameters(    0.713383, 0.0389179  );
  //  //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
  //  //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
  //  //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
  //  //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
  //  fWgtPA->SetParameters(    0.424411, 0.0966854  );
  //  fWgtPA_Ap->SetParameters( 0.487855, 0.0966854  );
  //  fWgtPA_Bp->SetParameters( 0.424411, 0.10523    );
  //  fWgtPA_Am->SetParameters( 0.360966, 0.0966854  );
  //  fWgtPA_Bm->SetParameters( 0.424411, 0.0881409  );
  //}
  //else if(Cat_==2){
  //  fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
  //  fWgtPP->SetParameters(    0.572979, 0.0694361  );
  //  fWgtPP_Ap->SetParameters( 0.596206, 0.0694361  );
  //  fWgtPP_Bp->SetParameters( 0.572979, 0.0727404  );
  //  fWgtPP_Am->SetParameters( 0.549751, 0.0694361  );
  //  fWgtPP_Bm->SetParameters( 0.572979, 0.0661317  );
  //  fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
  //  //fWgtPA->SetParameters(    0.713383, 0.0389179  );
  //  //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
  //  //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
  //  //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
  //  //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
  //  fWgtPA->SetParameters(    0.698849, 0.0322766  );
  //  fWgtPA_Ap->SetParameters( 0.854923, 0.0322766  );
  //  fWgtPA_Bp->SetParameters( 0.698849, 0.0443489  );
  //  fWgtPA_Am->SetParameters( 0.542774, 0.0322766  );
  //  fWgtPA_Bm->SetParameters( 0.698849, 0.0202042  );
  //}

  //NEW Nominal
  //if(Cat_==0){
  //  fWgtPP = new TF1("fWgtPP","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Cp = new TF1("fWgtPP_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Dp = new TF1("fWgtPP_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Cm = new TF1("fWgtPP_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP_Dm = new TF1("fWgtPP_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPP->SetParameters(    215.047,  1.91119, 26.8299, -4.58216  );
  //  fWgtPP_Ap->SetParameters( 278.596,  1.91119, 26.8299, -4.58216  );
  //  fWgtPP_Bp->SetParameters( 215.047,  11.3537, 26.8299, -4.58216  );
  //  fWgtPP_Cp->SetParameters( 215.047,  1.91119, 29.1364, -4.58216  );
  //  fWgtPP_Dp->SetParameters( 215.047,  1.91119, 26.8299, -4.05929  );
  //  fWgtPP_Am->SetParameters( 151.499,  1.91119, 26.8299, -4.58216  );
  //  fWgtPP_Bm->SetParameters( 215.047, -7.53128, 26.8299, -4.58216  );
  //  fWgtPP_Cm->SetParameters( 215.047,  1.91119, 24.5235, -4.58216  );
  //  fWgtPP_Dm->SetParameters( 215.047,  1.91119, 26.8299, -5.10502  );

  //  fWgtPA = new TF1("fWgtPA","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Cp = new TF1("fWgtPA_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Dp = new TF1("fWgtPA_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Cm = new TF1("fWgtPA_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA_Dm = new TF1("fWgtPA_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
  //  fWgtPA->SetParameters(    236.326, -48.7544,  46.716, -5.47027  );
  //  fWgtPA_Ap->SetParameters( 308.106, -48.7544,  46.716, -5.47027  );
  //  fWgtPA_Bp->SetParameters( 236.326,  -14.631,  46.716, -5.47027  );
  //  fWgtPA_Cp->SetParameters( 236.326, -48.7544, 52.1208, -5.47027  );
  //  fWgtPA_Dp->SetParameters( 236.326, -48.7544,  46.716, -4.4605   );
  //  fWgtPA_Am->SetParameters( 164.546, -48.7544,  46.716, -5.47027  );
  //  fWgtPA_Bm->SetParameters( 236.326, -82.8779,  46.716, -5.47027  );
  //  fWgtPA_Cm->SetParameters( 236.326, -48.7544, 41.3112, -5.47027  );
  //  fWgtPA_Dm->SetParameters( 236.326, -48.7544,  46.716, -6.48005  );
  //}
  //else if(Cat_==1){
  //  fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
  //  fWgtPP->SetParameters(    0.586834, 0.0715173 );
  //  fWgtPP_Ap->SetParameters( 0.617484, 0.0715173 );
  //  fWgtPP_Bp->SetParameters( 0.586834, 0.0768735 );
  //  fWgtPP_Am->SetParameters( 0.556185, 0.0715173 );
  //  fWgtPP_Bm->SetParameters( 0.586834, 0.0661611 );
  //  fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
  //  //fWgtPA->SetParameters(    0.713383, 0.0389179  );
  //  //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
  //  //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
  //  //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
  //  //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
  //  fWgtPA->SetParameters(    0.445835, 0.0918013 );
  //  fWgtPA_Ap->SetParameters( 0.512795, 0.0918013 );
  //  fWgtPA_Bp->SetParameters( 0.445835, 0.1002    );
  //  fWgtPA_Am->SetParameters( 0.378874, 0.0918013 );
  //  fWgtPA_Bm->SetParameters( 0.445835, 0.0834029 );
  //}
  //else if(Cat_==2){
  //  fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
  //  fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
  //  fWgtPP->SetParameters(    0.390621, 0.110652  );
  //  fWgtPP_Ap->SetParameters( 0.437638, 0.110652  );
  //  fWgtPP_Bp->SetParameters( 0.390621, 0.120475  );
  //  fWgtPP_Am->SetParameters( 0.343603, 0.110652  );
  //  fWgtPP_Bm->SetParameters( 0.390621, 0.100829  );
  //  fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
  //  fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
  //  fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
  //  //fWgtPA->SetParameters(    0.713383, 0.0389179  );
  //  //fWgtPA_Ap->SetParameters( 0.782035, 0.0389179  );
  //  //fWgtPA_Bp->SetParameters( 0.713383, 0.0454362  );
  //  //fWgtPA_Am->SetParameters( 0.644731, 0.0389179  );
  //  //fWgtPA_Bm->SetParameters( 0.713383, 0.0323997  );
  //  fWgtPA->SetParameters(    0.735233, 0.0297727 );
  //  fWgtPA_Ap->SetParameters( 0.897287, 0.0297727 );
  //  fWgtPA_Bp->SetParameters( 0.735233, 0.0415512 );
  //  fWgtPA_Am->SetParameters( 0.573179, 0.0297727 );
  //  fWgtPA_Bm->SetParameters( 0.735233, 0.0179941 );
  //}
  //NEW Nominal 2018_03_20
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
    fWgtPP->SetParameters(     250.913, -2.62827, 27.2277, -4.74677 );
    fWgtPP_Ap->SetParameters(  296.118, -2.62827, 27.2277, -4.74677 );
    fWgtPP_Bp->SetParameters(  250.913, 4.36927, 27.2277, -4.74677  );
    fWgtPP_Cp->SetParameters(  250.913, -2.62827, 28.4258, -4.74677 );
    fWgtPP_Dp->SetParameters(  250.913, -2.62827, 27.2277, -4.40976 );
    fWgtPP_Am->SetParameters(  205.708, -2.62827, 27.2277, -4.74677 );
    fWgtPP_Bm->SetParameters(  250.913, -9.62582, 27.2277, -4.74677 );
    fWgtPP_Cm->SetParameters(  250.913, -2.62827, 26.0297, -4.74677 );
    fWgtPP_Dm->SetParameters(  250.913, -2.62827, 27.2277, -5.08378 );

    fWgtPA = new TF1("fWgtPA","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Cp = new TF1("fWgtPA_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Dp = new TF1("fWgtPA_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Cm = new TF1("fWgtPA_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA_Dm = new TF1("fWgtPA_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPA->SetParameters(     522.141, 46.1862, 62.6194, -8.20691 );
    fWgtPA_Ap->SetParameters(  766.345, 46.1862, 62.6194, -8.20691 );
    fWgtPA_Bp->SetParameters(  522.141, 119.858, 62.6194, -8.20691 );
    fWgtPA_Cp->SetParameters(  522.141, 46.1862, 71.1234, -8.20691 );
    fWgtPA_Dp->SetParameters(  522.141, 46.1862, 62.6194, -6.71281 );
    fWgtPA_Am->SetParameters(  277.937, 46.1862, 62.6194, -8.20691 );
    fWgtPA_Bm->SetParameters(  522.141, -27.486, 62.6194, -8.20691 );
    fWgtPA_Cm->SetParameters(  522.141, 46.1862, 54.1154, -8.20691 );
    fWgtPA_Dm->SetParameters(  522.141, 46.1862, 62.6194, -9.70101 );

    fWgtPAcross = new TF1("fWgtPAcross","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Ap = new TF1("fWgtPAcross_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Bp = new TF1("fWgtPAcross_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Cp = new TF1("fWgtPAcross_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Dp = new TF1("fWgtPAcross_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Am = new TF1("fWgtPAcross_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Bm = new TF1("fWgtPAcross_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Cm = new TF1("fWgtPAcross_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross_Dm = new TF1("fWgtPAcross_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
    fWgtPAcross->SetParameters(    433.169, 24.0137, 58.1574, -7.54821  );
    fWgtPAcross_Ap->SetParameters( 612.963, 24.0137, 58.1574, -7.54821  );
    fWgtPAcross_Bp->SetParameters( 433.169, 81.9038, 58.1574, -7.54821  );
    fWgtPAcross_Cp->SetParameters( 433.169, 24.0137, 65.2252, -7.54821  );
    fWgtPAcross_Dp->SetParameters( 433.169, 24.0137, 58.1574, -6.28171  );
    fWgtPAcross_Am->SetParameters( 253.375, 24.0137, 58.1574, -7.54821  );
    fWgtPAcross_Bm->SetParameters( 433.169, -33.8764, 58.1574, -7.54821 );
    fWgtPAcross_Cm->SetParameters( 433.169, 24.0137, 51.0896, -7.54821  );
    fWgtPAcross_Dm->SetParameters( 433.169, 24.0137, 58.1574, -8.81472  );
  }
  else if(Cat_==1){
    fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
    fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
    fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
    fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
    fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
    fWgtPP->SetParameters(    0.598258, 0.0648555  );
    fWgtPP_Ap->SetParameters( 0.621414, 0.0648555  );
    fWgtPP_Bp->SetParameters( 0.598258, 0.0677232  );
    fWgtPP_Am->SetParameters( 0.575102, 0.0648555  );
    fWgtPP_Bm->SetParameters( 0.598258, 0.0619879  );
    fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
    fWgtPA->SetParameters(    0.405682, 0.0994858 );
    fWgtPA_Ap->SetParameters( 0.468708, 0.0994858 );
    fWgtPA_Bp->SetParameters( 0.405682, 0.108031  );
    fWgtPA_Am->SetParameters( 0.342655, 0.0994858 );
    fWgtPA_Bm->SetParameters( 0.405682, 0.0909409 );
    fWgtPAcross = new TF1("fWgtPAcross",    "( [0] + [1]*x  )",0,30);
    fWgtPAcross_Ap = new TF1("fWgtPAcross_Ap","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Bp = new TF1("fWgtPAcross_Bp","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Am = new TF1("fWgtPAcross_Am","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Bm = new TF1("fWgtPAcross_Bm","( [0] + [1]*x  )",0,30);
    fWgtPAcross->SetParameters(    0.421508, 0.0975407);
    fWgtPAcross_Ap->SetParameters( 0.481827, 0.0975407);
    fWgtPAcross_Bp->SetParameters( 0.421508, 0.105621 );
    fWgtPAcross_Am->SetParameters( 0.361188, 0.0975407);
    fWgtPAcross_Bm->SetParameters( 0.421508, 0.0894604);
  }
  else if(Cat_==2){
    fWgtPP = new TF1("fWgtPP",    "( [0] + [1]*x  )",0,30);
    fWgtPP_Ap = new TF1("fWgtPP_Ap","( [0] + [1]*x  )",0,30);
    fWgtPP_Bp = new TF1("fWgtPP_Bp","( [0] + [1]*x  )",0,30);
    fWgtPP_Am = new TF1("fWgtPP_Am","( [0] + [1]*x  )",0,30);
    fWgtPP_Bm = new TF1("fWgtPP_Bm","( [0] + [1]*x  )",0,30);
    fWgtPP->SetParameters(    0.453252, 0.091353  );
    fWgtPP_Ap->SetParameters( 0.489889, 0.091353  );
    fWgtPP_Bp->SetParameters( 0.453252, 0.0961323 );
    fWgtPP_Am->SetParameters( 0.416615, 0.091353  );
    fWgtPP_Bm->SetParameters( 0.453252, 0.0865736 );
    fWgtPA = new TF1("fWgtPA",    "( [0] + [1]*x  )",0,30);
    fWgtPA_Ap = new TF1("fWgtPA_Ap","( [0] + [1]*x  )",0,30);
    fWgtPA_Bp = new TF1("fWgtPA_Bp","( [0] + [1]*x  )",0,30);
    fWgtPA_Am = new TF1("fWgtPA_Am","( [0] + [1]*x  )",0,30);
    fWgtPA_Bm = new TF1("fWgtPA_Bm","( [0] + [1]*x  )",0,30);
    fWgtPA->SetParameters(    0.664791, 0.0370995 );
    fWgtPA_Ap->SetParameters( 0.831036, 0.0370995 );
    fWgtPA_Bp->SetParameters( 0.664791, 0.0494552 );
    fWgtPA_Am->SetParameters( 0.498545, 0.0370995 );
    fWgtPA_Bm->SetParameters( 0.664791, 0.0247437 );
    fWgtPAcross = new TF1("fWgtPAcross",    "( [0] + [1]*x  )",0,30);
    fWgtPAcross_Ap = new TF1("fWgtPAcross_Ap","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Bp = new TF1("fWgtPAcross_Bp","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Am = new TF1("fWgtPAcross_Am","( [0] + [1]*x  )",0,30);
    fWgtPAcross_Bm = new TF1("fWgtPAcross_Bm","( [0] + [1]*x  )",0,30);
    fWgtPAcross->SetParameters(    0.71427, 0.0319624  );
    fWgtPAcross_Ap->SetParameters( 0.876864, 0.0319624 );
    fWgtPAcross_Bp->SetParameters( 0.71427, 0.0439052  );
    fWgtPAcross_Am->SetParameters( 0.551677, 0.0319624 );
    fWgtPAcross_Bm->SetParameters( 0.71427, 0.0200195  );
  }

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample
  int nPtBins = 1;    double* ptBin = NULL;    
  int nYBins = 1;    double *yBin = NULL;   
  int nYBins3s = 1;    double *yBin3s = NULL;   
  int nYBinsPP = 1;    double *yBinPP = NULL; 
  int nYBins23 = 1;    double *yBin23 = NULL; 
  int nYBinsCrPA = 1; double *yBinCrPA = NULL; 
  int nYBinsCrI = 1;  double *yBinCrI = NULL;
  int nYBinsI   = 1;  double *yBinI = NULL;
  //integrated
  const int nPtBinsi  = 1;  double ptBini[nPtBinsi+1]   = {0, 30};
  const int nYBinsi   = 1;  double yBini[nYBinsi+1]     = {-1.93, 1.93};
  const int nYBinsiPA   = 1;  double yBiniPA[nYBinsiPA+1]     = {-2.87, 1.93};
  //pT
  const int nPtBins1  = 6;  double ptBin1[nPtBins1+1]   = {0, 2, 4, 6, 9, 12, 30};// 1S
  const int nPtBins2  = 3;  double ptBin2[nPtBins2+1]   = {0, 4, 9, 30};          // 2S
  const int nPtBins3  = 2;  double ptBin3[nPtBins3+1]   = {0, 6, 30};             // 3S
  //rapidity
  const int nYBins1   = 8;  double yBin1[nYBins1+1]     = {-1.93, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.93};// 1S
  const int nYBins2   = 4;  double yBin2[nYBins2+1]     = {-1.93, -0.8, 0, 0.8, 1.93};                      // 2S
  const int nYBins3   = 2;  double yBin3[nYBins3+1]     = {-1.93, 0, 1.93};                                 // 3S
  //pA Cross_section
  const int nYBinsCri = 1;  double yBinCri[nYBinsCri+1] = {-2.87, 1.93};                                           // integrated
  const int nYBinsCr1 = 9; double yBinCr1[nYBinsCr1+1] = {-2.87, -1.93, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.93};// 1S
  const int nYBinsCr2 = 5;  double yBinCr2[nYBinsCr2+1] = {-2.87, -1.93, -0.8, 0, 0.8, 1.93};                      // 2S
  const int nYBinsCr3 = 3;  double yBinCr3[nYBinsCr3+1] = {-2.87, -1.93, 0, 1.93};                                 // 3S
  //for pp
  const int nYBinsPPi = 1;  double yBinPPi[nYBinsPPi+1] = {0, 1.93};                // integrated
  const int nYBinsPP1 = 4;  double yBinPP1[nYBinsPP1+1] = {0, 0.4, 0.8, 1.2, 1.93}; // 1S Bin  
  const int nYBinsPP2 = 3;  double yBinPP2[nYBinsPP2+1] = {0, 0.8, 1.93, 2.4};           // 2S Bin
  const int nYBinsPP3 = 2;  double yBinPP3[nYBinsPP3+1] = {0, 1.93, 2.4};                // 3S Bin
  const int nYBins231  = 6; double yBin231[nYBins231+1] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4}; // 1S Bin  
  const int nYBins232  = 3; double yBin232[nYBins232+1] = {0, 0.8, 1.6, 2.4}; // 1S Bin  
  const int nYBins233  = 2; double yBin233[nYBins233+1] = {0, 1.2, 2.4}; // 1S Bin  

  TFile *in; // input skimed files
  if(Cat_ == 0) in = new TFile("skimedForAcc_MC_Ups1S_20170808.root");
  if(Cat_ == 1) in = new TFile("skimedForAcc_MC_Ups2S_20170808.root");
  if(Cat_ == 2) in = new TFile("skimedForAcc_MC_Ups3S_20170808.root");

  TFile *out; // define output file
  if(Cat_ == 0) out = new TFile(Form("%d/acceptance_wgt_1S_%d_2Dplot.root",date, date),"RECREATE");
  if(Cat_ == 1) out = new TFile(Form("%d/acceptance_wgt_2S_%d_2Dplot.root",date, date),"RECREATE");
  if(Cat_ == 2) out = new TFile(Form("%d/acceptance_wgt_3S_%d_2Dplot.root",date, date),"RECREATE");

  // 1 : denominator, 2 : numerator 
  // for RpA
  TH1D *hIntAccPPNoW = new TH1D("hIntAccPPNoW",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPPNoW1 = new TH1D("hIntAccPPNoW1",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPPNoW2 = new TH1D("hIntAccPPNoW2",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPANoW = new TH1D("hIntAccPANoW",";;",nYBinsi,yBini);
  TH1D *hIntAccPANoW1 = new TH1D("hIntAccPANoW1",";;",nYBinsi,yBini);
  TH1D *hIntAccPANoW2 = new TH1D("hIntAccPANoW2",";;",nYBinsi,yBini);
  TH1D *hIntAccPP = new TH1D("hIntAccPP",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPP_1 = new TH1D("hIntAccPP_1",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPP_2 = new TH1D("hIntAccPP_2",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntAccPPPt1_1 = new TH1D("hIntAccPPPt1_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPPPt1_2 = new TH1D("hIntAccPPPt1_2",";;",nYBinsi,yBini);
  TH1D *hIntAccPPPt2_1 = new TH1D("hIntAccPPPt2_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPPPt2_2 = new TH1D("hIntAccPPPt2_2",";;",nYBinsi,yBini);
  TH1D *hIntAccPA = new TH1D("hIntAccPA",";;",nYBinsi,yBini);
  TH1D *hIntAccPA_1 = new TH1D("hIntAccPA_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPA_2 = new TH1D("hIntAccPA_2",";;",nYBinsi,yBini);
  TH1D *hIntAccPAPt1_1 = new TH1D("hIntAccPAPt1_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPAPt1_2 = new TH1D("hIntAccPAPt1_2",";;",nYBinsi,yBini);
  TH1D *hIntAccPAPt2_1 = new TH1D("hIntAccPAPt2_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPAPt2_2 = new TH1D("hIntAccPAPt2_2",";;",nYBinsi,yBini);
  TH1D *hIntAccPANoW_1 = new TH1D("hIntAccPANoW_1",";;",nYBinsi,yBini);
  TH1D *hIntAccPANoW_2 = new TH1D("hIntAccPANoW_2",";;",nYBinsi,yBini);
  TH1D *hIntAccCross_1 = new TH1D("hIntAccCross_1",";;",nYBinsiPA,yBiniPA);
  TH1D *hIntAccCross_2 = new TH1D("hIntAccCross_2",";;",nYBinsiPA,yBiniPA);
  TH1D *hIntSysAccPP=new TH1D("hIntSysAccPP",";y_{CM};",nYBinsi,yBini);
  TH1D *hIntSysAccPA=new TH1D("hIntSysAccPA",";y_{CM};",nYBinsi,yBini);
  //for PA cross section
  TH1D *hIntAccCrossNoW  = new TH1D("hIntAccCrossNoW",";;",nYBinsCri,yBinCri);
  TH1D *hIntAccCrossNoW1 = new TH1D("hIntAccCrossNoW1",";;",nYBinsCri,yBinCri);
  TH1D *hIntAccCrossNoW2 = new TH1D("hIntAccCrossNoW2",";;",nYBinsCri,yBinCri);
  TH1D *hIntAccCrPA     = new TH1D("hIntAccCrPA",";;",nYBinsCri,yBinCri);
  TH1D *hIntAccCrPA1    = new TH1D("hIntAccCrPA1",";;",nYBinsCri,yBinCri);
  TH1D *hIntAccCrPA2    = new TH1D("hIntAccCrPA2",";;",nYBinsCri,yBinCri);

  //2D acceptance
  TH2F *hIntAccPP2D = new TH2F("hIntAccPP2D",";;",144,-2.4,2.4, 120,0,30);
  TH2F *hIntAccPP2D_1 = new TH2F("hIntAccPP2D_1",";y_{CM};p_{T} (GeV/c)",144,-2.4,2.4, 120,0,30);
  TH2F *hIntAccPP2D_2 = new TH2F("hIntAccPP2D_2",";y_{CM};p_{T} (GeV/c)",144,-2.4,2.4, 120,0,30);
  TH2F *hIntAccPA2D = new TH2F("hIntAccPA2D",";;",144,-2.4,2.4, 120,0,30);
  TH2F *hIntAccPA2D_1 = new TH2F("hIntAccPA2D_1",";y_{CM};p_{T} (GeV/c)",144,-2.4,2.4, 120,0,30);
  TH2F *hIntAccPA2D_2 = new TH2F("hIntAccPA2D_2",";y_{CM};p_{T} (GeV/c)",144,-2.4,2.4, 120,0,30);
  hIntAccPP2D->Sumw2();
  hIntAccPP2D_1->Sumw2();
  hIntAccPP2D_2->Sumw2();

  hIntAccPA2D->Sumw2();
  hIntAccPA2D_1->Sumw2();
  hIntAccPA2D_2->Sumw2();

  hIntAccPPNoW->Sumw2();
  hIntAccPPNoW1->Sumw2();
  hIntAccPPNoW2->Sumw2();
  hIntAccPANoW->Sumw2();
  hIntAccPANoW1->Sumw2();
  hIntAccPANoW2->Sumw2();
  hIntAccPP->Sumw2();
  hIntAccPP_1->Sumw2();
  hIntAccPP_2->Sumw2();
  hIntAccPPPt1_1->Sumw2();
  hIntAccPPPt1_2->Sumw2();
  hIntAccPPPt2_1->Sumw2();
  hIntAccPPPt2_2->Sumw2();
  hIntAccPA->Sumw2();
  hIntAccPA_1->Sumw2();
  hIntAccPA_2->Sumw2();
  hIntAccPAPt1_1->Sumw2();
  hIntAccPAPt1_2->Sumw2();
  hIntAccPAPt2_1->Sumw2();
  hIntAccPAPt2_2->Sumw2();
  hIntAccPANoW_1->Sumw2();
  hIntAccPANoW_2->Sumw2();
  hIntAccCross_1->Sumw2();
  hIntAccCross_2->Sumw2();

  hIntAccCrossNoW->Sumw2();
  hIntAccCrossNoW1->Sumw2();
  hIntAccCrossNoW2->Sumw2();
  hIntAccCrPA->Sumw2();
  hIntAccCrPA1->Sumw2();
  hIntAccCrPA2->Sumw2();

  if(Cat_ == 0) {
    nPtBins = nPtBins1; ptBin = ptBin1;
    nYBins = nYBins1; yBin = yBin1;
    nYBinsPP = nYBinsPP1; yBinPP = yBinPP1;
    nYBins23 = nYBins231; yBin23 = yBin231;
    nYBinsCrPA = nYBinsCr1; yBinCrPA = yBinCr1;
    nYBins3s = nYBins3; yBin3s = yBin3;
  }else if(Cat_ == 1){
    nPtBins = nPtBins2; ptBin = ptBin2;
    nYBins = nYBins2; yBin = yBin2;
    nYBinsPP = nYBinsPP2; yBinPP = yBinPP2;
    nYBins23 = nYBins232; yBin23 = yBin232;
    nYBinsCrPA = nYBinsCr2; yBinCrPA = yBinCr2;
    nYBins3s = nYBins3; yBin3s = yBin3;
  }else if(Cat_ == 2){
    nPtBins = nPtBins3; ptBin = ptBin3;
    nYBins = nYBins3; yBin = yBin3;
    nYBinsPP = nYBinsPP3; yBinPP = yBinPP3;
    nYBins23 = nYBins233; yBin23 = yBin233;
    nYBinsCrPA = nYBinsCr3; yBinCrPA = yBinCr3;
    nYBins3s = nYBins3; yBin3s = yBin3;
  }

  //for pp
  TH1D* hptSysAccPP=new TH1D("hptSysAccPP",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptSysAccPPRap1=new TH1D("hptSysAccPPRap1",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptSysAccPPRap2=new TH1D("hptSysAccPPRap2",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hrapSysAccPP=new TH1D("hrapSysAccPP",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPPPt1=new TH1D("hrapSysAccPPPt1",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPPPt2=new TH1D("hrapSysAccPPPt2",";y_{CM};",nYBins,yBin);
  //for pa
  TH1D* hptSysAccCross=new TH1D("hptSysAccCross",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hrapSysAccCross=new TH1D("hrapSysAccCross",";p_{T}(GeV/c);",nYBinsCrPA,yBinCrPA);
  TH1D* hrapSysAccCrossPt1=new TH1D("hrapSysAccCrossPt1",";y_{CM};",nYBinsCrPA,yBinCrPA);
  TH1D* hrapSysAccCrossPt2=new TH1D("hrapSysAccCrossPt2",";y_{CM};",nYBinsCrPA,yBinCrPA);
  //for Rpa
  TH1D* hptSysAccPA=new TH1D("hptSysAccPA",";y_{CM};",nYBins,yBin);
  TH1D* hptSysAccPARap1=new TH1D("hptSysAccPARap1",";y_{CM};",nYBins,yBin);
  TH1D* hptSysAccPARap2=new TH1D("hptSysAccPARap2",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPA=new TH1D("hrapSysAccPA",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPAPt1=new TH1D("hrapSysAccPAPt1",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPAPt2=new TH1D("hrapSysAccPAPt2",";y_{CM};",nYBins,yBin);
  TH1D* hrapSysAccPA2bin=new TH1D("hrapSysAccPA2bin",";y_{CM};",nYBins3s,yBin3s);
  TH1D* hintSysPA=new TH1D("hintSysPA",";y_{CM};",nYBins,yBin);

  TH1D *hptAccNoW = new TH1D("hptAccNoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPNoW = new TH1D("hptAccPPNoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPNoW1 = new TH1D("hptAccPPNoW1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPNoW2 = new TH1D("hptAccPPNoW2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPNoW2Cp = new TH1D("hptAccPPNoW2Cp",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPANoW = new TH1D("hptAccPANoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPANoW1 = new TH1D("hptAccPANoW1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPANoW2 = new TH1D("hptAccPANoW2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPANoW2Cp = new TH1D("hptAccPANoW2Cp",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccCrossNoW = new TH1D("hptAccCrossNoW",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccCrossNoW1 = new TH1D("hptAccCrossNoW1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccCrossNoW2 = new TH1D("hptAccCrossNoW2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap1 = new TH1D("hptAccPPRap1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap2 = new TH1D("hptAccPPRap2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap1_1 = new TH1D("hptAccPPRap1_1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap1_2 = new TH1D("hptAccPPRap1_2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap2_1 = new TH1D("hptAccPPRap2_1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPPRap2_2 = new TH1D("hptAccPPRap2_2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPP1 = new TH1D("hptAccPP1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPP2 = new TH1D("hptAccPP2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPA1 = new TH1D("hptAccPA1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPA2 = new TH1D("hptAccPA2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap1 = new TH1D("hptAccPARap1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap2 = new TH1D("hptAccPARap2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap1_1 = new TH1D("hptAccPARap1_1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap1_2 = new TH1D("hptAccPARap1_2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap2_1 = new TH1D("hptAccPARap2_1",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccPARap2_2 = new TH1D("hptAccPARap2_2",";p_{T} (GeV/c);Acceptance",nPtBins,ptBin);

  TH1D *hptAccCross_1 = new TH1D("hptAccCross_1",";p_{T}(GeV/c);Acceptance",nPtBins,ptBin);
  TH1D *hptAccCross_2 = new TH1D("hptAccCross_2",";p_{T}(GeV/c);Acceptance",nPtBins,ptBin);

  hptAccPPNoW->Sumw2();
  hptAccPPNoW1->Sumw2();
  hptAccPPNoW2->Sumw2();
  hptAccPANoW->Sumw2();
  hptAccPANoW1->Sumw2();
  hptAccPANoW2->Sumw2();
  hptAccPANoW2Cp->Sumw2();
  hptAccCrossNoW->Sumw2();
  hptAccCrossNoW1->Sumw2();
  hptAccCrossNoW2->Sumw2();
  hptAccPPRap1->Sumw2();
  hptAccPPRap2->Sumw2();
  hptAccPPRap1_1->Sumw2();
  hptAccPPRap1_2->Sumw2();
  hptAccPPRap2_1->Sumw2();
  hptAccPPRap2_2->Sumw2();
  hptAccPP1->Sumw2();
  hptAccPP2->Sumw2();
  hptAccPA1->Sumw2();
  hptAccPA2->Sumw2();
  hptAccPARap1->Sumw2();
  hptAccPARap2->Sumw2();
  hptAccPARap1_1->Sumw2();
  hptAccPARap1_2->Sumw2();
  hptAccPARap2_1->Sumw2();
  hptAccPARap2_2->Sumw2();

  TH1D *hptAccPP_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPP_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hptAccPPRap1_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap1_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hptAccPPRap2_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPPRap2_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hptAccPA_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPA_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hptAccPARap1_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap1_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hptAccPARap2_Ap1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Bp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Cp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Dp1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Am1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Bm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Cm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Dm1 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Ap2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Bp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Cp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Dp2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Am2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Bm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Cm2 = (TH1D*)hptAccNoW->Clone();
  TH1D *hptAccPARap2_Dm2 = (TH1D*)hptAccNoW->Clone();

  TH1D *hrapAccNoW = new TH1D("hrapAccNoW",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPNoW = new TH1D("hrapAccPPNoW",";|y_{CM}|;Acceptance",nYBins,yBin);
  TH1D *hrapAccPPNoW1 = new TH1D("hrapAccPPNoW1",";|y_{CM}|;Acceptance",nYBins,yBin);
  TH1D *hrapAccPPNoW2 = new TH1D("hrapAccPPNoW2",";|y_{CM}|;Acceptance",nYBins,yBin);
  TH1D *hrapAccPPNoW2Cp = new TH1D("hrapAccPPNoW2Cp",";|y_{CM}|;Acceptance",nYBins,yBin);
  TH1D *hrapAccPANoW = new TH1D("hrapAccPANoW",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPA2binNoW = new TH1D("hrapAccPA2binNoW",";y_{CM};Acceptance",nYBins3s,yBin3s);
  TH1D *hrapAccCrossNoW2Cp = new TH1D("hrapAccCrossNoW2Cp",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccPANoW1 = new TH1D("hrapAccPANoW1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPANoW2 = new TH1D("hrapAccPANoW2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPANoW2Cp = new TH1D("hrapAccPANoW2Cp",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccCrossNoW = new TH1D("hrapAccCrossNoW",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossNoW1 = new TH1D("hrapAccCrossNoW1",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossNoW2 = new TH1D("hrapAccCrossNoW2",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);

  TH1D *hrapAccPPPt1 = new TH1D("hrapAccPPPt1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPPt2 = new TH1D("hrapAccPPPt2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPPt1_1 = new TH1D("hrapAccPPPt1_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPPt1_2 = new TH1D("hrapAccPPPt1_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPPt2_1 = new TH1D("hrapAccPPPt2_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPPt2_2 = new TH1D("hrapAccPPPt2_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPP_1 = new TH1D("hrapAccPP_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPP_2 = new TH1D("hrapAccPP_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPA_1 = new TH1D("hrapAccPA_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPA_2 = new TH1D("hrapAccPA_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt1 = new TH1D("hrapAccPAPt1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt2 = new TH1D("hrapAccPAPt2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt1_1 = new TH1D("hrapAccPAPt1_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt1_2 = new TH1D("hrapAccPAPt1_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt2_1 = new TH1D("hrapAccPAPt2_1",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPAPt2_2 = new TH1D("hrapAccPAPt2_2",";y_{CM};Acceptance",nYBins,yBin);
  TH1D *hrapAccPPAN_1 = new TH1D("hrapAccPPAN_1",";|Y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1D *hrapAccPPAN_2 = new TH1D("hrapAccPPAN_2",";|Y_{CM}|;Acceptance",nYBinsPP,yBinPP);
  TH1D *hrapAcc23 = new TH1D("hrapAcc23",";|Y_{CM}|;Acceptance",nYBins23,yBin23);

  TH1D *hrapAccPA2bin_1 = new TH1D("hrapAccPA2bin_1",";y_{CM};Acceptance",nYBins3s,yBin3s);
  TH1D *hrapAccPA2bin_2 = new TH1D("hrapAccPA2bin_2",";y_{CM};Acceptance",nYBins3s,yBin3s);
  TH1D *hrapAccCross_1 = new TH1D("hrapAccCross_1",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCross_2 = new TH1D("hrapAccCross_2",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossNoW_1 = new TH1D("hrapAccCrossNoW_1",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossNoW_2 = new TH1D("hrapAccCrossNoW_2",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossPt1_1 = new TH1D("hrapAccCrossPt1_1",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossPt1_2 = new TH1D("hrapAccCrossPt1_2",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossPt2_1 = new TH1D("hrapAccCrossPt2_1",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);
  TH1D *hrapAccCrossPt2_2 = new TH1D("hrapAccCrossPt2_2",";y_{CM};Acceptance",nYBinsCrPA,yBinCrPA);

  TH1D *hraptest1 = new TH1D("hraptest1",";|Y_{CM}|;Acceptance",nYBins23,yBin23);
  TH1D *hraptest2 = new TH1D("hraptest2",";|Y_{CM}|;Acceptance",nYBins,yBin);

  hrapAccPPNoW->Sumw2();
  hrapAccPPNoW1->Sumw2();
  hrapAccPPNoW2->Sumw2();
  hrapAccPPNoW2Cp->Sumw2();
  hrapAccPANoW->Sumw2();
  hrapAccPA2binNoW->Sumw2();
  hrapAccCrossNoW->Sumw2();
  hrapAccPANoW1->Sumw2();
  hrapAccPANoW2->Sumw2();
  hrapAccPANoW2Cp->Sumw2();
  hrapAccCrossNoW1->Sumw2();
  hrapAccCrossNoW2->Sumw2();
  hrapAccPPPt1->Sumw2();
  hrapAccPPPt2->Sumw2();
  hrapAccPPPt1_1->Sumw2();
  hrapAccPPPt1_2->Sumw2();
  hrapAccPPPt2_1->Sumw2();
  hrapAccPPPt2_2->Sumw2();
  hrapAccPP_1->Sumw2();
  hrapAccPP_2->Sumw2();
  hrapAccPPAN_1->Sumw2();
  hrapAccPPAN_2->Sumw2();
  hrapAccPA2bin_1->Sumw2();
  hrapAccPA2bin_2->Sumw2();
  hrapAccPA_1->Sumw2();
  hrapAccPA_2->Sumw2();
  hrapAccPAPt1->Sumw2();
  hrapAccPAPt2->Sumw2();
  hrapAccPAPt1_1->Sumw2();
  hrapAccPAPt1_2->Sumw2();
  hrapAccPAPt2_1->Sumw2();
  hrapAccPAPt2_2->Sumw2();

  hptAccCross_1->Sumw2();
  hptAccCross_2->Sumw2();

  hrapAccCross_1->Sumw2();
  hrapAccCross_2->Sumw2();
  hrapAccCrossNoW_1->Sumw2();
  hrapAccCrossNoW_2->Sumw2();
  hrapAccCrossPt1_1->Sumw2();
  hrapAccCrossPt1_2->Sumw2();
  hrapAccCrossPt2_1->Sumw2();
  hrapAccCrossPt2_2->Sumw2();
  hraptest1->Sumw2();
  hraptest2->Sumw2();

  TH1D *hIntAccPP_Ap1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Bp1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Cp1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Dp1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Am1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Bm1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Cm1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Dm1 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Ap2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Bp2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Cp2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Dp2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Am2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Bm2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Cm2 = (TH1D*)hIntAccPPNoW->Clone();
  TH1D *hIntAccPP_Dm2 = (TH1D*)hIntAccPPNoW->Clone();

  TH1D *hrapAccPP_Ap1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Bp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Cp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Dp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Am1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Bm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Cm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Dm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Ap2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Bp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Cp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Dp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Am2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Bm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Cm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPP_Dm2 = (TH1D*)hrapAccPPNoW->Clone();

  TH1D *hrapAccPPPt1_Ap1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Bp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Cp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Dp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Am1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Bm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Cm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Dm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Ap2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Bp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Cp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Dp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Am2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Bm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Cm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt1_Dm2 = (TH1D*)hrapAccPPNoW->Clone();

  TH1D *hrapAccPPPt2_Ap1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Bp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Cp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Dp1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Am1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Bm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Cm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Dm1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Ap2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Bp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Cp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Dp2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Am2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Bm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Cm2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *hrapAccPPPt2_Dm2 = (TH1D*)hrapAccPPNoW->Clone();

  TH1D *hIntAccPA_Ap1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Bp1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Cp1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Dp1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Am1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Bm1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Cm1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Dm1 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Ap2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Bp2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Cp2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Dp2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Am2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Bm2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Cm2 = (TH1D*)hIntAccPANoW->Clone();
  TH1D *hIntAccPA_Dm2 = (TH1D*)hIntAccPANoW->Clone();

  TH1D *hrapAccPA2bin_Ap1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Bp1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Cp1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Dp1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Am1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Bm1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Cm1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Dm1 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Ap2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Bp2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Cp2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Dp2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Am2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Bm2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Cm2 = (TH1D*)hrapAccPA2binNoW->Clone();
  TH1D *hrapAccPA2bin_Dm2 = (TH1D*)hrapAccPA2binNoW->Clone();

  TH1D *hrapAccPAPt1_Ap1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Bp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Cp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Dp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Am1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Bm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Cm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Dm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Ap2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Bp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Cp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Dp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Am2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Bm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Cm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt1_Dm2 = (TH1D*)hrapAccPANoW->Clone();

  TH1D *hrapAccPAPt2_Ap1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Bp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Cp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Dp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Am1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Bm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Cm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Dm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Ap2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Bp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Cp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Dp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Am2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Bm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Cm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPAPt2_Dm2 = (TH1D*)hrapAccPANoW->Clone();

  TH1D *hrapAccPA_Ap1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Bp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Cp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Dp1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Am1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Bm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Cm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Dm1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Ap2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Bp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Cp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Dp2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Am2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Bm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Cm2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *hrapAccPA_Dm2 = (TH1D*)hrapAccPANoW->Clone();

  TH1D *hptAccCross_Ap1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Bp1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Cp1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Dp1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Am1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Bm1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Cm1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Dm1 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Ap2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Bp2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Cp2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Dp2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Am2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Bm2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Cm2 = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *hptAccCross_Dm2 = (TH1D*)hptAccCrossNoW->Clone();

  TH1D *hrapAccCross_Ap1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Bp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Cp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Dp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Am1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Bm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Cm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Dm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Ap2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Bp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Cp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Dp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Am2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Bm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Cm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCross_Dm2 = (TH1D*)hrapAccCrossNoW->Clone();

  TH1D *hrapAccCrossPt1_Ap1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Bp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Cp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Dp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Am1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Bm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Cm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Dm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Ap2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Bp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Cp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Dp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Am2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Bm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Cm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt1_Dm2 = (TH1D*)hrapAccCrossNoW->Clone();

  TH1D *hrapAccCrossPt2_Ap1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Bp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Cp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Dp1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Am1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Bm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Cm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Dm1 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Ap2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Bp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Cp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Dp2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Am2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Bm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Cm2 = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *hrapAccCrossPt2_Dm2 = (TH1D*)hrapAccCrossNoW->Clone();

  TH1D *max_ptPPRap1 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *min_ptPPRap1 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *max_ptPPRap2 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *min_ptPPRap2 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *max_ptPARap1 = (TH1D*)hptAccPANoW->Clone();
  TH1D *min_ptPARap1 = (TH1D*)hptAccPANoW->Clone();
  TH1D *max_ptPARap2 = (TH1D*)hptAccPANoW->Clone();
  TH1D *min_ptPARap2 = (TH1D*)hptAccPANoW->Clone();
  TH1D *max_ptPP = (TH1D*)hptAccPPNoW->Clone();
  TH1D *min_ptPP = (TH1D*)hptAccPPNoW->Clone();
  TH1D *max_ptPA = (TH1D*)hptAccPANoW->Clone();
  TH1D *min_ptPA = (TH1D*)hptAccPANoW->Clone();

  TH1D *max_rapPPPt1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *min_rapPPPt1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *max_rapPPPt2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *min_rapPPPt2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *max_rapPP = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *min_rapPP = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *max_rapPA = (TH1D*)hrapAccPANoW->Clone();
  TH1D *min_rapPA = (TH1D*)hrapAccPANoW->Clone();
  TH1D *max_rapPAPt1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *min_rapPAPt1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *max_rapPAPt2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *min_rapPAPt2 = (TH1D*)hrapAccPANoW->Clone();

  TH1D *max_rapCross = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *min_rapCross = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *max_ptCross = (TH1D*)hptAccCrossNoW->Clone();
  TH1D *min_ptCross = (TH1D*)hptAccCrossNoW->Clone();

  TH1D *MD_ptPP = (TH1D*)hptAccPPNoW->Clone();
  TH1D *MD_ptPA = (TH1D*)hptAccPANoW->Clone();
  TH1D *MD_rapPP = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *MD_rapPA = (TH1D*)hrapAccPANoW->Clone();
  TH1D *MD_ptPPRap1 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *MD_ptPPRap2 = (TH1D*)hptAccPPNoW->Clone();
  TH1D *MD_ptPARap1 = (TH1D*)hptAccPANoW->Clone();
  TH1D *MD_ptPARap2 = (TH1D*)hptAccPANoW->Clone();
  TH1D *MD_rapPPPt1 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *MD_rapPPPt2 = (TH1D*)hrapAccPPNoW->Clone();
  TH1D *MD_rapPAPt1 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *MD_rapPAPt2 = (TH1D*)hrapAccPANoW->Clone();
  TH1D *MD_rapCross = (TH1D*)hrapAccCrossNoW->Clone();
  TH1D *MD_ptCross = (TH1D*)hptAccCrossNoW->Clone();

  //getting trees from skimed files
  float mass = 0.0, pt = 0.0, phi = 0.0, y = 0.0, eta = 0.0;
  float pt1 = 0.0, phi1 = 0.0, eta1 = 0.0;
  float pt2 = 0.0, phi2 = 0.0, eta2 = 0.0;

  TBranch *b_mass;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_phi;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_phi1;
  TBranch *b_phi2;

  TTree *mm = (TTree*)in->Get("mmGen");
  mm->SetBranchAddress("mass", &mass, &b_mass);
  mm->SetBranchAddress("pt", &pt, &b_pt);
  mm->SetBranchAddress("y", &y, &b_y);
  mm->SetBranchAddress("phi", &phi, &b_phi);
  mm->SetBranchAddress("eta", &eta, &b_eta);
  mm->SetBranchAddress("pt1", &pt1, &b_pt1);
  mm->SetBranchAddress("pt2", &pt2, &b_pt2);
  mm->SetBranchAddress("eta1", &eta1, &b_eta1);
  mm->SetBranchAddress("eta2", &eta2, &b_eta2);
  mm->SetBranchAddress("phi1", &phi1, &b_phi1);
  mm->SetBranchAddress("phi2", &phi2, &b_phi2);

  int nEntries = mm->GetEntries();
  //nEntries = 100;

  int cnt1 = 0; // just for counting for denominator for pp
  int cnt2 = 0; // just for counting for numerator for pp
  int cnt3 = 0; // just for counting for denominator for PA
  int cnt4 = 0; // just for counting for numerator for PA
  int cnt5 = 0; // just for counting for denominator for Cr PA
  int cnt6 = 0; // just for counting for numerator for Cr PA

  float boost = 0.47;

  for(int i = 0; i < nEntries; i++){
    mm->GetEntry(i);
    if (i%300000==0) cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    //cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    if( y>-2.4 && y<2.4 && pt<30.0 ){
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
      hrapAccPPAN_1->Fill( fabs(y), ppwgt);
      hIntAccPP2D_1->Fill ( y, pt, ppwgt );
      hIntAccPP_1->Fill( y, ppwgt );
      hIntAccPP_Ap1->Fill( y, ppwgtAp );
      hIntAccPP_Bp1->Fill( y, ppwgtBp );
      hIntAccPP_Am1->Fill( y, ppwgtAm );
      hIntAccPP_Bm1->Fill( y, ppwgtBm );
      hrapAccPP_1->Fill( y, ppwgt);
      //denominator for systematics with variations
      hrapAccPP_Ap1->Fill( y, ppwgtAp );
      hrapAccPP_Bp1->Fill( y, ppwgtBp );
      hrapAccPP_Am1->Fill( y, ppwgtAm );
      hrapAccPP_Bm1->Fill( y, ppwgtBm );
      if(Cat_==0){
        hIntAccPP_Cp1->Fill( y, ppwgtCp );
        hIntAccPP_Dp1->Fill( y, ppwgtDp );
        hIntAccPP_Cm1->Fill( y, ppwgtCm );
        hIntAccPP_Dm1->Fill( y, ppwgtDm );
        hrapAccPP_Cp1->Fill( y, ppwgtCp );
        hrapAccPP_Dp1->Fill( y, ppwgtDp );
        hrapAccPP_Cm1->Fill( y, ppwgtCm );
        hrapAccPP_Dm1->Fill( y, ppwgtDm );
      }
      if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
        hrapAccPPAN_2->Fill( fabs(y),ppwgt );
        hIntAccPP2D_2->Fill( y, pt, ppwgt );
        hIntAccPP_2->Fill( y, ppwgt );
        hIntAccPP_Ap2->Fill( y, ppwgtAp );
        hIntAccPP_Bp2->Fill( y, ppwgtBp );
        hIntAccPP_Am2->Fill( y, ppwgtAm );
        hIntAccPP_Bm2->Fill( y, ppwgtBm );
        hrapAccPP_2->Fill( y, ppwgt);
        //numerator for systematics with variations 
        hrapAccPP_Ap2->Fill( y, ppwgtAp );
        hrapAccPP_Bp2->Fill( y, ppwgtBp );
        hrapAccPP_Am2->Fill( y, ppwgtAm );
        hrapAccPP_Bm2->Fill( y, ppwgtBm );
        if(Cat_==0){
          hIntAccPP_Cp2->Fill( y, ppwgtCp );
          hIntAccPP_Dp2->Fill( y, ppwgtDp );
          hIntAccPP_Cm2->Fill( y, ppwgtCm );
          hIntAccPP_Dm2->Fill( y, ppwgtDm );
          hrapAccPP_Cp2->Fill( y, ppwgtCp );
          hrapAccPP_Dp2->Fill( y, ppwgtDp );
          hrapAccPP_Cm2->Fill( y, ppwgtCm );
          hrapAccPP_Dm2->Fill( y, ppwgtDm );
        }
      }
      if(pt<6){
        hrapAccPPPt1_1->Fill( y, ppwgt);
        hIntAccPPPt1_1->Fill( y,ppwgt);
        //denominator for systematics with variations
        hrapAccPPPt1_Ap1->Fill( y, ppwgtAp );
        hrapAccPPPt1_Bp1->Fill( y, ppwgtBp );
        hrapAccPPPt1_Am1->Fill( y, ppwgtAm );
        hrapAccPPPt1_Bm1->Fill( y, ppwgtBm );
        if(Cat_==0){
          hrapAccPPPt1_Cp1->Fill( y, ppwgtCp );
          hrapAccPPPt1_Dp1->Fill( y, ppwgtDp );
          hrapAccPPPt1_Cm1->Fill( y, ppwgtCm );
          hrapAccPPPt1_Dm1->Fill( y, ppwgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hrapAccPPPt1_2->Fill( y,ppwgt);
          hIntAccPPPt1_2->Fill( y,ppwgt);
          //numerator for systematics with variations 
          hrapAccPPPt1_Ap2->Fill( y, ppwgtAp );
          hrapAccPPPt1_Bp2->Fill( y, ppwgtBp );
          hrapAccPPPt1_Am2->Fill( y, ppwgtAm );
          hrapAccPPPt1_Bm2->Fill( y, ppwgtBm );
          if(Cat_==0){
            hrapAccPPPt1_Cp2->Fill( y, ppwgtCp );
            hrapAccPPPt1_Dp2->Fill( y, ppwgtDp );
            hrapAccPPPt1_Cm2->Fill( y, ppwgtCm );
            hrapAccPPPt1_Dm2->Fill( y, ppwgtDm );
          }
        }
      }
      if(pt>6 && pt<30){
        hrapAccPPPt2_1->Fill( y, ppwgt);
        hIntAccPPPt2_1->Fill( y,ppwgt);
        //denominator for systematics with variations
        hrapAccPPPt2_Ap1->Fill( y, ppwgtAp );
        hrapAccPPPt2_Bp1->Fill( y, ppwgtBp );
        hrapAccPPPt2_Am1->Fill( y, ppwgtAm );
        hrapAccPPPt2_Bm1->Fill( y, ppwgtBm );
        if(Cat_==0){
          hrapAccPPPt2_Cp1->Fill( y, ppwgtCp );
          hrapAccPPPt2_Dp1->Fill( y, ppwgtDp );
          hrapAccPPPt2_Cm1->Fill( y, ppwgtCm );
          hrapAccPPPt2_Dm1->Fill( y, ppwgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hIntAccPPPt2_2->Fill( y,ppwgt);
          hrapAccPPPt2_2->Fill( y,ppwgt);
          //numerator for systematics with variations 
          hrapAccPPPt2_Ap2->Fill( y, ppwgtAp );
          hrapAccPPPt2_Bp2->Fill( y, ppwgtBp );
          hrapAccPPPt2_Am2->Fill( y, ppwgtAm );
          hrapAccPPPt2_Bm2->Fill( y, ppwgtBm );
          if(Cat_==0){
            hrapAccPPPt2_Cp2->Fill( y, ppwgtCp );
            hrapAccPPPt2_Dp2->Fill( y, ppwgtDp );
            hrapAccPPPt2_Cm2->Fill( y, ppwgtCm );
            hrapAccPPPt2_Dm2->Fill( y, ppwgtDm );
          }
        }
      }

      if(-1.93<y && y<1.93){
        hptAccPP1->Fill( pt, ppwgt );
        //denominator for systematics with variations
        hptAccPP_Ap1->Fill( pt, ppwgtAp );
        hptAccPP_Bp1->Fill( pt, ppwgtBp );
        hptAccPP_Am1->Fill( pt, ppwgtAm );
        hptAccPP_Bm1->Fill( pt, ppwgtBm );

        if(Cat_==0){
          hptAccPP_Cp1->Fill( pt, ppwgtCp );
          hptAccPP_Dp1->Fill( pt, ppwgtDp );
          hptAccPP_Cm1->Fill( pt, ppwgtCm );
          hptAccPP_Dm1->Fill( pt, ppwgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPP2->Fill( pt, ppwgt );
          hptAccPP_Ap2->Fill( pt, ppwgtAp );
          hptAccPP_Bp2->Fill( pt, ppwgtBp );
          hptAccPP_Am2->Fill( pt, ppwgtAm );
          hptAccPP_Bm2->Fill( pt, ppwgtBm );
          if(Cat_==0){
            hptAccPP_Cp2->Fill( pt, ppwgtCp );
            hptAccPP_Dp2->Fill( pt, ppwgtDp );
            hptAccPP_Cm2->Fill( pt, ppwgtCm );
            hptAccPP_Dm2->Fill( pt, ppwgtDm );
          }
        }     
      }

      if(-1.93<y && y<0){
        hptAccPPRap1_1->Fill( pt, ppwgt);
        //denominator for systematics with variations
        hptAccPPRap1_Ap1->Fill( pt, ppwgtAp );
        hptAccPPRap1_Bp1->Fill( pt, ppwgtBp );
        hptAccPPRap1_Am1->Fill( pt, ppwgtAm );
        hptAccPPRap1_Bm1->Fill( pt, ppwgtBm );
        if(Cat_==0){
          hptAccPPRap1_Cp1->Fill( pt, ppwgtCp );
          hptAccPPRap1_Dp1->Fill( pt, ppwgtDp );
          hptAccPPRap1_Cm1->Fill( pt, ppwgtCm );
          hptAccPPRap1_Dm1->Fill( pt, ppwgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPPRap1_2->Fill(pt,ppwgt);
          //numerator for systematics with variations 
          //
          hptAccPPRap1_Ap2->Fill( pt, ppwgtAp );
          hptAccPPRap1_Bp2->Fill( pt, ppwgtBp );
          hptAccPPRap1_Am2->Fill( pt, ppwgtAm );
          hptAccPPRap1_Bm2->Fill( pt, ppwgtBm );
          if(Cat_==0){
            hptAccPPRap1_Cp2->Fill( pt, ppwgtCp );
            hptAccPPRap1_Dp2->Fill( pt, ppwgtDp );
            hptAccPPRap1_Cm2->Fill( pt, ppwgtCm );
            hptAccPPRap1_Dm2->Fill( pt, ppwgtDm );
          }
        }
      }
      if(0<y && y<1.93){
        hptAccPPRap2_1->Fill( pt, ppwgt);
        //denominator for systematics with variations
        hptAccPPRap2_Ap1->Fill( pt, ppwgtAp );
        hptAccPPRap2_Bp1->Fill( pt, ppwgtBp );
        hptAccPPRap2_Am1->Fill( pt, ppwgtAm );
        hptAccPPRap2_Bm1->Fill( pt, ppwgtBm );
        if(Cat_==0){
          hptAccPPRap2_Cp1->Fill( pt, ppwgtCp );
          hptAccPPRap2_Dp1->Fill( pt, ppwgtDp );
          hptAccPPRap2_Cm1->Fill( pt, ppwgtCm );
          hptAccPPRap2_Dm1->Fill( pt, ppwgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPPRap2_2->Fill(pt,ppwgt);
          //numerator for systematics with variations 
          //
          hptAccPPRap2_Ap2->Fill( pt, ppwgtAp );
          hptAccPPRap2_Bp2->Fill( pt, ppwgtBp );
          hptAccPPRap2_Am2->Fill( pt, ppwgtAm );
          hptAccPPRap2_Bm2->Fill( pt, ppwgtBm );
          if(Cat_==0){
            hptAccPPRap2_Cp2->Fill( pt, ppwgtCp );
            hptAccPPRap2_Dp2->Fill( pt, ppwgtDp );
            hptAccPPRap2_Cm2->Fill( pt, ppwgtCm );
            hptAccPPRap2_Dm2->Fill( pt, ppwgtDm );
          }
        }
      }
    //End of PP loop  
    }
    //Start pPb loop
    if( y>-2.4 && y<2.4 ){
      double pawgt = fWgtPA->Eval(pt); // apply weighting factor from functions above-defined for pa
      double pawgtAp = fWgtPA_Ap->Eval(pt);
      double pawgtBp = fWgtPA_Bp->Eval(pt);
      double pawgtAm = fWgtPA_Am->Eval(pt);
      double pawgtBm = fWgtPA_Bm->Eval(pt);
      double pawgtCp;
      double pawgtDp;
      double pawgtCm;
      double pawgtDm;
      double pacrosswgt = fWgtPAcross->Eval(pt); // apply weighting factor from functions above-defined for pa
      double pacrosswgtAp = fWgtPAcross_Ap->Eval(pt);
      double pacrosswgtBp = fWgtPAcross_Bp->Eval(pt);
      double pacrosswgtAm = fWgtPAcross_Am->Eval(pt);
      double pacrosswgtBm = fWgtPAcross_Bm->Eval(pt);
      double pacrosswgtCp;
      double pacrosswgtDp;
      double pacrosswgtCm;
      double pacrosswgtDm;
      if(Cat_==0){
        pawgtCp = fWgtPA_Cp->Eval(pt);
        pawgtDp = fWgtPA_Dp->Eval(pt);
        pawgtCm = fWgtPA_Cm->Eval(pt);
        pawgtDm = fWgtPA_Dm->Eval(pt);
        pacrosswgtCp = fWgtPAcross_Cp->Eval(pt);
        pacrosswgtDp = fWgtPAcross_Dp->Eval(pt);
        pacrosswgtCm = fWgtPAcross_Cm->Eval(pt);
        pacrosswgtDm = fWgtPAcross_Dm->Eval(pt);
      }
      hIntAccPA2D_1->Fill ( y-boost, pt, pawgt);
      hIntAccPA_1->Fill ( y-boost, pawgt);
      hIntAccPANoW_1->Fill ( y-boost );
      hIntAccCross_1->Fill ( y-boost, pacrosswgt);
      hIntAccPA_Ap1->Fill( y-boost, pawgtAp );
      hIntAccPA_Bp1->Fill( y-boost, pawgtBp );
      hIntAccPA_Am1->Fill( y-boost, pawgtAm );
      hIntAccPA_Bm1->Fill( y-boost, pawgtBm );
      hptAccCross_1->Fill( pt, pacrosswgt );
      hptAccCross_Ap1->Fill( pt, pacrosswgtAp );
      hptAccCross_Bp1->Fill( pt, pacrosswgtBp );
      hptAccCross_Am1->Fill( pt, pacrosswgtAm );
      hptAccCross_Bm1->Fill( pt, pacrosswgtBm );
      hrapAccCross_1->Fill( y-boost, pawgt );
      if(Cat_==0 && -1.465<y&&y<2.395){
      hrapAccCrossNoW_1->Fill( y-boost );}
      hrapAccCross_Ap1->Fill( y-boost, pacrosswgtAp );
      hrapAccCross_Bp1->Fill( y-boost, pacrosswgtBp );
      hrapAccCross_Am1->Fill( y-boost, pacrosswgtAm );
      hrapAccCross_Bm1->Fill( y-boost, pacrosswgtBm );
      hrapAccPA2bin_1->Fill( y-boost, pawgt );
      hrapAccPA2bin_Ap1->Fill( y-boost, pawgtAp );
      hrapAccPA2bin_Bp1->Fill( y-boost, pawgtBp );
      hrapAccPA2bin_Am1->Fill( y-boost, pawgtAm );
      hrapAccPA2bin_Bm1->Fill( y-boost, pawgtBm );
      hrapAccPA_1->Fill( y-boost, pawgt );
      hrapAccPA_Ap1->Fill( y-boost, pawgtAp );
      hrapAccPA_Bp1->Fill( y-boost, pawgtBp );
      hrapAccPA_Am1->Fill( y-boost, pawgtAm );
      hrapAccPA_Bm1->Fill( y-boost, pawgtBm );
      if(Cat_==0){
        hIntAccPA_Cp1->Fill( y-boost, pawgtCp );
        hIntAccPA_Dp1->Fill( y-boost, pawgtDp );
        hIntAccPA_Cm1->Fill( y-boost, pawgtCm );
        hIntAccPA_Dm1->Fill( y-boost, pawgtDm );
        hrapAccCross_Cp1->Fill( y-boost, pacrosswgtCp );
        hrapAccCross_Dp1->Fill( y-boost, pacrosswgtDp );
        hrapAccCross_Cm1->Fill( y-boost, pacrosswgtCm );
        hrapAccCross_Dm1->Fill( y-boost, pacrosswgtDm );
        hptAccCross_Cp1->Fill( pt, pacrosswgtCp );
        hptAccCross_Dp1->Fill( pt, pacrosswgtDp );
        hptAccCross_Cm1->Fill( pt, pacrosswgtCm );
        hptAccCross_Dm1->Fill( pt, pacrosswgtDm );
        hrapAccPA2bin_Cp1->Fill( y-boost, pawgtCp );
        hrapAccPA2bin_Dp1->Fill( y-boost, pawgtDp );
        hrapAccPA2bin_Cm1->Fill( y-boost, pawgtCm );
        hrapAccPA2bin_Dm1->Fill( y-boost, pawgtDm );
        hrapAccPA_Cp1->Fill( y-boost, pawgtCp );
        hrapAccPA_Dp1->Fill( y-boost, pawgtDp );
        hrapAccPA_Cm1->Fill( y-boost, pawgtCm );
        hrapAccPA_Dm1->Fill( y-boost, pawgtDm );
      }
      if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
        hIntAccPA2D_2->Fill ( y-boost, pt, pawgt);
        hIntAccPA_2->Fill ( y-boost, pawgt);
        if(Cat_==0 && -1.465<y&&y<2.395 && -1.465<eta1 &&eta1<2.395 && -1.465<eta2 &&eta2<2.395){
        hIntAccPANoW_2->Fill ( y-boost );}
        hIntAccCross_2->Fill ( y-boost, pacrosswgt);
        hIntAccPA_Ap2->Fill( y-boost, pawgtAp );
        hIntAccPA_Bp2->Fill( y-boost, pawgtBp );
        hIntAccPA_Am2->Fill( y-boost, pawgtAm );
        hIntAccPA_Bm2->Fill( y-boost, pawgtBm );
        hptAccCross_2->Fill( pt, pacrosswgt );
        hptAccCross_Ap2->Fill( pt, pacrosswgtAp );
        hptAccCross_Bp2->Fill( pt, pacrosswgtBp );
        hptAccCross_Am2->Fill( pt, pacrosswgtAm );
        hptAccCross_Bm2->Fill( pt, pacrosswgtBm );
        hrapAccCross_2->Fill( y-boost, pacrosswgt );
        if(Cat_==0 && 9<mass&&mass<10){
        hrapAccCrossNoW_2->Fill( y-boost );}
        hrapAccCross_Ap2->Fill( y-boost, pacrosswgtAp );
        hrapAccCross_Bp2->Fill( y-boost, pacrosswgtBp );
        hrapAccCross_Am2->Fill( y-boost, pacrosswgtAm );
        hrapAccCross_Bm2->Fill( y-boost, pacrosswgtBm );
        hrapAccPA2bin_2->Fill( y-boost, pawgt );
        hrapAccPA2bin_Ap2->Fill( y-boost, pawgtAp );
        hrapAccPA2bin_Bp2->Fill( y-boost, pawgtBp );
        hrapAccPA2bin_Am2->Fill( y-boost, pawgtAm );
        hrapAccPA2bin_Bm2->Fill( y-boost, pawgtAm );
        hrapAccPA_2->Fill( y-boost, pawgt );
        hrapAccPA_Ap2->Fill( y-boost, pawgtAp );
        hrapAccPA_Bp2->Fill( y-boost, pawgtBp );
        hrapAccPA_Am2->Fill( y-boost, pawgtAm );
        hrapAccPA_Bm2->Fill( y-boost, pawgtAm );
        if(Cat_==0){
          hIntAccPA_Cp2->Fill( y-boost, pawgtCp );
          hIntAccPA_Dp2->Fill( y-boost, pawgtDp );
          hIntAccPA_Cm2->Fill( y-boost, pawgtCm );
          hIntAccPA_Dm2->Fill( y-boost, pawgtDm );
          hptAccCross_Cp2->Fill( pt, pacrosswgtCp );
          hptAccCross_Dp2->Fill( pt, pacrosswgtDp );
          hptAccCross_Cm2->Fill( pt, pacrosswgtCm );
          hptAccCross_Dm2->Fill( pt, pacrosswgtDm );
          hrapAccCross_Cp2->Fill( y-boost, pacrosswgtCp );
          hrapAccCross_Dp2->Fill( y-boost, pacrosswgtDp );
          hrapAccCross_Cm2->Fill( y-boost, pacrosswgtCm );
          hrapAccCross_Dm2->Fill( y-boost, pacrosswgtDm );
          hrapAccPA2bin_Cp2->Fill( y-boost, pawgtCp );
          hrapAccPA2bin_Dp2->Fill( y-boost, pawgtDp );
          hrapAccPA2bin_Cm2->Fill( y-boost, pawgtCm );
          hrapAccPA2bin_Dm2->Fill( y-boost, pawgtDm );
          hrapAccPA_Cp2->Fill( y-boost, pawgtCp );
          hrapAccPA_Dp2->Fill( y-boost, pawgtDp );
          hrapAccPA_Cm2->Fill( y-boost, pawgtCm );
          hrapAccPA_Dm2->Fill( y-boost, pawgtDm );
        }
      }

      if(pt<6){
        hrapAccCrossPt1_1->Fill ( y-boost, pacrosswgt);
        //denominator for systematics with variations
        hrapAccCrossPt1_Ap1->Fill( y-boost, pacrosswgtAp );
        hrapAccCrossPt1_Bp1->Fill( y-boost, pacrosswgtBp );
        hrapAccCrossPt1_Am1->Fill( y-boost, pacrosswgtAm );
        hrapAccCrossPt1_Bm1->Fill( y-boost, pacrosswgtBm );
        if(Cat_==0){
          hrapAccCrossPt1_Cp1->Fill( y-boost, pacrosswgtCp );
          hrapAccCrossPt1_Dp1->Fill( y-boost, pacrosswgtDp );
          hrapAccCrossPt1_Cm1->Fill( y-boost, pacrosswgtCm );
          hrapAccCrossPt1_Dm1->Fill( y-boost, pacrosswgtDm );
        }
        if( -1.93+boost<y ){
          hIntAccPAPt1_1->Fill( y-boost,pawgt);
          hrapAccPAPt1_1->Fill( y-boost, pawgt);
          hrapAccPAPt1_Ap1->Fill( y-boost, pawgtAp );
          hrapAccPAPt1_Bp1->Fill( y-boost, pawgtBp );
          hrapAccPAPt1_Am1->Fill( y-boost, pawgtAm );
          hrapAccPAPt1_Bm1->Fill( y-boost, pawgtBm );
          if(Cat_==0){
            hrapAccPAPt1_Cp1->Fill( y-boost, pawgtCp );
            hrapAccPAPt1_Dp1->Fill( y-boost, pawgtDp );
            hrapAccPAPt1_Cm1->Fill( y-boost, pawgtCm );
            hrapAccPAPt1_Dm1->Fill( y-boost, pawgtDm );
          }
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hrapAccCrossPt1_2->Fill ( y-boost, pacrosswgt);
          //numerator for systematics with variations 
          hrapAccCrossPt1_Ap2->Fill( y-boost, pacrosswgtAp );
          hrapAccCrossPt1_Bp2->Fill( y-boost, pacrosswgtBp );
          hrapAccCrossPt1_Am2->Fill( y-boost, pacrosswgtAm );
          hrapAccCrossPt1_Bm2->Fill( y-boost, pacrosswgtBm );
          if(Cat_==0){
            hrapAccCrossPt1_Cp2->Fill( y-boost, pacrosswgtCp );
            hrapAccCrossPt1_Dp2->Fill( y-boost, pacrosswgtDp );
            hrapAccCrossPt1_Cm2->Fill( y-boost, pacrosswgtCm );
            hrapAccCrossPt1_Dm2->Fill( y-boost, pacrosswgtDm );
          }
          if( -1.93+boost<y  ){
            hIntAccPAPt1_2->Fill( y-boost,pawgt);
            hrapAccPAPt1_2->Fill( y-boost,pawgt);
            hrapAccPAPt1_Ap2->Fill( y-boost, pawgtAp );
            hrapAccPAPt1_Bp2->Fill( y-boost, pawgtBp );
            hrapAccPAPt1_Am2->Fill( y-boost, pawgtAm );
            hrapAccPAPt1_Bm2->Fill( y-boost, pawgtBm );
            if(Cat_==0){
              hrapAccPAPt1_Cp2->Fill( y-boost, pawgtCp );
              hrapAccPAPt1_Dp2->Fill( y-boost, pawgtDp );
              hrapAccPAPt1_Cm2->Fill( y-boost, pawgtCm );
              hrapAccPAPt1_Dm2->Fill( y-boost, pawgtDm );
            }
          }
        }
      }
      if(pt>6 && pt<30){
        hrapAccCrossPt2_1->Fill( y-boost, pacrosswgt);
        //denominator for systematics with variations
        hrapAccCrossPt2_Ap1->Fill( y-boost, pacrosswgtAp );
        hrapAccCrossPt2_Bp1->Fill( y-boost, pacrosswgtBp );
        hrapAccCrossPt2_Am1->Fill( y-boost, pacrosswgtAm );
        hrapAccCrossPt2_Bm1->Fill( y-boost, pacrosswgtBm );
        if(Cat_==0){
          hrapAccCrossPt2_Cp1->Fill( y-boost, pacrosswgtCp );
          hrapAccCrossPt2_Dp1->Fill( y-boost, pacrosswgtDp );
          hrapAccCrossPt2_Cm1->Fill( y-boost, pacrosswgtCm );
          hrapAccCrossPt2_Dm1->Fill( y-boost, pacrosswgtDm );
        }
        if( -1.93+boost<y ){
          hIntAccPAPt2_1->Fill( y-boost,pawgt);
          hrapAccPAPt2_1->Fill( y-boost, pawgt);
          hrapAccPAPt2_Ap1->Fill( y-boost, pawgtAp );
          hrapAccPAPt2_Bp1->Fill( y-boost, pawgtBp );
          hrapAccPAPt2_Am1->Fill( y-boost, pawgtAm );
          hrapAccPAPt2_Bm1->Fill( y-boost, pawgtBm );
          if(Cat_==0){
            hrapAccPAPt2_Cp1->Fill( y-boost, pawgtCp );
            hrapAccPAPt2_Dp1->Fill( y-boost, pawgtDp );
            hrapAccPAPt2_Cm1->Fill( y-boost, pawgtCm );
            hrapAccPAPt2_Dm1->Fill( y-boost, pawgtDm );
          }
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hrapAccCrossPt2_2->Fill( y-boost,pacrosswgt);
          //numerator for systematics with variations 
          hrapAccCrossPt2_Ap2->Fill( y-boost, pacrosswgtAp );
          hrapAccCrossPt2_Bp2->Fill( y-boost, pacrosswgtBp );
          hrapAccCrossPt2_Am2->Fill( y-boost, pacrosswgtAm );
          hrapAccCrossPt2_Bm2->Fill( y-boost, pacrosswgtBm );
          if(Cat_==0){
            hrapAccCrossPt2_Cp2->Fill( y-boost, pacrosswgtCp );
            hrapAccCrossPt2_Dp2->Fill( y-boost, pacrosswgtDp );
            hrapAccCrossPt2_Cm2->Fill( y-boost, pacrosswgtCm );
            hrapAccCrossPt2_Dm2->Fill( y-boost, pacrosswgtDm );
          }
          if(y>-1.93+boost){
            hIntAccPAPt2_2->Fill( y-boost,pawgt);
            hrapAccPAPt2_2->Fill( y-boost,pawgt);
            hrapAccPAPt2_Ap2->Fill( y-boost, pawgtAp );
            hrapAccPAPt2_Bp2->Fill( y-boost, pawgtBp );
            hrapAccPAPt2_Am2->Fill( y-boost, pawgtAm );
            hrapAccPAPt2_Bm2->Fill( y-boost, pawgtBm );
            if(Cat_==0){
              hrapAccPAPt2_Cp2->Fill( y-boost, pawgtCp );
              hrapAccPAPt2_Dp2->Fill( y-boost, pawgtDp );
              hrapAccPAPt2_Cm2->Fill( y-boost, pawgtCm );
              hrapAccPAPt2_Dm2->Fill( y-boost, pawgtDm );
            }
          }
        }
      }
      if(y>-1.93+boost){
        hptAccPA1->Fill( pt, pawgt );
        hptAccPA_Ap1->Fill( pt, pawgtAp );
        hptAccPA_Bp1->Fill( pt, pawgtBp );
        hptAccPA_Am1->Fill( pt, pawgtAm );
        hptAccPA_Bm1->Fill( pt, pawgtBm );
        if(Cat_==0){
          hptAccPA_Cp1->Fill( pt, pawgtAp );
          hptAccPA_Dp1->Fill( pt, pawgtBp );
          hptAccPA_Cm1->Fill( pt, pawgtAm );
          hptAccPA_Dm1->Fill( pt, pawgtBm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPA2->Fill( pt, pawgt );
          hptAccPA_Ap2->Fill( pt, pawgtAp );
          hptAccPA_Bp2->Fill( pt, pawgtBp );
          hptAccPA_Am2->Fill( pt, pawgtAm );
          hptAccPA_Bm2->Fill( pt, pawgtBm );
          if(Cat_==0){
            hptAccPA_Cp2->Fill( pt, pawgtAp );
            hptAccPA_Dp2->Fill( pt, pawgtBp );
            hptAccPA_Cm2->Fill( pt, pawgtAm );
            hptAccPA_Dm2->Fill( pt, pawgtBm );         
          }
        }
      }

      if( -1.93+boost<y && y<0+boost ){
        hptAccPARap1_1->Fill( pt, pawgt );
        //denominator for systematics with variations
        hptAccPARap1_Ap1->Fill( pt, pawgtAp );
        hptAccPARap1_Bp1->Fill( pt, pawgtBp );
        hptAccPARap1_Am1->Fill( pt, pawgtAm );
        hptAccPARap1_Bm1->Fill( pt, pawgtBm );
        hraptest1->Fill( y-boost, pawgt );
        if(Cat_==0){
          hptAccPARap1_Cp1->Fill( pt, pawgtCp );
          hptAccPARap1_Dp1->Fill( pt, pawgtDp );
          hptAccPARap1_Cm1->Fill( pt, pawgtCm );
          hptAccPARap1_Dm1->Fill( pt, pawgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPARap1_2->Fill(pt,pawgt);
          //numerator for systematics with variations 
          hptAccPARap1_Ap2->Fill( pt, pawgtAp );
          hptAccPARap1_Bp2->Fill( pt, pawgtBp );
          hptAccPARap1_Am2->Fill( pt, pawgtAm );
          hptAccPARap1_Bm2->Fill( pt, pawgtBm );
          hraptest2->Fill( y-boost, pawgt );
          if(Cat_==0){
            hptAccPARap1_Cp2->Fill( pt, pawgtCp );
            hptAccPARap1_Dp2->Fill( pt, pawgtDp );
            hptAccPARap1_Cm2->Fill( pt, pawgtCm );
            hptAccPARap1_Dm2->Fill( pt, pawgtDm );
          }
        }
      }
      if( 0+boost<y && y<1.93+boost ){
        hptAccPARap2_1->Fill( pt, pawgt);
        //denominator for systematics with variations
        hptAccPARap2_Ap1->Fill( pt, pawgtAp );
        hptAccPARap2_Bp1->Fill( pt, pawgtBp );
        hptAccPARap2_Am1->Fill( pt, pawgtAm );
        hptAccPARap2_Bm1->Fill( pt, pawgtBm );
        if(Cat_==0){
          hptAccPARap2_Cp1->Fill( pt, pawgtCp );
          hptAccPARap2_Dp1->Fill( pt, pawgtDp );
          hptAccPARap2_Cm1->Fill( pt, pawgtCm );
          hptAccPARap2_Dm1->Fill( pt, pawgtDm );
        }
        if(pt1>4.0 && pt2>4.0 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
          hptAccPARap2_2->Fill(pt,pawgt);
          //numerator for systematics with variations 
          //
          hptAccPARap2_Ap2->Fill( pt, pawgtAp );
          hptAccPARap2_Bp2->Fill( pt, pawgtBp );
          hptAccPARap2_Am2->Fill( pt, pawgtAm );
          hptAccPARap2_Bm2->Fill( pt, pawgtBm );
          if(Cat_==0){
            hptAccPARap2_Cp2->Fill( pt, pawgtCp );
            hptAccPARap2_Dp2->Fill( pt, pawgtDp );
            hptAccPARap2_Cm2->Fill( pt, pawgtCm );
            hptAccPARap2_Dm2->Fill( pt, pawgtDm );
          }
        }
      }
    }
  }
  cout<<"cnt1 : "<<cnt1<<", cnt2 : "<<cnt2<<", cnt3 : "<<cnt3<<", cnt4 : "<<cnt4<<", cnt5 : "<<cnt5<<", cnt6 : "<<cnt6<<endl;
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
  TH1D * hptAccPPRap1_2Cp = (TH1D*)hptAccPPRap1_2->Clone();
  TH1D * hptAccPPRap2_2Cp = (TH1D*)hptAccPPRap2_2->Clone();
  TH1D * hrapAccPPPt1_2Cp = (TH1D*)hrapAccPPPt1_2->Clone();
  TH1D * hrapAccPPPt2_2Cp = (TH1D*)hrapAccPPPt2_2->Clone();
  TH1D * hrapAccPP_2Cp = (TH1D*)hrapAccPP_2->Clone();
  TH1D * hrapAccPPAN_2Cp = (TH1D*)hrapAccPPAN_2->Clone();

  TH1D * hptAccPP_2Cp = (TH1D*)hptAccPP2->Clone();
  TH1D * hptAccPA_2Cp = (TH1D*)hptAccPA2->Clone();
  TH1D * hptAccPARap1_2Cp = (TH1D*)hptAccPARap1_2->Clone();
  TH1D * hptAccPARap2_2Cp = (TH1D*)hptAccPARap2_2->Clone();
  TH1D * hrapAccPAPt1_2Cp = (TH1D*)hrapAccPAPt1_2->Clone();
  TH1D * hrapAccPAPt2_2Cp = (TH1D*)hrapAccPAPt2_2->Clone();

  TH1D * hrapAccCross_2Cp = (TH1D*)hrapAccCross_2->Clone();
  TH1D * hrapAccCrossPt1_2Cp = (TH1D*)hrapAccCrossPt1_2->Clone();
  TH1D * hrapAccCrossPt2_2Cp = (TH1D*)hrapAccCrossPt2_2->Clone();
  TH1D * hptAccCross_2Cp = (TH1D*)hptAccCross_2->Clone();
  TH1D * hrapAccPA2bin_2Cp = (TH1D*)hrapAccPA2bin_2->Clone();
  TH1D * hrapAccPA_2Cp = (TH1D*)hrapAccPA_2->Clone();

  hptAccPPRap1_2Cp->Divide(hptAccPPRap1_1);
  hptAccPPRap2_2Cp->Divide(hptAccPPRap2_1);
  hrapAccPPPt1_2Cp->Divide(hrapAccPPPt1_1);
  hrapAccPPPt2_2Cp->Divide(hrapAccPPPt2_1);
  hrapAccPP_2Cp->Divide(hrapAccPP_1);
  hrapAccPPAN_2Cp->Divide(hrapAccPPAN_1);

  cout<<hrapAccPPPt1_2->GetBinContent(1)<<" / "<<hrapAccPPPt1_1->GetBinContent(1)<<" is "<<hrapAccPPPt1_2->GetBinContent(1)/hrapAccPPPt1_1->GetBinContent(1)<<endl;
  cout<<hrapAccPPPt1_2->GetBinContent(2)<<" / "<<hrapAccPPPt1_1->GetBinContent(2)<<" is "<<hrapAccPPPt1_2->GetBinContent(2)/hrapAccPPPt1_1->GetBinContent(2)<<endl;

  if(Cat_!=2){
    cout<<hrapAccPPPt1_2->GetBinContent(3)<<" / "<<hrapAccPPPt1_1->GetBinContent(3)<<" is "<<hrapAccPPPt1_2->GetBinContent(3)/hrapAccPPPt1_1->GetBinContent(3)<<endl;
    cout<<hrapAccPPPt1_2->GetBinContent(4)<<" / "<<hrapAccPPPt1_1->GetBinContent(4)<<" is "<<hrapAccPPPt1_2->GetBinContent(4)/hrapAccPPPt1_1->GetBinContent(4)<<endl;
  }
  if(Cat_==0){
    cout<<hrapAccPPPt1_2->GetBinContent(5)<<" / "<<hrapAccPPPt1_1->GetBinContent(5)<<" is "<<hrapAccPPPt1_2->GetBinContent(5)/hrapAccPPPt1_1->GetBinContent(5)<<endl;
    cout<<hrapAccPPPt1_2->GetBinContent(6)<<" / "<<hrapAccPPPt1_1->GetBinContent(6)<<" is "<<hrapAccPPPt1_2->GetBinContent(6)/hrapAccPPPt1_1->GetBinContent(6)<<endl;
    cout<<hrapAccPPPt1_2->GetBinContent(7)<<" / "<<hrapAccPPPt1_1->GetBinContent(7)<<" is "<<hrapAccPPPt1_2->GetBinContent(7)/hrapAccPPPt1_1->GetBinContent(7)<<endl;
    cout<<hrapAccPPPt1_2->GetBinContent(8)<<" / "<<hrapAccPPPt1_1->GetBinContent(8)<<" is "<<hrapAccPPPt1_2->GetBinContent(8)/hrapAccPPPt1_1->GetBinContent(8)<<endl;
  }

  hIntAccPP_2->Divide(hIntAccPP_1);
  hIntAccPPPt1_2->Divide(hIntAccPPPt1_1);
  hIntAccPPPt2_2->Divide(hIntAccPPPt2_1);
  hIntAccPA_2->Divide(hIntAccPA_1);
  hIntAccPAPt1_2->Divide(hIntAccPAPt1_1);
  hIntAccPAPt2_2->Divide(hIntAccPAPt2_1);
  hIntAccPANoW_2->Divide(hIntAccPANoW_1);

  hIntAccCross_2->Divide(hIntAccCross_1);
  hptAccPP_2Cp->Divide(hptAccPP1);
  hptAccPA_2Cp->Divide(hptAccPA1);
  hptAccPARap1_2Cp->Divide(hptAccPARap1_1);
  hptAccPARap2_2Cp->Divide(hptAccPARap2_1);
  hrapAccPAPt1_2Cp->Divide(hrapAccPAPt1_1);
  hrapAccPAPt2_2Cp->Divide(hrapAccPAPt2_1);
  hrapAccCrossPt1_2Cp->Divide(hrapAccCrossPt1_1);
  hrapAccCrossPt2_2Cp->Divide(hrapAccCrossPt2_1);
  hrapAccCross_2Cp->Divide(hrapAccCross_1);
  hrapAccCrossNoW_2->Divide(hrapAccCrossNoW_1);
  hraptest2->Divide(hraptest1);
  hptAccCross_2Cp->Divide(hptAccCross_1);
  hrapAccPA2bin_2Cp->Divide(hrapAccPA2bin_1);
  hrapAccPA_2Cp->Divide(hrapAccPA_1);
  //// for systematics with variation by parameters error
  hptAccPPRap1_Ap2->Divide(hptAccPPRap1_Ap1);
  hptAccPPRap1_Bp2->Divide(hptAccPPRap1_Bp1);
  hptAccPPRap1_Am2->Divide(hptAccPPRap1_Am1);
  hptAccPPRap1_Bm2->Divide(hptAccPPRap1_Bm1);
  hrapAccPPPt1_Ap2->Divide(hrapAccPPPt1_Ap1);
  hrapAccPPPt1_Bp2->Divide(hrapAccPPPt1_Bp1);
  hrapAccPPPt1_Am2->Divide(hrapAccPPPt1_Am1);
  hrapAccPPPt1_Bm2->Divide(hrapAccPPPt1_Bm1);
  hptAccPPRap2_Ap2->Divide(hptAccPPRap2_Ap1);
  hptAccPPRap2_Bp2->Divide(hptAccPPRap2_Bp1);
  hptAccPPRap2_Am2->Divide(hptAccPPRap2_Am1);
  hptAccPPRap2_Bm2->Divide(hptAccPPRap2_Bm1);
  hrapAccPPPt2_Ap2->Divide(hrapAccPPPt2_Ap1);
  hrapAccPPPt2_Bp2->Divide(hrapAccPPPt2_Bp1);
  hrapAccPPPt2_Am2->Divide(hrapAccPPPt2_Am1);
  hrapAccPPPt2_Bm2->Divide(hrapAccPPPt2_Bm1);
  hrapAccPP_Ap2->Divide(hrapAccPP_Ap1);
  hrapAccPP_Bp2->Divide(hrapAccPP_Bp1);
  hrapAccPP_Am2->Divide(hrapAccPP_Am1);
  hrapAccPP_Bm2->Divide(hrapAccPP_Bm1);
  hIntAccPP_Ap2->Divide(hIntAccPP_Ap1);
  hIntAccPP_Bp2->Divide(hIntAccPP_Bp1);
  hIntAccPP_Am2->Divide(hIntAccPP_Am1);
  hIntAccPP_Bm2->Divide(hIntAccPP_Bm1);
  if(Cat_==0){
    hIntAccPP_Cp2->Divide(hIntAccPP_Cp1);
    hIntAccPP_Dp2->Divide(hIntAccPP_Dp1);
    hIntAccPP_Cm2->Divide(hIntAccPP_Cm1);
    hIntAccPP_Dm2->Divide(hIntAccPP_Dm1);
    hptAccPPRap1_Cp2->Divide(hptAccPPRap1_Cp1);
    hptAccPPRap1_Dp2->Divide(hptAccPPRap1_Dp1);
    hptAccPPRap1_Cm2->Divide(hptAccPPRap1_Cm1);
    hptAccPPRap1_Dm2->Divide(hptAccPPRap1_Dm1);
    hrapAccPP_Cp2->Divide(hrapAccPP_Cp1);
    hrapAccPP_Dp2->Divide(hrapAccPP_Dp1);
    hrapAccPP_Cm2->Divide(hrapAccPP_Cm1);
    hrapAccPP_Dm2->Divide(hrapAccPP_Dm1);
    hrapAccPPPt1_Cp2->Divide(hrapAccPPPt1_Cp1);
    hrapAccPPPt1_Dp2->Divide(hrapAccPPPt1_Dp1);
    hrapAccPPPt1_Cm2->Divide(hrapAccPPPt1_Cm1);
    hrapAccPPPt1_Dm2->Divide(hrapAccPPPt1_Dm1);
    hptAccPPRap2_Cp2->Divide(hptAccPPRap2_Cp1);
    hptAccPPRap2_Dp2->Divide(hptAccPPRap2_Dp1);
    hptAccPPRap2_Cm2->Divide(hptAccPPRap2_Cm1);
    hptAccPPRap2_Dm2->Divide(hptAccPPRap2_Dm1);
    hrapAccPPPt2_Cp2->Divide(hrapAccPPPt2_Cp1);
    hrapAccPPPt2_Dp2->Divide(hrapAccPPPt2_Dp1);
    hrapAccPPPt2_Cm2->Divide(hrapAccPPPt2_Cm1);
    hrapAccPPPt2_Dm2->Divide(hrapAccPPPt2_Dm1);
  }
  hptAccPP_Ap2->Divide(hptAccPP_Ap1);
  hptAccPP_Bp2->Divide(hptAccPP_Bp1);
  hptAccPP_Am2->Divide(hptAccPP_Am1);
  hptAccPP_Bm2->Divide(hptAccPP_Bm1);
  hptAccPA_Ap2->Divide(hptAccPA_Ap1);
  hptAccPA_Bp2->Divide(hptAccPA_Bp1);
  hptAccPA_Am2->Divide(hptAccPA_Am1);
  hptAccPA_Bm2->Divide(hptAccPA_Bm1);
  hptAccPARap1_Ap2->Divide(hptAccPARap1_Ap1);
  hptAccPARap1_Bp2->Divide(hptAccPARap1_Bp1);
  hptAccPARap1_Am2->Divide(hptAccPARap1_Am1);
  hptAccPARap1_Bm2->Divide(hptAccPARap1_Bm1);
  hptAccPARap2_Ap2->Divide(hptAccPARap2_Ap1);
  hptAccPARap2_Bp2->Divide(hptAccPARap2_Bp1);
  hptAccPARap2_Am2->Divide(hptAccPARap2_Am1);
  hptAccPARap2_Bm2->Divide(hptAccPARap2_Bm1);
  hrapAccPAPt1_Ap2->Divide(hrapAccPAPt1_Ap1);
  hrapAccPAPt1_Bp2->Divide(hrapAccPAPt1_Bp1);
  hrapAccPAPt1_Am2->Divide(hrapAccPAPt1_Am1);
  hrapAccPAPt1_Bm2->Divide(hrapAccPAPt1_Bm1);
  hrapAccPAPt2_Ap2->Divide(hrapAccPAPt2_Ap1);
  hrapAccPAPt2_Bp2->Divide(hrapAccPAPt2_Bp1);
  hrapAccPAPt2_Am2->Divide(hrapAccPAPt2_Am1);
  hrapAccPAPt2_Bm2->Divide(hrapAccPAPt2_Bm1);
  hrapAccPA_Ap2->Divide(hrapAccPA_Ap1);
  hrapAccPA_Bp2->Divide(hrapAccPA_Bp1);
  hrapAccPA_Am2->Divide(hrapAccPA_Am1);
  hrapAccPA_Bm2->Divide(hrapAccPA_Bm1);

  hrapAccPA2bin_Ap2->Divide(hrapAccPA2bin_Ap1);
  hrapAccPA2bin_Bp2->Divide(hrapAccPA2bin_Bp1);
  hrapAccPA2bin_Am2->Divide(hrapAccPA2bin_Am1);
  hrapAccPA2bin_Bm2->Divide(hrapAccPA2bin_Bm1);

  hptAccCross_Ap2->Divide(hptAccCross_Ap1);
  hptAccCross_Bp2->Divide(hptAccCross_Bp1);
  hptAccCross_Am2->Divide(hptAccCross_Am1);
  hptAccCross_Bm2->Divide(hptAccCross_Bm1);

  hrapAccCrossPt1_Ap2->Divide(hrapAccCrossPt1_Ap1);
  hrapAccCrossPt1_Bp2->Divide(hrapAccCrossPt1_Bp1);
  hrapAccCrossPt1_Am2->Divide(hrapAccCrossPt1_Am1);
  hrapAccCrossPt1_Bm2->Divide(hrapAccCrossPt1_Bm1);
  hrapAccCrossPt2_Ap2->Divide(hrapAccCrossPt2_Ap1);
  hrapAccCrossPt2_Bp2->Divide(hrapAccCrossPt2_Bp1);
  hrapAccCrossPt2_Am2->Divide(hrapAccCrossPt2_Am1);
  hrapAccCrossPt2_Bm2->Divide(hrapAccCrossPt2_Bm1);

  hrapAccCross_Ap2->Divide(hrapAccCross_Ap1);
  hrapAccCross_Bp2->Divide(hrapAccCross_Bp1);
  hrapAccCross_Am2->Divide(hrapAccCross_Am1);
  hrapAccCross_Bm2->Divide(hrapAccCross_Bm1);

  hIntAccPA_Ap2->Divide(hIntAccPA_Ap1);
  hIntAccPA_Bp2->Divide(hIntAccPA_Bp1);
  hIntAccPA_Am2->Divide(hIntAccPA_Am1);
  hIntAccPA_Bm2->Divide(hIntAccPA_Bm1);
  if(Cat_==0){
    hIntAccPA_Cp2->Divide(hIntAccPA_Cp1);
    hIntAccPA_Dp2->Divide(hIntAccPA_Dp1);
    hIntAccPA_Cm2->Divide(hIntAccPA_Cm1);
    hIntAccPA_Dm2->Divide(hIntAccPA_Dm1);
    hptAccPP_Cp2->Divide(hptAccPP_Cp1);
    hptAccPP_Dp2->Divide(hptAccPP_Dp1);
    hptAccPP_Cm2->Divide(hptAccPP_Cm1);
    hptAccPP_Dm2->Divide(hptAccPP_Dm1);
    hptAccPA_Cp2->Divide(hptAccPA_Cp1);
    hptAccPA_Dp2->Divide(hptAccPA_Dp1);
    hptAccPA_Cm2->Divide(hptAccPA_Cm1);
    hptAccPA_Dm2->Divide(hptAccPA_Dm1);
    hptAccPARap1_Cp2->Divide(hptAccPARap1_Cp1);
    hptAccPARap1_Dp2->Divide(hptAccPARap1_Dp1);
    hptAccPARap1_Cm2->Divide(hptAccPARap1_Cm1);
    hptAccPARap1_Dm2->Divide(hptAccPARap1_Dm1);
    hptAccPARap2_Cp2->Divide(hptAccPARap2_Cp1);
    hptAccPARap2_Dp2->Divide(hptAccPARap2_Dp1);
    hptAccPARap2_Cm2->Divide(hptAccPARap2_Cm1);
    hptAccPARap2_Dm2->Divide(hptAccPARap2_Dm1);
    hrapAccPA_Cp2->Divide(hrapAccPA_Cp1);
    hrapAccPA_Dp2->Divide(hrapAccPA_Dp1);
    hrapAccPA_Cm2->Divide(hrapAccPA_Cm1);
    hrapAccPA_Dm2->Divide(hrapAccPA_Dm1);
    hrapAccPAPt1_Cp2->Divide(hrapAccPAPt1_Cp1);
    hrapAccPAPt1_Dp2->Divide(hrapAccPAPt1_Dp1);
    hrapAccPAPt1_Cm2->Divide(hrapAccPAPt1_Cm1);
    hrapAccPAPt1_Dm2->Divide(hrapAccPAPt1_Dm1);
    hrapAccPAPt2_Cp2->Divide(hrapAccPAPt2_Cp1);
    hrapAccPAPt2_Dp2->Divide(hrapAccPAPt2_Dp1);
    hrapAccPAPt2_Cm2->Divide(hrapAccPAPt2_Cm1);
    hrapAccPAPt2_Dm2->Divide(hrapAccPAPt2_Dm1);
    hptAccCross_Cp2->Divide(hptAccCross_Cp1);
    hptAccCross_Dp2->Divide(hptAccCross_Dp1);
    hptAccCross_Cm2->Divide(hptAccCross_Cm1);
    hptAccCross_Dm2->Divide(hptAccCross_Dm1);
    hrapAccPA2bin_Cp2->Divide(hrapAccPA2bin_Cp1);
    hrapAccPA2bin_Dp2->Divide(hrapAccPA2bin_Dp1);
    hrapAccPA2bin_Cm2->Divide(hrapAccPA2bin_Cm1);
    hrapAccPA2bin_Dm2->Divide(hrapAccPA2bin_Dm1);

    hrapAccCross_Cp2->Divide(hrapAccCross_Cp1);
    hrapAccCross_Dp2->Divide(hrapAccCross_Dp1);
    hrapAccCross_Cm2->Divide(hrapAccCross_Cm1);
    hrapAccCross_Dm2->Divide(hrapAccCross_Dm1);
    hrapAccCrossPt1_Cp2->Divide(hrapAccCrossPt1_Cp1);
    hrapAccCrossPt1_Dp2->Divide(hrapAccCrossPt1_Dp1);
    hrapAccCrossPt1_Cm2->Divide(hrapAccCrossPt1_Cm1);
    hrapAccCrossPt1_Dm2->Divide(hrapAccCrossPt1_Dm1);
    hrapAccCrossPt2_Cp2->Divide(hrapAccCrossPt2_Cp1);
    hrapAccCrossPt2_Dp2->Divide(hrapAccCrossPt2_Dp1);
    hrapAccCrossPt2_Cm2->Divide(hrapAccCrossPt2_Cm1);
    hrapAccCrossPt2_Dm2->Divide(hrapAccCrossPt2_Dm1);
  }
  hIntAccPP_Ap2->SetName(  "hIntAccPP_Sys_Ap"  );
  hIntAccPP_Bp2->SetName(  "hIntAccPP_Sys_Bp"  );
  hIntAccPP_Am2->SetName(  "hIntAccPP_Sys_Am"  );
  hIntAccPP_Bm2->SetName(  "hIntAccPP_Sys_Bm"  );
  hptAccPPRap1_Ap2->SetName(  "hptAccPPRap1_Sys_Ap"  );
  hptAccPPRap1_Bp2->SetName(  "hptAccPPRap1_Sys_Bp"  );
  hptAccPPRap1_Am2->SetName(  "hptAccPPRap1_Sys_Am"  );
  hptAccPPRap1_Bm2->SetName(  "hptAccPPRap1_Sys_Bm"  );
  hrapAccPPPt1_Ap2->SetName(  "hrapAccPPPt1_Sys_Ap"  );
  hrapAccPPPt1_Bp2->SetName(  "hrapAccPPPt1_Sys_Bp"  );
  hrapAccPPPt1_Am2->SetName(  "hrapAccPPPt1_Sys_Am"  );
  hrapAccPPPt1_Bm2->SetName(  "hrapAccPPPt1_Sys_Bm"  );
  hptAccPPRap2_Ap2->SetName(  "hptAccPPRap2_Sys_Ap"  );
  hptAccPPRap2_Bp2->SetName(  "hptAccPPRap2_Sys_Bp"  );
  hptAccPPRap2_Am2->SetName(  "hptAccPPRap2_Sys_Am"  );
  hptAccPPRap2_Bm2->SetName(  "hptAccPPRap2_Sys_Bm"  );
  hrapAccPPPt2_Ap2->SetName(  "hrapAccPPPt2_Sys_Ap"  );
  hrapAccPPPt2_Bp2->SetName(  "hrapAccPPPt2_Sys_Bp"  );
  hrapAccPPPt2_Am2->SetName(  "hrapAccPPPt2_Sys_Am"  );
  hrapAccPPPt2_Bm2->SetName(  "hrapAccPPPt2_Sys_Bm"  );
  if(Cat_==0){
    hptAccPPRap1_Dp2->SetName(  "hptAccPPRap1_Sys_Dp"  );
    hptAccPPRap1_Cm2->SetName(  "hptAccPPRap1_Sys_Cm"  );
    hptAccPPRap1_Dm2->SetName(  "hptAccPPRap1_Sys_Dm"  );
    hrapAccPPPt1_Cp2->SetName(  "hrapAccPPPt1_Sys_Cp"  );
    hrapAccPPPt1_Dp2->SetName(  "hrapAccPPPt1_Sys_Dp"  );
    hrapAccPPPt1_Cm2->SetName(  "hrapAccPPPt1_Sys_Cm"  );
    hrapAccPPPt1_Dm2->SetName(  "hrapAccPPPt1_Sys_Dm"  );
    hptAccPPRap2_Cp2->SetName(  "hptAccPPRap2_Sys_Cp"  );
    hptAccPPRap2_Dp2->SetName(  "hptAccPPRap2_Sys_Dp"  );
    hptAccPPRap2_Cm2->SetName(  "hptAccPPRap2_Sys_Cm"  );
    hptAccPPRap2_Dm2->SetName(  "hptAccPPRap2_Sys_Dm"  );
    hrapAccPPPt2_Cp2->SetName(  "hrapAccPPPt2_Sys_Cp"  );
    hrapAccPPPt2_Dp2->SetName(  "hrapAccPPPt2_Sys_Dp"  );
    hrapAccPPPt2_Cm2->SetName(  "hrapAccPPPt2_Sys_Cm"  );
    hrapAccPPPt2_Dm2->SetName(  "hrapAccPPPt2_Sys_Dm"  );
  }
  hIntAccPA_Ap2->SetName(  "hIntAccPA_Sys_Ap"  );
  hIntAccPA_Bp2->SetName(  "hIntAccPA_Sys_Bp"  );
  hIntAccPA_Am2->SetName(  "hIntAccPA_Sys_Am"  );
  hIntAccPA_Bm2->SetName(  "hIntAccPA_Sys_Bm"  );
  hptAccPARap1_Ap2->SetName(  "hptAccPARap1_Sys_Ap"  );
  hptAccPARap1_Bp2->SetName(  "hptAccPARap1_Sys_Bp"  );
  hptAccPARap1_Am2->SetName(  "hptAccPARap1_Sys_Am"  );
  hptAccPARap1_Bm2->SetName(  "hptAccPARap1_Sys_Bm"  );
  hptAccPARap2_Ap2->SetName(  "hptAccPARap2_Sys_Ap"  );
  hptAccPARap2_Bp2->SetName(  "hptAccPARap2_Sys_Bp"  );
  hptAccPARap2_Am2->SetName(  "hptAccPARap2_Sys_Am"  );
  hptAccPARap2_Bm2->SetName(  "hptAccPARap2_Sys_Bm"  );
  hrapAccPAPt1_Ap2->SetName(  "hrapAccPAPt1_Sys_Ap"  );
  hrapAccPAPt1_Bp2->SetName(  "hrapAccPAPt1_Sys_Bp"  );
  hrapAccPAPt1_Am2->SetName(  "hrapAccPAPt1_Sys_Am"  );
  hrapAccPAPt1_Bm2->SetName(  "hrapAccPAPt1_Sys_Bm"  );
  hrapAccPAPt2_Ap2->SetName(  "hrapAccPAPt2_Sys_Ap"  );
  hrapAccPAPt2_Bp2->SetName(  "hrapAccPAPt2_Sys_Bp"  );
  hrapAccPAPt2_Am2->SetName(  "hrapAccPAPt2_Sys_Am"  );
  hrapAccPAPt2_Bm2->SetName(  "hrapAccPAPt2_Sys_Bm"  );
  hrapAccPA_Ap2->SetName(  "hrapAccPA_Sys_Ap"  );
  hrapAccPA_Bp2->SetName(  "hrapAccPA_Sys_Bp"  );
  hrapAccPA_Am2->SetName(  "hrapAccPA_Sys_Am"  );
  hrapAccPA_Bm2->SetName(  "hrapAccPA_Sys_Bm"  );
  hrapAccCrossPt1_Ap2->SetName(  "hrapAccCrossPt1_Sys_Ap"  );
  hrapAccCrossPt1_Bp2->SetName(  "hrapAccCrossPt1_Sys_Bp"  );
  hrapAccCrossPt1_Am2->SetName(  "hrapAccCrossPt1_Sys_Am"  );
  hrapAccCrossPt1_Bm2->SetName(  "hrapAccCrossPt1_Sys_Bm"  );
  hrapAccCrossPt2_Ap2->SetName(  "hrapAccCrossPt2_Sys_Ap"  );
  hrapAccCrossPt2_Bp2->SetName(  "hrapAccCrossPt2_Sys_Bp"  );
  hrapAccCrossPt2_Am2->SetName(  "hrapAccCrossPt2_Sys_Am"  );
  hrapAccCrossPt2_Bm2->SetName(  "hrapAccCrossPt2_Sys_Bm"  );
  if(Cat_==0){
    hIntAccPA_Cp2->SetName(  "hIntAccPA_Sys_Cp"  );
    hIntAccPA_Dp2->SetName(  "hIntAccPA_Sys_Dp"  );
    hIntAccPA_Cm2->SetName(  "hIntAccPA_Sys_Cm"  );
    hIntAccPA_Cp2->SetName(  "hIntAccPA_Sys_Cp"  );
    hptAccPARap1_Cp2->SetName(  "hptAccPARap1_Sys_Cp"  );
    hptAccPARap1_Dp2->SetName(  "hptAccPARap1_Sys_Dp"  );
    hptAccPARap1_Cm2->SetName(  "hptAccPARap1_Sys_Cm"  );
    hptAccPARap1_Dm2->SetName(  "hptAccPARap1_Sys_Dm"  );
    hptAccPARap2_Cp2->SetName(  "hptAccPARap2_Sys_Cp"  );
    hptAccPARap2_Dp2->SetName(  "hptAccPARap2_Sys_Dp"  );
    hptAccPARap2_Cm2->SetName(  "hptAccPARap2_Sys_Cm"  );
    hptAccPARap2_Dm2->SetName(  "hptAccPARap2_Sys_Dm"  );
    hrapAccPAPt1_Cp2->SetName(  "hrapAccPAPt1_Sys_Cp"  );
    hrapAccPAPt1_Dp2->SetName(  "hrapAccPAPt1_Sys_Dp"  );
    hrapAccPAPt1_Cm2->SetName(  "hrapAccPAPt1_Sys_Cm"  );
    hrapAccPAPt1_Dm2->SetName(  "hrapAccPAPt1_Sys_Dm"  );
    hrapAccPAPt2_Cp2->SetName(  "hrapAccPAPt2_Sys_Cp"  );
    hrapAccPAPt2_Dp2->SetName(  "hrapAccPAPt2_Sys_Dp"  );
    hrapAccPAPt2_Cm2->SetName(  "hrapAccPAPt2_Sys_Cm"  );
    hrapAccPAPt2_Dm2->SetName(  "hrapAccPAPt2_Sys_Dm"  );
    hrapAccPA_Cp2->SetName(  "hrapAccPA_Sys_Cp"  );
    hrapAccPA_Dp2->SetName(  "hrapAccPA_Sys_Dp"  );
    hrapAccPA_Cm2->SetName(  "hrapAccPA_Sys_Cm"  );
    hrapAccPA_Dm2->SetName(  "hrapAccPA_Sys_Dm"  );
    hrapAccCrossPt1_Cp2->SetName(  "hrapAccCrossPt1_Sys_Cp"  );
    hrapAccCrossPt1_Dp2->SetName(  "hrapAccCrossPt1_Sys_Dp"  );
    hrapAccCrossPt1_Cm2->SetName(  "hrapAccCrossPt1_Sys_Cm"  );
    hrapAccCrossPt1_Dm2->SetName(  "hrapAccCrossPt1_Sys_Dm"  );
    hrapAccCrossPt2_Cp2->SetName(  "hrapAccCrossPt2_Sys_Cp"  );
    hrapAccCrossPt2_Dp2->SetName(  "hrapAccCrossPt2_Sys_Dp"  );
    hrapAccCrossPt2_Cm2->SetName(  "hrapAccCrossPt2_Sys_Cm"  );
    hrapAccCrossPt2_Dm2->SetName(  "hrapAccCrossPt2_Sys_Dm"  );
  }

  if(Cat_ == 0) cout<<" %%%%% Upsilon 1S %%%%% "<<endl;
  if(Cat_ == 1) cout<<" %%%%% Upsilon 2S %%%%% "<<endl;
  if(Cat_ == 2) cout<<" %%%%% Upsilon 3S %%%%% "<<endl;

  cout<<""<<endl;
  cout<<"%%% Pt  %%%"<<endl;
  cout<<"Pt_PP_Rap1 :: "<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    //  cout<<"No Weighted pp :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPPNoW2Cp->GetBinError(i+1))<<endl;
    //  cout<<"No Weighted PA :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPANoW2Cp->GetBinError(i+1))<<endl;
    cout<<i+1<<" :: "<<Form("%0.3f",hptAccPPRap1_2Cp->GetBinContent(i+1))<<" & "<<Form("%0.4f",hptAccPARap1_2Cp->GetBinContent(i+1))<<" & "<<hptAccPPRap1_2Cp->GetBinContent(i+1)/hptAccPARap1_2Cp->GetBinContent(i+1)<<endl;
  }
  cout<<"Pt_PP_Rap2 :: "<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    //  cout<<"No Weighted pp :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPPNoW2Cp->GetBinContent(i+1))<<endl;
    //  cout<<"No Weighted PA :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPANoW2Cp->GetBinContent(i+1))<<endl;
    cout<<i+1<<" :: "<<Form("%0.3f",hptAccPPRap2_2Cp->GetBinContent(i+1))<<" & "<<Form("%0.4f",hptAccPARap2_2Cp->GetBinContent(i+1))<<" & "<<hptAccPPRap2_2Cp->GetBinContent(i+1)/hptAccPARap2_2Cp->GetBinContent(i+1)<<endl;
  }
  cout<<"Pt PA :: -1.93 < y < 1.93 :: "<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    //  cout<<"No Weighted pp :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPPNoW2Cp->GetBinError(i+1))<<endl;
    //  cout<<"No Weighted PA :: "<<i+1<<" :: "<<Form("%0.3f",hptAccPANoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPANoW2Cp->GetBinError(i+1))<<endl;
    cout<<i+1<<" :: "<<Form("%0.3f",hptAccPP_2Cp->GetBinContent(i+1))<<" & "<<Form("%0.4f",hptAccPA_2Cp->GetBinContent(i+1))<<" & "<<hptAccPP_2Cp->GetBinContent(i+1)/hptAccPA_2Cp->GetBinContent(i+1)<<endl;
  }
  cout<<""<<endl;
  cout<<"%%% Rapidity  %%%"<<endl;
  cout<<"---===Acc in rapidity bins for pp===---"<<endl;
  cout<<"Rap_Pt1 :: "<<endl;
  cout<<hrapAccCrossPt1_2Cp->GetBinContent(1)<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    //  cout<<"No Weighted pp :: "<<i+1<<" "<<Form("%0.3f",hrapAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPPNoW2Cp->GetBinContent(i+1))<<endl;
    cout<<i+1<<" :: "<<Form("%0.3f",hrapAccPPPt1_2Cp->GetBinContent(i+1))<<" & "<<Form("%0.4f",hrapAccPAPt1_2Cp->GetBinContent(i+1))<<" & "<<hrapAccPPPt1_2Cp->GetBinContent(i+1)/hrapAccPAPt1_2Cp->GetBinContent(i+1)<<endl;
  }
  cout<<"Rap_Pt2 :: "<<endl;
  cout<<hrapAccCrossPt2_2Cp->GetBinContent(1)<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    //  cout<<"No Weighted pp :: "<<i+1<<" "<<Form("%0.3f",hrapAccPPNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPPNoW2Cp->GetBinContent(i+1))<<endl;
    cout<<i+1<<" :: "<<Form("%0.3f",hrapAccPPPt2_2Cp->GetBinContent(i+1))<<" & "<<Form("%0.4f",hrapAccPAPt2_2Cp->GetBinContent(i+1))<<" & "<<hrapAccPPPt2_2Cp->GetBinContent(i+1)/hrapAccPAPt2_2Cp->GetBinContent(i+1)<<endl;
  }

  cout<<""<<endl;
  //for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
  //  cout<<"%%% Systematics :: Rapidity %%%"<<endl;
  //  cout<<"pp "<<Form("%0.3f",fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPNoW2Cp->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100)<<endl;
  //  cout<<"PA "<<Form("%0.3f",fabs(hrapAccPA2Cp->GetBinContent(i+1)-hrapAccPANoW2Cp->GetBinContent(i+1))/hrapAccPA2Cp->GetBinContent(i+1)*100)<<endl;
  //}

  // Calculation of Systematics for variation by parameters error
  cout<<"\n %%% Calculation of Systematics for variation by error of params up and down %%% "<<endl;
  cout<<"\n %%% Systematics :: Int Bins %%%"<<endl;
  for(int i = 0; i < hIntAccPP->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Ap2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1),fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Am2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1));
    double sys2 = max(fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Bp2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1),fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Bm2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1));
    double sys5 = max(fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Ap2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1),fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Am2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1));
    double sys6 = max(fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Bp2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1),fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Bm2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Cp2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1),fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Cm2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1));
      double sys4 = max(fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Dp2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1),fabs(hIntAccPP_2->GetBinContent(i+1)-hIntAccPP_Dm2->GetBinContent(i+1))/hIntAccPP_2->GetBinContent(i+1));
      double sys7 = max(fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Cp2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1),fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Cm2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1));
      double sys8 = max(fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Dp2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1),fabs(hIntAccPA_2->GetBinContent(i+1)-hIntAccPA_Dm2->GetBinContent(i+1))/hIntAccPA_2->GetBinContent(i+1));
      double Sys_Int_PP = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_Int_PA = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"PP Int A : "<<Form("%0.3f",sys1*100)<<", PP Int B : "<<Form("%0.3f",sys2*100)<<", PP Int C : "<<Form("%0.3f",sys3*100)<<", PP Int D : "<<Form("%0.3f",sys4*100)<<", PP Int total : "<<Form("%0.3f",abs(Sys_Int_PP*100))<<endl;
      cout<<"PA Int A : "<<Form("%0.3f",sys5*100)<<", PA Int B : "<<Form("%0.3f",sys6*100)<<", PA Int C : "<<Form("%0.3f",sys7*100)<<", PA Int D : "<<Form("%0.3f",sys8*100)<<", PA Int total : "<<Form("%0.3f",abs(Sys_Int_PA*100))<<endl;
      hIntSysAccPP->SetBinContent(i+1,Sys_Int_PP);
      hIntSysAccPA->SetBinContent(i+1,Sys_Int_PA);
    }
    else if(Cat_!=0){
      double Sys_Int_PP = max(sys1,sys2);
      double Sys_Int_PA = max(sys5,sys6);
      cout<<"PP Int A : "<<Form("%0.3f",sys1*100)<<", PP Int B : "<<Form("%0.3f",sys2*100)<<", PP Int total : "<<Form("%0.3f",abs(Sys_Int_PP*100))<<endl;
      cout<<"PA Int A : "<<Form("%0.3f",sys5*100)<<", PA Int B : "<<Form("%0.3f",sys6*100)<<", PA Int total : "<<Form("%0.3f",abs(Sys_Int_PA*100))<<endl;
      hIntSysAccPP->SetBinContent(i+1,Sys_Int_PP);
      hIntSysAccPA->SetBinContent(i+1,Sys_Int_PA);
    }
  }


  cout<<"\n %%% Systematics :: Pt Bins %%%"<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Ap2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1),fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Am2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1));
    double sys2 = max(fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Bp2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1),fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Bm2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1));
    double sys5 = max(fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Ap2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1),fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Am2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Bp2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1),fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Bm2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Cp2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1),fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Cm2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1));
      double sys4 = max(fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Dp2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1),fabs(hptAccPP_2Cp->GetBinContent(i+1)-hptAccPP_Dm2->GetBinContent(i+1))/hptAccPP_2Cp->GetBinContent(i+1));
      double sys7 = max(fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Cp2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1),fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Cm2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Dp2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1),fabs(hptAccPA_2Cp->GetBinContent(i+1)-hptAccPA_Dm2->GetBinContent(i+1))/hptAccPA_2Cp->GetBinContent(i+1));
      double Sys_pt_PP = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_pt_PA = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"PP  A : "<<Form("%0.3f",sys1*100)<<", PP  B : "<<Form("%0.3f",sys2*100)<<", PP  C : "<<Form("%0.3f",sys3*100)<<", PP  D : "<<Form("%0.3f",sys4*100)<<", PP  total : "<<Form("%0.3f",abs(Sys_pt_PP*100))<<endl;
      cout<<"PA  A : "<<Form("%0.3f",sys5*100)<<", PA  B : "<<Form("%0.3f",sys6*100)<<", PA  C : "<<Form("%0.3f",sys7*100)<<", PA  D : "<<Form("%0.3f",sys8*100)<<", PA  total : "<<Form("%0.3f",abs(Sys_pt_PA*100))<<endl;
      hptSysAccPP->SetBinContent(i+1,Sys_pt_PP);
      hptSysAccPA->SetBinContent(i+1,Sys_pt_PA);
    }
    else if(Cat_!=0){
      double Sys_pt_PP = max(sys1,sys2);
      double Sys_pt_PA = max(sys5,sys6);
      cout<<"PP  A : "<<Form("%0.3f",sys1*100)<<", PP  B : "<<Form("%0.3f",sys2*100)<<", PP  total : "<<Form("%0.3f",abs(Sys_pt_PP*100))<<endl;
      cout<<"PA  A : "<<Form("%0.3f",sys5*100)<<", PA  B : "<<Form("%0.3f",sys6*100)<<", PA  total : "<<Form("%0.3f",abs(Sys_pt_PA*100))<<endl;
      hptSysAccPP->SetBinContent(i+1,Sys_pt_PP);
      hptSysAccPA->SetBinContent(i+1,Sys_pt_PA);
    }
  }
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Ap2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1),fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Am2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1));
    double sys2 = max(fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Bp2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1),fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Bm2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1));
    double sys5 = max(fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Ap2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1),fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Am2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Bp2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1),fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Bm2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Cp2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1),fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Cm2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1));
      double sys4 = max(fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Dp2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1),fabs(hptAccPPRap1_2Cp->GetBinContent(i+1)-hptAccPPRap1_Dm2->GetBinContent(i+1))/hptAccPPRap1_2Cp->GetBinContent(i+1));
      double sys7 = max(fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Cp2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1),fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Cm2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Dp2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1),fabs(hptAccPARap1_2Cp->GetBinContent(i+1)-hptAccPARap1_Dm2->GetBinContent(i+1))/hptAccPARap1_2Cp->GetBinContent(i+1));
      double Sys_pt_PPy1 = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_pt_PAy1 = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"PP y1 A : "<<Form("%0.3f",sys1*100)<<", PP y1 B : "<<Form("%0.3f",sys2*100)<<", PP y1 C : "<<Form("%0.3f",sys3*100)<<", PP y1 D : "<<Form("%0.3f",sys4*100)<<", PP y1 total : "<<Form("%0.3f",abs(Sys_pt_PPy1*100))<<endl;
      cout<<"PA y1 A : "<<Form("%0.3f",sys5*100)<<", PA y1 B : "<<Form("%0.3f",sys6*100)<<", PA y1 C : "<<Form("%0.3f",sys7*100)<<", PA y1 D : "<<Form("%0.3f",sys8*100)<<", PA y1 total : "<<Form("%0.3f",abs(Sys_pt_PAy1*100))<<endl;
      hptSysAccPPRap1->SetBinContent(i+1,Sys_pt_PPy1);
      hptSysAccPARap1->SetBinContent(i+1,Sys_pt_PAy1);
    }
    else if(Cat_!=0){
      double Sys_pt_PPy1 = max(sys1,sys2);
      double Sys_pt_PAy1 = max(sys5,sys6);
      cout<<"PP y1 A : "<<Form("%0.3f",sys1*100)<<", PP y1 B : "<<Form("%0.3f",sys2*100)<<", PP y1 total : "<<Form("%0.3f",abs(Sys_pt_PPy1*100))<<endl;
      cout<<"PA y1 A : "<<Form("%0.3f",sys5*100)<<", PA y1 B : "<<Form("%0.3f",sys6*100)<<", PA y1 total : "<<Form("%0.3f",abs(Sys_pt_PAy1*100))<<endl;
      hptSysAccPPRap1->SetBinContent(i+1,Sys_pt_PPy1);
      hptSysAccPARap1->SetBinContent(i+1,Sys_pt_PAy1);
    }
  }
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Ap2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1),fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Am2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1));
    double sys2 = max(fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Bp2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1),fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Bm2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1));
    double sys5 = max(fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Ap2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1),fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Am2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Bp2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1),fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Bm2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Cp2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1),fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Cm2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1));
      double sys4 = max(fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Dp2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1),fabs(hptAccPPRap2_2Cp->GetBinContent(i+1)-hptAccPPRap2_Dm2->GetBinContent(i+1))/hptAccPPRap2_2Cp->GetBinContent(i+1));
      double sys7 = max(fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Cp2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1),fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Cm2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Dp2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1),fabs(hptAccPARap2_2Cp->GetBinContent(i+1)-hptAccPARap2_Dm2->GetBinContent(i+1))/hptAccPARap2_2Cp->GetBinContent(i+1));
      double Sys_pt_PPy2 = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_pt_PAy2 = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"PP y2 A : "<<Form("%0.3f",sys1*100)<<", PP y2 B : "<<Form("%0.3f",sys2*100)<<", PP y2 C : "<<Form("%0.3f",sys3*100)<<", PP y2 D : "<<Form("%0.3f",sys4*100)<<", PP y2 total : "<<Form("%0.3f",abs(Sys_pt_PPy2*100))<<endl;
      cout<<"PA y2 A : "<<Form("%0.3f",sys5*100)<<", PA y2 B : "<<Form("%0.3f",sys6*100)<<", PA y2 C : "<<Form("%0.3f",sys7*100)<<", PA y2 D : "<<Form("%0.3f",sys8*100)<<", PA y2 total : "<<Form("%0.3f",abs(Sys_pt_PAy2*100))<<endl;
      hptSysAccPPRap2->SetBinContent(i+1,Sys_pt_PPy2);
      hptSysAccPARap2->SetBinContent(i+1,Sys_pt_PAy2);
    }
    else if(Cat_!=0){
      double Sys_pt_PPy2 = max(sys1,sys2);
      double Sys_pt_PAy2 = max(sys5,sys6);
      cout<<"PP y2 A : "<<Form("%0.3f",sys1*100)<<", PP y2 B : "<<Form("%0.3f",sys2*100)<<", PP y2 total : "<<Form("%0.3f",abs(Sys_pt_PPy2*100))<<endl;
      cout<<"PA y2 A : "<<Form("%0.3f",sys5*100)<<", PA y2 B : "<<Form("%0.3f",sys6*100)<<", PA y2 total : "<<Form("%0.3f",abs(Sys_pt_PAy2*100))<<endl;
      hptSysAccPPRap2->SetBinContent(i+1,Sys_pt_PPy2);
      hptSysAccPARap2->SetBinContent(i+1,Sys_pt_PAy2);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hptAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"check for pt1 A value :: total : "<<hptAccCross_2Cp->GetBinContent(i+1)<<" Ap : "<<hptAccCross_Ap2->GetBinContent(i+1)<<" Am : "<<hptAccCross_Am2->GetBinContent(i+1)<<endl;
    cout<<"check for pt1 B value :: total : "<<hptAccCross_2Cp->GetBinContent(i+1)<<" Bp : "<<hptAccCross_Bp2->GetBinContent(i+1)<<" Bm : "<<hptAccCross_Bm2->GetBinContent(i+1)<<endl;
    double sys5 = max(fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Ap2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1),fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Am2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Bp2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1),fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Bm2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys7 = max(fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Cp2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1),fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Cm2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Dp2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1),fabs(hptAccCross_2Cp->GetBinContent(i+1)-hptAccCross_Dm2->GetBinContent(i+1))/hptAccCross_2Cp->GetBinContent(i+1));
      double Sys_pt_Cross = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"Cross  A : "<<Form("%0.3f",sys5*100)<<", Cross  B : "<<Form("%0.3f",sys6*100)<<", Cross  C : "<<Form("%0.3f",sys7*100)<<", Cross  D : "<<Form("%0.3f",sys8*100)<<", Cross  total : "<<Form("%0.3f",abs(Sys_pt_Cross*100))<<endl;
      hptSysAccCross->SetBinContent(i+1,Sys_pt_Cross);
    }
    else if(Cat_!=0){
      double Sys_pt_Cross = max(sys5,sys6);
      cout<<"Cross  A : "<<Form("%0.3f",sys5*100)<<", Cross  B : "<<Form("%0.3f",sys6*100)<<", Cross  total : "<<Form("%0.3f",abs(Sys_pt_Cross*100))<<endl;
      hptSysAccCross->SetBinContent(i+1,Sys_pt_Cross);
    }
  }
  /////////////////////////////////////////////////////////////////
  cout<<"\n %%% Systematics :: Rapidity Bins %%%"<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"check for pt1 A value :: total : "<<hrapAccPPPt1_2Cp->GetBinContent(i+1)<<" Ap : "<<hrapAccPPPt1_Ap2->GetBinContent(i+1)<<" Am : "<<hrapAccPPPt1_Am2->GetBinContent(i+1)<<endl;
    cout<<"check for pt1 B value :: total : "<<hrapAccPPPt1_2Cp->GetBinContent(i+1)<<" Bp : "<<hrapAccPPPt1_Bp2->GetBinContent(i+1)<<" Bm : "<<hrapAccPPPt1_Bm2->GetBinContent(i+1)<<endl;
    double sys1 = max(fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Ap2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Am2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1));
    double sys2 = max(fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Bp2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Bm2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1));
    double sys5 = max(fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Ap2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Am2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Bp2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Bm2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1));
    double sys9 = max(fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Ap2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1),fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Am2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1));
    double sys10 = max(fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Bp2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1),fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Bm2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Cp2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Cm2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1));
      double sys4 = max(fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Dp2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt1_2Cp->GetBinContent(i+1)-hrapAccPPPt1_Dm2->GetBinContent(i+1))/hrapAccPPPt1_2Cp->GetBinContent(i+1));
      double sys7 = max(fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Cp2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Cm2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Dp2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt1_2Cp->GetBinContent(i+1)-hrapAccPAPt1_Dm2->GetBinContent(i+1))/hrapAccPAPt1_2Cp->GetBinContent(i+1));
      double sys11 = max(fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Cp2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1),fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Cm2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1));
      double sys12 = max(fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Dp2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1),fabs(hrapAccPP_2Cp->GetBinContent(i+1)-hrapAccPP_Dm2->GetBinContent(i+1))/hrapAccPP_2Cp->GetBinContent(i+1));
      double Sys_rap_PPpt1 = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_rap_PApt1 = max(sys5,max(sys6,max(sys7,sys8)));
      double Sys_rap_PP = max(sys9,max(sys10,max(sys11,sys12)));
      cout<<"PP pt1 A : "<<Form("%0.3f",sys1*100)<<", PP pt1 B : "<<Form("%0.3f",sys2*100)<<", PP pt1 C : "<<Form("%0.3f",sys3*100)<<", PP pt1 D : "<<Form("%0.3f",sys4*100)<<", PP pt1 total : "<<Form("%0.3f",abs(Sys_rap_PPpt1*100))<<endl;
      cout<<"PA pt1 A : "<<Form("%0.3f",sys5*100)<<", PA pt1 B : "<<Form("%0.3f",sys6*100)<<", PA pt1 C : "<<Form("%0.3f",sys7*100)<<", PA pt1 D : "<<Form("%0.3f",sys8*100)<<", PA pt1 total : "<<Form("%0.3f",abs(Sys_rap_PApt1*100))<<endl;
      hrapSysAccPPPt1->SetBinContent(i+1,Sys_rap_PPpt1);
      hrapSysAccPAPt1->SetBinContent(i+1,Sys_rap_PApt1);
      hrapSysAccPP->SetBinContent(i+1,Sys_rap_PP);
    }
    else if(Cat_!=0){
      double Sys_rap_PPpt1 = max(sys1,sys2);
      double Sys_rap_PApt1 = max(sys5,sys6);
      double Sys_rap_PP = max(sys9,sys10);
      cout<<"PP pt1 A : "<<Form("%0.3f",sys1*100)<<", PP pt1 B : "<<Form("%0.3f",sys2*100)<<", PP pt1 total : "<<Form("%0.3f",abs(Sys_rap_PPpt1*100))<<endl;
      cout<<"PA pt1 A : "<<Form("%0.3f",sys5*100)<<", PA pt1 B : "<<Form("%0.3f",sys6*100)<<", PA pt1 total : "<<Form("%0.3f",abs(Sys_rap_PApt1*100))<<endl;
      hrapSysAccPPPt1->SetBinContent(i+1,Sys_rap_PPpt1);
      hrapSysAccPAPt1->SetBinContent(i+1,Sys_rap_PApt1);
      hrapSysAccPP->SetBinContent(i+1,Sys_rap_PP);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hrapAccPPNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Ap2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Am2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1));
    double sys2 = max(fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Bp2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Bm2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1));
    double sys5 = max(fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Ap2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Am2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Bp2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Bm2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1));
    double sys9 = max(fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Ap2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1),fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Am2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1));
    double sys10 = max(fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Bp2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1),fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Bm2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys3 = max(fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Cp2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Cm2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1));
      double sys4 = max(fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Dp2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPPPt2_2Cp->GetBinContent(i+1)-hrapAccPPPt2_Dm2->GetBinContent(i+1))/hrapAccPPPt2_2Cp->GetBinContent(i+1));
      double sys7 = max(fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Cp2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Cm2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Dp2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1),fabs(hrapAccPAPt2_2Cp->GetBinContent(i+1)-hrapAccPAPt2_Dm2->GetBinContent(i+1))/hrapAccPAPt2_2Cp->GetBinContent(i+1));
      double sys11 = max(fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Cp2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1),fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Cm2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1));
      double sys12 = max(fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Dp2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1),fabs(hrapAccPA_2Cp->GetBinContent(i+1)-hrapAccPA_Dm2->GetBinContent(i+1))/hrapAccPA_2Cp->GetBinContent(i+1));
      double Sys_rap_PPpt2 = max(sys1,max(sys2,max(sys3,sys4)));
      double Sys_rap_PApt2 = max(sys5,max(sys6,max(sys7,sys8)));
      double Sys_rap_PA = max(sys9,max(sys10,max(sys11,sys12)));
      cout<<"PP pt2 A : "<<Form("%0.3f",sys1*100)<<", PP pt2 B : "<<Form("%0.3f",sys2*100)<<", PP pt2 C : "<<Form("%0.3f",sys3*100)<<", PP pt2 D : "<<Form("%0.3f",sys4*100)<<", PP pt2 total : "<<Form("%0.3f",abs(Sys_rap_PPpt2*100))<<endl;
      cout<<"PA pt2 A : "<<Form("%0.3f",sys5*100)<<", PA pt2 B : "<<Form("%0.3f",sys6*100)<<", PA pt2 C : "<<Form("%0.3f",sys7*100)<<", PA pt2 D : "<<Form("%0.3f",sys8*100)<<", PA pt2 total : "<<Form("%0.3f",abs(Sys_rap_PApt2*100))<<endl;
      hrapSysAccPPPt2->SetBinContent(i+1,Sys_rap_PPpt2);
      hrapSysAccPAPt2->SetBinContent(i+1,Sys_rap_PApt2);
      hrapSysAccPA->SetBinContent(i+1,Sys_rap_PA);
    }
    else if(Cat_!=0){
      double Sys_rap_PPpt2 = max(sys1,sys2);
      double Sys_rap_PApt2 = max(sys5,sys6);
      double Sys_rap_PA = max(sys9,sys10);
      cout<<"PP pt2 A : "<<Form("%0.3f",sys1*100)<<", PP pt2 B : "<<Form("%0.3f",sys2*100)<<", PP pt2 total : "<<Form("%0.3f",abs(Sys_rap_PPpt2*100))<<endl;
      cout<<"PA pt2 A : "<<Form("%0.3f",sys5*100)<<", PA pt2 B : "<<Form("%0.3f",sys6*100)<<", PA pt2 total : "<<Form("%0.3f",abs(Sys_rap_PApt2*100))<<endl;
      hrapSysAccPPPt2->SetBinContent(i+1,Sys_rap_PPpt2);
      hrapSysAccPAPt2->SetBinContent(i+1,Sys_rap_PApt2);
      hrapSysAccPA->SetBinContent(i+1,Sys_rap_PA);
    }
  }

  cout<<""<<endl;
  for(int i = 0; i < hrapAccPA2binNoW->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"check for pt1 A value :: total : "<<hrapAccPA2bin_2Cp->GetBinContent(i+1)<<" Ap : "<<hrapAccPA2bin_Ap2->GetBinContent(i+1)<<" Am : "<<hrapAccPA2bin_Am2->GetBinContent(i+1)<<endl;
    cout<<"check for pt1 B value :: total : "<<hrapAccPA2bin_2Cp->GetBinContent(i+1)<<" Bp : "<<hrapAccPA2bin_Bp2->GetBinContent(i+1)<<" Bm : "<<hrapAccPA2bin_Bm2->GetBinContent(i+1)<<endl;
    double sys5 = max(fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Ap2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1),fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Am2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Bp2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1),fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Bm2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys7 = max(fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Cp2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1),fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Cm2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Dp2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1),fabs(hrapAccPA2bin_2Cp->GetBinContent(i+1)-hrapAccPA2bin_Dm2->GetBinContent(i+1))/hrapAccPA2bin_2Cp->GetBinContent(i+1));
      double Sys_rap_PA2bin = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"Cross pt1 A : "<<Form("%0.3f",sys5*100)<<", Cross pt1 B : "<<Form("%0.3f",sys6*100)<<", Cross pt1 C : "<<Form("%0.3f",sys7*100)<<", Cross pt1 D : "<<Form("%0.3f",sys8*100)<<", Cross pt1 total : "<<Form("%0.3f",abs(Sys_rap_PA2bin*100))<<endl;
      hrapSysAccPA2bin->SetBinContent(i+1,Sys_rap_PA2bin);
    }
    else if(Cat_!=0){
      double Sys_rap_PA2bin = max(sys5,sys6);
      cout<<"Cross pt1 A : "<<Form("%0.3f",sys5*100)<<", Cross pt1 B : "<<Form("%0.3f",sys6*100)<<", Cross pt1 total : "<<Form("%0.3f",abs(Sys_rap_PA2bin*100))<<endl;
      hrapSysAccPA2bin->SetBinContent(i+1,Sys_rap_PA2bin);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hrapAccCrossNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"check for  A value :: total : "<<hrapAccCross_2Cp->GetBinContent(i+1)<<" Ap : "<<hrapAccCross_Ap2->GetBinContent(i+1)<<" Am : "<<hrapAccCross_Am2->GetBinContent(i+1)<<endl;
    cout<<"check for  B value :: total : "<<hrapAccCross_2Cp->GetBinContent(i+1)<<" Bp : "<<hrapAccCross_Bp2->GetBinContent(i+1)<<" Bm : "<<hrapAccCross_Bm2->GetBinContent(i+1)<<endl;
    double sys5 = max(fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Ap2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1),fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Am2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Bp2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1),fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Bm2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys7 = max(fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Cp2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1),fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Cm2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Dp2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1),fabs(hrapAccCross_2Cp->GetBinContent(i+1)-hrapAccCross_Dm2->GetBinContent(i+1))/hrapAccCross_2Cp->GetBinContent(i+1));
      double Sys_rap_Cross = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"Cross  A : "<<Form("%0.3f",sys5*100)<<", Cross  B : "<<Form("%0.3f",sys6*100)<<", Cross  C : "<<Form("%0.3f",sys7*100)<<", Cross  D : "<<Form("%0.3f",sys8*100)<<", Cross  total : "<<Form("%0.3f",abs(Sys_rap_Cross*100))<<endl;
      hrapSysAccCross->SetBinContent(i+1,Sys_rap_Cross);
    }
    else if(Cat_!=0){
      double Sys_rap_Cross = max(sys5,sys6);
      cout<<"Cross  A : "<<Form("%0.3f",sys5*100)<<", Cross  B : "<<Form("%0.3f",sys6*100)<<", Cross  total : "<<Form("%0.3f",abs(Sys_rap_Cross*100))<<endl;
      hrapSysAccCross->SetBinContent(i+1,Sys_rap_Cross);
    }
  }
  cout<<""<<endl;
  for(int i = 0; i < hrapAccCrossNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"check for pt1 A value :: total : "<<hrapAccCrossPt1_2Cp->GetBinContent(i+1)<<" Ap : "<<hrapAccCrossPt1_Ap2->GetBinContent(i+1)<<" Am : "<<hrapAccCrossPt1_Am2->GetBinContent(i+1)<<endl;
    cout<<"check for pt1 B value :: total : "<<hrapAccCrossPt1_2Cp->GetBinContent(i+1)<<" Bp : "<<hrapAccCrossPt1_Bp2->GetBinContent(i+1)<<" Bm : "<<hrapAccCrossPt1_Bm2->GetBinContent(i+1)<<endl;
    double sys5 = max(fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Ap2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Am2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Bp2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Bm2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys7 = max(fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Cp2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Cm2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Dp2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt1_2Cp->GetBinContent(i+1)-hrapAccCrossPt1_Dm2->GetBinContent(i+1))/hrapAccCrossPt1_2Cp->GetBinContent(i+1));
      double Sys_rap_Crosspt1 = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"Cross pt1 A : "<<Form("%0.3f",sys5*100)<<", Cross pt1 B : "<<Form("%0.3f",sys6*100)<<", Cross pt1 C : "<<Form("%0.3f",sys7*100)<<", Cross pt1 D : "<<Form("%0.3f",sys8*100)<<", Cross pt1 total : "<<Form("%0.3f",abs(Sys_rap_Crosspt1*100))<<endl;
      hrapSysAccCrossPt1->SetBinContent(i+1,Sys_rap_Crosspt1);
    }
    else if(Cat_!=0){
      double Sys_rap_Crosspt1 = max(sys5,sys6);
      cout<<"Cross pt1 A : "<<Form("%0.3f",sys5*100)<<", Cross pt1 B : "<<Form("%0.3f",sys6*100)<<", Cross pt1 total : "<<Form("%0.3f",abs(Sys_rap_Crosspt1*100))<<endl;
      hrapSysAccCrossPt1->SetBinContent(i+1,Sys_rap_Crosspt1);
    }
  }
  for(int i = 0; i < hrapAccCrossNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys5 = max(fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Ap2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Am2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1));
    double sys6 = max(fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Bp2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Bm2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1));
    if(Cat_==0){
      double sys7 = max(fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Cp2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Cm2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1));
      double sys8 = max(fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Dp2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1),fabs(hrapAccCrossPt2_2Cp->GetBinContent(i+1)-hrapAccCrossPt2_Dm2->GetBinContent(i+1))/hrapAccCrossPt2_2Cp->GetBinContent(i+1));
      double Sys_rap_Crosspt2 = max(sys5,max(sys6,max(sys7,sys8)));
      cout<<"Cross pt2 A : "<<Form("%0.3f",sys5*100)<<", Cross pt2 B : "<<Form("%0.3f",sys6*100)<<", Cross pt2 C : "<<Form("%0.3f",sys7*100)<<", Cross pt2 D : "<<Form("%0.3f",sys8*100)<<", Cross pt2 total : "<<Form("%0.3f",abs(Sys_rap_Crosspt2*100))<<endl;
      hrapSysAccCrossPt2->SetBinContent(i+1,Sys_rap_Crosspt2);
    }
    else if(Cat_!=0){
      double Sys_rap_Crosspt2 = max(sys5,sys6);
      cout<<"Cross pt2 A : "<<Form("%0.3f",sys5*100)<<", Cross pt2 B : "<<Form("%0.3f",sys6*100)<<", Cross pt2 total : "<<Form("%0.3f",abs(Sys_rap_Crosspt2*100))<<endl;
      hrapSysAccCrossPt2->SetBinContent(i+1,Sys_rap_Crosspt2);
    }
  }
  //////////////////////////////////Make Systematics/////////////////////////////////
  TFile *fout = new TFile(Form("%d/sys_acceptance_ups%dS_%d.root",date,Cat_+1,date),"recreate");

  fout->cd();

  TH1D * hrapSysAccPPPt1_Clone = (TH1D*)hrapSysAccPPPt1->Clone();
  TH1D * hrapSysAccPPPt2_Clone = (TH1D*)hrapSysAccPPPt2->Clone();
  if(Cat_==0){
    hrapSysAccPPPt1_Clone->SetName("hrapSysAccPPPt1_023");
    hrapSysAccPPPt1_Clone->SetBinContent( 1, hrapSysAccPPPt1->GetBinContent(1) );
    hrapSysAccPPPt1_Clone->SetBinContent( 2, hrapSysAccPPPt1->GetBinContent(2) );
    hrapSysAccPPPt1_Clone->SetBinContent( 3, hrapSysAccPPPt1->GetBinContent(3) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 4, hrapSysAccPPPt1->GetBinContent(4) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 5, hrapSysAccPPPt1->GetBinContent(5) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 6, hrapSysAccPPPt1->GetBinContent(6) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 7, hrapSysAccPPPt1->GetBinContent(7) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 8, hrapSysAccPPPt1->GetBinContent(8) ); 
    hrapSysAccPPPt1_Clone->SetBinContent( 9, hrapSysAccPPPt1->GetBinContent(9) ); 
    hrapSysAccPPPt2_Clone->SetName("hrapSysAccPPPt2_023");
    hrapSysAccPPPt2_Clone->SetBinContent( 1, hrapSysAccPPPt2->GetBinContent(1) );
    hrapSysAccPPPt2_Clone->SetBinContent( 2, hrapSysAccPPPt2->GetBinContent(2) );
    hrapSysAccPPPt2_Clone->SetBinContent( 3, hrapSysAccPPPt2->GetBinContent(3) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 4, hrapSysAccPPPt2->GetBinContent(4) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 5, hrapSysAccPPPt2->GetBinContent(5) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 6, hrapSysAccPPPt2->GetBinContent(6) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 7, hrapSysAccPPPt2->GetBinContent(7) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 8, hrapSysAccPPPt2->GetBinContent(8) ); 
    hrapSysAccPPPt2_Clone->SetBinContent( 9, hrapSysAccPPPt2->GetBinContent(9) ); 
  }
  if(Cat_==1){
    hrapSysAccPPPt1_Clone->SetName("hrapSysAccPPPt1_023");
    hrapSysAccPPPt1_Clone->SetBinContent( 1, hrapSysAccPPPt1->GetBinContent(1) );
    hrapSysAccPPPt1_Clone->SetBinContent( 2, hrapSysAccPPPt1->GetBinContent(2) );
    hrapSysAccPPPt1_Clone->SetBinContent( 3, hrapSysAccPPPt1->GetBinContent(3) );
    hrapSysAccPPPt1_Clone->SetBinContent( 4, hrapSysAccPPPt1->GetBinContent(4) );
    hrapSysAccPPPt1_Clone->SetBinContent( 5, hrapSysAccPPPt1->GetBinContent(5) );
    hrapSysAccPPPt2_Clone->SetName("hrapSysAccPPPt2_023");
    hrapSysAccPPPt2_Clone->SetBinContent( 1, hrapSysAccPPPt2->GetBinContent(1) );
    hrapSysAccPPPt2_Clone->SetBinContent( 2, hrapSysAccPPPt2->GetBinContent(2) );
    hrapSysAccPPPt2_Clone->SetBinContent( 3, hrapSysAccPPPt2->GetBinContent(3) );
    hrapSysAccPPPt2_Clone->SetBinContent( 4, hrapSysAccPPPt2->GetBinContent(4) );
    hrapSysAccPPPt2_Clone->SetBinContent( 5, hrapSysAccPPPt2->GetBinContent(5) );
  }
  if(Cat_==2){
    hrapSysAccPPPt1_Clone->SetName("hrapSysAccPPPt1_023");
    hrapSysAccPPPt1_Clone->SetBinContent( 1, hrapSysAccPPPt1->GetBinContent(1) );
    hrapSysAccPPPt1_Clone->SetBinContent( 2, hrapSysAccPPPt1->GetBinContent(2) );
    hrapSysAccPPPt1_Clone->SetBinContent( 3, hrapSysAccPPPt1->GetBinContent(3) );
    hrapSysAccPPPt2_Clone->SetName("hrapSysAccPPPt2_023");
    hrapSysAccPPPt2_Clone->SetBinContent( 1, hrapSysAccPPPt2->GetBinContent(1) );
    hrapSysAccPPPt2_Clone->SetBinContent( 2, hrapSysAccPPPt2->GetBinContent(2) );
    hrapSysAccPPPt2_Clone->SetBinContent( 3, hrapSysAccPPPt2->GetBinContent(3) );
  }
  hIntSysAccPP->Write();
  hptSysAccPPRap1->Write();
  hptSysAccPPRap2->Write();
  hptSysAccPARap1->Write();
  hptSysAccPARap2->Write();
  hptSysAccPP->Write();
  hrapSysAccPPPt1->Write();
  hrapSysAccPPPt2->Write();
  hrapSysAccPP->Write();
  hIntSysAccPA->Write();
  hrapSysAccPA->Write();
  hrapSysAccPAPt1->Write();
  hrapSysAccPAPt2->Write();
  hrapSysAccPA2bin->Write();
  hptSysAccCross->Write();
  hrapSysAccCross->Write();
  hrapSysAccCrossPt1->Write();
  hrapSysAccCrossPt2->Write();
  hptSysAccPA->Write();
  //if(Cat_==0){
  //  hrapSysXsPA->Write();
  //}
  //hIntSysPA->Write();

  fout->Close();

  out->cd();
  hptAccPPRap1_1->Write();
  hptAccPPRap2_1->Write();
  hptAccPARap1_1->Write();
  hptAccPARap2_1->Write();
  hrapAccPPPt1_1->Write();
  hrapAccPPPt2_1->Write();
  hrapAccPAPt1_1->Write();
  hrapAccPAPt2_1->Write();
  hrapAccPA_1->Write();
  hrapAccPA_2->Write();
  hrapAccCrossPt1_1->Write();
  hrapAccCrossPt2_1->Write();

  hptAccPPRap1_2->Write();
  hptAccPPRap2_2->Write();
  hptAccPARap1_2->Write();
  hptAccPARap2_2->Write();
  hrapAccPPPt1_2->Write();
  hrapAccPPPt2_2->Write();
  hrapAccPAPt1_2->Write();
  hrapAccPAPt2_2->Write();
  hptAccCross_2->Write();
  hrapAccPA2bin_2->Write();
  hrapAccCrossPt1_2->Write();
  hrapAccCrossPt2_2->Write();

  if(Cat_ == 0) {
    hptAccPP_2Cp->SetName("hptAccPP_1S");
    hptAccPA_2Cp->SetName("hptAccPA_1S");
    hptAccPPRap1_2Cp->SetName("hptAccPPRap1_1S");
    hptAccPPRap2_2Cp->SetName("hptAccPPRap2_1S");
    hrapAccPPPt1_2Cp->SetName("hrapAccPPPt1_1S");
    hrapAccPPPt2_2Cp->SetName("hrapAccPPPt2_1S");
    hrapAccPP_2Cp->SetName("hrapAccPP_1S");
    hptAccPARap1_2Cp->SetName("hptAccPARap1_1S");
    hptAccPARap2_2Cp->SetName("hptAccPARap2_1S");
    hrapAccPA2bin_2Cp->SetName("hrapAccPA2bin_1S");
    hrapAccPA_2Cp->SetName("hrapAccPA_1S");
    hrapAccPAPt1_2Cp->SetName("hrapAccPAPt1_1S");
    hrapAccPAPt2_2Cp->SetName("hrapAccPAPt2_1S");
    hptAccCross_2Cp->SetName("hptAccCross_1S");
    hrapAccCross_2Cp->SetName("hrapAccCross_1S");
    hrapAccCrossPt1_2Cp->SetName("hrapAccCrossPt1_1S");
    hrapAccCrossPt2_2Cp->SetName("hrapAccCrossPt2_1S");
    hIntAccPP_2->SetName("hIntAccPP_1S");
    hIntAccPA_2->SetName("hIntAccPA_1S");
  }

  if(Cat_ == 1) {
    hptAccPP_2Cp->SetName("hptAccPP_2S");
    hptAccPA_2Cp->SetName("hptAccPA_2S");
    hptAccPPRap1_2Cp->SetName("hptAccPPRap1_2S");
    hptAccPPRap2_2Cp->SetName("hptAccPPRap2_2S");
    hptAccPARap1_2Cp->SetName("hptAccPARap1_2S");
    hptAccPARap2_2Cp->SetName("hptAccPARap2_2S");
    hrapAccPA2bin_2Cp->SetName("hrapAccPA2bin_2S");
    hrapAccPPPt1_2Cp->SetName("hrapAccPPPt1_2S");
    hrapAccPPPt2_2Cp->SetName("hrapAccPPPt2_2S");
    hrapAccPP_2Cp->SetName("hrapAccPP_2S");
    hrapAccPA_2Cp->SetName("hrapAccPA_2S");
    hrapAccPAPt1_2Cp->SetName("hrapAccPAPt1_2S");
    hrapAccPAPt2_2Cp->SetName("hrapAccPAPt2_2S");
    hptAccCross_2Cp->SetName("hptAccCross_2S");
    hrapAccCross_2Cp->SetName("hrapAccCross_2S");
    hrapAccCrossPt1_2Cp->SetName("hrapAccCrossPt1_2S");
    hrapAccCrossPt2_2Cp->SetName("hrapAccCrossPt2_2S");
    hIntAccPP_2->SetName("hIntAccPP_2S");
    hIntAccPA_2->SetName("hIntAccPA_2S");
  }

  if(Cat_ == 2) {
    hptAccPP_2Cp->SetName("hptAccPP_3S");
    hptAccPA_2Cp->SetName("hptAccPA_3S");
    hptAccPPRap1_2Cp->SetName("hptAccPPRap1_3S");
    hptAccPPRap2_2Cp->SetName("hptAccPPRap2_3S");
    hptAccPARap1_2Cp->SetName("hptAccPARap1_3S");
    hptAccPARap2_2Cp->SetName("hptAccPARap2_3S");
    hrapAccPPPt1_2Cp->SetName("hrapAccPPPt1_3S");
    hrapAccPPPt2_2Cp->SetName("hrapAccPPPt2_3S");
    hrapAccPP_2Cp->SetName("hrapAccPP_3S");
    hrapAccPA_2Cp->SetName("hrapAccPA_3S");
    hrapAccPAPt1_2Cp->SetName("hrapAccPAPt1_3S");
    hrapAccPAPt2_2Cp->SetName("hrapAccPAPt2_3S");
    hrapAccPA2bin_2Cp->SetName("hrapAccPA2bin_3S");
    hptAccCross_2Cp->SetName("hptAccCross_3S");
    hrapAccCross_2Cp->SetName("hrapAccCross_3S");
    hrapAccCrossPt1_2Cp->SetName("hrapAccCrossPt1_3S");
    hrapAccCrossPt2_2Cp->SetName("hrapAccCrossPt2_3S");
    hIntAccPP_2->SetName("hIntAccPP_3S");
    hIntAccPA_2->SetName("hIntAccPA_3S");
  }
  hIntAccPP_2->Write();
  hIntAccPA_2->Write();
  hIntAccCross_2->Write();
  hptAccPPRap1_2Cp->Write();
  hptAccPPRap2_2Cp->Write();
  hptAccPARap1_2Cp->Write();
  hptAccPARap2_2Cp->Write();
  hrapAccPA2bin_2Cp->Write(); 
  hrapAccPAPt1_2Cp->Write(); 
  hrapAccPAPt2_2Cp->Write(); 
  hrapAccPP_2Cp->Write(); 
  hptAccCross_2Cp->Write();
  hrapAccCrossPt1_2Cp->Write();
  hrapAccCrossPt2_2Cp->Write();
  //to make the acc fixed from HIN-16-023
  TH1D * hrapAcc23_Clone = (TH1D*)hrapAcc23->Clone();//make it same with rapidity bins in pA
  if(Cat_==0){
    hrapAcc23->SetName("hrapAccPPPt1_1S_023");
    hrapAcc23->SetBinContent( 1, 2.494152e-01 );
    hrapAcc23->SetBinContent( 2, 2.480554e-01 );
    hrapAcc23->SetBinContent( 3, 2.475569e-01 );
    hrapAcc23->SetBinContent( 4, 2.460665e-01 );
    hrapAcc23->SetBinContent( 5, 2.211435e-01 );
    hrapAcc23->SetBinContent( 6, 9.887582e-02 );
    hrapAcc23->Write();
  }
  if(Cat_==1){
    hrapAcc23->SetName("hrapAccPPPt1_2S_023");
    hrapAcc23->SetBinContent( 1, 3.095539e-01 );
    hrapAcc23->SetBinContent( 2, 3.073197e-01 );
    hrapAcc23->SetBinContent( 3, 1.978147e-01 );
    hrapAcc23->Write();
  }
  if(Cat_==2){
    hrapAcc23->SetName("hrapAccPPPt1_3S_023");
    hrapAcc23->SetBinContent( 1, 3.478020e-01 );
    hrapAcc23->SetBinContent( 2, 2.651429e-01 );
    hrapAcc23->Write();
  }

  //hptAccPPRap1_Ap2->Write();
  //hptAccPPRap1_Bp2->Write();
  //hptAccPPRap1_Am2->Write();
  //hptAccPPRap1_Bm2->Write();
  //hrapAccPPPt1_Ap2->Write();
  //hrapAccPPPt1_Bp2->Write();
  //hrapAccPPPt1_Am2->Write();
  //hrapAccPPPt1_Bm2->Write();
  //hptAccPPRap2_Ap2->Write();
  //hptAccPPRap2_Bp2->Write();
  //hptAccPPRap2_Am2->Write();
  //hptAccPPRap2_Bm2->Write();
  //hrapAccPPPt2_Ap2->Write();
  //hrapAccPPPt2_Bp2->Write();
  //hrapAccPPPt2_Am2->Write();
  //hrapAccPPPt2_Bm2->Write();
  //if(Cat_==0){
  //  hptAccPPRap1_Cp2->Write();
  //  hptAccPPRap1_Dp2->Write();
  //  hptAccPPRap1_Cm2->Write();
  //  hptAccPPRap1_Dm2->Write();
  //  hrapAccPPPt1_Cp2->Write();
  //  hrapAccPPPt1_Dp2->Write();
  //  hrapAccPPPt1_Cm2->Write();
  //  hrapAccPPPt1_Dm2->Write();
  //  hptAccPPRap2_Cp2->Write();
  //  hptAccPPRap2_Dp2->Write();
  //  hptAccPPRap2_Cm2->Write();
  //  hptAccPPRap2_Dm2->Write();
  //  hrapAccPPPt2_Cp2->Write();
  //  hrapAccPPPt2_Dp2->Write();
  //  hrapAccPPPt2_Cm2->Write();
  //  hrapAccPPPt2_Dm2->Write();
  //}
  //hptAccPARap1_Ap2->Write();
  //hptAccPARap1_Bp2->Write();
  //hptAccPARap1_Am2->Write();
  //hptAccPARap1_Bm2->Write();
  //hptAccPARap2_Ap2->Write();
  //hptAccPARap2_Bp2->Write();
  //hptAccPARap2_Am2->Write();
  //hptAccPARap2_Bm2->Write();
  //hrapAccPAPt1_Ap2->Write();
  //hrapAccPAPt1_Bp2->Write();
  //hrapAccPAPt1_Am2->Write();
  //hrapAccPAPt1_Bm2->Write();
  //hrapAccPAPt2_Ap2->Write();
  //hrapAccPAPt2_Bp2->Write();
  //hrapAccPAPt2_Am2->Write();
  //hrapAccPAPt2_Bm2->Write();
  //hrapAccCrossPt1_Ap2->Write();
  //hrapAccCrossPt1_Bp2->Write();
  //hrapAccCrossPt1_Am2->Write();
  //hrapAccCrossPt1_Bm2->Write();
  //hrapAccCrossPt2_Ap2->Write();
  //hrapAccCrossPt2_Bp2->Write();
  //hrapAccCrossPt2_Am2->Write();
  //hrapAccCrossPt2_Bm2->Write();
  //if(Cat_==0){
  //  hptAccPARap1_Cp2->Write();
  //  hptAccPARap1_Dp2->Write();
  //  hptAccPARap1_Cm2->Write();
  //  hptAccPARap1_Dm2->Write();
  //  hptAccPARap2_Cp2->Write();
  //  hptAccPARap2_Dp2->Write();
  //  hptAccPARap2_Cm2->Write();
  //  hptAccPARap2_Dm2->Write();
  //  hrapAccPAPt1_Cp2->Write();
  //  hrapAccPAPt1_Dp2->Write();
  //  hrapAccPAPt1_Cm2->Write();
  //  hrapAccPAPt1_Dm2->Write();
  //  hrapAccPAPt2_Cp2->Write();
  //  hrapAccPAPt2_Dp2->Write();
  //  hrapAccPAPt2_Cm2->Write();
  //  hrapAccPAPt2_Dm2->Write();
  //  hrapAccCrossPt1_Cp2->Write();
  //  hrapAccCrossPt1_Dp2->Write();
  //  hrapAccCrossPt1_Cm2->Write();
  //  hrapAccCrossPt1_Dm2->Write();
  //  hrapAccCrossPt2_Cp2->Write();
  //  hrapAccCrossPt2_Dp2->Write();
  //  hrapAccCrossPt2_Cm2->Write();
  //  hrapAccCrossPt2_Dm2->Write();
  //}

  cout<<"CHECK"<<endl;

  TLegend *leg11 = new TLegend(0.15,0.62,0.83,0.89);
  leg11->AddEntry(hptAccPP_2Cp,Form("Nominal(pp %dS, -1.93<y_{CM}<1.93)", Cat_+1),"l");
  leg11->AddEntry(max_ptPP,Form("Largest(pp %dS)", Cat_+1),"lp");
  leg11->AddEntry(min_ptPP,Form("Smallest(pp %dS)", Cat_+1),"lp");
  leg11->SetLineColor(kWhite);
  leg11->SetFillStyle(0);
  leg11->SetTextSize(0.045);

  TLegend *leg12 = new TLegend(0.13,0.62,0.83,0.89);
  leg12->AddEntry(hptAccPA_2Cp,Form("Nominal(pA %dS, -1.93<y_{CM}<1.93)", Cat_+1),"l");
  leg12->AddEntry(max_ptPA,Form("Largest(pA %dS)", Cat_+1),"p");
  leg12->AddEntry(min_ptPA,Form("Smallest(pA %dS)", Cat_+1),"p");
  leg12->SetLineColor(kWhite);
  leg12->SetFillStyle(0);
  leg12->SetTextSize(0.045);

  TLegend *leg13 = new TLegend(0.15,0.62,0.83,0.89);
  leg13->AddEntry(hrapAccPP_2Cp,Form("Nominal(pp %dS, 0<pt<30)", Cat_+1),"l");
  leg13->AddEntry(max_rapPP,Form("Largest(pp %dS)", Cat_+1),"p");
  leg13->AddEntry(min_rapPP,Form("Smallest(pp %dS)", Cat_+1),"p");
  leg13->SetLineColor(kWhite);
  leg13->SetFillStyle(0);
  leg13->SetTextSize(0.045);

  TLegend *leg14 = new TLegend(0.15,0.62,0.83,0.89);
  leg14->AddEntry(hrapAccPA_2Cp,Form("Nominal(pA %dS, 0<pt<30)", Cat_+1),"l");
  leg14->AddEntry(max_rapPA,Form("Largest(pA %dS)", Cat_+1),"p");
  leg14->AddEntry(min_rapPA,Form("Smallest(pA %dS)", Cat_+1),"p");
  //leg14->AddEntry(hrapAccPPAN_2Cp,Form("HIN-18-005(pp %dS)", Cat_+1),"l");
  //leg14->AddEntry(hrapAcc23,Form("HIN-16-023(pp %dS)", Cat_+1),"p");
  leg14->SetLineColor(kWhite);
  leg14->SetFillStyle(0);
  leg14->SetTextSize(0.045);

  TLegend *leg17 = new TLegend(0.15,0.62,0.83,0.89);
  leg17->AddEntry(hrapAccPPAN_2Cp,Form("HIN-18-005(pp %dS)", Cat_+1),"l");
  leg17->AddEntry(hrapAcc23,Form("HIN-16-023(pp %dS)", Cat_+1),"p");
  leg17->SetLineColor(kWhite);
  leg17->SetFillStyle(0);
  leg17->SetTextSize(0.045);

  TLegend *leg1 = new TLegend(0.15,0.62,0.83,0.89);
  leg1->AddEntry(hptAccPPRap1_2Cp,Form("Nominal(pp %dS, -1.93<y_{CM}<0)", Cat_+1),"l");
  leg1->AddEntry(max_ptPPRap1,Form("Largest(pp %dS)", Cat_+1),"p");
  leg1->AddEntry(min_ptPPRap1,Form("Smallest(pp %dS)", Cat_+1),"p");
  leg1->SetLineColor(kWhite);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.045);

  TLegend *leg2 = new TLegend(0.15,0.62,0.83,0.89);
  leg2->AddEntry(hptAccPPRap2_2Cp,Form("Nominal(pp %dS, 0<y_{CM}<1.93)", Cat_+1),"l");
  leg2->AddEntry(max_ptPPRap2,Form("Largest(pp %dS)", Cat_+1),"p");
  leg2->AddEntry(min_ptPPRap2,Form("Smallest(pp %dS)", Cat_+1),"p");
  leg2->SetLineColor(kWhite);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.045);

  TLegend *leg3 = new TLegend(0.13,0.62,0.83,0.89);
  leg3->AddEntry(hptAccPARap1_2Cp,Form("Nominal(pA %dS, -1.93<y_{CM}<0)", Cat_+1),"l");
  leg3->AddEntry(max_ptPARap1,Form("Largest(pA %dS)", Cat_+1),"p");
  leg3->AddEntry(min_ptPARap1,Form("Smallest(pA %dS)", Cat_+1),"p");
  leg3->SetLineColor(kWhite);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.045);

  TLegend *leg4 = new TLegend(0.13,0.62,0.83,0.89);
  leg4->AddEntry(hptAccPARap2_2Cp,Form("Nominal(pA %dS, 0<y_{CM}<1.93)", Cat_+1),"l");
  leg4->AddEntry(max_ptPARap2,Form("Largest(pA %dS)", Cat_+1),"p");
  leg4->AddEntry(min_ptPARap2,Form("Smallest(pA %dS)", Cat_+1),"p");
  leg4->SetLineColor(kWhite);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.045);

  TLegend *leg5 = new TLegend(0.13,0.62,0.83,0.89);
  leg5->AddEntry(hrapAccPPPt1_2Cp,Form("Nominal(pp %dS, 0<pt<6GeV/c)",Cat_+1),"l");
  leg5->AddEntry(max_rapPPPt1,Form("Largest(pp %dS)",Cat_+1),"p");
  leg5->AddEntry(min_rapPPPt1,Form("Smallest(pp %dS)",Cat_+1),"p");
  leg5->SetLineColor(kWhite);
  leg5->SetFillStyle(0);
  leg5->SetTextSize(0.045);

  TLegend *leg6 = new TLegend(0.13,0.62,0.83,0.89);
  leg6->AddEntry(hrapAccPPPt2_2Cp,Form("Nominal(pp %dS, 6<pt<30GeV/c)",Cat_+1),"l");
  leg6->AddEntry(max_rapPPPt2,Form("Largest(pp %dS)",Cat_+1),"p");
  leg6->AddEntry(min_rapPPPt2,Form("Smallest(pp %dS)",Cat_+1),"p");
  leg6->SetLineColor(kWhite);
  leg6->SetFillStyle(0);
  leg6->SetTextSize(0.045);

  TLegend *leg7 = new TLegend(0.13,0.62,0.83,0.89);
  leg7->AddEntry(hrapAccPAPt1_2Cp,Form("Nominal(pA %dS, 0<pt<6GeV/c)",Cat_+1),"l");
  leg7->AddEntry(max_rapPAPt1,Form("Largest(pA %dS)",Cat_+1),"p");
  leg7->AddEntry(min_rapPAPt1,Form("Smallest(pA %dS)",Cat_+1),"p");
  leg7->SetLineColor(kWhite);
  leg7->SetFillStyle(0);
  leg7->SetTextSize(0.045);

  TLegend *leg8 = new TLegend(0.13,0.62,0.83,0.89);
  leg8->AddEntry(hrapAccPAPt2_2Cp,Form("Nominal(pA %dS, 6<pt<30GeV/c)",Cat_+1),"l");
  leg8->AddEntry(max_rapPAPt2,Form("Largest(pA %dS)",Cat_+1),"p");
  leg8->AddEntry(min_rapPAPt2,Form("Smallest(pA %dS)",Cat_+1),"p");
  leg8->SetLineColor(kWhite);
  leg8->SetFillStyle(0);
  leg8->SetTextSize(0.045);

  TLegend *leg9 = new TLegend(0.13,0.64,0.83,0.89);
  leg9->AddEntry(hrapAccCross_2Cp,Form("Nominal(pA %dS, 0<pt<30GeV/c)",Cat_+1),"l");
  leg9->AddEntry(max_rapCross,Form("Largest(pA %dS)",Cat_+1),"p");
  leg9->AddEntry(min_rapCross,Form("Smallest(pA %dS)",Cat_+1),"p");
  leg9->SetLineColor(kWhite);
  leg9->SetFillStyle(0);
  leg9->SetTextSize(0.045);

  TLegend *leg10 = new TLegend(0.13,0.62,0.83,0.89);
  leg10->AddEntry(hptAccCross_2Cp,Form("Nominal(pA %dS, -2.87<y_{CM}<1.93)",Cat_+1),"l");
  leg10->AddEntry(max_ptCross,Form("Largest(pA %dS)",Cat_+1),"p");
  leg10->AddEntry(min_ptCross,Form("Smallest(pA %dS)",Cat_+1),"p");
  leg10->SetLineColor(kWhite);
  leg10->SetFillStyle(0);
  leg10->SetTextSize(0.045);

  cout<<"CHECK"<<endl;
  for(int i=0; i<hptAccPPRap1_1->GetXaxis()->GetNbins(); i++){
    if(Cat_==0){
      double acc1 = max(hptAccPPRap1_Ap2->GetBinContent(i+1),max(hptAccPPRap1_Bp2->GetBinContent(i+1),max(hptAccPPRap1_Cp2->GetBinContent(i+1),hptAccPPRap1_Dp2->GetBinContent(i+1))));
      double acc2 = min(hptAccPPRap1_Am2->GetBinContent(i+1),min(hptAccPPRap1_Bm2->GetBinContent(i+1),min(hptAccPPRap1_Cm2->GetBinContent(i+1),hptAccPPRap1_Dm2->GetBinContent(i+1))));
      double acc3 = max(hptAccPPRap2_Ap2->GetBinContent(i+1),max(hptAccPPRap2_Bp2->GetBinContent(i+1),max(hptAccPPRap2_Cp2->GetBinContent(i+1),hptAccPPRap2_Dp2->GetBinContent(i+1))));
      double acc4 = min(hptAccPPRap2_Am2->GetBinContent(i+1),min(hptAccPPRap2_Bm2->GetBinContent(i+1),min(hptAccPPRap2_Cm2->GetBinContent(i+1),hptAccPPRap2_Dm2->GetBinContent(i+1))));
      double acc5 = max(hptAccPARap1_Ap2->GetBinContent(i+1),max(hptAccPARap1_Bp2->GetBinContent(i+1),max(hptAccPARap1_Cp2->GetBinContent(i+1),hptAccPARap1_Dp2->GetBinContent(i+1))));
      double acc6 = min(hptAccPARap1_Am2->GetBinContent(i+1),min(hptAccPARap1_Bm2->GetBinContent(i+1),min(hptAccPARap1_Cm2->GetBinContent(i+1),hptAccPARap1_Dm2->GetBinContent(i+1))));
      double acc7 = max(hptAccPARap2_Ap2->GetBinContent(i+1),max(hptAccPARap2_Bp2->GetBinContent(i+1),max(hptAccPARap2_Cp2->GetBinContent(i+1),hptAccPARap2_Dp2->GetBinContent(i+1))));
      double acc8 = min(hptAccPARap2_Am2->GetBinContent(i+1),min(hptAccPARap2_Bm2->GetBinContent(i+1),min(hptAccPARap2_Cm2->GetBinContent(i+1),hptAccPARap2_Dm2->GetBinContent(i+1))));
      double acc9 =  max(hptAccPP_Ap2->GetBinContent(i+1),max(hptAccPP_Bp2->GetBinContent(i+1),max(hptAccPP_Cp2->GetBinContent(i+1),hptAccPP_Dp2->GetBinContent(i+1))));
      double acc10 = min(hptAccPP_Am2->GetBinContent(i+1),min(hptAccPP_Bm2->GetBinContent(i+1),min(hptAccPP_Cm2->GetBinContent(i+1),hptAccPP_Dm2->GetBinContent(i+1))));
      double acc11 = max(hptAccPA_Ap2->GetBinContent(i+1),max(hptAccPA_Bp2->GetBinContent(i+1),max(hptAccPA_Cp2->GetBinContent(i+1),hptAccPA_Dp2->GetBinContent(i+1))));
      double acc12 = min(hptAccPA_Am2->GetBinContent(i+1),min(hptAccPA_Bm2->GetBinContent(i+1),min(hptAccPA_Cm2->GetBinContent(i+1),hptAccPA_Dm2->GetBinContent(i+1))));
      cout<<"################ "<<i+1<<" ################ : "<< acc1 <<" "<<endl;
      cout<<"################ "<<i+1<<" ################ : "<< acc3 <<" "<<endl;
      max_ptPPRap1->SetBinContent(i+1,acc1);
      min_ptPPRap1->SetBinContent(i+1,acc2);
      getMaxDev(hptAccPPRap1_2Cp, max_ptPPRap1, min_ptPPRap1, MD_ptPPRap1);
      max_ptPPRap2->SetBinContent(i+1,acc3);
      min_ptPPRap2->SetBinContent(i+1,acc4);
      getMaxDev(hptAccPPRap2_2Cp, max_ptPPRap2, min_ptPPRap2, MD_ptPPRap2);
      max_ptPARap1->SetBinContent(i+1,acc5);
      min_ptPARap1->SetBinContent(i+1,acc6);
      getMaxDev(hptAccPARap1_2Cp, max_ptPARap1, min_ptPARap1, MD_ptPARap1);
      max_ptPARap2->SetBinContent(i+1,acc7);
      min_ptPARap2->SetBinContent(i+1,acc8);
      getMaxDev(hptAccPARap2_2Cp, max_ptPARap2, min_ptPARap2, MD_ptPARap2);
      max_ptPP->SetBinContent(i+1,acc9);
      min_ptPP->SetBinContent(i+1,acc10);
      getMaxDev(hptAccPPRap1_2Cp, max_ptPP, min_ptPP, MD_ptPP);
      max_ptPA->SetBinContent(i+1,acc11);
      min_ptPA->SetBinContent(i+1,acc12);
      getMaxDev(hptAccPARap2_2Cp, max_ptPA, min_ptPA, MD_ptPA);
    }
    if(Cat_!=0){
      double acc1 = max(hptAccPPRap1_Ap2->GetBinContent(i+1),hptAccPPRap1_Bp2->GetBinContent(i+1));
      double acc2 = min(hptAccPPRap1_Am2->GetBinContent(i+1),hptAccPPRap1_Bm2->GetBinContent(i+1));
      double acc3 = max(hptAccPPRap2_Ap2->GetBinContent(i+1),hptAccPPRap2_Bp2->GetBinContent(i+1));
      double acc4 = min(hptAccPPRap2_Am2->GetBinContent(i+1),hptAccPPRap2_Bm2->GetBinContent(i+1));
      double acc5 = max(hptAccPARap1_Ap2->GetBinContent(i+1),hptAccPARap1_Bp2->GetBinContent(i+1));
      double acc6 = min(hptAccPARap1_Am2->GetBinContent(i+1),hptAccPARap1_Bm2->GetBinContent(i+1));
      double acc7 = max(hptAccPARap2_Ap2->GetBinContent(i+1),hptAccPARap2_Bp2->GetBinContent(i+1));
      double acc8 = min(hptAccPARap2_Am2->GetBinContent(i+1),hptAccPARap2_Bm2->GetBinContent(i+1));
      double acc9 = max(hptAccPP_Ap2->GetBinContent(i+1),hptAccPP_Bp2->GetBinContent(i+1));
      double acc10 = min(hptAccPP_Am2->GetBinContent(i+1),hptAccPP_Bm2->GetBinContent(i+1));
      double acc11 = max(hptAccPA_Ap2->GetBinContent(i+1),hptAccPA_Bp2->GetBinContent(i+1));
      double acc12 = min(hptAccPA_Am2->GetBinContent(i+1),hptAccPA_Bm2->GetBinContent(i+1));
      max_ptPPRap1->SetBinContent(i+1,acc1);
      min_ptPPRap1->SetBinContent(i+1,acc2);
      getMaxDev(hptAccPPRap1_2Cp, max_ptPPRap1, min_ptPPRap1, MD_ptPPRap1);
      max_ptPPRap2->SetBinContent(i+1,acc3);
      min_ptPPRap2->SetBinContent(i+1,acc4);
      getMaxDev(hptAccPPRap2_2Cp, max_ptPPRap2, min_ptPPRap2, MD_ptPPRap2);
      max_ptPARap1->SetBinContent(i+1,acc5);
      min_ptPARap1->SetBinContent(i+1,acc6);
      getMaxDev(hptAccPARap1_2Cp, max_ptPARap1, min_ptPARap1, MD_ptPARap1);
      max_ptPARap2->SetBinContent(i+1,acc7);
      min_ptPARap2->SetBinContent(i+1,acc8);
      getMaxDev(hptAccPARap2_2Cp, max_ptPARap2, min_ptPARap2, MD_ptPARap2);
      max_ptPP->SetBinContent(i+1,acc9);
      min_ptPP->SetBinContent(i+1,acc10);
      getMaxDev(hptAccPP_2Cp, max_ptPP, min_ptPP, MD_ptPP);
      max_ptPA->SetBinContent(i+1,acc11);
      min_ptPA->SetBinContent(i+1,acc12);
      getMaxDev(hptAccPA_2Cp, max_ptPA, min_ptPA, MD_ptPA);
    }
  }
  for(int i=0; i<hrapAccPPPt1_1->GetXaxis()->GetNbins(); i++){
    if(Cat_==0){
      double acc1 = max(hrapAccPPPt1_Ap2->GetBinContent(i+1),max(hrapAccPPPt1_Bp2->GetBinContent(i+1),max(hrapAccPPPt1_Cp2->GetBinContent(i+1),hrapAccPPPt1_Dp2->GetBinContent(i+1))));
      double acc2 = min(hrapAccPPPt1_Am2->GetBinContent(i+1),min(hrapAccPPPt1_Bm2->GetBinContent(i+1),min(hrapAccPPPt1_Cm2->GetBinContent(i+1),hrapAccPPPt1_Dm2->GetBinContent(i+1))));
      double acc3 = max(hrapAccPPPt2_Ap2->GetBinContent(i+1),max(hrapAccPPPt2_Bp2->GetBinContent(i+1),max(hrapAccPPPt2_Cp2->GetBinContent(i+1),hrapAccPPPt2_Dp2->GetBinContent(i+1))));
      double acc4 = min(hrapAccPPPt2_Am2->GetBinContent(i+1),min(hrapAccPPPt2_Bm2->GetBinContent(i+1),min(hrapAccPPPt2_Cm2->GetBinContent(i+1),hrapAccPPPt2_Dm2->GetBinContent(i+1))));
      double acc5 = max(hrapAccPAPt1_Ap2->GetBinContent(i+1),max(hrapAccPAPt1_Bp2->GetBinContent(i+1),max(hrapAccPAPt1_Cp2->GetBinContent(i+1),hrapAccPAPt1_Dp2->GetBinContent(i+1))));
      double acc6 = min(hrapAccPAPt1_Am2->GetBinContent(i+1),min(hrapAccPAPt1_Bm2->GetBinContent(i+1),min(hrapAccPAPt1_Cm2->GetBinContent(i+1),hrapAccPAPt1_Dm2->GetBinContent(i+1))));
      double acc7 = max(hrapAccPAPt2_Ap2->GetBinContent(i+1),max(hrapAccPAPt2_Bp2->GetBinContent(i+1),max(hrapAccPAPt2_Cp2->GetBinContent(i+1),hrapAccPAPt2_Dp2->GetBinContent(i+1))));
      double acc8 = min(hrapAccPAPt2_Am2->GetBinContent(i+1),min(hrapAccPAPt2_Bm2->GetBinContent(i+1),min(hrapAccPAPt2_Cm2->GetBinContent(i+1),hrapAccPAPt2_Dm2->GetBinContent(i+1))));
      double acc9 = max(hrapAccPP_Ap2->GetBinContent(i+1),max(hrapAccPP_Bp2->GetBinContent(i+1),max(hrapAccPP_Cp2->GetBinContent(i+1),hrapAccPP_Dp2->GetBinContent(i+1))));
      double acc10 = min(hrapAccPP_Am2->GetBinContent(i+1),min(hrapAccPP_Bm2->GetBinContent(i+1),min(hrapAccPP_Cm2->GetBinContent(i+1),hrapAccPP_Dm2->GetBinContent(i+1))));
      double acc11 = max(hrapAccPA_Ap2->GetBinContent(i+1),max(hrapAccPA_Bp2->GetBinContent(i+1),max(hrapAccPA_Cp2->GetBinContent(i+1),hrapAccPA_Dp2->GetBinContent(i+1))));
      double acc12 = min(hrapAccPA_Am2->GetBinContent(i+1),min(hrapAccPA_Bm2->GetBinContent(i+1),min(hrapAccPA_Cm2->GetBinContent(i+1),hrapAccPA_Dm2->GetBinContent(i+1))));
      cout<<"################ "<<i+1<<" ################ : "<< acc1 <<" "<<endl;
      cout<<"################ "<<i+1<<" ################ : "<< acc3 <<" "<<endl;
      max_rapPPPt1->SetBinContent(i+1,acc1);
      min_rapPPPt1->SetBinContent(i+1,acc2);
      getMaxDev(hrapAccPPPt1_2Cp, max_rapPPPt1, min_rapPPPt1, MD_rapPPPt1);
      max_rapPPPt2->SetBinContent(i+1,acc3);
      min_rapPPPt2->SetBinContent(i+1,acc4);
      getMaxDev(hrapAccPPPt2_2Cp, max_rapPPPt2, min_rapPPPt2, MD_rapPPPt2);
      max_rapPAPt1->SetBinContent(i+1,acc5);
      min_rapPAPt1->SetBinContent(i+1,acc6);
      getMaxDev(hrapAccPAPt1_2Cp, max_rapPAPt1, min_rapPAPt1, MD_rapPAPt1);
      max_rapPAPt2->SetBinContent(i+1,acc7);
      min_rapPAPt2->SetBinContent(i+1,acc8);
      getMaxDev(hrapAccPAPt2_2Cp, max_rapPAPt2, min_rapPAPt2, MD_rapPAPt2);
      max_rapPP->SetBinContent(i+1,acc9);
      min_rapPP->SetBinContent(i+1,acc10);
      getMaxDev(hrapAccPP_2Cp, max_rapPP, min_rapPP, MD_rapPP);
      max_rapPA->SetBinContent(i+1,acc11);
      min_rapPA->SetBinContent(i+1,acc12);
      getMaxDev(hrapAccPA_2Cp, max_rapPA, min_rapPA, MD_rapPA);
    }
    if(Cat_!=0){
      double acc1 = max(hrapAccPPPt1_Ap2->GetBinContent(i+1),hrapAccPPPt1_Bp2->GetBinContent(i+1));
      double acc2 = min(hrapAccPPPt1_Am2->GetBinContent(i+1),hrapAccPPPt1_Bm2->GetBinContent(i+1));
      double acc3 = max(hrapAccPPPt2_Ap2->GetBinContent(i+1),hrapAccPPPt2_Bp2->GetBinContent(i+1));
      double acc4 = min(hrapAccPPPt2_Am2->GetBinContent(i+1),hrapAccPPPt2_Bm2->GetBinContent(i+1));
      double acc5 = max(hrapAccPAPt1_Ap2->GetBinContent(i+1),hrapAccPAPt1_Bp2->GetBinContent(i+1));
      double acc6 = min(hrapAccPAPt1_Am2->GetBinContent(i+1),hrapAccPAPt1_Bm2->GetBinContent(i+1));
      double acc7 = max(hrapAccPAPt2_Ap2->GetBinContent(i+1),hrapAccPAPt2_Bp2->GetBinContent(i+1));
      double acc8 = min(hrapAccPAPt2_Am2->GetBinContent(i+1),hrapAccPAPt2_Bm2->GetBinContent(i+1));
      double acc9 = max(hrapAccPP_Ap2->GetBinContent(i+1),hrapAccPP_Bp2->GetBinContent(i+1));
      double acc10 = min(hrapAccPP_Am2->GetBinContent(i+1),hrapAccPP_Bm2->GetBinContent(i+1));
      double acc11 = max(hrapAccPA_Ap2->GetBinContent(i+1),hrapAccPA_Bp2->GetBinContent(i+1));
      double acc12 = min(hrapAccPA_Am2->GetBinContent(i+1),hrapAccPA_Bm2->GetBinContent(i+1));
      cout<<"################ "<<i+1<<" ################ : "<< acc1 <<" "<<endl;
      cout<<"################ "<<i+1<<" ################ : "<< acc3 <<" "<<endl;
      max_rapPPPt1->SetBinContent(i+1,acc1);
      min_rapPPPt1->SetBinContent(i+1,acc2);
      getMaxDev(hrapAccPPPt1_2Cp, max_rapPPPt1, min_rapPPPt1, MD_rapPPPt1);
      max_rapPPPt2->SetBinContent(i+1,acc3);
      min_rapPPPt2->SetBinContent(i+1,acc4);
      getMaxDev(hrapAccPPPt2_2Cp, max_rapPPPt2, min_rapPPPt2, MD_rapPPPt2);
      max_rapPAPt1->SetBinContent(i+1,acc5);
      min_rapPAPt1->SetBinContent(i+1,acc6);
      getMaxDev(hrapAccPAPt1_2Cp, max_rapPAPt1, min_rapPAPt1, MD_rapPAPt1);
      max_rapPAPt2->SetBinContent(i+1,acc7);
      min_rapPAPt2->SetBinContent(i+1,acc8);
      getMaxDev(hrapAccPAPt2_2Cp, max_rapPAPt2, min_rapPAPt2, MD_rapPAPt2);
      max_rapPP->SetBinContent(i+1,acc9);
      min_rapPP->SetBinContent(i+1,acc10);
      getMaxDev(hrapAccPP_2Cp, max_rapPP, min_rapPP, MD_rapPP);
      max_rapPA->SetBinContent(i+1,acc11);
      min_rapPA->SetBinContent(i+1,acc12);
      getMaxDev(hrapAccPA_2Cp, max_rapPA, min_rapPA, MD_rapPA);
    }
  }
  cout<<""<<endl;
  for(int i=0; i<hrapAccCross_1->GetXaxis()->GetNbins(); i++){
    if(Cat_==0){
      double acc1 = max(hrapAccCross_Ap2->GetBinContent(i+1),max(hrapAccCross_Bp2->GetBinContent(i+1),max(hrapAccCross_Cp2->GetBinContent(i+1),hrapAccCross_Dp2->GetBinContent(i+1))));
      double acc2 = min(hrapAccCross_Am2->GetBinContent(i+1),min(hrapAccCross_Bm2->GetBinContent(i+1),min(hrapAccCross_Cm2->GetBinContent(i+1),hrapAccCross_Dm2->GetBinContent(i+1))));
      max_rapCross->SetBinContent(i+1,acc1);
      min_rapCross->SetBinContent(i+1,acc2);
      getMaxDev(hrapAccCross_2Cp, max_rapCross, min_rapCross, MD_rapCross);
      cout<<"Max Deviation : "<<MD_rapCross->GetBinContent(i+1)<<"ref : "<<hrapAccCross_2Cp->GetBinContent(i+1)<<", "<<hrapAccCross_Bp2->GetBinContent(i+1)<<endl;
    }
    if(Cat_!=0){
      double acc1 = max(hrapAccCross_Ap2->GetBinContent(i+1),hrapAccCross_Bp2->GetBinContent(i+1));
      double acc2 = min(hrapAccCross_Am2->GetBinContent(i+1),hrapAccCross_Bm2->GetBinContent(i+1));
      max_rapCross->SetBinContent(i+1,acc1);
      min_rapCross->SetBinContent(i+1,acc2);
      getMaxDev(hrapAccCross_2Cp, max_rapCross, min_rapCross, MD_rapCross);
    }
  }
  cout<<""<<endl;
  for(int i=0; i<hptAccCross_1->GetXaxis()->GetNbins(); i++){
    if(Cat_==0){
      double acc1 = max(hptAccCross_Ap2->GetBinContent(i+1),max(hptAccCross_Bp2->GetBinContent(i+1),max(hptAccCross_Cp2->GetBinContent(i+1),hptAccCross_Dp2->GetBinContent(i+1))));
      double acc2 = min(hptAccCross_Am2->GetBinContent(i+1),min(hptAccCross_Bm2->GetBinContent(i+1),min(hptAccCross_Cm2->GetBinContent(i+1),hptAccCross_Dm2->GetBinContent(i+1))));
      max_ptCross->SetBinContent(i+1,acc1);
      min_ptCross->SetBinContent(i+1,acc2);
      getMaxDev(hptAccCross_2Cp, max_ptCross, min_ptCross, MD_ptCross);
      cout<<"Max Deviation : "<<MD_ptCross->GetBinContent(i+1)<<"ref : "<<hptAccCross_2Cp->GetBinContent(i+1)<<", "<<hptAccCross_Bp2->GetBinContent(i+1)<<endl;
    }
    if(Cat_!=0){
      double acc1 = max(hptAccCross_Ap2->GetBinContent(i+1),hptAccCross_Bp2->GetBinContent(i+1));
      double acc2 = min(hptAccCross_Am2->GetBinContent(i+1),hptAccCross_Bm2->GetBinContent(i+1));
      max_ptCross->SetBinContent(i+1,acc1);
      min_ptCross->SetBinContent(i+1,acc2);
      getMaxDev(hptAccCross_2Cp, max_ptCross, min_ptCross, MD_ptCross);
    }
  }

  ///////////////////////////////////////////////////////////////////////////Draw//////////////////////////////////////////////////////////////////////////////
  hptAccPP_2Cp->SetLineColor(kGreen+2);
  hrapAccPP_2Cp->SetLineColor(kGreen+2);
  hptAccPA_2Cp->SetLineColor(kGreen+2);
  hrapAccPA_2Cp->SetLineColor(kGreen+2);
  hrapAccPPPt1_2Cp->SetLineColor(kGreen+2);
  hrapAccPPPt2_2Cp->SetLineColor(kGreen+2);
  hrapAccPAPt1_2Cp->SetLineColor(kGreen+2);
  hrapAccPAPt2_2Cp->SetLineColor(kGreen+2);
  hptAccPPRap1_2Cp->SetLineColor(kGreen+2);
  hptAccPPRap2_2Cp->SetLineColor(kGreen+2);
  hptAccPARap1_2Cp->SetLineColor(kGreen+2);
  hptAccPARap2_2Cp->SetLineColor(kGreen+2);
  hrapAccCross_2Cp->SetLineColor(kGreen+2);
  hptAccCross_2Cp->SetLineColor(kGreen+2);

  TCanvas* c11 =  new TCanvas("c11","",400, 400);
  c11->cd();
  //hrapAccCross_1->Draw();
  hptAccPP_2Cp->Draw("HIST SAME");
  max_ptPP->Draw("a,p,same");
  min_ptPP->Draw("a,p,same");
  //hrapAccPA1_4->Draw("same");
  //easyLeg(leg9);
  //cleverRange(hptAccPP_2Cp);
  leg11->Draw("same");
  //if(Cat_==0){
  hptAccPP_2Cp->SetAxisRange(0.,1.05,"Y");
  hptAccPP_2Cp->SetTitleSize(0.047,"Y");
  hptAccPP_2Cp->SetTitleSize(0.047,"X");//}
//  jumSun(0,hIntAccPP_2->GetBinContent(1),30,hIntAccPP_2->GetBinContent(1));
//if(Cat_!=0){
//  hptAccPP_2Cp->SetAxisRange(0.,0.8,"Y");
//  hptAccPP_2Cp->SetTitleSize(0.047,"Y");
//  hptAccPP_2Cp->SetTitleSize(0.047,"X");}

ColMar( max_ptPP, kRed+2, 4);
ColMar( min_ptPP, kBlue+2, 3);

TCanvas* c12 =  new TCanvas("c12","",400, 400);
c12->cd();
//hrapAccCross_1->Draw();
hptAccPA_2Cp->Draw("HIST SAME");
max_ptPA->Draw("a,p,same");
min_ptPA->Draw("a,p,same");
//hptAccPA1_4->Draw("same");
//easyLeg(leg9);
//cleverRange(hptAccPA_2Cp);
leg12->Draw("same");
//if(Cat_==0){
hptAccPA_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccPA_2Cp->SetTitleSize(0.047,"Y");
hptAccPA_2Cp->SetTitleSize(0.047,"X");//}
//jumSun(0,hIntAccPA_2->GetBinContent(1),30,hIntAccPA_2->GetBinContent(1));
//if(Cat_!=0){
//hptAccPA_2Cp->SetAxisRange(0.,0.8,"Y");
//hptAccPA_2Cp->SetTitleSize(0.047,"Y");
//hptAccPA_2Cp->SetTitleSize(0.047,"X");}

ColMar( max_ptPA, kRed+2, 4);
ColMar( min_ptPA, kBlue+2, 3);

TCanvas* c13 =  new TCanvas("c13","",400, 400);
c13->cd();
//hrapAccCross_1->Draw();
hrapAccPP_2Cp->Draw("HIST SAME");
//hrapAcc23->Draw("same");
max_rapPP->Draw("a,p,same");
min_rapPP->Draw("a,p,same");
//hrapAccPA1_4->Draw("same");
//easyLeg(leg9);
//cleverRange(hrapAccPP_2Cp);
leg13->Draw("same");
//if(Cat_==0){
hrapAccPP_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPP_2Cp->SetTitleSize(0.047,"Y");
hrapAccPP_2Cp->SetTitleSize(0.047,"X");//}
jumSun(-1.93,hIntAccPP_2->GetBinContent(1),1.93,hIntAccPP_2->GetBinContent(1));
//if(Cat_!=0){
//hrapAccPP_2Cp->SetAxisRange(0.,0.8,"Y");
//hrapAccPP_2Cp->SetTitleSize(0.047,"Y");
//hrapAccPP_2Cp->SetTitleSize(0.047,"X");}

ColMar( max_rapPP, kRed+2, 4);
ColMar( min_rapPP, kBlue+2, 3);

TCanvas* c14 =  new TCanvas("c14","",400, 400);
c14->cd();
hrapAccPA_2Cp->Draw();
//hrapAcc23->Draw("a,p,same");
max_rapPA->Draw("a,p,same");
min_rapPA->Draw("a,p,same");
//hrapAccPPAN_2Cp->SetLineColor(kRed+2);
//hrapAccPPAN_2Cp->SetMarkerColor(kRed+2);
//hrapAcc23->SetLineColor(kBlue+2);
//hrapAcc23->SetMarkerColor(kBlue+2);
//ColMar(hrapAccPA_2Cp, kRed+2, 4);
//ColMar(hrapAcc23, kBlue+2, 3);
leg14->Draw("same");
hrapAccPA_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPA_2Cp->SetTitleSize(0.047,"Y");
hrapAccPA_2Cp->SetTitleSize(0.047,"X");
jumSun(-1.93,hIntAccPA_2->GetBinContent(1),1.93,hIntAccPA_2->GetBinContent(1));
ColMar( max_rapPA, kRed+2, 4);
ColMar( min_rapPA, kBlue+2, 3);

TCanvas* c17 =  new TCanvas("c17","",400, 400);
c17->cd();
hrapAccPPAN_2Cp->Draw();
hrapAcc23->Draw("a,p,same");
//hrapAccPPAN_2Cp->SetLineColor(kRed+2);
//hrapAccPPAN_2Cp->SetMarkerColor(kRed+2);
//hrapAcc23->SetLineColor(kBlue+2);
//hrapAcc23->SetMarkerColor(kBlue+2);
ColMar(hrapAccPPAN_2Cp, kRed+2, 4);
ColMar(hrapAcc23, kBlue+2, 3);
leg17->Draw("same");
hrapAccPPAN_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPPAN_2Cp->SetTitleSize(0.047,"Y");
hrapAccPPAN_2Cp->SetTitleSize(0.047,"X");

TCanvas* c1 =  new TCanvas("c1","",400, 400);
c1->cd();
//hptAccPPRap1->Draw();
hptAccPPRap1_2Cp->Draw("hist same");
max_ptPPRap1->Draw("a,p,same");
min_ptPPRap1->Draw("a,p,same");
hptAccPPRap1_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccPPRap1_2Cp->SetTitleSize(0.047,"Y");
hptAccPPRap1_2Cp->SetTitleSize(0.047,"X");
leg1->Draw("same");
cout<<"check"<<endl;

TCanvas* c2 =  new TCanvas("c2","",400, 400);
c2->cd();
//hptAccPPRap2->Draw();
hptAccPPRap2_2Cp->Draw("HIST SAME");
max_ptPPRap2->Draw("a,p,same");
min_ptPPRap2->Draw("a,p,same");
//hrapAccPP1_4->Draw("same");
hptAccPPRap2_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccPPRap2_2Cp->SetTitleSize(0.047,"Y");
hptAccPPRap2_2Cp->SetTitleSize(0.047,"X");
leg2->Draw("same");

TCanvas* c3 =  new TCanvas("c3","",400, 400);
c3->cd();
//hptAccPARap1->Draw();
hptAccPARap1_2Cp->Draw("HIST SAME");
max_ptPARap1->Draw("a,p,same");
min_ptPARap1->Draw("a,p,same");
//hptAccPA1_4->Draw("same");
hptAccPARap1_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccPARap1_2Cp->SetTitleSize(0.047,"Y");
hptAccPARap1_2Cp->SetTitleSize(0.047,"X");
leg3->Draw("same");

TCanvas* c4 =  new TCanvas("c4","",400, 400);
c4->cd();
//hptAccPARap2->Draw();
hptAccPARap2_2Cp->Draw("HIST SAME");
max_ptPARap2->Draw("a,p,same");
min_ptPARap2->Draw("a,p,same");
//hrapAccPA1_4->Draw("same");
hptAccPARap2_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccPARap2_2Cp->SetTitleSize(0.047,"Y");
hptAccPARap2_2Cp->SetTitleSize(0.047,"X");
leg4->Draw("same");

ColMar(max_ptPPRap1, kRed+2, 4);
ColMar(max_ptPPRap2, kRed+2, 4);
ColMar(max_ptPARap1, kRed+2, 4);
ColMar(max_ptPARap2, kRed+2, 4);
ColMar(min_ptPPRap1, kBlue+2, 3);
ColMar(min_ptPPRap2, kBlue+2, 3);
ColMar(min_ptPARap1, kBlue+2, 3);
ColMar(min_ptPARap2, kBlue+2, 3);

TCanvas* c5 =  new TCanvas("c5","",400, 400);
c5->cd();
//hrapAccPPPt1->Draw();
hrapAccPPPt1_2Cp->Draw("hist same");
max_rapPPPt1->Draw("a,p,same");
min_rapPPPt1->Draw("a,p,same");
//if(Cat_==0){
hrapAccPPPt1_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPPPt1_2Cp->SetTitleSize(0.047,"Y");
hrapAccPPPt1_2Cp->SetTitleSize(0.047,"X");//}
jumSun(-1.93,hIntAccPPPt1_2->GetBinContent(1),1.93,hIntAccPPPt1_2->GetBinContent(1));
//if(Cat_!=0){
//hrapAccPPPt1_2Cp->SetAxisRange(0.,0.65,"Y");
//hrapAccPPPt1_2Cp->SetTitleSize(0.047,"Y");
//hrapAccPPPt1_2Cp->SetTitleSize(0.047,"X");}
leg5->Draw("same");
cout<<"check"<<endl;
c5->Update();

TCanvas* c6 =  new TCanvas("c6","",400, 400);
c6->cd();
//hrapAccPPPt2->Draw();
hrapAccPPPt2_2Cp->Draw("HIST SAME");
max_rapPPPt2->Draw("a,p,same");
min_rapPPPt2->Draw("a,p,same");
//hrapAccPP1_4->Draw("same");
hrapAccPPPt2_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPPPt2_2Cp->SetTitleSize(0.047,"Y");
hrapAccPPPt2_2Cp->SetTitleSize(0.047,"X");
jumSun(-1.93,hIntAccPPPt2_2->GetBinContent(1),1.93,hIntAccPPPt2_2->GetBinContent(1));
leg6->Draw("same");

TCanvas* c7 =  new TCanvas("c7","",400, 400);
c7->cd();
//hrapAccPAPt1->Draw();
hrapAccPAPt1_2Cp->Draw("HIST SAME");
max_rapPAPt1->Draw("a,p,same");
min_rapPAPt1->Draw("a,p,same");
//hrapAccPA1_4->Draw("same");
//if(Cat_==0){
hrapAccPAPt1_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPAPt1_2Cp->SetTitleSize(0.047,"Y");
hrapAccPAPt1_2Cp->SetTitleSize(0.047,"X");//}
jumSun(-1.93,hIntAccPAPt1_2->GetBinContent(1),1.93,hIntAccPAPt1_2->GetBinContent(1));
//if(Cat_!=0){
//hrapAccPAPt1_2Cp->SetAxisRange(0.,0.65,"Y");
//hrapAccPAPt1_2Cp->SetTitleSize(0.047,"Y");
//hrapAccPAPt1_2Cp->SetTitleSize(0.047,"X");}
leg7->Draw("same");

TCanvas* c8 =  new TCanvas("c8","",400, 400);
c8->cd();
//hrapAccPAPt2->Draw();
hrapAccPAPt2_2Cp->Draw("HIST SAME");
max_rapPAPt2->Draw("a,p,same");
min_rapPAPt2->Draw("a,p,same");
//hrapAccPA1_4->Draw("same");
hrapAccPAPt2_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccPAPt2_2Cp->SetTitleSize(0.047,"Y");
hrapAccPAPt2_2Cp->SetTitleSize(0.047,"X");
jumSun(-1.93,hIntAccPAPt2_2->GetBinContent(1),1.93,hIntAccPAPt2_2->GetBinContent(1));
leg8->Draw("same");

ColMar(max_rapPPPt1, kRed+2, 4);
ColMar(max_rapPPPt2, kRed+2, 4);
ColMar(max_rapPAPt1, kRed+2, 4);
ColMar(max_rapPAPt2, kRed+2, 4);
ColMar(min_rapPPPt1, kBlue+2, 3);
ColMar(min_rapPPPt2, kBlue+2, 3);
ColMar(min_rapPAPt1, kBlue+2, 3);
ColMar(min_rapPAPt2, kBlue+2, 3);

TCanvas* c9 =  new TCanvas("c9","",400, 400);
c9->cd();
//hrapAccCross_1->Draw();
hrapAccCross_2Cp->Draw("HIST SAME");
max_rapCross->Draw("a,p,same");
min_rapCross->Draw("a,p,same");
//hrapAccPA1_4->Draw("same");
//easyLeg(leg9);
//cleverRange(hrapAccCross_2Cp);
leg9->Draw("same");
//if(Cat_==0){
hrapAccCross_2Cp->SetAxisRange(0.,1.05,"Y");
hrapAccCross_2Cp->SetTitleSize(0.047,"Y");
hrapAccCross_2Cp->SetTitleSize(0.047,"X");//}
jumSun(-2.87,hIntAccCross_2->GetBinContent(1),1.93,hIntAccCross_2->GetBinContent(1));
//if(Cat_!=0){
//hrapAccCross_2Cp->SetAxisRange(0.,0.65,"Y");
//hrapAccCross_2Cp->SetTitleSize(0.047,"Y");
//hrapAccCross_2Cp->SetTitleSize(0.047,"X");}

TCanvas* c10 =  new TCanvas("c10","",400, 400);
c10->cd();
//hrapAccCross_1->Draw();
hptAccCross_2Cp->Draw("HIST SAME");
max_ptCross->Draw("a,p,same");
min_ptCross->Draw("a,p,same");
//hptAccPA1_4->Draw("same");
//easyLeg(leg9);
//cleverRange(hptAccCross_2Cp);
leg10->Draw("same");
//if(Cat_==0){
hptAccCross_2Cp->SetAxisRange(0.,1.05,"Y");
hptAccCross_2Cp->SetTitleSize(0.047,"Y");
hptAccCross_2Cp->SetTitleSize(0.047,"X");//}
//jumSun(0,hIntAccCross_2->GetBinContent(1),30,hIntAccCross_2->GetBinContent(1));
//if(Cat_!=0){
//hptAccCross_2Cp->SetAxisRange(0.,0.8,"Y");
//hptAccCross_2Cp->SetTitleSize(0.047,"Y");
//hptAccCross_2Cp->SetTitleSize(0.047,"X");}

ColMar(max_rapCross, kRed+2, 4);
ColMar(min_rapCross, kBlue+2, 3);
ColMar(max_ptCross, kRed+2, 4);
ColMar(min_ptCross, kBlue+2, 3);

c11->SaveAs(Form("%d/acceptance_pt_PP_%dS_%d.pdf", date, Cat_+1,date));
c12->SaveAs(Form("%d/acceptance_pt_PA_%dS_%d.pdf",date, Cat_+1,date));
c13->SaveAs(Form("%d/acceptance_rap_PP_%dS_%d.pdf", date, Cat_+1,date));
c14->SaveAs(Form("%d/acceptance_rap_PA_%dS_%d.pdf", date, Cat_+1,date));
c17->SaveAs(Form("%d/acceptance_rap_PP_vs_23_%dS_%d.pdf", date, Cat_+1,date));
c1->SaveAs(Form("%d/acceptance_pt_PP_Rap1_%dS_%d.pdf", date, Cat_+1,date));
c2->SaveAs(Form("%d/acceptance_pt_PP_Rap2_%dS_%d.pdf", date, Cat_+1,date));
c3->SaveAs(Form("%d/acceptance_pt_PA_Rap1_%dS_%d.pdf",date, Cat_+1,date));
c4->SaveAs(Form("%d/acceptance_pt_PA_Rap2_%dS_%d.pdf",date, Cat_+1,date));
c5->SaveAs(Form("%d/acceptance_rap_PP_Pt1_%dS_%d.pdf", date, Cat_+1,date));
c6->SaveAs(Form("%d/acceptance_rap_PP_Pt2_%dS_%d.pdf", date, Cat_+1,date));
c7->SaveAs(Form("%d/acceptance_rap_PA_pt1_%dS_%d.pdf",date, Cat_+1,date));
c8->SaveAs(Form("%d/acceptance_rap_PA_pt2_%dS_%d.pdf",date, Cat_+1,date));
c9->SaveAs(Form("%d/acceptance_rap_PA_Cross_%dS_%d.pdf",date, Cat_+1,date));
c10->SaveAs(Form("%d/acceptance_pt_PA_Cross_%dS_%d.pdf",date, Cat_+1,date));

out->Write();

if(Cat_==1){
TH1D * h2DPP = (TH1D*)hIntAccPP2D_2->Clone();
TH1D * h2DPA = (TH1D*)hIntAccPA2D_2->Clone();
TCanvas* c15 =  new TCanvas("c15","",1200, 600);
c15->Divide(3,2);
c15->cd(1);
hIntAccPP2D_1->Draw("colz");
hIntAccPP2D_1->SetTitleSize(0.048,"Y");
hIntAccPP2D_1->SetTitleSize(0.048,"X");
c15->cd(2);
hIntAccPP2D_2->Draw("colz");
hIntAccPP2D_2->SetTitleSize(0.048,"Y");
hIntAccPP2D_2->SetTitleSize(0.048,"X");
c15->cd(3);
h2DPP->Divide(hIntAccPP2D_1);
h2DPP->Draw("colz");
h2DPP->SetTitleSize(0.048,"Y");
h2DPP->SetTitleSize(0.048,"X");
c15->cd(4);
hIntAccPA2D_1->Draw("colz");
hIntAccPA2D_1->SetTitleSize(0.048,"Y");
hIntAccPA2D_1->SetTitleSize(0.048,"X");
c15->cd(5);
hIntAccPA2D_2->Draw("colz");
hIntAccPA2D_2->SetTitleSize(0.048,"Y");
hIntAccPA2D_2->SetTitleSize(0.048,"X");
c15->cd(6);
h2DPA->Divide(hIntAccPA2D_1);
h2DPA->Draw("colz");
h2DPA->SetTitleSize(0.048,"Y");
h2DPA->SetTitleSize(0.048,"X");
c15->SaveAs(Form("%d/acceptance_2D_%dS_%d.pdf",date,Cat_+1,date));
}
//TCanvas* c15 =  new TCanvas("c15","",600, 600);
//c15->cd();
//hIntAccPA2D_2->Divide(hIntAccPA2D_1);
//hIntAccPA2D_2->Draw("colz");
//c15->SaveAs(Form("acceptance_2D_PA_%dS_%d.pdf",Cat_+1,date));

if(Cat_==0){
cout<<"Integrated bin pPb"<<endl;
cout<<"No weighted Acc : "<<hIntAccPANoW_2->GetBinContent(1)<<endl;
cout<<"Weighted Acc : "<<hIntAccPA_2->GetBinContent(1)<<endl;
cout<<""<<endl;
}

if(Cat_==0){
for(int i=1; i<=hrapAccCross_2Cp->GetNbinsX(); i++){
  cout<<"Acc Cross "<<i<<": "<<hrapAccCross_2Cp->GetBinContent(i)<<" , NoW : "<<hrapAccCrossNoW_2->GetBinContent(i)<<endl;
}
TCanvas* c16 =  new TCanvas("c16","",600, 600);
TH1D * hrap13 = (TH1D*)hIntAccPA_2->Clone();
hrap13->SetBinContent(1, 0.210);
hrapAccCrossNoW_2->Draw();
hrapAccCross_2Cp->Draw("same");
hrap13->Draw("same");
hIntAccPA_2->Draw("same");
hIntAccPANoW_2->Draw("same");
hrapAccPPAN_2Cp->Draw("same");
hraptest2->Draw("same");
hrapAccCrossNoW_2->SetAxisRange(0.,0.8,"Y");
ColMar(hrap13, kRed+2, 23);
ColMar(hrapAccPPAN_2Cp, kRed+2, 26);
ColMar(hrapAccCross_2Cp, kGreen+2, 20);
ColMar(hrapAccCrossNoW_2, kBlue+2, 25);
ColMar(hIntAccPA_2, kGreen+2, 20);
ColMar(hIntAccPANoW_2, kBlue+2, 25);
}

}

double getMaxDev(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4)
{
  double MaxDev;
  int imax=h1->GetNbinsX();

  for ( int i=0 ; i<=imax+1; i++)
  {
    if( fabs(h1->GetBinContent(i)-h2->GetBinContent(i)) > fabs(h1->GetBinContent(i)-h3->GetBinContent(i)) ){
      MaxDev = h2->GetBinContent(i);
    }
    else MaxDev = h3->GetBinContent(i);
    //  cout<<"TEST"<<i<<" : "<<MaxDev<<endl;
    h4->SetBinContent(i,MaxDev);
  }
  return MaxDev;
}

void ColMar(TH1D* h1, Color_t kBlack, int a=0)
{
  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(a);
  return ColMar;
}
