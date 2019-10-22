#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
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
double Acc_NoW[10][10] = {0.0};
double Acc_NoW_Err[10][10] = {0.0};

double Acc_pp[10][10] = {0.0};
double Acc_pp_Err[10][10] = {0.0};

double Acc_AA[10][10] = {0.0};
double Acc_AA_Err[10][10] = {0.0};

TH1D *hAccIntNoW;
TH1D *hAccIntpp;
TH1D *hAccIntAA;

TH1D *hAccPtNoW;
TH1D *hAccPtpp;
TH1D *hAccPtAA;

TH1D *hAccYNoW;
TH1D *hAccYpp;
TH1D *hAccYAA;

// New weighting function 20161207
TF1* fWgtPP1 = new TF1("fWgtPP1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtPP2 = new TF1("fWgtPP2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");
TF1* fWgtAA1 = new TF1("fWgtAA1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtAA2 = new TF1("fWgtAA2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");

TH1D *hPadPt = new TH1D("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,30.0);
TH1D *hPadY = new TH1D("hPadY",";|y|;Acceptance",10,0.0,2.4);

TLatex *lt1 = new TLatex();

float getAcceptanceSingleBin(int state=k1S, int icoll=kNoWeight, TString kineCut="", int varianceID=0);
void getAcceptance(int doDraw = 0, int state = k1S) { // doDraw == 1 : draw, == 0 : no draw
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

  // Integrated bin : 
  const int nPtBinsInt = 1;   double ptBinInt[nPtBinsInt+1] = {0,30};  const int nYBinsInt = 1;   double yBinInt[nYBinsInt+1] = {0,2.4}; // integrated
  
  // A : pT dependence 
  int nPtBinsA;    
  double* ptBinA; 
  int nYBinsA;    
  double* yBinA; 
  // B : y dependence 
  int nPtBinsB;    
  double* ptBinB; 
  int nYBinsB;    
  double* yBinB; 
  
  if ( state == 1 ) {
    // pt depdence bin
    nPtBinsA = nPtBins1s;    ptBinA = ptBin1s;
    nYBinsA = 1;   yBinA = yBinInt;
    // y dependence bin 
    nPtBinsB = 1;  ptBinB = ptBinInt ; 
    nYBinsB =  nYBins1S; yBinB = yBin1S;
  }
  else if ( state == 2 ) {
    nPtBinsA = nPtBins2s;    ptBinA = ptBin2s;
    nYBinsA = 1;   yBinA = yBinInt;

    nPtBinsB = 1;   ptBinB = ptBinInt ; 
    nYBinsB =  nYBins2S; yBinB = yBin2S;
  }
  else if ( state == 3 ) {
    nPtBinsA = nPtBins3s;    ptBinA = ptBin3s;
    nYBinsA = 1;   yBinA = yBinInt;

    nPtBinsB = 1;   ptBinB = ptBinInt ; 
    nYBinsB =  nYBins3S; yBinB = yBin3S;
  }
  TFile *out = new TFile(Form("acceptance_wgt_final_XXXXXXX_%dS.root",state),"RECREATE");
  //TFile *out = new TFile(Form("acceptance_bin%d_wgt.root",bin),"RECREATE");

  
  // Integrated bin :
  hAccIntNoW = new TH1D("hAccIntNoW",";p_{T} (GeV/c);Acceptance of #Upsilon()",nPtBinsInt,ptBinInt);
  hAccIntpp = new TH1D("hAccIntpp",";p_{T} (GeV/c);Acceptance of #Upsilon()",nPtBinsInt,ptBinInt);
  hAccIntAA = new TH1D("hAccIntAA",";p_{T} (GeV/c);Acceptance of #Upsilon()",nPtBinsInt,ptBinInt);
  hAccIntNoW->Sumw2();
  hAccIntpp->Sumw2();
  hAccIntAA->Sumw2();
  hAccIntNoW->SetMarkerStyle(20);hAccIntNoW->SetMarkerColor(kBlue+2);
  hAccIntpp->SetMarkerStyle(20);hAccIntpp->SetMarkerColor(kBlue+2);
  hAccIntAA->SetMarkerStyle(20);hAccIntAA->SetMarkerColor(kBlue+2);
  

  // pT dependence : Plan A
  hAccPtNoW = new TH1D("hAccPtNoW",";p_{T} (GeV/c);Acceptance of #Upsilon",nPtBinsA,ptBinA);
  hAccPtpp = new TH1D("hAccPtpp",";p_{T} (GeV/c);Acceptance of #Upsilon",nPtBinsA,ptBinA);
  hAccPtAA = new TH1D("hAccPtAA",";p_{T} (GeV/c);Acceptance of #Upsilon",nPtBinsA,ptBinA);
  hAccPtNoW->Sumw2();
  hAccPtpp->Sumw2();
  hAccPtAA->Sumw2();
  hAccPtNoW->SetMarkerStyle(20);hAccPtNoW->SetMarkerColor(kBlue+2);
  hAccPtpp->SetMarkerStyle(20);hAccPtpp->SetMarkerColor(kBlue+2);
  hAccPtAA->SetMarkerStyle(20);hAccPtAA->SetMarkerColor(kBlue+2);
  
  // y dependdence : Plan B
  hAccYNoW = new TH1D("hAccYNoW",";|y|;Acceptance of #Upsilon",nYBinsB,yBinB);
  hAccYpp = new TH1D("hAccYpp",";|y|;Acceptance of #Upsilon",nYBinsB,yBinB);
  hAccYAA = new TH1D("hAccYAA",";|y|;Acceptance of #Upsilon",nYBinsB,yBinB);
  hAccYNoW->Sumw2();
  hAccYpp->Sumw2();
  hAccYAA->Sumw2();
  hAccYNoW->SetMarkerStyle(20);hAccYNoW->SetMarkerColor(kBlue+2);
  hAccYpp->SetMarkerStyle(20);hAccYpp->SetMarkerColor(kBlue+2);
  hAccYAA->SetMarkerStyle(20);hAccYAA->SetMarkerColor(kBlue+2);
  
    
  // Integrated Bin
  
  for ( int ipt = 1 ; ipt<=nPtBinsInt ; ipt++ ) { 
    for ( int iy =1 ; iy<=nYBinsInt ; iy++) {
      TString kRange = Form("pt>%0.1f && pt<%0.1f && abs(y)>%0.1f && abs(y)<%0.1f",(float)ptBinInt[ipt-1],(float)ptBinInt[ipt], (float)yBinInt[iy-1], (float)yBinInt[iy] ) ;
      cout << "*===*===*===*===*  " << kRange << "  *===*===*===*===*" << endl;
      float y_noRwt = getAcceptanceSingleBin( 1, kNoWeight, kRange,kNoVar);
      float y_ppRwt = getAcceptanceSingleBin( 1, kPP, kRange,kNoVar);
      float y_aaRwt = getAcceptanceSingleBin( 1, kAA, kRange,kNoVar);

      float err_ppPt = 0 ;  // tentatively stopped this function
      float err_ppY = 0; 
      float err_aaPt = 0;
      float err_aaY = 0;

      //      float err_ppPt = max ( abs(y_ppRwt_ptPlus/y_ppRwt-1), abs(y_ppRwt_ptMinus/y_ppRwt-1) ) ;
      //      float err_ppY = max ( abs(y_ppRwt_yPlus/y_ppRwt-1), abs(y_ppRwt_yMinus/y_ppRwt-1) ) ;
      //      float err_aaPt = max ( abs(y_aaRwt_ptPlus/y_aaRwt-1), abs(y_aaRwt_ptMinus/y_aaRwt-1) ) ;
      //      float err_aaY = max ( abs(y_aaRwt_yPlus/y_aaRwt-1), abs(y_aaRwt_yMinus/y_aaRwt-1) ) ;
      
      // dmoon added systematcis : unweighted - weighted
      float err_pp = abs(y_noRwt - y_ppRwt)/y_ppRwt;
      float err_aa = abs(y_noRwt - y_aaRwt)/y_aaRwt;    
      
      Acc_NoW[ipt][iy] = y_noRwt;
      Acc_NoW_Err[ipt][iy] = err_ppPt;
      
      Acc_pp[ipt][iy] = y_ppRwt;
      Acc_pp_Err[ipt][iy] = err_pp;
      
      Acc_AA[ipt][iy] = y_aaRwt;
      Acc_AA_Err[ipt][iy] = err_aa;    

      
      cout << "Un-weighed acceptance :  = " << y_noRwt << endl;
      cout << "Re-weighed Acceptance : " << "pp = " << y_ppRwt << ",  PbPb = " << y_aaRwt << endl;
      
      cout << int(1000*y_noRwt)*0.001 << " & " << int(1000*y_ppRwt)*0.001<< " & " << int(1000*y_aaRwt)*0.001<< " & " ;
      cout <<""<<endl;
      
      cout << int(1000*err_ppPt)*0.1 << "\% & " << int(1000*err_aaPt)*0.1 << "\% & " ;
      cout << int(1000*err_ppY)*0.1 << "\% & " << int(1000*err_aaY)*0.1 << "\% & " ;
      cout <<""<<endl;
      
      cout<<"##### Systematics #####"<<endl;
      cout<<"pp "<<  "    AA "<<  endl;
      cout<< int(1000*err_pp)*0.1 << "\% & " << int(1000*err_aa)*0.1 << "\% & ";
      //cout << int(1000*sqrt( err_ppPt*err_ppPt +  err2s_ppPt* err2s_ppPt + err_ppY*err_ppY + err2s_ppY*err2s_ppY + err_aaPt*err_aaPt +  err2s_aaPt* err2s_aaPt + err_aaY*err_aaY + err2s_aaY*err2s_aaY ) ) *0.1 << "\%" << endl;
      cout<<""<<endl;
      
      hAccIntNoW->SetBinContent(ipt,Acc_NoW[ipt][iy]);
      hAccIntpp->SetBinContent(ipt,Acc_pp[ipt][iy]);
      hAccIntAA->SetBinContent(ipt,Acc_AA[ipt][iy]);
      
      hAccIntNoW->SetBinError(ipt,Acc_NoW_Err[ipt][iy]);
      hAccIntpp->SetBinError(ipt,Acc_pp_Err[ipt][iy]);
      hAccIntAA->SetBinError(ipt,Acc_AA_Err[ipt][iy]);
    } // iy
  }    // ipt

    // pt dependence 
    for ( int ipt = 1 ; ipt<=nPtBinsA ; ipt++ ) { 
      for ( int iy =1 ; iy<=nYBinsA ; iy++) {
	TString kRange = Form("pt>%0.1f && pt<%0.1f && abs(y)>%0.1f && abs(y)<%0.1f",(float)ptBinA[ipt-1],(float)ptBinA[ipt], (float)yBinA[iy-1], (float)yBinA[iy] ) ;
	cout << "*===*===*===*===*  " << kRange << "  *===*===*===*===*" << endl;
	float y_noRwt = getAcceptanceSingleBin( 1, kNoWeight, kRange,kNoVar);
	float y_ppRwt = getAcceptanceSingleBin( 1, kPP, kRange,kNoVar);
	float y_aaRwt = getAcceptanceSingleBin( 1, kAA, kRange,kNoVar);
	
	float err_ppPt = 0 ;  // tentatively stopped this function
	float err_ppY = 0; 
	float err_aaPt = 0;
	float err_aaY = 0;
	//	float err_ppPt = max ( abs(y_ppRwt_ptPlus/y_ppRwt-1), abs(y_ppRwt_ptMinus/y_ppRwt-1) ) ;
	//	float err_ppY = max ( abs(y_ppRwt_yPlus/y_ppRwt-1), abs(y_ppRwt_yMinus/y_ppRwt-1) ) ;
	//      float err_aaPt = max ( abs(y_aaRwt_ptPlus/y_aaRwt-1), abs(y_aaRwt_ptMinus/y_aaRwt-1) ) ;
	//	float err_aaY = max ( abs(y_aaRwt_yPlus/y_aaRwt-1), abs(y_aaRwt_yMinus/y_aaRwt-1) ) ;
	
	// dmoon added systematcis : unweighted - weighted
	float err_pp = abs(y_noRwt - y_ppRwt)/y_ppRwt;
	float err_aa = abs(y_noRwt - y_aaRwt)/y_aaRwt;
	
	Acc_NoW[ipt][iy] = y_noRwt;
	Acc_NoW_Err[ipt][iy] = err_ppPt;
	
	Acc_pp[ipt][iy] = y_ppRwt;
	Acc_pp_Err[ipt][iy] = err_pp;
	
	Acc_AA[ipt][iy] = y_aaRwt;
	Acc_AA_Err[ipt][iy] = err_aa;
	
	cout << "Un-weighed acceptance :  = " << y_noRwt << endl;
	cout << "Re-weighed Acceptance : " << "pp = " << y_ppRwt << ",  PbPb = " << y_aaRwt << endl;
      
	cout << int(1000*y_noRwt)*0.001 << " & " << int(1000*y_ppRwt)*0.001<< " & " << int(1000*y_aaRwt)*0.001<< " & " ;
	cout <<""<<endl;
	
	cout << int(1000*err_ppPt)*0.1 << "\% & " << int(1000*err_aaPt)*0.1 << "\% & " ;
	cout << int(1000*err_ppY)*0.1 << "\% & " << int(1000*err_aaY)*0.1 << "\% & " ;
	cout <<""<<endl;
      
	cout<<"##### Systematics #####"<<endl;
	cout<<"pp "<<  "    AA "<<  endl;
	cout<< int(1000*err_pp)*0.1 << "\% & " << int(1000*err_aa)*0.1 << "\% & ";
	//cout << int(1000*sqrt( err_ppPt*err_ppPt +  err2s_ppPt* err2s_ppPt + err_ppY*err_ppY + err2s_ppY*err2s_ppY + err_aaPt*err_aaPt +  err2s_aaPt* err2s_aaPt + err_aaY*err_aaY + err2s_aaY*err2s_aaY ) ) *0.1 << "\%" << endl;
	cout<<""<<endl;
	
	hAccPtNoW->SetBinContent(ipt,Acc_NoW[ipt][iy]);
	hAccPtpp->SetBinContent(ipt,Acc_pp[ipt][iy]);
	hAccPtAA->SetBinContent(ipt,Acc_AA[ipt][iy]);
	hAccPtNoW->SetBinError(ipt,Acc_NoW_Err[ipt][iy]);
	hAccPtpp->SetBinError(ipt,Acc_pp_Err[ipt][iy]);
	hAccPtAA->SetBinError(ipt,Acc_AA_Err[ipt][iy]);
      } // iy
    }   // ipt
    
    
    // y dependence 
    for ( int ipt = 1 ; ipt<=nPtBinsB ; ipt++ ) { 
      for ( int iy =1 ; iy<=nYBinsB ; iy++) {
	TString kRange = Form("pt>%0.1f && pt<%0.1f && abs(y)>%0.1f && abs(y)<%0.1f",(float)ptBinB[ipt-1],(float)ptBinB[ipt], (float)yBinB[iy-1], (float)yBinB[iy] ) ;
	cout << "*===*===*===*===*  " << kRange << "  *===*===*===*===*" << endl;
	float y_noRwt = getAcceptanceSingleBin( 1, kNoWeight, kRange,kNoVar);
	float y_ppRwt = getAcceptanceSingleBin( 1, kPP, kRange,kNoVar);
	float y_aaRwt = getAcceptanceSingleBin( 1, kAA, kRange,kNoVar);
	
	float err_ppPt = 0;
	float err_ppY = 0;
	float err_aaPt =0;
	float err_aaY = 0;

	//	float err_ppPt = max ( abs(y_ppRwt_ptPlus/y_ppRwt-1), abs(y_ppRwt_ptMinus/y_ppRwt-1) ) ;
	//	float err_ppY = max ( abs(y_ppRwt_yPlus/y_ppRwt-1), abs(y_ppRwt_yMinus/y_ppRwt-1) ) ;
	//	float err_aaPt = max ( abs(y_aaRwt_ptPlus/y_aaRwt-1), abs(y_aaRwt_ptMinus/y_aaRwt-1) ) ;
	//	float err_aaY = max ( abs(y_aaRwt_yPlus/y_aaRwt-1), abs(y_aaRwt_yMinus/y_aaRwt-1) ) ;
	
	// dmoon added systematcis : unweighted - weighted
	float err_pp = abs(y_noRwt - y_ppRwt)/y_ppRwt;
	float err_aa = abs(y_noRwt - y_aaRwt)/y_aaRwt;
	
	Acc_NoW[ipt][iy] = y_noRwt;
	Acc_NoW_Err[ipt][iy] = err_ppPt;
	
	Acc_pp[ipt][iy] = y_ppRwt;
	Acc_pp_Err[ipt][iy] = err_pp;
	
	Acc_AA[ipt][iy] = y_aaRwt;
	Acc_AA_Err[ipt][iy] = err_aa;
	
	cout << "Un-weighed acceptance :  = " << y_noRwt << endl;
	cout << "Re-weighed Acceptance : " << "pp = " << y_ppRwt << ",  PbPb = " << y_aaRwt << endl;
      
	cout << int(1000*y_noRwt)*0.001 << " & " << int(1000*y_ppRwt)*0.001<< " & " << int(1000*y_aaRwt)*0.001<< " & " ;
	cout <<""<<endl;
	
	cout << int(1000*err_ppPt)*0.1 << "\% & " << int(1000*err_aaPt)*0.1 << "\% & " ;
	cout << int(1000*err_ppY)*0.1 << "\% & " << int(1000*err_aaY)*0.1 << "\% & " ;
	cout <<""<<endl;
      
	cout<<"##### Systematics #####"<<endl;
	cout<<"pp "<<  "    AA "<<  endl;
	cout<< int(1000*err_pp)*0.1 << "\% & " << int(1000*err_aa)*0.1 << "\% & ";
	//cout << int(1000*sqrt( err_ppPt*err_ppPt +  err2s_ppPt* err2s_ppPt + err_ppY*err_ppY + err2s_ppY*err2s_ppY + err_aaPt*err_aaPt +  err2s_aaPt* err2s_aaPt + err_aaY*err_aaY + err2s_aaY*err2s_aaY ) ) *0.1 << "\%" << endl;
	cout<<""<<endl;
	
	hAccYNoW->SetBinContent(iy,Acc_NoW[ipt][iy]);
	hAccYpp->SetBinContent(iy,Acc_pp[ipt][iy]);
	hAccYAA->SetBinContent(iy,Acc_AA[ipt][iy]);
	hAccYNoW->SetBinError(iy,Acc_NoW_Err[ipt][iy]);
	hAccYpp->SetBinError(iy,Acc_pp_Err[ipt][iy]);
	hAccYAA->SetBinError(iy,Acc_AA_Err[ipt][iy]);
      } // iy
    } // ipt
  
    
    out->cd();
    
    hAccIntNoW->SetName(Form("hIntAccNoW%dS",state));
    hAccIntpp->SetName(Form("hIntAccPP%dS",state));
    hAccIntAA->SetName(Form("hIntAccAA%dS",state));
    hAccIntNoW->Write();
    hAccIntpp->Write();
    hAccIntAA->Write();
    
    hAccPtNoW->SetName(Form("hptAccNoW%dS",state));
    hAccPtpp->SetName(Form("hptAccPP%dS",state));
    hAccPtAA->SetName(Form("hptAccAA%dS",state));
    hAccPtNoW->Write();
    hAccPtpp->Write();
    hAccPtAA->Write();
    hAccYNoW->SetName(Form("hrapAccNoW%dS",state));
    hAccYpp->SetName(Form("hrapAccPP%dS",state));
    hAccYAA->SetName(Form("hrapAccAA%dS",state));
    hAccYNoW->Write();
    hAccYpp->Write();
    hAccYAA->Write();
    
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


