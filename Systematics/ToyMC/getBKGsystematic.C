#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "../../commonUtility.h"
#include "../../cutsAndBin.h"

using namespace std;

void getBKGsystematic(int state)
{
  
  gStyle->SetOptStat(0);

  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double *yBin;

  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
    nCentBins = nCentBins1s;  centBin = centBin1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
    nYBins = nYBins2S;  yBin = yBin2S;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
    nYBins = nYBins3S;  yBin = yBin3S;
  }
 
  TFile *file_AA_cent[nCentBins];

  TFile *file_PP_Pt[nPtBins];
  TFile *file_AA_Pt[nPtBins];

  TFile *file_PP_Rap[nYBins];
  TFile *file_AA_Rap[nYBins];

  TFile *file_PP_Int;
  TFile *file_AA_Int;
 
  for(icent=0;icent<nCentBins;icent++)
  {
    file_AA_cent[icent] = new TFile(Form("addFiles/AA_fit_pt0.0-30.0_rap_0.0-2.4_cent%d-%d_Gen1000000_input3_useCentBkg1_nToys1_all.root",centBin[icent],centBin[icent+1]),"read");
  }

  for(ipt=0;ipt<nPtBins;ipt++)
  {
    file_PP_Pt[ipt] = new TFile(Form("addFiles/PP_fit_pt%.1f-%.1f_rap0.0-2.4_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root",ptBin[ipt],ptBin[ipt+1]),"read");
    file_AA_Pt[ipt] = new TFile(Form("addFiles/AA_fit_pt%.1f-%.1f_rap0.0-2.4_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root",ptBin[ipt],ptBin[ipt+1]),"read");
  }

  for(irap=0;irap<nYBins;irap++)
  {
    file_PP_Rap[irap] = new TFile(Form("addFiles/PP_fit_pt0.0-30.0_rap%.1f-%.1f_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root",yBin[irap],yBin[irap+1]),"read");
    file_AA_Rap[irap] = new TFile(Form("addFiles/AA_fit_pt0.0-30.0_rap%.1f-%.1f_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root",yBin[irap],yBin[irap+1]),"read");
  }

  file_PP_Int = new TFile("addFiles/PP_fit_pt0.0-30.0_rap0.0-2.4_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root","read");
  file_AA_Int = new TFile("addFiles/AA_fit_pt0.0-30.0_rap0.0-2.4_cent0-200_Gen1000000_input3_useCentBkg1_nToys1_all.root","read");

  
  TCanvas *c_PP_Pt[3];
  TCanvas *c_AA_Pt[3];

  TCanvas *c_PP_Rap[3];
  TCanvas *c_AA_Rap[3];

  TCanvas *c_AA_Cent[3];

  TCanvas *c_PP_int[3];
  TCanvas *c_AA_int[3];

  int canv_pt_dev1, canv_pt_dev2;
  int canv_rap_dev1, canv_rap_dev2;
  int canv_cent_dev1, canv_cent_dev2;

  if(state==1) {canv_pt_dev1 = 2; canv_pt_dev2 = 3; canv_rap_dev1 = 2; canv_rap_dev2 = 3; canv_cent_dev1 = 3; canv_cent_dev2 = 3;}
  else if(state==2) {canv_pt_dev1 = 1; canv_pt_dev2 = 3; canv_rap_dev1 = 1; canv_rap_dev2 = 2; canv_cent_dev1 = 3; canv_cent_dev2 = 3;}
  else if(state==3) {canv_pt_dev1 = 1; canv_pt_dev2 = 3; canv_rap_dev1 = 1; canv_rap_dev2 = 2; canv_cent_dev1 = 2; canv_cent_dev2 = 2;}



  for(int i=0;i<3;i++)
  {
    c_PP_Pt[i] = new TCanvas(Form("c_PP_Pt_%d",i),"",700,700);
    c_AA_Pt[i] = new TCanvas(Form("c_PP_Pt_%d",i),"",700,700);
    c_PP_Rap[i] = new TCanvas(Form("c_PP_Rap%d",i),"",700,700);
    c_AA_Rap[i] = new TCanvas(Form("c_AA_Rap%d",i),"",700,700);
    c_PP_int[i] = new TCanvas(Form("c_PP_int%d",i),"",700,700);
    c_AA_int[i] = new TCanvas(Form("c_AA_int%d",i),"",700,700);
    c_AA_Cent[i] = new TCanvas(Form("c_AA_Cent%d",i),"",700,700);

    c_PP_Pt[i] -> Divide(canv_pt_dev1,canv_pt_dev2);
    c_AA_Pt[i] -> Divide(canv_pt_dev1,canv_pt_dev2);
    c_PP_Rap[i] -> Divide(canv_rap_dev1,canv_rap_dev2);
    c_AA_Rap[i] -> Divide(canv_rap_dev1,canv_rap_dev2);
    c_AA_Cent[i] -> Divide(canv_cent_dev1,canv_cent_dev2);
  }

  TH1D *hPP_Pt_Nom[nPtBins];
  TH1D *hPP_Pt_Alt[nPtBins];
  TH1D *hPP_Pt_Dev[nPtBins];

  TH1D *hAA_Pt_Nom[nPtBins];
  TH1D *hAA_Pt_Alt[nPtBins];
  TH1D *hAA_Pt_Dev[nPtBins];
  
  TH1D *hPP_Rap_Nom[nYBins];
  TH1D *hPP_Rap_Alt[nYBins];
  TH1D *hPP_Rap_Dev[nYBins];
  
  TH1D *hAA_Rap_Nom[nYBins];
  TH1D *hAA_Rap_Alt[nYBins];
  TH1D *hAA_Rap_Dev[nYBins];

  TH1D *hAA_cent_Nom[nCentBins];
  TH1D *hAA_cent_Alt[nCentBins];
  TH1D *hAA_cent_Dev[nCentBins];

  TH1D *hPP_int_Nom;
  TH1D *hPP_int_Alt;
  TH1D *hPP_int_Dev;
  TH1D *hAA_int_Nom;
  TH1D *hAA_int_Alt;
  TH1D *hAA_int_Dev; 

  TF1 *f1 = new TF1("f1","gaus",0,100);

  TString perc="%";

  for(int ipt=0;ipt<nPtBins;ipt++)
  {
    hPP_Pt_Nom[ipt] = (TH1D*) file_PP_Pt[ipt]->Get(Form("hNom_%dS",state));
    hPP_Pt_Alt[ipt] = (TH1D*) file_PP_Pt[ipt]->Get(Form("hAlt_%dS",state));
    hPP_Pt_Dev[ipt] = (TH1D*) file_PP_Pt[ipt]->Get(Form("hDev_%dS",state));
    hAA_Pt_Nom[ipt] = (TH1D*) file_AA_Pt[ipt]->Get(Form("hNom_%dS",state));
    hAA_Pt_Alt[ipt] = (TH1D*) file_AA_Pt[ipt]->Get(Form("hAlt_%dS",state));
    hAA_Pt_Dev[ipt] = (TH1D*) file_AA_Pt[ipt]->Get(Form("hDev_%dS",state));
  
    hPP_Pt_Nom[ipt]->GetXaxis()->CenterTitle();
    hPP_Pt_Alt[ipt]->GetXaxis()->CenterTitle();
    hPP_Pt_Dev[ipt]->GetXaxis()->CenterTitle();
    hAA_Pt_Nom[ipt]->GetXaxis()->CenterTitle();
    hAA_Pt_Alt[ipt]->GetXaxis()->CenterTitle();
    hAA_Pt_Dev[ipt]->GetXaxis()->CenterTitle();
    
    hPP_Pt_Nom[ipt] -> SetTitle(Form("PP N(%dS) nominal 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    hPP_Pt_Alt[ipt] -> SetTitle(Form("PP N(%dS) alternative 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    hPP_Pt_Dev[ipt] -> SetTitle(Form("PP N(%dS) deviation 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    hAA_Pt_Nom[ipt] -> SetTitle(Form("AA N(%dS) nominal 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    hAA_Pt_Alt[ipt] -> SetTitle(Form("AA N(%dS) alternative 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    hAA_Pt_Dev[ipt] -> SetTitle(Form("AA N(%dS) deviation 1M events, 100toys, %.1f < p_{T} < %.1f GeV/c",state,ptBin[i], ptBin[i+1]));
    
    hPP_Pt_Nom[ipt]->GetXaxis()->SetTitleOffset(1.4);
    hPP_Pt_Alt[ipt]->GetXaxis()->SetTitleOffset(1.4);
    hPP_Pt_Dev[ipt]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Pt_Nom[ipt]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Pt_Alt[ipt]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Pt_Dev[ipt]->GetXaxis()->SetTitleOffset(1.4);
    
    hPP_Pt_Nom[ipt]->GetYaxis()->SetTitleOffset(1.4);
    hPP_Pt_Alt[ipt]->GetYaxis()->SetTitleOffset(1.4);
    hPP_Pt_Dev[ipt]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Pt_Nom[ipt]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Pt_Alt[ipt]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Pt_Dev[ipt]->GetYaxis()->SetTitleOffset(1.4);
  }

  for(int irap=0;irap<nYBins;irap++)
  {
    hPP_Rap_Nom[irap] = (TH1D*) file_PP_Rap[irap] -> Get(Form("hNom_%dS",state));
    hPP_Rap_Alt[irap] = (TH1D*) file_PP_Rap[irap] -> Get(Form("hAlt_%dS",state));
    hPP_Rap_Dev[irap] = (TH1D*) file_PP_Rap[irap] -> Get(Form("hDev_%dS",state));
    hAA_Rap_Nom[irap] = (TH1D*) file_AA_Rap[irap] -> Get(Form("hNom_%dS",state));
    hAA_Rap_Alt[irap] = (TH1D*) file_AA_Rap[irap] -> Get(Form("hAlt_%dS",state));
    hAA_Rap_Dev[irap] = (TH1D*) file_AA_Rap[irap] -> Get(Form("hDev_%dS",state));
  
    hPP_Rap_Nom[irap]->GetXaxis()->CenterTitle();
    hPP_Rap_Alt[irap]->GetXaxis()->CenterTitle();
    hPP_Rap_Dev[irap]->GetXaxis()->CenterTitle();
    hAA_Rap_Nom[irap]->GetXaxis()->CenterTitle();
    hAA_Rap_Alt[irap]->GetXaxis()->CenterTitle();
    hAA_Rap_Dev[irap]->GetXaxis()->CenterTitle();
    
    hPP_Rap_Nom[irap] -> SetTitle(Form("PP N(%dS) nominal 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    hPP_Rap_Alt[irap] -> SetTitle(Form("PP N(%dS) alternative 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    hPP_Rap_Dev[irap] -> SetTitle(Form("PP N(%dS) deviation 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    hAA_Rap_Nom[irap] -> SetTitle(Form("AA N(%dS) nominal 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    hAA_Rap_Alt[irap] -> SetTitle(Form("AA N(%dS) alternative 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    hAA_Rap_Dev[irap] -> SetTitle(Form("AA N(%dS) deviation 1M events, 100toys, %.1f < |y| < %.1f GeV/c",state,yBin[i], yBin[i+1]));
    
    hPP_Rap_Nom[irap]->GetXaxis()->SetTitleOffset(1.4);
    hPP_Rap_Alt[irap]->GetXaxis()->SetTitleOffset(1.4);
    hPP_Rap_Dev[irap]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Rap_Nom[irap]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Rap_Alt[irap]->GetXaxis()->SetTitleOffset(1.4);
    hAA_Rap_Dev[irap]->GetXaxis()->SetTitleOffset(1.4);
    
    hPP_Rap_Nom[irap]->GetYaxis()->SetTitleOffset(1.4);
    hPP_Rap_Alt[irap]->GetYaxis()->SetTitleOffset(1.4);
    hPP_Rap_Dev[irap]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Rap_Nom[irap]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Rap_Alt[irap]->GetYaxis()->SetTitleOffset(1.4);
    hAA_Rap_Dev[irap]->GetYaxis()->SetTitleOffset(1.4);
  }

  for(int icent=0;icent<nCentBins;icent++)
  {
    hAA_cent_Nom[icent] = (TH1D*) file_AA_cent[icent] -> Get(Form("hNom_%dS",state));
    hAA_cent_Alt[icent] = (TH1D*) file_AA_cent[icent] -> Get(Form("hAlt_%dS",state));
    hAA_cent_Dev[icent] = (TH1D*) file_AA_cent[icent] -> Get(Form("hDev_%dS",state));
  
    hPP_cent_Nom[icent]->GetXaxis()->CenterTitle();
    hPP_cent_Alt[icent]->GetXaxis()->CenterTitle();
    hPP_cent_Dev[icent]->GetXaxis()->CenterTitle();
    hAA_cent_Nom[icent]->GetXaxis()->CenterTitle();
    hAA_cent_Alt[icent]->GetXaxis()->CenterTitle();
    hAA_cent_Dev[icent]->GetXaxis()->CenterTitle();
    
    hPP_cent_Nom[icent] -> SetTitle(Form("PP N(%dS) nominal 1M events, 100toys, cent. %d-%d GeV/c",state,centBin[i],centBin[i+1]));
    hPP_cent_Alt[icent] -> SetTitle(Form("PP N(%dS) alternative 1M events, 100toys,cent. %d-%d GeV/c",state,centBin[i],centBin[i+1]));
    hPP_cent_Dev[icent] -> SetTitle(Form("PP N(%dS) deviation 1M events, 100toys, cent. %d-%d GeV/c",state,centBin[i],centBin[i+1])); 
    hAA_cent_Nom[icent] -> SetTitle(Form("AA N(%dS) nominal 1M events, 100toys, cent. %d-%d GeV/c",state,centBin[i],centBin[i+1])); 
    hAA_cent_Alt[icent] -> SetTitle(Form("AA N(%dS) alternative 1M events, 100toys, cent. %d-%d GeV/c",state,centBin[i],centBin[i+1])); 
    hAA_cent_Dev[icent] -> SetTitle(Form("AA N(%dS) deviation 1M events, 100toys, cent. %d-%d GeV/c",state,centBin[i],centBin[i+1]));
    
    hPP_cent_Nom[icent]->GetXaxis()->SetTitleOffset(1.4);
    hPP_cent_Alt[icent]->GetXaxis()->SetTitleOffset(1.4);
    hPP_cent_Dev[icent]->GetXaxis()->SetTitleOffset(1.4);
    hAA_cent_Nom[icent]->GetXaxis()->SetTitleOffset(1.4);
    hAA_cent_Alt[icent]->GetXaxis()->SetTitleOffset(1.4);
    hAA_cent_Dev[icent]->GetXaxis()->SetTitleOffset(1.4);
    
    hPP_cent_Nom[icent]->GetYaxis()->SetTitleOffset(1.4);
    hPP_cent_Alt[icent]->GetYaxis()->SetTitleOffset(1.4);
    hPP_cent_Dev[icent]->GetYaxis()->SetTitleOffset(1.4);
    hAA_cent_Nom[icent]->GetYaxis()->SetTitleOffset(1.4);
    hAA_cent_Alt[icent]->GetYaxis()->SetTitleOffset(1.4);
    hAA_cent_Dev[icent]->GetYaxis()->SetTitleOffset(1.4);
  }
  
  hPP_int_Nom = (TH1D*) file_PP_Int -> Get(Form("hNom_%dS",state)); 
  hPP_int_Alt = (TH1D*) file_PP_Int -> Get(Form("hAlt_%dS",state)); 
  hPP_int_Dev = (TH1D*) file_PP_Int -> Get(Form("hDev_%dS",state)); 
  hAA_int_Nom = (TH1D*) file_AA_Int -> Get(Form("hNom_%dS",state));
  hAA_int_Alt = (TH1D*) file_AA_Int -> Get(Form("hAlt_%dS",state));
  hAA_int_Dev = (TH1D*) file_AA_Int -> Get(Form("hDev_%dS",state));
  
  hPP_cent_Nom->GetXaxis()->CenterTitle();
  hPP_cent_Alt->GetXaxis()->CenterTitle();
  hPP_cent_Dev->GetXaxis()->CenterTitle();
  hAA_cent_Nom->GetXaxis()->CenterTitle();
  hAA_cent_Alt->GetXaxis()->CenterTitle();
  hAA_cent_Dev->GetXaxis()->CenterTitle();

  hPP_cent_Nom -> SetTitle(Form("PP N(%dS) nominal 1M events, 100toys, integrated",state));
  hPP_cent_Alt -> SetTitle(Form("PP N(%dS) alternative 1M events, 100toys, integrated",state));
  hPP_cent_Dev -> SetTitle(Form("PP N(%dS) deviation 1M events, 100toys, integrated",state)); 
  hAA_cent_Nom -> SetTitle(Form("AA N(%dS) nominal 1M events, 100toys, integrated",state)); 
  hAA_cent_Alt -> SetTitle(Form("AA N(%dS) alternative 1M events, 100toys, integrated",state)); 
  hAA_cent_Dev -> SetTitle(Form("AA N(%dS) deviation 1M events, 100toys, integrated",state));

  hPP_cent_Nom->GetXaxis()->SetTitleOffset(1.4);
  hPP_cent_Alt->GetXaxis()->SetTitleOffset(1.4);
  hPP_cent_Dev->GetXaxis()->SetTitleOffset(1.4);
  hAA_cent_Nom->GetXaxis()->SetTitleOffset(1.4);
  hAA_cent_Alt->GetXaxis()->SetTitleOffset(1.4);
  hAA_cent_Dev->GetXaxis()->SetTitleOffset(1.4);

  hPP_cent_Nom->GetYaxis()->SetTitleOffset(1.4);
  hPP_cent_Alt->GetYaxis()->SetTitleOffset(1.4);
  hPP_cent_Dev->GetYaxis()->SetTitleOffset(1.4);
  hAA_cent_Nom->GetYaxis()->SetTitleOffset(1.4);
  hAA_cent_Alt->GetYaxis()->SetTitleOffset(1.4);
  hAA_cent_Dev->GetYaxis()->SetTitleOffset(1.4);

    
  
  double nSig_PP_Pt_nom[nPtBins] = {0.};
  double nSig_PP_Pt_alt[nPtBins] = {0.};
  double nSig_AA_Pt_nom[nPtBins] = {0.};
  double nSig_AA_Pt_alt[nPtBins] = {0.};
  
  double nSig_PP_Rap_nom[nYBins] = {0.};
  double nSig_PP_Rap_alt[nYBins] = {0.};
  double nSig_AA_Rap_nom[nYBins] = {0.};
  double nSig_AA_Rap_alt[nYBins] = {0.};

  double nSig_AA_cent_nom[nCentBins] = {0.};
  double nSig_AA_cent_alt[nCentBins] = {0.};
  
  double nSig_PP_int_nom = 0.;
  double nSig_PP_int_alt = 0.;
  double nSig_AA_int_nom = 0.;
  double nSig_AA_int_alt = 0.;
  
  double dev_PP_Pt[nPtBins] = {0.};  
  double dev_PP_Rap[nYBins] = {0.};  
  double dev_PP_int = 0.;  

  double dev_AA_Pt[nPtBins] = {0.};  
  double dev_AA_Rap[nYBins] = {0.};  
  double dev_AA_cent[nCentBins] = {0.};  
  double dev_AA_int = 0.;  

  double dev_Pt[nPtBins] = {0.}; 
  double dev_Rap[nYBins] = {0.}; 
  double dev_cent[nCentBins] = {0.}; 
  double dev_Int = 0.; 
  
  double rms_PP_Pt[nPtBins] = {0.};
  double rms_AA_Pt[nPtBins] = {0.};
  double rms_PP_Rap[nYBins] = {0.};
  double rms_AA_Rap[nYBins] = {0.};
  double rms_AA_cent[nCentBins] = {0.};
  double rms_PP_int = 0.;
  double rms_AA_int = 0.;


  //FitToData Values
  double rap_FTD_SR_nom_pp[2]={0.326937,0.320951};
  double rap_FTD_SR_alt_pp[2]={0.330472,0.328755};
  double rap_FTD_SR_dev_pp[2]={1.08103,2.43177};
  double rap_FTD_SR_nom_aa[2]={0.11099,0.0638839};
  double rap_FTD_SR_alt_aa[2]={0.115006,0.0718792};
  double rap_FTD_SR_dev_aa[2]={2.51013,12.5154};

  double int_FTD_SR_nom_pp = 0.325764;
  double int_FTD_SR_alt_pp = 0.32841;
  double int_FTD_SR_dev_pp = 0.812282;
  double int_FTD_SR_nom_aa = 0.10002;
  double int_FTD_SR_alt_aa = 0.100965;
  double int_FTD_SR_dev_aa = 0.945045;

  double pt_FTD_SR_nom_pp[3]={0.305984,0.313963,0.383257};
  double pt_FTD_SR_alt_pp[3]={0.326476,0.31364,0.381771};
  double pt_FTD_SR_dev_pp[3]={6.69699,0.102687,0.387764};
  double pt_FTD_SR_nom_aa[3]={0.0960655,0.0771805,0.146024};
  double pt_FTD_SR_alt_aa[3]={0.129783,0.0604429,0.133573};
  double pt_FTD_SR_dev_aa[3]={35.0986,21.6863,8.52615};

  double cent_FTD_SR_nom_aa[9]={0.0341179,0.121653,0.0977837,0.12798,0.126112,0.138736,0.111563,0.165479,0.254205};
  double cent_FTD_SR_alt_aa[9]={0.0213719,0.12329,0.0998268,0.131644,0.133368,0.144165,0.123964,0.150344,0.265456};
  double cent_FTD_SR_dev_aa[9]={37.3587,1.34595,2.08934,2.86332,5.7537,3.91269,11.1161,9.14564,4.42562};


  TCanvas *c_unc = new TCanvas("c_unc","",700,700);
  TGraphErrors *gr_Pt = new TGraphErrors();
  TGraphErrors *gr_Rap = new TGraphErrors();
  TGraphErrors *gr_Cent = new TGraphErrors();

  int nBins = 0;

  for(int ipt=0; ipt<nPtBins; ipt++)
  {
    c_PP_Pt[0] -> cd(ipt+1);
    hPP_Pt_Nom[ipt]->Draw();
    drawText("PP nominal pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nSig_PP_Pt_nom[ipt] = hPP_Pt_Nom[ipt]->GetMean();

    c_PP_Pt[1] -> cd(ipt+1);
    hPP_Pt_Alt[ipt]->Draw();
    drawText("PP alternative pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nSig_PP_Pt_alt[ipt] = hPP_Pt_Nom[ipt]->GetMean();
    
    c_PP_Pt[2] -> cd(ipt+1);
    hPP_Pt_Dev[ipt]->Draw();
    drawText("PP Dev.",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nBins = hPP_Pt_Dev[ipt]->GetNBinsX();
    double avg = 0;
    double hcontent =0;
    for(int ibin = 1; ibin<=nBins;ibin++)
    {
      avg += hPP_Pt_Dev[ipt]->GetBinContent(ibin)*TMath::Abs(hPP_Pt_Dev[ipt]->GetBinCenter(ibin));
      hcontent += hPP_Pt_Dev[ipt]->GetBinContent(ibin);
    }
    avg = avg/hcontent;
    dev_PP_Pt[ipt] = avg;
    
    for(int ibin=1; ibin<=nBins;ibin++)
    {
      rms_PP_Pt[ipt] += (avg-hPP_Pt_Dev[ipt]->GetBinContent(ibin))*(avg-hPP_Pt_Dev[ipt]->GetBinContent(ibin));
    }
    rms_PP_Pt[ipt] = TMath::Sqrt(rms_PP_Pt[ipt]/hcontent);
    
    c_AA_Pt[0] -> cd(ipt+1);
    hAA_Pt_Nom[ipt]->Draw();
    drawText("AA nominal pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nSig_AA_Pt_nom[ipt] = hAA_Pt_Nom[ipt]->GetMean();
    
    c_AA_Pt[1] -> cd(ipt+1);
    hAA_Pt_Alt[ipt]->Draw();
    drawText("AA alternative pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nSig_AA_Pt_alt[ipt] = hAA_Pt_Nom[ipt]->GetMean();
    
    
    c_AA_Pt[2] -> cd(ipt+1);
    hAA_Pt_Dev[ipt]->Draw();
    drawText("AA Dev.",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[ipt],ptBin[ipt+1]));
    gPad->Update();
    nBins = hPP_Pt_Dev[ipt]->GetNBinsX();
    avg=0;
    hcontent=0;
    for(int ibin = 1; ibin<=nBins;ibin++)
    {
      avg += hAA_Pt_Dev[ipt]->GetBinContent(ibin)*TMath::Abs(hAA_Pt_Dev[ipt]->GetBinCenter(ibin));
      hcontent += hAA_Pt_Dev[ipt]->GetBinContent(ibin);
    }
    avg = avg/hcontent;
    dev_AA_Pt[ipt] = avg;
    dev_AA_Pt[ipt] = hAA_Pt_Dev[ipt] -> GetMean();
    
    for(int ibin=1; ibin<=nBins;ibin++)
    {
      rms_AA_Pt[ipt] += (avg-hAA_Pt_Dev[ipt]->GetBinContent(ibin))*(avg-hAA_Pt_Dev[ipt]->GetBinContent(ibin));
    }
    rms_AA_Pt[ipt] = TMath::Sqrt(rms_AA_Pt[ipt]/hcontent);
    
    dev_Pt[ipt] = ( (100+dev_AA_Pt[ipt]) / (100+dev_PP_Pt[ipt]) - 1) * 100;
    rms_Pt[ipt] = (100+dev_AA_Pt[ipt]) / (100+dev_PP_Pt[ipt]) * TMath::Sqrt(TMath::Power(rms_AA_Pt[ipt]/dev_AA_Pt[ipt],2)+TMath::Power(rms_PP_Pt[ipt]/dev_PP_Pt[ipt],2)) * 100; 

    gr_Pt -> SetPoint(ipt,(ptBin[ipt]+ptBin[ipt+1])/2,dev_Pt[ipt]);
    gr_Pt -> SetPointError(ipt,0,rms_Pt[ipt]);
  }


  for(int irap=0; irap<nYBins; irap++)
  {
    c_PP_Rap[0] -> cd(irap+1);
    hPP_Rap_Nom[irap]->Draw();
    drawText("PP nominal pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nSig_PP_Rap_nom[irap] = hPP_Rap_Nom[irap]->GetMean();

    c_PP_Rap[1] -> cd(irap+1);
    hPP_Rap_Alt[irap]->Draw();
    drawText("PP alternative pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nSig_PP_Rap_alt[irap] = hPP_Rap_Nom[irap]->GetMean();
    
    c_PP_Rap[2] -> cd(irap+1);
    hPP_Rap_Dev[irap]->Draw();
    drawText("PP Dev.",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nBins = hPP_Rap_Dev[irap]->GetNBinsX();
    double avg = 0;
    double hcontent =0;
    for(int ibin = 1; ibin<=nBins;ibin++)
    {
      avg += hPP_Rap_Dev[irap]->GetBinContent(ibin)*TMath::Abs(hPP_Rap_Dev[irap]->GetBinCenter(ibin));
      hcontent += hPP_Rap_Dev[irap]->GetBinContent(ibin);
    }
    avg = avg/hcontent;
    dev_PP_Rap[irap] = avg;
    
    for(int ibin=1; ibin<=nBins;ibin++)
    {
      rms_PP_Rap[irap] += (avg-hPP_Rap_Dev[irap]->GetBinContent(ibin))*(avg-hPP_Rap_Dev[irap]->GetBinContent(ibin));
    }
    rms_PP_Rap[irap] = TMath::Sqrt(rms_PP_Rap[irap]/hcontent);
    
    c_AA_Rap[0] -> cd(irap+1);
    hAA_Rap_Nom[irap]->Draw();
    drawText("AA nominal pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nSig_AA_Rap_nom[irap] = hAA_Rap_Nom[irap]->GetMean();
    
    c_AA_Rap[1] -> cd(irap+1);
    hAA_Rap_Alt[irap]->Draw();
    drawText("AA alternative pdf",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nSig_AA_Rap_alt[irap] = hAA_Rap_Nom[irap]->GetMean();
    
    
    c_AA_Rap[2] -> cd(irap+1);
    hAA_Rap_Dev[irap]->Draw();
    drawText("AA Dev.",0.51,0.84,1,13);
    drawText("1M event generation",0.51,0.77,1,13);
    drawText(Form("%.1f < p_{T} < %.1f",ptBin[irap],ptBin[irap+1]));
    gPad->Update();
    nBins = hPP_Rap_Dev[irap]->GetNBinsX();
    avg=0;
    hcontent=0;
    for(int ibin = 1; ibin<=nBins;ibin++)
    {
      avg += hAA_Rap_Dev[irap]->GetBinContent(ibin)*TMath::Abs(hAA_Rap_Dev[irap]->GetBinCenter(ibin));
      hcontent += hAA_Rap_Dev[irap]->GetBinContent(ibin);
    }
    avg = avg/hcontent;
    dev_AA_Rap[irap] = avg;
    dev_AA_Rap[irap] = hAA_Rap_Dev[irap] -> GetMean();
    
    for(int ibin=1; ibin<=nBins;ibin++)
    {
      rms_AA_Rap[irap] += (avg-hAA_Rap_Dev[irap]->GetBinContent(ibin))*(avg-hAA_Rap_Dev[irap]->GetBinContent(ibin));
    }
    rms_AA_Rap[irap] = TMath::Sqrt(rms_AA_Rap[irap]/hcontent);
    
    dev_Rap[irap] = ( (100+dev_AA_Rap[irap]) / (100+dev_PP_Rap[irap]) - 1) * 100;
    rms_Rap[irap] = (100+dev_AA_Rap[irap]) / (100+dev_PP_Rap[irap]) * TMath::Sqrt(TMath::Power(rms_AA_Rap[irap]/dev_AA_Rap[irap],2)+TMath::Power(rms_PP_Rap[irap]/dev_PP_Rap[irap],2)) * 100; 

    gr_Rap -> SetPoint(irap,(ptBin[irap]+ptBin[irap+1])/2,dev_Rap[irap]);
    gr_Rap -> SetPointError(irap,0,rms_Rap[irap]);
  }

  

  //Rapidity and Integrated Bins
  for(int i=0;i<2;i++)
  {
    cpp_y1[0]->cd(i+1);
    hist3_pp_y1[i] -> Draw();
    drawText("PP nominal pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_nom_pp[i],0,rap_FTD_SR_nom_pp[i],gPad->GetFrame()->GetY2(),2,2);
    cpp_y1[1]->cd(i+1);
    hist3_pp_y2[i] -> Draw();
    drawText("PP 4th poly pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_alt_pp[i],0,rap_FTD_SR_alt_pp[i],gPad->GetFrame()->GetY2(),2,2);
    cpp_y1[2]->cd(i+1);
    hist3_pp_y3[i] -> Draw();
    drawText("PP",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.49,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_dev_pp[i],0,rap_FTD_SR_dev_pp[i],gPad->GetFrame()->GetY2(),2,2);
    cAA_y1[0]->cd(i+1);
    hist3_AA_y1[i] -> Draw();
    drawText("AA nominal pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_nom_aa[i],0,rap_FTD_SR_nom_aa[i],gPad->GetFrame()->GetY2(),2,2);
    cAA_y1[1]->cd(i+1);
    hist3_AA_y2[i] -> Draw();
    drawText("AA 4th poly pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_alt_aa[i],0,rap_FTD_SR_alt_aa[i],gPad->GetFrame()->GetY2(),2,2);
    cAA_y1[2]->cd(i+1);
    hist3_AA_y3[i] -> Draw();
    drawText("PbPb",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.49,0.71,1,13);
    drawText(Form("%.1f<|y|<%.1f",yBin[i],yBin[i+1]),0.49,0.64,1,13);
    gPad->Update();
    jumSun(rap_FTD_SR_dev_aa[i],0,rap_FTD_SR_dev_aa[i],gPad->GetFrame()->GetY2(),2,2);
     
    cint[0]->cd(i+1);
    hist3_int3[0] -> GetXaxis()->SetRange(0,300);
    hist3_int1[i]->Draw();

    if(i==0)
    {
      drawText("PP nominal pdf",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("pp integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_nom_pp,0,int_FTD_SR_nom_pp,gPad->GetFrame()->GetY2(),2,2);
    }
    else if(i==1)
    {
      drawText("AA nominal pdf",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("PbPb integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_nom_aa,0,int_FTD_SR_nom_aa,gPad->GetFrame()->GetY2(),2,2);
    }
    cint[1]->cd(i+1);
    hist3_int2[i]->Draw();
    
    if(i==0)
    {
      drawText("PP 4th poly pdf",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("pp integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_alt_pp,0,int_FTD_SR_alt_pp,gPad->GetFrame()->GetY2(),2,2);
    }
    else if(i==1)
    {
      drawText("AA 4th poly pdf",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("PbPb integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_alt_aa,0,int_FTD_SR_alt_aa,gPad->GetFrame()->GetY2(),2,2);
    }
    cint[2]->cd(i+1);
    hist3_int3[i]->Draw();
    if(i==0)
    {
      drawText("PP",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("pp integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_dev_pp,0,int_FTD_SR_dev_pp,gPad->GetFrame()->GetY2(),2,2);
    }
    else if(i==1)
    {
      drawText("PbPb",0.49,0.84,1,13);
      drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
      drawText("1000 toys",0.49,0.71,1,13);
      drawText("PbPb integrated",0.49,0.64,1,13);
      gPad->Update();
      jumSun(int_FTD_SR_dev_aa,0,int_FTD_SR_dev_aa,gPad->GetFrame()->GetY2(),2,2);
    }

    mean_int[i] = hist3_int3[i]->GetMean();
    rms_int[i] = hist3_int3[i]->GetRMS();
    rms_int[i] = rms_int[i]/TMath::Sqrt(1000);

    if(hist3_int1[i]->GetMean()>hist3_int2[i]->GetMean()) mean_int[i] = -mean_int[i];

    meany1_1[i] = hist3_pp_y1[i] -> GetMean();
    meany1_2[i] = hist3_pp_y2[i] -> GetMean();
    meany2_1[i] = hist3_AA_y1[i] -> GetMean();
    meany2_2[i] = hist3_AA_y2[i] -> GetMean();

    meany1[i] = hist3_pp_y3[i] -> GetMean();
    rmsy1[i] = hist3_pp_y3[i] -> GetRMS();
    rmsy1[i] = rmsy1[i]/TMath::Sqrt(1000);
    meany2[i] = hist3_AA_y3[i] -> GetMean();
    rmsy2[i] = hist3_AA_y3[i] -> GetRMS();
    rmsy2[i] = rmsy2[i]/TMath::Sqrt(1000);

    if(meany1_1[i]>meany1_2[i]) meany1[i] = -meany1[i];
    if(meany2_1[i]>meany2_2[i]) meany2[i] = -meany1[i];

  
    meany[i] = ((100+meany2[i])/(100+meany1[i]) - 1)*100;

    cout << endl;
    cout << endl;
    cout << "y bin : " << yBin[i] << " - " << yBin[i+1] << endl;
    cout << "Unc from pp : " << meany1[i] << endl;
    cout << "error : " << rmsy1[i] << endl;
    cout << endl;
    cout << "Unc from AA : " << meany2[i] << endl;
    cout << "error : " << rmsy2[i] << endl;
    cout << endl;
    cout << "Total Unc : " << meany[i] << endl;
    cout << endl;
    cout << "Integrated bin : " << i << endl;
    cout << "Unc : " << mean_int[i] << endl;
    cout << "error : " << rms_int[i] << endl;
    cout << endl;

  }

  for(int i=0;i<3;i++)
  {
    cpp_pT1[0]->cd(i+1);
    hist3_pp_pT1[i]->Draw();
    drawText("PP nominal pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<pT<%.1f",ptBin[i],ptBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_nom_pp[i],0,pt_FTD_SR_nom_pp[i],gPad->GetFrame()->GetY2(),2,2);
    
    cpp_pT1[1]->cd(i+1);
    hist3_pp_pT2[i] -> Draw();
    drawText("PP 4th poly",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<pT<%.1f",ptBin[i],ptBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_alt_pp[i],0,pt_FTD_SR_alt_pp[i],gPad->GetFrame()->GetY2(),2,2);
    
    cpp_pT1[2]->cd(i+1);
    hist3_pp_pT3[i]->Rebin(4);
    hist3_pp_pT3[i]->GetXaxis()->SetRange(0,60);
    hist3_pp_pT3[i]->Draw();
    drawText("PP",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.49,0.71,1,13);
    drawText(Form("%.1f<pT<%.1f GeV/c",ptBin[i],ptBin[i+1]),0.49,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_dev_pp[i],0,pt_FTD_SR_dev_pp[i],gPad->GetFrame()->GetY2(),2,2);

    meanpt1_1[i] = hist3_pp_pT1[i]->GetMean();
    meanpt1_2[i] = hist3_pp_pT2[i]->GetMean();

    meanpt1[i] = hist3_pp_pT3[i]->GetMean();  
    if(meanpt1_1[i]>meanpt1_2[i]) meanpt1[i] = -meanpt1[i];
    rmspt1[i] = hist3_pp_pT3[i]->GetRMS();  

    cAA_pT1[0]->cd(i+1);
    hist3_AA_pT1[i] -> Draw();
    drawText("AA nominal pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<pT<%.1f",ptBin[i],ptBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_nom_aa[i],0,pt_FTD_SR_nom_aa[i],gPad->GetFrame()->GetY2(),2,2);
    
    cAA_pT1[1]->cd(i+1);
    hist3_AA_pT2[i] -> Draw();
    drawText("AA 4th poly",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%.1f<pT<%.1f",ptBin[i],ptBin[i+1]),0.51,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_alt_aa[i],0,pt_FTD_SR_alt_aa[i],gPad->GetFrame()->GetY2(),2,2);
    
    cAA_pT1[2]->cd(i+1);
    hist3_AA_pT3[i]->Rebin(4);
    hist3_AA_pT3[i]->Draw();
    drawText("PbPb",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.49,0.71,1,13);
    drawText(Form("%.1f<p_{T}<%.1f GeV/c",ptBin[i],ptBin[i+1]),0.49,0.64,1,13);
    gPad->Update();
    jumSun(pt_FTD_SR_dev_aa[i],0,pt_FTD_SR_dev_aa[i],gPad->GetFrame()->GetY2(),2,2);

    meanpt2_1[i] = hist3_AA_pT1[i]->GetMean();
    meanpt2_2[i] = hist3_AA_pT2[i]->GetMean();

    meanpt2[i] = hist3_AA_pT3[i]->GetMean();  
    if(meanpt2_1[i]>meanpt2_2[i]) meanpt2[i] = -meanpt2[i];
    rmspt2[i] = hist3_AA_pT3[i]->GetRMS();  

    
    meanpt[i] = ((100+meanpt2[i])/(100+meanpt1[i]) - 1)*100;
    rmspt1[i]=rmspt1[i]/TMath::Sqrt(1000);
    rmspt2[i]=rmspt2[i]/TMath::Sqrt(1000);

    cout << endl;
    cout << endl;
    cout << "pT bin : " << ptBin[i] << " - " << ptBin[i+1] << " GeV/c" << endl;
    cout << "Unc from pp : " << meanpt1[i] << endl;
    cout << "error : " << rmspt1[i] << endl;
    cout << endl;
    cout << "Unc from AA : " << meanpt2[i] << endl;
    cout << "error : " << rmspt2[i] << endl;
    cout << endl;
    cout << "Total Unc : " << meanpt[i] << endl;
    cout << endl;

  }

  for(int i=0;i<9;i++)
  {
    c1[0]->cd(i+1);
    hist3_1[i]->Draw();
    drawText("AA Nominal pdf",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.49,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.49,0.64,1,13);
    gPad->Update();
    jumSun(cent_FTD_SR_nom_aa[i],0,cent_FTD_SR_nom_aa[i],gPad->GetFrame()->GetY2(),2,2);
    c1[1]->cd(i+1);
    hist3_2[i]->Draw();
    drawText("AA 4th order poly",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.51,0.64,1,13);
    gPad->Update();
    jumSun(cent_FTD_SR_alt_aa[i],0,cent_FTD_SR_alt_aa[i],gPad->GetFrame()->GetY2(),2,2);
    c1[2]->cd(i+1);
    hist3_3[i]->Rebin(4);
    hist3_3[i]->Draw();
    drawText("PbPb",0.49,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.49,0.77,1,13);
    drawText("1000 toys",0.49,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.49,0.64,1,13);
    gPad->Update();
    jumSun(cent_FTD_SR_dev_aa[i],0,cent_FTD_SR_dev_aa[i],gPad->GetFrame()->GetY2(),2,2);
    
    mean1[i]=hist3_3[i]->GetMean();
    rms1[i]=hist3_3[i]->GetRMS();

//    histV3_3[i]->Fit("f1");
//    value = f1->GetParameter(1);
/*    c2[0]->cd(i+1); 
    hist4_1[i]->Draw();
    drawText("AA Nominal pdf",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.51,0.64,1,13);
    c2[1]->cd(i+1); 
    hist4_2[i]->Draw();
    drawText("AA Nom+Exp",0.54,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.54,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.51,0.64,1,13);
    c2[2]->cd(i+1); 
    hist4_3[i]->Rebin(4);
    hist4_3[i]->Draw();
    drawText("AA Dev. Nom & Nom+Exp",0.51,0.84,1,13);
    drawText(Form("%s",ToySel.Data()),0.51,0.77,1,13);
    drawText("1000 toys",0.51,0.71,1,13);
    drawText(Form("%d - %d %s",CentBin[i]/2,CentBin[i+1]/2,perc.Data()),0.51,0.64,1,13);
*/

    if(hist3_1[i]->GetMean()>hist3_2[i]->GetMean()) mean1[i] = -mean1[i];

    rms1[i] = rms1[i]/TMath::Sqrt(1000);
    rms2[i] = rms2[i]/TMath::Sqrt(1000);
    //    hist4_3[i]->Fit("f1");
//    value1 = f1->GetParameter(1);
//    rms = TMath::Sqrt((value*value+value1*value1)/2);
    
    mean[i] = TMath::Sqrt((mean1[i]*mean1[i]+mean2[i]*mean2[i])/2);
    sigmaA[i] = TMath::Sqrt(2*(mean1[i]*mean1[i]*rms1[i]*rms1[i]+mean2[i]*mean2[i]*rms2[i]*rms2[i]));
    rms[i]=0.5*TMath::Power(mean[i],-0.5)*sigmaA[i];

    cout << endl;
    cout << Form("FINAL Nom & 4th mean %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << mean1[i] << endl;
    cout << Form("FINAL Nom & 4thrms %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << rms1[i] << endl;
    cout << endl;
    cout << Form("FINAL Nom & nomexp mean %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << mean2[i] << endl;
    cout << Form("FINAL Nom & nomexp rms %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << rms2[i] << endl;
    cout << endl;
    cout << Form("FINAL mean %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << mean[i] << endl;
    cout << Form("FINAL rms %d-%d %s : ",CentBin[i],CentBin[i+1],perc.Data())  << rms[i] << endl;
    cout << endl;
    
    c1[0]->Update();    
    c1[1]->Update();    
    c1[2]->Update();    
  }


  c1[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_NomSingle_Cent_100K_1000Toys.png");
  c1[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_AltSingle_Cent_100K_1000Toys.png");
  c1[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_Dev_Cent_100K_1000Toys.png");
  cint[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/4thPoly_NomSingle_int_100K_1000Toys.png");
  cint[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/4thPoly_AltSingle_int_100K_1000Toys.png");
  cint[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/4thPoly_Dev_int_100K_1000Toys.png");

  cpp_y1[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_NomSingle_rap_100K_1000Toys.png");
  cpp_y1[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_AltSingle_rap_100K_1000Toys.png");
  cpp_y1[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_Dev_rap_100K_1000Toys.png");
  cAA_y1[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_NomSingle_rap_100K_1000Toys.png");
  cAA_y1[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_AltSingle_rap_100K_1000Toys.png");
  cAA_y1[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_Dev_rap_100K_1000Toys.png");
  
  cpp_pT1[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_NomSingle_pt_100K_1000Toys.png");
  cpp_pT1[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_AltSingle_pt_100K_1000Toys.png");
  cpp_pT1[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/PP_4thPoly_Dev_pt_100K_1000Toys.png");
  cAA_pT1[0]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_NomSingle_pt_100K_1000Toys.png");
  cAA_pT1[1]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_AltSingle_pt_100K_1000Toys.png");
  cAA_pT1[2]->SaveAs("ffsys_BkgVar/finalbkg_Gen100K_2S/TOYMC_2S/plot/AA_4thPoly_Dev_pt_100K_1000Toys.png");
}

