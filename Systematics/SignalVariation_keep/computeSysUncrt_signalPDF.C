#include "../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../cutsAndBin.h"
#include <TGraphErrors.h>
using namespace std;

const int nVar = 5 ;     // 0 : nominal value

double getRMS(TString fileName = "", int state=1);
//double getRMS31(TString fileName = "");
double getUncVar(TString ppfilename = "", 
    TString pbpbfilename = "", int varnum = 1);

void computeSysUncrt_signalPDF(int state=3)
{ 
  double pp ;  // placeholder for pp
  double aa ; // placeholder for PbPb
  double finalUnc; 
  
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double* yBin;

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

  TH1D *hptAA  = new TH1D("hptAA", " uncertainty vs p_{T};p_{T}(GeV/c);uncertainty",nPtBins,ptBin);
  TH1D *hrapAA = new TH1D("hrapAA", "uncertainty vs rap;y;uncertainty",nYBins,yBin);
  TH1D *hcentAA= new TH1D("hcentAA","uncertainty vs cent;cent;uncertainty",nCentBins,centBin);
  TH1D *hIntAA = new TH1D("hIntAA", "uncertainty vs integrated bin;cent;uncertainty",0,0,200);
  TH1D *hptPP  = new TH1D("hptPP",  "uncertainty vs p_{T};p_{T}(GeV/c);uncertainty",nPtBins,ptBin);
  TH1D *hrapPP = new TH1D("hrapPP", "uncertainty vs rap;y;uncertainty",nYBins,yBin);
  TH1D *hIntPP = new TH1D("hIntPP", "uncertainty vs integrated bin;cent;uncertainty",0,0,200);

  TFile *fout=new TFile(Form("sys_signalPDFVariaion_%ds.root", state),"recreate");

  cout << endl;
  cout << " pT DEPENDENCE===================== " << endl;
  for ( int ipt = 1 ; ipt<= nPtBins ; ipt++ ) {
    cout << Form("pT : %.1f - %.1f GeV/c", (float)ptBin[ipt-1], ptBin[ipt]) << endl;
    cout << " pp :   " << endl; 
    pp = getRMS(Form("fitresults_upsilon_DoubleCB_PP_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0", ptBin[ipt-1], ptBin[ipt]), state); // pp
    cout << " AA :   " << endl; 
    aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI", ptBin[ipt-1], ptBin[ipt] ),state ); // PbPb
    finalUnc = sqrt( pp*pp + aa*aa );
    cout << "Unc. of double ratio  " << 1.000000* double(100*finalUnc) << "%%" << endl;
    cout << endl;
    hptPP->SetBinContent(ipt,pp);
    hptAA->SetBinContent(ipt,aa);
  }
    fout->cd();
    hptAA->Write();
    hptPP->Write();

  cout << endl;
  cout << " y DEPENDENCE===================== " << endl;
  for ( int iy = 1 ; iy<= nYBins ; iy++ ) {
    cout << Form("y : %.1f - %.1f", (float)yBin[iy-1], yBin[iy]) << endl;
    cout << " pp :   " << endl; 
    pp = getRMS(Form("fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0", yBin[iy-1], yBin[iy]),state ); // pp
    cout << " AA :   " << endl; 
    aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI", yBin[iy-1], yBin[iy]),state ); // PbPb
    finalUnc = sqrt( pp*pp + aa*aa );
    cout << "Unc. of double ratio  " << 1.000000* double(100*finalUnc) << "%%" << endl;
    hrapPP->SetBinContent(iy,pp);
    hrapAA->SetBinContent(iy,aa);
  }
    fout->cd();
    hrapAA->Write();
    hrapPP->Write();

  // centrality dependence :
  cout << endl;
  cout << " CENTRALITY DEPENDENCE===================== " << endl;
  cout << " pp :   " << endl; 
  pp = getRMS("fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0",state); // pp.  pT/y integrated
  hIntPP->SetBinContent(1,pp);
  fout->cd();
  hIntPP->Write();
  cout << " AA :   " << endl; 
  for ( int icent = 0 ; icent<= nCentBins ; icent++ ) 
  {
    if(icent == 0) cout << "Centrality 0 - 200" << endl;
    else cout << Form("Centrality  %0.f - %0.f", (double)centBin[icent-1], (double)centBin[icent]) << endl;
    if(icent == 0) {
    aa = getRMS("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI",state);
    hIntAA->SetBinContent(1,aa);
    fout->cd();
    hIntAA->Write();
    }
    else 
    {
      if((double)centBin[icent] <= 60) aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality%0.f-%0.f_dphiEp_0.00PI_100.00PI", (double)centBin[icent-1], (double)centBin[icent]),state );
      else if((double)centBin[icent] > 60) aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality%0.f-%0.f_dphiEp_0.00PI_100.00PI", (double)centBin[icent-1], (double)centBin[icent]),state );
    hcentAA->SetBinContent(icent,aa);
    }
    finalUnc = sqrt( pp*pp + aa*aa );
    cout << "Unc. of double ratio  " << 0.00001* int(10000000*finalUnc) << "%%" << endl;
    cout << endl;
  }
    fout->cd();
    hcentAA->Write();
}

double getUncVar(TString ppfileName, TString pbpbfileName, int varnum)
{
  double r21pp[2];
  double relVarpp;   //pp relative variance
  double r21pbpb[2];
  double relVarpbpb;   //pbpb relative variance
  double varnum_unc;
  
  //nominal value
  TString finpp_nominal = Form("%s_0.root", ppfileName.Data() );
  TFile* fpp_nominal = new TFile(finpp_nominal.Data() );
  TH1D* hpp_nominal = (TH1D*)fpp_nominal->Get("fitResults");
  TString finpbpb_nominal = Form("%s_0.root", pbpbfileName.Data() );
  TFile* fpbpb_nominal = new TFile(finpbpb_nominal.Data() );
  TH1D* hpbpb_nominal = (TH1D*)fpbpb_nominal->Get("fitResults");


  TString finpp = Form("%s_%d.root", ppfileName.Data(), varnum);
  TFile* fpp = new TFile(finpp.Data() );
  TH1D* hpp = (TH1D*)fpp->Get("fitResults");
  

  TString finpbpb = Form("%s_%d.root", pbpbfileName.Data(), varnum);
  TFile* fpbpb = new TFile(finpbpb.Data() );
  TH1D* hpbpb = (TH1D*)fpbpb->Get("fitResults");
  
  r21pp[0] = (double)hpp_nominal->GetBinContent(4);  
  r21pbpb[0] = (double)hpbpb_nominal->GetBinContent(4);  
  r21pp[1] = (double)hpp->GetBinContent(4);
  r21pbpb[1] = (double)hpbpb->GetBinContent(4);
  relVarpp = r21pp[1]/r21pp[0] -1 ;
  relVarpbpb = r21pbpb[1]/r21pbpb[0] -1 ;
  varnum_unc = TMath::Sqrt(relVarpp*relVarpp + relVarpbpb*relVarpbpb);
  
  cout << " varnum_unc : " << varnum_unc * 100 << endl;
  return varnum_unc;
}

double getRMS(TString fileName , int state) 
{
  double yield[nVar+1];
  double relVar[nVar+1];   // relative variance
  double quadSum;
  TString fin;
  cout << "Var1 - Var5 : " ;
  for ( int ivar=0 ; ivar<=nVar ; ivar++) 
  { 
    if(ivar ==0) fin = Form("%s.root",fileName.Data());
    else fin = Form("Sys_SignalVar_%s_%d.root", fileName.Data(), ivar);
    TFile* f1 = new TFile(fin.Data() );
    TH1D* h = (TH1D*)f1->Get("fitResults");
    yield[ivar] = (double)h->GetBinContent(state);
    relVar[ivar] = yield[ivar]/yield[0] -1 ;
    if ( ivar>0&&ivar<5) 
      cout << " & " << Form("%.6f",100*double(relVar[ivar])) << "\\\%";
      //cout << " & " << 0.1 * double(1000*relVar[ivar]) << "%%";
    if ( ivar==5) 
      cout << " & \\multicolumn{1}{c|}{" << Form("%.6f",100*double(relVar[ivar])) << "\\\%}";
  }
  
  quadSum=0;
  for ( int ivar=0 ; ivar<=nVar ; ivar++) 
  {
    quadSum = quadSum + relVar[ivar]*relVar[ivar] ;
  }

  quadSum = quadSum/(double)nVar;
  double theRMS = sqrt(quadSum);
  cout << " & " << Form("%.6f",100*double(theRMS)) << " \\\\"<< endl;
  return theRMS;
}

