#include "../../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../../cutsAndBin.h"
#include <TGraphErrors.h>
using namespace std;

const int nVar = 1 ;     // 0 : nominal value

double getRMS(TString fileName = "DR_SignalSys/", int state=1);
double getFinalUnc(double aa = 1, double pp = 1);

void computeSysUncrt_bkgPDF(int state=1)
{ 
  double pp ;  // placeholder for pp
  double aa ; // placeholder for PbPb
  double finalUnc; 
  
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
  
  TH1D *hUncPtPP = new TH1D("hptPP",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapPP = new TH1D("hrapPP",";Uncertainty;|y|",nYBins,yBin);
  TH1D *hUncIntPP = new TH1D("hIntPP",";Uncertainty;Centrality",1,0,200);
  
  TH1D *hUncCentAA = new TH1D("hcentAA","; Uncertainty ; Centrality Bins", nCentBins, centBin);
  TH1D *hUncPtAA = new TH1D("hptAA",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapAA = new TH1D("hrapAA",";Uncertainty;|y|",nYBins,yBin);
  TH1D *hUncIntAA = new TH1D("hIntAA",";Uncertainty;Centrality",1,0,200);
  
  TH1D *hUncIntAAoPP = new TH1D("hIntAAoPP",";Uncertainty;Centrality",1,0,200);
  TH1D *hUncCentAAoPP = new TH1D("hcentAAoPP","; Uncertainty ; Centrality Bins", nCentBins, centBin);
  TH1D *hUncPtAAoPP = new TH1D("hptAAoPP",";Uncertainty;p_{T} (GeV/c)",nPtBins,ptBin);
  TH1D *hUncRapAAoPP = new TH1D("hrapAAoPP",";Uncertainty;|y|",nYBins,yBin);


  cout << "state : " << state << endl;
  // centrality dependence :
  cout << endl;
  cout << " CENTRALITY DEPENDENCE===================== " << endl;
  cout << " pp :   " << endl; 
  pp = getRMS("fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0",state); // pp.  pT/y integrated
  hUncIntPP -> SetBinContent(1,pp);
  cout << " AA :   " << endl; 
  for ( int icent = 0 ; icent<= nCentBins ; icent++ ) 
  {
    cout << " icent : " << icent << endl;
    if(icent == 0) cout << "Centrality 0 - 200" << endl;
    else cout << Form("Centrality  %.0f - %.0f", centBin[icent-1], centBin[icent]) << endl;
    if(icent == 0){ 
      aa = getRMS("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI",state);
    }
    else 
    {
      if(centBin[icent] <= 60) aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality%.0f-%.0f_dphiEp_0.00PI_100.00PI", centBin[icent-1], centBin[icent]),state );
      else if(centBin[icent] > 60) aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0_centrality%.0f-%.0f_dphiEp_0.00PI_100.00PI", centBin[icent-1], centBin[icent]),state );
    }
    finalUnc = getFinalUnc(pp,aa);
    cout << "finalUnc: " << finalUnc << "     cent : " << centBin[icent-1] << " - " << centBin[icent] << endl;
    if(icent == 0) {
      hUncIntAA -> SetBinContent(1,aa);
      hUncIntAAoPP -> SetBinContent(1,finalUnc);
    } 
    else if(icent != 0) {
      hUncCentAA -> SetBinContent(icent, aa);
      hUncCentAAoPP -> SetBinContent(icent,finalUnc);
    }
    cout << "Unc. of double ratio  " << 0.00001* double(10000000*finalUnc) << "%%" << endl;
    cout << endl;
  }

  cout << endl;
  cout << " pT DEPENDENCE===================== " << endl;
  cout << " pp :   " << endl; 
  cout << " AA :   " << endl; 
  for ( int ipt = 1 ; ipt<= nPtBins ; ipt++ ) {
    cout << Form("pT : %.1f - %.1f GeV/c", (float)ptBin[ipt-1], ptBin[ipt]) << endl;
    pp = getRMS(Form("fitresults_upsilon_DoubleCB_PP_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0", ptBin[ipt-1], ptBin[ipt]), state); // pp
    aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt%.1f-%.1f_y0.0-2.4_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI", ptBin[ipt-1], ptBin[ipt] ),state ); // PbPb
    finalUnc = getFinalUnc(pp,aa);
    hUncPtPP->SetBinContent(ipt,pp);
    hUncPtAA->SetBinContent(ipt,aa);
    hUncPtAAoPP->SetBinContent(ipt,finalUnc);
    cout << "Unc. of double ratio  " << 1.000000* double(100*finalUnc) << "%%" << endl;
    cout << endl;
  }

  cout << endl;
  cout << " y DEPENDENCE===================== " << endl;
  cout << " pp :   " << endl; 
  cout << " AA :   " << endl; 
  for ( int iy = 1 ; iy<= nYBins ; iy++ ) {
    cout << Form("y : %.1f - %.1f", (float)yBin[iy-1], yBin[iy]) << endl;
    pp = getRMS(Form("fitresults_upsilon_DoubleCB_PP_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0", yBin[iy-1], yBin[iy]),state ); // pp
    aa = getRMS(Form("fitresults_upsilon_DoubleCB_AA_DATA_pt0.0-30.0_y%.1f-%.1f_muPt4.0_centrality0-200_dphiEp_0.00PI_100.00PI", yBin[iy-1], yBin[iy]),state ); // PbPb
    finalUnc = getFinalUnc(pp,aa);
    hUncRapPP->SetBinContent(iy,pp);
    hUncRapAA->SetBinContent(iy,aa);
    hUncRapAAoPP->SetBinContent(iy,finalUnc);
    cout << "Unc. of double ratio  " << 1.000000* double(100*finalUnc) << "%%" << endl;
  }

  TFile *wf = new TFile(Form("sys_bkgPDFVariaion_linear_%ds.root",state),"recreate");

  wf->cd();
  hUncRapPP->Write();
  hUncRapAA->Write();
  hUncRapAAoPP->Write();
  hUncPtPP->Write();
  hUncPtAA->Write();
  hUncPtAAoPP->Write();
  hUncIntAA->Write();
  hUncIntPP->Write();
  hUncIntAAoPP->Write();
  hUncCentAA->Write();
  hUncCentAAoPP->Write();

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
    if(ivar ==0) fin = Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/NomPlot/%s.root",fileName.Data());
    else fin = Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/Systematic/BkgVariation/bkglinear/Sys_BkgVar_%s.root", fileName.Data());
    TFile* f1 = new TFile(fin.Data() );
    TH1D* h = (TH1D*)f1->Get("fitResults");
    yield[ivar] = (double)h->GetBinContent(state);
    relVar[ivar] = yield[ivar]/yield[0] -1 ;
    if ( ivar>0) 
      cout << " & " << Form("%.6f",100*double(relVar[ivar])) << "%%";
      //cout << " & " << 0.1 * double(1000*relVar[ivar]) << "%%";
  }
  
  quadSum=0;
  for ( int ivar=0 ; ivar<=nVar ; ivar++) 
  {
    quadSum = quadSum + relVar[ivar]*relVar[ivar] ;
  }

  quadSum = quadSum/(double)nVar;
  double theRMS = sqrt(quadSum);
  cout << " & " << Form("%.6f",100*double(theRMS)) << endl;
  return relVar[1];
}

double getFinalUnc(double pp, double aa)
{
  double FinalUnc;

  FinalUnc = ((1+aa)/(1+pp)-1);
  return FinalUnc;
}
  



