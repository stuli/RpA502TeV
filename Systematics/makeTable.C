#include <iostream>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include <sstream>
#include "TROOT.h"
#include "TSystem.h"

using namespace std;

void makeTable(int state = 2)
{
  TFile *wf = new TFile(Form("CBGaus_Variation/sys_CBGausVariaion_%ds.root",state),"read");
  TH1D* hIntPP = (TH1D*) wf->Get("hIntPP");
  TH1D* hIntAA = (TH1D*) wf->Get("hIntAA");
  TH1D* hptPP = (TH1D*) wf->Get("hptPP");
  TH1D* hptAA = (TH1D*) wf->Get("hptAA");
  TH1D* hrapPP = (TH1D*) wf->Get("hrapPP");
  TH1D* hrapAA = (TH1D*) wf->Get("hrapAA");
  TH1D* hcentAA = (TH1D*) wf->Get("hcentAA");
  TH1D* hIntAAoPP = (TH1D*) wf->Get("hIntAAoPP");
  TH1D* hptAAoPP = (TH1D*) wf->Get("hptAAoPP");
  TH1D* hrapAAoPP = (TH1D*) wf->Get("hrapAAoPP");
  TH1D* hcentAAoPP = (TH1D*) wf->Get("hcentAAoPP");

  if(state == 1){
    
cout << "                                    & pp                    & PbPb  & PbPb/pp hline " << endl;
cout << "Inegrated                           & " << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "                  & " << Form("%.2f",1.00*double(hIntAA->GetBinContent(1)*100)) << " & " << Form("%.3f",1.00*double(hIntAAoPP->GetBinContent(1))) << "     hline" << endl;
cout << "pt textless 2 GeVc                & " << Form("%.2f",1.00*double(hptPP->GetBinContent(1)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "2 textless pt textless 4 GeVc     & " << Form("%.2f",1.00*double(hptPP->GetBinContent(2)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "4 textless pt textless 6 GeVc       & " << Form("%.2f",1.00*double(hptPP->GetBinContent(3)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(3)*100)) << "   hline" << endl;
cout << "6 textless pt textless 9 GeVc      & " << Form("%.2f",1.00*double(hptPP->GetBinContent(4)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(4)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(4)*100)) << "   hline" << endl;
cout << "9 textless pt textless 12 GeVc     & " << Form("%.2f",1.00*double(hptPP->GetBinContent(5)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(5)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(5)*100)) << "   hline" << endl;
cout << "12 textless pt textless 30 GeVc     & " << Form("%.2f",1.00*double(hptPP->GetBinContent(6)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(6)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(6)*100)) << "   hline" << endl;
cout << "$|y|$ textless 0.4                  & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(1)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "0.4 textless $|y|$ textless 0.8     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(2)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "0.8 textless $|y|$ textless 1.2     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(3)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(3)*100)) << "   hline" << endl;
cout << "1.2 textless $|y|$ textless 1.6     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(4)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(4)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(4)*100)) << "   hline" << endl;
cout << "1.6 textless $|y|$ textless 2.0     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(5)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(5)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(5)*100)) << "   hline" << endl;
cout << "2.0 textless $|y|$ textless 2.4     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(6)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(6)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(6)*100)) << "   hline" << endl;
cout << "0-5                                 & multirow{9}{*}{" << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "}  & " << Form("%.2f",1.00*double(hcentAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(1)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "5-10                                &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(2)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "10-20                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(3)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "20-30                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(4)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(4)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "30-40                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(5)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(5)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "40-50                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(6)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(6)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "50-60                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(7)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(7)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "60-70                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(8)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(8)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "70-100                              &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(9)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(9)*100)) << "  hline" << endl;
  }
  
  else if(state == 2){
    
cout << "                                    & pp                    & PbPb  & PbPb/pp hline " << endl;
cout << "Inegrated                           & " << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "                  & " << Form("%.2f",1.00*double(hIntAA->GetBinContent(1)*100)) << " & " << Form("%.3f",1.00*double(hIntAAoPP->GetBinContent(1))) << "     hline" << endl;
cout << "pt textless 4 GeVc                  & " << Form("%.2f",1.00*double(hptPP->GetBinContent(1)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "4 textless pt textless 9 GeVc      & " << Form("%.2f",1.00*double(hptPP->GetBinContent(2)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "9 textless pt textless 30 GeVc     & " << Form("%.2f",1.00*double(hptPP->GetBinContent(3)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(3)*100)) << "   hline" << endl;
cout << "$|y|$ textless 0.8                  & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(1)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "0.8 textless $|y|$ textless 1.6                  & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(2)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "1.6 textless $|y|$ textless 2.4     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(3)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(3)*100)) << "   hline" << endl;
cout << "0-5                                 & multirow{9}{*}{" << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "}  & " << Form("%.2f",1.00*double(hcentAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(1)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "5-10                                &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(2)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "10-20                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(3)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(3)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "20-30                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(4)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(4)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "30-40                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(5)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(5)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "40-50                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(6)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(6)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "50-60                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(7)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(7)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "60-70                               &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(8)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(8)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "70-100                              &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(9)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(9)*100)) << "  hline" << endl;
  }
  else if(state == 3){
    
cout << "                                    & pp                    & PbPb  & PbPb/pp hline " << endl;
cout << "Inegrated                           & " << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "                  & " << Form("%.2f",1.00*double(hIntAA->GetBinContent(1)*100)) << " & " << Form("%.3f",1.00*double(hIntAAoPP->GetBinContent(1))) << "     hline" << endl;
cout << "pt textless 6 GeVc                  & " << Form("%.2f",1.00*double(hptPP->GetBinContent(1)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "6 textless pt textless 30 GeVc      & " << Form("%.2f",1.00*double(hptPP->GetBinContent(2)*100)) << "                 & " << Form("%.2f",1.00*double(hptAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hptAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "$|y|$ textless 1.2                  & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(1)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(1)*100)) << "   hline" << endl;
cout << "1.2 textless $|y|$ textless 2.4     & " << Form("%.2f",1.00*double(hrapPP->GetBinContent(2)*100)) << "                & " << Form("%.2f",1.00*double(hrapAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hrapAAoPP->GetBinContent(2)*100)) << "   hline" << endl;
cout << "0-30                                 & multirow{9}{*}{" << Form("%.2f",1.00*double(hIntPP->GetBinContent(1)*100)) << "}  & " << Form("%.2f",1.00*double(hcentAA->GetBinContent(1)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(1)*100)) << "  cline{1-1} cline{3-4}" << endl;
cout << "30-100                              &                      &  " << Form("%.2f",1.00*double(hcentAA->GetBinContent(2)*100)) << " & " << Form("%.2f",1.00*double(hcentAAoPP->GetBinContent(2)*100)) << "  hline" << endl;
  }
}
