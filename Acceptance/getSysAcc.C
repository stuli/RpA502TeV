#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");

void getSysAcc(int state = 1) { 
  
  TH1::SetDefaultSumw2();

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

  TFile* f0 = new TFile("acceptance.root");
  TH1D* hptAccPP = (TH1D*)f0->Get(Form("hptAccPP%dS",state));
  TH1D* hrapAccPP = (TH1D*)f0->Get(Form("hrapAccPP%dS",state));
  TH1D* hcentAccPP = (TH1D*)f0->Get(Form("hIntAccPP%dS",state));  // Integrated bin for pp.
  TH1D* hptAccAA = (TH1D*)f0->Get(Form("hptAccAA%dS",state));
  TH1D* hrapAccAA = (TH1D*)f0->Get(Form("hrapAccAA%dS",state));
  TH1D* hcentAccAA_int = (TH1D*)f0->Get(Form("hIntAccAA%dS",state));
  TH1D* hcentAccAA = new TH1D(Form("hcentAccAA%dS",state), "",nCentBins, centBin);
  for ( int icent = 1 ; icent<= nCentBins ;icent++) {
    hcentAccAA->SetBinContent( icent ,hcentAccAA_int->GetBinContent( 1 ) );
    hcentAccAA->SetBinError( icent ,hcentAccAA_int->GetBinError( 1 ) );
  }

  TH1D* hptAccNoW = (TH1D*)f0->Get(Form("hptAccNoW%dS",state));
  TH1D* hrapAccNoW = (TH1D*)f0->Get(Form("hrapAccNoW%dS",state));
  TH1D* hcentAccNoW_int = (TH1D*)f0->Get(Form("hIntAccNoW%dS",state));  // Integrated bin for NoW.
  TH1D* hcentAccNoW = new TH1D(Form("hcentAccNoW%dS",state), "",nCentBins, centBin);
  for ( int icent = 1 ; icent<= nCentBins ;icent++) {
    hcentAccNoW->SetBinContent( icent ,hcentAccNoW_int->GetBinContent( 1 ) );
    hcentAccNoW->SetBinError( icent ,hcentAccNoW_int->GetBinError( 1 ) );
  }
  
  TH1D* hptSysPP = (TH1D*)hptAccPP->Clone("hptSysPP");
  TH1D* hrapSysPP = (TH1D*)hrapAccPP->Clone("hrapSysPP");
  TH1D* hcentSysPP = (TH1D*)hcentAccPP->Clone("hcentSysPP");

  TH1D* hptSysAA = (TH1D*)hptAccAA->Clone("hptSysAA");
  TH1D* hrapSysAA = (TH1D*)hrapAccAA->Clone("hrapSysAA");
  TH1D* hcentSysAA = (TH1D*)hcentAccAA->Clone("hcentSysAA");
  TH1D* hcentSysAA_int = (TH1D*)hcentAccAA_int->Clone("hcentSysAA_int");

  hptSysPP->Add(hptAccNoW, -1);
  hrapSysPP->Add(hrapAccNoW, -1);
  hcentSysPP->Add(hcentAccNoW_int, -1);
  hptSysAA->Add(hptAccNoW, -1);
  hrapSysAA->Add(hrapAccNoW, -1);
  hcentSysAA->Add(hcentAccNoW, -1);
  hcentSysAA_int->Add(hcentAccNoW_int, -1);

  cout << " here1 " << endl;
  hptSysPP->Divide(hptAccNoW);
  cout << " here2 " << endl;
  hrapSysPP->Divide(hrapAccNoW);
  cout << " here3 " << endl;
  hcentSysPP->Divide(hcentAccNoW_int);
  cout << " here4 " << endl;
  hptSysAA->Divide(hptAccNoW);
  cout << " here5 " << endl;
  hrapSysAA->Divide(hrapAccNoW);
  cout << " here6 " << endl;
  hcentSysAA->Divide(hcentAccNoW);
  cout << " here7 " << endl;
  hcentSysAA_int->Divide(hcentAccNoW_int);
  cout << " here8 " << endl;

  hptSysPP->SetAxisRange(-0.1,0.1,"Y");;   
  hptSysAA->SetAxisRange(-0.1,0.1,"Y");;
  hrapSysPP->SetAxisRange(-0.1,0.1,"Y");;   
  hrapSysAA->SetAxisRange(-0.1,0.1,"Y");;
  hcentSysAA_int->SetAxisRange(-0.1,0.1,"Y");;  
  hcentSysAA->SetAxisRange(-0.1,0.1,"Y");;
  hcentSysPP->SetAxisRange(-0.1,0.1,"Y");;

  
  TCanvas* c = new TCanvas("c1","",800,800);
  c->Divide(2,4);
  c->cd(1) ; hptSysPP->Draw();   c->cd(2) ; hptSysAA->Draw();
  c->cd(3) ; hrapSysPP->Draw();   c->cd(4) ; hrapSysAA->Draw();
  c->cd(5) ; hcentSysAA_int->Draw();   c->cd(6) ; hcentSysPP->Draw();   
  c->cd(7) ; hcentSysAA->Draw();
 
  c->SaveAs(Form("sys_acceptance_ups%d.pdf",state));

  TFile *fout = new TFile(Form("sys_acceptance_ups%d.root",state),"recreate");
  hptSysPP->Write();   
  hptSysAA->Write();
  hrapSysPP->Write();   
  hrapSysAA->Write();
  hcentSysAA_int->Write();  
  hcentSysAA->Write();
  hcentSysPP->Write();
    
  fout->Close();

}

