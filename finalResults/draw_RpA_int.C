#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
#include "../commonUtility.h"
using namespace std;

TString ResultDir  = "nominalFits";
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0) ;

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 

void draw_RpA_int() //1 or 2 (1S or 2S)
{
  TH1::SetDefaultSumw2();
  setTDRStyle();
  writeExtraText = false;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 5; // 0: 15001, 1: ours
  double xmax = 2.85;
//  double relsys = 0.1;
  
  TFile* fEff1s = new TFile("../Efficiency/RootFiles/EffNomCor_SysRpA_1S.root");
  TFile* fAcc1s = new TFile("../Acceptance/20180724/acceptance_wgt_1S_20180724_2Dplot.root");
  TH1D* hEff1s = (TH1D*) fEff1s->Get("EffNomRatInt");
  TH1D* hAcc1spp = (TH1D*) fAcc1s->Get("hIntAccPP_1S");
  TH1D* hAcc1spa = (TH1D*) fAcc1s->Get("hIntAccPA_1S");
  double eff1s = hEff1s->GetBinContent(1);
  double acc1s = hAcc1spp->GetBinContent(1)/hAcc1spa->GetBinContent(1);

  TFile* fEff2s = new TFile("../Efficiency/RootFiles/EffNomCor_SysRpA_2S.root");
  TFile* fAcc2s = new TFile("../Acceptance/20180724/acceptance_wgt_2S_20180724_2Dplot.root");
  TH1D* hEff2s = (TH1D*) fEff2s->Get("EffNomRatInt");
  TH1D* hAcc2spp = (TH1D*) fAcc2s->Get("hIntAccPP_2S");
  TH1D* hAcc2spa = (TH1D*) fAcc2s->Get("hIntAccPA_2S");
  double eff2s = hEff2s->GetBinContent(1);
  double acc2s = hAcc2spp->GetBinContent(1)/hAcc2spa->GetBinContent(1);

  TFile* fEff3s = new TFile("../Efficiency/RootFiles/EffNomCor_SysRpA_3S.root");
  TFile* fAcc3s = new TFile("../Acceptance/20180724/acceptance_wgt_3S_20180724_2Dplot.root");
  TH1D* hEff3s = (TH1D*) fEff3s->Get("EffNomRatInt");
  TH1D* hAcc3spp = (TH1D*) fAcc3s->Get("hIntAccPP_3S");
  TH1D* hAcc3spa = (TH1D*) fAcc3s->Get("hIntAccPA_3S");
  double eff3s = hEff3s->GetBinContent(1);
  double acc3s = hAcc3spp->GetBinContent(1)/hAcc3spa->GetBinContent(1);

  
  valErr yieldPP;
  valErr yieldPA;
  yieldPP = getYield(1, kPPDATA, 0,30, 0, 1.93, 0,200,0,100); 
  yieldPA = getYield(1, kPADATA, 0,30, -1.93, 1.93, 0,200,0,100);

  TH1D* hpp = new TH1D("hpp_1s",";;;",1,0,100);  
  TH1D* hpa = new TH1D("hpa_1s",";;;",1,0,100);  
  hpp->SetBinContent(1,yieldPP.val);
  hpp->SetBinError(1,yieldPP.err);
  hpa->SetBinContent(1,yieldPA.val);
  hpa->SetBinError(1,yieldPA.err);
  hpa->Divide(hpp);
  double rpa_corr;
  rpa_corr = eff1s*acc1s*lumi_pp*1000/(lumi_pa*208);
  hpa->Scale(rpa_corr);
  double rpa_1s = hpa->GetBinContent(1);   
  double rpa_1s_err = hpa->GetBinError(1);   
  hpp->Reset(); hpa->Reset(); rpa_corr = 0.;

  yieldPP = getYield(2, kPPDATA, 0,30, 0, 1.93, 0,200,0,100); 
  yieldPA = getYield(2, kPADATA, 0,30, -1.93, 1.93, 0,200,0,100);
  hpp->SetBinContent(1,yieldPP.val);
  hpp->SetBinError(1,yieldPP.err);
  hpa->SetBinContent(1,yieldPA.val);
  hpa->SetBinError(1,yieldPA.err);
  hpa->Divide(hpp);
  rpa_corr = eff2s*acc2s*lumi_pp*1000/(lumi_pa*208);
  hpa->Scale(rpa_corr);
  double rpa_2s = hpa->GetBinContent(1);   
  double rpa_2s_err = hpa->GetBinError(1);   
  hpp->Reset(); hpa->Reset(); rpa_corr = 0.;
  
  yieldPP = getYield(3, kPPDATA, 0,30, 0, 1.93, 0,200,0,100); 
  yieldPA = getYield(3, kPADATA, 0,30, -1.93, 1.93, 0,200,0,100);
  hpp->SetBinContent(1,yieldPP.val);
  hpp->SetBinError(1,yieldPP.err);
  hpa->SetBinContent(1,yieldPA.val);
  hpa->SetBinError(1,yieldPA.err);
  hpa->Divide(hpp);
  rpa_corr = eff3s*acc3s*lumi_pp*1000/(lumi_pa*208);
  hpa->Scale(rpa_corr);
  double rpa_3s = hpa->GetBinContent(1);   
  double rpa_3s_err = hpa->GetBinError(1);   
  hpp->Reset(); hpa->Reset(); rpa_corr = 0.;
  
  TFile* fsys1s = new TFile("../Systematics/mergedSys_ups1s.root");
  TH1D* hInt1s = (TH1D*) fsys1s->Get("hintRPA_merged");
  TFile* fsys2s = new TFile("../Systematics/mergedSys_ups2s.root");
  TH1D* hInt2s = (TH1D*) fsys2s->Get("hintRPA_merged");
  TFile* fsys3s = new TFile("../Systematics/mergedSys_ups3s.root");
  TH1D* hInt3s = (TH1D*) fsys3s->Get("hintRPA_merged");
  
  double eysys1s = hInt1s->GetBinContent(1);
  double eysys2s = hInt2s->GetBinContent(1);
  double eysys3s = hInt3s->GetBinContent(1);


  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3.,12.};

  double exsys = 0.05;
  double exsys_align = 0.0;
  
  cout << "rpa_1s : " << rpa_1s << endl; 
  const int cn_1s =  3;
  double cpx_1s[cn_1s] =  {0.55-exsys_align, 1.469-exsys_align, 2.4-exsys_align};
  double cpx_1s_exsys[cn_1s] = {0.55+exsys_align, 1.469+exsys_align, 2.4+exsys_align};
  double cpy_1s[cn_1s] =  {rpa_1s,rpa_2s,rpa_3s}; 
  double cex_1s[cn_1s] =  {0., 0., 0};
  double cey_1s[cn_1s] =  {rpa_1s_err,rpa_2s_err,rpa_3s_err};
  double cexsys_1s[cn_1s] =  {exsys, exsys, exsys};
  double ceysys_1s[cn_1s] =  {rpa_1s*eysys1s,rpa_2s*eysys2s,rpa_3s*eysys3s};
  


  eysys1s = TMath::Sqrt(eysys1s*eysys1s+lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  eysys2s = TMath::Sqrt(eysys2s*eysys2s+lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  eysys3s = TMath::Sqrt(eysys3s*eysys3s+lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);


  cout << "$\\pt$, $y_{CM}$ integrated & " << Form("%.3f",rpa_1s) << " & " << Form("%.3f",rpa_1s_err) << " & " << Form("%.3f",eysys1s) << " \\\\ " << endl; 
  cout << "$\\pt$, $y_{CM}$ integrated & " << Form("%.3f",rpa_2s) << " & " << Form("%.3f",rpa_2s_err) << " & " << Form("%.3f",eysys2s) << " \\\\ " << endl; 
  cout << "$\\pt$, $y_{CM}$ integrated & " << Form("%.3f",rpa_3s) << " & " << Form("%.3f",rpa_3s_err) << " & " << Form("%.3f",eysys3s) << " \\\\ " << endl; 

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphAsymmErrors* gRAA_sys[nfile];
  
  gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
  gRAA_sys[0] = new TGraphAsymmErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, cexsys_1s, ceysys_1s, ceysys_1s); 
  

  //// graph style 
  SetGraphStyle(gRAA[0], 0, 0); 
  SetGraphStyleSys(gRAA_sys[0], 0); 
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(31); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.038);
  
  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.052);
  
  //// legend
//  TLegend *leg= new TLegend(0.804, 0.67, 0.994, 0.81);
  TLegend *leg= new TLegend(0.22, 0.81, 0.465, 0.9);
  SetLegendStyle(leg);
  leg->SetTextSize(0.042);
  leg->SetTextFont(22);
  leg -> SetHeader("");
  //leg -> SetHeader("#Upsilon's");
//  leg -> AddEntry(gRAA[0],"#Upsilon's","lp");
  leg -> AddEntry(gRAA[0],"#Upsilon(nS) R_{pPb}","lp");

  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.046);
  header->SetTextFont(62);
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.4); //1.3
 
  for(int i=0;i<1;i++)
  {
    gRAA[i]->GetXaxis()->SetBinLabel(10,"");
    gRAA_sys[i]->GetXaxis()->SetBinLabel(10,"");
  }
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14); //0.1
  gPad->SetTopMargin(0.067);
  gRAA_sys[0]->Draw("A5");
  gRAA[0]->Draw("P");
  leg->Draw("same");


  dashedLine(0,1.,xmax,1.,1,1);
  //// draw text
//  double sz_init = 0.87; double sz_step = 0.0535;
  double sz_init = 0.87; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.24, sz_init, "p_{T}^{#Upsilon} < 30 GeV");
//  globtex->DrawLatex(0.24, sz_init-sz_step-0.007, "|y_{CM}^{#Upsilon}| < 1.93");
  globtex->DrawLatex(0.92, sz_init-sz_step-0.0165, "p_{T}^{#Upsilon} < 30 GeV/c");
  globtex->DrawLatex(0.92, sz_init-sz_step*2.1-0.0165, "|y_{CM}^{#Upsilon}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
/*  globtex_label->DrawLatex(0.265, sz_init-sz_step*15.34, "#Upsilon(1S)");
  globtex_label->DrawLatex(0.522, sz_init-sz_step*15.34, "#Upsilon(2S)");
  globtex_label->DrawLatex(0.782, sz_init-sz_step*15.34, "#Upsilon(3S)");
// */
  globtex_label->DrawLatex(0.243, sz_init-sz_step*14.75, "#Upsilon(1S)");
  globtex_label->DrawLatex(0.505, sz_init-sz_step*14.75, "#Upsilon(2S)");
  globtex_label->DrawLatex(0.782, sz_init-sz_step*14.75, "#Upsilon(3S)");

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.1;
  TBox *globalUncBox = new TBox(0,1-sys_global_y_Lo,sys_global_x,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox -> SetLineWidth(2);
  globalUncBox -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("plots/RpPb_int.pdf");
  c1->SaveAs("plots/RpPb_int.png");

/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_pt.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nfile; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

  valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh,int cLow, int cHigh,   float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut,cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TFile* inf = new TFile(Form("../NominalFitResult/jaredFit/NominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("NominalFitResult/jaredFit/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/afs/cern.ch/work/j/jaebeom/private/Analysis/RpA502TeV/NominalFitResult/jaredFit/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //Santona
  //TFile* inf = new TFile(Form("/home/stuli/CMSResearch/UpsilonAna_Run1/GitRepo/WorkingPlots/FitResults/nomfitresults_upsilon_%s.root",kineLabel.Data()));
  //Jaebeom
  //TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_%s.root",kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/samba/UpsilonAnalysis/fitResultFiles/mcFit_MuPt4_2016_11_04/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

void stripErrorBars( TH1* h, double defaultErr  ) {
  
  for ( int i=1;  i<= h->GetNbinsX() ; i++) {
    h->SetBinError( i, defaultErr);
  }
}
