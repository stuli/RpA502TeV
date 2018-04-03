#include "../SONGKYO.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_raaCent.C"
#include "../../cutsAndBin.h"
using namespace std;

TString ResultDir  = "nominalFits";
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0,int cLow=0, int cHigh=0, float dphiEp2Low=0, float dphiEp2High=0) ;

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 
void compare_RFB_ALICE(int istate=1) //1 or 2 (1S or 2S)
{
  TH1::SetDefaultSumw2();
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 2; // 0: ALICE, 1: ATLAS, 2: LHCb, 3: CMS
  double xmin = 0;
  double xmax = 4.3;
//  double relsys = 0.1;
  
  double exsys_1s[1] =  {0.965};
  double exsys_2s[5] =  {0.235, 1.13/2, 0.4,0.4,1.13/2};
  double exsys_3s[3] =  {0.235, 1.93/2, 1.93/2};

  //// LHCb values
  const int cn_lhcb =  1;
  double cpx_lhcb[cn_lhcb] =  {3.25};
  double cpy_lhcb[cn_lhcb] =  {0.75};
  double cex_lhcb[cn_lhcb] =  {0.};
  double cey_lhcb[cn_lhcb] =  {0.16};
  double cexsys_lhcb[cn_lhcb] =  {0.75};
  double ceysys_lhcb[cn_lhcb] = {0.08};

  valErr yieldPP;
  valErr yieldPA;
  yieldPP = getYield(1, kPADATA, 0,30, -1.93, 0, 0,200,0,100); 
  yieldPA = getYield(1, kPADATA, 0,30, 0, 1.93, 0,200,0,100);

  TFile *f_acc = new TFile(Form("../../Acceptance/20180305/acceptance_wgt_%dS_20180305_2Dplot.root",istate),"read");
  TFile *f_eff = new TFile(Form("../../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_SysRFB_%dS.root",istate),"read");

  TH1D* hRFB1 = new TH1D("hRFB",";;;",1,0,1.93);
  hRFB1->SetBinContent(1,yieldPA.val);
  hRFB1->SetBinError(1,yieldPA.err);
  TH1D* hRFB2 = new TH1D("hRFB",";;;",1,0,1.93);
  hRFB2->SetBinContent(1,yieldPP.val);
  hRFB2->SetBinError(1,yieldPP.err);
  hRFB1->Divide(hRFB2);

  TH1D* hAcc_div = (TH1D*) hRFB2->Clone("hAcc_div"); hAcc_div->Reset();
  TH1D* hEff_div = (TH1D*) hRFB2->Clone("hAcc_div"); hEff_div->Reset();
  TH1D* hAcc_r = (TH1D*) f_acc->Get("hrapAccPA2bin_1S");
  TH1D* hEff_r = (TH1D*) f_eff->Get("EffNomRatIntRFB");
  double corr_acc = hAcc_r->GetBinContent(1)/hAcc_r->GetBinContent(2);
  double corr_eff = hEff_r->GetBinContent(1);
  for(int i=1;i<=hAcc_div->GetNbinsX();i++){
    hAcc_div->SetBinContent(i,corr_acc);
    hEff_div->SetBinContent(i,corr_eff);
  }
  hRFB1->Multiply(hAcc_div);
  hRFB1->Multiply(hEff_div);
  

  // eff 
  TFile* feff = new TFile("../../CrossChecks/efficiency_Santona/ForAnalysisNote/RootFiles/EffNomCor_SysRFB_1S.root");
  TH1D* heff = (TH1D*)feff->Get("EffSysIntRFB");
  // acceptance
  TFile* facc = new TFile("../../Acceptance/20180305/sys_acceptance_ups1S_20180305.root");
  TH1D* hrapPA_rfbint = (TH1D*) facc->Get("hrapSysAccPA2bin");
  double sys1 = hrapPA_rfbint->GetBinContent(1);
  double sys2 = hrapPA_rfbint->GetBinContent(2);
  double sysf = TMath::Abs(1-(1+sys1)/(1+sys2));
  //bkg
  TFile *fbkg = new TFile("../../Systematics/Graham_BkgVariation/BkgPdfSystematics.root");
  TH1D* hbkg = (TH1D*)fbkg->Get("hInt1S_Rfbdiff");
  //sig
  TFile *fsig = new TFile("../../Systematics/Jared_SignalShapeVariation/ErrorEstimates/SysSig1s.root");
  TH1D* hsig = (TH1D*)fsig->Get("hSysSigRFBIntActivity");
  
  double sys0 = heff->GetBinContent(1);
  double sys3 = hbkg->GetBinContent(1);
  double sys4 = hsig->GetBinContent(1);
  double sys_tot = TMath::Sqrt(sys0*sys0+sysf*sysf+sys3*sys3+sys4*sys4);


  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphErrors* gRAA_sys[nfile];
  //// LHC exp
    gRAA[0] = new TGraphErrors(cn_lhcb, cpx_lhcb, cpy_lhcb, cex_lhcb, cey_lhcb); 
    gRAA_sys[0] = new TGraphErrors(cn_lhcb, cpx_lhcb, cpy_lhcb, cexsys_lhcb, ceysys_lhcb); 

  //// 2) ours
  gRAA[1]= new TGraphErrors(hRFB1);
  gRAA_sys[1] = new TGraphErrors(hRFB1);

  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  int npoint = 1; 
  if (npoint != gRAA[1]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
    gRAA[1]->GetPoint(ipt, pxtmp, pytmp);
    cout << "pxtmp : " << pxtmp << endl;
    extmp=gRAA[1]->GetErrorX(ipt);
    eytmp=gRAA[1]->GetErrorY(ipt);
    relsys=sys_tot;
    // 1) remove ex from gRAA
    gRAA[1]->SetPointError(ipt, 0, eytmp);
    // 2) set ey for gRAA_sys (assign 10% temporarily)
    //gRAA_sys[1]->SetPointError(ipt, extmp, pytmp*relsys);
    gRAA_sys[1]->SetPointError(ipt, exsys_1s[ipt],pytmp*relsys);
  }
  ////////////////////////////////////////////////////////////////
  
  //// graph style 
  SetGraphStyle_comp(gRAA[0], 1, 1); 
  SetGraphStyleSys_comp(gRAA_sys[0], 1); 
  SetGraphStyle(gRAA[1], 0, 0); 
  SetGraphStyleSys(gRAA_sys[1], 0); 

  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  //// legend
  TLegend *leg= new TLegend(0.20, 0.70, 0.56, 0.91);
  SetLegendStyle(leg);
  leg -> SetTextSize(0.034);
  leg -> AddEntry(gRAA[1],Form("#varUpsilon(%dS), p_{T}^{#varUpsilon} < 30 GeV/c",istate),"lp");
  leg -> AddEntry(gRAA[0],Form("LHCb #varUpsilon(%dS), p^{#varUpsilon}_{T} < 15 GeV/c",istate),"lp");
  
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("|y_{CM}^{#varUpsilon}|");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetTitleOffset(1.0);
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{FB}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0,xmax);
  gRAA_sys[0]->GetXaxis()->SetNdivisions(505);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.8);
  /// for rap
 
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nfile; is++){
    if ( is==0) gRAA_sys[is]->Draw("A5");
    else gRAA_sys[is]->Draw("5");
    gRAA[is]->Draw("P");
	}
  dashedLine(xmin,1.,xmax,1.,1,1);
  leg->Draw("same");
//  leg_comp->Draw("same");

  //// draw text
  double sz_init = 0.927; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "|y|^{#mu#mu} < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  
  double TAA_unc_Global_Hi = 0.068;
  double TAA_unc_Global_Lo = 0.072;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.4;
  //double sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_aa*lumi_unc_aa);
  double sys_global_y_alice = 0.016;
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("../plots/%dS_comp_RFB_vs_rap_LHCb.pdf",istate));
  c1->SaveAs(Form("../plots/%dS_comp_RFB_vs_rap_LHCb.png",istate));

/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_rap.root", "RECREATE");
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
  TFile* inf = new TFile(Form("../../NominalFitResult/jaredFit/NominalFits/nomfitresults_upsilon_%s.root",kineLabel.Data()));
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
