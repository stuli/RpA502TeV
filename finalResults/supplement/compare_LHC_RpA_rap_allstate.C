#include "../SONGKYO.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_raaCent.C"
#include "../../cutsAndBin.h"

void compare_LHC_RpA_rap_allstate(int istate=1) //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 3; // 0: ALICE, 1: ATLAS, 2: LHCb
  double xmin = -4.7;
  double xmax = 5.3;
//  double relsys = 0.1;
  
  double exsys_1s[8] =  {0.73/2, 0.2, 0.2, 0.2, 0.2,0.2,0.2, 0.73/2};
  double exsys_2s[4] =  {1.13/2, 0.4,0.4,1.13/2};
  double exsys_3s[2] =  {1.93/2, 1.93/2};

  //// ALICE values
  const int cn_alice =  4;
  double cpx_alice[cn_alice] =  {-3.995,-3.245,2.495,3.245};
  double cpy_alice[cn_alice] =  {0.78,0.91,0.66,0.86};
  double cex_alice[cn_alice] =  {0., 0., 0., 0.};
  double cey_alice[cn_alice] =  {0.15,0.15,0.09,0.17};
  double cexsys_alice[cn_alice] =  {0.465,0.285,0.465,0.285};
  double ceysys_alice[cn_alice] = {0.13,0.15,0.08,0.13};

  //// ATLAS values
  const int cn_atlas =  5;
  double cpx_atlas[cn_atlas] =  {-1.75,-1.125,-0.375,0.375,1.125};
  double cpy_atlas[cn_atlas] =  {0.743866, 0.710037, 0.704833, 0.830855, 0.876952};
  double cex_atlas[cn_atlas] =  {0., 0., 0., 0., 0.};
  double cey_atlas[cn_atlas] =  {0.070632, 0.0520446, 0.0501859, 0.0576208, 0.063197};
  double cexsys_atlas[cn_atlas] =  {0.25,0.375,0.375,0.375,0.375};
  double ceysys_atlas[cn_atlas] = {0.0910781, 0.0576208, 0.0743494, 0.0966543, 0.063197};

  //// LHCb values
  const int cn_lhcb =  2;
  double cpx_lhcb[cn_lhcb] =  {-3.25,3.25};
  double cpy_lhcb[cn_lhcb] =  {1.21, 0.9};
  double cex_lhcb[cn_lhcb] =  {0., 0.};
  double cey_lhcb[cn_alice] =  {0.23,0.1};
  double cexsys_lhcb[cn_lhcb] =  {0.75,0.75};
  double ceysys_lhcb[cn_lhcb] = {0.1,0.08};


  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphErrors* gRAA_sys[nfile];
	TGraphErrors* gRAAn[nfile];
	TGraphErrors* gRAA_sysn[nfile];
  //// LHC exp
    gRAA[1] = new TGraphErrors(cn_alice, cpx_alice, cpy_alice, cex_alice, cey_alice); 
    gRAA_sys[1] = new TGraphErrors(cn_alice, cpx_alice, cpy_alice, cexsys_alice, ceysys_alice); 
    gRAA[0] = new TGraphErrors(cn_atlas, cpx_atlas, cpy_atlas, cex_atlas, cey_atlas); 
    gRAA_sys[0] = new TGraphErrors(cn_atlas, cpx_atlas, cpy_atlas, cexsys_atlas, ceysys_atlas); 
    gRAA[2] = new TGraphErrors(cn_lhcb, cpx_lhcb, cpy_lhcb, cex_lhcb, cey_lhcb); 
    gRAA_sys[2] = new TGraphErrors(cn_lhcb, cpx_lhcb, cpy_lhcb, cexsys_lhcb, ceysys_lhcb); 

  //// 2) ours
  TFile* fIn[3];
  TFile* fInSys[3];
  TH1D* hSys[3];
  int npoint[3];
  for(int i=0;i<3;i++){
    fIn[i] = new TFile(Form("../Ups_%d_1D.root",i+1),"READ");
    gRAAn[i]=(TGraphErrors*)fIn[i]->Get("gRPA_rap");
    gRAA_sysn[i]=(TGraphErrors*)fIn[i]->Get("gRPA_rap");
    fInSys[i] = new TFile(Form("../../Systematics/mergedSys_ups%ds.root",i+1),"READ");
    hSys[i] = (TH1D*)fInSys[i]->Get("hrapRPA_merged");
    npoint[i] = hSys[i]->GetSize()-2;
    cout << "*** Y("<<i+1<<") : # of point = " << npoint[i] << endl;
  }
  //// read input file : syst. 

  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  for(int i=0;i<3;i++){
    if (npoint[i] != gRAAn[i]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[i]; ipt++) {
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRAAn[i]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAAn[i]->GetErrorX(ipt);
      eytmp=gRAAn[i]->GetErrorY(ipt);
      relsys=hSys[i]->GetBinContent(ipt+1);
      // 1) remove ex from gRAA
      gRAAn[i]->SetPointError(ipt, 0, eytmp);
      // 2) set ey for gRAA_sys (assign 10% temporarily)
      //gRAA_sys[1]->SetPointError(ipt, extmp, pytmp*relsys);
      if(i==0) gRAA_sysn[i]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if(i==1) gRAA_sysn[i]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else if(i==2) gRAA_sysn[i]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
    }
  } 
  ////////////////////////////////////////////////////////////////
  
  //remove 1st point

  //// graph style
  for(int i=0;i<3;i++){ 
    SetGraphStyle_comp(gRAA[i], i, i); 
    SetGraphStyleSys_comp(gRAA_sys[i], i); 
    SetGraphStyle(gRAAn[i], i, i); 
    SetGraphStyleSys(gRAA_sysn[i], i); 
  }
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.03);
  
  //// legend
  TLegend *leg= new TLegend(0.20, 0.71, 0.56, 0.88);
  SetLegendStyle(leg);
  leg -> SetTextSize(0.033);
  leg -> AddEntry(gRAA[1],Form("ALICE #varUpsilon(%dS), p^{#varUpsilon}_{T} > 0 GeV/c",istate),"lp");
  leg -> AddEntry(gRAA[0],Form("ATLAS #varUpsilon(%dS), p^{#varUpsilon}_{T} < 40 GeV/c",istate),"lp");
  leg -> AddEntry(gRAA[2],Form("LHCb #varUpsilon(%dS), p^{#varUpsilon}_{T} < 15 GeV/c",istate),"lp");
  
  TLegend *leg1= new TLegend(0.75, 0.62, 0.96, 0.78);
  SetLegendStyle(leg1);
  leg1 -> SetTextSize(0.033);
  for(int i=0;i<3;i++){
  leg1 -> AddEntry(gRAAn[i],Form("#varUpsilon(%dS)",i+1),"lp");
  }
  
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("y_{CM}^{#varUpsilon}");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetTitleOffset(1.0);
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRAA_sys[0]->GetXaxis()->SetNdivisions(505);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(2.1);
  /// for rap
 
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nfile; is++){
    if ( is==0) { gRAA_sys[is]->Draw("A5"); gRAA_sysn[is]->Draw("5");}
    else { gRAA_sys[is]->Draw("5");  gRAA_sysn[is]->Draw("5");}
    gRAA[is]->Draw("P");
    gRAAn[is]->Draw("P");
	}
  dashedLine(xmin,1.,xmax,1.,1,1);
  leg->Draw("same");
  leg1->Draw("same");
//  leg_comp->Draw("same");

  //// draw text
  double sz_init = 0.593; double sz_step = 0.0535;
  globtex->DrawLatex(0.75, sz_init, "p_{T}^{#varUpsilon} < 30 GeV/c");
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
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.1);
  globalUncBox -> SetLineWidth(2);
  globalUncBox -> Draw("l same");
  
  TBox *globalUncBox_alice = new TBox(xmax-sys_global_x*2,1-sys_global_y_alice,xmax-sys_global_x-0.003,1+sys_global_y_alice);
  globalUncBox_alice -> SetLineColor(kOrange+3);
  globalUncBox_alice -> SetFillColorAlpha(kOrange-7,0.6);
  globalUncBox_alice -> SetLineWidth(2);
  globalUncBox_alice -> Draw("l same");
  
  double sys_global_y_atlas = 0.0892193;
  TBox *globalUncBox_atlas = new TBox(xmax-sys_global_x*3,1-sys_global_y_atlas,xmax-sys_global_x*2-0.005,1+sys_global_y_atlas);
  globalUncBox_atlas -> SetLineColor(kViolet+3);
  globalUncBox_atlas -> SetFillColorAlpha(kViolet-5,0.4);
  globalUncBox_atlas -> SetLineWidth(2);
  globalUncBox_atlas -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("../plots/%dS_comp_RpA_vs_rap_LHC_all.pdf",istate));
  c1->SaveAs(Form("../plots/%dS_comp_RpA_vs_rap_LHC_all.png",istate));

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

