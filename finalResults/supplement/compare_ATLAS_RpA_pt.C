#include "../SONGKYO.h"
#include "../tdrstyle.C"
#include "../CMS_lumi_raaCent.C"
#include "../../cutsAndBin.h"

void compare_ATLAS_RpA_pt(int istate=1) //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 2; // 0: 15001, 1: ours
  double xmin = 0;
  double xmax = 46;
//  double relsys = 0.1;
  
  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3., 12.};

  //// ALICE values
  const int cn_1s =  8;
  double cpx_1s[cn_1s] =  {0.75,2.25,4,6,8.5,12,17,30};
  double cpy_1s[cn_1s] =  {0.671747,0.726022,0.792937,0.778067,0.846840,0.848327,1.006320,1.144240};
  double cex_1s[cn_1s] =  {0., 0., 0., 0.,0,0,0,0};
  double cey_1s[cn_1s] =  {0.0743494, 0.0650558, 0.0613383, 0.0613383, 0.0613383, 0.0613383, 0.0855019, 0.148699};
  double cexsys_1s[cn_1s] =  {0.75, 0.75, 1, 1, 1.5, 2, 3, 10};
  double ceysys_1s[cn_1s] = {0.115242, 0.0780669, 0.107807, 0.0762082, 0.070632, 0.0464684, 0.0613383, 0.0483271};
//  for (int it=0; it < cn_1s ; it++) {
//    ceysys_1s[it] = TMath::Sqrt(ceysys_1s_1[it]*ceysys_1s_1[it]+ceysys_1s_2[it]*ceysys_1s_2[it]);
//  }

  const int cn_2s =  2;
  double cpx_2s[cn_2s] =  {0.6, 1.8};
  double cpy_2s[cn_2s] =  {0.116, 0.138};
  double cex_2s[cn_2s] =  {0., 0.};
  double cey_2s[cn_2s] =  {0.034, 0.049};
  double cexsys_2s[cn_2s] =  {0.6, 0.6};
  double ceysys_2s_1[cn_2s] =  {0.036, 0.073};
  double ceysys_2s_2[cn_2s] =  {0.002, 0.002};
  double ceysys_2s[cn_2s] = {0.039, 0.064};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphErrors* gRAA_sys[nfile];
  //// ALICE
  if (istate==1) {
    gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
    gRAA_sys[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, ceysys_1s); 
  }
  else {
    gRAA[0] = new TGraphErrors(cn_2s, cpx_2s, cpy_2s, cex_2s, cey_2s); 
    gRAA_sys[0] = new TGraphErrors(cn_2s, cpx_2s, cpy_2s, cexsys_2s, ceysys_2s); 
  } 
  //// 2) ours
  TFile* fIn = new TFile(Form("../Ups_%d_1D.root",istate),"READ");
  gRAA[1]=(TGraphErrors*)fIn->Get("gRPA_pt");
  gRAA_sys[1]=(TGraphErrors*)fIn->Get("gRPA_pt");
  //// read input file : syst. 
  TFile* fInSys = new TFile(Form("../../Systematics/mergedSys_ups%ds.root",istate),"READ");
  TH1D* hSys = (TH1D*)fInSys->Get("hptRPA_merged");
  int npoint = hSys->GetSize()-2;
  cout << "*** Y("<<istate<<") : # of point = " << npoint << endl;
  
  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  if (npoint != gRAA[1]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
    gRAA[1]->GetPoint(ipt, pxtmp, pytmp);
    extmp=gRAA[1]->GetErrorX(ipt);
    eytmp=gRAA[1]->GetErrorY(ipt);
    relsys=hSys->GetBinContent(ipt+1);
    // 1) remove ex from gRAA
    gRAA[1]->SetPointError(ipt, 0, eytmp);
    // 2) set ey for gRAA_sys (assign 10% temporarily)
    //gRAA_sys[1]->SetPointError(ipt, extmp, pytmp*relsys);
    if (istate==1) gRAA_sys[1]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
    else gRAA_sys[1]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
  }
  
  ////////////////////////////////////////////////////////////////
  
  //remove 1st point
  
  //// graph style 
  SetGraphStyle(gRAA[0], 1, 1); 
  SetGraphStyleSys(gRAA_sys[0], 1); 
  SetGraphStyle(gRAA[1], 0, 0); 
  SetGraphStyleSys(gRAA_sys[1], 0); 
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  //// legend
  TLegend *leg= new TLegend(0.36, 0.28, 0.87, 0.44);
  SetLegendStyle(leg);
  leg -> SetTextSize(0.038);
  leg -> AddEntry(gRAA[1],Form("#varUpsilon(%dS), |y^{#varUpsilon}_{CM}| < 1.93",istate),"lp");
  leg -> AddEntry(gRAA[0],Form("ATLAS #varUpsilon(%dS), -2 < y^{#varUpsilon}_{CM} < 1.5 ",istate),"lp");
  
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetTitleOffset(1.0);
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRAA_sys[0]->GetXaxis()->SetNdivisions(505);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.5);
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
  double sys_global_x = 1.4;
  //double sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_aa*lumi_unc_aa);
  double sys_global_y_15001 = 0.0892193;
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kRed-2);
  globalUncBox -> SetFillColorAlpha(kPink-6,0.4);
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");
  
  TBox *globalUncBox_15001 = new TBox(xmax-sys_global_x*2,1-sys_global_y_15001,xmax-sys_global_x-0.001,1+sys_global_y_15001);
  globalUncBox_15001 -> SetLineColor(kBlue-3);
  globalUncBox_15001 -> SetFillColorAlpha(kBlue-3,0.4);
  globalUncBox_15001 -> SetLineWidth(1);
  globalUncBox_15001 -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("../plots/%dS_comp_RpA_vs_pt_ATLAS.pdf",istate));
  c1->SaveAs(Form("../plots/%dS_comp_RpA_vs_pt_ATLAS.png",istate));

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

