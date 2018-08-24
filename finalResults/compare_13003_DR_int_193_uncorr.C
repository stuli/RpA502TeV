#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"

void compare_13003_DR_int_193_uncorr() //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
 
  TH1::SetDefaultSumw2(); 
  const int nfile = 5; // 0: 15001, 1: ours
  double xmax = 2.85;
//  double relsys = 0.1;
  
  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3.,12.};

  double exsys = 0.05;
  double exsys_align = 0.075;

  double x_loc1 = 0.65;
  double x_loc2 = 2.08;
    
  //// 15001 values
  const int cn_1s =  2;
  double cpx_1s[cn_1s] =  {x_loc1-exsys_align, x_loc2-exsys_align};
  double cpx_1s_exsys[cn_1s] = {x_loc1+exsys_align, x_loc2+exsys_align};
  double cpx_1s_exsys_193[cn_1s] = {x_loc1+exsys_align*3, x_loc2+exsys_align*3};
  double cpy_1s[cn_1s] =  {0.83*0.779/0.796,0.71*0.625/0.631}; 
  double cex_1s[cn_1s] =  {0., 0};
  double cey_1s[cn_1s] =  {0.05, 0.08};
  double cexsys_1s[cn_1s] =  {exsys, exsys};
  double ceysys_1s[cn_1s] =  {0.05,0.09};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphAsymmErrors* gRAA_sys[nfile];
  //// 1) 15001
  gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
  gRAA_sys[0] = new TGraphAsymmErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, cexsys_1s, ceysys_1s, ceysys_1s); 
  //// 2) ours
  TFile* fIn_pp = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PP_DATA_pt0.0-30.0_y0.00-1.93_muPt4.0.root","READ");
  TFile* fIn_pPb = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu2.4/FitResults/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PA_DATA_pt0.0-30.0_y-1.93-1.93_muPt4.0.root","READ");
  TFile* fIn_pp_193 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu1.93/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PP_DATA_pt0.0-30.0_y-1.93-1.93_muPt4.0.root","READ");
  TFile* fIn_pPb_193 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_RpA/UpsilonpPb5TeV/RpA5.02TeV/Fitting/AllParmFree_SingleMu1.93/AllParmFree_fitresults_upsilon_DoubleCB_5TeV_PA_DATA_pt0.0-30.0_y-1.93-1.93_muPt4.0.root","READ");
  TH1D* h_pp = (TH1D*) fIn_pp -> Get("fitResults");
  TH1D* h_pPb = (TH1D*) fIn_pPb -> Get("fitResults");
  TH1D* h_pp_193 = (TH1D*) fIn_pp_193 -> Get("fitResults");
  TH1D* h_pPb_193 = (TH1D*) fIn_pPb_193 -> Get("fitResults");

  TH1D* h_pPb_div = (TH1D*) h_pPb -> Clone("fitResults_clone");
  TH1D* h_pp_div = (TH1D*) h_pp -> Clone("fitResults_clone");
  TH1D* h_pPb_div_193 = (TH1D*) h_pPb_193 -> Clone("fitResults_clone");
  TH1D* h_pp_div_193 = (TH1D*) h_pp_193 -> Clone("fitResults_clone");
  for(int i=0;i<3;i++)
  {
    h_pPb_div -> SetBinContent(i+1,h_pPb->GetBinContent(1));
    h_pPb_div -> SetBinError(i+1,h_pPb->GetBinError(1));
    h_pp_div -> SetBinContent(i+1,h_pp->GetBinContent(1));
    h_pp_div -> SetBinError(i+1,h_pp->GetBinError(1));
    h_pPb_div_193 -> SetBinContent(i+1,h_pPb_193->GetBinContent(1));
    h_pPb_div_193 -> SetBinError(i+1,h_pPb_193->GetBinError(1));
    h_pp_div_193 -> SetBinContent(i+1,h_pp_193->GetBinContent(1));
    h_pp_div_193 -> SetBinError(i+1,h_pp_193->GetBinError(1));
  }

  h_pPb->Divide(h_pPb_div);
  h_pp->Divide(h_pp_div);
  h_pPb->Divide(h_pp);

  h_pPb_193->Divide(h_pPb_div_193);
  h_pp_193->Divide(h_pp_div_193);
  h_pPb_193->Divide(h_pp_193);


  gRAA[1]= new TGraphErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cey_1s);
  gRAA_sys[1]= new TGraphAsymmErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cex_1s, cey_1s, cey_1s); 
  gRAA[2]= new TGraphErrors(cn_1s, cpx_1s_exsys_193, cpy_1s, cex_1s, cey_1s);
  gRAA_sys[2]= new TGraphAsymmErrors(cn_1s, cpx_1s_exsys_193, cpy_1s, cex_1s, cex_1s, cey_1s, cey_1s); 

  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Hi, relsys_Lo;
  double pxtmp_193, pytmp_193, extmp_193, eytmp_193;
  double relsys_Hi_193;
  int npoint = gRAA[0]->GetN();
  if (npoint !=cn_1s) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) 
  {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Hi=0; relsys_Lo=0;
    pxtmp = cpx_1s_exsys[ipt];
    pxtmp_193 = cpx_1s_exsys_193[ipt];
    extmp=exsys;
    extmp_193=exsys;
    relsys_Hi=0.05;
    relsys_Hi_193=0.05;
    pytmp = h_pPb->GetBinContent(ipt+2);
    pytmp_193 = h_pPb_193->GetBinContent(ipt+2);
    eytmp = h_pPb->GetBinError(ipt+2);
    eytmp_193 = h_pPb_193->GetBinError(ipt+2);
    gRAA[1]->SetPoint(ipt, pxtmp, pytmp);
    gRAA_sys[1]->SetPoint(ipt, pxtmp, pytmp);
    gRAA[1]->SetPointError(ipt, 0, eytmp);
    gRAA_sys[1]->SetPointError(ipt, extmp, extmp, pytmp*relsys_Hi, pytmp*relsys_Hi);

    gRAA[2]->SetPoint(ipt, pxtmp_193, pytmp_193);
    gRAA_sys[2]->SetPoint(ipt, pxtmp_193, pytmp_193);
    gRAA[2]->SetPointError(ipt, 0, eytmp_193);
    gRAA_sys[2]->SetPointError(ipt, extmp_193, extmp_193, pytmp_193*relsys_Hi_193, pytmp_193*relsys_Hi_193);
  }
 
  //// graph style 
  SetGraphStyle(gRAA[0], 1, 1); 
  SetGraphStyleSys(gRAA_sys[0], 1); 
  SetGraphStyle(gRAA[1], 0, 0); 
  SetGraphStyleSys(gRAA_sys[1], 0); 
  SetGraphStyle(gRAA[2], 3, 3); 
  SetGraphStyleSys(gRAA_sys[2], 3); 
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.047);
  
  //// legend
  TLegend *leg= new TLegend(0.511, 0.63, 0.801, 0.80);
  SetLegendStyle(leg);
  leg -> SetHeader("");
  //leg -> SetHeader("#Upsilon's");
  leg -> AddEntry(gRAA[0],"HIN-13-003  |#eta^{#mu}| < 1.93","lp");
  leg -> AddEntry(gRAA[1],"New Result  |#eta^{#mu}| < 2.4","lp");
  leg -> AddEntry(gRAA[2],"New Result  |#eta^{#mu}| < 1.93","lp");

  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.046);
  header->SetTextFont(62);
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("[#Upsilon(nS)/#Upsilon(1S)]_{pPb}/[#Upsilon(nS)/#Upsilon(1S)]_{pp}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.6);
 
  for(int i=0;i<=2;i++)
  {
    gRAA[i]->GetXaxis()->SetBinLabel(10,"");
    gRAA_sys[i]->GetXaxis()->SetBinLabel(10,"");
  }
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.1);
  gPad->SetTopMargin(0.067);
  gRAA_sys[0]->Draw("A5");
  gRAA[0]->Draw("P");
  gRAA_sys[1]->Draw("5");
  gRAA[1]->Draw("P");
  gRAA_sys[2]->Draw("5");
  gRAA[2]->Draw("P");
  dashedLine(0.,1.,xmax,1.,1,1);
  leg->Draw();

  //// draw text
  double sz_init = 0.87; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu#mu} < 30 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step-0.007, "|y^{#mu#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2-0.007, "|#eta^{#mu}| < 2.4");
  globtex_label->DrawLatex(0.248, sz_init-sz_step*15.24, "#Upsilon(2S)/#Upsilon(1S)");
  globtex_label->DrawLatex(0.652, sz_init-sz_step*15.24, "#Upsilon(3S)/#Upsilon(1S)");
  
  double TAA_unc_Global_Hi = 0.068;
  double TAA_unc_Global_Lo = 0.072;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Hi*TAA_unc_Global_Hi+nMB_unc*nMB_unc);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Lo*TAA_unc_Global_Lo+nMB_unc*nMB_unc);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.8;
  //double sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_aa*lumi_unc_aa);
  double sys_global_y_15001 = 0.079;
  TBox *globalUncBox = new TBox(xmax-sys_global_x*2,1-sys_global_y_Lo,xmax-sys_global_x-0.05,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kRed-2);
  globalUncBox -> SetFillColorAlpha(kPink-6,0.6);
  globalUncBox -> SetLineWidth(1);
  //globalUncBox -> Draw("l same");
  
  TBox *globalUncBox_15001 = new TBox(xmax-sys_global_x,1-sys_global_y_15001,xmax,1+sys_global_y_15001);
  globalUncBox_15001 -> SetLineColor(kBlue-3);
  globalUncBox_15001 -> SetFillColorAlpha(kBlue-3,0.6);
  globalUncBox_15001 -> SetLineWidth(1);
  //globalUncBox_15001 -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("Comp13003_DR_int_asym_193_uncorr.pdf");
  c1->SaveAs("Comp13003_DR_int_asym_193_uncorr.png");

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

