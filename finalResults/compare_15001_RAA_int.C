#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"

void compare_15001_RAA_int() //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 101; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 5; // 0: 15001, 1: ours
  double xmax = 2.85;
//  double relsys = 0.1;
  
  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3.,12.};

  double exsys = 0.05;
  double exsys_align = 0.075;

  //// 15001 values
  const int cn_1s =  3;
  double cpx_1s[cn_1s] =  {0.51-exsys_align, 1.425-exsys_align, 2.3-exsys_align};
  double cpx_1s_exsys[cn_1s] = {0.51+exsys_align, 1.425+exsys_align, 2.3+exsys_align};
  double cpy_1s[cn_1s] =  {0.453, 0.119, 0.145}; 
  double cex_1s[cn_1s] =  {0., 0., 0};
  double cey_1s[cn_1s] =  {0.014, 0.028, 0.031};
  double cexsys_1s[cn_1s] =  {exsys, exsys, exsys};
  double ceysys_1s[cn_1s] =  {0.046, 0.015, 0.058};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphAsymmErrors* gRAA_sys[nfile];
  //// 1) 15001
  gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
  gRAA_sys[0] = new TGraphAsymmErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, cexsys_1s, ceysys_1s, ceysys_1s); 
  //// 2) ours
  TFile* fIn_1S = new TFile("Ups_1_RAA.root","READ");
  TFile* fIn_2S = new TFile("Ups_2_RAA.root","READ");
  TFile* fIn_3S = new TFile("Ups_3_RAA.root","READ");
  gRAA[1]=(TGraphErrors*)fIn_1S->Get("gRAA_int");
  gRAA_sys[1]= new TGraphAsymmErrors();
  gRAA[2]=(TGraphErrors*)fIn_2S->Get("gRAA_int");
  gRAA_sys[2]= new TGraphAsymmErrors();
  gRAA[3]=(TGraphErrors*)fIn_3S->Get("gRAA_int");
  gRAA_sys[3]= new TGraphAsymmErrors();
  gRAA[4]= new TGraphErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cey_1s);
  gRAA_sys[4]= new TGraphAsymmErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cex_1s, cey_1s, cey_1s);
  //// read input file : syst.
  TFile* fInSys_1S_Hi = new TFile("../Systematic/mergedSys_ups1s_asymHi.root","READ");
  TFile* fInSys_2S_Hi = new TFile("../Systematic/mergedSys_ups2s_asymHi.root","READ");
  TFile* fInSys_3S_Hi = new TFile("../Systematic/mergedSys_ups3s_asymHi.root","READ");
  TFile* fInSys_1S_Lo = new TFile("../Systematic/mergedSys_ups1s_asymLo.root","READ");
  TFile* fInSys_2S_Lo = new TFile("../Systematic/mergedSys_ups2s_asymLo.root","READ");
  TFile* fInSys_3S_Lo = new TFile("../Systematic/mergedSys_ups3s_asymLo.root","READ");
  TH1D* hSys_1S_Hi = (TH1D*)fInSys_1S_Hi->Get("hintRAA_merged");
  TH1D* hSys_2S_Hi = (TH1D*)fInSys_2S_Hi->Get("hintRAA_merged");
  TH1D* hSys_3S_Hi = (TH1D*)fInSys_3S_Hi->Get("hintRAA_merged");
  TH1D* hSys_1S_Lo = (TH1D*)fInSys_1S_Lo->Get("hintRAA_merged");
  TH1D* hSys_2S_Lo = (TH1D*)fInSys_2S_Lo->Get("hintRAA_merged");
  TH1D* hSys_3S_Lo = (TH1D*)fInSys_3S_Lo->Get("hintRAA_merged");
  
  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Hi, relsys_Lo;
  npoint = gRAA[0]->GetN();
  if (npoint !=cn_1s) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) 
  {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Hi=0; relsys_Lo;
    gRAA[ipt+1]->GetPoint(0, pxtmp, pytmp);
    extmp=gRAA[ipt+1]->GetErrorX(0);
    eytmp=gRAA[ipt+1]->GetErrorY(0);
    relsys_Hi=hSys_1S_Hi->GetBinContent(1);
    relsys_Hi=TMath::Sqrt(relsys_Hi*relsys_Hi + lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    relsys_Lo=hSys_1S_Lo->GetBinContent(1);
    relsys_Lo=TMath::Sqrt(relsys_Lo*relsys_Lo + lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    gRAA[4]->SetPoint(ipt, cpx_1s_exsys[ipt], pytmp);
    gRAA_sys[4]->SetPoint(ipt, cpx_1s_exsys[ipt], pytmp);
    gRAA[4]->SetPointError(ipt, 0, eytmp);
    gRAA_sys[4]->SetPointError(ipt, exsys, exsys, pytmp*relsys_Lo, pytmp*relsys_Hi);
  }
 

  gRAA[0]->SetPoint(2,-100,-100); 
  gRAA_sys[0]->SetPoint(2,-100,-100); 
  gRAA[4]->SetPoint(2,-100,-100); 
  gRAA_sys[4]->SetPoint(2,-100,-100); 
  ////////////////////////////////////////////////////////////////


  //******************************************
  // Upper Limit
  //******************************************
  
  const int numComp = 2;
  double lower95_int[numComp] = {0,lower95_cint}; 
  double upper95_int[numComp] = {0.145,upper95_cint};
  
  double align_upper = 0.558;
  double exsys_upper = 0.15;

   TArrow *arr95per_int[numComp];
   for(int icomp = 0; icomp<numComp; icomp++)
   {
     arr95per_int[icomp] = new TArrow((cpx_1s[cn_1s-1]+cpx_1s[cn_1s-2])/2-exsys/2+exsys_upper*icomp+align_upper,lower95_int[icomp],(cpx_1s[cn_1s-1]+cpx_1s[cn_1s-2])/2-exsys/2+exsys_upper*icomp+align_upper,upper95_int[icomp],0.027,"<-|");
     arr95per_int[icomp]->SetLineWidth(2);
     if(icomp==0) arr95per_int[icomp]->SetLineColor(kBlue-3);
     else if(icomp!=0) arr95per_int[icomp]->SetLineColor(kPink-6);
   }


  //// graph style 
  SetGraphStyle(gRAA[0], 1, 1); 
  SetGraphStyleSys(gRAA_sys[0], 1); 
  SetGraphStyle(gRAA[4], 0, 0); 
  SetGraphStyleSys(gRAA_sys[4], 0); 
  
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
  TLegend *leg= new TLegend(0.574, 0.57, 0.974, 0.74);
  SetLegendStyle(leg);
  leg -> SetHeader("");
  //leg -> SetHeader("#Upsilon's");
  leg -> AddEntry(gRAA[0],"#sqrt{s_{NN}} = 2.76 TeV","lp");
  leg -> AddEntry(gRAA[4],"#sqrt{s_{NN}} = 5.02 TeV","lp");

  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.046);
  header->SetTextFont(62);
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.);
 
  for(int i=0;i<=4;i++)
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
  gRAA_sys[4]->Draw("5");
  gRAA[4]->Draw("P");
  arr95per_int[0]->Draw();
  arr95per_int[1]->Draw();
  dashedLine(0.,1.,xmax,1.,1,1);
  leg->Draw();

  //// draw text
  double sz_init = 0.87; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu#mu} < 30 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step-0.007, "|y^{#mu#mu}| < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");
  globtex_label->DrawLatex(0.243, sz_init-sz_step*15.24, "#Upsilon(1S)");
  globtex_label->DrawLatex(0.505, sz_init-sz_step*15.24, "#Upsilon(2S)");
  globtex_label->DrawLatex(0.782, sz_init-sz_step*15.24, "#Upsilon(3S)");
  
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
  c1->SaveAs("comp15001_RAA_int_asym.pdf");
  c1->SaveAs("comp15001_RAA_int_asym.png");

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

