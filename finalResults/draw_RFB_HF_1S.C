#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"

void draw_RFB_HF_1S(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 3; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 4; 
  double xmin = 0;
  double xmax = 120;
//  double relsys = 0.1;

  double ex_range = 15/10.;
  double exsys_1[4] =  {ex_range,ex_range,ex_range,ex_range};
  double exsys_2[4] =  {ex_range,ex_range,ex_range,ex_range};
  double exsys_3[4] =  {ex_range,ex_range,ex_range,ex_range};
  double exsys_4[4] =  {ex_range,ex_range,ex_range,ex_range};

  double px_1[4] = {15-15/8*4,45-15/8*4,75-15/8*4,105-15/8*4};
  double px_2[4] = {15-15/8*1,45-15/8*1,75-15/8*1,105-15/8*1};
  double px_3[4] = {15+15/8*1,45+15/8*1,75+15/8*1,105+15/8*1};
  double px_4[4] = {15+15/8*4,45+15/8*4,75+15/8*4,105+15/8*4};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn;
  fIn = new TFile("Ups_1_RFB_HF.root","READ");
	TGraphErrors* gRPA_read[nState];
	TGraphErrors* gRPA[nState];
	TGraphErrors* gRPA_sys[nState];

  for (int is=0; is<nState; is++){
    gRPA[is]=(TGraphErrors*)fIn->Get(Form("gRFB_%d",is+1));
    gRPA_sys[is]=(TGraphErrors*)fIn->Get(Form("gRFB_%d",is+1));
    //cout << "gRPA["<<is<<"] = " <<gRPA[is]->GetPoint(ipt, pxtmp, pytmp) << endl;
  }
  //// read input file : syst.
  TFile* fInSys;
  fInSys = new TFile("../Systematics/mergedSys_ups1s_rfb.root","READ");
  TH1D* hSys[nState];
  int npoint[nState] = {4,4,4,4};
  for (int is=0; is<nState; is++){
    hSys[is] =(TH1D*)fInSys->Get(Form("hHFmerged%d",is+1));
//    npoint[is] = hSys[is]->GetSize()-2;
//    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
//    if (npoint[is] != gRPA[is]->GetN()) {cout << "Error!! data file and syst. file have dAifferent binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA[is]->GetPoint(ipt, pxtmp, pytmp); 
      relsys=hSys[is]->GetBinContent(ipt+1);
      extmp=gRPA[is]->GetErrorX(ipt);
      eytmp=gRPA[is]->GetErrorY(ipt);
      cout << ipt <<"th bin RFB value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      //// 1) remove ex from gRAA
      if(is==0){
        gRPA[is]->SetPoint(ipt,px_1[ipt],pytmp);
        gRPA_sys[is]->SetPoint(ipt,px_1[ipt],pytmp);
        gRPA[is]->SetPointError(ipt, 0, eytmp);
        gRPA_sys[is]->SetPointError(ipt, exsys_1[ipt], pytmp*relsys);
      }
      else if(is==1){
        gRPA[is]->SetPoint(ipt,px_2[ipt],pytmp);
        gRPA_sys[is]->SetPoint(ipt,px_2[ipt],pytmp);
        gRPA[is]->SetPointError(ipt, 0, eytmp);
        gRPA_sys[is]->SetPointError(ipt, exsys_2[ipt], pytmp*relsys);
      }
      else if(is==2){
        gRPA[is]->SetPoint(ipt,px_3[ipt],pytmp);
        gRPA_sys[is]->SetPoint(ipt,px_3[ipt],pytmp);
        gRPA[is]->SetPointError(ipt, 0, eytmp);
        gRPA_sys[is]->SetPointError(ipt, exsys_3[ipt], pytmp*relsys);
      }
      else if(is==3){
        gRPA[is]->SetPoint(ipt,px_4[ipt],pytmp);
        gRPA_sys[is]->SetPoint(ipt,px_4[ipt],pytmp);
        gRPA[is]->SetPointError(ipt, 0, eytmp);
        gRPA_sys[is]->SetPointError(ipt, exsys_4[ipt], pytmp*relsys);
      }
    }
  }
 
  // No upper limit
  ////////////////////////////////////////////////////////////////
/*  int ulstate = 2; //3S
  static const int n3s = 2;
  double boxw = 0.6; // for syst. box (vs cent)
  double lower68[n3s] = {lower68_pt1,lower68_pt2};
  double upper68[n3s] = {upper68_pt1,upper68_pt2};
  double lower95[n3s] = {lower95_pt1,lower95_pt2};
  double upper95[n3s] = {upper95_pt1,upper95_pt2};
//  if (n3s != npoint[ulstate]) {cout<<"ERROR!! # of bins for UL is wrong!!"<<endl;return;} 

  //// --- vs centrality
  TBox *box68per[n3s];
  TArrow *arr95per[n3s];
  for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
    //lower68=0; upper68=0; lower95=0; upper95=0; 
//    gRPA[ulstate]->GetPoint(ipt, pxtmp, pytmp);
    box68per[ipt] = new TBox(pxtmp-boxw,lower68[ipt],pxtmp+boxw,upper68[ipt]);
    arr95per[ipt] = new TArrow(pxtmp,lower95[ipt],pxtmp,upper95[ipt],0.027,"<-|"); //95%
    box68per[ipt]->SetLineColor(kGreen+3);
    box68per[ipt]->SetFillColorAlpha(kGreen-6,0.5);
    box68per[ipt]->SetLineWidth(1);
    arr95per[ipt]->SetLineColor(kGreen+2);
    arr95per[ipt]->SetLineWidth(2);
  }
// */

  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle(gRPA[is], is, is); 
    SetGraphStyleSys(gRPA_sys[is], is); 
	}
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  //// legend
  //// axis et. al
  gRPA_sys[0]->GetXaxis()->SetTitle("E_{T}^{HF |#eta|>4} [GeV]");
  gRPA_sys[0]->GetXaxis()->SetTitleOffset(1.1);
  gRPA_sys[0]->GetXaxis()->CenterTitle();
  gRPA_sys[0]->GetYaxis()->SetTitle("R_{FB}");
  gRPA_sys[0]->GetYaxis()->CenterTitle();
  gRPA_sys[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRPA_sys[0]->SetMinimum(0.0);
  gRPA_sys[0]->SetMaximum(1.8);
  gRPA_sys[0]->GetXaxis()->SetNdivisions(505);
  gRPA_sys[0]->GetXaxis()->SetBinLabel(1,"");
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nState; is++){
    if ( is==0) {gRPA_sys[is]->Draw("A5");}
/*    else if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++){
        box68per[ipt]->Draw("l");
      }
      gRPA_sys[is]->Draw("5");
    }
// */
    else {gRPA_sys[is]->Draw("5");}
  }
  
  for(int is=0;is<nState;is++){
/*    if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++) {
        arr95per[ipt]->Draw();
      }
      gRPA[is]->Draw("P");
    }
// */
//    else {gRPA[is]->Draw("P");}
    gRPA[is]->Draw("P");
  }
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg1= new TLegend(0.574, 0.19, 0.759, 0.366);
  SetLegendStyle(leg1);
  TLegend *leg= new TLegend(0.214, 0.19, 0.399, 0.366);
  SetLegendStyle(leg);
  TLegend *leg_up= new TLegend(0.57, 0.50, 0.78, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(16.,0.532,16.,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

    leg -> AddEntry(gRPA[0],"0 < |y_{CM}| < 0.4","lp");
    leg -> AddEntry(gRPA[1],"0.4 < |y_{CM}| < 0.8","lp");
    leg1 -> AddEntry(gRPA[2],"0.8 < |y_{CM}| < 1.2","lp");
    leg1 -> AddEntry(gRPA[3],"1.2 < |y_{CM}| < 1.93","lp");
    //TLegendEntry *ent=leg_up->AddEntry("ent"," #Upsilon(3S) 68\% CL","f");
    //ent->SetLineColor(kGreen+3);
    //ent->SetFillColorAlpha(kGreen-6,0.5);
    //ent->SetFillStyle(1001);
    //ent=leg_up->AddEntry("ent"," #Upsilon(3S) 95\% CL","f");
    //ent->SetLineColor(kWhite);
//    leg_up->SetTextSize(0.03);
    leg->Draw("same");
    leg1->Draw("same");
    //leg_up->Draw("same");


  //// draw text
  double sz_init = 0.925; double sz_step = 0.0675;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#varUpsilon} < 30 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "#varUpsilon(1S)");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "|y^{#mu#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  
  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.042);
  globtex_label->DrawLatex(0.223, sz_init-sz_step*11.56, "0-15");
  globtex_label->DrawLatex(0.409, sz_init-sz_step*11.56, "15-22");
  globtex_label->DrawLatex(0.618, sz_init-sz_step*11.56, "22-30");
  globtex_label->DrawLatex(0.805, sz_init-sz_step*11.56, "30-120");

  onSun(30,0,30,0.06,1,1);
  onSun(60,0,60,0.06,1,1);
  onSun(90,0,90,0.06,1,1);
  //Global Unc.
 
  double TAA_unc_Global_Hi = 0.068;
  double TAA_unc_Global_Lo = 0.072;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 2.4;
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox -> SetLineWidth(1);
  //globalUncBox -> Draw("l same");
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("plots/RFB_vs_hf_1S.pdf");
  c1->SaveAs("plots/RFB_vs_hf_1S.png");

	return;

} // end of main func.

