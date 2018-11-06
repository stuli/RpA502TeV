#include "JaebeomStyle.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
void draw_RpA_1D(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = false;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmin_ = 0;
  double xmax_ = 30;
  double xmin = -1.93;
  double xmax = 1.93;
//  double relsys = 0.1;

  double exsys_1s_[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s_[3] =  {2., 2.5, 10.5};
  double exsys_3s_[2] =  {3., 12.};

  double exsys_1s[8] =  {0.73/2, 0.2, 0.2, 0.2, 0.2,0.2,0.2, 0.73/2};
  double exsys_2s[4] =  {1.13/2, 0.4,0.4,1.13/2};
  double exsys_3s[2] =  {1.93/2, 1.93/2};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRPA_l[nState];
	TGraphErrors* gRPA_h[nState];
	TGraphErrors* gRPA_d[nState];
	TGraphErrors* gRPA_sys_l[nState];
	TGraphErrors* gRPA_sys_h[nState];
	TGraphErrors* gRPA_sys_d[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_1D.root",is+1),"READ");
    gRPA_l[is]=(TGraphErrors*)fIn[is]->Get("gRPA_rap");
    gRPA_h[is]=(TGraphErrors*)fIn[is]->Get("gRPA_pt");
    gRPA_sys_l[is]=(TGraphErrors*)fIn[is]->Get("gRPA_rap");
    gRPA_sys_h[is]=(TGraphErrors*)fIn[is]->Get("gRPA_pt");
  }

  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys_h[nState];
  TH1D* hSys_l[nState];
  int npoint_l[nState];
  int npoint_h[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematics/mergedSys_ups%ds.root",is+1),"READ");
    hSys_l[is]=(TH1D*)fInSys[is]->Get("hrapRPA_merged");
    hSys_h[is]=(TH1D*)fInSys[is]->Get("hptRPA_merged");
    npoint_l[is] = hSys_l[is]->GetSize()-2;
    npoint_h[is] = hSys_h[is]->GetSize()-2;
    cout << "*** Y("<<is+1<<") rap : # of point = " << npoint_l[is] << endl;
    cout << "*** Y("<<is+1<<") pt : # of point = " << npoint_h[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint_l[is] != gRPA_l[is]->GetN()) {cout << "Error!! data file and syst. file have dAifferent binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint_l[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_l[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_l[is]->GetErrorX(ipt);
      eytmp=gRPA_l[is]->GetErrorY(ipt);
      relsys=hSys_l[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin low pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      gRPA_l[is]->SetPointError(ipt, 0, eytmp);
      if (is==0) gRPA_sys_l[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys_l[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRPA_sys_l[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
    }
    for (int ipt=0; ipt<npoint_h[is];ipt++){
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_h[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_h[is]->GetErrorX(ipt);
      eytmp=gRPA_h[is]->GetErrorY(ipt);
      relsys=hSys_h[is]->GetBinContent(ipt+1);
      cout << endl;
      cout << ipt <<"th bin high pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      gRPA_h[is]->SetPointError(ipt, 0, eytmp);
      if (is==0) gRPA_sys_h[is]->SetPointError(ipt, exsys_1s_[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys_h[is]->SetPointError(ipt, exsys_2s_[ipt], pytmp*relsys);
      else gRPA_sys_h[is]->SetPointError(ipt, exsys_3s_[ipt], pytmp*relsys);
    }
  }
 
  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle2(gRPA_l[is], is, is); 
    SetGraphStyleSys2(gRPA_sys_l[is], is); 
    SetGraphStyle2(gRPA_h[is], is, is); 
    SetGraphStyleSys2(gRPA_sys_h[is], is); 
	}
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
//  globtex->SetTextAlign(12); //left-center
//  globtex->SetTextAlign(11); //left-bottom
  globtex->SetTextAlign(31); //right-bottom
  globtex->SetTextFont(42);
//  globtex->SetTextSize(0.035);
  globtex->SetTextSize(0.040);
  
  //// legend
  //// axis et. al
  gRPA_sys_l[0]->GetXaxis()->SetTitle("y_{CM}^{#Upsilon}");
  gRPA_sys_l[0]->GetXaxis()->CenterTitle();
  gRPA_sys_l[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRPA_sys_l[0]->GetYaxis()->CenterTitle();
  gRPA_sys_l[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRPA_sys_l[0]->SetMinimum(0.0);
  gRPA_sys_l[0]->SetMaximum(1.7);
  gRPA_sys_l[0]->GetXaxis()->SetNdivisions(505);
  gRPA_sys_h[0]->GetXaxis()->SetTitle("p_{T}^{#Upsilon} (GeV/c)");
  gRPA_sys_h[0]->GetXaxis()->CenterTitle();
  gRPA_sys_h[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRPA_sys_h[0]->GetYaxis()->CenterTitle();
  gRPA_sys_h[0]->GetXaxis()->SetLimits(xmin_,xmax_);
  gRPA_sys_h[0]->SetMinimum(0.0);
  gRPA_sys_h[0]->SetMaximum(1.7);
  gRPA_sys_h[0]->GetXaxis()->SetNdivisions(505);
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nState; is++){
    if ( is==0) {gRPA_sys_l[is]->Draw("A5");}
    else {gRPA_sys_l[is]->Draw("5");}
  }
  
  for(int is=0;is<nState;is++){
    gRPA_l[is]->Draw("P");
  }
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.22, 0.68, 0.465, 0.876);
  SetLegendStyle(leg);
//  leg->SetTextSize(0.036);
  leg->SetTextSize(0.042);
  leg->SetTextFont(22);
  TLegend *leg_up= new TLegend(0.57, 0.50, 0.78, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(16.,0.532,16.,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
//      leg -> AddEntry(gRPA_l[is],Form(" #varUpsilon(%dS)",is+1),"lp");
      leg -> AddEntry(gRPA_l[is],Form(" #Upsilon(%dS)",is+1),"lp");
      leg -> Draw("same");
    }
  }


  //// draw text
//  double sz_init = 0.925; double sz_step = 0.1975;
  double sz_init = 1.1; double sz_step = 0.3;

//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
//  globtex->DrawLatex(0.717, sz_init-sz_step, "p_{T}^{#Upsilon} < 30 GeV/c");
  globtex->DrawLatex(0.92, sz_init-sz_step, "p_{T}^{#Upsilon} < 30 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  

  //Global Unc.
 
  double TAA_unc_Global_Hi = 0.068;
  double TAA_unc_Global_Lo = 0.072;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = .2;
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");

  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("plots/RpA_vs_rap_1D.pdf");
  c1->SaveAs("plots/RpA_vs_rap_1D.png");


  //// draw  
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nState; is++){
    if ( is==0) {gRPA_sys_h[is]->Draw("A5");}
    else {gRPA_sys_h[is]->Draw("5");}
  }
  
  for(int is=0;is<nState;is++){
    gRPA_h[is]->Draw("P");
  }
  
  dashedLine(xmin_,1.,xmax_,1.,1,1);
  TLegend *leg_h= new TLegend(0.22, 0.68, 0.465, 0.876);
  SetLegendStyle(leg_h);
//  leg_h->SetTextSize(0.036);
  leg_h->SetTextSize(0.042);
  leg_h->SetTextFont(22); //22);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
//      leg_h -> AddEntry(gRPA_h[is],Form(" #varUpsilon(%dS)",is+1),"lp");
	leg_h -> AddEntry(gRPA_h[is],Form(" #Upsilon(%dS)",is+1),"lp");
	leg_h -> Draw("same");
    }
  }


  //// draw text
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
//  globtex->DrawLatex(0.7614, sz_init-sz_step, "|y_{CM}^{#varUpsilon}| < 1.93");
//  globtex->DrawLatex(0.7614, sz_init-sz_step, "|y_{CM}^{#Upsilon}| < 1.93");
  globtex->DrawLatex(0.92, sz_init-sz_step, "|y_{CM}^{#Upsilon}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  

  double sys_global_x_ = 1.6;
  TBox *globalUncBox_ = new TBox(xmax_-sys_global_x_,1-sys_global_y_Lo,xmax_,1+sys_global_y_Hi);
  globalUncBox_ -> SetLineColor(kBlack);
  globalUncBox_ -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox_ -> SetLineWidth(1);
  globalUncBox_ -> Draw("l same");
  
  CMS_lumi_raaCent( c2, iPeriod, iPos );

	c2->Update();
  c2->SaveAs("plots/RpA_vs_pT_1D.pdf");
  c2->SaveAs("plots/RpA_vs_pt_1D.png");
  
  for (int is=0; is<nState; is++){
    double val[npoint_l[is]]; double val_stat[npoint_l[is]]; double val_sys[npoint_l[is]];
    for (int ipt=0; ipt< npoint_l[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_l[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRPA_l[is]->GetErrorX(ipt);
      eytmp=gRPA_l[is]->GetErrorY(ipt);
      relsys=hSys_l[is]->GetBinContent(ipt+1);
      val[ipt] = pytmp; val_stat[ipt] = eytmp; val_sys[ipt] = pytmp*relsys;
    }
    if(is==0){
      cout << "$-1.93  < \\rap <  -1.20$ & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$-1.20  < \\rap <  -0.80$ & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
      cout << "$-0.80  < \\rap <  -0.40$ & " << Form("%.3f",val[2])  << " & " << Form("%.3f",val_stat[2]) << " & " << Form("%.3f",val_sys[2]) << " \\\\ " << endl;
      cout << "$-0.40  < \\rap <   0.00$ & " << Form("%.3f",val[3])  << " & " << Form("%.3f",val_stat[3]) << " & " << Form("%.3f",val_sys[3]) << " \\\\ " << endl;
      cout << "$0.00  < \\rap <  0.40$ & "  << Form("%.3f",val[4])  << " & " << Form("%.3f",val_stat[4]) << " & " << Form("%.3f",val_sys[4]) << " \\\\ " << endl;
      cout << "$0.40  < \\rap <  0.80$ & " << Form("%.3f",val[5])  << " & " << Form("%.3f",val_stat[5]) << " & " << Form("%.3f",val_sys[5]) << " \\\\ " << endl;
      cout << "$0.80  < \\rap <  1.20$ & " << Form("%.3f",val[6])  << " & " << Form("%.3f",val_stat[6]) << " & " << Form("%.3f",val_sys[6]) << " \\\\ " << endl;
      cout << "$1.20  < \\rap <  1.93$ & " << Form("%.2f",val[7])  << " & " << Form("%.2f",val_stat[7]) << " & " << Form("%.2f",val_sys[7]) << " \\\\ " << endl;
    }
    else if(is==1){
      cout << "$-1.93  < \\rap <  -0.80$ & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$-0.80  < \\rap <   0.00$ & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
      cout << "$0.00  < \\rap <  0.80$ & " << Form("%.3f",val[2])  << " & " << Form("%.3f",val_stat[2]) << " & " << Form("%.3f",val_sys[2]) << " \\\\ " << endl;
      cout << "$0.80  < \\rap <  1.93$ & " << Form("%.2f",val[3])  << " & " << Form("%.2f",val_stat[3]) << " & " << Form("%.2f",val_sys[3]) << " \\\\ " << endl;
    }
    else if(is==2){
      cout << "$-1.93  < \\rap <   0.00$ & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$0.00  < \\rap <   1.93$ & " << Form("%.2f",val[1])  << " & " << Form("%.2f",val_stat[1]) << " & " << Form("%.2f",val_sys[1]) << " \\\\ " << endl;
    }
  }

  for (int is=0; is<nState; is++){
    double val[npoint_h[is]]; double val_stat[npoint_h[is]]; double val_sys[npoint_h[is]];
    for (int ipt=0; ipt< npoint_h[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_h[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRPA_h[is]->GetErrorX(ipt);
      eytmp=gRPA_h[is]->GetErrorY(ipt);
      relsys=hSys_h[is]->GetBinContent(ipt+1);
      val[ipt] = pytmp; val_stat[ipt] = eytmp; val_sys[ipt] = pytmp*relsys;
    }
    if(is==0){
      cout << "$\\pt < 2$ \\GeVc & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$2 < \\pt < 4$ \\GeVc & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
      cout << "$4 < \\pt < 6$ \\GeVc & " << Form("%.3f",val[2])  << " & " << Form("%.3f",val_stat[2]) << " & " << Form("%.3f",val_sys[2]) << " \\\\ " << endl;
      cout << "$6 < \\pt < 9$ \\GeVc & " << Form("%.3f",val[3])  << " & " << Form("%.3f",val_stat[3]) << " & " << Form("%.3f",val_sys[3]) << " \\\\ " << endl;
      cout << "$9 < \\pt < 12$ \\GeVc & " << Form("%.3f",val[4])  << " & " << Form("%.3f",val_stat[4]) << " & " << Form("%.3f",val_sys[4]) << " \\\\ " << endl;
      cout << "$12 < \\pt < 30$ \\GeVc & " << Form("%.3f",val[5])  << " & " << Form("%.3f",val_stat[5]) << " & " << Form("%.3f",val_sys[5]) << " \\\\ " << endl;
    }
    else if(is==1){
      cout << "$\\pt < 4$ \\GeVc & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$4 < \\pt < 9$ \\GeVc & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
      cout << "$9 < \\pt < 30$ \\GeVc & " << Form("%.3f",val[2])  << " & " << Form("%.3f",val_stat[2]) << " & " << Form("%.3f",val_sys[2]) << " \\\\ " << endl;
    }
    else if(is==2){
      cout << "$\\pt < 6$ \\GeVc & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$6 < \\pt < 30$ \\GeVc & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
    }
  }

  return;

} // end of main func.

