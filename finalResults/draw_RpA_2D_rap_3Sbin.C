#include "JaebeomStyle.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
void draw_RpA_2D_rap_3Sbin(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = false;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmin = -1.93;
  double xmax = 1.93;
//  double relsys = 0.1;
/*
  double exsys_1s[2] =  {1.93/2, 1.93/2};
  double exsys_2s[2] =  {1.93/2, 1.93/2};
  double exsys_3s[2] =  {1.93/2, 1.93/2};
*/
  double xrange_width = 0.012*(xmax-xmin); //0.038; //0.035;
  double exsys_1s[2] =  {xrange_width,xrange_width};
  double exsys_2s[2] =  {xrange_width,xrange_width};
  double exsys_3s[2] =  {xrange_width,xrange_width};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn_l[nState];
  TFile* fIn_h[nState];
	TGraphErrors* gRPA_l[nState];
	TGraphErrors* gRPA_h[nState];
	TGraphErrors* gRPA_sys_l[nState];
	TGraphErrors* gRPA_sys_h[nState];
  for (int is=0; is<nState; is++){
  	fIn_l[is] = new TFile(Form("Ups_%d_RPA_2D_rap_3Sbin.root",is+1),"READ");
    gRPA_l[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_rap_low");
    gRPA_h[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_rap_high");
    gRPA_sys_l[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_rap_low");
    gRPA_sys_h[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_rap_high");
  }

  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys_h[nState];
  TH1D* hSys_l[nState];
  int npoint[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematics/rap2D_3Sbin_mergedSys_ups%ds.root",is+1),"READ");
    hSys_l[is]=(TH1D*)fInSys[is]->Get("hrapRPA_merged1");
    hSys_h[is]=(TH1D*)fInSys[is]->Get("hrapRPA_merged2");
    npoint[is] = hSys_l[is]->GetSize()-2;
    if(npoint[is] != hSys_h[is]->GetSize()-2) {cout << "Inconsistent number of bins! " << endl; return; }
    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  double pxp[2] = {-1.93/2, 1.93/2};
  double pxpshift[3] = {-0.1, 0, 0.1};

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRPA_l[is]->GetN()) {cout << "Error!! data file and syst. file have dAifferent binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_l[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_l[is]->GetErrorX(ipt);
      eytmp=gRPA_l[is]->GetErrorY(ipt);
      relsys=hSys_l[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin low pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      gRPA_l[is]->SetPoint(ipt,pxp[ipt]+pxpshift[is],pytmp);
      gRPA_sys_l[is]->SetPoint(ipt,pxp[ipt]+pxpshift[is],pytmp);
      gRPA_l[is]->SetPointError(ipt, 0, eytmp);
      if (is==0) gRPA_sys_l[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys_l[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRPA_sys_l[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
 
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_h[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_h[is]->GetErrorX(ipt);
      eytmp=gRPA_h[is]->GetErrorY(ipt);
      relsys=hSys_h[is]->GetBinContent(ipt+1);
      cout << endl;
      cout << ipt <<"th bin high pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      gRPA_h[is]->SetPoint(ipt,pxp[ipt]+pxpshift[is],pytmp);
      gRPA_sys_h[is]->SetPoint(ipt,pxp[ipt]+pxpshift[is],pytmp);
      gRPA_h[is]->SetPointError(ipt, 0, eytmp);
      if (is==0) gRPA_sys_h[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys_h[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRPA_sys_h[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
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
  globtex->SetTextFont(22);
//  globtex->SetTextFont(42);
//  globtex->SetTextSize(0.035);
//  globtex->SetTextSize(0.038);
  globtex->SetTextSize(0.042);

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
  gRPA_sys_h[0]->GetXaxis()->SetTitle("y_{CM}^{#Upsilon}");
  gRPA_sys_h[0]->GetXaxis()->CenterTitle();
  gRPA_sys_h[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRPA_sys_h[0]->GetYaxis()->CenterTitle();
  gRPA_sys_h[0]->GetXaxis()->SetLimits(xmin,xmax);
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
  

/*  for (int is=1; is<2; is++){
    gRPA_sys_l[is]->GetYaxis()->SetRangeUser(0,1.7);
    gRPA_sys_l[is]->GetXaxis()->SetRangeUser(-1.93,1.93);
    gRPA_sys_l[is]->Draw("A5");
    gRPA_l[is]->Draw("P");
  }
  */
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
      leg -> AddEntry(gRPA_l[is],Form(" #Upsilon(%dS)",is+1),"lp");
      leg -> Draw("same");
    }
  }


  //// draw text
//  double sz_init = 0.925; double sz_step = 0.1975;
  double sz_init = 1.1; double sz_step = 0.3;

//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
//  globtex->DrawLatex(0.745, sz_init-sz_step, "p_{T}^{#Upsilon} < 6 GeV/c");
  globtex->DrawLatex(0.92, sz_init-sz_step, "p_{T}^{#Upsilon} < 6 GeV/c");
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
  c1->SaveAs("plots/RpA_vs_rap_2D_3Sbin_lowPt.pdf");
  c1->SaveAs("plots/RpA_vs_rap_2D_3Sbin_lowPt.png");


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
  
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg_h= new TLegend(0.22, 0.68, 0.465, 0.876);
  SetLegendStyle(leg_h);
//  leg_h->SetTextSize(0.036);
  leg_h->SetTextSize(0.042);
  leg_h->SetTextFont(22);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
      leg_h -> AddEntry(gRPA_h[is],Form(" #Upsilon(%dS)",is+1),"lp");
      leg_h -> Draw("same");
    }
  }


  //// draw text
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
//  globtex->DrawLatex(0.68, sz_init-sz_step, "6 < p_{T}^{#Upsilon} < 30 GeV/c");
  globtex->DrawLatex(0.92, sz_init-sz_step, "6 < p_{T}^{#Upsilon} < 30 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  

  globalUncBox -> Draw("l same");
  
  CMS_lumi_raaCent( c2, iPeriod, iPos );

	c2->Update();
  c2->SaveAs("plots/RpA_vs_rap_2D_3Sbin_highPt.pdf");
  c2->SaveAs("plots/RpA_vs_rap_2D_3Sbin_highPt.png");
  
  for (int is=0; is<nState; is++){
    double val[npoint[is]]; double val_stat[npoint[is]]; double val_sys[npoint[is]];
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_l[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRPA_l[is]->GetErrorX(ipt);
      eytmp=gRPA_l[is]->GetErrorY(ipt);
      relsys=hSys_l[is]->GetBinContent(ipt+1);
      val[ipt] = pytmp; val_stat[ipt] = eytmp; val_sys[ipt] = pytmp*relsys;
    }
      cout << "$-1.93  < y_{CM} <   0.00$ & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$ 0.00  < y_{CM} <   1.93$ & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
  }

  cout << endl;
  cout << endl;
  for (int is=0; is<nState; is++){
    double val[npoint[is]]; double val_stat[npoint[is]]; double val_sys[npoint[is]];
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_h[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRPA_h[is]->GetErrorX(ipt);
      eytmp=gRPA_h[is]->GetErrorY(ipt);
      relsys=hSys_h[is]->GetBinContent(ipt+1);
      val[ipt] = pytmp; val_stat[ipt] = eytmp; val_sys[ipt] = pytmp*relsys;
    }
      cout << "$-1.93  < y_{CM} <   0.00$ & " << Form("%.3f",val[0])  << " & " << Form("%.3f",val_stat[0]) << " & " << Form("%.3f",val_sys[0]) << " \\\\ " << endl;
      cout << "$ 0.00  < y_{CM} <   1.93$ & " << Form("%.3f",val[1])  << " & " << Form("%.3f",val_stat[1]) << " & " << Form("%.3f",val_sys[1]) << " \\\\ " << endl;
  }


	return;

} // end of main func.

