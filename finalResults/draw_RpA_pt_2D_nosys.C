#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
void draw_RpA_pt_2D_nosys(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmin = 0;
  double xmax = 30;
//  double relsys = 0.1;

  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3., 12.};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn_l[nState];
  TFile* fIn_h[nState];
	TGraphErrors* gRPA_l[nState];
	TGraphErrors* gRPA_h[nState];
	TGraphErrors* gRPA_sys_l[nState];
	TGraphErrors* gRPA_sys_h[nState];
  for (int is=0; is<nState; is++){
  	fIn_l[is] = new TFile(Form("Ups_%d_RPA_lowPt.root",is+1),"READ");
  	fIn_h[is] = new TFile(Form("Ups_%d_RPA_highPt.root",is+1),"READ");
    gRPA_l[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_pt");
    gRPA_h[is]=(TGraphErrors*)fIn_h[is]->Get("gRPA_pt");
    gRPA_sys_l[is]=(TGraphErrors*)fIn_l[is]->Get("gRPA_pt");
    gRPA_sys_h[is]=(TGraphErrors*)fIn_h[is]->Get("gRPA_pt");
  }

  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys[nState];
  int npoint[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematics/mergedSys_ups%ds.root",is+1),"READ");
    hSys[is]=(TH1D*)fInSys[is]->Get("hptRPA_merged");
    npoint[is] = hSys[is]->GetSize()-2;
    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRPA_l[is]->GetN()) {cout << "Error!! data file and syst. file have dAifferent binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_l[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_l[is]->GetErrorX(ipt);
      eytmp=gRPA_l[is]->GetErrorY(ipt);
      relsys=0.0;
      cout << ipt <<"th bin low pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      gRPA_l[is]->SetPointError(ipt, 0, eytmp);
      if (is==0) gRPA_sys_l[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys_l[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRPA_sys_l[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
 
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA_h[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA_h[is]->GetErrorX(ipt);
      eytmp=gRPA_h[is]->GetErrorY(ipt);
      relsys=0.0;
      cout << endl;
      cout << ipt <<"th bin high pt RPA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
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
    SetGraphStyle2(gRPA_h[is], 3+is, 3+is); 
    SetGraphStyleSys2(gRPA_sys_h[is], 3+is); 
	}
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.033);
  
  //// legend
  //// axis et. al
  gRPA_sys_l[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gRPA_sys_l[0]->GetXaxis()->CenterTitle();
  gRPA_sys_l[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRPA_sys_l[0]->GetYaxis()->CenterTitle();
  gRPA_sys_l[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRPA_sys_l[0]->SetMinimum(0.0);
  gRPA_sys_l[0]->SetMaximum(1.7);
  gRPA_sys_l[0]->GetXaxis()->SetNdivisions(505);
  gRPA_sys_h[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
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
    if ( is==0) {gRPA_sys_l[is]->Draw("A5");gRPA_sys_h[is]->Draw("5");}
/*    else if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++){
        box68per[ipt]->Draw("l");
      }
      gRPA_sys[is]->Draw("5");
    }
// */
    else {gRPA_sys_l[is]->Draw("5");gRPA_sys_h[is]->Draw("5");}
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
    gRPA_l[is]->Draw("P"); gRPA_h[is]->Draw("P");
  }
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.20, 0.62, 0.405, 0.896);
  SetLegendStyle(leg);
  leg->SetTextSize(0.034);
  TLegend *leg_up= new TLegend(0.57, 0.50, 0.78, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(16.,0.532,16.,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
      leg -> AddEntry(gRPA_l[is],Form(" #Upsilon(%dS), 0<y_{CM}^{#varUpsilon}<1.93",is+1),"lp");
      leg -> AddEntry(gRPA_h[is],Form(" #Upsilon(%dS), -1.93<y_{CM}^{#varUpsilon}<0",is+1),"lp");
      leg -> Draw("same");
    }
  }


  //// draw text
  double sz_init = 0.925; double sz_step = 0.1975;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
  globtex->DrawLatex(0.73, sz_init-sz_step, "|y_{CM}^{#varUpsilon}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  

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
  globalUncBox -> Draw("l same");

  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("RpA_vs_pt_2D_nosys.pdf");
  c1->SaveAs("RpA_vs_pt_2D_nosys.png");

	return;

} // end of main func.

