#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"

void draw_RpA_rap(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 502; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmin = -2;
  double xmax = 2;
//  double relsys = 0.1;

  double exsys_1s[8] =  {0.73/2, 0.2, 0.2, 0.2, 0.2,0.2,0.2, 0.73/2};
  double exsys_2s[4] =  {1.13/2, 0.4,0.4,1.13/2};
  double exsys_3s[2] =  {1.93/2, 1.93/2};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRPA[nState];
	TGraphErrors* gRPA_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_RPA.root",is+1),"READ");
    gRPA[is]=(TGraphErrors*)fIn[is]->Get("gRPA_rap");
    gRPA_sys[is]=(TGraphErrors*)fIn[is]->Get("gRPA_rap");
    //cout << "gRAA["<<is<<"] = " <<gRAA[is] << endl;
  }
  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys[nState];
  int npoint[nState]={8,4,2};/*
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematic/mergedSys_ups%ds.root",is+1),"READ");
    hSys[is]=(TH1D*)fInSys[is]->Get("hptRAA_merged");
    npoint[is] = hSys[is]->GetSize()-2;
    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  */
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRPA[is]->GetN()) {cout << "Error!! data file and syst. file have dAifferent binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRPA[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRPA[is]->GetErrorX(ipt);
      eytmp=gRPA[is]->GetErrorY(ipt);
      relsys=0.00;
      //relsys=0.05;
      //relsys=hSys[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      //// 1) remove ex from gRAA
      gRPA[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys
      //gRAA_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      if (is==0) gRPA_sys[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRPA_sys[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRPA_sys[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
    }
  }
 
  ////////////////////////////////////////////////////////////////
  int ulstate = 2; //3S
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
  gRPA_sys[0]->GetXaxis()->SetTitle("y_{CM}");
  gRPA_sys[0]->GetXaxis()->CenterTitle();
  gRPA_sys[0]->GetYaxis()->SetTitle("R_{pPb}");
  gRPA_sys[0]->GetYaxis()->CenterTitle();
  gRPA_sys[0]->GetXaxis()->SetLimits(xmin,xmax);
  gRPA_sys[0]->SetMinimum(0.0);
  gRPA_sys[0]->SetMaximum(1.5);


  gRPA_sys[0]->GetXaxis()->SetNdivisions(505);
  if (isArrow == true){
        gRPA_sys[2]->SetPoint(0,-10,-10);
        gRPA_sys[2]->SetPointError(0,0,0);
        gRPA_sys[2]->SetPoint(1,-11,-11);
        gRPA_sys[2]->SetPointError(1,0,0);
        gRPA_sys[2]->SetPoint(2,-12,-12);
        gRPA_sys[2]->SetPointError(2,0,0);
        gRPA[2]->SetPoint(0,-10,-10);
        gRPA[2]->SetPointError(0,0,0);
        gRPA[2]->SetPoint(1,-11,-11);
        gRPA[2]->SetPointError(1,0,0);
        gRPA[2]->SetPoint(2,-12,-12);
        gRPA[2]->SetPointError(2,0,0);
        gRPA_sys[2]->GetHistogram()->GetXaxis()->SetLimits(0,30);
        gRPA_sys[2]->GetHistogram()->GetXaxis()->SetRangeUser(0,30);
        gRPA_sys[2]->SetMinimum(0.0);
        gRPA_sys[2]->SetMaximum(1.3);
        gRPA[2]->GetHistogram()->GetXaxis()->SetRangeUser(0,30);
        gRPA[2]->GetHistogram()->GetXaxis()->SetLimits(0,30);
        gRPA[2]->SetMinimum(0.0);
        gRPA[2]->SetMaximum(1.3);
      }
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nState; is++){
    if ( is==0) {gRPA_sys[is]->Draw("A5");}
    else if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++){
        box68per[ipt]->Draw("l");
      }
      gRPA_sys[is]->Draw("5");
    }
    else {gRPA_sys[is]->Draw("5");}
  }
  
  for(int is=0;is<nState;is++){
    if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++) {
        arr95per[ipt]->Draw();
      }
      gRPA[is]->Draw("P");
    }
    else {gRPA[is]->Draw("P");}
  }
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.50, 0.72, 0.705, 0.896);
  SetLegendStyle(leg);
  TLegend *leg_up= new TLegend(0.57, 0.50, 0.78, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(16.,0.532,16.,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
      leg -> AddEntry(gRPA[is],Form(" #Upsilon(%dS)",is+1),"lp");
      leg -> Draw("same");
    }
  }
  else {
    leg -> AddEntry(gRPA[0]," #Upsilon(1S)","lp");
    leg -> AddEntry(gRPA[1]," #Upsilon(2S)","lp");
    leg -> AddEntry(gRPA[2]," #Upsilon(3S)","lp");
    TLegendEntry *ent=leg_up->AddEntry("ent"," #Upsilon(3S) 68\% CL","f");
    ent->SetLineColor(kGreen+3);
    ent->SetFillColorAlpha(kGreen-6,0.5);
    ent->SetFillStyle(1001);
    ent=leg_up->AddEntry("ent"," #Upsilon(3S) 95\% CL","f");
    ent->SetLineColor(kWhite);
//    leg_up->SetTextSize(0.03);
    leg->Draw("same");
    leg_up->Draw("same");
    arrLeg->Draw();
  }


  //// draw text
  double sz_init = 0.925; double sz_step = 0.0675;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T}^{#mu#mu} < 30 GeV/c");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "|y^{#mu#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 1.93");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");

  

  //Global Unc.
 
  double TAA_unc_Global_Hi = 0.068;
  double TAA_unc_Global_Lo = 0.072;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_pa*lumi_unc_pa);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.4;
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("RpA_vs_rap_isArrow%d_asym.pdf",(int)isArrow));
  c1->SaveAs(Form("RpA_vs_rap_isArrow%d_asym.png",(int)isArrow));

	return;

} // end of main func.

