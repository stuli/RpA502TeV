#include "JaebeomStyle.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"
void theory_comp_Ferreiro_RpA_1D_rap(int drawState=2)
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
    npoint_l[is] = hSys_l[is]->GetSize()-2;
    cout << "*** Y("<<is+1<<") rap : # of point = " << npoint_l[is] << endl;
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
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  
  //// legend
  //// axis et. al
  if(drawState==0){
    gRPA_sys_l[0]->GetXaxis()->SetTitle("y_{CM}^{#varUpsilon}");
    gRPA_sys_l[0]->GetXaxis()->CenterTitle();
    gRPA_sys_l[0]->GetYaxis()->SetTitle("R_{pPb}");
    gRPA_sys_l[0]->GetYaxis()->CenterTitle();
    gRPA_sys_l[0]->GetXaxis()->SetLimits(xmin,xmax);
    gRPA_sys_l[0]->SetMinimum(0.0);
    gRPA_sys_l[0]->SetMaximum(1.7);
    gRPA_sys_l[0]->GetXaxis()->SetNdivisions(505);
  }
  else{
    gRPA_sys_l[drawState-1]->GetXaxis()->SetTitle("y_{CM}^{#varUpsilon}");
    gRPA_sys_l[drawState-1]->GetXaxis()->CenterTitle();
    gRPA_sys_l[drawState-1]->GetYaxis()->SetTitle("R_{pPb}");
    gRPA_sys_l[drawState-1]->GetYaxis()->CenterTitle();
    gRPA_sys_l[drawState-1]->GetXaxis()->SetLimits(xmin,xmax);
    gRPA_sys_l[drawState-1]->SetMinimum(0.0);
    gRPA_sys_l[drawState-1]->SetMaximum(1.7);
    gRPA_sys_l[drawState-1]->GetXaxis()->SetNdivisions(505);
  }
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  if(drawState==0){
    for (int is=0; is<nState; is++){
      if ( is==0) {gRPA_sys_l[is]->Draw("A5");}
      else {gRPA_sys_l[is]->Draw("5");}
    }

    for(int is=0;is<nState;is++){
      gRPA_l[is]->Draw("P");
    }
  }
  else{
      gRPA_sys_l[drawState-1]->Draw("A5");
      gRPA_l[drawState-1]->Draw("P");
  }
  
  dashedLine(xmin,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.22, 0.68, 0.465, 0.876);
  SetLegendStyle(leg);
  leg->SetTextSize(0.036);
  leg->SetTextFont(22);
  TLegend *leg_up= new TLegend(0.57, 0.50, 0.78, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(16.,0.532,16.,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);


  //// draw text
  double sz_init = 0.925; double sz_step = 0.1975;
  globtex->DrawLatex(0.717, sz_init-sz_step, "p_{T}^{#varUpsilon} < 30 GeV/c");

  

  //Global Unc.
 
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


  TFile *f_e = new TFile("Theory/Elena_RpPb_graph.root");
  TFile *f_ep = new TFile("Theory/Elena_RpPb_graph_eps09.root");

  Int_t sh_color[] = { kOrange+3, kBlue+3, kGreen+4 }; 
  Int_t sh_colorep[] = { kPink+7, kBlue+3, kGreen+4 }; //{ kMagenta-7, kBlue+3, kGreen+4 }; 
  TGraph *gsh[nState];
  for(int i=0;i<nState;i++){
    gsh[i] = (TGraph*) f_e -> Get(Form("RpA_%ds_rap_shade",i+1));
    gsh[i] -> SetLineStyle(7);
    gsh[i] -> SetLineWidth(2);
    gsh[i] -> SetLineColor(sh_color[0]);
    //gsh[i] -> SetFillColor(sh_color[0]);
    gsh[i] -> SetFillColor(0);
    gsh[i] -> SetFillStyle(0);
  }

  if(drawState==0){
    for (int is=0; is<nState; is++){
      gsh[is]->Draw("L");
      gsh[is]->Draw("f same");
      leg -> AddEntry(gRPA_l[is],Form(" #varUpsilon(%dS)",is+1),"lp");
    }
  }
  else{
    //gsh[drawState-1] -> Draw("L");
    //gsh[drawState-1] -> Draw("f same");
    gsh[drawState-1] -> Draw("same");  
    leg -> AddEntry(gRPA_l[drawState-1],"Data","lp");
  }

  TGraph *gsh_eps[nState];
  for(int i=0;i<nState;i++){
    gsh_eps[i] = (TGraph*) f_ep -> Get(Form("RpA_%ds_rap_shade",i+1));
    gsh_eps[i] -> SetLineStyle(3);
    gsh_eps[i] -> SetLineWidth(2);
    gsh_eps[i] -> SetLineColor(sh_colorep[0]);
    //gsh_eps[i] -> SetFillColor(sh_colorep[0]);
    gsh_eps[i] -> SetFillColor(0);
    gsh_eps[i] -> SetFillStyle(0);
  }

  if(drawState==0){
    for (int is=0; is<nState; is++){
      gsh_eps[is]->Draw("L");
      gsh_eps[is]->Draw("f same");
    }
  }
  else{
    //gsh_eps[drawState-1] -> Draw("L");  
    //gsh_eps[drawState-1] -> Draw("f same");
    gsh_eps[drawState-1] -> Draw("same");
  }

  leg -> Draw("same");

  //drawText("E. Ferreiro and J. Lansberg",0.387,0.88,1,19);
  //drawText("nPDF+Comover",0.387,0.84,1,19);
  drawText(Form(" #varUpsilon(%dS)",drawState),0.219,0.82,1,20);
//  drawText("nCTEQ15+Comover",0.387,0.84,1,19);
//  drawText("EPS09+Comover",0.387,0.80,1,19);
  double leg1_x1 = 0.36;
  double leg1_x2 = 0.67;
  double leg1_y1 = 0.70;
  double leg1_y2 = 0.886;

  if(drawState!=0) leg1_y1 = 0.75;
  TLegend *leg1= new TLegend(leg1_x1,leg1_y1,leg1_x2,leg1_y2);
  SetLegendStyle(leg1);
  leg1->SetTextSize(0.034);
//  leg1->SetHeader("Ferreiro, nCTEQ15+comovers","");
  for(int i=0;i<nState;i++){
    if(drawState==0){leg1->AddEntry(gsh[i],Form("#varUpsilon(%dS)",i+1),"f");}
    else{
      if(i==0){
        leg1->AddEntry(gsh[drawState-1],"CIM + nCTEQ15","f");  
        leg1->AddEntry(gsh_eps[drawState-1],"CIM + EPS09 LO","f"); 
      }
    }
  }
  leg1 -> Draw("same"); 

  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("plots/Theory_Ferreiro_vs_rap_1D_drawState%d.pdf",drawState));
  c1->SaveAs(Form("plots/Theory_Ferreiro_vs_rap_1D_drawState%d.png",drawState));


  return;

} // end of main func.

