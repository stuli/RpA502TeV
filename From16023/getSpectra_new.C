#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "cutsAndBin.h"
#include "multiTreeUtil.h"
using namespace std;

TString ResultDir  = "nominalFits";


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 
valErr getYield(int state=0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0, 	float dphiEp2Low=0,  float dphiEp2High=0) ;

double getScale(int fTAA = 1, double* TAA =  TAA1s, double* centBin = centBin1s, int nCentBins=0);

void stripErrorBars( TH1* h =0, double defaultErr = 0 ); 

void getSpectra_new(int state = 2 ) {  

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  int nPtBins=0;
  int nYBins=0;
  double* ptBin;
  double* yBin;
  int nCentBins=0;
  double* centBin;
  double* nPart;  // In order from peripheral to central 
  double* nColl;  // In order from central to peripheral 
  double* TAA;
  //  double nColl1s[nCentBins] = {1819,1432,1005,606,349,186,90.7,40.1,7.67};
  //  double nPart1s[nCentBins] = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};



  if ( state == 1 ) { 
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;    yBin = yBin1S; 
    nCentBins = nCentBins1s;  centBin = centBin1s; nPart = nPart1s; nColl = nColl1s; TAA = TAA1s;
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nYBins = nYBins2S;    yBin = yBin2S; 
    nCentBins = nCentBins2s;  centBin = centBin2s; nPart = nPart2s; nColl = nColl2s; TAA = TAA2s;
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nYBins = nYBins3S;    yBin = yBin3S; 
    nCentBins = nCentBins3s;  centBin = centBin3s; nPart = nPart3s; nColl = nColl3s; TAA = TAA3s;
  }
  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];  
  double yMin = yBin[0];    double yMax = yBin[nYBins];  
  double centMin = centBin[0];    double centMax = centBin[nCentBins];  

  TH1D* hrapEffAA;
  TH1D* hrapEffPP;
  TH1D* hrapAccAA;
  TH1D* hrapAccPP;
  TH1D* hrapSigAA;
  TH1D* hrapSigPP;

  TH1D* hptEffAA;
  TH1D* hptEffPP;
  TH1D* hptAccAA;
  TH1D* hptAccPP;
  TH1D* hptSigAA;
  TH1D* hptSigPP;

  TH1D* hcentEffAA;
  TH1D* hcentEffPP;
  TH1D* hcentSigAA;
  TH1D* hcentSigPP;  // There is only one bin for pp. 

  TH1D* hIntAccAA;
  TH1D* hIntAccPP;

  // ##################################
  // ~*~*~*~*~* Acceptance ~*~*~*~*~*~*
  // ##################################

  TFile* infacc = new TFile(Form("acceptance/acceptance_wgt_norm_%dS.root",state));
  //TFile* infacc = new TFile(Form("compareDataMc/20170405_AccCleaned/20170106/acceptance_wgt_%dS_20170529.root",state));
  hrapAccAA  = (TH1D*)infacc->Get(Form("hrapAccAA%dS",state));
  hrapAccPP  = (TH1D*)infacc->Get(Form("hrapAccPP%dS",state));
  hptAccAA  = (TH1D*) infacc->Get(Form("hptAccAA%dS",state));
  hptAccPP  = (TH1D*) infacc->Get(Form("hptAccPP%dS",state));
  hIntAccAA = (TH1D*) infacc->Get(Form("hIntAccAA%dS",state));
  hIntAccPP = (TH1D*) infacc->Get(Form("hIntAccPP%dS",state));

  // temporary fix for acceptance uncertainty bar 
  stripErrorBars(hrapAccPP); 
  stripErrorBars(hrapAccAA); 
  stripErrorBars(hptAccPP); 
  stripErrorBars(hptAccAA); 
  stripErrorBars(hIntAccAA); 
  stripErrorBars(hIntAccPP); 



  // ##################################
  // ~*~*~*~*~* Rapidity ~*~*~*~*~*~*~*
  // ##################################
  //TFile* inf = new TFile(Form("efficiency/efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnpWeight1_tnpIdx0.root",state));
  TFile* inf = new TFile(Form("efficiency/efficiencyTable/efficiency_ups%ds_useDataPtWeight1_tnp_trgId0_trkId0_muId-100_staId-100.root",state));
  hrapEffAA  = (TH1D*)inf->Get("hrapEffAA");
  hrapEffPP  = (TH1D*)inf->Get("hrapEffPP");
  stripErrorBars(hrapEffAA);
  stripErrorBars(hrapEffPP);

  TCanvas* c_rap_eff =  new TCanvas("c_rap_eff","",400,400);
  c_rap_eff->cd();
  cleverRange(hrapEffAA, 1.3, 0);
  handsomeTH1(hrapEffAA,2);
  handsomeTH1(hrapEffPP,2);
  hrapEffAA->Draw();
  hrapEffPP->SetMarkerStyle(24);
  hrapEffPP->Draw("same");
  
  // ############################
  // ~*~*~*~*~* Pt ~*~*~*~*~*~*~*
  // ############################
  hptEffAA  = (TH1D*)inf->Get("hptEffAA");
  hptEffPP  = (TH1D*)inf->Get("hptEffPP");
  stripErrorBars(hptEffAA);
  stripErrorBars(hptEffPP);

  TCanvas* c1 =  new TCanvas("c1","",400,400);
  c1->cd();
  cleverRange(hptEffAA, 1.3, 0);
  handsomeTH1(hptEffAA,2);
  handsomeTH1(hptEffPP,2);
  hptEffAA->Draw();
  hptEffPP->SetMarkerStyle(24);
  hptEffPP->Draw("same");
  
  
  // ############################
  // ~*~*~*~*~* Cent ~*~*~*~*~*~*
  // ############################
  hcentEffAA = (TH1D*) inf -> Get("hcentEffAA");
  hcentEffPP = (TH1D*) inf -> Get("hcentintEffPP");
  stripErrorBars(hcentEffAA);
  stripErrorBars(hcentEffPP);


  hcentSigAA = (TH1D*) hcentEffAA -> Clone("hcentSigAA_cent");
  hcentSigPP = (TH1D*) hcentEffPP -> Clone("hcentintSigPP");
  hcentSigAA -> Reset(); 
  hcentSigPP -> Reset(); 
  
  TH1D *hSetBin;
  hSetBin = (TH1D*) hcentSigAA ->Clone("hSetBin");
  TH1D *hcentEffAA_int = (TH1D*) inf -> Get("hcentintEffAA");
  TH1D *hcentSigAA_int = (TH1D*) hcentEffAA_int -> Clone("hcentSigAA_int");
  hcentSigAA_int->Reset();



  // signals :
  TH1D* hRAAraw_rap;   // w/o efficiency correction
  TH1D* hRAA_rap;   // w/ efficiency correction
  TH1D* hRAAraw_pt;   // w/o efficiency correction
  TH1D* hRAA_pt;   // w efficiency correction
  TH1D* hRAAraw_cent; // w/o efficiency correction
  TH1D* hRAA_cent;   // w efficiency correction

  //***Rapidity***
  hrapSigAA = (TH1D*)  hrapEffAA->Clone("hrapSigAA");
  hrapSigPP = (TH1D*)  hrapEffPP->Clone("hrapSigPP");
  hrapSigAA->Reset();
  hrapSigPP->Reset();
    for ( int iy = 1 ; iy<= nYBins ; iy++) {
      valErr yieldPP = getYield(state, kPPDATA, ptMin, ptMax, yBin[iy-1], yBin[iy], centMin, centMax, 0, 100);
      valErr yieldAA = getYield(state, kAADATA, ptMin, ptMax, yBin[iy-1], yBin[iy], centMin, centMax, 0, 100);
      hrapSigAA->SetBinContent( iy, yieldAA.val ) ;
      hrapSigAA->SetBinError( iy, yieldAA.err ) ;
      hrapSigPP->SetBinContent( iy, yieldPP.val ) ;
      hrapSigPP->SetBinError( iy, yieldPP.err ) ;
    }
    
  
  //*****RAA Rap ******
  TCanvas* cRAA_rap =  new TCanvas("cRAA_rap","",400,400);
  hRAAraw_rap = (TH1D*)hrapSigAA->Clone("hRAAraw_rap");
  hRAAraw_rap->Divide( hrapSigPP );
  double scale_rap = getScale(nCentBins+1, TAA, centBin, nCentBins);

//  hRAAraw_rap->Scale( 26000000. / 351 ) ;   // pp : 26pb-1,  PbPb : 351 microBarn-1
//  hRAAraw_rap->Scale( 1./ (208.*208) );
  
  hRAAraw_rap->Scale(scale_rap);
  hRAAraw_rap->SetAxisRange(0,1.2,"Y");
  hRAAraw_rap->SetYTitle("R_{AA} (efficiency UNcorrected)");
  hRAAraw_rap->Draw();
  jumSun(0,1,30,1);
  
  TCanvas* cRAA_rap_effcor =  new TCanvas("cRAA_rap_effcor","",400,400);
  hRAA_rap = (TH1D*)hRAAraw_rap->Clone("hRAA_rap");

  TH1D* relativeEff_rap = (TH1D*)hrapEffAA->Clone("relEffAA_rap");
  relativeEff_rap->Divide(hrapEffPP);
  TH1D* relativeAcc_rap = (TH1D*)hrapAccAA->Clone("relAccAA_rap");
  relativeAcc_rap->Divide(hrapAccPP);
  hRAA_rap->Divide( relativeEff_rap ) ;  // efficiency correction
  hRAA_rap->Divide( relativeAcc_rap ) ;  // acceptance correction
  hRAA_rap->SetAxisRange(0,1.2,"Y");
  hRAA_rap->SetYTitle("R_{AA}");
  hRAA_rap->Draw();
  jumSun(0,1,2.4,1);

  
  cRAA_rap->SaveAs(Form("raa_vs_rap_%ds.pdf",state));

  //***Pt***
  TCanvas* c2 =  new TCanvas("c2","",400,400);
  hptSigAA = (TH1D*)  hptEffAA->Clone("hptSigAA");
  hptSigPP = (TH1D*)  hptEffPP->Clone("hptSigPP");
  hptSigAA->Reset();
  hptSigPP->Reset();  hptSigPP->SetYTitle("Raw yields");
    for ( int ipt = 1 ; ipt<= nPtBins ; ipt++) {
      valErr yieldPP = getYield(state, kPPDATA, ptBin[ipt-1], ptBin[ipt], yMin, yMax , centMin, centMax , 0, 100);
      valErr yieldAA = getYield(state, kAADATA, ptBin[ipt-1], ptBin[ipt], yMin, yMax , centMin, centMax , 0, 100);
      hptSigAA->SetBinContent( ipt, yieldAA.val ) ;
      hptSigAA->SetBinError( ipt, yieldAA.err ) ;
      hptSigPP->SetBinContent( ipt, yieldPP.val ) ;
      hptSigPP->SetBinError( ipt, yieldPP.err ) ;
    }
    
  //*****Pt Yield***** 
  hptSigPP->SetAxisRange(10,1e5,"Y");
  hptSigAA->SetAxisRange(10,1e5,"Y");
  //    cleverRange(hptSigPP[iy], 1.3, 1);
  handsomeTH1(hptSigPP,2);
  handsomeTH1(hptSigAA,2);
  hptSigPP->SetMarkerStyle(24);
  
  hptSigPP->Draw();
  gPad->SetLogy();
  hptSigAA->Draw("same");

  TCanvas* ccsPt = new TCanvas("crossSectionPt","",400,400);
  TH1D* hcsAA_pt = (TH1D*)hptSigAA->Clone("hcsAA_pt");
  TH1D* hcsPP_pt = (TH1D*)hptSigPP->Clone("hcsPP_pt");
  TH1D* hcsAA_rap = (TH1D*)hrapSigAA->Clone("hcsAA_rap");
  TH1D* hcsPP_rap = (TH1D*)hrapSigPP->Clone("hcsPP_rap");
  
  hcsAA_pt->Divide(hptEffAA); // efficiency 
  hcsPP_pt->Divide(hptEffPP);
  hcsAA_rap->Divide(hrapEffAA); 
  hcsPP_rap->Divide(hrapEffPP);
  cout << " here1 " << endl;
  hcsAA_pt->Divide(hptAccAA);  // acceptance
  cout << " here2 " << endl;
  hcsPP_pt->Divide(hptAccPP);
  cout << " here3 " << endl;
  hcsAA_rap->Divide(hrapAccAA);
  hcsPP_rap->Divide(hrapAccPP);

  
  // d_sigma 
  hcsPP_pt->Scale( 1. / 28000. )  ;     // pp : 25.8pb-1 = 25800 nb-1
  hcsPP_rap->Scale( 1. / 28000. )  ;     // pp : 25.8pb-1 = 25800 nb-1
  
  hcsAA_pt->Scale(1000000./TAA[nCentBins] );  // TAA : 5.607mb-1 = 5.607e-6 nb-1
  hcsAA_rap->Scale(1000000./TAA[nCentBins] ); // TAA : 5.607mb-1 = 5.607e-6 nb-1
  hcsAA_pt->Scale(1./NumberOfMBColl);
  hcsAA_rap->Scale(1./NumberOfMBColl);

//  hcsAA_pt->Scale( 1000. / 351. ) ;     // PbPb : 351 microBarn-1 = 0.351 nb-1
//  hcsAA_pt->Scale( 1./(208. * 208) );
//  hcsAA_rap->Scale( 1000. / 351. ) ;     // PbPb : 351 microBarn-1 = 0.351 nb-1
//  hcsAA_rap->Scale( 1./(208. * 208) );
  // d_sigma/dpT and d_sigma/dy
  TH1ScaleByWidth(hcsAA_pt);
  TH1ScaleByWidth(hcsPP_pt);
  cout << "noproblem" << endl;
  TH1ScaleByWidth(hcsAA_rap);
  TH1ScaleByWidth(hcsPP_rap);

  // pT cr:oss-section is normalized by delta y 
  hcsAA_pt->Scale( 0.5 / yMax );   // 1 / 4.8
  hcsPP_pt->Scale( 0.5 / yMax );   // 1 / 4.8
  // y cross-section.  We measured it in |y| bin
  hcsAA_rap->Scale( 0.5);
  hcsPP_rap->Scale( 0.5);
  
  handsomeTH1(hcsAA_pt,2);
  handsomeTH1(hcsPP_pt,1);
  hcsAA_pt->SetYTitle("B #times #frac{d#sigma}{dp_{T}} #frac{1}{#Deltay} [nb/(GeV/c)]");
  hcsAA_pt->SetAxisRange(  hcsAA_pt->GetBinContent(nPtBins) * 0.1,  hcsPP_pt->GetBinContent(1) * 10, "Y");
  hcsAA_pt->Draw();
  hcsPP_pt->Draw("same");
//  gPad->SetLogy();
  cout << "noproblem" << endl;
  
  TLegend* legCS = new TLegend(0.5046176,0.65,0.9,0.9,NULL,"brNDC");
  easyLeg(legCS,(Form("#Upsilon(%dS),  |y| < 2.4",state)) );
  legCS->AddEntry(hcsAA_pt, "PbPb #times A^{2}");
  legCS->AddEntry(hcsPP_pt, "pp");
  legCS->Draw();

  ccsPt->SaveAs(Form("crossSection_pt_%ds.pdf",state));
  TCanvas* ccsRap = new TCanvas("crossSectionRap","",400,400);
  handsomeTH1(hcsAA_rap,2);
  handsomeTH1(hcsPP_rap,1);
  hcsAA_rap->SetYTitle("B #times d#sigma/dy [nb]");
  //  hcsAA_rap->SetAxisRange(  hcsAA_rap->GetBinContent(nYBins) * 0.5,  hcsPP_rap->GetBinContent(1) * 1.5, "Y");
  hcsAA_rap->SetAxisRange(  0,  hcsPP_rap->GetBinContent(1) * 2.0, "Y");
  hcsAA_rap->Draw();
  hcsPP_rap->Draw("same");
  
  
  TLegend* legCSrap = new TLegend(0.5046176,0.65,0.9,0.9,NULL,"brNDC");
  easyLeg(legCSrap,(Form("#Upsilon(%dS)",state)) );
  legCSrap->AddEntry(hcsAA_pt, "PbPb #times A^{2}");
  legCSrap->AddEntry(hcsPP_pt, "pp");
  legCSrap->Draw();
  ccsRap->SaveAs(Form("crossSection_rap_%ds.pdf",state));



  TCanvas* cptRAA1 =  new TCanvas("cRAA_ptUnCorr","",400,400);
  cptRAA1->cd();
  hRAAraw_pt = (TH1D*)hptSigAA->Clone("hRAAraw_pt");
  hRAAraw_pt->Divide( hptSigPP );

  double scale_pt = getScale(nCentBins+1, TAA, centBin, nCentBins);
//  hRAAraw_pt->Scale( 26000000. / 351 ) ;   // pp : 26pb-1,  PbPb : 351 microBarn-1
//  hRAAraw_pt->Scale( 1./ (208.*208) );
  
  hRAAraw_pt->Scale(scale_pt);
  hRAAraw_pt->SetAxisRange(0,1.2,"Y");
  hRAAraw_pt->SetYTitle("R_{AA} (efficiency UNcorrected)");
  hRAAraw_pt->Draw();
  jumSun(0,1,30,1);
  
  TCanvas* cPtRAA =  new TCanvas("cRAA_pt","",400,400);
  hRAA_pt = (TH1D*)hRAAraw_pt->Clone("raa_vs_pt");

  TH1D* relativeEff = (TH1D*)hptEffAA->Clone("relEffAA");
  relativeEff->Divide(hptEffPP);
  TH1D* relativeAcc = (TH1D*)hptAccAA->Clone("relAccAA");
  relativeAcc->Divide(hptAccPP);
  hRAA_pt->Divide( relativeEff ) ; //efficiency correction
  cout << " here4 " << endl;
  hRAA_pt->Divide( relativeAcc ) ; // acceptance correction
  cout << " here5 " << endl;
  hRAA_pt->SetAxisRange(0,1.2,"Y");
  hRAA_pt->SetYTitle("R_{AA}");
  hRAA_pt->Draw();
  jumSun(0,1,30,1);
  
  cPtRAA->SaveAs(Form("raa_vs_pt_%ds.pdf",state));
  
 

  //******RAA vs Centrality *****
  TCanvas* cRAACent =  new TCanvas("cRAACent","",400,400);
  for(int icent=1; icent<=nCentBins;icent++)
  {
    valErr yCentAA = getYield(state,kAADATA,ptMin,ptMax,yMin,yMax, (int) centBin[icent-1],(int) centBin[icent],0,100); 
    valErr yCentPP = getYield(state,kPPDATA,ptMin,ptMax,yMin,yMax, 0, 200, 0,100);
    hSetBin -> SetBinContent(icent,(double)((centBin[icent]-centBin[icent-1])*nColl[icent-1]));
    hcentSigAA -> SetBinContent(icent,yCentAA.val);
    hcentSigAA -> SetBinError(icent,yCentAA.err);
  }
  valErr yCentAA_int = getYield(state,kAADATA,0,30,0,2.4,0,200,0,100);
  valErr yCentPP = getYield(state,kPPDATA,0,30,0,2.4,0,200,0,100);
  hcentSigAA_int -> SetBinContent(1,yCentAA_int.val);
  hcentSigAA_int -> SetBinError(1,yCentAA_int.err);
  hcentSigPP -> SetBinContent(1,yCentPP.val);
  hcentSigPP -> SetBinError(1,yCentPP.err);

  //****RAA CENTRALITY 0-100%****
  TH1D *hRAA_int;
  hRAA_int = (TH1D*) hcentSigAA_int->Clone("hRAAraw_centint");
  hRAA_int ->Divide(hcentSigPP);
  double scale_cent_int = getScale(nCentBins+1, TAA, centBin, nCentBins);
//  hRAA_int ->Scale( 26000000. / 351 ) ;   // pp : 26pb-1,  PbPb : 351 microBarn-1
//  hRAA_int ->Scale( 1./(208*208));
  hRAA_int -> Scale(scale_cent_int);

  // Ratio of efficiency 
  hcentEffAA_int->Divide( hcentEffPP );

  hRAA_int ->Divide(hcentEffAA_int);
  double AccIntPP = hIntAccPP->GetBinContent(1);
  double AccIntAA = hIntAccAA->GetBinContent(1);
  hRAA_int -> Scale(AccIntPP/AccIntAA);
  hRAA_int ->SetAxisRange(0,1.6,"Y");
  hRAA_int ->SetYTitle("R_{AA}");
  
  TGraphErrors *gRAA_int = new TGraphErrors();
  gRAA_int->SetName("gRAA_int");
  gRAA_int->SetTitle("raa_integrated");
  gRAA_int->SetPoint(0,1,hRAA_int->GetBinContent(1));
  gRAA_int->SetPointError(0,0,hRAA_int->GetBinError(1));



  
  hRAAraw_cent = (TH1D*) hcentSigAA->Clone("hRAAraw_cent");
  TH1D* htempPP = (TH1D*)hRAAraw_cent->Clone("htempPP_cent"); 
  TH1D* hRAAScaleCent = (TH1D*)hRAAraw_cent->Clone("hRAAScaleCent"); 
  htempPP->Reset();
  hRAAScaleCent->Reset();
  double scale_cent[nCentBins];
  for(int icent=1; icent<=nCentBins;icent++) {
    htempPP->SetBinContent( icent, yCentPP.val );
    htempPP->SetBinError( icent,0 );
    scale_cent[icent-1] = getScale(icent, TAA, centBin, nCentBins);
    hRAAScaleCent->SetBinContent(icent, 1./scale_cent[icent-1]);
  }
 
  cout << "PP int stat : " << yCentPP.err/yCentPP.val << endl;
  hRAAraw_cent->Divide(htempPP);
  hRAAraw_cent->Divide(hRAAScaleCent); 

//  hRAAraw_cent ->Divide(hSetBin);
//  hRAAraw_cent ->Scale( 26000000. / 351 ) ;   // pp : 26pb-1,  PbPb : 351 microBarn-1
//  hRAAraw_cent ->Scale( 1./(208*208));
//  hRAAraw_cent ->Scale(200.*392.);
  hRAAraw_cent ->SetAxisRange(0,1.2,"Y");
  hRAAraw_cent ->SetYTitle("R_{AA} (efficiency UNcorrected)");
  hRAAraw_cent ->Draw();
  jumSun(0,1,200,1);
  
   
  TCanvas* cRAACentEffCor =  new TCanvas("cRAACentEffCor","",400,400);
  //  TGraphErrors *gre1 = new TGraphErrors(nCentBins);
  //  gre1->SetName("raa_vs_npart");
  //  gre1->SetTitle("Graph");
  
  //  Int_t cii;      // for color index setting
  //  TColor *color; // for color definition with alpha
  //  cii = TColor::GetColor("#6699ff");
  //  gre1->SetFillColor(cii);
  //  gre1->SetMarkerStyle(20);
  //  for(int ibin=0;ibin<nCentBins;ibin++){
  //  gre1->SetPoint(ibin,nPart[ibin],hRAAraw_cent->GetBinContent(nCentBins-ibin));
  //  gre1->SetPointError(ibin,10,hRAAraw_cent->GetBinError(nCentBins-ibin));}
  //  gre1->Draw();
  
  for(int ibin = 1; ibin<nCentBins+1; ibin++)
    {
      cout << "yield at " << ibin<<"th bin: "<< hRAAraw_cent->GetBinContent(ibin) << endl;
    }
  
  
  hRAA_cent = (TH1D*) hRAAraw_cent ->Clone("hRAAraw_cent_final");
  TH1D* relativeEff_cent = (TH1D*) hcentEffAA -> Clone("relativeEff_cent");
  // Efficiency ratio : 
  double ppIntEff = hcentEffPP->GetBinContent ( 1 ) ; 
  for(int icent=1; icent<=nCentBins;icent++) {
    relativeEff_cent -> SetBinContent(  icent,   hcentEffAA->GetBinContent( icent ) / ppIntEff ) ;
  }
  
  hRAA_cent -> Divide(relativeEff_cent);
  hRAA_cent -> Scale(AccIntPP/AccIntAA);
  hRAA_cent -> SetAxisRange(0,1.2,"Y");
  hRAA_cent -> SetTitle("R_{AA}");
  
  TGraphErrors *gRAA_cent = new TGraphErrors(nCentBins);
  gRAA_cent->SetName("gRAA_cent");
  gRAA_cent->SetTitle("raa_vs_npart");

  Int_t ci;
  ci = TColor::GetColor("#6699ff");
  gRAA_cent->SetFillColor(ci);
  gRAA_cent->SetMarkerStyle(10);
  for(int ibin=0;ibin<nCentBins;ibin++){
  gRAA_cent->SetPoint(ibin,nPart[ibin],hRAA_cent->GetBinContent(nCentBins-ibin));
  gRAA_cent->SetPointError(ibin,0,hRAA_cent->GetBinError(nCentBins-ibin));}

  
  TPad *padl = new TPad("padl","padl", 0, 0., 0.9, 1);
  TPad *padr = new TPad("padr","padr", 0.9, 0., 1, 1);
//  padl->SetBottomMargin(0);
//  padr->SetBottomMargin(0);

  cRAACentEffCor->cd();  
  padl->Draw();
  padl->cd();
  TH1D *htemp = new TH1D("htemp",";N_{Part};RAA",420,0,420);
  TH1D *htempFull = new TH1D("htempfull","",1,0,2);
  handsomeTG1(gRAA_cent,2);
  htemp->SetAxisRange(0,1.6,"Y");
  htemp->SetYTitle("R_{AA}");
  handsomeTH1(htemp);
  htemp->DrawCopy();

  gRAA_cent->Draw("P same");

  drawText(Form("#Upsilon(%dS),  |y| < 2.4",state),0.45,0.87,1,16);
  drawText("p_{T}^{#mu} > 4 GeV", 0.45,0.80,1,16);
  jumSun(0,1,420,1);
  
  padr->SetFrameBorderMode(0);
  padr->SetBorderMode(0);
  padr->SetBorderSize(0);
  padr->SetTicks(0,0);
  padr->Draw();
  padr->cd();
  handsomeTG1(gRAA_int,2);
  htempFull->SetAxisRange(0,1.6,"Y");
  htempFull->GetXaxis()->SetLabelOffset(999);
  htempFull->GetXaxis()->SetLabelSize(0);
  htempFull->GetYaxis()->SetTickLength(0.);
  htempFull->GetXaxis()->SetTickLength(0.);
  handsomeTH1(htempFull);
  htempFull->DrawCopy();
  gRAA_int->GetHistogram()->GetXaxis()->SetLabelOffset(999);
  gRAA_int->GetHistogram()->GetXaxis()->SetLabelSize(0);
  gRAA_int->GetHistogram()->GetYaxis()->SetLimits(0,1.6);
  gRAA_int->GetHistogram()->GetYaxis()->SetRangeUser(0,1.6);
  gRAA_int->GetHistogram()->GetXaxis()->SetLimits(0,2);
  gRAA_int->GetHistogram()->GetXaxis()->SetRangeUser(0,2);
  gRAA_int->SetFillColor(ci);
  gRAA_int->SetMarkerStyle(20);
  gRAA_int->Draw("P same");

  cRAACentEffCor->SaveAs(Form("raa_vs_cent_%ds.pdf",state));

  //// Save TGraphErrors as root files
  // (convert TH1D to TGraph) 
  TGraphErrors *gRAA_pt = new TGraphErrors(hRAA_pt);
  gRAA_pt->SetName("gRAA_pt");
  gRAA_pt->SetTitle("raa_vs_pt");
  TGraphErrors *gRAA_rap = new TGraphErrors(hRAA_rap);
  gRAA_rap->SetName("gRAA_rap");
  gRAA_rap->SetTitle("raa_vs_rap");


  TGraphErrors *gCSPP_pt = new TGraphErrors(hcsPP_pt);
  gCSPP_pt->SetName("gCrossSection_pt_PP");
  gCSPP_pt->SetTitle("cs_vs_pt_PP");
  TGraphErrors *gCSAA_pt = new TGraphErrors(hcsAA_pt);
  gCSAA_pt->SetName("gCrossSection_pt_AA");
  gCSAA_pt->SetTitle("cs_vs_pt_AA");
  TGraphErrors *gCSAA_rap = new TGraphErrors(hcsAA_rap);
  gCSAA_rap->SetName("gCrossSection_rap_AA");
  gCSAA_rap->SetTitle("cs_vs_rap_AA");
  TGraphErrors *gCSPP_rap = new TGraphErrors(hcsPP_rap);
  gCSPP_rap->SetName("gCrossSection_rap_PP");
  gCSPP_rap->SetTitle("cs_vs_rap_PP");

  
  TFile *wf = new TFile(Form("finalResults/Ups_%d_RAA.root",state),"recreate");
  gRAA_int->Write(); // integrated
  gRAA_cent->Write(); // vs npart
  gRAA_pt->Write(); // vs pt
  gRAA_rap->Write(); //vs rap
  gCSPP_pt->Write();
  gCSAA_pt->Write();
  gCSPP_rap->Write();
  gCSAA_rap->Write();

  wf->Close();
  
}

valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh,
		float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow, yHigh, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High) ;
  TString SignalCB = "Double";
  TFile* inf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_RAA/upsilonRAA5TeV/fitResults/Final_NomResult_170124/PAS_fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  //TFile* inf = new TFile(Form("/home/samba/UpsilonAnalysis/fitResultFiles/mcFit_MuPt4_2016_11_04/fitresults_upsilon_%sCB_%s.root",SignalCB.Data(),kineLabel.Data()));
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret; 
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  cout << kineLabel << ": " << ret.val << " +/- " << ret.err << endl; 
  return ret;
}

void stripErrorBars( TH1* h, double defaultErr  ) {
  
  for ( int i=1;  i<= h->GetNbinsX() ; i++) {
    h->SetBinError( i, defaultErr);
  }
}

double getScale(int fTAA, double* TAA, double* centBin, int nCentBins)
{
  double flumi_;
  if(fTAA == nCentBins+1) flumi_ = 368;
  else if(centBin[fTAA-1]>=60 && centBin[fTAA-1]<120 && fTAA !=nCentBins+1) flumi_ = 464;
  else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 464;
  //else if(centBin[fTAA-1]>=120 && fTAA!=nCentBins+1) flumi_ = 334.82249848;
  else if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1) flumi_ = 368;
  double nMBColl = NumberOfMBColl;
  if(centBin[fTAA-1]<60 && fTAA !=nCentBins+1)  nMBColl = NumberOfMBColl;
  else if(centBin[fTAA-1]>=60 && fTAA !=nCentBins+1) nMBColl = NumberOfMBColl1;
  double scaleFactor;
  if(fTAA == nCentBins+1) scaleFactor = 28000000000000./(nMBColl*TAA[fTAA-1]*1000);
  else scaleFactor = 28000000000000./(nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200.*TAA[fTAA-1]*1000);
  
  if(fTAA!=nCentBins+1){
  cout << endl;
  cout << "icent : " << centBin[fTAA-1] << " - " << centBin[fTAA] << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "nMBColl*(centBin[fTAA]-centBin[fTAA-1])/100. : " << nMBColl*(centBin[fTAA]-centBin[fTAA-1])/200. << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }

  else if(fTAA==nCentBins+1){
  cout << endl;
  cout << "nMBColl : " << nMBColl << endl;
  cout << "TAA : " << TAA[fTAA-1] << endl;
  cout << "flumi : " << flumi_ << endl;
  cout << "scale : " << scaleFactor << endl;
  cout << endl;
  }
  return scaleFactor;
}

