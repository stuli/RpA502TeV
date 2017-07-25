#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
#include "../efficiency/tnp_weight.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");



void getMcSepctra(int state = 1, bool useDataWeight=true, bool useTnpWeight=true, int tnpIdx=0) {  // 1S, 2S, 3S
  TH1::SetDefaultSumw2();

TH1D* ncoll1SBin = new TH1D("ncoll1sbin","", nCentBins1s, centBin1s );
for ( int ii = 1 ; ii <= nCentBins1s ; ii++ )
  {
    ncoll1SBin->SetBinContent( ii , nColl1s[ii-1] ) ;
  }

 TString outputDirName = "mcSpectra";

  
  float massLow = 8;
  float massHigh = 13;
  
  int nPtBins=0;
  double* ptBin;
  int nCentBins=0;
  double* centBin;
  int nYBins=0;
  double *yBin;
  if ( state == 1 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
    nCentBins = nCentBins1s;  centBin = centBin1s;
  }
  else if ( state == 2 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nCentBins = nCentBins1s;  centBin = centBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
  }
  else if ( state == 3 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nCentBins = nCentBins1s;  centBin = centBin1s;
    nYBins = nYBins1S;  yBin = yBin1S;
  }

  
  double ptMin = ptBin[0];    double ptMax = ptBin[nPtBins];
  double yMin = yBin[0];    double yMax = yBin[nYBins];
  

  TH1D* hptGenAA;
  TH1D* hptRecoAA;
  TH1D* hptGenPP;
  TH1D* hptRecoPP;

  TH1D* hrapGenAA;
  TH1D* hrapRecoAA;
  TH1D* hrapGenPP;
  TH1D* hrapRecoPP;

  TH1D* hcentGenAA;
  TH1D* hcentRecoAA;

  TH1D* hcentintGenAA;
  TH1D* hcentintRecoAA;
  TH1D* hcentintGenPP;
  TH1D* hcentintRecoPP;
  //*~*~*~* for integrated bins *~*~*~*

  hcentintGenAA = new TH1D("hcentintGenAA","",1,0,200);
  hcentintRecoAA = (TH1D*) hcentintGenAA->Clone("hcentintRecoAA");
  hcentintGenPP = (TH1D*) hcentintGenAA->Clone("hcentintGenPP");
  hcentintRecoPP = (TH1D*) hcentintGenAA->Clone("hcentintRecoPP");
  //*~*~*~* for pt bins *~*~*~*
  hptGenAA = new TH1D("hptGenAA","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  hptRecoAA = (TH1D*)hptGenAA->Clone("hptRecoAA");
  hptGenPP = (TH1D*)hptGenAA->Clone("hptGenPP");
  hptRecoPP = (TH1D*)hptGenAA->Clone("hptRecoPP");
  //*~*~*~* for rapidity bins *~*~*~*
  hrapGenAA = new TH1D("hrapGenAA","; |y| ; ", nYBins, yBin);
  hrapRecoAA = (TH1D*)hrapGenAA->Clone("hrapRecoAA");
  hrapGenPP = (TH1D*)hrapGenAA->Clone("hrapGenPP");
  hrapRecoPP = (TH1D*)hrapGenAA->Clone("hrapRecoPP");
  //*~*~*~* for centrality bins *~*~*~*
  hcentGenAA = new TH1D("hcentGenAA","; centrality bin ; ", nCentBins, centBin);
  hcentRecoAA = (TH1D*)hcentGenAA->Clone("hcentRecoAA");

  
  TString fnameAA;
  TString fnamePP;
  if ( state == 1) { 
    fnamePP = "../skimmedFilesWeight2/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281226_.root"; 
    fnameAA = "../skimmedFilesWeight2/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281233_.root";  
  }
  else if ( state == 2 ) { 
    fnamePP = "../skimmedFilesWeight2/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281228_.root";
    fnameAA = "../skimmedFilesWeight2/yskimAA_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281234_.root";
  }
  else if ( state == 3 ) { 
    fnamePP = "../skimmedFilesWeight2/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281230_.root"; 
    fnameAA = "../skimmedFilesWeight2/yskimAA_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20161281235_.root"; 
  }
  
  TChain *mmAA = new TChain("mm");
  mmAA->Add(fnameAA);
  
  DiMuon dmAA;
  TBranch        *b_dmAA;
  mmAA->SetBranchAddress("mm",&dmAA,&b_dmAA);
  
  for(int iev=0; iev<mmAA->GetEntries() ; ++iev)  
    {
      mmAA->GetEntry(iev);
      //Accpetnace 
      if (  !( (dmAA.pt1 > glbMuPtCut) && (dmAA.pt2 > glbMuPtCut) && ( fabs(dmAA.eta1) < 2.4 ) &&  ( fabs(dmAA.eta2) <2.4 ) ) )
	continue;
      
      if (  !( (dmAA.mass > massLow) 
	       && (dmAA.mass < massHigh) 
	       && ( dmAA.pt > ptMin) 
	       && ( dmAA.pt < ptMax) 
	       && ( fabs(dmAA.y) > yMin)
	       && ( fabs(dmAA.y) < yMax)	      )
	    )
	continue;
      
      float ptWeight = dmAA.weight0; 
      if (useDataWeight) ptWeight = dmAA.weight;
      float ncollWeight =  ncoll1SBin->GetBinContent(  ncoll1SBin->FindBin( dmAA.cBin ) ) ;
      float tnpWeight = 1;
      if (useTnpWeight) { 
	if ( tnpIdx == 200 ) {  
          tnpWeight = tnp_weight_muid_pbpb(dmAA.pt1, dmAA.eta1) *  tnp_weight_muid_pbpb(dmAA.pt2, dmAA.eta2) ;
	}
	else if ( tnpIdx == 300 ) { 
          tnpWeight = tnp_weight_sta_pp(dmAA.pt1, dmAA.eta1) *  tnp_weight_sta_pp(dmAA.pt2, dmAA.eta2) ;
	}
	else   {
	  tnpWeight = tnp_weight_trg_pbpb(dmAA.pt1, dmAA.eta1, tnpIdx) *  tnp_weight_trg_pbpb(dmAA.pt2, dmAA.eta2, tnpIdx) ; 
	}
      }
      hcentintRecoAA->Fill( dmAA.cBin, ptWeight*tnpWeight*ncollWeight);
      hcentRecoAA->Fill   ( dmAA.cBin, ptWeight*tnpWeight*ncollWeight);
      hptRecoAA->Fill     ( dmAA.pt,   ptWeight*tnpWeight*ncollWeight);
      hrapRecoAA->Fill    ( dmAA.y,    ptWeight*tnpWeight*ncollWeight);
    }
  
  
  TChain *mmGenAA = new TChain("mmGen");   
  mmGenAA->Add(fnameAA);
  DiMuon dmGenAA; 
  TBranch        *b_dmGenAA;
  mmGenAA->SetBranchAddress("mmGen",&dmGenAA,&b_dmGenAA);
  
  for(int iev=0; iev<mmGenAA->GetEntries() ; ++iev)  
    {
      mmGenAA->GetEntry(iev);
      //Accpetnace 
      if (  !( (dmGenAA.pt1 > glbMuPtCut) && (dmGenAA.pt2 > glbMuPtCut) && ( fabs(dmGenAA.eta1) < 2.4 ) &&  ( fabs(dmGenAA.eta2) <2.4 ) ) )
	continue;
      
      if (  !( (dmGenAA.mass > massLow) 
	       && (dmGenAA.mass < massHigh) 
	       && ( dmGenAA.pt > ptMin) 
	       && ( dmGenAA.pt < ptMax) 
	       && ( fabs(dmGenAA.y) > yMin)
	       && ( fabs(dmGenAA.y) < yMax)	      )
	    )
	continue;
      
      float ptWeight = dmGenAA.weight0;
      if (useDataWeight) ptWeight = dmGenAA.weight;
      float ncollWeight =  ncoll1SBin->GetBinContent(  ncoll1SBin->FindBin( dmGenAA.cBin ) ) ;
      hcentintGenAA->Fill( dmGenAA.cBin, ptWeight*ncollWeight);
      hcentGenAA->Fill   ( dmGenAA.cBin, ptWeight*ncollWeight);
      hptGenAA->Fill     ( dmGenAA.pt,   ptWeight*ncollWeight);
      hrapGenAA->Fill    ( dmGenAA.y,    ptWeight*ncollWeight);
    }

  // pp

  TChain *mmPP = new TChain("mm");
  mmPP->Add(fnamePP);
  DiMuon dmPP;
  TBranch        *b_dmPP;
  mmPP->SetBranchAddress("mm",&dmPP,&b_dmPP);
  
  for(int iev=0; iev<mmPP->GetEntries() ; ++iev)  
    {
      mmPP->GetEntry(iev);
      //Accpetnace
      if (  !( (dmPP.pt1 > glbMuPtCut) && (dmPP.pt2 > glbMuPtCut) && ( fabs(dmPP.eta1) < 2.4 ) &&  ( fabs(dmPP.eta2) <2.4 ) ) )
        continue;

      if (  !( (dmPP.mass > massLow)
               && (dmPP.mass < massHigh)
               && ( dmPP.pt > ptMin)
               && ( dmPP.pt < ptMax)
               && ( fabs(dmPP.y) > yMin)
               && ( fabs(dmPP.y) < yMax)              )
            )
        continue;
      
      float ptWeight = dmPP.weight0;
      if (useDataWeight) ptWeight = dmPP.weight;
      float tnpWeight = 1; 
      if (useTnpWeight)   {
	if ( tnpIdx == 200 ) {  
          tnpWeight = tnp_weight_muid_pp(dmPP.pt1, dmPP.eta1) *  tnp_weight_muid_pp(dmPP.pt2, dmPP.eta2) ;
	}
	else if ( tnpIdx == 300 ) { 
          tnpWeight = tnp_weight_sta_pp(dmPP.pt1, dmPP.eta1) *  tnp_weight_sta_pp(dmPP.pt2, dmPP.eta2) ;
	}
	else   {
	  tnpWeight = tnp_weight_trg_pbpb(dmPP.pt1, dmPP.eta1, tnpIdx) *  tnp_weight_trg_pbpb(dmPP.pt2, dmPP.eta2, tnpIdx) ; 
	}
      }
      
      hcentintRecoPP->Fill( dmPP.cBin, ptWeight*tnpWeight);
      hptRecoPP->Fill     ( dmPP.pt,   ptWeight*tnpWeight);
      hrapRecoPP->Fill    ( dmPP.y,    ptWeight*tnpWeight);
    }
  
  
  TChain *mmGenPP = new TChain("mmGen");
  mmGenPP->Add(fnamePP);
  
  DiMuon dmGenPP;
  TBranch        *b_dmGenPP;
  mmGenPP->SetBranchAddress("mmGen",&dmGenPP,&b_dmGenPP);
  
  for(int iev=0; iev<mmGenPP->GetEntries() ; ++iev)  
    {
      mmGenPP->GetEntry(iev);
      //Accpetnace
      if (  !( (dmGenPP.pt1 > glbMuPtCut) && (dmGenPP.pt2 > glbMuPtCut) && ( fabs(dmGenPP.eta1) < 2.4 ) &&  ( fabs(dmGenPP.eta2) <2.4 ) ) )
	continue;

      if (  !( (dmGenPP.mass > massLow)
               && (dmGenPP.mass < massHigh)
               && ( dmGenPP.pt > ptMin)
	       && ( dmGenPP.pt < ptMax)
               && ( fabs(dmGenPP.y) > yMin)
               && ( fabs(dmGenPP.y) < yMax)              )
            )
        continue;
      
      float ptWeight = dmGenPP.weight0;
      if (useDataWeight) ptWeight = dmGenPP.weight;
      hcentintGenPP->Fill( dmGenPP.cBin, ptWeight);
      hptGenPP->Fill     ( dmGenPP.pt,   ptWeight);
      hrapGenPP->Fill    ( dmGenPP.y,    ptWeight);
    }
  
  
  // DRAW! 
  // Ncoll
  TCanvas* cnoll = new TCanvas("cnoll","",400,400);
  ncoll1SBin->Draw();

  // Pt Bins
  TCanvas* c_pt =  new TCanvas("c_pt","",600,600);
  c_pt->Divide(2,1);

  c_pt->cd(1);
  handsomeTH1(hptGenAA,1);
  handsomeTH1(hptRecoAA,2);
  cleverRange(hptGenAA, 1.3, 0);
  hptGenAA->SetYTitle("dN/dp_{T}");
  hptGenAA->Draw("hist");
  hptRecoAA->Draw("same");
  drawText("PbPb", 0.3, 0.8, 1, 15);
  drawText("|y| < 2.4", 0.3, 0.4, 1, 17);

  c_pt->cd(2);
  handsomeTH1(hptGenPP,1);
  handsomeTH1(hptRecoPP,2);
  cleverRange(hptGenPP, 1.3, 0);
  hptGenPP->SetYTitle("dN/dp_{T}");
  hptGenPP->Draw("hist");
  hptRecoPP->Draw("same");
  drawText("pp", 0.3, 0.8, 1, 15);
  drawText("|y| < 2.4", 0.3, 0.4, 1, 17);



  
  // Rap Bins
  TCanvas* c_rap =  new TCanvas("c_rap","",600,600);
  c_rap->Divide(2,1);

  c_rap->cd(1);
  handsomeTH1(hrapGenAA,1);
  handsomeTH1(hrapRecoAA,2);
  cleverRange(hrapGenAA, 1.3, 0);
  hrapGenAA->SetYTitle("dN/dp_{T}");
  hrapGenAA->Draw("hist");
  hrapRecoAA->Draw("same");
  drawText("PbPb", 0.3, 0.8, 1, 15);
  drawText("0 < p_{T} < 30 GeV/c", 0.3, 0.4, 1, 17);

  c_rap->cd(2);
  handsomeTH1(hrapGenPP,1);
  handsomeTH1(hrapRecoPP,2);
  cleverRange(hrapGenPP, 1.3, 0);
  hrapGenPP->SetYTitle("dN/dp_{T}");
  hrapGenPP->Draw("hist");
  hrapRecoPP->Draw("same");
  drawText("pp", 0.3, 0.8, 1, 15);
  drawText("0 < p_{T} < 30 GeV/c", 0.3, 0.4, 1, 17);


  // Centrality Bins
  TCanvas* c_cent =  new TCanvas("c_cent","",400,400);

  handsomeTH1(hcentGenAA,1);
  handsomeTH1(hcentRecoAA,2);
  handsomeTH1(hcentintGenAA,1);
  handsomeTH1(hcentintRecoAA,2);
  cleverRange(hcentGenAA, 1.3, 0);
  cleverRange(hcentintGenAA, 1.3, 0);
  hcentGenAA->SetYTitle("dN/dN_{part}");
  hcentintGenAA->SetYTitle("dN/dN_{part}");
  hcentGenAA->Draw("hist");
  hcentRecoAA->Draw("same");
  drawText("PbPb", 0.3, 0.8, 1, 15);
  drawText("|y| < 2.4", 0.3, 0.4, 1, 17);

  // *************************
  // ****** Efficiency *******
  // *************************

  // Efficiency Pt
  TCanvas* c_eff_pt =  new TCanvas("c_eff_pt","",400,400);
  TH1D* hptEffAA;
  TH1D* hptEffPP;
  c_eff_pt->cd();
  hptEffAA = (TH1D*)hptRecoAA->Clone("hptEffAA");
  hptEffAA ->Divide(hptGenAA);
  hptEffAA ->SetAxisRange(0,1.2,"Y");
  hptEffAA->SetYTitle("efficiency");
  hptEffAA->Draw();
  hptEffPP = (TH1D*)hptRecoPP->Clone("hptEffPP");
  hptEffPP->Divide(hptGenPP);
  hptEffPP->SetAxisRange(0,1.2,"Y");
  hptEffPP->SetYTitle("efficiency");
  hptEffPP->SetMarkerStyle(24);
  hptEffPP->Draw("same");
  TLegend* leg2 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg2,"");
  leg2->AddEntry(hptEffAA, "PbPb (0-100%)");
  leg2->AddEntry(hptEffPP, "pp");
  leg2->Draw();
  drawText(Form("#Upsilon(%dS),  p_{T}^{#mu} > 4GeV/c",state),0.25,0.87,1,15);
  jumSun(0,1,30,1);

  //Print the results for AN table 
  for ( int ii = 1 ; ii<= nPtBins ; ii++)   {
    if ( state == 1 ) {  
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << "$ \\GeVc &" <<  int(hptEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hptEffAA->GetBinContent(ii) *1000) / 1000. << " & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {  
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << "$ \\GeVc & & & " <<  int(hptEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hptEffAA->GetBinContent(ii) *1000) / 1000. << " & &   \\\\ " << endl;
    }    
    if ( state == 3 ) {  
      cout << "$" << ptBin[ii-1] << " < \\pt < " << ptBin[ii] << "$ \\GeVc & & & & & " <<  int(hptEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hptEffAA->GetBinContent(ii) *1000) / 1000. << " \\\\ " << endl;
    }    
  }
  // Efficiency Rap
  TCanvas* c_eff_rap =  new TCanvas("c_eff_rap","",400,400);
  TH1D* hrapEffAA;
  TH1D* hrapEffPP;
  c_eff_rap->cd();
  hrapEffAA = (TH1D*)hrapRecoAA->Clone("hrapEffAA");
  hrapEffAA ->Divide(hrapGenAA);
  hrapEffAA ->SetAxisRange(0,1.2,"Y");
  hrapEffAA->SetYTitle("efficiency");
  hrapEffAA->Draw();
  hrapEffPP = (TH1D*)hrapRecoPP->Clone("hrapEffPP");
  hrapEffPP->Divide(hrapGenPP);
  hrapEffPP->SetAxisRange(0,1.2,"Y");
  hrapEffPP->SetYTitle("efficiency");
  hrapEffPP->SetMarkerStyle(24);
  hrapEffPP->Draw("same");
  TLegend* leg3 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg3,"p_{T} < 30 GEV/c");
  leg3->AddEntry(hrapEffAA, "PbPb (0-100%)");
  leg3->AddEntry(hrapEffPP, "pp");
  leg3->Draw();
  jumSun(0,1,30,1);
  drawText(Form("#Upsilon(%dS),  p_{T}^{#mu} > 4GeV/c",state),0.25,0.87,1,15);

  // Print the results for the table in for AN
  for ( int ii = 1 ; ii<= nYBins ; ii++)   {
    if ( state == 1 ) {  
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ &" <<  int(hrapEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hrapEffAA->GetBinContent(ii) *1000) / 1000. << " & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {  
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ & & & " <<  int(hrapEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hrapEffAA->GetBinContent(ii) *1000) / 1000. << " & &   \\\\ " << endl;
    }    
    if ( state == 3 ) {  
      cout << "$" << yBin[ii-1] << " < |y| < " << yBin[ii] << "$ & & & & & " <<  int(hrapEffPP->GetBinContent(ii) *1000) / 1000. << " & " <<  int(hrapEffAA->GetBinContent(ii) *1000) / 1000. << " \\\\ " << endl;
    }    
  }

  
  // Centrality Efficiency
  TCanvas* c_eff_cent =  new TCanvas("c_eff_cent","",400,400);
  TH1D* hcentEffAA;
  TH1D* hcentintEffPP;
  TH1D* hcentintEffAA;
  hcentintEffAA = (TH1D*) hcentintRecoAA->Clone("hcentintEffAA");
  hcentintEffAA -> Divide(hcentintGenAA);
  hcentintEffAA -> SetAxisRange(0,1.2,"Y");
  hcentintEffAA -> SetYTitle("efficiency");
  hcentintEffPP = (TH1D*)hcentintRecoPP->Clone("hcentintEffPP");
  hcentintEffPP ->Divide(hcentintGenPP);
  hcentintEffPP ->SetMarkerStyle(24);
  hcentintEffPP ->SetLineStyle(2);
  hcentEffAA = (TH1D*)hcentRecoAA->Clone("hcentEffAA");
  hcentEffAA ->Divide(hcentGenAA);
  hcentEffAA ->SetAxisRange(0,1.2,"Y");
  hcentEffAA ->SetYTitle("efficiency");
  hcentEffAA ->Draw();
  hcentintEffAA -> Draw("same hist");
  hcentintEffPP ->Draw("same hist");

  TLegend* leg4 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg4,"|y| < 2.4");
  leg4->AddEntry(hcentintEffPP, "pp","l");
  leg4->AddEntry(hcentintEffAA, "PbPb (0-100%)","l");
  leg4->AddEntry(hcentEffAA, "PbPb","pl");
  leg4->Draw();
  drawText(Form("#Upsilon(%dS),  p_{T}^{#mu} > 4GeV/c",state),0.25,0.87,1,15);
  jumSun(0,1,200,1);


  // Print the results for the table in for AN
  for ( int ii = 1 ; ii<= nCentBins ; ii++)   {
    if ( state == 1 ) {  
      cout << "$" << centBin[ii-1]/2 << "\\% -- " << centBin[ii]/2 << "\\% $ &   & " <<  int(hcentEffAA->GetBinContent(ii) *1000) / 1000. << " & & & &   \\\\ " << endl;
    }
    if ( state == 2 ) {  
      cout << "$" << centBin[ii-1]/2 << "\\% -- " << centBin[ii]/2 << "\\% $ & & &   & " <<  int(hcentEffAA->GetBinContent(ii) *1000) / 1000. << " & &   \\\\ " << endl;
    }    
    if ( state == 3 ) {  
      cout << "$" << centBin[ii-1]/2 << "\\% -- " << centBin[ii]/2 << "\\% $ & & & & &   & " <<  int(hcentEffAA->GetBinContent(ii) *1000) / 1000. << " \\\\ " << endl;
    }    
  }

  // Integrated bin:
  cout << " pp = " << int(hcentintEffPP->GetBinContent(1)*1000)/1000. << ",  PbPb = " << int(hcentintEffAA->GetBinContent(1)*1000)/1000. << endl; 

  
  TFile* fout = new TFile(Form("%s/efficiency_ups%ds_useDataPtWeight%d_tnpWeight%d_tnpIdx%d.root",outputDirName.Data(), state,useDataWeight,useTnpWeight,tnpIdx),"recreate");

  c_eff_pt->SaveAs(Form("%s/eff_vs_pt_%ds_useDataPtWeight%d_tnpWeight%d_tnpIdx%d.pdf",outputDirName.Data(), state,useDataWeight,useTnpWeight,tnpIdx)) ;
  c_eff_rap->SaveAs(Form("%s/eff_vs_rap_%ds_useDataPtWeight%d_tnpWeight%d_tnpIdx%d.pdf",outputDirName.Data(), state,useDataWeight,useTnpWeight,tnpIdx)) ;
  c_eff_cent->SaveAs(Form("%s/eff_vs_cent_%ds_useDataPtWeight%d_tnpWeight%d_tnpIdx%d.pdf",outputDirName.Data(), state,useDataWeight,useTnpWeight,tnpIdx)) ;
  
  hptGenPP->Write();
  hptRecoPP->Write();
  hptGenAA->Write();
  hptRecoAA->Write();
  hrapGenPP->Write();
  hrapRecoPP->Write();
  hrapGenAA->Write();
  hrapRecoAA->Write();
  hptEffPP->Write();
  hptEffAA->Write();
  hrapEffPP->Write();
  hrapEffAA->Write();
  hcentintEffPP->Write();
  hcentEffAA->Write();
  hcentintEffAA->Write();
  
  fout->Close();


  
  
}
