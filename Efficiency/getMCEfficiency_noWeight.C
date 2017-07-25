#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBin.h"
#include "../multiTreeUtil.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");

void setupMultiTreeTool( multiTreeUtil* mt=0, int UpsState=kPPMCUps1S, bool isGen=true) ;
void getMCEfficiency_noWeight(int state = 1) {  // 1S, 2S, 3S
  
  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample

  

  
  TCut accCut = Form("(pt1>%f) && (pt2>%f)", (float)glbMuPtCut, (float)glbMuPtCut);
  multiTreeUtil* genAA = new multiTreeUtil();
  multiTreeUtil* recoAA = new multiTreeUtil();
  multiTreeUtil* genPP = new multiTreeUtil();
  multiTreeUtil* recoPP = new multiTreeUtil();

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
    setupMultiTreeTool(genAA, kAAMCUps1S, true);  // isGen = 1
    setupMultiTreeTool(recoAA, kAAMCUps1S, false);  
    setupMultiTreeTool(genPP, kPPMCUps1S, true);  // isGen = 1
    setupMultiTreeTool(recoPP, kPPMCUps1S, false);  
  }
  else if ( state == 2 ) { 
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nCentBins = nCentBins2s;  centBin = centBin2s;
    nYBins = nYBins2S;  yBin = yBin2S;
    setupMultiTreeTool(genAA, kAAMCUps2S, true);  // isGen = 1
    setupMultiTreeTool(recoAA, kAAMCUps2S, false);  
    setupMultiTreeTool(genPP, kPPMCUps2S, true);  // isGen = 1
    setupMultiTreeTool(recoPP, kPPMCUps2S, false);  
  }
  else if ( state == 3 ) { 
    nPtBins = nPtBins3s;    ptBin = ptBin3s;
    nCentBins = nCentBins3s;  centBin = centBin3s;
    nYBins = nYBins3S;  yBin = yBin3S;
    setupMultiTreeTool(genAA, kAAMCUps3S, true);  // isGen = 1
    setupMultiTreeTool(recoAA, kAAMCUps3S, false);  
    setupMultiTreeTool(genPP, kPPMCUps3S, true);  // isGen = 1
    setupMultiTreeTool(recoPP, kPPMCUps3S, false);  
  }


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
  TH1D* hcentGenPP;
  TH1D* hcentRecoPP;

  TH1D* hcentintGenAA;
  TH1D* hcentintRecoAA;
  TH1D* hcentintGenPP;
  TH1D* hcentintRecoPP;
  
  bool doCent = true; 
  TCut PtCut = " pt >= 0 && pt < 30 ";
  TCut CbinCut = " cBin >= 0 && cBin <= 200 ";
  TCut yCut = " abs(y) >=0 && abs(y) < 2.4 ";

  
  //*~*~*~* for integrated bins *~*~*~*

  hcentintGenAA = new TH1D("hcentintGenAA","",1,0,200);
  hcentintRecoAA = (TH1D*) hcentintGenAA->Clone("hcentintRecoAA");
  hcentintGenPP = (TH1D*) hcentintGenAA->Clone("hcentintGenPP");
  hcentintRecoPP = (TH1D*) hcentintGenAA->Clone("hcentintRecoPP");

  genAA->Draw2(hcentintGenAA,"cBin", accCut && PtCut && CbinCut && yCut,"weight");
  recoAA->Draw2(hcentintRecoAA,"cBin", accCut && PtCut && CbinCut && yCut,"weight");
  genPP->Draw2(hcentintGenPP,"cBin", accCut && PtCut && CbinCut && yCut,"weight");
  recoPP->Draw2(hcentintRecoPP,"cBin", accCut && PtCut && CbinCut && yCut,"weight");

  
  //*~*~*~* for pt bins *~*~*~*
  
  hptGenAA = new TH1D("hptGenAA","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  hptRecoAA = (TH1D*)hptGenAA->Clone("hptRecoAA");
  hptGenPP = (TH1D*)hptGenAA->Clone("hptGenPP");
  hptRecoPP = (TH1D*)hptGenAA->Clone("hptRecoPP");

  genAA->Draw2( hptGenAA, "pt", accCut && yCut && CbinCut, "weight" ) ;
  recoAA->Draw2( hptRecoAA, "pt", accCut && yCut && CbinCut,"weight" ) ;
  genPP->Draw2( hptGenPP, "pt", accCut && yCut && CbinCut,"weight") ;
  recoPP->Draw2( hptRecoPP, "pt", accCut && yCut && CbinCut,"weight") ;


  //*~*~*~* for rapidity bins *~*~*~* 

  hrapGenAA = new TH1D("hrapGenAA","; |y| ; ", nYBins, yBin);
  hrapRecoAA = (TH1D*)hrapGenAA->Clone("hrapRecoAA");
  hrapGenPP = (TH1D*)hrapGenAA->Clone("hrapGenPP");
  hrapRecoPP = (TH1D*)hrapGenAA->Clone("hrapRecoPP");

  genAA->Draw2( hrapGenAA, "y",accCut && PtCut && CbinCut, "weight");
  recoAA->Draw2( hrapRecoAA, "y",accCut && PtCut && CbinCut, "weight");
  genPP->Draw2( hrapGenPP, "y",accCut && PtCut && CbinCut, "weight");
  recoPP->Draw2( hrapRecoPP, "y",accCut && PtCut && CbinCut, "weight");


  //*~*~*~* for centrality bins *~*~*~* 
  if ( doCent ) 
  { 
    hcentGenAA = new TH1D("hcentGenAA","; centrality bin ; ", nCentBins, centBin);
    hcentRecoAA = (TH1D*)hcentGenAA->Clone("hcentRecoAA");
    hcentGenPP = (TH1D*)hcentintGenAA->Clone("hcentGenPP");
    hcentRecoPP = (TH1D*)hcentintGenAA->Clone("hcentRecoPP");

    genAA->Draw2( hcentGenAA, "cBin", accCut && yCut && PtCut, "weight" ) ;
    recoAA->Draw2( hcentRecoAA, "cBin", accCut && yCut && PtCut,"weight" ) ;

    genPP->Draw2( hcentGenPP, "cBin", accCut && yCut && PtCut, "weight" ) ;
    recoPP->Draw2( hcentRecoPP, "cBin", accCut && yCut && PtCut,"weight" ) ;
  }


  // *************************
  // ****** Gen vs Reco ******
  // *************************
  
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
  TCanvas* c_cent =  new TCanvas("c_cent","",600,600);
  c_cent->Divide(2,1);
  
  c_cent->cd(1);
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
    
  c_cent->cd(2);
  handsomeTH1(hcentGenPP,1);
  handsomeTH1(hcentRecoPP,2);
  handsomeTH1(hcentintGenPP,1);
  handsomeTH1(hcentintRecoPP,2);
  cleverRange(hcentGenPP, 1.3, 0);
  cleverRange(hcentintGenPP, 1.3, 0);
  hcentGenPP->SetYTitle("dN/dN_{part}");
  hcentintGenPP->SetYTitle("dN/dN_{part}");
  hcentGenPP->Draw("hist");
  hcentRecoPP->Draw("same");
  drawText("pp", 0.3, 0.8, 1, 15);
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
  drawText("Accepted muon p_{T} > 4GeV/c",0.25,0.3,1,15); 
  drawText("|y|<2.4",0.25,0.5,1,15); 
  jumSun(0,1,30,1);


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
  easyLeg(leg3,"|y| < 2.4");
  leg3->AddEntry(hrapEffAA, "PbPb (0-100%)");
  leg3->AddEntry(hrapEffPP, "pp");
  leg3->Draw();
  drawText("Accepted muon p_{T} > 4GeV/c",0.25,0.2,1,15); 
  jumSun(0,1,30,1);

 
  // Centrality Efficiency 
  TCanvas* c_eff_cent =  new TCanvas("c_eff_cent","",400,400);
  TH1D* hcentEffAA;
  TH1D* hcentEffPP;
  TH1D* hcentEffAA_int;
  hcentEffAA_int = (TH1D*) hcentintRecoAA->Clone("hcentEffAA_int");
  hcentEffAA_int -> Divide(hcentintGenAA);
  hcentEffAA_int -> SetAxisRange(0,1.2,"Y");
  hcentEffAA_int -> SetYTitle("efficiency");
  hcentEffPP = (TH1D*)hcentRecoPP->Clone("hcentEffPP");
  hcentEffPP ->Divide(hcentGenPP);
  hcentEffPP ->SetMarkerStyle(24);
  hcentEffPP ->SetLineStyle(2);
  hcentEffAA = (TH1D*)hcentRecoAA->Clone("hcentEffAA");
  hcentEffAA ->Divide(hcentGenAA);
  hcentEffAA ->SetAxisRange(0,1.2,"Y");
  hcentEffAA ->SetYTitle("efficiency");
  hcentEffAA ->Draw();
  hcentEffAA_int -> Draw("same hist");
  hcentEffPP ->Draw("same hist");

  TLegend* leg4 = new TLegend(0.4046176,0.3500982,0.8492568,0.5304435,NULL,"brNDC");
  easyLeg(leg4,"|y| < 2.4");
  leg4->AddEntry(hcentEffPP, "pp","l");
  leg4->AddEntry(hcentEffAA_int, "PbPb (0-100%)","l");
  leg4->AddEntry(hcentEffAA, "PbPb","pl");
  leg4->Draw();
  drawText("Accepted muon p_{T} > 4GeV/c",0.25,0.2,1,15); 
  jumSun(0,1,200,1);

  
  TFile *fout = new TFile(Form("efficiency_ups%ds_MC_noWeight.root",state),"recreate");
  hptGenPP->Write();
  hptRecoPP->Write();
  hptEffPP->Write();
  hptGenAA->Write();
  hptRecoAA->Write();
  hptEffAA->Write();
  
  hrapGenPP->Write();
  hrapRecoPP->Write();
  hrapEffPP->Write();
  hrapGenAA->Write();
  hrapRecoAA->Write();
  hrapEffAA->Write();

  hcentEffPP->Write();
  hcentEffAA->Write();
  hcentEffAA_int->Write();
  fout->Close();


}

void setupMultiTreeTool( multiTreeUtil* mt, int UpsState, bool isGen) { 
  TString treeName = "mm";
  if ( isGen) treeName = "mmGen";

  if ( UpsState == kPPMCUps1S) {    
    //    mt->addFile("../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168142115_3c54df0419c4813e2d7256dc8952ac699405d027.root",treeName,""); }// pT weighted 
    mt->addFile(
"../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121653_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root"
,treeName,""); } 
  else if ( UpsState == kPPMCUps2S) { 
    mt->addFile(
		"../skimmedFiles/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121655_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root"
		,treeName,"");  }
  else if ( UpsState == kPPMCUps3S) { 
    mt->addFile(
		"../skimmedFiles/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121658_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root"
		,treeName,"");  }
  
  else if ( UpsState == kAAMCUps1S) { 
    mt->addFile(
		"../skimmedFiles/yskimAA_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016812170_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root",
		treeName,"");  } // Weighted
  else if ( UpsState == kAAMCUps2S) { 
    mt->addFile("../skimmedFiles/yskimAA_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016812171_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root",treeName,"");
  }
  else if ( UpsState == kAAMCUps3S) {
    mt->addFile("../skimmedFiles/yskimAA_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_2016812172_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root"
		,treeName,"");
  }
  
}


