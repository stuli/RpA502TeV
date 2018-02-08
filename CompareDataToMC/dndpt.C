#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "cutsAndBin.h"
#include <fstream>
//#include "SimplePtFit.C"
using namespace std;

valErr getYield(int state= 0, int collId=0, float ptLow=0, float ptHigh=0, float yLow=0, float yHigh=0, int cLow=0, int cHigh=0,   float dphiEp2Low=0,  float dphiEp2High=0) ;
void dndpt(int state= 1, int collId= kPPDATA, double getyieldfromPAS= 1, int WRITE=1) {

  cout<<"dN/dp_{T} starting "<< collId <<" state: "<< state <<"S, "<<"GetYield? : "<<getyieldfromPAS<<endl;

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(1);

  TH1::SetDefaultSumw2();

  TString fcollId;
  if(collId == kPPDATA) fcollId = "PP";
  else if(collId == kPADATA) fcollId = "PA";

  //// modify by hand according to the pt range of the sample
  int nPtBins=0;
  double* ptBin;
  int nPtBinsMC=0;
  double* ptBinMC;
  int nYBins=0;
  double* yBin;

  if ( state == 1 ) {
    nPtBins = nPtBins1s;    ptBin = ptBin1s;
    nPtBinsMC = nPtBins1sMC;    ptBinMC = ptBin1sMC;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  else if ( state == 2 ) {
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nPtBinsMC = nPtBins1sMC;    ptBinMC = ptBin1sMC;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  else if ( state == 3 ) {
    nPtBins = nPtBins2s;    ptBin = ptBin2s;
    nPtBinsMC = nPtBins1sMC;    ptBinMC = ptBin1sMC;
    nYBins = nYBins1S;    yBin = yBin1S;
  }
  //For 2D plot :
  float ptLo = -1.93; float ptHi = 1.93;
  float yLo  = -1.93; float  yHi = 1.93;

  // Get MC :
  float massLow = 8; float massHigh = 13;
  double ptMin = ptBinMC[0]; double ptMax = ptBinMC[nPtBinsMC];
  double yMin = yBin[0];     double yMax = yBin[nYBins];

  TH1D* hptData=new TH1D("hptData",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData1=new TH1D("hptData1",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData2=new TH1D("hptData2",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData3=new TH1D("hptData3",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData4=new TH1D("hptData4",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData5=new TH1D("hptData5",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData6=new TH1D("hptData6",";p_{T}(GeV/c);",nPtBins,ptBin);
  TH1D* hptData7=new TH1D("hptData7",";|y|;",nYBins,yBin);

  TH1D *DataFit1 = new TH1D("DataFit1","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *DataFit2 = new TH1D("DataFit2","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *DataFit3 = new TH1D("DataFit3","; p_{T} (GeV/c) ; ", nPtBins, ptBin);

  TH1D *WeightFactor = new TH1D("WeightFactor","; p_{T} (GeV/c) ; ", nPtBinsMC, ptBinMC);
  TH1D *WeightFactor_fit = new TH1D("WeightFactor_fit","; p_{T} (GeV/c) ; ", nPtBinsMC, ptBinMC);

  TH1D *hptMC  = new TH1D("hptMC","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC1 = new TH1D("hptMC1","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC2 = new TH1D("hptMC2","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC3 = new TH1D("hptMC3","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC4 = new TH1D("hptMC4","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC5 = new TH1D("hptMC5","; p_{T} (GeV/c) ; ", nPtBins, ptBin);
  TH1D *hptMC6 = new TH1D("hptMC6","; p_{T} (GeV/c) ; ",nPtBins,ptBin);
  TH1D *hptMC7 = new TH1D("hptMC7","; |y| ; ", nYBins, yBin);
  TH1D *MCfit  = new TH1D("MCfit","; p_{T} (GeV/c) ; ", nPtBinsMC, ptBinMC);

  TChain *mm = new TChain("mm");

  if(state==1){
    if(collId==kPPDATA){
      mm->Add("../skimmedFiles/yskimPP_MC_Ups1S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121653_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root");
    }
    else if(collId==kPADATA){
      mm->Add("../skimmedFiles/yskimPAMC_OpSign_20170727.root");
    }
  }
  else if(state==2){
    if(collId==kPPDATA){
      mm->Add("../skimmedFiles/yskimPP_MC_Ups2S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121655_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root");
    }
    else if(collId==kPADATA){
      mm->Add("../skimmedFiles/yskimPA2S_OpSign_20178161019_unIdentified.root");
    }
  }
  else if(state==3){
    if(collId==kPPDATA){
      mm->Add("../skimmedFiles/yskimPP_MC_Ups3S_Trig-L1DoubleMu0_OpSign_EP-OppositeHF_20168121658_ce8f82aaf8e612dcd1c7c161216161d988fbf9a6.root");
    }
    else if(collId==kPADATA){
      mm->Add("../skimmedFiles/yskimPA3S_L1DoubleMu0PD_Trig-L1DoubleMuOpen2016_OpSign_20178161028_unIdentified.root");
    }
  }

  DiMuon  dm;
  TBranch *b_dm;

  mm->SetBranchAddress("mm",&dm,&b_dm);

  for(int iev=0; iev<mm->GetEntries() ; ++iev)
  {
    mm->GetEntry(iev);

    if(collId==kPPDATA){
      if ( !(dm.pt1>4 && dm.pt2>4 && dm.y>yLo && dm.y<yHi && (fabs(dm.eta1)<2.4) &&  (fabs(dm.eta2)<2.4)))
      continue;
    }
    else if(collId==kPADATA){
      if ( !(dm.pt1>4 && dm.pt2>4 && dm.y>(yLo-0.47) && dm.y<(yHi-0.47) && (fabs(dm.eta1)<2.4) &&  (fabs(dm.eta2)<2.4)))
      continue;
    }
    //if (  !( (dm.mass > massLow)
    //      && (dm.mass < massHigh)
    //      && ( dm.pt > ptMin)
    //      && ( dm.pt < ptMax)
    //      //&& ( fabs(dm.y) > yMin)
    //      //&& ( fabs(dm.y) < yMax) 
    //      && ( dm.y > yMin)
    //      && ( dm.y < yMax) 
    //      )
    //   )
    //  continue;

    //if( (    dm.pt1>4 
    //      && dm.pt2>4 
    //      && (fabs(dm.eta1)<2.4) 
    //      && (fabs(dm.eta2)<2.4) 
    //      && (dm.mass > massLow)  
    //      && (dm.mass < massHigh)
    //      && ( dm.pt > ptMin)
    //      && ( dm.pt < ptMax)
    //      && ( fabs(dm.y) > yMin)
    //      && ( fabs(dm.y) < yMax)
    //      )  
    //    ){
   
    //if(dm.weight!=0){ 
    //cout<<dm.weight<<endl;
    //}
    if(collId==kPPDATA){
    hptMC->Fill      ( dm.pt );
    hptMC1->Fill     ( dm.pt );
    hptMC2->Fill     ( dm.pt );
    hptMC3->Fill     ( dm.pt );
    hptMC6->Fill     ( dm.pt );
    }
    else if(collId==kPADATA){
    hptMC->Fill      ( dm.pt );
    hptMC1->Fill     ( dm.pt );
    hptMC2->Fill     ( dm.pt );
    hptMC3->Fill     ( dm.pt );
    hptMC6->Fill     ( dm.pt );
    }
  }

  //------------------------------------------ Get Data :
  //FROM FINAL RESULTS
  if(getyieldfromPAS==0){
    double pxtmp, pytmp, extmp, eytmp;
    TFile* inf1 = new TFile(Form("../finalResults/Ups_%d_RAA.root",state));
    TGraphErrors* gCrossSection;
    if(collId==kPPDATA){
      gCrossSection = (TGraphErrors*)inf1->Get("gCrossSection_pt_PP");}
    else if(collId==kPADATA){
      gCrossSection = (TGraphErrors*)inf1->Get("gCrossSection_pt_PA");}
    cout<<"Loop Start (Get Data)"<<endl;
    for ( int ipt = 0 ; ipt<nPtBins ; ipt++) {
      gCrossSection->GetPoint(ipt, pxtmp, pytmp);
      extmp=gCrossSection->GetErrorX(ipt);
      eytmp=gCrossSection->GetErrorY(ipt);
      hptData->SetBinContent(ipt+1,pytmp);
      hptData->SetBinError(ipt+1,eytmp);
      hptData1->SetBinContent(ipt+1,pytmp);
      hptData1->SetBinError(ipt+1,eytmp);
      hptData2->SetBinContent(ipt+1,pytmp);
      hptData2->SetBinError(ipt+1,eytmp);
      hptData3->SetBinContent(ipt+1,pytmp);
      hptData3->SetBinError(ipt+1,eytmp);
      hptData4->SetBinContent(ipt+1,pytmp);
      hptData4->SetBinError(ipt+1,eytmp);
      hptData5->SetBinContent(ipt+1,pytmp);
      hptData5->SetBinError(ipt+1,eytmp);
      hptData6->SetBinContent(ipt+1,pytmp);
      hptData6->SetBinError(ipt+1,eytmp);
      cout<<"[ 1S-"<<ipt+1<<" ] :  "<<pytmp<<" +/- "<<eytmp<<endl;
    }
  }
  if(getyieldfromPAS==1){
    for(int ipt=1;ipt<=nPtBins;ipt++)
    {
      if(collId==kPPDATA){
        valErr yieldPP;
        yieldPP = getYield(state, 0, ptBin[ipt-1],ptBin[ipt], yLo, yHi, 0,200,0,100);
        hptData->SetBinContent(ipt,yieldPP.val);
        hptData->SetBinError(ipt,yieldPP.err);
        hptData1->SetBinContent(ipt,yieldPP.val);
        hptData1->SetBinError(ipt,yieldPP.err);
        //hptData2->SetBinContent(ipt,yieldPP.val);
        //hptData2->SetBinError(ipt,yieldPP.err);
        //hptData3->SetBinContent(ipt,yieldPP.val);
        //hptData3->SetBinError(ipt,yieldPP.err);
        hptData2->SetBinContent(ipt,(yieldPP.val)+(yieldPP.err));
        hptData3->SetBinContent(ipt,(yieldPP.val)-(yieldPP.err));
        hptData4->SetBinContent(ipt,yieldPP.val);
        hptData4->SetBinError(ipt,yieldPP.err);
        hptData5->SetBinContent(ipt,yieldPP.val);
        hptData5->SetBinError(ipt,yieldPP.err);
        hptData6->SetBinContent(ipt,yieldPP.val);
        hptData6->SetBinError(ipt,yieldPP.err);
        cout << Form(" PP %dS yield, pt  ",state) << ptBin[ipt-1] << " - " << ptBin[ipt] << " : " << yieldPP.val <<" +/- "<< yieldPP.err 
          <<", error : " << 100*yieldPP.err/yieldPP.val <<"%"<< endl;
      }
      else if(collId==kPADATA){
        valErr yieldPA;
        yieldPA = getYield(state, 1, ptBin[ipt-1],ptBin[ipt], yLo-0.47, yHi-0.47, 0,200,0,100);
        hptData->SetBinContent(ipt,yieldPA.val);
        hptData->SetBinError(ipt,yieldPA.err);
        hptData1->SetBinContent(ipt,yieldPA.val);
        hptData1->SetBinError(ipt,yieldPA.err);
        //hptData2->SetBinContent(ipt,yieldPA.val);
        //hptData2->SetBinError(ipt,yieldPA.err);
        //hptData3->SetBinContent(ipt,yieldPA.val);
        //hptData3->SetBinError(ipt,yieldPA.err);
        hptData2->SetBinContent(ipt,(yieldPA.val)+(yieldPA.err));
        hptData3->SetBinContent(ipt,(yieldPA.val)-(yieldPA.err));
        hptData4->SetBinContent(ipt,yieldPA.val);
        hptData4->SetBinError(ipt,yieldPA.err);
        hptData5->SetBinContent(ipt,yieldPA.val);
        hptData5->SetBinError(ipt,yieldPA.err);
        hptData6->SetBinContent(ipt,yieldPA.val);
        hptData6->SetBinError(ipt,yieldPA.err);
        cout << Form(" PA %dS yield, pt  ",state) << ptBin[ipt-1] << " - " << ptBin[ipt] << " : " << yieldPA.val <<" +/- "<< yieldPA.err 
          << ", error : " << 100*yieldPA.err/yieldPA.val <<"%"<< endl;
      }
    }
  }
  ///////////////////////////////////Normalization///////////////////////////////////
  hptMC->Scale(1./hptMC->Integral());
  hptMC1->Scale(1./hptMC1->Integral());
  hptMC2->Scale(1./hptMC2->Integral());
  hptMC3->Scale(1./hptMC3->Integral());
  hptMC4->Scale(1./hptMC4->Integral());
  hptMC5->Scale(1./hptMC5->Integral());
  hptMC6->Scale(1./hptMC6->Integral());
  hptMC7->Scale(1./hptMC7->Integral());

  hptData->Scale(1./hptData->Integral());
  hptData1->Scale(1./hptData1->Integral());
  hptData2->Scale(1./hptData2->Integral());
  hptData3->Scale(1./hptData3->Integral());
  hptData4->Scale(1./hptData4->Integral());
  hptData5->Scale(1./hptData5->Integral());
  hptData6->Scale(1./hptData6->Integral());

  //TH1ScaleByWidth(hptMC);   
  //TH1ScaleByWidth(hptMC1); 
  //TH1ScaleByWidth(hptMC2); 
  //TH1ScaleByWidth(hptMC3); 
  //TH1ScaleByWidth(hptMC4); 
  //TH1ScaleByWidth(hptMC5); 
  //TH1ScaleByWidth(hptMC6); 
  //TH1ScaleByWidth(hptMC7); 

  //TH1ScaleByWidth(hptData); 
  //TH1ScaleByWidth(hptData1);  
  //TH1ScaleByWidth(hptData2);  
  //TH1ScaleByWidth(hptData3);
  //TH1ScaleByWidth(hptData4);  
  //TH1ScaleByWidth(hptData5);  
  //TH1ScaleByWidth(hptData6);

  handsomeTH1(hptMC,1);     
  handsomeTH1(hptMC1,1);   
  handsomeTH1(hptData,1);   
  handsomeTH1(hptData1,1);  
  handsomeTH1(hptData2,1); 
  handsomeTH1(hptData3,1);
  handsomeTH1(hptData4,2);  
  handsomeTH1(hptData5,1); 
  handsomeTH1(hptData6,1);

  ///////////////////////////////////////////////////Fit Function///////////////////////////
  TF1* fitmc1;
  TF1* fitdata1;
  TF1* fitdata2;
  TF1* fitdata3;
  TF1* fitRatio1;
  TF1* fitRatio2;
  TF1* fitRatio3;
  TF1* fitRatio_Ap;
  TF1* fitRatio_Bp;
  TF1* fitRatio_Cp;
  TF1* fitRatio_Dp;
  TF1* fitRatio_Am;
  TF1* fitRatio_Bm;
  TF1* fitRatio_Cm;
  TF1* fitRatio_Dm;
  //1.Tsallis
  if( state == 1 ){
    if(collId==kPPDATA){
      fitmc1 = new TF1("fitmc1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[0]*x*TMath::Exp([1]*x+[2])",0,30);
      fitdata2 = new TF1("fitdata2","[0]*x*TMath::Exp([1]*x+[2])",0,30);
      fitdata3 = new TF1("fitdata3","[0]*x*TMath::Exp([1]*x+[2])",0,30);
      //fitRatio1 = new TF1("fitRatio1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))/([2]-1)*([3]-2)/([2]*[3]*([2]*[3] + ([2]-2)*[2])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([2]*[3])),-[2])",0,30);
      //fitRatio1 = new TF1("fitRatio1","([0] + [1]*x + [2]*x*x +[3]*x*x*x +[5]*x*x*x*x ) / ( (x-[4]) * (x-[4]) *(x-[4]) *(x-[6]) )",0,30);
      //fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x +[4]*x*x*x) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Cp = new TF1("fitRatio_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Dp = new TF1("fitRatio_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Cm = new TF1("fitRatio_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Dm = new TF1("fitRatio_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);

      fitRatio1->SetParameters( 2.55074e+02, -9.34016e+01, 4.42256e+01, -4.81048e+00 );
      //fitRatio1->SetParLimits( 0, 2.55074e+02+6.56645e+01,2.55074e+02-6.56645e+01);
    }
    else if(collId==kPADATA){
      fitmc1 = new TF1("fitmc1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*0.92 * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*0.92*x*0.92))/([0]*[1])),-[0]))",0,30);
      fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*1.18 * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*1.18*x*1.18))/([0]*[1])),-[0]))",0,30);
      //fitRatio1 = new TF1("fitRatio1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([0]*[1])),-[0]))/([2]-1)*([3]-2)/([2]*[3]*([2]*[3] + ([2]-2)*[2])) * x * TMath::Power(( 1+ (TMath::Sqrt(89.4916 + x*x))/([2]*[3])),-[2])",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Cp = new TF1("fitRatio_Cp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Dp = new TF1("fitRatio_Dp","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Cm = new TF1("fitRatio_Cm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio_Dm = new TF1("fitRatio_Dm","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);

      fitRatio1->SetParameters( 2.57960e+02,  -9.62878e+01,  4.38407e+01,  -4.73982e+00);

      //fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x )",0,30);
      //fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x )",0,30);
      //fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x )",0,30);
      //fitRatio1->SetParameters( 2.57960e+02,  -9.62878e+01);
    }
  }
  if( state == 2){
    if(collId==kPPDATA){
      fitmc1 = new TF1("fitmc1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*0.90 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*0.90*x*0.90))/([0]*[1])),-[0]))",0,30);
      fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*1.10 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*1.10*x*1.10))/([0]*[1])),-[0]))",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x )",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x )",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x )",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( ( [0] + [1]*x )  )",0,30);

    }
    else if(collId==kPADATA){
      fitmc1 = new TF1("fitmc1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*0.90 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*0.90*x*0.90))/([0]*[1])),-[0]))",0,30);
      fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*1.10 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*1.10*x*1.10))/([0]*[1])),-[0]))",0,30);
      //FitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      //FitRatio2 = new TF1("fitRatio2","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      //FitRatio3 = new TF1("fitRatio3","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x )  ",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x )  ",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x )  ",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( ( [0] + [1]*x )  )",0,30);

    }
  } 
  if( state == 3){
    if(collId==kPPDATA){
      fitmc1 = new TF1("fitmc1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*0.90 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*0.90*x*0.90))/([0]*[1])),-[0]))",0,30);
      fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*1.10 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*1.10*x*1.10))/([0]*[1])),-[0]))",0,30);
      //fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x* TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      //fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      //fitRatio1 = new TF1("fitRatio1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))/([2]-1)*([3]-2)/([2]*[3]*([2]*[3] + ([2]-2)*[2])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([2]*[3])),-[2])",0,30);
      //fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])  )",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x )",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x )",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x )",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( ( [0] + [1]*x )  )",0,30);
    }
    else if(collId==kPADATA){
      fitmc1 = new TF1("fitmc1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata1 = new TF1("fitdata1","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*0.90 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*0.90*x*0.90))/([0]*[1])),-[0]))",0,30);
      fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x*1.10 * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*1.10*x*1.10))/([0]*[1])),-[0]))",0,30);
      //fitdata2 = new TF1("fitdata2","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x* TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      //fitdata3 = new TF1("fitdata3","[2]*(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))",0,30);
      //fitRatio1 = new TF1("fitRatio1","(([0]-1)*([0]-2)/([0]*[1]*([0]*[1] + ([0]-2)*[0])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([0]*[1])),-[0]))/([2]-1)*([3]-2)/([2]*[3]*([2]*[3] + ([2]-2)*[2])) * x * TMath::Power(( 1+ (TMath::Sqrt(100.460529 + x*x))/([2]*[3])),-[2])",0,30);
      fitRatio1 = new TF1("fitRatio1","( [0] + [1]*x )",0,30);
      fitRatio2 = new TF1("fitRatio2","( [0] + [1]*x )",0,30);
      fitRatio3 = new TF1("fitRatio3","( [0] + [1]*x )",0,30);

      fitRatio_Ap = new TF1("fitRatio_Ap","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bp = new TF1("fitRatio_Bp","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Am = new TF1("fitRatio_Am","( ( [0] + [1]*x )  )",0,30);
      fitRatio_Bm = new TF1("fitRatio_Bm","( ( [0] + [1]*x )  )",0,30);
    }
  } 


  //////////////////////////////////////////////HISTOGRAM/////////////////////////////////////////////////////////
  //
  //TCanvas* c2 =  new TCanvas("c2","",800, 800);
  //c2->cd();

  TLegend *leg1 = new TLegend(0.65,0.75,0.85,0.85);
  leg1->AddEntry(hptData,Form("Data(%s %dS)",fcollId.Data(), state),"p");
  leg1->AddEntry(hptMC,Form("MC(%s %dS)",fcollId.Data(), state),"l");
  leg1->SetLineColor(kWhite);

  TLegend *leg2 = new TLegend(0.5,0.65,0.87,0.85);
  leg2->AddEntry(fitRatio1,Form("reweighting factor using fit(%s %dS)",fcollId.Data(), state),"l");
  //leg2->AddEntry(hptData1,Form("reweighting factor using bin by bin(%s %dS)",fcollId.Data(), state),"lp");
  leg2->AddEntry(fitRatio2,Form("reweighting factor using varied fit1(%s %dS)",fcollId.Data(), state),"l");
  leg2->AddEntry(fitRatio3,Form("reweighting factor using varied fit2(%s %dS)",fcollId.Data(), state),"l");
  leg2->SetLineColor(kWhite);

  TCanvas* c1 =  new TCanvas("c1","",800, 600);
  c1->Divide(1,2);
  c1->cd(1);
  hptData->Draw();
  hptMC->Draw("same hist");
  leg1->Draw("same");
  if(state==1){
  hptData->SetAxisRange(0.,0.5,"Y");
  }
  else if(state==2||state==3){
  hptData->SetAxisRange(0.,1.,"Y");
  }
  c1->cd(2);
  hptData1->Divide(hptMC1);
  hptData2->Divide(hptMC2);
  hptData3->Divide(hptMC3);
  hptData1->Fit(fitRatio1,"IE","",0,30);
  ofstream of1(Form("params_results_%s_%dS.txt",fcollId.Data(),state));
  TFitResultPtr r = hptData1->Fit(fitRatio1,"S");
  r.Get()->Print(cout);
    if(state==1){
  fitRatio2->FixParameter(0,fitRatio1->GetParameter(0)+fitRatio1->GetParError(0));
  fitRatio2->FixParameter(1,fitRatio1->GetParameter(1)+fitRatio1->GetParError(1));
  fitRatio2->FixParameter(2,fitRatio1->GetParameter(2)+fitRatio1->GetParError(2));
  fitRatio2->FixParameter(3,fitRatio1->GetParameter(3)+fitRatio1->GetParError(3));
  fitRatio3->FixParameter(0,fitRatio1->GetParameter(0)-fitRatio1->GetParError(0));
  fitRatio3->FixParameter(1,fitRatio1->GetParameter(1)-fitRatio1->GetParError(1));
  fitRatio3->FixParameter(2,fitRatio1->GetParameter(2)-fitRatio1->GetParError(2));
  fitRatio3->FixParameter(3,fitRatio1->GetParameter(3)-fitRatio1->GetParError(3));
  ////////////////////////To get variation////////////////////////////////////////////
  fitRatio_Ap->FixParameter(0,fitRatio1->GetParameter(0)+fitRatio1->GetParError(0));//
  fitRatio_Ap->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Ap->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Ap->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Am->FixParameter(0,fitRatio1->GetParameter(0)-fitRatio1->GetParError(0));//
  fitRatio_Am->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Am->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Am->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Bp->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Bp->FixParameter(1,fitRatio1->GetParameter(1)+fitRatio1->GetParError(1));//
  fitRatio_Bp->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Bp->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Bm->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Bm->FixParameter(1,fitRatio1->GetParameter(1)-fitRatio1->GetParError(1));//
  fitRatio_Bm->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Bm->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Cp->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Cp->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Cp->FixParameter(2,fitRatio1->GetParameter(2)+fitRatio1->GetParError(2));//
  fitRatio_Cp->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Cm->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Cm->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Cm->FixParameter(2,fitRatio1->GetParameter(2)-fitRatio1->GetParError(2));//
  fitRatio_Cm->FixParameter(3,fitRatio1->GetParameter(3));                          //
  fitRatio_Dp->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Dp->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Dp->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Dp->FixParameter(3,fitRatio1->GetParameter(3)+fitRatio1->GetParError(3));//
  fitRatio_Dm->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Dm->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Dm->FixParameter(2,fitRatio1->GetParameter(2));                          //
  fitRatio_Dm->FixParameter(3,fitRatio1->GetParameter(3)-fitRatio1->GetParError(3));//
  ////////////////////////////////////////////////////////////////////////////////////
  of1<<"\n"<<"Fit Parameter as Nomonal"<<endl;
  of1<<fitRatio1->GetParameter(0)<<", "<<fitRatio1->GetParameter(1)<<", "<<fitRatio1->GetParameter(2)<<", "<<fitRatio1->GetParameter(3)<<"\n"<<endl;
  of1<<"Fit Parameter + variation"<<endl;
  of1<<fitRatio2->GetParameter(0)<<", "<<fitRatio2->GetParameter(1)<<", "<<fitRatio2->GetParameter(2)<<", "<<fitRatio2->GetParameter(3)<<"\n"<<endl;
  of1<<"Fit Parameter - variation"<<endl;
  of1<<fitRatio3->GetParameter(0)<<", "<<fitRatio3->GetParameter(1)<<", "<<fitRatio3->GetParameter(2)<<", "<<fitRatio3->GetParameter(3)<<endl;
  }
  else if(state!=1){
  fitRatio2->FixParameter(0,fitRatio1->GetParameter(0)+fitRatio1->GetParError(0));
  fitRatio2->FixParameter(1,fitRatio1->GetParameter(1)+fitRatio1->GetParError(1));
  fitRatio3->FixParameter(0,fitRatio1->GetParameter(0)-fitRatio1->GetParError(0));
  fitRatio3->FixParameter(1,fitRatio1->GetParameter(1)-fitRatio1->GetParError(1));
  ////////////////////////To get variation////////////////////////////////////////////
  fitRatio_Ap->FixParameter(0,fitRatio1->GetParameter(0)+fitRatio1->GetParError(0));//
  fitRatio_Ap->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Am->FixParameter(0,fitRatio1->GetParameter(0)-fitRatio1->GetParError(0));//
  fitRatio_Am->FixParameter(1,fitRatio1->GetParameter(1));                          //
  fitRatio_Bp->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Bp->FixParameter(1,fitRatio1->GetParameter(1)+fitRatio1->GetParError(1));//
  fitRatio_Bm->FixParameter(0,fitRatio1->GetParameter(0));                          //
  fitRatio_Bm->FixParameter(1,fitRatio1->GetParameter(1)-fitRatio1->GetParError(1));//
  ////////////////////////////////////////////////////////////////////////////////////
  of1<<"\n"<<"Fit Parameter as Nomonal"<<endl;
  of1<<fitRatio1->GetParameter(0)<<", "<<fitRatio1->GetParameter(1)<<endl;
  of1<<"Fit Parameter + variation"<<endl;
  of1<<fitRatio2->GetParameter(0)<<", "<<fitRatio2->GetParameter(1)<<endl;
  of1<<"Fit Parameter - variation"<<endl;
  of1<<fitRatio3->GetParameter(0)<<", "<<fitRatio3->GetParameter(1)<<endl;
  }
  
  hptData1->Draw();
  hptData2->SetLineColor(kGreen+2);
  hptData3->SetLineColor(kRed+2);
  fitRatio1->SetLineColor(kGreen+2);
  fitRatio2->SetLineColor(kRed+2);
  fitRatio3->SetLineColor(kBlue+2);
  fitRatio1->Draw("same");
  fitRatio2->Draw("same");
  fitRatio3->Draw("same");
  hptData1->SetAxisRange(0.,5.,"Y");
  leg2->Draw("same");
  jumSun(0,1,30,1);
  //fitRatio1->Draw();
//  TCanvas* c1 =  new TCanvas("c1","",1200, 800);
//  c1->Divide(3,2);
//  //MC plot
//  c1->cd(1) ;
//  //hptMC->Fit(fitmc1,"I");
//  hptMC->Draw();
//  gPad->SetLogy();
//  //Data plot
//  cout<<"FIND!!!!!"<<endl;
//  c1->cd(2);
//  gPad->SetLogy();
//  hptData->Fit(fitdata1,"I");
//  hptData->Draw();
//  fitdata2->FixParameter(0,fitdata1->GetParameter(0));
//  fitdata2->FixParameter(1,fitdata1->GetParameter(1));
//  fitdata2->FixParameter(2,fitdata1->GetParameter(2));
//  fitdata3->FixParameter(0,fitdata1->GetParameter(0));
//  fitdata3->FixParameter(1,fitdata1->GetParameter(1));
//  fitdata3->FixParameter(2,fitdata1->GetParameter(2));
//  fitdata2->Draw("same");
//  fitdata3->Draw("same");
//  //hptData->SetAxisRange(0.,0.3,"Y");
//  //Weighting Factor
//  //handsomeTH1(WeightFactor1,1);   handsomeTH1(WeightFactor2,1); handsomeTH1(WeightFactor3,1);
//
//  c1->cd(3);
//  cout<<"FIND!!!!!"<<endl;
//  float fittmp1, fittmp2, fittmp3;
//  float mctmp1, mctmp2, mctmp3;
//  for(int i =0; i<nPtBinsMC; ++i){
//    fittmp1=fitdata1->Eval(hptMC1->GetBinCenter(i+1));
//    fittmp2=fitdata2->Eval(hptMC1->GetBinCenter(i+1));
//    fittmp3=fitdata3->Eval(hptMC1->GetBinCenter(i+1));
//    //mctmp1=fitmc1->Eval(hptMC1->GetBinCenter(i+1));
//    //fittmp1=fitdata1->Integral(ptBinMC[i],ptBinMC[i+1]);
//    //fittmp2=fitdata2->Integral(ptBinMC[i],ptBinMC[i+1]);
//    //fittmp3=fitdata3->Integral(ptBinMC[i],ptBinMC[i+1]);
//    WeightFactor1->SetBinContent(i+1,fittmp1);
//    WeightFactor2->SetBinContent(i+1,fittmp2);
//    WeightFactor3->SetBinContent(i+1,fittmp3);
//    //WeightDeno1->SetBinContent(i+1,mctmp1);
//  }
//  //WeightFactor1->Divide(WeightDeno1);
//  WeightFactor1->Divide(hptMC1);
//  WeightFactor2->Divide(hptMC2);
//  WeightFactor3->Divide(hptMC3);
//  WeightFactor1->SetMarkerColor(kBlack);     WeightFactor2->SetMarkerColor(kRed+1);     WeightFactor3->SetMarkerColor(kGreen+2);
//  WeightFactor1->SetMarkerStyle(1);          WeightFactor2->SetMarkerStyle(2);          WeightFactor3->SetMarkerStyle(2);
//  WeightFactor1->Draw();
//  WeightFactor2->Draw("same");
//  WeightFactor3->Draw("same");
//  WeightFactor1->SetAxisRange(0.,2.,"Y");
//  if(state==2){
//    WeightFactor1->SetAxisRange(0.,6.,"Y");}
//  //hptData6->Divide(hptMC6);
//  //hptData6->Draw();
//  jumSun(0,1,30,1);
//  //for(int i =0; i<nPtBins; ++i){
//  //WeightFactor1->SetBinContent( i+1, hptData6->GetBinContent(i+1));
//  //}
//  //WeightFactor1->Draw("same");
//  //Fit-Data Ratio
//  c1->cd(4);
//  hptMC4->Draw("hist");
//  hptData4->Draw("same");
//  hptMC4->SetAxisRange(0.,0.165,"Y");
//  jumSun(0,1,30,1);
//
//  c1->cd(5);
//  hptMC5->Draw("hist");
//  //hptData4->Draw("same");
//  hptData5->Draw("same");
//  hptMC5->SetAxisRange(0.,0.165,"Y");
//  jumSun(0,1,30,1);
//
//  c1->cd(6);
//  float fitval1, fitval2, fitval3;
//  float fitint1, fitint2, fitint3;
//  handsomeTH1(DataFit1,1);  handsomeTH1(DataFit2,1); handsomeTH1(DataFit3,1);
//  for(int i =0; i<nPtBins; i++){
//    fitval1=fitdata1->Eval(hptData1->GetBinCenter(i+1));
//    fitval2=fitdata2->Eval(hptData2->GetBinCenter(i+1));
//    fitval3=fitdata3->Eval(hptData3->GetBinCenter(i+1));
//    fitint1=fitdata1->Integral(ptBin[i],ptBin[i+1]);
//    fitint2=fitdata2->Integral(ptBin[i],ptBin[i+1]);
//    fitint3=fitdata3->Integral(ptBin[i],ptBin[i+1]);
//    DataFit1->SetBinContent(i+1,fitint1);
//    DataFit2->SetBinContent(i+1,fitint2);
//    DataFit3->SetBinContent(i+1,fitint3);
//  }
//  TH1ScaleByWidth(DataFit1);  TH1ScaleByWidth(DataFit2);  TH1ScaleByWidth(DataFit3);
//  for(int i=1; i<=nPtBins; i++){
//    cout<<"  ["<<ptBin[i-1]<<"-"<<ptBin[i]<<"] Data :"<<hptData1->GetBinContent(i)<<" Fit Int : "<<DataFit1->GetBinContent(i)
//      <<" , Ratio : "<<DataFit1->GetBinContent(i)/hptData1->GetBinContent(i)<<endl;
//    //", Fit val : "<<fitval1<<endl;
//  }
//  DataFit1->SetMarkerColor(kBlack); DataFit2->SetMarkerColor(kRed+1); DataFit3->SetMarkerColor(kGreen+2);
//  DataFit1->SetLineColor(kBlack); DataFit2->SetLineColor(kRed+1); DataFit3->SetLineColor(kGreen+2);
//  DataFit1->Divide(hptData1);
//  DataFit2->Divide(hptData2);
//  DataFit3->Divide(hptData3);
//  DataFit1->Draw();
//  DataFit2->Draw("same");
//  DataFit3->Draw("same");
//  DataFit1->SetAxisRange(0.,2.,"Y");
//  jumSun(0,1,30,1);
//
//
//
  int date = 20180131;
  c1->SaveAs(Form("dNdpt_plot_%s_%dS_%d.pdf",fcollId.Data(),state,date));
  c1->SaveAs(Form("dNdpt_plot_%s_%dS_%d.png",fcollId.Data(),state,date));

  if(WRITE==1){
    if(collId==kPPDATA){
      if(state==1){
        TFile *fUp1Spp = new TFile(Form("WeightedFcN_fit/ratioDataMC_PP_DATA_1s_%d.root",date),"RECREATE");
        fUp1Spp->cd();
      }
      else if(state==2){
        TFile *fUp2Spp = new TFile(Form("WeightedFcN_fit/ratioDataMC_PP_DATA_2s_%d.root",date),"RECREATE");
        fUp2Spp->cd();
      }
      else if(state==3){
        TFile *fUp2Spp = new TFile(Form("WeightedFcN_fit/ratioDataMC_PP_DATA_3s_%d.root",date),"RECREATE");
        fUp2Spp->cd();
      }
      hptData1->SetName("WeightFactor");
      hptData1->Write();
      fitRatio1->SetName("dataMC_Ratio_norm");
      fitRatio1->Write();
      fitRatio2->SetName("dataMC_Ratio_All_p");
      fitRatio2->Write();
      fitRatio3->SetName("dataMC_Ratio_All_m");
      fitRatio3->Write();
      fitRatio_Ap->SetName("dataMC_Ratio_Ap");
      fitRatio_Ap->Write();
      fitRatio_Bp->SetName("dataMC_Ratio_Bp");
      fitRatio_Bp->Write();
      fitRatio_Am->SetName("dataMC_Ratio_Am");
      fitRatio_Am->Write();
      fitRatio_Bm->SetName("dataMC_Ratio_Bm");
      fitRatio_Bm->Write();
      if(state==1){
      fitRatio_Cp->SetName("dataMC_Ratio_Cp");
      fitRatio_Cp->Write();
      fitRatio_Dp->SetName("dataMC_Ratio_Dp");
      fitRatio_Dp->Write();
      fitRatio_Cm->SetName("dataMC_Ratio_Cm");
      fitRatio_Cm->Write();
      fitRatio_Dm->SetName("dataMC_Ratio_Dm");
      fitRatio_Dm->Write();
      }
    }
    else if(collId==kPADATA){
      if(state==1){
        TFile *fUp1Spb = new TFile(Form("WeightedFcN_fit/ratioDataMC_PA_DATA_1s_%d.root",date),"RECREATE");
        fUp1Spb->cd();
      }
      else if(state==2){
        TFile *fUp2Spb = new TFile(Form("WeightedFcN_fit/ratioDataMC_PA_DATA_2s_%d.root",date),"RECREATE");
        fUp2Spb->cd();
      }
      else if(state==3){
        TFile *fUp2Spb = new TFile(Form("WeightedFcN_fit/ratioDataMC_PA_DATA_3s_%d.root",date),"RECREATE");
        fUp2Spb->cd();
      }
      hptData1->SetName("WeightFactor");
      hptData1->Write();
      fitRatio1->SetName("dataMC_Ratio_norm");
      fitRatio1->Write();
      fitRatio2->SetName("dataMC_Ratio_All_p");
      fitRatio2->Write();
      fitRatio3->SetName("dataMC_Ratio_All_m");
      fitRatio3->Write();
      fitRatio_Ap->SetName("dataMC_Ratio_Ap");
      fitRatio_Ap->Write();
      fitRatio_Bp->SetName("dataMC_Ratio_Bp");
      fitRatio_Bp->Write();
      fitRatio_Am->SetName("dataMC_Ratio_Am");
      fitRatio_Am->Write();
      fitRatio_Bm->SetName("dataMC_Ratio_Bm");
      fitRatio_Bm->Write();
      if(state==1){
      fitRatio_Cp->SetName("dataMC_Ratio_Cp");
      fitRatio_Cp->Write();
      fitRatio_Dp->SetName("dataMC_Ratio_Dp");
      fitRatio_Dp->Write();
      fitRatio_Cm->SetName("dataMC_Ratio_Cm");
      fitRatio_Cm->Write();
      fitRatio_Dm->SetName("dataMC_Ratio_Dm");
      fitRatio_Dm->Write();
      }
    }
  }
}

//Get Yield
valErr getYield(int state, int collId, float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh, float dphiEp2Low,  float dphiEp2High) {
  TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow-0.47, yHigh-0.47, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  //TString kineLabel = getKineLabel (collId, ptLow, ptHigh, yLow-0.47, yHigh-0.47, glbMuPtCut, cLow, cHigh, dphiEp2Low, dphiEp2High);
  TString SignalCB = "Double";
  //TFile* inf = new TFile(Form("/afs/cern.ch/work/g/goni/Analysis/Upsilon/upslion_RpA/src/Usercode/UpsilonpPb5TeV/Fitting/FixedParm/TEST/results/Fixed_fitresults_upsilon_%sCB_5TeV_%s.root",SignalCB.Data(),kineLabel.Data())); //Fixed Parameter
  //TFile* inf = new TFile(Form("/afs/cern.ch/work/g/goni/Analysis/Upsilon/upslion_RpA/src/Usercode/UpsilonpPb5TeV/Fitting/AllParmFree/results/nominal_170808/AllParmFree_fitresults_upsilon_%sCB_5TeV_%s.root",SignalCB.Data(),kineLabel.Data())); //Free Parameter
  TFile* inf = new TFile(Form("../Fitting/AllParmFree/Nominal/AllParmFree_fitresults_upsilon_%sCB_5TeV_%s.root",SignalCB.Data(),kineLabel.Data())); //Free Parameter
  cout<<kineLabel.Data()<<endl;
  RooWorkspace* ws = (RooWorkspace*)inf->Get("workspace");
  TH1D* fitResults = (TH1D*)inf->Get("fitResults");
  valErr ret;
  ret.val = fitResults->GetBinContent(state);
  ret.err = fitResults->GetBinError(state);
  //cout << kineLabel << ": " << " & " << ret.val << " $\pm$ " << ret.err << " & " <<ws->var("nBkg")->getVal() << " $\pm$ "<< ws->var("nBkg")->getError() << "\\\\" << endl;
  return ret;
}

