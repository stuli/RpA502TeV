#ifndef CutAndBinCollection_C
#define CutAndBinCollection_C

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>


float glbMuPtCut = 4; // for acceptance
const int nPtBins1s  = 5;   double ptBin1s[nPtBins1s+1] = {0,2.5,5,8,15,30};
const int nPtBins2s  = 3;   double ptBin2s[nPtBins2s+1] = {0,5,15,30};
const int nPtBins3s  = 3;   double ptBin3s[nPtBins3s+1] = {0,5,15,30};

const int nYBins1S  = 6;   double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 2;   double yBin2S[nYBins2S+1] ={0, 1.2, 2.4};
const int nYBins3S  = 2;   double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};

const int nYBins  = 2;   double yBin[nYBins+1] ={0, 1.2, 2.4}; // for event reweighting

const int nCentBins1s  = 9;   double centBin1s[nCentBins1s+1] = {0,10,20,40,60,80,100,120,140,200}; 
const int nCentBins2s  = 9;   double centBin2s[nCentBins2s+1] = {0,10,20,40,60,80,100,120,140,200};
const int nCentBins3s  = 4;   double centBin3s[nCentBins3s+1] = {0,20,60,100,200};
double nPart1s[nCentBins1s]   = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4}; // HIN-16-008 paper
double nPart2s[nCentBins2s]   = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4};
double nPart3s[nCentBins3s]   = {21.9, 109.2, 226.7, 358.8}; 
double nColl1s[nCentBins1s]   = {1819,1432,1005,606,349,186,90.7,40.1,7.67}; 
double nColl2s[nCentBins2s]   = {1819,1432,1005,606,349,186,90.7,40.1,7.67}; 
double nColl3s[nCentBins3s]   = {1626, 805, 267, 30.8};

double TAA1s[nCentBins1s+1] = {25.98, 20.46, 14.35, 8.66, 4.978, 2.66, 1.296, 0.5729, 0.1095, 5.607};
double TAA2s[nCentBins2s+1] = {25.98, 20.46, 14.35, 8.66, 4.978, 2.66, 1.296, 0.5729, 0.1095, 5.607};
double TAA3s[nCentBins3s+1] = {23.22, 11.51, 3.819, 0.4395, 5.607};

const double inel_cross_PbPb = 7716;

//const int nCentBins2s  = 4;   double centBin2s[nCentBins2s+1] = {0,20,60,100,200};


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };

struct valErr { float val, err; } ; 

int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;
int kPPMCUps1S = 8 ;
int kPPMCUps2S = 9 ;
int kPPMCUps3S = 10 ;
int kAAMCUps1S = 11 ;
int kAAMCUps2S = 12 ;
int kAAMCUps3S = 13 ;
int kPPAADATASIMUL = 20 ; // 2 and 0 simultaneous fit
int kPPAADATAPeriSIMUL = 60 ; // 6 and 0 simultaneous fit

TString getCollID( int collid ) {
  if ( collid == kPPDATA ) return "PP_DATA";
  else if ( collid == kPADATA ) return "PA_DATA";
  else if ( collid == kAADATA ) return "AA_DATA";
  else if ( collid == kPPMC ) return "PP_MC";
  else if ( collid == kPAMC ) return "PA_MC";
  else if ( collid == kAAMC ) return "AA_MC";
//  else if ( collid == kAADATAPeri ) return "AA_DATA_PeriL1";
  else if ( collid == kAADATAPeri ) return "AA_DATA";
  else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
  else if ( collid == kPPMCUps1S ) return "PP_MC_Ups1S";
  else if ( collid == kPPMCUps2S ) return "PP_MC_Ups2S";
  else if ( collid == kPPMCUps3S ) return "PP_MC_Ups3S";
  else if ( collid == kAAMCUps1S ) return "AA_MC_Ups1S";
  else if ( collid == kAAMCUps2S ) return "AA_MC_Ups2S";
  else if ( collid == kAAMCUps3S ) return "AA_MC_Ups3S";
  else if ( collid == kPPAADATASIMUL ) return "PP_AA_DATA_SIMUL";
  else if ( collid == kPPAADATAPeriSIMUL ) return "PP_AA_DATA_PeriL1_SIMUL";

  else return "none";
}

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;


TString getEPSel( int eventPln) {
  if ( eventPln == kEPl2HF)  return "BothHFs";
  else if ( eventPln == kEPOppositeHF ) return "OppositeHF" ;
  else if ( eventPln == kEPSameSideHF ) return "SameSideHF" ;
  else return "none";
}


int kSoftMuCut = 0;
int kHighPtMuCut = 0;




class DiMuon {
 public:
 DiMuon() :
  run(0),   lumi(0), event(0), cBin(0), ep2(0), dphiEp2(0),
    vz(-99),  mass(-1), pt(-1), y(999), phi(999), eta(999),
    pt1(-1), eta1(-1), phi1(-1),        
    pt2(-1), eta2(-1), phi2(-1), weight0(0), weight(0),       
    oniaIndex(-1), softFlag(0), highPtFlag(0)
    {}
  
  int run;
  int lumi;
  int event;
  int cBin;
  float ep2;
  float dphiEp2;
  float vz;
  float mass;
  float pt;
  float y;
  float phi;    
  float eta;
  float pt1; 
  float eta1;
  float phi1;
  float pt2;
  float eta2;
  float phi2;    
  float weight0;
  float weight;
  int oniaIndex;
  int softFlag;
  int highPtFlag;
  
  void clear() {
    run = -99;  lumi=-99; event=-99; cBin=-99; ep2=-99, dphiEp2=-99; 
    vz=-99;     mass = -99; pt=-99; y=-99; phi=-99; eta=-99;      
    pt1=-99; eta1=-99; phi1=-99; pt2=-99; eta2=-99; phi2=-99; weight0=-99, weight=-99;
    oniaIndex=-1; softFlag=-1; highPtFlag=-1; 
  }

};
TString branchString = "run/I:lumi:event:cBin:ep2/F:dphiEp2:vz:mass:pt:y:phi:eta:pt1:eta1:phi1:pt2:eta2:phi2:weight0:weight:oniaIndex/I:softFlag:highPtFlag";



// Upsilon nominal bins
const int nPtBinsUps = 2;   double ptBinUps[nPtBinsUps+1] = {0, 5,     100};
const int  nYBinsUps = 2;   double yBinUps[nYBinsUps+1] =   {0, 1.2,   2.4};
const int nPBinsUps  = 3;   double pBinUps[nPBinsUps+1] =   {0, 0.167, 0.333,  0.5};


TString getKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh, float dphiEp2Low, float dphiEp2High) {
  TString kineLabel = Form("%s_pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, (float)muPtCut_) ;
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATAPeri ) || ( collId == kAADATACentL3) || (collId == kAAMCUps1S) || ( collId==kAAMCUps2S) || (collId == kAAMCUps3S) || (collId == kPPAADATASIMUL) || (collId == kPPAADATAPeriSIMUL))
    kineLabel = kineLabel+ Form("_centrality%d-%d_dphiEp_%.2fPI_%.2fPI",(int)cLow, (int)cHigh, (float)dphiEp2Low, (float)dphiEp2High ) ;
  return kineLabel;
}

#endif
