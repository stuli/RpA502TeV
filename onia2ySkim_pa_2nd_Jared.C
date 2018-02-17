#include <ctime>

#include "cutsAndBin.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include <TLorentzVector.h>
#include "TriggerManipulation.h" 
#include "commonUtility.h"
static const long MAXTREESIZE = 10000000000;

TString getDayAndTime();

Double_t fTsallis1SR(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y1S));
  Double_t mT = TMath::Sqrt(pdgMass.Y1S*pdgMass.Y1S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y1S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = c*xx*pow;


  Double_t c1 = (fpar[2]-1)*(fpar[2]-2)/(fpar[2]*fpar[3]*(fpar[2]*fpar[3]+(fpar[2]-2)*pdgMass.Y1S));
  Double_t mT1 = TMath::Sqrt(pdgMass.Y1S*pdgMass.Y1S+xx*xx);
  Double_t pow1 = TMath::Power((1+(mT1-pdgMass.Y1S)/(fpar[2]*fpar[3])),-fpar[2]);

  Double_t f1 = c1*xx*pow1;

  Double_t fr = f/f1;

                   
  return fr;
}

Double_t fTsallis2SR(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  Double_t c = (fpar[0]-1)*(fpar[0]-2)/(fpar[0]*fpar[1]*(fpar[0]*fpar[1]+(fpar[0]-2)*pdgMass.Y2S));
  Double_t mT = TMath::Sqrt(pdgMass.Y2S*pdgMass.Y2S+xx*xx);
  Double_t pow = TMath::Power((1+(mT-pdgMass.Y2S)/(fpar[0]*fpar[1])),-fpar[0]);

  Double_t f = c*xx*pow;


  Double_t c1 = (fpar[2]-1)*(fpar[2]-2)/(fpar[2]*fpar[3]*(fpar[2]*fpar[3]+(fpar[2]-2)*pdgMass.Y2S));
  Double_t mT1 = TMath::Sqrt(pdgMass.Y2S*pdgMass.Y2S+xx*xx);
  Double_t pow1 = TMath::Power((1+(mT1-pdgMass.Y2S)/(fpar[2]*fpar[3])),-fpar[2]);

  Double_t f1 = c1*xx*pow1;

  Double_t fr = f/f1;

                   
  return fr;
}


void onia2ySkim_pa_2nd_Jared( int nevt = -1,
		 int fileID = kPADATA,
		 int trigId=kL1DoubleMuOpen2016, 
		 int epSelection = kEPOppositeHF, 
		 bool saveTracks=false, 
		 TString skimVersion="unIdentified", 
		 bool DiMuSign = false
		 ) 
{
  using namespace std;
  
  bool isMC = false; 
  if ( (fileID == kPPMC) || (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) || (fileID == kAAMC) || (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )
    isMC = true;

  TChain *mytree = new TChain("myTree");

  TString fname;  TString fname1;  TString fname2;  TString fname3;  TString fname4;   TString fname5;

  const int nFiles = 6;
  double fileBin[nFiles+1] = {0,3,6,9,12,15,9999};
  const int nFiles3S = 4;
  double fileBin3S[nFiles3S+1] = {0,3,6,9,9999};
  TH1D* hWeight;
  if ( fileID == kAAMCUps3S )   {  
    hWeight = new TH1D("hWeight","hWeight",nFiles3S, fileBin3S); 
  }
  else { 
    hWeight = new TH1D("hWeight","hWeight",nFiles,   fileBin  );
  }
  
  TFile *inf_func;
  if(fileID == kPPMCUps1S) inf_func = new TFile("compareDataMc/ratioDataMC_PP_DATA_1sState.root","read");
  else if(fileID == kPPMCUps2S) inf_func = new TFile("compareDataMc/ratioDataMC_PP_DATA_2sState.root","read");
  else if(fileID == kPPMCUps3S) inf_func = new TFile("compareDataMc/ratioDataMC_PP_DATA_2sState.root","read");
  else if(fileID == kAAMCUps1S) inf_func = new TFile("compareDataMc/ratioDataMC_PA_DATA_1sState.root","read");
  else if(fileID == kAAMCUps2S) inf_func = new TFile("compareDataMc/ratioDataMC_PA_DATA_2sState.root","read");
  else if(fileID == kAAMCUps3S) inf_func = new TFile("compareDataMc/ratioDataMC_PA_DATA_2sState.root","read");

  TF1* fWgt = new TF1("fWgt","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )");

  if( fileID == kPPMCUps1S ) fWgt->SetParameters( 4.31213e+02,1.68660e+01,3.05209e+01,-5.91966e+00);
  else if( (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) ) fWgt->SetParameters( 3.31946e+03, 9.42285e+02, 2.04657e+02, -1.86672e+01);
  else if( fileID == kAAMCUps1S) fWgt->SetParameters(3.15619e+02,-8.65884e+01,4.65292e+01,-5.34732e+00);
  else if( (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) ) fWgt->SetParameters(2.28610e+02, 4.80360e+01, 5.29009e+01, -6.88728e+00);

  TF1* wFunc[nYBins+1];
  wFunc[1] = new TF1("weightCurve_1s",fTsallis1SR,0,30,4);
  wFunc[2] = new TF1("weightCurve_2s",fTsallis2SR,0,30,4);
  if ( (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) )  { 
    wFunc[1]->SetParameters( 0.988141, 3.0971, 1.81891, 10.0239);
    wFunc[2]->SetParameters(11.518, 7.53196, 2.38444, 2.68481);
  }
  else if ( (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )  { 
    wFunc[1]->SetParameters( 1.0001, 5.1, 2.0024, 12.4243);
    wFunc[2]->SetParameters( 3.46994, 11.8612, 2.10006, 3.25859);
  }
  
  //wFunc[1] = (TF1*) inf_func -> Get("dataMcRatio");
  //wFunc[1]  = new TF1("weightCurve_1s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2]))))",0,30); 
  //wFunc[2]  = new TF1("weightCurve_2s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2]))))",0,30); 
  
  if (fileID == kPPDATA) {
    fname = "/home/samba/OniaTree/Onia5TeV/ppData/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root";
    mytree->Add(fname.Data());
  }
  if (fileID == kPADATA) {
    //fname = "/home/samba.old/UpsilonAnalysis/tempfiles/Upsilon_pPb/data/RD2013_pa_2nd_run_merged.root"; // in lxplus
    fname = "/media/jared/Acer/Users/Jared/Desktop/Ubuntu_Overflow/RD2013_pa_2nd_run_merged.root";
    mytree->Add(fname.Data());
  }
  else if  (fileID == kAADATA) { 
    fname = "/home/samba/OniaTree/Onia5TeV/PbPbData/OniaTree_HIOniaL1DoubleMu0ABCD_HIRun2015-PromptReco-v1_Run_262620_263757.root";
    mytree->Add(fname.Data());
  }
  else if  (fileID == kAADATAPeri) { 
    fname = "/home/samba/OniaTree/Onia5TeV/PbPbData/OniaTree_HIOniaPeripheral30100_HIRun2015-PromptReco-v1_Run_262620_263757.root";
    mytree->Add(fname.Data());
  } 
  else if  (fileID == kAADATACentL3) {
    fname = "/home/samba/OniaTree/Onia5TeV/PbPbData/OniaTree_HIOniaCentral30L2L3_HIRun2015-PromptReco-v1_Run_262548_263757.root";  // Jan 29th
    mytree->Add(fname.Data());
  }
  else if  (fileID == kPPMCUps1S){
    fname = "/home/samba/OniaTree/Onia5TeV/ppOfficialMC/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root";
    mytree->Add(fname.Data());
  }
  else if  (fileID == kPPMCUps2S) {
    fname = "/home/samba/OniaTree/Onia5TeV/ppOfficialMC/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root";  
    mytree->Add(fname.Data());
  }
  else if  (fileID == kPPMCUps3S) {
    fname = "/home/samba/OniaTree/Onia5TeV/ppOfficialMC/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root";  
    mytree->Add(fname.Data());
  }
  else if  (fileID == kAAMCUps1S) {
    fname  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname1 = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname2 = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname3 = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname4 = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname5 = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups1SMM_ptUps_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    mytree->Add(fname.Data());    mytree->Add(fname1.Data());    mytree->Add(fname2.Data());    mytree->Add(fname3.Data());    mytree->Add(fname4.Data());     mytree->Add(fname5.Data());
    hWeight->SetBinContent(1,  3.10497);
    hWeight->SetBinContent(2,  4.11498);
    hWeight->SetBinContent(3,  2.2579);
    hWeight->SetBinContent(4,  1.2591);
    hWeight->SetBinContent(5,  0.567094);
    hWeight->SetBinContent(6,  0.783399);
  }
  else if  (fileID == kAAMCUps2S) {
    fname   = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname1  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname2  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname3  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname4  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname5  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups2SMM_ptUps2S_15_inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    mytree->Add(fname.Data());    mytree->Add(fname1.Data());    mytree->Add(fname2.Data());    mytree->Add(fname3.Data());    mytree->Add(fname4.Data());     mytree->Add(fname5.Data());
    hWeight->SetBinContent(1,  5.89168);
    hWeight->SetBinContent(2,  9.08207);
    hWeight->SetBinContent(3,  3.106);
    hWeight->SetBinContent(4,  1.10018);
    hWeight->SetBinContent(5,  0.534916);
    hWeight->SetBinContent(6,  0.776183);
  }
  else if  (fileID == kAAMCUps3S) {
    fname   = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups3SMM_ptUps3S_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname1  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups3SMM_ptUps3S_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname2  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups3SMM_ptUps3S_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    fname3  = "/home/samba/OniaTree/Onia5TeV/PbPbOfficialMC/OniaTree_Pythia8_Ups3SMM_ptUps3S_09_inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root";
    mytree->Add(fname.Data());    mytree->Add(fname1.Data());    mytree->Add(fname2.Data());    mytree->Add(fname3.Data());    
    hWeight->SetBinContent(1,  6.86815);
    hWeight->SetBinContent(2,  8.29618);
    hWeight->SetBinContent(3,  6.75153);
    hWeight->SetBinContent(4,  5.48684);
  }



  cout << endl << "*==*==*==*==*==*==*==* INPUT FILE *==*==*==*==*==*==*==*==*" << endl;

  TFile *f1 = new TFile(Form("%s",fname.Data()),"read");
  if (f1->IsZombie()) { cout << "*** KYO : CANNOT open the root file!! Macro terminated ***" << endl; return;} 
  cout <<"* file ::" << fname << endl;
  
  
  // Same or Opposite sign event
  TString fdimusign;
  if(!DiMuSign) fdimusign = "OpSign";
  else if(DiMuSign) fdimusign = "SSign";

  cout << endl;
  cout << "*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*" << endl;
  cout << " Sign of selecting dimuons : " << fdimusign.Data() << endl;
  cout << endl;

  

  // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
  TFile* newfile;
  if (fileID == kPPDATA) {
    newfile = new TFile(Form("Jared_yskimPP_L1DoubleMu0PD_%s_%s_%s.root",fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  if (fileID == kPADATA) {
    newfile = new TFile(Form("Jared_yskimPA2nd_%s_%s_%s.root",fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
   
   
  // import the tree to the RooDataSet
  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Ntracks;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_QQ_trig[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[200];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Ntracks;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_QQ_mupl_isHighPurity[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_isHighPurity[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_isHighPurity;   //!
  TBranch        *b_Reco_QQ_mumi_isHighPurity;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_isHighPurity", Reco_QQ_mupl_isHighPurity, &b_Reco_QQ_mupl_isHighPurity);
  mytree->SetBranchAddress("Reco_QQ_mumi_isHighPurity", Reco_QQ_mumi_isHighPurity, &b_Reco_QQ_mumi_isHighPurity);

  Bool_t          Reco_QQ_mupl_TrkMuArb[200];
  Bool_t          Reco_QQ_mumi_TrkMuArb[200];
  TBranch        *b_Reco_QQ_mupl_TrkMuArb;   //!
  TBranch        *b_Reco_QQ_mumi_TrkMuArb;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
  mytree->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);


  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  //  mytree->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  //  mytree->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_QQ_mupl_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_normChi2_global[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
  TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  mytree->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_QQ_mupl_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nMuValHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
  TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_QQ_mupl_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_StationsMatched[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
  TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  mytree->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_QQ_mupl_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dxy[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dxy;   //!
  TBranch        *b_Reco_QQ_mumi_dxy;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  mytree->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_QQ_mupl_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dz[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[200];   //[Reco_mu_size]
  Float_t         SumET_HFplusEta4;   //
  Float_t         SumET_HFminusEta4;   //
  Float_t         SumET_HFplus;   //
  Float_t         SumET_HFminus;   //
  Float_t         SumET_HF;   //
  TBranch        *b_Reco_QQ_mupl_dz;   //!
  TBranch        *b_Reco_QQ_mumi_dz;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  TBranch        *b_SumET_HFplusEta4;   //!
  TBranch        *b_SumET_HFminusEta4;   //!
  TBranch        *b_SumET_HFplus;   //!
  TBranch        *b_SumET_HFminus;   //!
  TBranch        *b_SumET_HF;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  mytree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  mytree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  mytree->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
  mytree->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
  mytree->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
  mytree->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  Int_t           Reco_QQ_mupl_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_QQ_mupl_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_mu_TMOneStaTight[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
  TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  mytree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_QQ_mupl_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nPixWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_QQ_mupl_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  Int_t           Reco_QQ_mumi_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  Float_t         Reco_QQ_mupl_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  Float_t         Reco_QQ_mumi_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  dimuon tree 
  ////////////////////////////////////////////////////////////////////////
  DiMuon dm;
  TTree *mmTree = new TTree("mm","diMuonPairs");
  TTree *hfTree = new TTree("hf","hf");
  mmTree->SetMaxTreeSize(MAXTREESIZE);
  mmTree->Branch("mm",&dm.run,branchString.Data());
  hfTree->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
  hfTree->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
  hfTree->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
  hfTree->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
  hfTree->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
  hfTree->Branch("Ntracks",&Ntracks,"Ntracks/I");
  
  /////////////////////////////////////////
  ////// Gen QQ 
  /////////////////////////////////////////
  Int_t           Gen_QQ_size;
  Int_t           Gen_QQ_type[200];   //[Gen_QQ_size]
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_type;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_mupl_4mom;   //!
  TBranch        *b_Gen_QQ_mumi_4mom;   //!
  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;

  if (isMC) { 
    mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    mytree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
    mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    mytree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
    mytree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
  }
  
  DiMuon dmGen;
  TTree *mmGenTree;
  if (isMC)  {
    mmGenTree = new TTree("mmGen","Gen Di-muon Pairs");
    mmGenTree->SetMaxTreeSize(MAXTREESIZE);
    mmGenTree->Branch("mmGen",&dmGen.run,branchString.Data());
  }
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","pt weight", 0, 10000,"");
  RooRealVar* hfpluseta4 = new RooRealVar("hfpluseta4","HF pluseta4", 0, 2000,"GeV/c^{2}");
  RooRealVar* hfminuseta4 = new RooRealVar("hfminuseta4","HF minuseta4", 0, 2000,"GeV/c^{2}");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var, *hfpluseta4, *hfminuseta4);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  Event selection tree 
  ////////////////////////////////////////////////////////////////////////
  TH1D* hEvent = new TH1D("hFilter","",20 ,0.5, 21.5);

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;
  
  
  int HLTPASS = 0;
  int DIMUTRIGPASS = 0;
  int DIMUIDPASS = 0;
  int DIMUPASS_all = 0;
  
  int DIMUPASS_CENT1 = 0;
  int DIMUPASS_CENT2 = 0;
  int DIMUPASS_CENT3 = 0;
  int DIMUPASS_CENT4 = 0;
  int DIMUPASS_CENT5 = 0;
  int DIMUPASS_CENT6 = 0;
  int DIMUPASS_CENT7 = 0;
  int DIMUPASS_CENT8 = 0;
  int DIMUPASS_CENT9 = 0;

  double ch_pluseta4,ch_minuseta4,ch_plus,ch_minus;
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;
    
    ///////////////////////////
    ///// Call the values /////
    ///////////////////////////
    mytree->GetEntry(iev);
        
    hEvent->GetXaxis()->SetBinLabel(1,"Events total");     hEvent->Fill(1);  

    if ( !((HLTriggers&1)==1))  
      continue; // trigger selection
    
    dmGen.clear();
    /////////////////////////////////////////////////////
    /////////  gen QQ loop first
    /////////////////////////////////////////////////////
    if (isMC) 
    { 
      //dmGen.clear();
      dmGen.run = runNb;
      dmGen.lumi = LS ;
      dmGen.event = eventNb ;
      dmGen.vz = zVtx;
      dmGen.cBin = Centrality ;
      for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) 
      {
        JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
        mupl_Gen = (TLorentzVector*) Gen_QQ_mupl_4mom->At(irqq);
        mumi_Gen = (TLorentzVector*) Gen_QQ_mumi_4mom->At(irqq);
        dmGen.mass   = JP_Gen->M();
        dmGen.pt     = JP_Gen->Pt();
        dmGen.phi    = JP_Gen->Phi();
        dmGen.y      = JP_Gen->Rapidity();
        dmGen.eta      = JP_Gen->Eta();
        dmGen.pt1  = mupl_Gen->Pt();
        dmGen.eta1 = mupl_Gen->Eta();
        dmGen.phi1 = mupl_Gen->Phi();
        dmGen.pt2  = mumi_Gen->Pt();
        dmGen.eta2 = mumi_Gen->Eta();
        dmGen.phi2 = mumi_Gen->Phi();

        if ( (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )   {
          dmGen.weight0 = (float) hWeight->GetBinContent( hWeight->FindBin(dmGen.pt) ) ;
        }
        else  {
          dmGen.weight0 = 1;  
        }
        // MC pT weight :  
        if ( (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) ||  (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S))  {
          if ( (fileID == kPPMCUps1S) || (fileID == kAAMCUps1S)) {
            dmGen.weight = dmGen.weight0 * fWgt->Eval(dmGen.pt); 
            //dmGen.weight = dmGen.weight0 * wFunc[1]->Eval(dmGen.pt); 
          }	  
          else {
      	    dmGen.weight = dmGen.weight0 * fWgt->Eval(dmGen.pt); 
      	    //dmGen.weight = dmGen.weight0 * wFunc[2]->Eval(dmGen.pt); 
          }
        }		 

        dmGen.oniaIndex = irqq;
        mmGenTree->Fill();
        
          /* 
        massVar->setVal( (double)dmGen.mass ) ;
        ptVar->setVal(   (double)dmGen.pt   ) ;
        yVar->setVal(    (double)dmGen.y    ) ;
        pt1Var->setVal(  (double)dmGen.pt1  ) ;
        eta1Var->setVal( (double)dmGen.eta1 ) ;
        pt2Var->setVal(  (double)dmGen.pt2  ) ;
        eta2Var->setVal( (double)dmGen.eta2 ) ;
        ep2Var->setVal( (double)dmGen.ep2 ) ;
        cBinVar->setVal( (double)dmGen.cBin ) ;
        dataSet->add( *argSet);
        */
      } //end of Gen QQ tree
    } //end of isMC condition for Gen QQ tree
    

    HLTPASS++;

    dm.clear();
    dm.run = runNb;
    dm.lumi = LS ;
    dm.event = eventNb ;
    dm.vz = zVtx;
    dm.cBin = Centrality ;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      hEvent->GetXaxis()->SetBinLabel(4,"Di-muons Total");      hEvent->Fill(4);
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      
      if ( !((Reco_QQ_trig[irqq]&1)==1) )
        continue; // trigger selection
      hEvent->GetXaxis()->SetBinLabel(7,"Di-muons trig");      hEvent->Fill(7);

      DIMUTRIGPASS++;


      if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta() ) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta()) )   )
        continue;      
      hEvent->GetXaxis()->SetBinLabel(5,"Di-muons Accep");      hEvent->Fill(5);

      bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) && (Reco_QQ_mupl_TrkMuArb[irqq] == true) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mupl_nPixWMea[irqq] > 0) &&
          (Reco_QQ_mupl_dxy[irqq]<0.3) &&
          (Reco_QQ_mupl_dz[irqq]<20.)	 &&  (Reco_QQ_mupl_isHighPurity[irqq]==true) 
          ) ; 

      bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) && (Reco_QQ_mumi_TrkMuArb[irqq] = true) && 
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mumi_nPixWMea[irqq] > 0) &&
          (Reco_QQ_mumi_dxy[irqq]<0.3) &&
          (Reco_QQ_mumi_dz[irqq]<20.)  && (Reco_QQ_mumi_isHighPurity[irqq]==true)
          ) ; 

      bool muplHighPtCut = ( (Reco_QQ_mupl_nMuValHits[irqq]>0) &&
          (Reco_QQ_mupl_StationsMatched[irqq]>1) &&
          (Reco_QQ_mupl_ptErr_global[irqq]/mupl_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mupl_dxy[irqq]<0.2) &&
          (Reco_QQ_mupl_dz[irqq] <0.5) &&
          (Reco_QQ_mupl_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) 
          );
      bool mumiHighPtCut = ( (Reco_QQ_mumi_nMuValHits[irqq]>0) &&
          (Reco_QQ_mumi_StationsMatched[irqq]>1) &&
          (Reco_QQ_mumi_ptErr_global[irqq]/mumi_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mumi_dxy[irqq]<0.2) &&
          (Reco_QQ_mumi_dz[irqq] <0.5) &&
          (Reco_QQ_mumi_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) 
          );


      if ( !(muplSoft && mumiSoft) ) 
        continue;   
      hEvent->GetXaxis()->SetBinLabel(8,"Di-muons mu ID");      hEvent->Fill(8);

      DIMUIDPASS++;

      // vertex probablitiy cut:
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;
      hEvent->GetXaxis()->SetBinLabel(6,"Di-muons Vtx prob.");      hEvent->Fill(6);

      if( !DiMuSign )
      {
        if ( Reco_QQ_sign[irqq] != 0 ) continue;
      }
      else if(DiMuSign)
      {
        if( Reco_QQ_sign[irqq] == 0) continue;
      }

      hEvent->GetXaxis()->SetBinLabel(9,"Di-muoons charge sign");      hEvent->Fill(9);


      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();
      dm.phi    = JP_Reco->Phi();
      dm.y      = JP_Reco->Rapidity();
      dm.eta      = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();

      dm.oniaIndex = irqq;
      mmTree->Fill();

      massVar->setVal( (double)dm.mass ) ;
      ptVar->setVal(   (double)dm.pt   ) ;
      yVar->setVal(    (double)dm.y    ) ;
      pt1Var->setVal(  (double)dm.pt1  ) ;
      dphiEp2Var->setVal(   (double)dm.dphiEp2  ) ;
      eta1Var->setVal( (double)dm.eta1 ) ;
      pt2Var->setVal(  (double)dm.pt2  ) ;
      eta2Var->setVal( (double)dm.eta2 ) ;
      ep2Var->setVal( (double)dm.ep2 ) ;
      cBinVar->setVal( (double)dm.cBin ) ;
      evtWeight->setVal( (double)dm.weight ) ;
      hfpluseta4->setVal( (double) SumET_HFminusEta4) ;
      hfminuseta4->setVal( (double) SumET_HFplusEta4) ;
      dataSet->add( *argSet);
    } // end of dimuon loop
     
      ch_pluseta4 = SumET_HFplusEta4; 
      
      ch_minuseta4 = SumET_HFminusEta4; 
      ch_plus = SumET_HFplus; 
      ch_minus = SumET_HFminus;

      SumET_HFplusEta4 = ch_minuseta4; 
      SumET_HFminusEta4 = ch_pluseta4; 
      SumET_HFplus = ch_minus; 
      SumET_HFminus = ch_plus; 
      hfTree->Fill();
  } //end of event loop
  
  
  cout << endl;
  cout << endl;
  cout << "******************" << endl;
  cout << "**** SUMMARY *****" << endl;
  cout << endl;

  cout << "1. # of events passing HLT : " << HLTPASS << endl;
  cout << "2. # of dimuons passing trigger matching after pass 1 : " << DIMUTRIGPASS << endl;
  cout << "3. # of dimuons passing muon ID cut after pass 1-2 & acceptance cut : " << DIMUIDPASS << endl;
  cout << "******************" << endl;
  cout << "******************" << endl;
  cout << endl;
  cout << endl;

  dataSet->Write();
  mmTree->Write();  // Don't need to call Write() for trees
  hEvent->Write();
  hfTree->Write();
  newfile->Close();
} 



TString getDayAndTime() 
{ 
  time_t currentTime;
  struct tm *localTime;
  
  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time
  
  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%d%d%d%d%d",Year,Month,Day,Hour,Min);
}
