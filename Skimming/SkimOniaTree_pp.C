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

bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) ;
void SkimOniaTree_pp( int nevt = -1,
		 int fileID = kPPDATA,
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

  TChain *mytree = new TChain("hionia/myTree");
  /*  TChain *trkTree;
  if ( (fileID == kPPDATA) || (fileID == kPADATA) || (fileID == kPPMC) || (fileID == kPAMC) || (fileID == kPPMCUps1S) || (fileID==kPPMCUps2S) || (fileID == kPPMCUps3S) )
  trkTree = new TChain("ppTrack/trackTree");
  else 
  trkTree = new TChain("anaTrack/trackTree"); */

  TString fname;  TString fname1;  TString fname2;  TString fname3;  TString fname4;   TString fname5;

  const int nFiles = 6;
  double fileBin[nFiles+1] = {0,3,6,9,12,15,9999};
  const int nFiles3S = 4;
  double fileBin3S[nFiles3S+1] = {0,3,6,9,9999};

if (isMC) {
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
  else if(fileID == kAAMCUps1S) inf_func = new TFile("compareDataMc/ratioDataMC_AA_DATA_1sState.root","read");
  else if(fileID == kAAMCUps2S) inf_func = new TFile("compareDataMc/ratioDataMC_AA_DATA_2sState.root","read");
  else if(fileID == kAAMCUps3S) inf_func = new TFile("compareDataMc/ratioDataMC_AA_DATA_2sState.root","read");



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
}
  
  //wFunc[1] = (TF1*) inf_func -> Get("dataMcRatio");
  //wFunc[1]  = new TF1("weightCurve_1s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2]))))",0,30); 
  //wFunc[2]  = new TF1("weightCurve_2s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2]))))",0,30); 
  
  if (fileID == kPPDATA) {
    fname = "/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PP_Data/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root";
    mytree->Add(fname.Data());
  }
  if (fileID == kPADATA) {
    fname = "/afs/cern.ch/user/k/kimy/workDir/oniaTree/OniaTree_HIOnia_HIRun2016_ExpressStream_Run_285480-285539.root"; // in lxplus
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
  //  if (f1->IsZombie()) { cout << "*** KYO : CANNOT open the root file!! Macro terminated ***" << endl; return;} 
  cout <<"* file ::" << fname << endl;
  
  
  // Same or Opposite sign event
  TString fdimusign;
  if(!DiMuSign) fdimusign = "OpSign";
  else if(DiMuSign) fdimusign = "SSign";

  cout << endl;
  cout << "*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*" << endl;
  cout << " Sign of selecting dimuons : " << fdimusign.Data() << endl;
  cout << endl;

  

  // *==*==*==*==*==*==* Trigger selection *==*==*==*==*==*==* //
  TString trigName = getTrig(trigId);
  hltIndex hltBits = getTrigIndex(trigId, fname);
  
  // *==*==*==*==*==*==* Event Plane  *==*==*==*==*==*==* //
  TString epName = getEPSel(epSelection) ;
  cout << " Event Plane : " << epName << endl;
  
  
  // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
  TFile* newfile;
  if (fileID == kPPDATA) {
    newfile = new TFile(Form("SkimmedFiles_Santona/yskimPP_L1DoubleMu0PD_Trig-%s_%s_%s_%s.root",trigName.Data(), fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  if (fileID == kPADATA) {
    newfile = new TFile(Form("skimmedFiles/yskimPA_Trig-%s_%s_%s_%s.root",trigName.Data(), fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  else if (fileID == kAADATA) { 
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_L1DoubleMu0PD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if  (fileID == kAADATACentL3) {
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_CentralPD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if  (fileID == kAADATAPeri) { 
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_PeripheralPD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMC) {
    newfile = new TFile(Form("skimmedFiles/yskimPP_MC_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMC) {
    newfile = new TFile(Form("skimmedFiles/yskimPbPb_MC_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate");
  }
  else if (fileID == kPPMCUps1S) {
    newfile = new TFile(Form("skimmedFiles/yskimPP_MC_Ups1S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMCUps2S) {
    newfile = new TFile(Form("skimmedFiles/yskimPP_MC_Ups2S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMCUps3S) {
    newfile = new TFile(Form("skimmedFiles/yskimPP_MC_Ups3S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps1S) {
    newfile = new TFile(Form("skimmedFiles/yskimAA_MC_Ups1S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps2S) {
    newfile = new TFile(Form("skimmedFiles/yskimAA_MC_Ups2S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps3S) {
    newfile = new TFile(Form("skimmedFiles/yskimAA_MC_Ups3S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
   
   
  // import the tree to the RooDataSet
  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
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
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_QQ_mupl_highPurity[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_highPurity[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_highPurity;   //!
  TBranch        *b_Reco_QQ_mumi_highPurity;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
  mytree->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);


  
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
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
  TBranch        *b_Reco_QQ_mupl_dz;   //!
  TBranch        *b_Reco_QQ_mumi_dz;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  mytree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  mytree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
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
  
  
  


  
  /*  const int NMAXTRK = 10000;
  Int_t           nTrk;
  Float_t         trkPt[NMAXTRK];   //[nTrk]
  Float_t         trkPtError[NMAXTRK];   //[nTrk]
  UChar_t         trkNHit[NMAXTRK];   //[nTrk]
  Float_t         trkEta[NMAXTRK];   //[nTrk]
  Float_t         trkPhi[NMAXTRK];   //[nTrk]
  Int_t           trkCharge[NMAXTRK];   //[nTrk]
  Bool_t          highPurity[NMAXTRK];   //[nTrk]
  Bool_t          tight[NMAXTRK];   //[nTrk]
  Bool_t          loose[NMAXTRK];   //[nTrk]
  Float_t         trkChi2[NMAXTRK];   //[nTrk]
  UChar_t         trkNdof[NMAXTRK];   //[nTrk]
  Float_t         trkDxy1[NMAXTRK];   //[nTrk]
  Float_t         trkDxyError1[NMAXTRK];   //[nTrk]
  Float_t         trkDz1[NMAXTRK];   //[nTrk]
  Float_t         trkDzError1[NMAXTRK];   //[nTrk]
  Bool_t          trkFake[NMAXTRK];   //[nTrk]
  UChar_t         trkAlgo[NMAXTRK];   //[nTrk]
  Float_t         trkMVA[NMAXTRK];   //[nTrk]
  Float_t         pfEcal[NMAXTRK]; //[nTrk]
  Float_t         pfHcal[NMAXTRK]; //[nTrk]
  
  TBranch        *b_nTrk;   //!
  TBranch        *b_trkPt;   //!
  TBranch        *b_trkPtError;   //!
  TBranch        *b_trkNHit;   //!
  TBranch        *b_trkEta;   //!
  TBranch        *b_trkPhi;   //!
  TBranch        *b_trkCharge;   //!
  TBranch        *b_highPurity;   //!
  TBranch        *b_tight;   //!
  TBranch        *b_loose;   //!
  TBranch        *b_trkChi2;   //!
  TBranch        *b_trkNdof;   //!
  TBranch        *b_trkDxy1;   //!
  TBranch        *b_trkDxyError1;   //!
  TBranch        *b_trkDz1;   //!
  TBranch        *b_trkDzError1;   //!
  TBranch        *b_trkFake;   //!
  TBranch        *b_trkAlgo;   //!
  TBranch        *b_trkMVA;   //!
  TBranch        *b_pfEcal;   //!
  TBranch        *b_pfHcal;   //!
  if (saveTracks) { 
    trkTree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
    trkTree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
    trkTree->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
    trkTree->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
    trkTree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
    trkTree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
    trkTree->SetBranchAddress("trkCharge", trkCharge, &b_trkCharge);
    trkTree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
    trkTree->SetBranchAddress("tight", tight, &b_tight);
    trkTree->SetBranchAddress("loose", loose, &b_loose);
    trkTree->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
    trkTree->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
    trkTree->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
    trkTree->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
    trkTree->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
    trkTree->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
    trkTree->SetBranchAddress("trkFake", trkFake, &b_trkFake);
    trkTree->SetBranchAddress("trkAlgo", trkAlgo, &b_trkAlgo);
    trkTree->SetBranchAddress("trkMVA", trkMVA, &b_trkMVA);
    trkTree->SetBranchAddress("pfEcal", pfEcal, &b_pfEcal);
    trkTree->SetBranchAddress("pfHcal", pfHcal, &b_pfHcal);
    } */

  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  dimuon tree 
  ////////////////////////////////////////////////////////////////////////
  DiMuon dm;
  TTree *mmTree = new TTree("mm","diMuonPairs");
  mmTree->SetMaxTreeSize(MAXTREESIZE);
  mmTree->Branch("mm",&dm.run,branchString.Data());
  
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
  //  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *ep2Var, *pt1Var, *eta1Var, *pt2Var, *eta2Var);
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  if ( (fileID == kAAMC) || (fileID == kAADATA) || (fileID == kAADATAPeri) || (fileID == kAADATACentL3) || (fileID == kPAMC) || (fileID == kPADATA) || (fileID==kAAMCUps1S) || (fileID==kAAMCUps2S)|| (fileID==kAAMCUps3S) )
  {argSet->add(*cBinVar);argSet->add(*ep2Var); argSet->add(*dphiEp2Var);}

  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  Gen QQ tree        
  ////////////////////////////////////////////////////////////////////////
  DiMuon dmGen;
  TTree *mmGenTree;
  if (isMC)  {
    mmGenTree = new TTree("mmGen","Gen Di-muon Pairs");
    mmGenTree->SetMaxTreeSize(MAXTREESIZE);
    mmGenTree->Branch("mmGen",&dmGen.run,branchString.Data());
  }
    
  /*
////////////////////////////////////////////////////////////////////////
  //////////////////  Track tree        
  ////////////////////////////////////////////////////////////////////////
  TTree *trkOutTree = new TTree("track","trackTree");
  static const int MAXTRK  = 10000;   // This must be  enough.
  int nTrkOut;
  float trkOutPt[MAXTRK];
  float trkOutEta[MAXTRK];
  float trkOutPhi[MAXTRK];
  float trkOutDphi[MAXTRK];
  float trkOutDeta[MAXTRK];
  float trkOutDr[MAXTRK];
  float trkPiTripleMass[MAXTRK];
  float trkKaTripleMass[MAXTRK];
  float trkOutCorr[MAXTRK];
  trkOutTree->SetMaxTreeSize(MAXTREESIZE);
  trkOutTree->Branch("nTrack",&nTrkOut,"nTrack/I");
  trkOutTree->Branch("pt",trkOutPt,"pt[nTrack]/F");
  trkOutTree->Branch("eta",trkOutEta,"eta[nTrack]/F");
  trkOutTree->Branch("phi",trkOutPhi,"phi[nTrack]/F");
  trkOutTree->Branch("dphi",trkOutDphi,"dphi[nTrack]/F");
  trkOutTree->Branch("deta",trkOutDeta,"deta[nTrack]/F");
  trkOutTree->Branch("dr",trkOutDr,"dr[nTrack]/F");
  trkOutTree->Branch("corr",trkOutCorr,"corr[nTrack]/F");
  trkOutTree->Branch("massPiTriple",trkPiTripleMass,"massPiTriple[nTrack]/F");
  trkOutTree->Branch("massKaTriple",trkKaTripleMass,"massKaTriple[nTrack]/F");
  */  
  
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
    
    ///////////////////////////
    ///// Vertex cut //////////
    ///////////////////////////
    //    if (TMath::Abs(zVtx) > 15.) continue;
    //    hEvent->GetXaxis()->SetBinLabel(2,"Events z vertex < 15cm");    hEvent->Fill(2);
    
    // Trigger cut is now moved into the RECO_QQ loop..
   
    if ( isTrigMatched(hltBits,HLTriggers) == false ) 
      continue; // trigger selection

    HLTPASS++;

    hEvent->GetXaxis()->SetBinLabel(3,Form("Events Trig %s",trigName.Data()));    hEvent->Fill(3);

    dm.clear();
    dm.run = runNb;
    dm.lumi = LS ;
    dm.event = eventNb ;
    dm.vz = zVtx;
    dm.cBin = Centrality ;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      hEvent->GetXaxis()->SetBinLabel(4,"Di-muons Total");      hEvent->Fill(4);
      //struct Condition Jpsi_Reco; 
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      //       if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta(), mupl_Reco->P()) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta(), mumi_Reco->P())) )

      if (  isTrigMatched(hltBits,Reco_QQ_trig[irqq]) == false )
        continue; // trigger selection
      hEvent->GetXaxis()->SetBinLabel(7,"Di-muons trig");      hEvent->Fill(7);

      DIMUTRIGPASS++;


      if ( !  ( acceptance( mupl_Reco->Pt(), mupl_Reco->Eta() ) && acceptance( mumi_Reco->Pt(), mumi_Reco->Eta()) )   )
        continue;      
      hEvent->GetXaxis()->SetBinLabel(5,"Di-muons Accep");      hEvent->Fill(5);

      // muon id cut :  tight cut  pixel cut is missing!!
      // Reference : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Soft_Muon
      bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mupl_nPixWMea[irqq] > 0) &&
          (Reco_QQ_mupl_dxy[irqq]<0.3) &&
          (Reco_QQ_mupl_dz[irqq]<20.)   			//			 &&  (Reco_QQ_mupl_highPurity[irqq]==true) 
          ) ; 

      bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mumi_nPixWMea[irqq] > 0) &&
          (Reco_QQ_mumi_dxy[irqq]<0.3) &&
          (Reco_QQ_mumi_dz[irqq]<20.)  // && (Reco_QQ_mumi_highPurity[irqq]==true)
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

      dm.softFlag = 0;
      if (muplSoft && mumiSoft) 
        dm.softFlag = 1;

      dm.highPtFlag = 0;
      if ( muplHighPtCut && mumiHighPtCut ) 
        dm.highPtFlag = 1;

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

      dm.weight0 = dmGen.weight0;
      dm.weight = dmGen.weight;

      if      ( epSelection == kEPl2HF )            dm.ep2 = rpAng[8];  
      else if ( epSelection == kEPOppositeHF )      
      {
        if ( dm.y > 0 )     dm.ep2 = rpAng[6];  // [6] = HFm2;
        else                 dm.ep2 = rpAng[7];  // [7] = HFm2;
      }
      else if ( epSelection == kEPSameSideHF )      
      {
        if ( dm.y > 0 )     dm.ep2 = rpAng[7];  
        else                 dm.ep2 = rpAng[6]; 
      }

      dm.dphiEp2 = getDPHI( dm.phi, dm.ep2);
      
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
      dataSet->add( *argSet);
    } // end of dimuon loop

  } //end of event loop
  
  
  cout << endl;
  cout << endl;
  cout << "******************" << endl;
  cout << "**** SUMMARY *****" << endl;
  cout << endl;

  cout << "1. # of events passing HLT : " << HLTPASS << endl;
  cout << "2. # of dimuons passing trigger matching after pass 1 : " << DIMUTRIGPASS << endl;
  cout << "3. # of dimuons passing muon ID cut after pass 1-2 & acceptance cut : " << DIMUIDPASS << endl;
  cout << "Cent All : # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for all central events (0-100%) " << DIMUPASS_all << endl;
  cout << "Cent 1. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (0-5%)   : " << DIMUPASS_CENT1 << endl;
  cout << "Cent 2. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (5-10%)  : " << DIMUPASS_CENT2 << endl;
  cout << "Cent 3. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (10-20%) : " << DIMUPASS_CENT3 << endl;
  cout << "Cent 4. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (20-30%) : " << DIMUPASS_CENT4 << endl;
  cout << "Cent 5. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (30-40%) : " << DIMUPASS_CENT5 << endl;
  cout << "Cent 6. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (40-50%) : " << DIMUPASS_CENT6 << endl;
  cout << "Cent 7. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (50-60%) : " << DIMUPASS_CENT7 << endl;
  cout << "Cent 8. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (60-70%) : " << DIMUPASS_CENT8 << endl;
  cout << "Cent 9. # of dimuons passing mass range 7.5-14 after pass 1-3 & vertex probability cut for centrality (70-80%) : " << DIMUPASS_CENT9 << endl;
  cout << "******************" << endl;
  cout << "******************" << endl;
  cout << endl;
  cout << endl;

  dataSet->Write();
  mmTree->Write();  // Don't need to call Write() for trees
  hEvent->Write();
  if (isMC) mmGenTree->Write();
  //   if (saveTracks) trkOutTree->Write();
  newfile->Close();
} 



/*
      Jpsi_Reco.theMass = JP_Reco->M();
      Jpsi_Reco.theRapidity = JP_Reco->Rapidity(); // y_{lab}	
      Jpsi_Reco.thePt = JP_Reco->Pt();
      Jpsi_Reco.thePhi = JP_Reco->Phi();
      
      Jpsi_Reco.mupl_p = sqrt( (mupl_Reco->Px())*(mupl_Reco->Px()) + (mupl_Reco->Py())*(mupl_Reco->Py()) + (mupl_Reco->Pz())*(mupl_Reco->Pz()) );
      Jpsi_Reco.mumi_p = sqrt( (mumi_Reco->Px())*(mumi_Reco->Px()) + (mumi_Reco->Py())*(mumi_Reco->Py()) + (mumi_Reco->Pz())*(mumi_Reco->Pz()) );
      Jpsi_Reco.mupl_pt = mupl_Reco->Pt();
      Jpsi_Reco.mumi_pt = mumi_Reco->Pt();
      Jpsi_Reco.mupl_eta = mupl_Reco->Eta();
      Jpsi_Reco.mumi_eta = mumi_Reco->Eta();
      */    

/*  ////// Gen_QQ_size loop
    for (Int_t igqq=0; igqq<Gen_QQ_size; ++igqq) {
    //struct Condition Jpsi_Gen; 
    mupl_Gen = (TLorentzVector*) Gen_QQ_mupl_4mom->At(igqq);
    mumi_Gen = (TLorentzVector*) Gen_QQ_mumi_4mom->At(igqq);
    JP_Gen_tmp_qq = (TLorentzVector*) Gen_QQ_4mom->At(igqq); // Gen Jpsi (for filling NoCut)
    *JP_Gen_tmp = *mupl_Gen +  *mumi_Gen; // Gen dimuon pairs (for filling NoCut)
	  
	  // variables only used for Reco
	  Jpsi_Gen.HLTriggers = 0;
	  Jpsi_Gen.Reco_QQ_trig  = 0;
	  Jpsi_Gen.theSign = 0; //already +- pair
	  
	  Jpsi_Gen.theCentrality = 97.5; // for pp!
	  Jpsi_Gen.theType = Gen_QQ_type[igqq]; // PR or NP
	  Jpsi_Gen.mupl_p = sqrt( (mupl_Gen->Px())*(mupl_Gen->Px()) + (mupl_Gen->Py())*(mupl_Gen->Py()) + (mupl_Gen->Pz())*(mupl_Gen->Pz()) );
	  Jpsi_Gen.mumi_p = sqrt( (mumi_Gen->Px())*(mumi_Gen->Px()) + (mumi_Gen->Py())*(mumi_Gen->Py()) + (mumi_Gen->Pz())*(mumi_Gen->Pz()) );
	  Jpsi_Gen.mupl_pt = mupl_Gen->Pt();
	  Jpsi_Gen.mumi_pt = mumi_Gen->Pt();
	  Jpsi_Gen.mupl_eta = mupl_Gen->Eta();
	  Jpsi_Gen.mumi_eta = mumi_Gen->Eta();
	  
	  // --- cut01 : GEN - No cut (only DimuCut = +-pair from 443)
	  h2D_NoCut_Gen_pt_y->Fill(JP_Gen_tmp->Rapidity(),JP_Gen_tmp->Pt()); // Gen dimuon
	  h2D_NoCut_GenJpsi_pt_y->Fill(JP_Gen_tmp_qq->Rapidity(),JP_Gen_tmp_qq->Pt()); // Gen Jpsi
	  
	  // --- cut02 : GEN for denominator
	  bool yn_gen = false;
	  if ( kineCut(Jpsi_Gen.mupl_pt, Jpsi_Gen.mupl_eta, Jpsi_Gen.mupl_p) 
			&& kineCut(Jpsi_Gen.mumi_pt, Jpsi_Gen.mumi_eta, Jpsi_Gen.mumi_p)) {
			*JP_Gen = *mupl_Gen +  *mumi_Gen; // used for actual denominator ( GEN dimuon pairs)
			Jpsi_Gen.theMass = JP_Gen->M();
			Jpsi_Gen.theRapidity = JP_Gen->Rapidity();	
			Jpsi_Gen.thePt = JP_Gen->Pt();
			Jpsi_Gen.thePhi = JP_Gen->Phi();
			if ( massCut1(Jpsi_Gen.theMass)  
			&& minpt<=Jpsi_Gen.thePt && Jpsi_Gen.thePt <maxpt 
			&& minylab<=Jpsi_Gen.theRapidity && Jpsi_Gen.theRapidity <maxylab) {
			yn_gen = true;
			}
			}
			if (yn_gen == true){
			h2D_Den_pt_y->Fill(Jpsi_Gen.theRapidity,Jpsi_Gen.thePt,zWeight01*zWeight02);
			} // end of yn_gen
			} //end of Gen_QQ_size loop
      */



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


bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
  if (  (pt1-pt2)/pt1 > 0.02 ) 
    return false;
  if ( fabs( getDPHI(phi1, phi2) ) > 0.1 )
    return false;
  if ( fabs( eta1 - eta2 ) > 0.1 )
    return false;

  return true;
}
