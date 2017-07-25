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
bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) ;
void onia2mmNtuple( int nevt = -1,
		    int fileID = kPPMCUps1S,
		    int trigId=kL1DoubleMu0,
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
  else if(fileID == kAAMCUps1S) inf_func = new TFile("compareDataMc/ratioDataMC_AA_DATA_1sState.root","read");
  else if(fileID == kAAMCUps2S) inf_func = new TFile("compareDataMc/ratioDataMC_AA_DATA_2sState.root","read");

  TF1* wFunc[nYBins+1];
  wFunc[1] = (TF1*) inf_func -> Get("dataMcRatio");
  //wFunc[1]  = new TF1("weightCurve_1s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2]))))",0,30); 
  //wFunc[2]  = new TF1("weightCurve_2s","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2]))))",0,30); 
/*  if ( (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) )  { 
    wFunc[1]->SetParameters( 0.988141, 3.0971, 1.81891, 10.0239);
    wFunc[2]->SetParameters(11.518, 7.53196, 2.38444, 2.68481);
  }
  else if ( (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )  { 
    wFunc[1]->SetParameters( 1.0001, 5.1, 2.0024, 12.4243);
    wFunc[2]->SetParameters( 3.46994, 11.8612, 2.10006, 3.25859);
  }
  */
  
  if  (fileID == kPPMCUps1S){
    fname = "OniaTree/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root";
    mytree->Add(fname.Data());
  }
  else { 
    cout << "options other than kPPMC Y(1S) is not ready" << endl;
    return; 
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
    newfile = new TFile(Form("skimmedFiles/mmSkimPP_L1DoubleMu0PD_Trig-%s_%s_%s_%s.root",trigName.Data(), fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  if (fileID == kPADATA) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPA_Trig-%s_%s_%s_%s.root",trigName.Data(), fdimusign.Data(), getDayAndTime().Data(), skimVersion.Data() ),"recreate");   
  }
  else if (fileID == kAADATA) { 
    newfile = new TFile(Form("skimmedFiles/mmSkimPbPb_L1DoubleMu0PD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if  (fileID == kAADATACentL3) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPbPb_CentralPD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if  (fileID == kAADATAPeri) { 
    newfile = new TFile(Form("skimmedFiles/mmSkimPbPb_PeripheralPD_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMC) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPP_MC_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMC) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPbPb_MC_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate");
  }
  else if (fileID == kPPMCUps1S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPP_MC_Ups1S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMCUps2S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPP_MC_Ups2S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kPPMCUps3S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimPP_MC_Ups3S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps1S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimAA_MC_Ups1S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps2S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimAA_MC_Ups2S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
  }
  else if (fileID == kAAMCUps3S) {
    newfile = new TFile(Form("skimmedFiles/mmSkimAA_MC_Ups3S_Trig-%s_%s_EP-%s_%s_%s.root",trigName.Data(),fdimusign.Data(),epName.Data(), getDayAndTime().Data(), skimVersion.Data()),"recreate"); 
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
  
  
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  dimuon tree 
  ////////////////////////////////////////////////////////////////////////
  DiMuon dm;
  TTree *mmTree = new TTree("mm","diMuonPairs");
  mmTree->SetMaxTreeSize(MAXTREESIZE);
  mmTree->Branch("mm",&dm.run,branchString.Data());

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

      float GenMaxPt = -1;
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
	
	if (  !( (dmGen.pt1 > glbMuPtCut) && (dmGen.pt2 > glbMuPtCut) && ( fabs(dmGen.eta1) < 2.4 ) &&  ( fabs(dmGen.eta2) <2.4 ) ) )
	  continue;
	if ( (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S) )   {
	  dmGen.weight0 = (float) hWeight->GetBinContent( hWeight->FindBin(dmGen.pt) ) ;
	}
	else  {
	  dmGen.weight0 = 1;  
	}
	// MC pT weight :  
	if ( (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) ||  (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S))  {
	    dmGen.weight = dmGen.weight0 * wFunc[1]->Eval(dmGen.pt); 
	}		 
	
	dmGen.oniaIndex = irqq;
	
	if ( dmGen.pt < GenMaxPt )  
	  continue;  // select the highest pt gen muon pair only
	GenMaxPt = dmGen.pt;
      } //end of Gen QQ tree
      mmGenTree->Fill();
      


      // Now reco tree 	
      dm.clear();
      dm.run = runNb;
      dm.lumi = LS ;
      dm.event = eventNb ;
      dm.vz = zVtx;
      dm.cBin = Centrality ;
      float maxRecoPt = -1;
      for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
	{
	  JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
	  mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
	  mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
	  if ( isTrigMatched(hltBits,HLTriggers) == false ) 
	    continue; // trigger selection
	  if (  isTrigMatched(hltBits,Reco_QQ_trig[irqq]) == false )
	    continue; // trigger selection
	  if (  !( (mupl_Reco->Pt() > glbMuPtCut) && (mumi_Reco->Pt() > glbMuPtCut) && ( fabs(mupl_Reco->Eta()) < 2.4 ) &&  ( fabs(mumi_Reco->Eta()) <2.4 ) ) )
	    continue;

	  
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
      
	  if ( !(muplSoft && mumiSoft) ) 
	    continue;   
	  if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
	    continue;
	  if( !DiMuSign )
	    {
	      if ( Reco_QQ_sign[irqq] != 0 ) continue;
	    }
	  else if(DiMuSign)
	    {
	      if( Reco_QQ_sign[irqq] == 0) continue;
	    }
	  
	  dm.softFlag = 0;
	  if (muplSoft && mumiSoft) 
	    dm.softFlag = 1;
	  
	  dm.highPtFlag = 0;
	  
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
	  if ( (fileID == kPPMCUps1S) || (fileID == kPPMCUps2S) || (fileID == kPPMCUps3S) ||  (fileID == kAAMCUps1S) || (fileID == kAAMCUps2S) || (fileID == kAAMCUps3S)  )	{
	    dm.weight0 = dmGen.weight0 ; 
	    dm.weight = dmGen.weight ; 
	  }
	  else { 
	dm.weight0 = 1.0 ; 
	dm.weight = 1.0 ; 
	  }      
	  
	  dm.oniaIndex = irqq;
	  if ( dm.pt < maxRecoPt) 
	    continue;
	  maxRecoPt = dm.pt ;
	} // end of dimuon loop
      
      mmTree->Fill();
      
      
    } //end of isMC condition for Gen QQ tree
  } //end of event loop
  
  mmTree->Write();  // Don't need to call Write() for trees
  if (isMC) mmGenTree->Write();
  //   if (saveTracks) trkOutTree->Write();
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


bool isTrackMatched(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
  if (  (pt1-pt2)/pt1 > 0.02 ) 
    return false;
  if ( fabs( getDPHI(phi1, phi2) ) > 0.1 )
    return false;
  if ( fabs( eta1 - eta2 ) > 0.1 )
    return false;

  return true;
}
