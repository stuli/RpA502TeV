#include <ctime>
#include "TTree.h"

#include "../../cutsAndBin.h"
#include "RooRealVar.h"

#include "RooDataSet.h"
#include "RooGaussian.h"
#include <TLorentzVector.h>
#include "../../TriggerManipulation.h" 
#include "../../commonUtility.h"
static const long MAXTREESIZE = 10000000000;
// Make tree from onia for acceptance study
// Editor : Dongho Moon (dmoon@cern.ch)

TString getDayAndTime();
void makeOnia2AccTree(int iCat_ = 2) { // iCat_ = 0 (1S), 1 (2S), 2 (3S)

  using namespace std;

  TChain *mytree = new TChain("hionia/myTree");

  TString fname;

  // Input file pp MC Upsilon 1S, 2S, 3S
  if(iCat_ == 0) fname = "/afs/cern.ch/work/g/goni/Analysis/Upsilon/upslion_RpA/src/Usercode/UpsilonpPb5TeV/OniaTrees/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  if(iCat_ == 1) fname = "/afs/cern.ch/work/g/goni/Analysis/Upsilon/upslion_RpA/src/Usercode/UpsilonpPb5TeV/OniaTrees/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";
  if(iCat_ == 2) fname = "/afs/cern.ch/work/g/goni/Analysis/Upsilon/upslion_RpA/src/Usercode/UpsilonpPb5TeV/OniaTrees/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root";

  mytree->Add(fname.Data());

  cout << endl << "*==*==*==*==*==*==*==* INPUT FILE *==*==*==*==*==*==*==*==*" << endl;
  //  if (f1->IsZombie()) { cout << "*** KYO : CANNOT open the root file!! Macro terminated ***" << endl; return;} 
  cout <<"* file ::" << fname << endl;

  TString fout;
  cout << endl << "*==*==*==*==*==*==*==* OUTPUT FILE *==*==*==*==*==*==*==*==*" << endl;
  if(iCat_ == 0) fout = "skimedForAcc_MC_Ups1S_20170808.root";
  if(iCat_ == 1) fout = "skimedForAcc_MC_Ups2S_20170808.root";
  if(iCat_ == 2) fout = "skimedForAcc_MC_Ups3S_20170808.root";
  cout <<"* file ::" << fout << endl;
  TFile *out = new TFile(fout,"RECREATE");


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
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  mytree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);


  ////////////////////////////////////////////////////////////////////////
  //////////////////  Gen QQ tree        
  ////////////////////////////////////////////////////////////////////////
  TTree *mmGenTree;
  mmGenTree = new TTree("mmGen","Gen Di-muon Pairs");
  mmGenTree->SetMaxTreeSize(MAXTREESIZE);

  int genRun = 0, genEvt = 0;
  float mass = 0.0, pt = 0.0, phi = 0.0, y = 0.0, eta = 0.0;
  float pt1 = 0.0, phi1 = 0.0, eta1 = 0.0; 
  float pt2 = 0.0, phi2 = 0.0, eta2 = 0.0; 
  
  mmGenTree->Branch("mass",&mass,"mass/F");
  mmGenTree->Branch("pt",&pt,"pt/F");
  mmGenTree->Branch("phi",&phi,"phi/F");
  mmGenTree->Branch("y",&y,"y/F");
  mmGenTree->Branch("eta",&eta,"eta/F");
  mmGenTree->Branch("pt1",&pt1,"pt1/F");
  mmGenTree->Branch("eta1",&eta1,"eta1/F");
  mmGenTree->Branch("phi1",&phi1,"phi1/F");
  mmGenTree->Branch("pt2",&pt2,"pt2/F");
  mmGenTree->Branch("eta2",&eta2,"eta2/F");
  mmGenTree->Branch("phi2",&phi2,"phi2/F");

  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Gen = new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;
  int nevt = -1;

  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  cout<<"nevt : "<<nevt<<endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    ///////////////////////////
    ///// Call the values /////
    ///////////////////////////
    mytree->GetEntry(iev);


    /////////////////////////////////////////////////////
    /////////  gen QQ loop first
    /////////////////////////////////////////////////////
    for (Int_t irqq=0; irqq<Gen_QQ_size; ++irqq) 
    {
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(irqq);
      mupl_Gen = (TLorentzVector*) Gen_QQ_mupl_4mom->At(irqq);
      mumi_Gen = (TLorentzVector*) Gen_QQ_mumi_4mom->At(irqq);

      mass = JP_Gen->M();
      pt = JP_Gen->Pt();
      phi = JP_Gen->Phi();
      y = JP_Gen->Rapidity();
      eta = JP_Gen->Eta();

      pt1 = mupl_Gen->Pt();
      eta1 = mupl_Gen->Eta();
      phi1 = mupl_Gen->Phi();
      pt2 = mumi_Gen->Pt();
      eta2 = mumi_Gen->Eta();
      phi2 = mumi_Gen->Phi();
      mmGenTree->Fill();
     
    } //end of Gen QQ tree


  } //end of event loop


  out->Write();
  out->Close();
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

