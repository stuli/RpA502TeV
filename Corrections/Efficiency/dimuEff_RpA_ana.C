//////////////////////////////////////////////////////////////////////////////////////
// This code can be used to obtain the type I error rate (or Efficiency) of the CMS //
// detector when used to perform analysis on muon pairs (dimuons). It checks the    //
// triggering, muon identification, and reconstruction efficiencies using Monte     //
// Carlo (MC) simulations of muon pairs in the Upsilon mass range. The MC efficien- //
// cies are further corrected by data/MC ratios (scale factor) of the error rate of //
// single muons in the CMS detector. To obtain the scale factors, a data-driven     //
// method known as tag-and-probe (TnP) is used. This code can be used to determine  //
// the statistical uncertainty in the efficiency corrections, as well as systematic //
// uncertainties arising from the TnP method as well as reweighting the transverse  //
// momentum (pT) spectra, if such reweighting to the MC is required.                //
//                                                                                  //
// For the analysis HIN-18-005, this code is used to determine the efficiency       //
// corrections for pp and pPb cross section, RpA and differential RpA measurements. //
//////////////////////////////////////////////////////////////////////////////////////

#include "effCommon.h"
#include "tnp_weight_pp.h"

const double muonPtCut = 4.0; // For Upsilon candidates

// Set true for efficiency in differential RpA bins
bool isRpA2D = false;

// Set true if systematics due to pT reweighting is wanted
bool isPt_reweight_sys = false;

// Select whether to use the nominal TnP weights or to determine systematic uncertainties.
// For pp, binned is also a type of TnP systematic. There are 9 total: 
// (TnP up, TnP down) X (Trigger, Tracking, MuID, STA) & TnP binned X Trigger. 
// For pPb, the only TnP single muon efficiency is a combined muon id and triggering 
// efficiency, so only upward and downward variations are tested.
enum tnp_correction {Nominal, Systematics, Binned};
enum pp_tnp_types {Trigger, Tracking, Muid, Sta};

// Open and read the rootfiles containing the TnP scale factors for pPb. SFs for pp
// are contained in tnp_weight_pp.h
// pPb TnP nominal:
TFile* fTnp_pa = new TFile("output_official_5eta_cutG_all_nominal_v3.root","READ");
TF1* hTnp_pa_eta0_09 = (TF1*)fTnp_pa->Get("func_1");
TF1* hTnp_pa_eta09_12 = (TF1*)fTnp_pa->Get("func_2");
TF1* hTnp_pa_eta12_16 = (TF1*)fTnp_pa->Get("func_3");
TF1* hTnp_pa_eta16_21 = (TF1*)fTnp_pa->Get("func_4");
TF1* hTnp_pa_eta21_24 = (TF1*)fTnp_pa->Get("func_5");

// pPb TnP systematic variations:
TFile* fTnp_pa_sys_up = new TFile("pPb_official_error_ratio_max.root", "READ");
TF1* hTnp_sys_up_pa_eta0_09 =     (TF1*)   fTnp_pa_sys_up->Get("func_1");
TF1* hTnp_sys_up_pa_eta09_12 =    (TF1*)   fTnp_pa_sys_up->Get("func_2");
TF1* hTnp_sys_up_pa_eta12_16 =    (TF1*)   fTnp_pa_sys_up->Get("func_3");
TF1* hTnp_sys_up_pa_eta16_21 =    (TF1*)   fTnp_pa_sys_up->Get("func_4");
TF1* hTnp_sys_up_pa_eta21_24 =    (TF1*)   fTnp_pa_sys_up->Get("func_5");

TFile* fTnp_pa_sys_down = new TFile("pPb_official_error_ratio_min.root", "READ");
TF1* hTnp_sys_down_pa_eta0_09 =     (TF1*)   fTnp_pa_sys_down->Get("func_1");
TF1* hTnp_sys_down_pa_eta09_12 =    (TF1*)   fTnp_pa_sys_down->Get("func_2");
TF1* hTnp_sys_down_pa_eta12_16 =    (TF1*)   fTnp_pa_sys_down->Get("func_3");
TF1* hTnp_sys_down_pa_eta16_21 =    (TF1*)   fTnp_pa_sys_down->Get("func_4");
TF1* hTnp_sys_down_pa_eta21_24 =    (TF1*)   fTnp_pa_sys_down->Get("func_5");


// Function Initializations
bool PtCut(TLorentzVector* Muon);
bool IsAccept(TLorentzVector* Muon);
bool MassCut(TLorentzVector* DiMuon, double LowM, double HighM);
double PtReweight(TLorentzVector* DiMuon, TF1 *Pt_ReWeights);
// TnP weight functions
double weight_tp_pp(double pt, double eta);
double weight_tp_pp_binned(double pt, double eta);
double sys_SF_tp_pp_trigger(double pt, double eta, int idx_variation);
double sys_SF_tp_pp_tracking(double pt, double eta, int idx_variation);
double sys_SF_tp_pp_muid(double pt, double eta, int idx_variation);
double sys_SF_tp_pp_sta(double pt, double eta, int idx_variation);
double weight_tp_pPb(double mupt1,double mupt2,double mueta1, double mueta2);
double sys_SF_tp_pPb(double mupt1,double mupt2,double mueta1, double mueta2, bool isSysUp);


// Variable declarations
int  nPtBin;     
int  nRapBin;     
int  nIntBin;

double m1S_low = 8.0;
double m1S_high = 10.0;
double m2S_low = 8.563;
double m2S_high = 10.563;
double m3S_low = 8.895;
double m3S_high = 10.895;

int iPeriod = 5;
int iPos = 33;


/////////////// Start of main function ///////////////////////////////////////////////////////
void dimuEff_RpA_ana(
int oniaMode = VVV, //1 = 1S, 2 = 2S, 3 = 3S
bool ispPb = WWW //true = pPb and false = pp
){   
    // Set style convention used at CMS
	setTDRStyle();

	cout << Form("Calculating corrected efficiencies for %dS in %s collisions.", oniaMode, ispPb ? "pPb" : "pp") << endl;

    // Tag and Probe settings:
	// Tag and Probe correction type
	int tnp_type = tnp_choice;
    // Upward or downward fluctuation
    char fluctuation_dirrection[] = "up_or_down";
	// For tnp up and down systematics in pp collisions, also choose the type of tnp scale factors to be used 
	int pp_tnp_type;
    if(!ispPb) pp_tnp_type = pp_tnp_sys_choice;

	// Loading data trees for reference
	if(ispPb){
        TChain myTree_Data("myTree");
		myTree_Data.Add("../../Data/pPb_Data/RD2013_pa_1st_run_merged.root");
		cout<<"Number of records in Pb->P data tree = "<<myTree_Data.GetEntries()<<endl;
		myTree_Data.Add("../../Data/pPb_Data/RD2013_pa_2nd_run_merged.root");
		cout<<"Number of records in p->Pb data tree = "<<myTree_Data.GetEntries()<<endl;
	}
	else{
	        TChain myTree_Data("hionia/myTree");
		myTree_Data.Add("../../Data/pp_Data/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root");
    		cout<<"Number of records in pp data tree = "<<myTree_Data.GetEntries()<<endl;
	}

	// Loading pp collision MC trees for efficiency determination
	TChain *myTree_pp = new TChain("hionia/myTree");
        
    if ((oniaMode == 1) && !ispPb){
		myTree_pp->Add("../../Data/pp_MC/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
	}

    if ((oniaMode == 2) && !ispPb){
        myTree_pp->Add("../../Data/pp_MC/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
    }

    if ((oniaMode == 3) && !ispPb){
    	myTree_pp->Add("../../Data/pp_MC/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
    }

    // Loading pPb collision MC trees for efficiency determination
    TChain *myTree_pPb = new TChain("myTree");

	if ((oniaMode == 1) && ispPb){
	   myTree_pPb->Add("../../Data/pPb_MC/OniaTree_Ups1S_PA_MC_PbP_5TeV02.root");
	   // 1S p->Pb MC included below. Would need to change rapidity bins etc. to use.
	   // myTree_pPb->Add("../../Data/pPb_MC/OniaTree_Ups1S_PA_MC_pPb_5TeV02.root");
	}

    if ((oniaMode == 2) && ispPb){
    	myTree_pPb->Add("../../Data/pPb_MC/OniaTree_Ups2S_PA_MC_PbP_5TeV02.root");
    }

    if ((oniaMode == 3) && ispPb){
        myTree_pPb->Add("../../Data/pPb_MC/OniaTree_Ups3S_PA_MC_PbP_5TeV02.root");
    }

    TChain *myTree;
    if(!ispPb){
    	myTree = (TChain*)myTree_pp;
    }
	else{
        myTree = (TChain*)myTree_pPb;
	}
	cout<<"Number of records in MC tree = "<<myTree->GetEntries()<<endl;


	// Declaring variables to be read from the data	
	Float_t        muMiDxy;
	Float_t        muMiDz;
	Int_t          muMiNPxlLayers;
	Int_t          muMiNTrkLayers;
	Float_t        muPlDxy;
	Float_t        muPlDz;
	Int_t          muPlNPxlLayers;
	Int_t          muPlNTrkLayers;
	Float_t        vProb;

    Bool_t         muPlGoodMu;
	Bool_t         muMiGoodMu;

	Bool_t		   muPlHighPurity;
	Bool_t		   muPlTrkMuArb;
	Bool_t		   muPlTMOneStaTight;
	Bool_t         muMiHighPurity;
	Bool_t         muMiTrkMuArb;
	Bool_t         muMiTMOneStaTight;

	Int_t           Centrality;
    Int_t           HLTriggers;
	ULong64_t		HLTriggers_pp;
	Int_t           Reco_QQ_size;
	Int_t           Reco_QQ_sign[45];    
	TClonesArray    *Reco_QQ_4mom;
	TClonesArray    *Reco_QQ_mupl_4mom;
	TClonesArray    *Reco_QQ_mumi_4mom;
	Int_t           Reco_QQ_trig[45];
	ULong64_t       Reco_QQ_trig_pp[45]; 
	Float_t         Reco_QQ_VtxProb[45];    
	Int_t           Reco_QQ_mupl_nPixWMea[45];    
	Int_t           Reco_QQ_mumi_nPixWMea[45];    
	Int_t           Reco_QQ_mupl_nTrkWMea[45];    
	Int_t           Reco_QQ_mumi_nTrkWMea[45];    
	Float_t         Reco_QQ_mupl_dxy[45];    
	Float_t         Reco_QQ_mumi_dxy[45];    
	Float_t         Reco_QQ_mupl_dz[45];    
	Float_t         Reco_QQ_mumi_dz[45];    

	Bool_t          Reco_QQ_mupl_isGoodMuon[45];
	Bool_t          Reco_QQ_mumi_isGoodMuon[45];

    Bool_t        	Reco_QQ_mupl_isHighPurity[45];
    Bool_t        	Reco_QQ_mupl_TrkMuArb[45];
    Bool_t        	Reco_QQ_mupl_TMOneStaTight[45];
    Bool_t        	Reco_QQ_mumi_isHighPurity[45];
    Bool_t        	Reco_QQ_mumi_TrkMuArb[45];
    Bool_t        	Reco_QQ_mumi_TMOneStaTight[45];

	Int_t           Gen_QQ_size;
	Int_t           Gen_QQ_sign[45];    
	TClonesArray    *Gen_QQ_4mom;
	TClonesArray    *Gen_QQ_mupl_4mom;
	TClonesArray    *Gen_QQ_mumi_4mom;

    // Setting pointers for branches
	TBranch        	*b_Centrality;    
	TBranch        	*b_HLTriggers;
    TBranch         *b_HLTriggers_pp;    
	TBranch        	*b_Reco_QQ_size;    
	TBranch        	*b_Reco_QQ_sign;    
	TBranch        	*b_Reco_QQ_4mom;    
	TBranch        	*b_Reco_QQ_mupl_4mom;    
	TBranch        	*b_Reco_QQ_mumi_4mom;    
	TBranch        	*b_Reco_QQ_trig; 
    TBranch         *b_Reco_QQ_trig_pp;
	TBranch        	*b_Reco_QQ_VtxProb;    
	TBranch        	*b_Reco_QQ_mupl_nPixWMea;    
	TBranch        	*b_Reco_QQ_mumi_nPixWMea;    
	TBranch        	*b_Reco_QQ_mupl_nTrkWMea;    
	TBranch        	*b_Reco_QQ_mumi_nTrkWMea;    
	TBranch        	*b_Reco_QQ_mupl_dxy;    
	TBranch        	*b_Reco_QQ_mumi_dxy;    
	TBranch        	*b_Reco_QQ_mupl_dz;    
	TBranch        	*b_Reco_QQ_mumi_dz;    

	TBranch        	*b_Reco_QQ_mupl_isGoodMuon;
	TBranch        	*b_Reco_QQ_mumi_isGoodMuon;

	TBranch			*b_Reco_QQ_mupl_isHighPurity;
	TBranch			*b_Reco_QQ_mupl_TrkMuArb;
	TBranch			*b_Reco_QQ_mupl_TMOneStaTight;
	TBranch        	*b_Reco_QQ_mumi_isHighPurity;
	TBranch        	*b_Reco_QQ_mumi_TrkMuArb;
	TBranch        	*b_Reco_QQ_mumi_TMOneStaTight;

	TBranch        	*b_Gen_QQ_size;  
	TBranch        	*b_Gen_QQ_4mom;    
	TBranch        	*b_Gen_QQ_mupl_4mom;    
	TBranch        	*b_Gen_QQ_mumi_4mom;    


	// Grabbing data from trees. 
    // Set object pointers and Initialize.
	Reco_QQ_4mom = 0;
	Reco_QQ_mupl_4mom = 0;
	Reco_QQ_mumi_4mom = 0;

	Gen_QQ_4mom = 0;
	Gen_QQ_mupl_4mom = 0;
	Gen_QQ_mumi_4mom = 0;

	myTree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	if(ispPb) myTree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
    else myTree->SetBranchAddress("HLTriggers", &HLTriggers_pp, &b_HLTriggers_pp);
	myTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	myTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	myTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	myTree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
	myTree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
	if(ispPb) myTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
    else myTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig_pp, &b_Reco_QQ_trig_pp);
	myTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
	myTree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
	myTree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
	myTree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
	myTree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
	myTree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
	myTree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
	myTree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
	myTree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
	
	if(!ispPb) {
		myTree->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
		myTree->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
	}
	else {
		myTree->SetBranchAddress("Reco_QQ_mupl_isHighPurity", Reco_QQ_mupl_isHighPurity, &b_Reco_QQ_mupl_isHighPurity);
		myTree->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
		myTree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
        myTree->SetBranchAddress("Reco_QQ_mumi_isHighPurity", Reco_QQ_mumi_isHighPurity, &b_Reco_QQ_mumi_isHighPurity);
        myTree->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
        myTree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
	}

	myTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
	myTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
	myTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
	myTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

    // Activate branches
	myTree->SetBranchStatus("*", 0);

	myTree->SetBranchStatus("Centrality", 1);
	myTree->SetBranchStatus("HLTriggers", 1);
	myTree->SetBranchStatus("Reco_QQ_size", 1);
	myTree->SetBranchStatus("Reco_QQ_sign", 1);
	myTree->SetBranchStatus("Reco_QQ_4mom", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_4mom", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_4mom", 1);
    myTree->SetBranchStatus("Reco_QQ_trig", 1);
	myTree->SetBranchStatus("Reco_QQ_VtxProb", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_nPixWMea", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_nPixWMea", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_nTrkWMea", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_nTrkWMea", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_dxy", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_dxy", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_dz", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_dz", 1);

	if(!ispPb){
		myTree->SetBranchStatus("Reco_QQ_mupl_isGoodMuon", 1);
		myTree->SetBranchStatus("Reco_QQ_mumi_isGoodMuon", 1);
	}
	else{
		myTree->SetBranchStatus("Reco_QQ_mupl_isHighPurity", 1);
		myTree->SetBranchStatus("Reco_QQ_mupl_TrkMuArb", 1);
		myTree->SetBranchStatus("Reco_QQ_mupl_TMOneStaTight", 1);
		myTree->SetBranchStatus("Reco_QQ_mumi_isHighPurity", 1);
        myTree->SetBranchStatus("Reco_QQ_mumi_TrkMuArb", 1);
        myTree->SetBranchStatus("Reco_QQ_mumi_TMOneStaTight", 1);
	}

	myTree->SetBranchStatus("Gen_QQ_size", 1);
	myTree->SetBranchStatus("Gen_QQ_4mom", 1);
	myTree->SetBranchStatus("Gen_QQ_mupl_4mom", 1);
	myTree->SetBranchStatus("Gen_QQ_mumi_4mom", 1);


    // Declare and initialize constant variables
    const int nPtBins1s  = 6;  
    const int nPtBins2s  = 3;
    const int nPtBins3s  = 2;  
    // y bins including most bkwd for pPb
    const int nYBins1S  = 9; 
    const int nYBins2S  = 5; 
    const int nYBins3S  = 3; 
    // y bins for symm region only for pp
    const int nYBins1Spp  = 4; 
    const int nYBins2Spp  = 2; 
    const int nYBins3Spp  = 1; 
    
    std::vector<double> ptBinEdges;
    std::vector<double> ptBin;
    std::vector<double> rapBinEdges;
    std::vector<double> rapBin;
    
    // Declare the number of bins and assign bin edges for histograms
    if(isRpA2D){ 
        nPtBin = 2;
        ptBinEdges = {0.0,6.0,30.0};
        ptBin = {3.0,18.0};
        if(ispPb){
            nRapBin = 2;
            rapBinEdges = {-1.93, 0, 1.93};
            rapBin = {-0.965, 0.965};
        }
        else{
            nRapBin = 1;
            rapBinEdges = {0, 1.93};
            rapBin = {0.965};
        }
    }
    else{
        if(oniaMode ==1){
            nPtBin = nPtBins1s;
            if(ispPb){nRapBin = nYBins1S;}
            else{nRapBin = nYBins1Spp;}
        
            ptBinEdges = {0,2,4,6,9,12,30};
            ptBin = {1,3,5,7.5,10.5,21};
            if(ispPb){
                rapBinEdges = {-2.87, -1.93, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.93}; 
                rapBin = {-2.4, -1.565, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.565}; 
            }
            else{
                // Note: generalizable to 0 < |y| <2.4 for comparison to symmetric collision systems
                rapBinEdges = {0, 0.4, 0.8, 1.2, 1.93}; 
                    rapBin = {0.2, 0.6, 1.0, 1.565}; 
            }
        }
        if(oniaMode ==2){
            nPtBin = nPtBins2s;
            if(ispPb){nRapBin = nYBins2S;}
            else{nRapBin = nYBins2Spp;}
        
            ptBinEdges = {0,4,9,30};
            ptBin = {2,6.5,19.5};
            if(ispPb){
                rapBinEdges = {-2.87, -1.93, -0.8, 0, 0.8, 1.93}; 
                rapBin = {-2.4, -1.3565, -0.4, 0.4, 1.365}; 
            }
            else{
                rapBinEdges = {0, 0.8, 1.93};
                rapBin = {0.4, 1.365};
            }
        }
        if(oniaMode ==3){
            nPtBin = nPtBins3s;
            if(ispPb){nRapBin = nYBins3S;}
            else{nRapBin = nYBins3Spp;}
        
            ptBinEdges = {0.0,6.0,30.0};
            ptBin = {3.0,18.0};
            if(ispPb){
                rapBinEdges = {-2.87, -1.93, 0, 1.93}; 
                rapBin = {-2.4, -0.965, 0.965}; 
            }   
            else{
                rapBinEdges = {0, 1.93};
                rapBin = {0.965};   
            }
        }
    }

    // The pPb rapidity cuts are for Run 1 Only. We only have MC for run 1 for all 3 states.
    float rapLow = 0.0; // For pp XS
    float rapHigh = 1.93; // For pp XS //Possible to go up to 2.4
    float rapLowRpA = 0.0; // For pp RpA 
    float rapHighRpA = 1.93; // Same for pp and pPb, for RpA
    if(ispPb){
        rapLow = -2.87; // For pPb XS asymmetric
        rapHigh = 1.93; // For pPb XS
        rapLowRpA = -rapHighRpA; // For pPb RpA and XS symmetric
    }
    // For 2D RpA
    float rapLowRpANeg;
    float rapLowRpAPos;
    float rapHighRpANeg;
    float rapHighRpAPos;
    float lowpTLow;
    float highpTLow;
    float lowpTHigh;
    float highpTHigh;

    if(isRpA2D){
        rapLowRpANeg = 0.0;
        rapLowRpAPos = 0.0;
        rapHighRpANeg = 1.93;
        rapHighRpAPos = 1.93;
        if (ispPb) {
            rapLowRpANeg = -1.93;
            rapHighRpANeg = 0.0;
            rapLowRpAPos = 0.0;
            rapHighRpAPos = 1.93;
        }   
        lowpTLow = 0.0;
        lowpTHigh = 6.0;
        highpTLow = 6.0;
        highpTHigh = 30.0;
    }
  
    float       IntBin[1] = { 50 };
    float       IntBinEdges[2] = { 0, 100 };

    // Needed for TnP weighting and pT reweighting
    double      weighttp;   
    int         idx_nom = 0;
    int         idx_variation;

    float       ptReweight = 0.0;
    float       ptReweight_XS = 0.0;

    // Setting mass ranges for Upsilon states
    float       massLow = 0;
    float       massHigh = 0;

        if (oniaMode == 1){
            massLow = m1S_low;
            massHigh = m1S_high;
        }
        else if (oniaMode == 3){
            massLow = m3S_low;
            massHigh = m3S_high;
        }
        else{
            massLow = m2S_low;
            massHigh = m2S_high;
        }


	//convert bin and bin edge vectors to arrays to be used as parameter for declaring TH1Ds
	double* ptBinEdges_arr = &ptBinEdges[0];
	double* ptBin_arr = &ptBin[0];
	double* rapBinEdges_arr = &rapBinEdges[0];
	double* rapBin_arr = &rapBin[0];

	// Declaring histograms
	TH1D  *RecoEventsInt = new TH1D("RecoEventsInt", "Reconstructed", 1, IntBinEdges);
	TH1D  *GenEventsInt = new TH1D("GenEventsInt", "Generated", 1, IntBinEdges);
    TH1D  *RecoEventsIntRpA = new TH1D("RecoEventsIntRpA", "Reconstructed", 1, IntBinEdges);
    TH1D  *GenEventsIntRpA = new TH1D("GenEventsIntRpA", "Generated", 1, IntBinEdges);
	
	TH1D  *RecoEventsPt = new TH1D("RecoEventsPt", "Reconstructed", nPtBin, ptBinEdges_arr);
	TH1D  *GenEventsPt = new TH1D("GenEventsPt", "Generated", nPtBin, ptBinEdges_arr);
    TH1D  *RecoEventsPtRpA = new TH1D("RecoEventsPtRpA", "Reconstructed", nPtBin, ptBinEdges_arr);
    TH1D  *GenEventsPtRpA = new TH1D("GenEventsPtRpA", "Generated", nPtBin, ptBinEdges_arr);

	TH1D  *RecoEventsRap = new TH1D("RecoEventsRap", "Reconstructed", nRapBin, rapBinEdges_arr);
	TH1D  *GenEventsRap = new TH1D("GenEventsRap", "Generated", nRapBin, rapBinEdges_arr);
	TH1D  *RecoEventsRapRpA = new TH1D("RecoEventsRapRpA", "Reconstructed", nRapBin, rapBinEdges_arr);
	TH1D  *GenEventsRapRpA = new TH1D("GenEventsRapRpA", "Generated", nRapBin, rapBinEdges_arr);


    TH1D  *RecoEventsRapRpAlowpT = new TH1D("RecoEventsRapRpAlowpT", "Reconstructed", nRapBin, rapBinEdges_arr);
    TH1D  *GenEventsRapRpAlowpT = new TH1D("GenEventsRapRpAlowpT", "Generated", nRapBin, rapBinEdges_arr);

    TH1D  *RecoEventsRapRpAhighpT = new TH1D("RecoEventsRapRpAhighpT", "Reconstructed", nRapBin, rapBinEdges_arr);
    TH1D  *GenEventsRapRpAhighpT = new TH1D("GenEventsRapRpAhighpT", "Generated", nRapBin, rapBinEdges_arr);

    TH1D  *RecoEventsPtRpArapNeg = new TH1D("RecoEventsPtRpArapNeg", "Reconstructed", nPtBin, ptBinEdges_arr);
    TH1D  *GenEventsPtRpArapNeg = new TH1D("GenEventsPtRpArapNeg", "Generated", nPtBin, ptBinEdges_arr);

    TH1D  *RecoEventsPtRpArapPos = new TH1D("RecoEventsPtRpArapPos", "Reconstructed", nPtBin, ptBinEdges_arr);
    TH1D  *GenEventsPtRpArapPos = new TH1D("GenEventsPtRpArapPos", "Generated", nPtBin, ptBinEdges_arr);

    TH1D  *RecoEventsIntRpArapNeg = new TH1D("RecoEventsIntRpArapNeg", "Reconstructed", 1, IntBinEdges);
    TH1D  *GenEventsIntRpArapNeg = new TH1D("GenEventsIntRpArapNeg", "Generated", 1, IntBinEdges);

    TH1D  *RecoEventsIntRpArapPos = new TH1D("RecoEventsIntRpArapPos", "Reconstructed", 1, IntBinEdges);
    TH1D  *GenEventsIntRpArapPos = new TH1D("GenEventsIntRpArapPos", "Generated", 1, IntBinEdges);

	TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);

    // Store sum of square of weights before filling
	RecoEventsRapRpAlowpT->Sumw2();
	GenEventsRapRpAlowpT->Sumw2();
	RecoEventsRapRpAhighpT->Sumw2();
	GenEventsRapRpAhighpT->Sumw2();
	RecoEventsPtRpArapNeg->Sumw2();
	GenEventsPtRpArapNeg->Sumw2();
    RecoEventsPtRpArapPos->Sumw2();
    GenEventsPtRpArapPos->Sumw2();

    RecoEventsIntRpArapNeg->Sumw2();
    GenEventsIntRpArapNeg->Sumw2();
    RecoEventsIntRpArapPos->Sumw2();
    GenEventsIntRpArapPos->Sumw2();

	RecoEventsInt->Sumw2();
	GenEventsInt->Sumw2();
	RecoEventsPt->Sumw2();
	GenEventsPt->Sumw2();
    RecoEventsIntRpA->Sumw2();
    GenEventsIntRpA->Sumw2();
    RecoEventsPtRpA->Sumw2();
    GenEventsPtRpA->Sumw2();
	RecoEventsRap->Sumw2();
	GenEventsRap->Sumw2();
	hCrossCheck->Sumw2();


	// Get pT Reweight functions from rootfiles. 
    // Different reweighting functions for pp and pPb, as well as pPb symm and assym regions.
	const char *f_name;
	const char *f_name_XS;
	if(!ispPb){
		if(oniaMode == 1){
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_1s_Rpa_20181217.root";
		}else if(oniaMode ==2){
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_2s_Rpa_20181217.root";
		}else{
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_3s_Rpa_20181217.root";
		}
	}
	else{
		if(oniaMode == 1){
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_1s_Rpa_20181217.root";
			f_name_XS = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_1s_Cross_20181217.root";
		}else if(oniaMode ==2){
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_2s_Rpa_20181217.root";
			f_name_XS = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_2s_Cross_20181217.root";
		}else{
			f_name = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_3s_Rpa_20181217.root";
			f_name_XS = "../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_3s_Cross_20181217.root";
		}
	}

	TFile* PtReweightFunctions = new TFile(f_name, "Open");
    TF1* Pt_ReWeights = (TF1*)PtReweightFunctions->Get("dataMC_Ratio_norm");
	TFile* PtReweightFunctions_XS;
	TF1* Pt_ReWeights_XS;
	if(ispPb){
		PtReweightFunctions_XS = new TFile(f_name_XS, "Open");
		Pt_ReWeights_XS = (TF1*)PtReweightFunctions_XS->Get("dataMC_Ratio_norm");
	}


    // Get total number of event records in tree
	Long64_t nentries = myTree->GetEntries();

	// Start of calculation
	for (Long64_t jentry = 0; jentry < nentries; jentry++){

		myTree->GetEntry(jentry);
		if(jentry%100000 == 0){
			cout<<"--Processing Event: "<<jentry<<endl;
		}

		if (Gen_QQ_size > 0) {   // Analyzing events that have at least one generated dimuon for MC with some issue. 

			//Numerator Loop - RECO (Muon pairs reconstructed in the simulated detector)
			for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++){
				hCrossCheck->Fill(1);
				TLorentzVector *qq4mom = (TLorentzVector*)Reco_QQ_4mom->At(iQQ); // muon pair 4-momentum
				TLorentzVector *mumi4mom = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iQQ); // mu^{-} 4-momentum
				TLorentzVector *mupl4mom = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iQQ); // mu^{+} 4-momentum
	
				//--Muon identification (Muid) cuts for muon minus
				muMiDxy = Reco_QQ_mumi_dxy[iQQ]; // transverse distance from primary vertex
				muMiDz = Reco_QQ_mumi_dz[iQQ]; // longitudinal distance form primary vertex
				muMiNPxlLayers = Reco_QQ_mumi_nPixWMea[iQQ]; // number of pixel layers with hits
				muMiNTrkLayers = Reco_QQ_mumi_nTrkWMea[iQQ]; // number of tracker later with hits
				if(!ispPb) muMiGoodMu = Reco_QQ_mumi_isGoodMuon[iQQ]; // Conditions for being considered a good muon
				else {                                                
					muMiHighPurity = Reco_QQ_mumi_isHighPurity[iQQ];  // Different in pPb and pp collisions
					muMiTrkMuArb = Reco_QQ_mumi_TrkMuArb[iQQ];
					muMiTMOneStaTight = Reco_QQ_mumi_TMOneStaTight[iQQ];
					
					muMiGoodMu = muMiHighPurity && muMiTrkMuArb && muMiTMOneStaTight;
				}
				
				//--Muid cuts for muon plus
				muPlDxy = Reco_QQ_mupl_dxy[iQQ];
				muPlDz = Reco_QQ_mupl_dz[iQQ];
				muPlNPxlLayers = Reco_QQ_mupl_nPixWMea[iQQ];
				muPlNTrkLayers = Reco_QQ_mupl_nTrkWMea[iQQ];
				if(!ispPb) muPlGoodMu = Reco_QQ_mupl_isGoodMuon[iQQ];
                else { 
                    muPlHighPurity = Reco_QQ_mupl_isHighPurity[iQQ];
                    muPlTrkMuArb = Reco_QQ_mupl_TrkMuArb[iQQ];
                    muPlTMOneStaTight = Reco_QQ_mupl_TMOneStaTight[iQQ];

                    muPlGoodMu = muPlHighPurity && muPlTrkMuArb && muPlTMOneStaTight;
                }				

				// Vertex matching probability
				vProb = Reco_QQ_VtxProb[iQQ];
	               
                // About to apply selections (or "cuts") to data to maximize signal-to-background ratio.
                // Same selections are used in all components of the analysis.
				bool mupl_cut = 0;
				bool mumi_cut = 0;
				bool trigL1Dmu = 0;
				bool PtCutPass = 0;
				bool MassCutPass = 0;
				bool acceptMu = 0;
	
				//--Muon id selections for muon plus and muon minus.
	            if ( (muPlGoodMu == 1) && muPlNTrkLayers > 5 && muPlNPxlLayers > 0 && TMath::Abs(muPlDxy) < 0.3 && 
					TMath::Abs(muPlDz) < 20 && vProb > 0.01){ mupl_cut = 1; }
                if ( (muMiGoodMu == 1) && muMiNTrkLayers > 5 && muMiNPxlLayers > 0 && TMath::Abs(muMiDxy) < 0.3 && 
					TMath::Abs(muMiDz) < 20){ mumi_cut = 1; }
	
				// Check to see if muons satisfy pT selections and are in acceptance
				if (PtCut(mupl4mom) && PtCut(mumi4mom)){ PtCutPass = 1; }
				if (IsAccept(mupl4mom) && IsAccept(mumi4mom)){ acceptMu = 1; }
                // Check to see muon pair satisfies invariant mass requirement
				MassCutPass = MassCut(qq4mom, massLow, massHigh);
	
				// Check if trigger bit is matched to dimuon triggers
                if(ispPb){
				    if ((HLTriggers & 1) == 1 && (Reco_QQ_trig[iQQ] & 1) == 1) {trigL1Dmu = 1;}
				}else{
                    if ((HLTriggers_pp & 1) == 1 && (Reco_QQ_trig_pp[iQQ] & 1) == 1) {trigL1Dmu = 1;}
                }

				//getting reco pt and y
				float ptReco = 0;
				float rapReco = 0;
				float rapRecoCM = 0;
				ptReco = qq4mom->Pt();
	
				if(!ispPb){rapReco = TMath::Abs(qq4mom->Rapidity());
					rapRecoCM = rapReco;
				}
				else{rapReco = qq4mom->Rapidity();
					rapRecoCM = (-1.*rapReco)-0.47;  // Run 1 MC y_lab is flipped and then shifted by -0.47 to give y_CM.
				}

                // Preparing to get weights
                // TnP weights only needed for reco
                weighttp=1.0;
                float weight = 0;
                float weight_XS = 0;
                ptReweight = 0;
                ptReweight_XS = 0;
	
                // Getting pT reweights
				ptReweight = PtReweight(qq4mom, Pt_ReWeights);
				if(ispPb) ptReweight_XS = PtReweight(qq4mom, Pt_ReWeights_XS); 

				// Getting Tag and Probe single muon efficiency corrections as weights for dimuons
				if(!ispPb){
					switch(tnp_type){
						case Nominal:	
                            weighttp = weight_tp_pp(mupl4mom->Pt(),mupl4mom->Eta()) * \
			   					weight_tp_pp(mumi4mom->Pt(),mumi4mom->Eta());
			   				break;
                        case Systematics:   
                            if(strcmp(fluctuation_dirrection,"Up")==0) idx_variation = -1;
                            else if(strcmp(fluctuation_dirrection,"Down")==0) idx_variation = -2;
                            switch(pp_tnp_type){
                                case Trigger:   
                                    weighttp =  sys_SF_tp_pp_trigger(mupl4mom->Pt(), mupl4mom->Eta(), idx_variation) * \
                                         sys_SF_tp_pp_trigger(mumi4mom->Pt(), mumi4mom->Eta(), idx_variation);
                                     break;
                                 case Tracking:  
                                    weighttp = sys_SF_tp_pp_tracking(mupl4mom->Pt(), mupl4mom->Eta(), idx_variation) * \
                                         sys_SF_tp_pp_tracking(mumi4mom->Pt(), mumi4mom->Eta(), idx_variation);
                                     break;
                                 case Muid:      
                                    weighttp = sys_SF_tp_pp_muid(mupl4mom->Pt(), mupl4mom->Eta(), idx_variation) * \
                                         sys_SF_tp_pp_muid(mumi4mom->Pt(), mumi4mom->Eta(), idx_variation);
                                     break;
                                 case Sta:       
                                    weighttp =  sys_SF_tp_pp_sta(mupl4mom->Pt(), mupl4mom->Eta(), idx_variation) * \
                                         sys_SF_tp_pp_sta(mumi4mom->Pt(), mumi4mom->Eta(), idx_variation);
                                     break;
                             }
                        case Binned:	
                            weighttp = 	weight_tp_pp_binned(mupl4mom->Pt(),mupl4mom->Eta()) * \
                				weight_tp_pp_binned(mumi4mom->Pt(),mumi4mom->Eta());
                			break;
					}
				}
				else{
					switch(tnp_type){
						case Nominal:	    
                            weighttp = weight_tp_pPb(mupl4mom->Pt(),mumi4mom->Pt(),mupl4mom->Eta(), mumi4mom->Eta());
                            break;
                        case Systematics:   
                            if(strcmp(fluctuation_dirrection,"Up")==0) idx_variation = 1;
                            else if(strcmp(fluctuation_dirrection,"Down")==0) idx_variation = 0;
                            weighttp = sys_SF_tp_pPb(mupl4mom->Pt(),mumi4mom->Pt(),mupl4mom->Eta(), mumi4mom->Eta(),idx_variation);
                            break;
					}
				}
	
				// Nominal TnP weights.  
				weight = ptReweight * weighttp ; if(ispPb) weight_XS = ptReweight_XS * weighttp ;

                //For ptReweight Systematics. 
				if(isPt_reweight_sys && tnp_type==Nominal){
					weight = weighttp; weight_XS = weighttp;
				}

				// Typically not required, but if cross check is need:
                // No TnP correction (only in numerator)
				//weight = ptReweight ; if(ispPb) weight_XS = ptReweight_XS;

                // Checking to see if muon pair passes all requirements and are reconstructed
				bool recoPass = 0;
				if (Reco_QQ_sign[iQQ] == 0 && mupl_cut && mumi_cut && trigL1Dmu){ recoPass = 1; } 
	
				// Filling RecoEvent histograms if muon pair passes all above cuts
				if (recoPass == 1 && PtCutPass == 1 && MassCutPass == 1 && acceptMu == 1){
					// Filling histograms for differential RpA study
					if(isRpA2D){
						if ((rapLowRpANeg < rapRecoCM) && (rapRecoCM < rapHighRpANeg) && ptReco < 30 ){
                            RecoEventsPtRpArapNeg->Fill(ptReco, weight);
							RecoEventsIntRpArapNeg->Fill(Centrality/2., weight);
                        }
						if ((rapLowRpAPos < rapRecoCM) && (rapRecoCM < rapHighRpAPos) && ptReco < 30 ){
                            RecoEventsPtRpArapPos->Fill(ptReco, weight);
							RecoEventsIntRpArapPos->Fill(Centrality/2., weight);
                        }
						if ((rapLowRpANeg < rapRecoCM) && (rapRecoCM < rapHighRpAPos) && (lowpTLow < ptReco) \
								&& (ptReco < lowpTHigh) ){
                            RecoEventsRapRpAlowpT->Fill(rapRecoCM, weight);
                        }
						if ((rapLowRpANeg < rapRecoCM) && (rapRecoCM < rapHighRpAPos) && (highpTLow < ptReco) \
								&& (ptReco < highpTHigh) ){
                            RecoEventsRapRpAhighpT->Fill(rapRecoCM, weight);
                        }
					}
					// Filling histograms for RpA and cross section study
					else {
						// asymmetric y range
						if ((rapLow < rapRecoCM) && (rapRecoCM < rapHigh) && ptReco < 30){ 
                            if(ispPb){
                            	RecoEventsInt->Fill(Centrality/2., weight_XS);
                            	RecoEventsPt->Fill(ptReco, weight_XS);
                            	RecoEventsRap->Fill(rapRecoCM, weight_XS);
                            }
						    else{
						    	RecoEventsInt->Fill(Centrality/2., weight);
						    	RecoEventsPt->Fill(ptReco, weight);
						    	RecoEventsRap->Fill(rapRecoCM, weight);
						    }
                        }
                        // Additionally, for symmetric y range, for RpA and symm XS
                        if ((rapLowRpA < rapRecoCM) && (rapRecoCM < rapHighRpA) && ptReco < 30 ){ 
                            RecoEventsIntRpA->Fill(Centrality/2., weight);
                            RecoEventsPtRpA->Fill(ptReco, weight);
							RecoEventsRapRpA->Fill(rapRecoCM, weight);
                        }
					}	
				} // End Reco Histogram filling loop
			} // End Numerator loop
	
	
			// Denominator loop - GEN (Muon pairs generated by physics engine)
            // Here selections test that Upsilon candidate muon pairs were generated in the MC
			for (int iQQ = 0; iQQ < Gen_QQ_size; iQQ++){
	
				hCrossCheck->Fill(0);
				TLorentzVector *g_qq4mom = (TLorentzVector*)Gen_QQ_4mom->At(iQQ);
				TLorentzVector *g_mumi4mom = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iQQ);
				TLorentzVector *g_mupl4mom = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iQQ);
	
				bool acceptMu = 0;
				bool PtCutPass = 0;
				bool MassCutPass = 0;
	
				if (PtCut(g_mupl4mom) && PtCut(g_mumi4mom)){ PtCutPass = 1; }
				if (IsAccept(g_mupl4mom) && IsAccept(g_mumi4mom)){ acceptMu = 1; }
				MassCutPass = MassCut(g_qq4mom, massLow, massHigh);
	
				//getting gen pt and y 
				float ptGen = 0;
				float rapGen = 0;
				float rapGenCM = 0;
				ptGen = g_qq4mom->Pt();
	
				if(!ispPb){rapGen = TMath::Abs(g_qq4mom->Rapidity());
					rapGenCM = rapGen;
				}
				else{rapGen = g_qq4mom->Rapidity();
					rapGenCM = (-1.*rapGen)-0.47;
				}
	
                // Only pT reweights in generated muons/dimuons. No Tag and Probe.
                float weight = 0;
                float weight_XS = 0;
                ptReweight = 0;
                ptReweight_XS = 0;

				ptReweight = PtReweight(g_qq4mom, Pt_ReWeights);
				if(ispPb) ptReweight_XS = PtReweight(g_qq4mom, Pt_ReWeights_XS);
	
				// Nominal weights. 
				weight = ptReweight; if(ispPb) weight_XS = ptReweight_XS;

                // ptReweight systematics
                if(isPt_reweight_sys && tnp_type==Nominal){
				weight = 1.0; weight_XS = 1.0;
                }
	
				// Filling GenEvent histograms if generated muon pair passes all above cuts
				if (PtCutPass == 1 && MassCutPass == 1 && acceptMu == 1){
					// Filling histograms for differential RpA study
                    if(isRpA2D){
                        if ((rapLowRpANeg < rapGenCM) && (rapGenCM < rapHighRpANeg) && ptGen < 30 ){
                            GenEventsIntRpArapNeg->Fill(Centrality/2., weight);
                            GenEventsPtRpArapNeg->Fill(ptGen, weight);
                        }
                        if ((rapLowRpAPos < rapGenCM) && (rapGenCM < rapHighRpAPos) && ptGen < 30 ){
                            GenEventsIntRpArapPos->Fill(Centrality/2., weight);
                            GenEventsPtRpArapPos->Fill(ptGen, weight);
                        }
                        if ((rapLowRpA < rapGenCM) && (rapGenCM < rapHighRpA) && (ptGen > lowpTLow) \
								&& (ptGen < lowpTHigh) ){
                            GenEventsRapRpAlowpT->Fill(rapGenCM, weight);
                        }
                        if ((rapLowRpA < rapGenCM) && (rapGenCM < rapHighRpA) && (ptGen > highpTLow) \
								&& (ptGen < highpTHigh) ){
                            GenEventsRapRpAhighpT->Fill(rapGenCM, weight);
                        }
                    }
                    // Filling histograms for RpA and cross section study
                    else{
                    	// For asymmetric cross section
						if ((rapLow < rapGenCM) && (rapGenCM < rapHigh) && ptGen < 30 ){ 
							if(ispPb){
								GenEventsInt->Fill(Centrality/2., weight_XS);
								GenEventsPt->Fill(ptGen, weight_XS);
								GenEventsRap->Fill(rapGenCM, weight_XS);
							}
						    else{
						    	GenEventsInt->Fill(Centrality/2., weight);  
								GenEventsPt->Fill(ptGen, weight);  
								GenEventsRap->Fill(rapGenCM, weight);
							}
						}
						// Additionally, for RpA and symmetric cross section
						if ((rapLowRpA < rapGenCM) && (rapGenCM < rapHighRpA) && ptGen < 30 ){ 
							GenEventsIntRpA->Fill(Centrality/2., weight);
							GenEventsPtRpA->Fill(ptGen, weight);
							GenEventsRapRpA->Fill(rapGenCM, weight);
						}
					}
				}  // Filled Gen Histograms (Denominator)
			} // End Denominator/RECO Loop
		} // End Gen_QQ_Size > 0 Loop 
	} // End nentries (recorded events) loop


	// Plotting
	TGraphAsymmErrors *EffPtRpArapPos = new TGraphAsymmErrors(nPtBin);
	TGraphAsymmErrors *EffPtRpArapNeg = new TGraphAsymmErrors(nPtBin);
	TGraphAsymmErrors *EffIntRpArapPos = new TGraphAsymmErrors(1);
	TGraphAsymmErrors *EffIntRpArapNeg = new TGraphAsymmErrors(1);
	TGraphAsymmErrors *EffRapRpAlowpT = new TGraphAsymmErrors(nRapBin);
	TGraphAsymmErrors *EffRapRpAhighpT = new TGraphAsymmErrors(nRapBin);

	TGraphAsymmErrors *EffInt = new TGraphAsymmErrors(1);
	TGraphAsymmErrors *EffIntRpA = new TGraphAsymmErrors(1);
	TGraphAsymmErrors *EffPtRpA = new TGraphAsymmErrors(nPtBin);
	TGraphAsymmErrors *EffPt = new TGraphAsymmErrors(nPtBin);
	TGraphAsymmErrors *EffRap = new TGraphAsymmErrors(nRapBin);
	TGraphAsymmErrors *EffRapRpA = new TGraphAsymmErrors(nRapBin);

	TLatex tex_RpA2D;
	tex_RpA2D.SetTextAlign(12);
	tex_RpA2D.SetTextSize(0.038);

	// Declare graphs, set styles and make plots for differential RpA study
	if(isRpA2D){
	    EffIntRpArapPos->BayesDivide(RecoEventsIntRpArapPos, GenEventsIntRpArapPos);
	    EffIntRpArapPos->SetName("EffIntRpArapPos");
	    double IntValRpArapPos = EffIntRpArapPos->Eval(IntBin[0]);

	    EffIntRpArapNeg->BayesDivide(RecoEventsIntRpArapNeg, GenEventsIntRpArapNeg);
	    EffIntRpArapNeg->SetName("EffIntRpArapNeg");
	    double IntValRpArapNeg = EffIntRpArapNeg->Eval(IntBin[0]);

		// Lines used in pT plot for integrated efficiency in positive and negative y regions 
	    TLine *lIntRpArapPos = new TLine(0, IntValRpArapPos, 30, IntValRpArapPos);
	    lIntRpArapPos->SetLineStyle(2);   lIntRpArapPos->SetLineWidth(2);  lIntRpArapPos->SetLineColor(kBlue+2);

	    TLine *lIntRpArapNeg = new TLine(0, IntValRpArapNeg, 30, IntValRpArapNeg);
	    lIntRpArapNeg->SetLineStyle(2);   lIntRpArapNeg->SetLineWidth(2);  lIntRpArapNeg->SetLineColor(kBlue+2);

		TCanvas *c1 = new TCanvas("c1","c1",800,600);
		c1->SetRightMargin(1);
		c1->cd();

		EffPtRpArapPos->BayesDivide(RecoEventsPtRpArapPos, GenEventsPtRpArapPos);
		EffPtRpArapPos->SetName("EffPtRpArapPos");
		
		EffPtRpArapPos->SetMarkerSize(1.0);
		EffPtRpArapPos->SetMarkerColor(kRed);
		EffPtRpArapPos->SetMarkerStyle(20);
		
		EffPtRpArapPos->SetTitle("");
		EffPtRpArapPos->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		EffPtRpArapPos->GetXaxis()->SetTitle("p^{#varUpsilon}_{T} (GeV/c)");
		EffPtRpArapPos->GetYaxis()->SetRangeUser(0,1);
		EffPtRpArapPos->GetXaxis()->SetRangeUser(0.0, 30.0);
		EffPtRpArapPos->GetXaxis()->CenterTitle();
		EffPtRpArapPos->GetYaxis()->CenterTitle();
		EffPtRpArapPos->GetXaxis()->SetTitleOffset(0.92);
		EffPtRpArapPos->GetYaxis()->SetTitleOffset(0.92);
		EffPtRpArapPos->GetXaxis()->SetLabelSize(0.04);
		EffPtRpArapPos->GetYaxis()->SetLabelSize(0.04);
		
		EffPtRpArapPos->Draw("AP");
		lIntRpArapPos->Draw("sames");
		tex_RpA2D.DrawLatex(17, .55, ispPb ? "0 < y^{#varUpsilon}_{CM} < 1.93" : "|y^{#varUpsilon}_{CM}| < 1.93");
		CMS_lumi(c1,iPeriod, iPos);
		c1->Update();
		
		c1->SaveAs(Form("eff_XXXTAG/EfficiencyPtRpArapPos_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
		
		TCanvas *c2 = new TCanvas("c2","c2",800,600);
		c2->SetRightMargin(1);
		c2->cd();

		EffPtRpArapNeg->BayesDivide(RecoEventsPtRpArapNeg, GenEventsPtRpArapNeg);
		EffPtRpArapNeg->SetName("EffPtRpArapNeg");
		
		EffPtRpArapNeg->SetMarkerSize(1.0);
		EffPtRpArapNeg->SetMarkerColor(kRed);
		EffPtRpArapNeg->SetMarkerStyle(20);
		
		EffPtRpArapNeg->SetTitle("");
		EffPtRpArapNeg->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		EffPtRpArapNeg->GetXaxis()->SetTitle("p^{#varUpsilon}_{T} (GeV/c)");
		EffPtRpArapNeg->GetYaxis()->SetRangeUser(0,1);
		EffPtRpArapNeg->GetXaxis()->SetRangeUser(0.0, 30.0);
		EffPtRpArapNeg->GetXaxis()->CenterTitle();
		EffPtRpArapNeg->GetYaxis()->CenterTitle();
		EffPtRpArapNeg->GetXaxis()->SetTitleOffset(0.92);
		EffPtRpArapNeg->GetYaxis()->SetTitleOffset(0.92);
		EffPtRpArapNeg->GetXaxis()->SetLabelSize(0.04);
		EffPtRpArapNeg->GetYaxis()->SetLabelSize(0.04);
		
		EffPtRpArapNeg->Draw("AP");
		lIntRpArapNeg->Draw("sames");
		tex_RpA2D.DrawLatex(17, .55, ispPb ? "-1.93 < y^{#varUpsilon}_{CM} < 0" : "|y^{#varUpsilon}_{CM}| < 1.93");
		CMS_lumi(c2,iPeriod, iPos);
		c2->Update();
		
		c2->SaveAs(Form("eff_XXXTAG/EfficiencyPtRpArapNeg_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
		
		TCanvas *c3 = new TCanvas("c3","c3",800,600);
		c3->SetRightMargin(1);
		c3->cd();
		
		EffRapRpAlowpT->BayesDivide(RecoEventsRapRpAlowpT, GenEventsRapRpAlowpT);
		EffRapRpAlowpT->SetName("EffRapRpAlowpT");
		
		EffRapRpAlowpT->SetMarkerSize(1.0);
		EffRapRpAlowpT->SetMarkerColor(kRed);
		EffRapRpAlowpT->SetMarkerStyle(20);
		
		EffRapRpAlowpT->SetTitle("");
		EffRapRpAlowpT->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		if(ispPb){EffRapRpAlowpT->GetXaxis()->SetTitle("y^{#varUpsilon}_{CM}");}
		else{EffRapRpAlowpT->GetXaxis()->SetTitle("|y^{#varUpsilon}_{CM}|");}
		EffRapRpAlowpT->GetYaxis()->SetRangeUser(0,1);
		EffRapRpAlowpT->GetXaxis()->SetRangeUser(rapLowRpANeg,rapHighRpAPos);
		EffRapRpAlowpT->GetXaxis()->CenterTitle();
		EffRapRpAlowpT->GetYaxis()->CenterTitle();
		EffRapRpAlowpT->GetXaxis()->SetTitleOffset(0.92);
		EffRapRpAlowpT->GetYaxis()->SetTitleOffset(0.92);
		EffRapRpAlowpT->GetXaxis()->SetLabelSize(0.04);
		EffRapRpAlowpT->GetYaxis()->SetLabelSize(0.04);
		
		EffRapRpAlowpT->Draw("AP");
		tex_RpA2D.DrawLatex(0.3, .55, "p^{#varUpsilon}_{T} < 6 GeV");
		CMS_lumi(c3,iPeriod, iPos);
		c3->Update();
		
		c3->SaveAs(Form("eff_XXXTAG/EfficiencyRapRpAlowpT_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
		
		TCanvas *c4 = new TCanvas("c4","c4",800,600);
		c4->SetRightMargin(1);
		c4->cd();
		
		EffRapRpAhighpT->BayesDivide(RecoEventsRapRpAhighpT, GenEventsRapRpAhighpT);
		EffRapRpAhighpT->SetName("EffRapRpAhighpT");
		
		EffRapRpAhighpT->SetMarkerSize(1.0);
		EffRapRpAhighpT->SetMarkerColor(kRed);
		EffRapRpAhighpT->SetMarkerStyle(20);
		
		EffRapRpAhighpT->SetTitle("");
		EffRapRpAhighpT->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		if(ispPb){EffRapRpAhighpT->GetXaxis()->SetTitle("y^{#varUpsilon}_{CM}");}
		else{EffRapRpAhighpT->GetXaxis()->SetTitle("|y^{#varUpsilon}_{CM}|");}
		EffRapRpAhighpT->GetYaxis()->SetRangeUser(0,1);
		EffRapRpAhighpT->GetXaxis()->SetRangeUser(rapLowRpANeg,rapHighRpAPos);
		EffRapRpAhighpT->GetXaxis()->CenterTitle();
		EffRapRpAhighpT->GetYaxis()->CenterTitle();
		EffRapRpAhighpT->GetXaxis()->SetTitleOffset(0.92);
		EffRapRpAhighpT->GetYaxis()->SetTitleOffset(0.92);
		EffRapRpAhighpT->GetXaxis()->SetLabelSize(0.04);
		EffRapRpAhighpT->GetYaxis()->SetLabelSize(0.04);
		
		EffRapRpAhighpT->Draw("AP");
		tex_RpA2D.DrawLatex(0.3, .55, "6 GeV < p^{#varUpsilon}_{T} < 30 GeV");
		CMS_lumi(c4,iPeriod, iPos);
		c4->Update();
		
		c4->SaveAs(Form("eff_XXXTAG/EfficiencyRapRpAhighpT_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
	}

	// Declare graphs, set styles and make plots for RpA and Cross section study
	else {
		EffInt->BayesDivide(RecoEventsInt, GenEventsInt);
		EffInt->SetName("EffInt");
		double IntVal = EffInt->Eval(IntBin[0]);
		
		EffIntRpA->BayesDivide(RecoEventsIntRpA, GenEventsIntRpA);
		EffIntRpA->SetName("EffIntRpA");
		double IntValRpA = EffIntRpA->Eval(IntBin[0]);
		
		// Vertical line used in y plot to designate symmetric y region used tor RpA
		TLine *lylow;
		if(ispPb){
			lylow = new TLine(rapLowRpA, 0, rapLowRpA, 1);
			lylow->SetLineStyle(2);   lylow->SetLineWidth(2);  lylow->SetLineColor(kGreen+2);
		}

		// Line used in y plot for Integrated Efficiency value for entire y region available
		TLine *lInt = new TLine(rapLow, IntVal, rapHigh, IntVal);
		lInt->SetLineStyle(2);   lInt->SetLineWidth(2);  lInt->SetLineColor(kBlue+2); 
		
		// Line used in y plot for Integrated Efficiency value for symmetric y region
		TLine *lIntRpA = new TLine(rapLowRpA, IntValRpA, rapHighRpA, IntValRpA);
		lIntRpA->SetLineStyle(2);   lIntRpA->SetLineWidth(2);  lIntRpA->SetLineColor(kRed+2);
		
		// Line used in pT plot for Integrated Efficiency value for entire y region available
		TLine *lIntPt = new TLine(0, IntVal, 30, IntVal);
		lIntPt->SetLineStyle(2);   lIntPt->SetLineWidth(2);  lIntPt->SetLineColor(kBlue+2); 
		
		// Line used in pT plot for Integrated Efficiency value for symmetric y region used for RpA
		TLine *lIntPtRpA = new TLine(0, IntValRpA, 30, IntValRpA);
		lIntPtRpA->SetLineStyle(2);   lIntPtRpA->SetLineWidth(2);  lIntPtRpA->SetLineColor(kRed+2);
		
		TLegend *legyRpA = new TLegend(0.31,0.500,0.65,0.650);
		legyRpA->SetTextSize(0.038);
		legyRpA->AddEntry(lIntRpA,Form("%s, p^{#varUpsilon}_{T} < 30", Form("|y^{#varUpsilon}_{CM}| < %.2f", rapHighRpA)), "l");
		
		TLegend *legy = new TLegend(0.31,0.500,0.65,0.650);
		legy->SetTextSize(0.038);
		legy->AddEntry(lInt,Form("%s, p^{#varUpsilon}_{T} < 30", ispPb ? "-2.87 < y^{#varUpsilon}_{CM} < 1.93" \
					        : Form("|y^{#varUpsilon}_{CM}| < %.2f", rapHigh)),"l");

		TLegend *legpt = new TLegend(0.31,0.500,0.65,0.650);
		legpt->SetTextSize(0.038);
		legpt->AddEntry(lIntPt,Form("%s, p^{#varUpsilon}_{T} < 30", ispPb ? "-2.87 < y^{#varUpsilon}_{CM} < 1.93" \
				       		: Form("|y^{#varUpsilon}_{CM}| < %.2f", rapHigh)),"l");
		
		TLegend *legptRpA = new TLegend(0.31,0.500,0.65,0.650);
		legptRpA->SetTextSize(0.038);
		legptRpA->AddEntry(lIntPtRpA,Form("%s, p^{#varUpsilon}_{T} < 30", Form("|y^{#varUpsilon}_{CM}| < %.2f", rapHighRpA)), "l");
		
		//----------Pt---------//
		TCanvas *c2 = new TCanvas("c2","c2",800,600);
		c2->SetRightMargin(1);
		c2->cd();
		
		EffPt->BayesDivide(RecoEventsPt, GenEventsPt);
		EffPt->SetName("EffPt");
		
		EffPt->SetMarkerSize(1.0);
		EffPt->SetMarkerColor(kRed);
		EffPt->SetMarkerStyle(20);
		
		EffPt->SetTitle("");
		EffPt->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		EffPt->GetXaxis()->SetTitle("p^{#varUpsilon}_{T} (GeV/c)");
		EffPt->GetYaxis()->SetRangeUser(0,1);
		EffPt->GetXaxis()->SetRangeUser(0.0, 30.0);
		EffPt->GetXaxis()->CenterTitle();
		EffPt->GetYaxis()->CenterTitle();
		EffPt->GetXaxis()->SetTitleOffset(0.92);
		EffPt->GetYaxis()->SetTitleOffset(0.92);
		EffPt->GetXaxis()->SetLabelSize(0.04);
		EffPt->GetYaxis()->SetLabelSize(0.04);
		
		EffPt->Draw("AP");
		lIntPt->Draw("sames");
		legpt->Draw("sames");
		CMS_lumi(c2,iPeriod, iPos);
		c2->Update();
		
		c2->SaveAs(Form("eff_XXXTAG/EfficiencyPt_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
		
		//----------PtRpA----------//
		TCanvas *c1 = new TCanvas("c1","c1",800,600);
		c1->SetRightMargin(1);
		c1->cd();
		
		EffPtRpA->BayesDivide(RecoEventsPtRpA, GenEventsPtRpA);
		EffPtRpA->SetName("EffPtRpA");
		
		EffPtRpA->SetMarkerSize(1.0);
		EffPtRpA->SetMarkerColor(kRed);
		EffPtRpA->SetMarkerStyle(20);
		
		EffPtRpA->SetTitle("");
		EffPtRpA->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		EffPtRpA->GetXaxis()->SetTitle("p^{#varUpsilon}_{T} (GeV/c)");
		EffPtRpA->GetYaxis()->SetRangeUser(0,1);
		EffPtRpA->GetXaxis()->SetRangeUser(0.0, 30.0);
		EffPtRpA->GetXaxis()->CenterTitle();
		EffPtRpA->GetYaxis()->CenterTitle();
		EffPtRpA->GetXaxis()->SetTitleOffset(0.92);
		EffPtRpA->GetYaxis()->SetTitleOffset(0.92);
		EffPtRpA->GetXaxis()->SetLabelSize(0.04);
		EffPtRpA->GetYaxis()->SetLabelSize(0.04);
		
		EffPtRpA->Draw("AP");
		lIntPtRpA->Draw("sames");
		legptRpA->Draw("sames");
		CMS_lumi(c2,iPeriod, iPos);
		c1->Update();
		
		c1->SaveAs(Form("eff_XXXTAG/EfficiencyPtRpA_%dS_%s_TAG.png",oniaMode, ispPb ? "pPb" : "PP"));
		
		//------------Rap----------//
		TCanvas *c3 = new TCanvas("c3","c3",800,600);
		c3->SetRightMargin(1);
		c3->cd();
		
		EffRap->BayesDivide(RecoEventsRap, GenEventsRap);
		EffRap->SetName("EffRap");
		
		EffRap->SetMarkerSize(1.0);
		EffRap->SetMarkerColor(kRed);
		EffRap->SetMarkerStyle(20);
		
		EffRap->SetTitle("");
		EffRap->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
		if(ispPb){EffRap->GetXaxis()->SetTitle("y^{#varUpsilon}_{CM}");}
		else{EffRap->GetXaxis()->SetTitle("|y^{#varUpsilon}_{CM}|");}
		EffRap->GetYaxis()->SetRangeUser(0,1);
		EffRap->GetXaxis()->SetRangeUser(rapLow,rapHigh);
		EffRap->GetXaxis()->CenterTitle();
		EffRap->GetYaxis()->CenterTitle();
		EffRap->GetXaxis()->SetTitleOffset(0.92);
		EffRap->GetYaxis()->SetTitleOffset(0.92);
		EffRap->GetXaxis()->SetLabelSize(0.04);
		EffRap->GetYaxis()->SetLabelSize(0.04);
		
		EffRap->Draw("AP");
		lInt->Draw("sames");
		legy->Draw("sames");
		CMS_lumi(c3,iPeriod, iPos);
		c3->Update();
		
		c3->SaveAs(Form("eff_XXXTAG/EfficiencyRap_%dS_%s_TAG.png",oniaMode,ispPb ? "pPb" : "PP"));	

		//-----------Rap RpA---------//
	    TCanvas *c4 = new TCanvas("c4","c4",800,600);
	    c4->SetRightMargin(1);
	    c4->cd();

	    EffRapRpA->BayesDivide(RecoEventsRapRpA, GenEventsRapRpA);
	    EffRapRpA->SetName("EffRapRpA");

	    EffRapRpA->SetMarkerSize(1.0);
	    EffRapRpA->SetMarkerColor(kRed);
	    EffRapRpA->SetMarkerStyle(20);

	    EffRapRpA->SetTitle("");
	    EffRapRpA->GetYaxis()->SetTitle(Form("#varepsilon [#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
	    if(ispPb){EffRapRpA->GetXaxis()->SetTitle("y^{#varUpsilon}_{CM}");}
	    else{EffRapRpA->GetXaxis()->SetTitle("|y^{#varUpsilon}_{CM}|");}
	    EffRapRpA->GetYaxis()->SetRangeUser(0,1);
	    EffRapRpA->GetXaxis()->SetRangeUser(rapLowRpA,rapHigh);
	    EffRapRpA->GetXaxis()->CenterTitle();
	    EffRapRpA->GetYaxis()->CenterTitle();
	    EffRapRpA->GetXaxis()->SetTitleOffset(0.92);
	    EffRapRpA->GetYaxis()->SetTitleOffset(0.92);
	    EffRapRpA->GetXaxis()->SetLabelSize(0.04);
	    EffRapRpA->GetYaxis()->SetLabelSize(0.04);

	    EffRapRpA->Draw("AP");
	    lIntRpA->Draw("sames");
	    legyRpA->Draw("sames");
	    CMS_lumi(c4,iPeriod, iPos);
	    c4->Update();

	    c4->SaveAs(Form("eff_XXXTAG/EfficiencyRapRpA_%dS_%s_TAG.png",oniaMode,ispPb ? "pPb" : "PP"));
	}
	// Done making all efficiency plots

	// Writing efficiencies to root files
	TFile* MyFileEff;
	MyFileEff = new TFile(Form("eff_XXXTAG/Eff_%s_%dS_TAG.root","XXX",oniaMode), "Recreate");
	hCrossCheck->Write();

	// Differential RpA
	if(isRpA2D){
		RecoEventsPtRpArapPos->Write();
		RecoEventsPtRpArapNeg->Write();
		RecoEventsIntRpArapPos->Write();
		RecoEventsIntRpArapNeg->Write();
		RecoEventsRapRpAlowpT->Write();
		RecoEventsRapRpAhighpT->Write();

		GenEventsPtRpArapPos->Write();
		GenEventsPtRpArapNeg->Write();
		GenEventsIntRpArapPos->Write();
		GenEventsIntRpArapNeg->Write();
		GenEventsRapRpAlowpT->Write();
		GenEventsRapRpAhighpT->Write();

		EffPtRpArapPos->Write();
		EffPtRpArapNeg->Write();
		EffIntRpArapPos->Write();
		EffIntRpArapNeg->Write();
		EffRapRpAlowpT->Write();
		EffRapRpAhighpT->Write();
	}
	// RpA and cross section
	else {
		RecoEventsInt->Write();
		RecoEventsPt->Write();
		RecoEventsRap->Write();
		GenEventsInt->Write();
		GenEventsPt->Write();
		GenEventsRap->Write();

		EffPt->Write();
		EffRap->Write();
		EffInt->Write();

		RecoEventsIntRpA->Write();
		RecoEventsPtRpA->Write();
		GenEventsIntRpA->Write();
		GenEventsPtRpA->Write();

		EffPtRpA->Write();
		EffIntRpA->Write();
	}

	MyFileEff->Close();

	// Writing out efficiencies to terminal
	cout << Form(" %dS %s ",oniaMode, ispPb ? "pPb" : "pp") << endl;

	// Differential RpA
	if(isRpA2D){
    	cout << "2D RpA: " << endl;
        for (Int_t i = 0; i < (nPtBin); i++){
			cout << setprecision(3) << fixed << EffPtRpArapPos->Eval(ptBin_arr[i]) << endl;
		}

		cout << "" << endl;
        cout << setprecision(3) << fixed << EffIntRpArapPos->Eval(IntBin[0]) << endl;
		cout << "" << endl;
		cout << "" << endl;
		cout << "" << endl;

        for (Int_t i = 0; i < (nPtBin); i++){
			cout << setprecision(3) << fixed << EffPtRpArapNeg->Eval(ptBin_arr[i]) << endl;
		}

		cout << "" << endl;
		cout << setprecision(3) << fixed << EffIntRpArapNeg->Eval(IntBin[0]) << endl;
		cout << "" << endl;
		cout << "" << endl;
		cout << "" << endl;

        if(ispPb){
			for (Int_t i = 0; i < (nRapBin); i++){
				cout << setprecision(3) << fixed << EffRapRpAlowpT->Eval(rapBin_arr[i]) << endl;
			}
		}
	    else{
            for (Int_t i = (nRapBin-1) ; i >= 0; i--){
            	cout << setprecision(3) << fixed << EffRapRpAlowpT->Eval(rapBin_arr[i]) << endl;
            }
			for (Int_t i = 0; i < (nRapBin); i++){
				cout << setprecision(3) << fixed << EffRapRpAlowpT->Eval(rapBin_arr[i]) << endl;
			}
		}

		cout << "" << endl;
		cout << "" << endl;
		cout << "" << endl;

		if(ispPb){
			for (Int_t i = 0; i < (nRapBin); i++){
				cout << setprecision(3) << fixed << EffRapRpAhighpT->Eval(rapBin_arr[i]) << endl;
			}
		}
		else{
            for (Int_t i = (nRapBin-1) ; i >= 0; i--){
            	cout << setprecision(3) << fixed << EffRapRpAhighpT->Eval(rapBin_arr[i]) << endl;
            }
			for (Int_t i = 0; i < (nRapBin); i++){
				cout << setprecision(3) << fixed << EffRapRpAhighpT->Eval(rapBin_arr[i])  <<endl;
			}
		}
	}

	// RpA and cross section
	else{
		if(ispPb){
			cout << "Asymmetric y region: " << endl;
	       	for (Int_t i = 0; i < (nPtBin); i++){
				cout << setprecision(3) << fixed << EffPt->Eval(ptBin_arr[i]) << endl;
			}
			cout << "" << endl;
			for (Int_t i = 0; i < (nRapBin); i++){
				cout << setprecision(3) << fixed << EffRap->Eval(rapBin_arr[i]) << endl;
			}
			cout << "" << endl;
			cout << setprecision(3) << fixed << EffInt->Eval(IntBin[0]) << endl;
		}

		cout << "" << endl;
		cout << "Symmetric region: " << endl;
        for (Int_t i = 0; i < (nPtBin); i++){
			cout << setprecision(3) << fixed << EffPtRpA->Eval(ptBin_arr[i]) << endl;
        }
		cout << "" << endl;
		if(ispPb){
        	for (Int_t i = 1; i < (nRapBin); i++){	  
				cout << setprecision(3) << fixed << EffRap->Eval(rapBin_arr[i]) << endl;
        	}
		}
        else{
            for (Int_t i = (nRapBin-1) ; i >= 0; i--){ 
            	cout << setprecision(3) << fixed << EffRap->Eval(rapBin_arr[i]) << endl;
            }
        for (Int_t i = 0; i < (nRapBin); i++){       
				cout << setprecision(3) << fixed << EffRap->Eval(rapBin_arr[i]) << endl;
        	}
        }
		cout << "" << endl;
		cout << setprecision(3) << fixed << EffIntRpA->Eval(IntBin[0]) << endl;
	}

	// Close files containing weighting functions
    PtReweightFunctions->Close();

}  // End of main function



// Functions
bool PtCut(TLorentzVector* Muon){
        if (Muon->Pt() < muonPtCut){ return false; }
        else return true;
}

bool IsAccept(TLorentzVector* Muon){
    return ((Muon->Pt() > muonPtCut) && (fabs(Muon->Eta())<2.4));
}


bool MassCut(TLorentzVector* DiMuon, double LowM, double HighM){
        if (DiMuon->M() < LowM){ return false; }
        if (DiMuon->M() > HighM){ return false; }
        return true;
}

double PtReweight(TLorentzVector* DiMuon, TF1 *Pt_ReWeights){
        double pT = (DiMuon->Pt());
        return Pt_ReWeights->Eval(pT);
}

// TnP weight functions
// pp Nominal correction
double weight_tp_pp(double pt, double eta)
{
      double trg_SF = tnp_weight_trg_pp(pt, eta, 0);
      double trk_SF = tnp_weight_trk_pp(0);

      return trg_SF * trk_SF;
}
// pp Trigger binned
double weight_tp_pp_binned(double pt, double eta)
{
      double trg_SF = tnp_weight_trg_pp(pt, eta, -10);
      double trk_SF = tnp_weight_trk_pp(0);

      return trg_SF * trk_SF;
}
// pp Trigger Up and Down
double sys_SF_tp_pp_trigger(double pt, double eta, int idx_variation)
{
          double trg_sys_SF = tnp_weight_trg_pp(pt, eta, idx_variation);
          double trk_sys_SF = tnp_weight_trk_pp(0);

          return trg_sys_SF * trk_sys_SF ;
}
// pp Tracking Up and Down
double sys_SF_tp_pp_tracking(double pt, double eta, int idx_variation)
{
          double trg_sys_SF = tnp_weight_trg_pp(pt, eta, 0);
          double trk_sys_SF = tnp_weight_trk_pp(idx_variation);

          return trg_sys_SF * trk_sys_SF ;
}
// pp MuID Up and Down
double sys_SF_tp_pp_muid(double pt, double eta, int idx_variation)
{
      double trg_sys_SF = tnp_weight_trg_pp(pt, eta, 0);
      double trk_sys_SF = tnp_weight_trk_pp(0);
      double muid_sys_SF = tnp_weight_muid_pp(pt, eta, idx_variation);

      return trg_sys_SF * trk_sys_SF * muid_sys_SF ;
}
// pp STA Up and Down
double sys_SF_tp_pp_sta(double pt, double eta, int idx_variation)
{
      double trg_sys_SF = tnp_weight_trg_pp(pt, eta, 0);
      double trk_sys_SF = tnp_weight_trk_pp(0);
      double sta_sys_SF = tnp_weight_sta_pp(pt, eta, idx_variation);

      return trg_sys_SF * trk_sys_SF * sta_sys_SF ;
}
// pPb Nominal correction
double weight_tp_pPb(double mupt1,double mupt2,double mueta1, double mueta2)
{
            TF1* hw1;
        TF1* hw2;
    
        if (  TMath::Abs(mueta1) < 0.9 )      hw1 = hTnp_pa_eta0_09;
        else if ( TMath::Abs(mueta1) < 1.2 )  hw1 = hTnp_pa_eta09_12;
        else if ( TMath::Abs(mueta1) < 1.6 )  hw1 = hTnp_pa_eta12_16;
        else if ( TMath::Abs(mueta1) < 2.1 )  hw1 = hTnp_pa_eta16_21;
        else                                  hw1 = hTnp_pa_eta21_24;
        if (  TMath::Abs(mueta2) < 0.9 )      hw2 = hTnp_pa_eta0_09;
        else if ( TMath::Abs(mueta2) < 1.2 )  hw2 = hTnp_pa_eta09_12;
        else if ( TMath::Abs(mueta2) < 1.6 )  hw2 = hTnp_pa_eta12_16;
        else if ( TMath::Abs(mueta2) < 2.1 )  hw2 = hTnp_pa_eta16_21;
        else                                  hw2 = hTnp_pa_eta21_24;

        double tnpWeightMu1 = hw1->Eval(mupt1);
        double tnpWeightMu2 = hw2->Eval(mupt2);

        return tnpWeightMu1 * tnpWeightMu2;
}
// pPb systematic Up and Down (MuID+Trigger combined)
double sys_SF_tp_pPb(double mupt1,double mupt2,double mueta1, double mueta2, bool isSysUp)
{
    TF1* hw1;
    TF1* hw2;
        
    if(isSysUp){
                if (  TMath::Abs(mueta1) < 0.9 )      hw1 = hTnp_sys_up_pa_eta0_09;
                else if ( TMath::Abs(mueta1) < 1.2 )  hw1 = hTnp_sys_up_pa_eta09_12;
                else if ( TMath::Abs(mueta1) < 1.6 )  hw1 = hTnp_sys_up_pa_eta12_16;
                else if ( TMath::Abs(mueta1) < 2.1 )  hw1 = hTnp_sys_up_pa_eta16_21;
                else                                  hw1 = hTnp_sys_up_pa_eta21_24;
                if (  TMath::Abs(mueta2) < 0.9 )      hw2 = hTnp_sys_up_pa_eta0_09;
                else if ( TMath::Abs(mueta2) < 1.2 )  hw2 = hTnp_sys_up_pa_eta09_12;
                else if ( TMath::Abs(mueta2) < 1.6 )  hw2 = hTnp_sys_up_pa_eta12_16;
                else if ( TMath::Abs(mueta2) < 2.1 )  hw2 = hTnp_sys_up_pa_eta16_21;
                else                                  hw2 = hTnp_sys_up_pa_eta21_24;
    }
    else{
                if (  TMath::Abs(mueta1) < 0.9 )      hw1 = hTnp_sys_down_pa_eta0_09;
                else if ( TMath::Abs(mueta1) < 1.2 )  hw1 = hTnp_sys_down_pa_eta09_12;
                else if ( TMath::Abs(mueta1) < 1.6 )  hw1 = hTnp_sys_down_pa_eta12_16;
                else if ( TMath::Abs(mueta1) < 2.1 )  hw1 = hTnp_sys_down_pa_eta16_21;
                else                                  hw1 = hTnp_sys_down_pa_eta21_24;
                if (  TMath::Abs(mueta2) < 0.9 )      hw2 = hTnp_sys_down_pa_eta0_09;
                else if ( TMath::Abs(mueta2) < 1.2 )  hw2 = hTnp_sys_down_pa_eta09_12;
                else if ( TMath::Abs(mueta2) < 1.6 )  hw2 = hTnp_sys_down_pa_eta12_16;
                else if ( TMath::Abs(mueta2) < 2.1 )  hw2 = hTnp_sys_down_pa_eta16_21;
                else                                  hw2 = hTnp_sys_down_pa_eta21_24;      
    }

                double tnpSysWeightMu1 = hw1->Eval(mupt1);
                double tnpSysWeightMu2 = hw2->Eval(mupt2);

        return tnpSysWeightMu1 * tnpSysWeightMu2;   
}
