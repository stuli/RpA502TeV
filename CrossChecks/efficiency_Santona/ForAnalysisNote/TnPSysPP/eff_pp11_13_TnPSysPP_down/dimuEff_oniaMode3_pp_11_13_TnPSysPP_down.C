#include "effCommon.h"
#include "tnp_weight_pp.h"

const double muonPtCut = 4.0;

TFile* fTnp_pa = new TFile("output_official_5eta_cutG_all_nominal_v3.root","READ");
TF1* hTnp_pa_eta0_09 = (TF1*)fTnp_pa->Get("func_1");
TF1* hTnp_pa_eta09_12 = (TF1*)fTnp_pa->Get("func_2");	
TF1* hTnp_pa_eta12_16 = (TF1*)fTnp_pa->Get("func_3");
TF1* hTnp_pa_eta16_21 = (TF1*)fTnp_pa->Get("func_4");
TF1* hTnp_pa_eta21_24 = (TF1*)fTnp_pa->Get("func_5");


//Ratio Error Propogation
double RError(double A, double eA, double B, double eB){
	double f=A/B;
	double fA=eA/A;
	double fB=eB/B;
	double eR=  f*sqrt( (fA*fA + fB*fB )) ;
	return eR;
}

//Product Error Propogation
double PError(double A, double eA, double B, double eB){
	double f=A*B;
	double fA=eA/A;
	double fB=eB/B;
	double eR=  f*sqrt( (fA*fA + fB*fB )) ;
	return eR;
}


bool PtCut(TLorentzVector* Muon){
        if (Muon->Pt() < muonPtCut){ return false; }
        else return true;
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


double weight_tp_pp(double pt, double eta)
{
      double trg_SF = tnp_weight_trg_pp(pt, eta, 0);
      double trk_SF = tnp_weight_trk_pp(0);

      return trg_SF * trk_SF;
}

double sys_SF_tp_pp(double pt, double eta, int idx_variation)
{
      double trg_sys_SF = tnp_weight_trg_pp(pt, eta, idx_variation);
      double trk_sys_SF = tnp_weight_trk_pp(idx_variation);
      double muid_sys_SF = tnp_weight_muid_pp(pt, eta, idx_variation);
      double sta_sys_SF = tnp_weight_sta_pp(pt, eta, idx_variation);      

      return trg_sys_SF * trk_sys_SF * muid_sys_SF * sta_sys_SF ; 
}

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

void dimuEff_oniaMode3_pp_11_13_TnPSysPP_down(
	int oniaMode = 3, //1 = 1S, 2 = 2S, 3 = 3S
	bool ispPb = 0 //true = pPb and false = pp
	){   

	int idx_nom = 0;
	int idx_sys_up = -1;
	int idx_sys_down = -2;

	setTDRStyle();

	if(ispPb){
	        TChain myTree_Data("myTree");
    		myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_1st_run_merged.root");
    		cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
    		myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_2nd_run_merged.root");
    		cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
	}
	else{
	        TChain myTree_Data("hionia/myTree");
		myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PP_Data/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root");
    		cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
	}

	TChain *myTree_pp = new TChain("hionia/myTree");
        
        if ((oniaMode == 1) && !ispPb){
		myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
	}

        if ((oniaMode == 2) && !ispPb){
        	myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }

        if ((oniaMode == 3) && !ispPb){
        	myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }


        TChain *myTree_pPb = new TChain("myTree");

    	if ((oniaMode == 1) && ispPb){
		myTree_pPb->Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_MC_Ups1S_PA_5TeV02_WithFSR_tuneD6T.root");
	}

        if ((oniaMode == 2) && ispPb){
        	myTree_pPb->Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_MC_Ups2S_PA_5TeV02_WithFSR_tuneD6T.root");
        }

        if ((oniaMode == 3) && ispPb){
        	myTree_pPb->Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_Ups3S_PA_5TeV02.root");
        }

        TChain *myTree;
        if(!ispPb){
        	myTree = (TChain*)myTree_pp;
        }
	else{
	        myTree = (TChain*)myTree_pPb;
	}
    
		
	Float_t         muMiDxy;
	Float_t         muMiDz;
	Int_t           muMiNPxlLayers;
	Int_t           muMiNTrkLayers;
	Float_t         muPlDxy;
	Float_t         muPlDz;
	Int_t           muPlNPxlLayers;
	Int_t           muPlNTrkLayers;
	Float_t         vProb;

const int nPtBins1s  = 6;  // double ptBin1s[nPtBins1s+1] = {0,2.5,5,8,15,30};
const int nPtBins2s  = 3;  // double ptBin2s[nPtBins2s+1] = {0,5,15,30};
const int nPtBins3s  = 2;  //double ptBin3s[nPtBins3s+1] = {0,5,15,30};

const int nYBins1S  = 9;   //double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 5;   //double yBin2S[nYBins2S+1] ={0, 1.2, 2.4};
const int nYBins3S  = 3;   //double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};

const int nYBins1Spp  = 4;
const int nYBins2Spp  = 2;
const int nYBins3Spp  = 1;

std::vector<double> ptBinEdges;
std::vector<double> ptBin;
std::vector<double> rapBinEdges;
std::vector<double> rapBin;

//declare the number of bins and assign bin edges
if(oniaMode ==1){
	nPtBin = nPtBins1s;
	if(ispPb){nRapBin = nYBins1S;}
	else{nRapBin = nYBins1Spp;}

	ptBinEdges = {0,2,4,6,9,12,30};
	ptBin = {1,3,5,7.5,10.5,21};
	if(ispPb){
	rapBinEdges = {-2.4, -1.67, -1.27, -0.87, -0.47, -0.07, 0.33, 0.73, 1.46, 2.4};
	rapBin = {-2.035, -1.47, -1.07, -0.67, -0.27, 0.13, 0.53, 1.095, 1.93};
	}
	else{
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
	rapBinEdges = {-2.4, -1.27, -0.47, 0.33, 1.46, 2.4};
	rapBin = { -1.835, -0.6585, -0.07, 0.895, 1.93 };
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
	rapBinEdges = {-2.4, -0.47, 1.46, 2.4};
	rapBin = { -1.435, 0.495, 1.93 };
        }
        else{
        rapBinEdges = {0, 1.93};
        rapBin = {0.965};
        }
}
	// The pPb rapidity cuts are for Run 1 Only. We only have MC for run 1.
	float rapLow = 0.0;
	float rapHigh = 1.93;
        if(ispPb){rapLow = -2.4;
        rapHigh = 2.4;
        }
                
	float  ptReWeight;
	double weighttp;

	float           IntBin[1] = { 50 };
	float		IntBinEdges[2] = { 0, 100 };
	float         	ptReweight = 0.0;

	float 		massLow = 0;
	float 		massHigh = 0;

	Int_t           Centrality;
	ULong64_t       HLTriggers;
	Int_t           Reco_QQ_size;
	Int_t           Reco_QQ_sign[45];   //[Reco_QQ_size]
	TClonesArray    *Reco_QQ_4mom;
	TClonesArray    *Reco_QQ_mupl_4mom;
	TClonesArray    *Reco_QQ_mumi_4mom;
	ULong64_t       Reco_QQ_trig[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_VtxProb[45];   //[Reco_QQ_size]
	Int_t           Reco_QQ_mupl_nPixWMea[45];   //[Reco_QQ_size]
	Int_t           Reco_QQ_mumi_nPixWMea[45];   //[Reco_QQ_size]
	Int_t           Reco_QQ_mupl_nTrkWMea[45];   //[Reco_QQ_size]
	Int_t           Reco_QQ_mumi_nTrkWMea[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_mupl_dxy[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_mumi_dxy[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_mupl_dz[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_mumi_dz[45];   //[Reco_QQ_size]


	Int_t           Gen_QQ_size;
	Int_t           Gen_QQ_sign[45];   //[Gen_QQ_size]
	TClonesArray    *Gen_QQ_4mom;
	TClonesArray    *Gen_QQ_mupl_4mom;
	TClonesArray    *Gen_QQ_mumi_4mom;


	TBranch        *b_Centrality;   //!
	TBranch        *b_HLTriggers;   //!
	TBranch        *b_Reco_QQ_size;   //!
	TBranch        *b_Reco_QQ_sign;   //!
	TBranch        *b_Reco_QQ_4mom;   //!
	TBranch        *b_Reco_QQ_mupl_4mom;   //!
	TBranch        *b_Reco_QQ_mumi_4mom;   //!
	TBranch        *b_Reco_QQ_trig;   //!
	TBranch        *b_Reco_QQ_VtxProb;   //!
	TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
	TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
	TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
	TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
	TBranch        *b_Reco_QQ_mupl_dxy;   //!
	TBranch        *b_Reco_QQ_mumi_dxy;   //!
	TBranch        *b_Reco_QQ_mupl_dz;   //!
	TBranch        *b_Reco_QQ_mumi_dz;   //!


	TBranch        *b_Gen_QQ_size;   //
	TBranch        *b_Gen_QQ_4mom;   //!
	TBranch        *b_Gen_QQ_mupl_4mom;   //!
	TBranch        *b_Gen_QQ_mumi_4mom;   //!


	//Set object pointer, Initialize
	Reco_QQ_4mom = 0;
	Reco_QQ_mupl_4mom = 0;
	Reco_QQ_mumi_4mom = 0;

	Gen_QQ_4mom = 0;
	Gen_QQ_mupl_4mom = 0;
	Gen_QQ_mumi_4mom = 0;

	myTree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	myTree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	myTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	myTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	myTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	myTree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
	myTree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
	myTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	myTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
	myTree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
	myTree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
	myTree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
	myTree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
	myTree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
	myTree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
	myTree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
	myTree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);


	myTree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
	myTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
	myTree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
	myTree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

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


	myTree->SetBranchStatus("Gen_QQ_size", 1);
	myTree->SetBranchStatus("Gen_QQ_4mom", 1);
	myTree->SetBranchStatus("Gen_QQ_mupl_4mom", 1);
	myTree->SetBranchStatus("Gen_QQ_mumi_4mom", 1);


	//convert bin and bin edge vectors to arrays to be used as parameter when delcaring TH1Ds
	double* ptBinEdges_arr = &ptBinEdges[0];
	double* ptBin_arr = &ptBin[0];
	double* rapBinEdges_arr = &rapBinEdges[0];
	double* rapBin_arr = &rapBin[0];

	TH1D  *RecoEventsInt = new TH1D("RecoEventsInt", "Reconstructed", 1, IntBinEdges);
	TH1D  *GenEventsInt = new TH1D("GenEventsInt", "Generated", 1, IntBinEdges);
	
	TH1D  *RecoEventsPt = new TH1D("RecoEventsPt", "Reconstructed", nPtBin, ptBinEdges_arr);
	TH1D  *GenEventsPt = new TH1D("GenEventsPt", "Generated", nPtBin, ptBinEdges_arr);

	TH1D  *RecoEventsRap = new TH1D("RecoEventsRap", "Reconstructed", nRapBin, rapBinEdges_arr);
	TH1D  *GenEventsRap = new TH1D("GenEventsRap", "Generated", nRapBin, rapBinEdges_arr);

	TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);

	RecoEventsInt->Sumw2();
	GenEventsInt->Sumw2();
	RecoEventsPt->Sumw2();
	GenEventsPt->Sumw2();
	RecoEventsRap->Sumw2();
	GenEventsRap->Sumw2();
	hCrossCheck->Sumw2();


	// Get pT Reweight functions
	std::string fmode="1";

	const char *f_name;
	if(!ispPb){
		if(oniaMode == 1){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_1s_1108.root";
		}else if(oniaMode ==2){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_2s_1108.root";
		}else{
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_3s_1108.root";
		}
	}else{
		if(oniaMode == 1){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_1s_1108.root";
		}else if(oniaMode ==2){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_2s_1108.root";
		}else{
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_3s_1108.root";
		}
	}

	TFile* PtReweightFunctions = new TFile(f_name, "Open");
        TF1* Pt_ReWeights = (TF1*)PtReweightFunctions->Get("dataMC_Ratio_norm");

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

	Long64_t nentries = myTree->GetEntries();
	cout << nentries << endl;

	for (Long64_t jentry = 0; jentry < nentries; jentry++){

		myTree->GetEntry(jentry);
		if(jentry%100000 == 0){
			cout<<"--Processing Event: "<<jentry<<endl;
		}

		//Numerator Loop RECO
		for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++){
			hCrossCheck->Fill(1);
			TLorentzVector *qq4mom = (TLorentzVector*)Reco_QQ_4mom->At(iQQ);
			TLorentzVector *mumi4mom = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iQQ);
			TLorentzVector *mupl4mom = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iQQ);

			//--Muid cuts for muon minus
			muMiDxy = Reco_QQ_mumi_dxy[iQQ];
			muMiDz = Reco_QQ_mumi_dz[iQQ];
			muMiNPxlLayers = Reco_QQ_mumi_nPixWMea[iQQ];
			muMiNTrkLayers = Reco_QQ_mumi_nTrkWMea[iQQ];

			//--Muid cuts for muon plus
			muPlDxy = Reco_QQ_mupl_dxy[iQQ];
			muPlDz = Reco_QQ_mupl_dz[iQQ];
			muPlNPxlLayers = Reco_QQ_mupl_nPixWMea[iQQ];
			muPlNTrkLayers = Reco_QQ_mupl_nTrkWMea[iQQ];

			// Vertex matching probability
			vProb = Reco_QQ_VtxProb[iQQ];

			bool mupl_cut = 0;
			bool mumi_cut = 0;
			bool trigL1Dmu = 0;
			bool PtCutPass = 0;
			bool MassCutPass = 0;

			//--Muon id cuts
                        if ( muPlNTrkLayers > 5 && muPlNPxlLayers > 0 && TMath::Abs(muPlDxy) < 0.3 && TMath::Abs(muPlDz) < 20 && vProb > 0.01){ mupl_cut = 1; }
                        if ( muMiNTrkLayers > 5 && muMiNPxlLayers > 0 && TMath::Abs(muMiDxy) < 0.3 && TMath::Abs(muMiDz) < 20){ mumi_cut = 1; }

			//check mass and pT cuts (acceptance)
			if (PtCut(mupl4mom) && PtCut(mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(qq4mom, massLow, massHigh);

			//check if trigger bit is matched to dimuon
			if ((HLTriggers & 1) == 1 && (Reco_QQ_trig[iQQ] & 1) == 1){ trigL1Dmu = 1; }

			// TnP weights only needed for reco
			float weight = 0;
			ptReweight = 0;
			weighttp=1.0;

			//getting reco pt
			float ptReco = 0;
			float rapReco = 0;
			ptReco = qq4mom->Pt();


			if(!ispPb){rapReco = TMath::Abs(qq4mom->Rapidity());}
			else{rapReco = qq4mom->Rapidity();}

			ptReweight = PtReweight(qq4mom, Pt_ReWeights);

			if(!ispPb){
//				weighttp = weight_tp_pp(mupl4mom->Pt(),mupl4mom->Eta()) * weight_tp_pp(mumi4mom->Pt(),mumi4mom->Eta());
//				weighttp = sys_SF_tp_pp(mupl4mom->Pt(), mupl4mom->Eta(), idx_sys_up) * sys_SF_tp_pp(mumi4mom->Pt(), mumi4mom->Eta(), idx_sys_up);
                              weighttp = sys_SF_tp_pp(mupl4mom->Pt(), mupl4mom->Eta(), idx_sys_down) * sys_SF_tp_pp(mumi4mom->Pt(), mumi4mom->Eta(), idx_sys_down);
			}else{
				weighttp = weight_tp_pPb(mupl4mom->Pt(),mumi4mom->Pt(),mupl4mom->Eta(), mumi4mom->Eta());
			}

			weight = ptReweight * weighttp ;
//			weight = weighttp;

			bool recoPass = 0;

			if (Reco_QQ_sign[iQQ] == 0 && mupl_cut && mumi_cut && trigL1Dmu){ recoPass = 1; } // acceptMu

			//filling RecoEvent Histo if passing
			if ((rapLow < rapReco) && (rapReco < rapHigh) && ptReco < 30 && Centrality < 200){
				if (recoPass == 1 && PtCutPass == 1 && MassCutPass == 1){
					RecoEventsInt->Fill(Centrality/2., weight);
					RecoEventsPt->Fill(ptReco, weight);
					RecoEventsRap->Fill(rapReco, weight);
				}
			}
		}


		//Denominator loop  GEN
		for (int iQQ = 0; iQQ < Gen_QQ_size; iQQ++){

			hCrossCheck->Fill(0);
			TLorentzVector *g_qq4mom = (TLorentzVector*)Gen_QQ_4mom->At(iQQ);
			TLorentzVector *g_mumi4mom = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iQQ);
			TLorentzVector *g_mupl4mom = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iQQ);

			bool acceptMu = 0;
			bool PtCutPass = 0;
			bool MassCutPass = 0;

			//check if muons are in acceptance
			if (PtCut(g_mupl4mom) && PtCut(g_mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(g_qq4mom, massLow, massHigh);

			float weight = 0;
			ptReweight = 0;

			//getting a pt gen value 
			float ptGen = 0;
			float rapGen = 0;
			ptGen = g_qq4mom->Pt();

			if(!ispPb){rapGen = TMath::Abs(g_qq4mom->Rapidity());}
			else{rapGen = g_qq4mom->Rapidity();}

			ptReweight = PtReweight(g_qq4mom, Pt_ReWeights);
			weight = ptReweight;
//			weight = 1.0;

			//fill GenEvent Histo Denominator if passing 
			if ((rapLow < rapGen) && (rapGen < rapHigh) && ptGen < 30 && Centrality < 200 ){
				if (PtCutPass == 1 && MassCutPass == 1){  //acceptMu == 1
					GenEventsInt->Fill(Centrality/2., weight);
					GenEventsPt->Fill(ptGen, weight);
					GenEventsRap->Fill(rapGen, weight);
				}
			}
		}

	}


// Plotting
//----------Pt
TCanvas *c2 = new TCanvas("c2","c2",800,600);
c2->SetRightMargin(1);
c2->cd();

TGraphAsymmErrors *EffPt = new TGraphAsymmErrors(nPtBin);
EffPt->BayesDivide(RecoEventsPt, GenEventsPt);
EffPt->SetName("EffPt");

EffPt->SetMarkerSize(2.0);
EffPt->SetMarkerColor(kRed);
EffPt->SetMarkerStyle(20);

EffPt->SetTitle("");
EffPt->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
EffPt->GetXaxis()->SetTitle("p^{#mu^{+}#mu^{-}}_{T} (GeV/c)");
EffPt->GetYaxis()->SetRangeUser(0,1);
EffPt->GetXaxis()->SetRangeUser(0.0, 30.0);
EffPt->GetXaxis()->CenterTitle();
EffPt->GetYaxis()->CenterTitle();
EffPt->GetXaxis()->SetTitleOffset(0.9);
EffPt->GetYaxis()->SetTitleOffset(1);

EffPt->Draw("AP");
CMS_lumi(c2,iPeriod, iPos);
c2->Update();

c2->SaveAs(Form("eff_pp11_13_TnPSysPP_down/EfficiencyPt_%dS_%s_11_13_TnPSysPP_down.png",oniaMode, ispPb ? "pPb" : "PP"));

//------------Rap
TCanvas *c3 = new TCanvas("c3","c3",800,600);
c3->SetRightMargin(1);
c3->cd();

TGraphAsymmErrors *EffRap = new TGraphAsymmErrors(nRapBin);
EffRap->BayesDivide(RecoEventsRap, GenEventsRap);
EffRap->SetName("EffRap");

EffRap->SetMarkerSize(2.0);
EffRap->SetMarkerColor(kRed);
EffRap->SetMarkerStyle(20);

EffRap->SetTitle("");
EffRap->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
if(ispPb){EffRap->GetXaxis()->SetTitle("y^{#mu^{+}#mu^{-}}_{lab}");}
else{EffRap->GetXaxis()->SetTitle("|y|^{#mu^{+}#mu^{-}}_{lab}");}
//EffRap->GetXaxis()->SetTitle("y^{\mu^{+}\mu^{-}}_{CM}");
EffRap->GetYaxis()->SetRangeUser(0,1);
EffRap->GetXaxis()->SetRangeUser(rapLow,rapHigh);
EffRap->GetXaxis()->CenterTitle();
EffRap->GetYaxis()->CenterTitle();
EffRap->GetXaxis()->SetTitleOffset(0.9);
EffRap->GetYaxis()->SetTitleOffset(1);

EffRap->Draw("AP");
CMS_lumi(c3,iPeriod, iPos);
c3->Update();

c3->SaveAs(Form("eff_pp11_13_TnPSysPP_down/EfficiencyRap_%dS_%s_11_13_TnPSysPP_down.png",oniaMode,ispPb ? "pPb" : "PP"));

//------------Int
TCanvas *c4 = new TCanvas("c4","c4",800,600);
c4->SetRightMargin(1);
c4->cd();

TGraphAsymmErrors *EffInt = new TGraphAsymmErrors(1);
EffInt->BayesDivide(RecoEventsInt, GenEventsInt);
EffInt->SetName("EffInt");

EffInt->SetMarkerSize(2.0);
EffInt->SetMarkerColor(kRed);
EffInt->SetMarkerStyle(20);

EffInt->SetTitle("");
EffInt->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
EffInt->GetXaxis()->SetTitle("Integrated bin");
EffInt->GetYaxis()->SetRangeUser(0,1);
EffInt->GetXaxis()->SetRangeUser(0.0,100);
EffInt->GetXaxis()->CenterTitle();
EffInt->GetYaxis()->CenterTitle();
EffInt->GetXaxis()->SetTitleOffset(1);
EffInt->GetYaxis()->SetTitleOffset(1);

EffInt->Draw("AP");
CMS_lumi(c4,iPeriod, iPos);
c4->Update();

c4->SaveAs(Form("eff_pp11_13_TnPSysPP_down/EfficiencyInt_%dS_%s_11_13_TnPSysPP_down.png",oniaMode, ispPb ? "pPb" : "PP"));


// Writing efficiencies to file
TFile* MyFileEff;
MyFileEff = new TFile(Form("eff_pp11_13_TnPSysPP_down/Eff_%s_%dS_11_13_TnPSysPP_down.root","pp",oniaMode), "Recreate");
RecoEventsInt->Write();
RecoEventsPt->Write();
RecoEventsRap->Write();
GenEventsInt->Write();
GenEventsPt->Write();
GenEventsRap->Write();
hCrossCheck->Write();
EffPt->Write();
EffRap->Write();
EffInt->Write();
MyFileEff->Close();

	// Writing out efficiencies
        for (Int_t i = 0; i < (nPtBin); i++){
        cout << "Pt" << EffPt->Eval(ptBin_arr[i]) << " , - " << EffPt->GetErrorYlow(i) << " , + " << EffPt->GetErrorYhigh(i) << endl;
	}
        for (Int_t i = 0; i < (nRapBin); i++){
        cout << "Rapidity" << EffRap->Eval(rapBin_arr[i]) << " , - " << EffRap->GetErrorYlow(i) << " , + " << EffRap->GetErrorYhigh(i) << endl;
        }
        cout << "Integrated" << EffInt->Eval(IntBin[0]) << " , - " << EffInt->GetErrorYlow(0) << " , + " << EffInt->GetErrorYhigh(0) << endl;


        PtReweightFunctions->Close();

}  // end void

