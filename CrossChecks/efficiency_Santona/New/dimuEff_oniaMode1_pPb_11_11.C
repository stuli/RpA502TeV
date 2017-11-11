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
		// 20160626 (off8M, tagpt5)

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
		//cout<<tnpWeight1*tnpWeight2<<endl;
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

void dimuEff_oniaMode1_pPb_11_11(
	int oniaMode = 1, //1 = 1S, 2 = 2S, 3 = 3S
	bool ispPb = 1 //true = pPb and false = pp
	){   

	int idx_nom = 0;
	int idx_sys_up = -1;
	int idx_sys_down = -2;

	setTDRStyle();

	TChain myTree_Data("hionia/myTree");
	if(ispPb){
    		myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_1st_run_merged.root");
    		cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
    		myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_2nd_run_merged.root");
    		cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
	}
	else{
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
//	Bool_t          muMiGoodMu;
	Float_t         muPlDxy;
	Float_t         muPlDz;
	Int_t           muPlNPxlLayers;
	Int_t           muPlNTrkLayers;
//	Bool_t          muPlGoodMu;
	Float_t         vProb;

const int nPtBins1s  = 6;  // double ptBin1s[nPtBins1s+1] = {0,2.5,5,8,15,30};
const int nPtBins2s  = 3;  // double ptBin2s[nPtBins2s+1] = {0,5,15,30};
const int nPtBins3s  = 2;  //double ptBin3s[nPtBins3s+1] = {0,5,15,30};

const int nYBins1S  = 9;   //double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 5;   //double yBin2S[nYBins2S+1] ={0, 1.2, 2.4};
const int nYBins3S  = 3;   //double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};


std::vector<double> ptBinEdges;
std::vector<double> ptBin;
std::vector<double> rapBinEdges;
std::vector<double> rapBin;
//std::vector<double> IntBinEdges;
//std::vector<double> IntBin;

//declare the number of bins and assign bin edges
if(oniaMode ==1){
	nPtBin = nPtBins1s;
	nRapBin = nYBins1S;
//	nNtracksBin = nNtracksBins1s;
//	nSumET_HFBin = nSumET_HFBins1s;

	ptBinEdges = {0,2,4,6,9,12,30};
	ptBin = {1,3,5,7.5,10.5,21};
	rapBinEdges = {-2.4, -1.67, -1.27, -0.87, -0.47, -0.07, 0.33, 0.73, 1.46, 2.4};
	rapBin = {-2.035, -1.47, -1.07, -0.67, -0.27, 0.13, 0.53, 1.095, 1.93};
/*	NtracksBinEdges = { 0,10,15,20,27,36,200 };
	NtracksBin = { 5,12.5,23.5,31.5,118};
	SumET_HFBinEdges = { 0,9,13,18,24,32,200 };
	SumET_HFBin = { 4.5,11,15.5,21,28,116};
// */
}
if(oniaMode ==2){
	nPtBin = nPtBins2s;
	nRapBin = nYBins2S;
//	nNtracksBin = nNtracksBins2s3s;
//	nSumET_HFBin = nSumET_HFBins2s3s;

	ptBinEdges = {0,4,9,30};
	ptBin = {2,6.5,19.5};
	rapBinEdges = {-2.4, -1.27, -0.47, 0.33, 1.46, 2.4};
	rapBin = { -1.835, -0.6585, -0.07, 0.895, 1.93 };
/*	NtracksBinEdges = { 0,12,20,31,200 };
	NtracksBin = { 6, 16, 25.5, 115.5};
	SumET_HFBinEdges = { 0,11,18,28,200 };
	SumET_HFBin = { 5.5, 14.5,23,114};
// */
}
if(oniaMode ==3){
	nPtBin = nPtBins3s;
	nRapBin = nYBins3S;
//	nNtracksBin  = nNtracksBins2s3s;
//	nSumET_HFBin  = nSumET_HFBins2s3s;

	ptBinEdges = {0.0,6.0,30.0};
	ptBin = {3.0,18.0};
	rapBinEdges = {-2.4, -0.47, 1.46, 2.4};
	rapBin = { -1.435, 0.495, 1.93 };
/*	NtracksBinEdges = { 0,12,20,31,200 };
	NtracksBin = { 6, 16, 25.5, 115.5};
	SumET_HFBinEdges = { 0,11,18,28,200 };
	SumET_HFBin = { 5.5, 14.5,23,114};
// */
}
	// These rapidity cuts are for Run 1 Only. We only have MC for run 1.
        float rapLow = -2.4;
        float rapHigh = 2.4;
                        
	float  ptReWeight;
	double weighttp;
    	double weighttpsta;

	float           IntBin[1] = { 50 };
	float		IntBinEdges[2] = { 0, 100 };
	float         	ptReweight = 0.0;

	float 		massLow = 0;
	float 		massHigh = 0;

	Int_t           Centrality;
//	Int_t			Ntracks;
//	Float_t			SumET_HF;
	ULong64_t       HLTriggers;
	Int_t           Reco_QQ_size;
	Int_t           Reco_QQ_sign[45];   //[Reco_QQ_size]
	TClonesArray    *Reco_QQ_4mom;
	TClonesArray    *Reco_QQ_mupl_4mom;
	TClonesArray    *Reco_QQ_mumi_4mom;
	ULong64_t       Reco_QQ_trig[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_VtxProb[45];   //[Reco_QQ_size]
//	Bool_t          Reco_QQ_mupl_isGoodMuon[45];   //[Reco_QQ_size]
//	Bool_t          Reco_QQ_mumi_isGoodMuon[45];   //[Reco_QQ_size]
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
/*	Float_t         Gen_QQ_VtxProb[45];   //[Gen_QQ_size]
	Bool_t          Gen_QQ_mupl_isGoodMuon[45];   //[Gen_QQ_size]
	Bool_t          Gen_QQ_mumi_isGoodMuon[45];   //[Gen_QQ_size]
	Int_t           Gen_QQ_mupl_nPixWMea[45];   //[Gen_QQ_size]
	Int_t           Gen_QQ_mumi_nPixWMea[45];   //[Gen_QQ_size]
	Int_t           Gen_QQ_mupl_nTrkWMea[45];   //[Gen_QQ_size]
	Int_t           Gen_QQ_mumi_nTrkWMea[45];   //[Gen_QQ_size]
	Float_t         Gen_QQ_mupl_dxy[45];   //[Gen_QQ_size]
	Float_t         Gen_QQ_mumi_dxy[45];   //[Gen_QQ_size]
	Float_t         Gen_QQ_mupl_dz[45];   //[Gen_QQ_size]
	Float_t         Gen_QQ_mumi_dz[45];   //[Gen_QQ_size]
/ */


//	TBranch        *b_SumET_HF;   //!
//	TBranch        *b_Ntracks;   //!
	TBranch        *b_Centrality;   //!
	TBranch        *b_HLTriggers;   //!
	TBranch        *b_Reco_QQ_size;   //!
	TBranch        *b_Reco_QQ_sign;   //!
	TBranch        *b_Reco_QQ_4mom;   //!
	TBranch        *b_Reco_QQ_mupl_4mom;   //!
	TBranch        *b_Reco_QQ_mumi_4mom;   //!
	TBranch        *b_Reco_QQ_trig;   //!
	TBranch        *b_Reco_QQ_VtxProb;   //!
//	TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
//	TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
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

//	myTree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
//	myTree->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
	myTree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	myTree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	myTree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	myTree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	myTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	myTree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
	myTree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
	myTree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	myTree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
//	myTree->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
//	myTree->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
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

//	myTree->SetBranchStatus("SumET_HF", 1);
//	myTree->SetBranchStatus("Ntracks", 1);
	myTree->SetBranchStatus("Centrality", 1);
	myTree->SetBranchStatus("HLTriggers", 1);
	myTree->SetBranchStatus("Reco_QQ_size", 1);
	myTree->SetBranchStatus("Reco_QQ_sign", 1);
	myTree->SetBranchStatus("Reco_QQ_4mom", 1);
	myTree->SetBranchStatus("Reco_QQ_mupl_4mom", 1);
	myTree->SetBranchStatus("Reco_QQ_mumi_4mom", 1);
	myTree->SetBranchStatus("Reco_QQ_trig", 1);
	myTree->SetBranchStatus("Reco_QQ_VtxProb", 1);
//	myTree->SetBranchStatus("Reco_QQ_mupl_isGoodMuon", 1);
//	myTree->SetBranchStatus("Reco_QQ_mumi_isGoodMuon", 1);
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
/*	double* NtracksBinEdges_arr = &NtracksBinEdges[0];
	double* NtracksBin_arr = &NtracksBin[0];
	double* SumET_HFBinEdges_arr = &SumET_HFBinEdges[0];
	double* SumET_HFBin_arr = &SumET_HFBin[0];
// */

	TH1D  *RecoEventsInt = new TH1D("RecoEventsInt", "Reconstructed", 1, IntBinEdges);
	TH1D  *GenEventsInt = new TH1D("GenEventsInt", "Generated", 1, IntBinEdges);
	
	TH1D  *RecoEventsPt = new TH1D("RecoEventsPt", "Reconstructed", nPtBin, ptBinEdges_arr);
	TH1D  *GenEventsPt = new TH1D("GenEventsPt", "Generated", nPtBin, ptBinEdges_arr);

	TH1D  *RecoEventsRap = new TH1D("RecoEventsRap", "Reconstructed", nRapBin, rapBinEdges_arr);
	TH1D  *GenEventsRap = new TH1D("GenEventsRap", "Generated", nRapBin, rapBinEdges_arr);

	TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);

	//cout<<"STILL WORKING"<<endl;
	
	//TH1D  *SumET_HF_MC = new TH1D("SumET_HF_MC", "SumET_HF_MC", nSumET_HFBin, SumET_HFBinEdges_arr);
	//TH1D  *SumET_HF_Data = new TH1D("SumET_HF_Data", "SumET_HF_Data", nSumET_HFBin, SumET_HFBinEdges_arr);
	//TH1D  *RecoEventsNtracks = new TH1D("RecoEventsNtracks", "ReconstructedNtracks", nNtracksBin,NtracksBinEdges_arr);
	//TH1D  *GenEventsNtracks = new TH1D("GenEventsNtracks", "GeneratedNtracks", nNtracksBin, NtracksBinEdges_arr);

	//TH1D  *RecoEventsSumET_HF = new TH1D("RecoEventsSumET_HF", "ReconstructedSumET_HF", nSumET_HFBin, SumET_HFBinEdges_arr);
	//TH1D  *GenEventsSumET_HF = new TH1D("GenEventsSumET_HF", "GeneratedSumET_HF", nSumET_HFBin,SumET_HFBinEdges_arr);

	//TH1D  *Ntracks_MC = new TH1D("Ntracks_MC", "Ntracks_MC", nNtracksBin, NtracksBinEdges_arr);
	//TH1D  *Ntracks_Data = new TH1D("Ntracks_Data", "Ntracks_Data", nNtracksBin, NtracksBinEdges_arr);

	//TH1D  *SumET_HF_Weights = new TH1D("SumET_HF_Weights", "SumET_HF_Weights", nSumET_HFBin, SumET_HFBinEdges_arr);
	//TH1D  *Ntracks_Weights = new TH1D("Ntracks_Weights", "Ntracks_Weights", nNtracksBin, NtracksBinEdges_arr);
	
	//cout<<"STILL WORKING"<<endl;

	//myTree->Draw("Ntracks>>Ntracks_MC");
	//myTree_Data.Draw("Ntracks>>Ntracks_Data");
	//myTree->Draw("SumET_HF>>SumET_HF_MC");
	//myTree_Data.Draw("SumET_HF>>SumET_HF_Data");

	//SumET_HF_Weights->Sumw2();
	//Ntracks_Weights->Sumw2();
	
	/*SumET_HF_Weights->Divide(SumET_HF_Data,SumET_HF_MC);
	TCanvas *preCan1 = new TCanvas("preCan1","preCan1",800,600);
	SumET_HF_Weights->SetTitle("SumET_HF Weights");
	SumET_HF_Weights->GetXaxis()->SetTitle("E_{T}(MC)");
	SumET_HF_Weights->GetYaxis()->SetTitle("E_{T}(Data)/E_{T}(MC)");
	SumET_HF_Weights->Draw();

	TF1 *f_HFWeights = new TF1("f_HFWeights","[0]*TMath::Erf([1]*(x+[2]))+[3]",0,140);
	f_HFWeights->SetParameters(11,.025,-30,11);
	SumET_HF_Weights->Fit(f_HFWeights);
	f_HFWeights->Draw("SAME");

	preCan1->SaveAs(Form("eff_pp11_11/HFWeights_%dS_%s_11_11.png",oniaMode,"pp"));

	Ntracks_Weights->Divide(Ntracks_Data,Ntracks_MC);
	TCanvas *preCan2 = new TCanvas("preCan2","preCan2",800,600);
	Ntracks_Weights->SetTitle("Ntracks Weights");
	Ntracks_Weights->GetXaxis()->SetTitle("N_{tracks}(MC)");
	Ntracks_Weights->GetYaxis()->SetTitle("N_{tracks}(Data)/N_{tracks}(MC)");
	Ntracks_Weights->Draw();

	TF1 *f_Ntracks = new TF1("f_Ntracks","[0]*TMath::Erf([1]*(x+[2]))+[3]",0,200);
	f_Ntracks->SetParameters(12.5,.025,-60,12.5);
	Ntracks_Weights->Fit(f_Ntracks);
	f_Ntracks->Draw("SAME");

	preCan2->SaveAs(Form("eff_pp11_11/NtracksWeights_%dS_%s_11_11.png",oniaMode,"pp"));*/

	//RecoEventsNtracks->Sumw2();
	//GenEventsNtracks->Sumw2();
	//RecoEventsSumET_HF->Sumw2();
	//GenEventsSumET_HF->Sumw2();
	RecoEventsInt->Sumw2();
	GenEventsInt->Sumw2();
	RecoEventsPt->Sumw2();
	GenEventsPt->Sumw2();
	RecoEventsRap->Sumw2();
	GenEventsRap->Sumw2();
	hCrossCheck->Sumw2();

	std::string fmode="1";

	const char *f_name;
	if(!ispPb){
		if(oniaMode == 1){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_1s_20170816.root";
		}else if(oniaMode ==2){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_2s_20170816.root";
		}else{
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_3s_20170816.root";
		}
	}else{
		if(oniaMode == 1){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_1s_20170816.root";
		}else if(oniaMode ==2){
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_2s_20170816.root";
		}else{
			f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_3s_20170816.root";
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
//			muMiGoodMu = Reco_QQ_mumi_isGoodMuon[iQQ];

			//--Muid cuts for muon plus
			muPlDxy = Reco_QQ_mupl_dxy[iQQ];
			muPlDz = Reco_QQ_mupl_dz[iQQ];
			muPlNPxlLayers = Reco_QQ_mupl_nPixWMea[iQQ];
			muPlNTrkLayers = Reco_QQ_mupl_nTrkWMea[iQQ];
//			muPlGoodMu = Reco_QQ_mupl_isGoodMuon[iQQ];

			// Vertex matching probability
			vProb = Reco_QQ_VtxProb[iQQ];

			bool mupl_cut = 0;
			bool mumi_cut = 0;
			//bool acceptMu = 0;
			bool trigL1Dmu = 0;
			bool PtCutPass = 0;
			bool MassCutPass = 0;

			//--Muon id cuts
/*			if ((muPlGoodMu == 1) && muPlNTrkLayers > 5 && muPlNPxlLayers > 0 && TMath::Abs(muPlDxy) < 0.3 && TMath::Abs(muPlDz) < 20 && vProb > 0.01){ mupl_cut = 1; }
			if ((muMiGoodMu == 1) && muMiNTrkLayers > 5 && muMiNPxlLayers > 0 && TMath::Abs(muMiDxy) < 0.3 && TMath::Abs(muMiDz) < 20){ mumi_cut = 1; }
// */
                        if ( muPlNTrkLayers > 5 && muPlNPxlLayers > 0 && TMath::Abs(muPlDxy) < 0.3 && TMath::Abs(muPlDz) < 20 && vProb > 0.01){ mupl_cut = 1; }
                        if ( muMiNTrkLayers > 5 && muMiNPxlLayers > 0 && TMath::Abs(muMiDxy) < 0.3 && TMath::Abs(muMiDz) < 20){ mumi_cut = 1; }

			//check if muons are in acceptance
			//if (IsAccept(mupl4mom) && IsAccept(mumi4mom)){ acceptMu = 1; }
			if (PtCut(mupl4mom) && PtCut(mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(qq4mom, massLow, massHigh);

			//check if trigger bit is matched to dimuon
			if ((HLTriggers & 1) == 1 && (Reco_QQ_trig[iQQ] & 1) == 1){ trigL1Dmu = 1; }

			// TnP weights only needed for reco
			float weight = 0;
			ptReweight = 0;
			weighttp=1.;
         		weighttpsta=1.;

			//getting reco pt
			float ptReco = 0;
			float rapReco = 0;
			ptReco = qq4mom->Pt();


//			rapReco = TMath::Abs(qq4mom->Rapidity());
			rapReco = qq4mom->Rapidity();

			ptReweight = PtReweight(qq4mom, Pt_ReWeights);
			//cout<<ptReweight<<endl;

			if(!ispPb){
				weighttp=weight_tp_pp(mupl4mom->Pt(),mupl4mom->Eta())*weight_tp_pp(mumi4mom->Pt(),mumi4mom->Eta());
			}else{
				weighttp = weight_tp_pPb(mupl4mom->Pt(),mumi4mom->Pt(),mupl4mom->Eta(), mumi4mom->Eta());
			}

			weight = ptReweight * weighttp ;

			bool recoPass = 0;

			if (Reco_QQ_sign[iQQ] == 0 && mupl_cut && mumi_cut && trigL1Dmu){ recoPass = 1; } // acceptMu

			//filling RecoEvent Histo if passing
			if (rapLow < rapReco < rapHigh && ptReco < 30 && Centrality < 200){
				if (recoPass == 1 && PtCutPass == 1 && MassCutPass == 1){
					//RecoEventsNtracks->Fill(Ntracks, weight*sumET_HFWeight*weighttp);
					//RecoEventsSumET_HF->Fill(SumET_HF, weight*ntracksWeight*weighttp);
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
			//if (IsAccept(g_mupl4mom) && IsAccept(g_mumi4mom)){ acceptMu = 1; }
			if (PtCut(g_mupl4mom) && PtCut(g_mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(g_qq4mom, massLow, massHigh);

			float weight = 0;
			ptReweight = 0;

			//getting a pt gen value 
			float ptGen = 0;
			float rapGen = 0;
			ptGen = g_qq4mom->Pt();

//			rapGen = TMath::Abs(g_qq4mom->Rapidity());
			rapGen = g_qq4mom->Rapidity();

			ptReweight = PtReweight(g_qq4mom, Pt_ReWeights);
			//cout<<ptReweight<<endl;
			weight = ptReweight;

			//fill GenEvent Histo Denominator if passing 
			if (rapLow < rapGen < rapHigh && ptGen < 30 && Centrality < 200 ){
				if (PtCutPass == 1 && MassCutPass == 1){  //acceptMu == 1
					//GenEventsNtracks->Fill(Ntracks, weight*sumET_HFWeight);
					//GenEventsSumET_HF->Fill(SumET_HF, weight*ntracksWeight);
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
EffPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
EffPt->GetYaxis()->SetRangeUser(0,1);
EffPt->GetXaxis()->SetRangeUser(0.0, 30.0);
EffPt->GetXaxis()->CenterTitle();
EffPt->GetYaxis()->CenterTitle();
EffPt->GetXaxis()->SetTitleOffset(1);
EffPt->GetYaxis()->SetTitleOffset(1);

EffPt->Draw("AP");
CMS_lumi(c2,iPeriod, iPos);
c2->Update();

c2->SaveAs(Form("eff_pPb11_11/EfficiencyPt_%dS_%s_11_11.png",oniaMode, ispPb ? "pPb" : "PP"));

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
EffRap->GetXaxis()->SetTitle("y_{lab}");
EffRap->GetYaxis()->SetRangeUser(0,1);
EffRap->GetXaxis()->SetRangeUser(rapLow,rapHigh);
EffRap->GetXaxis()->CenterTitle();
EffRap->GetYaxis()->CenterTitle();
EffRap->GetXaxis()->SetTitleOffset(1);
EffRap->GetYaxis()->SetTitleOffset(1);

EffRap->Draw("AP");
CMS_lumi(c3,iPeriod, iPos);
c3->Update();

c3->SaveAs(Form("eff_pPb11_11/EfficiencyRap_%dS_%s_11_11.png",oniaMode,ispPb ? "pPb" : "PP"));

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

c4->SaveAs(Form("eff_pPb11_11/EfficiencyInt_%dS_%s_11_11.png",oniaMode, ispPb ? "pPb" : "PP"));


// Writing efficiencies to file
TFile* MyFileEff;
MyFileEff = new TFile(Form("eff_pPb11_11/Eff_%s_%dS_11_11.root","pPb",oniaMode), "Recreate");
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

        //for (Int_t i = 0; i < (nNtracksBin); i++){
        //cout << "Ntracks" << EffNtracks->Eval(NtracksBin_arr[i]) << " , - " << EffNtracks->GetErrorYlow(i) << " , + " << EffNtracks->GetErrorYhigh(i) << endl;
        //}
        //for (Int_t i = 0; i < (nSumET_HFBin); i++){
        //cout << "SumET_HF" << EffSumET_HF->Eval(SumET_HFBin_arr[i]) << " , - " << EffSumET_HF->GetErrorYlow(i) << " , + " << EffSumET_HF->GetErrorYhigh(i) << endl;
        //}

}  // end void




/*double FindNtracksWeight(int Bin) {
	const int nbins = 200;
	const double Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
	return Ncoll[Bin];
}*/
