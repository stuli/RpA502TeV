#include "effCommon.h"


double GetWeight(int numTree,int oniaMode);
bool IsAccept(TLorentzVector *Muon);
double FindCenWeight(int Bin);
double RError(double A, double eA, double B, double eB);
double PError(double A, double eA, double B, double eB);
bool PtCut(TLorentzVector* Muon);
bool MassCut(TLorentzVector* DiMuon, double LowM, double HighM);
double PtReweight(TLorentzVector* DiMuon, double coefficient, double constant);


int  nPtBin;     
int  nRapBin;     
int nCenBin;
int nNtracksBin;
int nSumET_HFBin;
const double muonPtCut = 4.0;

double m1S_low = 8.0;
double m1S_high = 10.0;
double m2S_low = 8.563;
double m2S_high = 10.563;
double m3S_low = 8.895;
double m3S_high = 10.895;


int iPeriod = 5;
int iPos = 33;

// Need to fix rap acceptance and cuts...

void dimuEff_oniaMode1_pPb_Ntrackseta_gt25(
	int oniaMode = 1, //1 = 1S, 2 = 2S, 3 = 3S
	bool ispPb = true, //true = pPb and false = pp
	int Ntracks_RapHigh = 2.5,
	int SumET_HF_RapHigh = 5,
	int SumET_HF_RapLow = 2.9
	){   

	setTDRStyle();

	TChain myTree_Data("myTree");
    myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_1st_run_merged.root");
    cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;
    myTree_Data.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_Data/RD2013_pa_2nd_run_merged.root");
    cout<<"Entries in Data Tree = "<<myTree_Data.GetEntries()<<endl;

	if (ispPb){
	TChain myTree("myTree");  
	}

        if ((oniaMode == 1) && ispPb){
	myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_MC_Ups1S_PA_5TeV02_WithFSR_tuneD6T.root");
	}

        if ((oniaMode == 2) && ispPb){
        myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_MC_Ups2S_PA_5TeV02_WithFSR_tuneD6T.root");
        }

        if ((oniaMode == 3) && ispPb){
        myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2013_5TeV02_pPb/pPb_MC/OniaTree_MC_Ups3S_PA_5TeV02_WithFSR_tuneD6T.root");
        }


        if (!ispPb){
        TChain myTree("hionia/myTree");
        }

	if (oniaMode == 2 && !ispPb){
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
	}

	if (oniaMode == 1 && !ispPb){
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
	}

        if (oniaMode == 3 && !ispPb){
                myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }
/*

	if (oniaMode == 2 && isPbPb){
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");   
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups2SMM_ptUps2S_15_inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	}

	if ((oniaMode == 1) && isPbPb){
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");   
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_09_12_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_12_15_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
		myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups1SMM_ptUps_15_30_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
	}

        if (oniaMode == 3 && isPbPb){
                myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups3SMM_ptUps3S_00_03_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
                myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups3SMM_ptUps3S_03_06_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
                myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups3SMM_ptUps3S_06_09_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");
                myTree.Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/PbPb_MC_Official/OniaTree_Pythia8_Ups3SMM_ptUps3S_09_inf_Hydjet_MB_HINPbPbWinter16DR-75X_mcRun2_HeavyIon_v13-v1.root");

        }
// */


	Float_t         muMiDxy;
	Float_t         muMiDz;
	Int_t           muMiNPxlLayers;
	Int_t           muMiNTrkLayers;
	Bool_t          muMiGoodMu;
	Float_t         muPlDxy;
	Float_t         muPlDz;
	Int_t           muPlNPxlLayers;
	Int_t           muPlNTrkLayers;
	Bool_t          muPlGoodMu;
	Float_t         vProb;

const int nPtBins1s  = 6;  // double ptBin1s[nPtBins1s+1] = {0,2.5,5,8,15,30};
const int nPtBins2s  = 3;  // double ptBin2s[nPtBins2s+1] = {0,5,15,30};
const int nPtBins3s  = 2;  //double ptBin3s[nPtBins3s+1] = {0,5,15,30};

const int nYBins1S  = 10;   //double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 6;   //double yBin2S[nYBins2S+1] ={0, 1.2, 2.4};
const int nYBins3S  = 4;   //double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};

const int nCenBins1s2s = 9;
const int nCenBins3s = 4;

const int nNtracksBins1s = 6;
const int nNtracksBins2s3s = 4;

const int nSumET_HFBins1s = 6;
const int nSumET_HFBins2s3s = 4;


/*   float CenBinEdges[nCenBins3s+1] = {0,10,30,50,100};
nPtBin = nPtBins3s;
float ptBin[nPtBins3s] = {2.5,10,22.5};
        float          ptBinEdges[nPtBins3s + 1] ={0,5,15,30};
nRapBin = nYBins3S;
float          rapBin[nYBins3S] = { 0.6, 1.8 };
float          rapBinEdges[nYBins3S + 1] = { 0, 1.2, 2.4 };
nCenBin  = nCenBins3s;
        float           CenBin[nCenBins3s] = { 5, 20, 40, 75};
// */

if (oniaMode == 1){
nPtBin = nPtBins1s;
cout << "This is fine" << endl;
float ptBin[nPtBins1s] = {1,3,5,7.5,10.5,21};
        float          ptBinEdges[nPtBins1s + 1] ={0,2,4,6,9,12,30};
nRapBin = nYBins1S;
if (ispPb){
float          rapBin[nYBins1S] = {-1.765, -1.4, -1.0,-.6,-.2,.2,.6, 1.0, 1.4, 1.765};
float          rapBinEdges[nYBins1S + 1] = {-1.93, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 1.93};
}
if (!ispPb){
float          rapBin[nYBins1S] = {-1.765, -1.4, -1.0,-.6,-.2,.2,.6, 1.0, 1.4, 1.765};
float          rapBinEdges[nYBins1S + 1] = {-1.93, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 1.93};
}
nCenBin = nCenBins1s2s;
        float           CenBin[nCenBins1s2s] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65, 85 };
        float           CenBinEdges[nCenBins1s2s + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 100 };
nNtracksBin = nNtracksBins1s;
        float           NtracksBin[nNtracksBins1s] = { 5,12.5,23.5,31.5,82};
        float           NtracksBinEdges[nNtracksBins1s + 1] = { 0,10,15,20,27,36,200 };
nSumET_HFBin = nSumET_HFBins1s;
        float           SumET_HFBin[nSumET_HFBins1s] = { 4.5,11,15.5,21,28,54};
        float           SumET_HFBinEdges[nSumET_HFBins1s + 1] = { 0,9,13,18,24,32,140 };
}
if (oniaMode == 2){
nPtBin = nPtBins2s;
float ptBin[nPtBins2s] = {2,6.5,19.5};
        float          ptBinEdges[nPtBins2s + 1] ={0,4,9,30};
nRapBin = nYBins2S;
if (ispPb){
float          rapBin[nYBins2S] = { -1.765,-1.2, -.4,.4,1.2, 1.765 };
float          rapBinEdges[nYBins2S + 1] = {-1.93, -1.6, -0.8, 0, 0.8, 1.6, 1.93};
}
if (!ispPb){
float          rapBin[nYBins2S] = { -1.765,-1.2, -.4,.4,1.2, 1.765 };
float          rapBinEdges[nYBins2S + 1] = {-1.93, -1.6, -0.8, 0, 0.8, 1.6, 1.93};
}
nCenBin = nCenBins1s2s;
        float           CenBin[nCenBins1s2s] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65, 85 };
        float           CenBinEdges[nCenBins1s2s + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 100 };  
nNtracksBin = nNtracksBins2s3s;
        float           NtracksBin[nNtracksBins2s3s] = { 6, 16, 25.5, 84.5};
        float           NtracksBinEdges[nNtracksBins2s3s + 1] = { 0,12,20,31,200 };
nSumET_HFBin = nSumET_HFBins2s3s;
        float           SumET_HFBin[nSumET_HFBins2s3s] = { 5.5, 14.5,23,56};
        float           SumET_HFBinEdges[nSumET_HFBins2s3s + 1] = { 0,11,18,28,140 };
} 
if (oniaMode == 3){
nPtBin = nPtBins3s;
float ptBin[nPtBins3s] = {3.0,18.0};
        float          ptBinEdges[nPtBins3s + 1] ={0.0,6.0,30.0};
nRapBin = nYBins3S;
if (ispPb){
float          rapBin[nYBins3S] = { -1.565, -.6,.6, 1.565 };
float          rapBinEdges[nYBins3S + 1] = {-1.93, -1.2, 0, 1.2, 1.93};
}
if (!ispPb){
float          rapBin[nYBins3S] = { -1.565, -.6,.6, 1.565 };
float          rapBinEdges[nYBins3S + 1] = {-1.93, -1.2, 0, 1.2, 1.93};
}
nCenBin  = nCenBins3s;
        float           CenBin[nCenBins3s] = { 5, 20, 40, 75};
   float CenBinEdges[nCenBins3s+1] = {0,10,30,50,100};
nNtracksBin  = nNtracksBins2s3s;
        float           NtracksBin[nNtracksBins3s] = { 6, 16, 25.5, 84.5};
   float NtracksBinEdges[nNtracksBins2s3s+1] = { 0,12,20,31,200 };
nSumET_HFBin  = nSumET_HFBins2s3s;
        float           SumET_HFBin[nSumET_HFBins2s3s] = { 5.5, 14.5,23,56};
   float SumET_HFBinEdges[nSumET_HFBins2s3s+1] = { 0,11,18,28,140 };
cout << "this is fine" << endl;
} // */

                        if (ispPb){float rapLow = -1.93;
                                float rapHigh = 1.93;
                        }
                        else {float rapLow = -1.93;
                                float rapHigh = 1.93;
                        }


	float          ptWeight;
	float          centWeight;
	float 		   ntracksWeight;
	float 		   sumET_HFWeight;
//	float          ptBin[nPtBin] = { 2.5, 8.5, 21 };
//	float          ptBinEdges[nPtBin + 1] = { 0, 5, 12, 30 };
//	float          rapBin[nRapBin] = { 0.6, 1.8 };
//	float          rapBinEdges[nRapBin + 1] = { 0, 1.2, 2.4 };

//	float           CenBin[nCenBin] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65, 85 };
//	float    	CenBinEdges[nCenBin + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 100 };

	float           IntBin[1] = { 100 };
	float		IntBinEdges[2] = { 0, 100 };
	float          ptReweight = 0.0;

	float 		massLow = 0;
	float 		massHigh = 0;

	Int_t           Centrality;
	Int_t			Ntracks;
	Float_t			SumET_HF;
	ULong64_t       HLTriggers;
	Int_t           Reco_QQ_size;
	Int_t           Reco_QQ_sign[45];   //[Reco_QQ_size]
	TClonesArray    *Reco_QQ_4mom;
	TClonesArray    *Reco_QQ_mupl_4mom;
	TClonesArray    *Reco_QQ_mumi_4mom;
	ULong64_t       Reco_QQ_trig[45];   //[Reco_QQ_size]
	Float_t         Reco_QQ_VtxProb[45];   //[Reco_QQ_size]
	Bool_t          Reco_QQ_mupl_isGoodMuon[45];   //[Reco_QQ_size]
	Bool_t          Reco_QQ_mumi_isGoodMuon[45];   //[Reco_QQ_size]
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
	Float_t         Gen_QQ_VtxProb[45];   //[Gen_QQ_size]
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



	TBranch        *b_SumET_HF;   //!
	TBranch        *b_Ntracks;   //!
	TBranch        *b_Centrality;   //!
	TBranch        *b_HLTriggers;   //!
	TBranch        *b_Reco_QQ_size;   //!
	TBranch        *b_Reco_QQ_sign;   //!
	TBranch        *b_Reco_QQ_4mom;   //!
	TBranch        *b_Reco_QQ_mupl_4mom;   //!
	TBranch        *b_Reco_QQ_mumi_4mom;   //!
	TBranch        *b_Reco_QQ_trig;   //!
	TBranch        *b_Reco_QQ_VtxProb;   //!
	TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
	TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
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

	myTree.SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
	myTree.SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
	myTree.SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	myTree.SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	myTree.SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	myTree.SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	myTree.SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	myTree.SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
	myTree.SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
	myTree.SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	myTree.SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
//	myTree.SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
//	myTree.SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
	myTree.SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
	myTree.SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
	myTree.SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
	myTree.SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
	myTree.SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
	myTree.SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
	myTree.SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
	myTree.SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);



	myTree.SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
	myTree.SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
	myTree.SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
	myTree.SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

	myTree.SetBranchStatus("*", 0);

	myTree.SetBranchStatus("SumET_HF", 1);
	myTree.SetBranchStatus("Ntracks", 1);
	myTree.SetBranchStatus("Centrality", 1);
	myTree.SetBranchStatus("HLTriggers", 1);
	myTree.SetBranchStatus("Reco_QQ_size", 1);
	myTree.SetBranchStatus("Reco_QQ_sign", 1);
	myTree.SetBranchStatus("Reco_QQ_4mom", 1);
	myTree.SetBranchStatus("Reco_QQ_mupl_4mom", 1);
	myTree.SetBranchStatus("Reco_QQ_mumi_4mom", 1);
	myTree.SetBranchStatus("Reco_QQ_trig", 1);
	myTree.SetBranchStatus("Reco_QQ_VtxProb", 1);
//	myTree.SetBranchStatus("Reco_QQ_mupl_isGoodMuon", 1);
//	myTree.SetBranchStatus("Reco_QQ_mumi_isGoodMuon", 1);
	myTree.SetBranchStatus("Reco_QQ_mupl_nPixWMea", 1);
	myTree.SetBranchStatus("Reco_QQ_mumi_nPixWMea", 1);
	myTree.SetBranchStatus("Reco_QQ_mupl_nTrkWMea", 1);
	myTree.SetBranchStatus("Reco_QQ_mumi_nTrkWMea", 1);
	myTree.SetBranchStatus("Reco_QQ_mupl_dxy", 1);
	myTree.SetBranchStatus("Reco_QQ_mumi_dxy", 1);
	myTree.SetBranchStatus("Reco_QQ_mupl_dz", 1);
	myTree.SetBranchStatus("Reco_QQ_mumi_dz", 1);



	myTree.SetBranchStatus("Gen_QQ_size", 1);
	myTree.SetBranchStatus("Gen_QQ_4mom", 1);
	myTree.SetBranchStatus("Gen_QQ_mupl_4mom", 1);
	myTree.SetBranchStatus("Gen_QQ_mumi_4mom", 1);



       	TH1D  *RecoEvents = new TH1D("RecoEvents", "Reconstructed", ispPb ? nCenBin : 1, ispPb ? CenBinEdges : IntBinEdges);
	TH1D  *GenEvents = new TH1D("GenEvents", "Generated", ispPb ? nCenBin : 1, ispPb ? CenBinEdges : IntBinEdges);

	TH1D  *RecoEventsNtracks = new TH1D("RecoEventsNtracks", "ReconstructedNtracks", ispPb ? nNtracksBin : 1, ispPb ? NtracksBinEdges : IntBinEdges);
	TH1D  *GenEventsNtracks = new TH1D("GenEventsNtracks", "GeneratedNtracks", ispPb ? nNtracksBin : 1, ispPb ? NtracksBinEdges : IntBinEdges);

	TH1D  *RecoEventsSumET_HF = new TH1D("RecoEventsSumET_HF", "ReconstructedSumET_HF", ispPb ? nSumET_HFBin : 1, ispPb ? SumET_HFBinEdges : IntBinEdges);
	TH1D  *GenEventsSumET_HF = new TH1D("GenEventsSumET_HF", "GeneratedSumET_HF", ispPb ? nSumET_HFBin : 1, ispPb ? SumET_HFBinEdges : IntBinEdges);

	TH1D  *RecoEventsInt = new TH1D("RecoEventsInt", "Reconstructed", 1, IntBinEdges);
	TH1D  *GenEventsInt = new TH1D("GenEventsInt", "Generated", 1, IntBinEdges);

cout << ptBinEdges << endl;

	TH1D  *RecoEventsPt = new TH1D("RecoEventsPt", "Reconstructed", nPtBin, ptBinEdges);
	TH1D  *GenEventsPt = new TH1D("GenEventsPt", "Generated", nPtBin, ptBinEdges);

	TH1D  *RecoEventsRap = new TH1D("RecoEventsRap", "Reconstructed", nRapBin, rapBinEdges);
	TH1D  *GenEventsRap = new TH1D("GenEventsRap", "Generated", nRapBin, rapBinEdges);

	TH1D  *hCentrality = new TH1D("hCentrality", "Centrality distribution", 202, -1, 201);
	TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);

	TH1D  *hRecoEventsD = new TH1D("hRecoEventsD", "Reconstructed", ispPb ? nCenBin : 1, ispPb ? CenBinEdges : IntBinEdges);
	TH1D  *hGenEventsD = new TH1D("hGenEventsD", "Generated", ispPb ? nCenBin : 1, ispPb ? CenBinEdges : IntBinEdges);
	//cout<<"STILL WORKING"<<endl;
	TH1D  *SumET_HF_MC = new TH1D("SumET_HF_MC", "SumET_HF_MC", ispPb ? nSumET_HFBin : 1, ispPb ? SumET_HFBinEdges : IntBinEdges);
	TH1D  *SumET_HF_Data = new TH1D("SumET_HF_Data", "SumET_HF_Data", ispPb ? nSumET_HFBin : 1, ispPb ? SumET_HFBinEdges : IntBinEdges);

	TH1D  *Ntracks_MC = new TH1D("Ntracks_MC", "Ntracks_MC", ispPb ? nNtracksBin : 1, ispPb ? NtracksBinEdges : IntBinEdges);
	TH1D  *Ntracks_Data = new TH1D("Ntracks_Data", "Ntracks_Data", ispPb ? nNtracksBin : 1, ispPb ? NtracksBinEdges : IntBinEdges);

	TH1D  *SumET_HF_Weights = new TH1D("SumET_HF_Weights", "SumET_HF_Weights", ispPb ? nSumET_HFBin : 1, ispPb ? SumET_HFBinEdges : IntBinEdges);
	TH1D  *Ntracks_Weights = new TH1D("Ntracks_Weights", "Ntracks_Weights", ispPb ? nNtracksBin : 1, ispPb ? NtracksBinEdges : IntBinEdges);
	//cout<<"STILL WORKING"<<endl;

	myTree.Draw("Ntracks>>Ntracks_MC");
	myTree_Data.Draw("Ntracks>>Ntracks_Data");
	myTree.Draw("SumET_HF>>SumET_HF_MC");
	myTree_Data.Draw("SumET_HF>>SumET_HF_Data");

	SumET_HF_Weights->Sumw2();
	Ntracks_Weights->Sumw2();
	
	SumET_HF_Weights->Divide(SumET_HF_Data,SumET_HF_MC);
	TCanvas *preCan2 = new TCanvas("preCan1","preCan1",800,600);
	SumET_HF_Weights->SetTitle("#Sigma E_{T}^{HF} Weights");
	SumET_HF_Weights->GetXaxis()->SetTitle("#Sigma E_{T}^{HF}(MC)");
	SumET_HF_Weights->GetYaxis()->SetTitle("#Sigma E_{T}^{HF}(Data)/#Sigma E_{T}^{HF}(MC)");
	SumET_HF_Weights->Draw();
	preCan1->SaveAs(Form("eff_pPb/HFWeights_%dS_%s_HFetagt25.png",oniaMode,ispPb ? "pPb" : "PP"));
	
	Ntracks_Weights->Divide(Ntracks_Data,Ntracks_MC);
	TCanvas *preCan2 = new TCanvas("preCan2","preCan2",800,600);
	Ntracks_Weights->SetTitle("Ntracks Weights");
	Ntracks_Weights->GetXaxis()->SetTitle("Ntracks(MC)");
	Ntracks_Weights->GetYaxis()->SetTitle("Ntracks(Data)/Ntracks(MC)");
	Ntracks_Weights->Draw();
	preCan2->SaveAs(Form("eff_pPb/NtracksWeights_%dS_%s_HFetagt25.png",oniaMode,ispPb ? "pPb" : "PP"));

	RecoEvents->Sumw2();
	GenEvents->Sumw2();
	RecoEventsNtracks->Sumw2();
	GenEventsNtracks->Sumw2();
	RecoEventsSumET_HF->Sumw2();
	GenEventsSumET_HF->Sumw2();
	RecoEventsInt->Sumw2();
	GenEventsInt->Sumw2();
	RecoEventsPt->Sumw2();
	GenEventsPt->Sumw2();
	RecoEventsRap->Sumw2();
	GenEventsRap->Sumw2();

	TF1* f1SAA;
	TF1* f2SAA;
	TF1* f1Spp;
	TF1* f2Spp;
	TFile* ReweightFunctions = new TFile("dNdpT_ratio_tsallis_June7.root", "Open");

	f1SAA = (TF1*)ReweightFunctions->Get("f1sraa_test");
	f2SAA = (TF1*)ReweightFunctions->Get("f2sraa_test");
	f1Spp = (TF1*)ReweightFunctions->Get("f1srpp_test");
	f2Spp = (TF1*)ReweightFunctions->Get("f2srpp_test");
// */

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

	Long64_t nentries = myTree.GetEntries();
	cout << nentries << endl;

	for (Long64_t jentry = 0; jentry < nentries; jentry++){
		myTree.GetEntry(jentry);
		if(jentry%1000 == 0){
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
			vProb = Reco_QQ_VtxProb[iQQ];

			bool mupl_cut = 0;
			bool mumi_cut = 0;
			bool acceptMu = 0;
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
			if (IsAccept(mupl4mom) && IsAccept(mumi4mom)){ acceptMu = 1; }
			if (PtCut(mupl4mom) && PtCut(mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(qq4mom, massLow, massHigh);

			//check if trigger bit is matched to dimuon
			if ((HLTriggers & 1) == 1 && (Reco_QQ_trig[iQQ] & 1) == 1){ trigL1Dmu = 1; }

			//weights only needed for pPb
			float weight = 0;
			ptWeight = 0;
			centWeight = 0;
			ntracksWeight=0;
			sumET_HFWeight= 0;
			centWeight = FindCenWeight(Centrality);
			ntracksWeight= FindNtracksWeight(Ntracks,Ntracks_Weights);
			sumET_HFWeight= FindSumET_HFWeight(SumET_HF,SumET_HF_Weights);
			ptReweight = 0;

			//getting reco pt
			float ptReco = 0;
			float rapReco = 0;
			ptReco = qq4mom->Pt();
//			rapReco = TMath::Abs(qq4mom->Rapidity());
			rapReco = qq4mom->Rapidity();

			//getting the tree weight from pt generated MC bins

			//reweight from dn/dpt function
			int tNum = -1;
			//total weighting factor
			if (ispPb){
				//if (oniaMode == 1){ ptReweight = (f1SAA->Eval(ptReco)); }
				//if (oniaMode == 2){ ptReweight = (f2SAA->Eval(ptReco)); }
				//if (oniaMode == 3){ ptReweight = 1;}
/*				tNum = myTree.GetTreeNumber();
				ptReweight = 1;
				ptWeight = GetWeight(tNum, oniaMode);
				weight = centWeight*ptWeight*ptReweight;
// */
				weight = 1;
			}
			else {
				if (oniaMode == 1){ ptReweight = (f1Spp->Eval(ptReco)); }
				if (oniaMode == 2){ ptReweight = (f2Spp->Eval(ptReco)); }
				if (oniaMode == 3){ ptReweight = 1;}
				//ptReweight = 1;
				weight = ptReweight;
			}

			bool recoPass = 0;

			if (Reco_QQ_sign[iQQ] == 0 && acceptMu && mupl_cut && mumi_cut && trigL1Dmu){ recoPass = 1; }


			//filling RecoEvent Histo if passing
			if (rapLow < rapReco < rapHigh && ptReco < 30 && Centrality < 200&&TMath::Abs(qq4mom->Eta())<2.5){
				if (recoPass == 1 && PtCutPass == 1 && MassCutPass == 1){
					RecoEvents->Fill(Centrality/2., weight);
					hRecoEventsD->Fill(Centrality/2., weight);
					RecoEventsNtracks->Fill(Ntracks, weight*sumET_HFWeight);
					RecoEventsSumET_HF->Fill(SumET_HF, weight*ntracksWeight);
					RecoEventsInt->Fill(Centrality/2., weight);
					RecoEventsPt->Fill(ptReco, weight*sumET_HFWeight*ntracksWeight);
					RecoEventsRap->Fill(rapReco, weight*sumET_HFWeight*ntracksWeight);
					hCentrality->Fill(Centrality, weight);
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
			if (IsAccept(g_mupl4mom) && IsAccept(g_mumi4mom)){ acceptMu = 1; }
			if (PtCut(g_mupl4mom) && PtCut(g_mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(g_qq4mom, massLow, massHigh);



			//weights only needed for pPb
			float weight = 0;
			ptWeight = 0;
			centWeight = 0;
			ntracksWeight=0;
			sumET_HFWeight= 0;
			centWeight = FindCenWeight(Centrality);
			ntracksWeight= FindNtracksWeight(Ntracks,Ntracks_Weights);
			sumET_HFWeight= FindSumET_HFWeight(SumET_HF,SumET_HF_Weights);
			ptReweight = 0;

			//getting a pt gen value 
			float ptGen = 0;
			float rapGen = 0;
			ptGen = g_qq4mom->Pt();
//			rapGen = TMath::Abs(g_qq4mom->Rapidity());
			rapGen = g_qq4mom->Rapidity();


			int tNum = -1;
			//getting the tree pt mc weighting from generation
			if (ispPb){
				//if (oniaMode == 1){ ptReweight = (f1SAA->Eval(ptGen)); }
				//if (oniaMode == 2){ ptReweight = (f2SAA->Eval(ptGen)); }
                                //if (oniaMode == 3){ ptReweight = 1;}
				//tNum = myTree.GetTreeNumber();
				//ptReweight = 1;
				//ptWeight = GetWeight(tNum, oniaMode);
				//weight = centWeight*ptWeight*ptReweight;
				weight = 1;
			}
			else{
				if (oniaMode == 1){ ptReweight = (f1Spp->Eval(ptGen)); }
				if (oniaMode == 2){ ptReweight = (f2Spp->Eval(ptGen)); }
                                if (oniaMode == 3){ ptReweight = 1;}
				//ptReweight = 1;
				weight = ptReweight;
			}

			//fill GenEvent Histo Denominator if passing 
			if (rapLow < rapGen < rapHigh && ptGen < 30 && Centrality < 200 &&TMath::Abs(g_qq4mom->Eta())<2.5){
				if (acceptMu == 1 && PtCutPass == 1 && MassCutPass == 1){
					GenEvents->Fill(Centrality/2., weight);
					GenEventsNtracks->Fill(Ntracks, weight*sumET_HFWeight);
					GenEventsSumET_HF->Fill(SumET_HF, weight*ntracksWeight);
					hGenEventsD->Fill(Centrality/2., weight);
					GenEventsInt->Fill(Centrality/2., weight);
					GenEventsPt->Fill(ptGen, weight*ntracksWeight*sumET_HFWeight);
					GenEventsRap->Fill(rapGen, weight*ntracksWeight*sumET_HFWeight);

				}
			}

		}


	}


if(!ispPb){
iPeriod = 6;
}
//------Cent---------       
//dividing the RecoEvents by GenEvents 
TGraphAsymmErrors *EffCent = new TGraphAsymmErrors(nCenBin);
EffCent->BayesDivide(RecoEvents, GenEvents);
EffCent->SetName("EffCent");

if(ispPb){
TCanvas *c1 = new TCanvas("c1","c1",800,600);
c1->SetRightMargin(1);
c1->cd();
EffCent->SetMarkerSize(2.0);
EffCent->SetMarkerColor(kRed);
EffCent->SetMarkerStyle(20);

EffCent->SetTitle("");
EffCent->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
EffCent->GetXaxis()->SetTitle(Form("%s",ispPb ? "Centrality" : "Integrated Bin"));
EffCent->GetYaxis()->SetRangeUser(0,1);
EffCent->GetXaxis()->SetRangeUser(0.0, 100.0);
EffCent->GetXaxis()->CenterTitle();
EffCent->GetYaxis()->CenterTitle();
EffCent->GetXaxis()->SetTitleOffset(1);
EffCent->GetYaxis()->SetTitleOffset(1);

EffCent->Draw("AP");
CMS_lumi(c1,iPeriod, iPos);
c1->Update();

c1->SaveAs(Form("eff_pPb/EfficiencyCent_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));
}

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

c2->SaveAs(Form("eff_pPb/EfficiencyPt_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));

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
EffRap->GetXaxis()->SetTitle("y");
EffRap->GetYaxis()->SetRangeUser(0,1);
EffRap->GetXaxis()->SetRangeUser(rapLow,rapHigh);
EffRap->GetXaxis()->CenterTitle();
EffRap->GetYaxis()->CenterTitle();
EffRap->GetXaxis()->SetTitleOffset(1);
EffRap->GetYaxis()->SetTitleOffset(1);

EffRap->Draw("AP");
CMS_lumi(c3,iPeriod, iPos);
c3->Update();

c3->SaveAs(Form("eff_pPb/EfficiencyRap_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));

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

c4->SaveAs(Form("eff_pPb/EfficiencyInt_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));

//------Ntracks---------       
//dividing the RecoEvents by GenEvents 
TGraphAsymmErrors *EffNtracks = new TGraphAsymmErrors(nNtracksBin);
EffNtracks->BayesDivide(RecoEventsNtracks, GenEventsNtracks);
EffNtracks->SetName("EffNtracks");

if(ispPb){
TCanvas *c5 = new TCanvas("c5","c5",800,600);
c5->SetRightMargin(1);
c5->cd();
EffNtracks->SetMarkerSize(2.0);
EffNtracks->SetMarkerColor(kRed);
EffNtracks->SetMarkerStyle(20);

EffNtracks->SetTitle("");
EffNtracks->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
EffNtracks->GetXaxis()->SetTitle(Form("%s",ispPb ? "Ntracks" : "Integrated Bin"));
EffNtracks->GetYaxis()->SetRangeUser(0,1);
EffNtracks->GetXaxis()->SetRangeUser(0.0, 200.0);
EffNtracks->GetXaxis()->CenterTitle();
EffNtracks->GetYaxis()->CenterTitle();
EffNtracks->GetXaxis()->SetTitleOffset(1);
EffNtracks->GetYaxis()->SetTitleOffset(1);

EffNtracks->Draw("AP");
CMS_lumi(c5,iPeriod, iPos);
c5->Update();

c5->SaveAs(Form("eff_pPb/EfficiencyNtracks_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));
}

//------SumET_HF---------       
//dividing the RecoEvents by GenEvents 
TGraphAsymmErrors *EffSumET_HF = new TGraphAsymmErrors(nSumET_HFBin);
EffSumET_HF->BayesDivide(RecoEventsSumET_HF, GenEventsSumET_HF);
EffSumET_HF->SetName("EffSumET_HF");

if(ispPb){
TCanvas *c6 = new TCanvas("c6","c6",800,600);
c6->SetRightMargin(1);
c6->cd();
EffSumET_HF->SetMarkerSize(2.0);
EffSumET_HF->SetMarkerColor(kRed);
EffSumET_HF->SetMarkerStyle(20);

EffSumET_HF->SetTitle("");
EffSumET_HF->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(%dS)]_{%s}",oniaMode, ispPb ? "pPb" : "PP"));
EffSumET_HF->GetXaxis()->SetTitle(Form("%s",ispPb ? "SumET_HF" : "Integrated Bin"));
EffSumET_HF->GetYaxis()->SetRangeUser(0,1);
EffSumET_HF->GetXaxis()->SetRangeUser(0.0, 140.0);
EffSumET_HF->GetXaxis()->CenterTitle();
EffSumET_HF->GetYaxis()->CenterTitle();
EffSumET_HF->GetXaxis()->SetTitleOffset(1);
EffSumET_HF->GetYaxis()->SetTitleOffset(1);

EffSumET_HF->Draw("AP");
CMS_lumi(c6,iPeriod, iPos);
c6->Update();

c6->SaveAs(Form("eff_pPb/EfficiencySumET_HF_%dS_%s_Ntrackseta_gt25.png",oniaMode,ispPb ? "pPb" : "PP"));
}

TFile* MyFileEff;
MyFileEff = new TFile(Form("eff_pPb/Eff_%s_%dS_Ntrackseta_gt25.root",ispPb ? "pPb" : "PP",oniaMode), "Recreate");
if(ispPb){
	EffSumET_HF->Write();
	EffNtracks->Write();
	EffCent->Write();
	GenEvents->Write();
	RecoEvents->Write();
}

hGenEventsD->Write();
hRecoEventsD->Write();
RecoEventsInt->Write();
RecoEventsPt->Write();
RecoEventsRap->Write();
GenEventsInt->Write();
GenEventsPt->Write();
GenEventsRap->Write();
hCentrality->Write();
hCrossCheck->Write();
EffPt->Write();
EffRap->Write();
EffInt->Write();
MyFileEff->Close();

// Writing out efficiencies
        for (Int_t i = 0; i < (nPtBin); i++){
        cout << "Pt" << EffPt->Eval(ptBin[i]) << " , - " << EffPt->GetErrorYlow(i) << " , + " << EffPt->GetErrorYhigh(i) << endl;
	}
        for (Int_t i = 0; i < (nRapBin); i++){
        cout << "Rapidity" << EffRap->Eval(rapBin[i]) << " , - " << EffRap->GetErrorYlow(i) << " , + " << EffRap->GetErrorYhigh(i) << endl;
        }
        for (Int_t i = 0; i < (nCenBin); i++){
        cout << "Centrality" << EffCent->Eval(CenBin[i]) << " , - " << EffCent->GetErrorYlow(i) << " , + " << EffCent->GetErrorYhigh(i) << endl;
        }
        for (Int_t i = 0; i < (nNtracksBin); i++){
        cout << "Ntracks" << EffNtracks->Eval(NtracksBin[i]) << " , - " << EffNtracks->GetErrorYlow(i) << " , + " << EffNtracks->GetErrorYhigh(i) << endl;
        }
        for (Int_t i = 0; i < (nSumET_HFBin); i++){
        cout << "SumET_HF" << EffSumET_HF->Eval(SumET_HFBin[i]) << " , - " << EffSumET_HF->GetErrorYlow(i) << " , + " << EffSumET_HF->GetErrorYhigh(i) << endl;
        }

        ReweightFunctions->Close();

}



//Returns a boolean for muon in acceptance
bool IsAccept(TLorentzVector *Muon){
	return (
			(( fabs(Muon->Eta())>=0.0 && fabs(Muon->Eta())<1.0 ) && Muon->Pt()>3.4) ||
			(( fabs(Muon->Eta())>=1.0 && fabs(Muon->Eta())<1.5 ) && Muon->Pt()>(5.8-2.4*fabs(Muon->Eta())) ) ||
			(( fabs(Muon->Eta())>=1.5 && fabs(Muon->Eta())<2.4 ) && Muon->Pt()>(3.4-0.78*fabs(Muon->Eta())) )
	       );
}


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

double PtReweight(TLorentzVector* DiMuon, double coefficient, double constant){
        double f = coefficient*(DiMuon->Pt()) + constant;
        return f;
}


double GetWeight(int numTree,int oniaMode){
  double weight1[6] = {3.10497,4.11498,2.2579,1.2591,0.567094,0.783399};
  double weight2[6] = {5.89168,9.08207,3.106,1.10018,0.534916,0.776183};
  double weight3[4] = {6.86815,8.29618,6.75153,5.48684};
  if(oniaMode ==1)return weight1[numTree];
  if(oniaMode ==2) return weight2[numTree];
  else return weight3[numTree];
}

//this re weights centrality dist. 
double FindCenWeight(int Bin) {
	const int nbins = 200;
	const double Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
	return Ncoll[Bin];
}
double FindNtracksWeight(int Ntracks, TH1D *Ntracks_Weights) {
	int nbin = Ntracks_Weights->GetXaxis()->FindBin(Ntracks);
	//cout<<"Ntracks: "<<Ntracks<<endl;
	//cout<<"Bin Number: "<<nbin<<endl;
	//cout<<Ntracks_Weights->GetBinContent(nbin)<<endl;
	return Ntracks_Weights->GetBinContent(nbin);
}
double FindSumET_HFWeight(int SumET_HF, TH1D *SumET_HF_Weights) {
	int nbin = SumET_HF_Weights->GetXaxis()->FindBin(SumET_HF);
	return SumET_HF_Weights->GetBinContent(nbin);
}
/*double FindNtracksWeight(int Bin) {
	const int nbins = 200;
	const double Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
	return Ncoll[Bin];
}*/
