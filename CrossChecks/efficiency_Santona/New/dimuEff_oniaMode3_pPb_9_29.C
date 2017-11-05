#include "effCommon.h"
#include "tnp_weight.h"
#include <algorithm>
const double muonPtCut = 4.0;
TFile* fTnp_pa_new = new TFile("output_official_5eta_cutG_all_nominal_v3.root","READ");
TF1* hTnp_pa_new_eta1 = (TF1*)fTnp_pa_new->Get("func_1");
TF1* hTnp_pa_new_eta2 = (TF1*)fTnp_pa_new->Get("func_2");	
		TF1* hTnp_pa_new_eta3 = (TF1*)fTnp_pa_new->Get("func_3");
		TF1* hTnp_pa_new_eta4 = (TF1*)fTnp_pa_new->Get("func_4");
		TF1* hTnp_pa_new_eta5 = (TF1*)fTnp_pa_new->Get("func_5");

//Returns a boolean for muon in acceptance for pPb!  Accounting for weird cut off at eta = 2.0
bool IsAccept(TLorentzVector *Muon){
	return (
			(( fabs(Muon->Eta())>=0.0 && fabs(Muon->Eta())<1.2 ) && Muon->Pt()>3.4) ||
			(( fabs(Muon->Eta())>=1.2 && fabs(Muon->Eta())<1.5 ) && Muon->Pt()>(10.6-6*fabs(Muon->Eta())) ) ||
			(( fabs(Muon->Eta())>=1.5 && fabs(Muon->Eta())<2.0 ) && Muon->Pt()>(3.4-1.2*fabs(Muon->Eta())) ) ||
			(( Muon->Eta()>-2.4 && Muon->Eta()<=-2.0 ) && Muon->Pt()>1.0 )
	       );
}

// For rapidity shifted pp!
bool IsAcceptpp(TLorentzVector *Muon){
	double BoostedEta = Muon->Eta() -0.47;
        return (
                        (( fabs(BoostedEta)>=0.0 && fabs(BoostedEta)<1.2 ) && Muon->Pt()>3.4) ||
                        (( fabs(BoostedEta)>=1.2 && fabs(BoostedEta)<1.5 ) && Muon->Pt()>(10.6-6*fabs(BoostedEta)) ) ||
                        (( fabs(BoostedEta)>=1.5 && fabs(BoostedEta)<2.0 ) && Muon->Pt()>(3.4-1.2*fabs(BoostedEta)) ) ||
                        (( BoostedEta>-2.4 && BoostedEta<=-2.0 ) && Muon->Pt()>1.0 )
                        //(( fabs(BoostedEta)>=0.0 && fabs(BoostedEta)<1.0 ) && Muon->Pt()>3.4) ||
                        //(( fabs(BoostedEta)>=1.0 && fabs(BoostedEta)<1.5 ) && Muon->Pt()>(5.8-2.4*fabs(BoostedEta)) ) ||
                        //(( fabs(BoostedEta)>=1.5 && fabs(BoostedEta)<2.4 ) && Muon->Pt()>(3.4-0.78*fabs(BoostedEta)) )
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

double PtReweight(TLorentzVector* DiMuon, TF1 *Pt_ReWeights){
        double pT = (DiMuon->Pt());
        return Pt_ReWeights->Eval(pT);
}

double weight_tp_pp(double pt, double eta, bool ispPb, int idx_variation)
{
   if (!ispPb)
   {
      return tnp_weight_muidtrg_pp(pt, eta, idx_variation);
      // if (fabs(eta)<1.6)
      //    return tnp_weight_pbpb_midrap(pt, idx_variation);
      // else
      //    return tnp_weight_pbpb_fwdrap(pt, idx_variation);
   }
}

double weight_tp_pPb(bool ispPb, int idx_variation, double mupt1,double mupt2,double mueta1, double mueta2,TFile* fTnp_pa_new)
{
   if (ispPb)
   {
      TF1* hw1;
		TF1* hw2;
		 // 20160626 (off8M, tagpt5)
		if (  TMath::Abs(mueta1) < 0.9 )      hw1 = hTnp_pa_new_eta1;
		else if ( TMath::Abs(mueta1) < 1.2 )  hw1 = hTnp_pa_new_eta2;
		else if ( TMath::Abs(mueta1) < 1.6 )  hw1 = hTnp_pa_new_eta3;
		else if ( TMath::Abs(mueta1) < 2.1 )  hw1 = hTnp_pa_new_eta4;
		else                                  hw1 = hTnp_pa_new_eta5;
		if (  TMath::Abs(mueta2) < 0.9 )      hw2 = hTnp_pa_new_eta1;
		else if ( TMath::Abs(mueta2) < 1.2 )  hw2 = hTnp_pa_new_eta2;
		else if ( TMath::Abs(mueta2) < 1.6 )  hw2 = hTnp_pa_new_eta3;
		else if ( TMath::Abs(mueta2) < 2.1 )  hw2 = hTnp_pa_new_eta4;
		else                                  hw2 = hTnp_pa_new_eta5;

		double tnpWeight1 = hw1->Eval(mupt1);
		double tnpWeight2 = hw2->Eval(mupt2);
		//cout<<tnpWeight1*tnpWeight2<<endl;
		return tnpWeight1 * tnpWeight2;
		//fTnp_pa_new->Close();
   }
}

double weight_tpsta(double pt, double eta, bool ispPb, int idx_variation){
   if (ispPb)
   {
      return 1;
   }
   else
   {
      return tnp_weight_sta_pp(pt, eta, idx_variation);
   }
}

int  nPtBin;     
int  nRapBin;     

double m1S_low = 8.0;
double m1S_high = 10.0;
double m2S_low = 8.563;
double m2S_high = 10.563;
double m3S_low = 8.895;
double m3S_high = 10.895;

int iPeriod = 5;
int iPos = 33;


void dimuEff_oniaMode3_pPb_9_29(
	int oniaMode = 3, //1 = 1S, 2 = 2S, 3 = 3S
	bool ispPb = 1 //true = pPb and false = pp
	){   
	int var_tp1 = 0;
	int var_tp2 = 0;

	setTDRStyle();

	TChain *myTree_pp = new TChain("hionia/myTree");
	// For using pp gen tree as denominator for all states.   
        if ((oniaMode == 1) && ispPb){
		myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups1SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
		}

        if ((oniaMode == 2) && ispPb){
        myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }

        if ((oniaMode == 3) && ispPb){
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

const int nYBins1S  = 8;   //double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 4;   //double yBin2S[nYBins2S+1] ={0, 1.2, 2.4};
const int nYBins3S  = 2;   //double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};

std::vector<double> ptBinEdges;
std::vector<double> ptBin;
std::vector<double> ptBinErr;
std::vector<double> rapBinEdges;
std::vector<double> rapBin;
//std::vector<double> rapBinErr;

//declare the number of bins and assign bin edges
if(oniaMode ==1){
	nPtBin = nPtBins1s;
	nRapBin = nYBins1S;

	ptBinEdges = {0,2,4,6,9,12,30};
	ptBin = {1,3,5,7.5,10.5,21};
	ptBinErr = {1,1,1,1.5,1.5,9};
        rapBinEdges = {-2.4, -1.67, -1.27, -0.87, -0.47, -0.07, 0.33, 0.73, 1.46};
        rapBin = {-2.035, -1.47, -1.07, -0.67, -0.27, 0.13, 0.53, 1.095};
}
if(oniaMode ==2){
	nPtBin = nPtBins2s;
	nRapBin = nYBins2S;

	ptBinEdges = {0,4,9,30};
	ptBin = {2,6.5,19.5};
	ptBinErr = {2,2.5,10.5};
        rapBinEdges = {-2.4, -1.27, -0.47, 0.33, 1.46};
        rapBin = { -1.835, -0.6585, -0.07, 0.895 };
}
if(oniaMode ==3){
	nPtBin = nPtBins3s;
	nRapBin = nYBins3S;

	ptBinEdges = {0.0,6.0,30.0};
	ptBin = {3.0,18.0};
	ptBinErr = {3.0,12.0};
        rapBinEdges = {-2.4, -0.47, 1.46};
        rapBin = { -1.435, 0.495 };
}


        float rapLow = -2.4;
        float rapHigh = 1.46;

	double weighttp;
	double weighttpsta;

//	float           IntBin[1] = { 100 };
//	float		IntBinEdges[2] = { 0, 100 };
	float           ptReweight = 0.0;

	float 		massLow = 0;
	float 		massHigh = 0;

//	Int_t           Centrality;
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
// */

//	TBranch        *b_Centrality;   //!
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
//	Reco_QQ_size = 0;
	Reco_QQ_4mom = 0;
	Reco_QQ_mupl_4mom = 0;
	Reco_QQ_mumi_4mom = 0;

//	Gen_QQ_size = 0;
	Gen_QQ_4mom = 0;
	Gen_QQ_mupl_4mom = 0;
	Gen_QQ_mumi_4mom = 0;

//	myTree_pPb->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
	myTree_pPb->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	myTree_pPb->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	myTree_pPb->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	myTree_pPb->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
	myTree_pPb->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	myTree_pPb->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
//	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
//	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
	myTree_pPb->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
	myTree_pPb->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);

	myTree_pp->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
	myTree_pp->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
	myTree_pp->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
	myTree_pp->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

	myTree_pPb->SetBranchStatus("*", 0);
        myTree_pp->SetBranchStatus("*", 0);

//	myTree_pPb->SetBranchStatus("Centrality", 1);
	myTree_pPb->SetBranchStatus("HLTriggers", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_size", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_sign", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_4mom", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_4mom", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_4mom", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_trig", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_VtxProb", 1);
//	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_isGoodMuon", 1);
//	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_isGoodMuon", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_nPixWMea", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_nPixWMea", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_nTrkWMea", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_nTrkWMea", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_dxy", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_dxy", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mupl_dz", 1);
	myTree_pPb->SetBranchStatus("Reco_QQ_mumi_dz", 1);

	myTree_pp->SetBranchStatus("Gen_QQ_size", 1);
	myTree_pp->SetBranchStatus("Gen_QQ_4mom", 1);
	myTree_pp->SetBranchStatus("Gen_QQ_mupl_4mom", 1);
	myTree_pp->SetBranchStatus("Gen_QQ_mumi_4mom", 1);

	//convert bin and bin edge vectors to arrays to be used as parameter when delcaring TH1Ds
	double* ptBinEdges_arr = &ptBinEdges[0];
	double* ptBin_arr = &ptBin[0];
	double* rapBinEdges_arr = &rapBinEdges[0];
	double* rapBin_arr = &rapBin[0];
	double* ptBinErr_arr = &ptBinErr[0];
/*	double* CenBinEdges_arr = &CenBinEdges[0];
	double* CenBin_arr = &CenBin[0];
// */
//	TH1D  *RecoEvents = new TH1D("RecoEvents", "Reconstructed", nCenBin, CenBinEdges_arr);
//	TH1D  *GenEvents = new TH1D("GenEvents", "Generated", nCenBin, CenBinEdges_arr);
//	TH1D  *RecoEventsInt = new TH1D("RecoEventsInt", "Reconstructed", 1, IntBinEdges);
//	TH1D  *GenEventsInt = new TH1D("GenEventsInt", "Generated", 1, IntBinEdges);
	
	TH1D  *RecoEventsPt = new TH1D("RecoEventsPt", "Reconstructed", nPtBin, ptBinEdges_arr);
	TH1D  *GenEventsPt = new TH1D("GenEventsPt", "Generated", nPtBin, ptBinEdges_arr);

	TH1D  *RecoEventsRap = new TH1D("RecoEventsRap", "Reconstructed", nRapBin, rapBinEdges_arr);
	TH1D  *GenEventsRap = new TH1D("GenEventsRap", "Generated", nRapBin, rapBinEdges_arr);

//	TH1D  *hCentrality = new TH1D("hCentrality", "Centrality distribution", 202, -1, 201);
	TH1D  *hCrossCheck = new TH1D("hCrossCheck", "Checking number of events", 2, 0, 2);

/*	TH1D  *hRecoEventsD = new TH1D("hRecoEventsD", "Reconstructed", nCenBin, CenBinEdges_arr);
	TH1D  *hGenEventsD = new TH1D("hGenEventsD", "Generated", nCenBin, CenBinEdges_arr);
// */
	RecoEventsPt->Sumw2();
	GenEventsPt->Sumw2();
	RecoEventsRap->Sumw2();
	GenEventsRap->Sumw2();
	hCrossCheck->Sumw2();

	std::string fmode="1";
	const char *f_name;
        const char *f_name_pp;

	//TF1* Pt_Weights = (TF1*)PtReweightFunctions->Get("dataMC_Ratio_norm");
/*	if(!ispPb){
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
// */

	f_name_pp = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PP_DATA_3s_20170816.root";
	f_name = "../../../CompareDataToMC/WeightedFcN_fit/ratioDataMC_PA_DATA_3s_20170816.root";

	TFile* PtReweightFunctions = new TFile(f_name, "Open");
        TFile* PtReweightFunctions_pp = new TFile(f_name_pp, "Open");

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

	TF1* Pt_ReWeights = (TF1*)PtReweightFunctions->Get("dataMC_Ratio_norm");
        TF1* Pt_ReWeights_pp = (TF1*)PtReweightFunctions_pp->Get("dataMC_Ratio_norm");

	Long64_t nentries_pPb = myTree_pPb->GetEntries();
	cout << nentries_pPb << endl;

	Long64_t nentries_pp = myTree_pp->GetEntries();
	cout << nentries_pp << endl;


	for (Long64_t jentry = 0; jentry < nentries_pPb; jentry++){
		myTree_pPb->GetEntry(jentry);
		if(jentry%100000 == 0){
			cout<<"--Processing Event: "<<jentry<<endl;
		}


		//Numerator Loop RECO
		for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++){

			//cout << " Reco QQ size = " << Reco_QQ_size << endl;
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

// */			//--Muon id cuts
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

			//tnp weights only needed for pPb
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

			//getting the tree weight from pt generated MC bins

				ptReweight = PtReweight(qq4mom, Pt_ReWeights);
				//cout<<ptReweight<<endl;
				weight = ptReweight;
				if(!ispPb){
					weighttp=weight_tp_pp(mupl4mom->Pt(),mupl4mom->Eta(),ispPb,var_tp1)*weight_tp_pp(mumi4mom->Pt(),mumi4mom->Eta(),ispPb,var_tp1);
				}else{
					weighttp = weight_tp_pPb(ispPb, var_tp1,mupl4mom->Pt(),mumi4mom->Pt(),mupl4mom->Eta(), mumi4mom->Eta(),fTnp_pa_new);
				}
         		weighttpsta=weight_tpsta(mupl4mom->Pt(),mupl4mom->Eta(),ispPb,var_tp2)*weight_tpsta(mumi4mom->Pt(),mumi4mom->Eta(),ispPb,var_tp2);
         		weighttp *= weighttpsta;
			bool recoPass = 0;

			if (Reco_QQ_sign[iQQ] == 0 && acceptMu && mupl_cut && mumi_cut && trigL1Dmu){ recoPass = 1; }

			//filling RecoEvent Histo if passing
			if (rapLow < rapReco < rapHigh && ptReco < 30){
				if (recoPass == 1 && PtCutPass == 1 && MassCutPass == 1){
					//RecoEventsInt->Fill(Centrality/2., weight*weighttp);
					RecoEventsPt->Fill(ptReco, weight*weighttp);
					RecoEventsRap->Fill(rapReco, weight*weighttp);
					//hCentrality->Fill(Centrality, weight*weighttp);
				}
			}
		}
	}

std::vector<Long64_t> range;
for (Long64_t i = 0; i < nentries_pPb; i++){
	if (i % 100000 == 0){ cout <<i<<endl;}
	range.push_back(i);
}
std::random_shuffle(range.begin(), range.end());
	
        for (Long64_t jentry2 = 0; jentry2 < nentries_pp; jentry2++){ // nentries_pPb for using same number of entries
                Long64_t jentry = range[jentry2];
		myTree_pp->GetEntry(jentry);
                if(jentry%100000 == 0){
                        cout<<"--Processing Event: "<<jentry<<endl;
               }
// */
		//cout <<  " Num loop done " << endl;
		//Denominator loop  GEN
		for (int iGQQ = 0; iGQQ < Gen_QQ_size; iGQQ++){

			//cout << " Gen QQ size = " << Gen_QQ_size << endl;
			hCrossCheck->Fill(0);
			TLorentzVector *g_qq4mom = (TLorentzVector*)Gen_QQ_4mom->At(iGQQ);
			TLorentzVector *g_mumi4mom = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGQQ);
			TLorentzVector *g_mupl4mom = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGQQ);

			bool acceptMu = 0;
			bool PtCutPass = 0;
			bool MassCutPass = 0;

			//check if muons are in acceptance
			if (IsAcceptpp(g_mupl4mom) && IsAcceptpp(g_mumi4mom)){ acceptMu = 1; }
			if (PtCut(g_mupl4mom) && PtCut(g_mumi4mom)){ PtCutPass = 1; }
			MassCutPass = MassCut(g_qq4mom, massLow, massHigh);

			float weight = 0;
			ptReweight = 0;

			float ptGen = 0;
			float rapGen = 0;
			ptGen = g_qq4mom->Pt();

//			rapGen = TMath::Abs(g_qq4mom->Rapidity());
			rapGen = g_qq4mom->Rapidity() - 0.47;
			//cout << "rapGen = " << rapGen << endl;

				ptReweight = PtReweight(g_qq4mom, Pt_ReWeights_pp);
				//cout<<ptReweight<<endl;
				weight = ptReweight;

			//fill GenEvent Histo Denominator if passing 
			if (rapLow < rapGen < rapHigh && ptGen < 30 ){
				if (acceptMu == 1 && PtCutPass == 1 && MassCutPass == 1){
					//GenEventsInt->Fill(Centrality/2., weight);
					GenEventsPt->Fill(ptGen, weight);
					GenEventsRap->Fill(rapGen, weight);
				}
			}
		}
	//cout << " Den loop done " << endl;
	}


// Plotting

// Checking the histograms that are being used to get the efficiencies
        TGraphAsymmErrors* TGNumPt = new TGraphAsymmErrors(nPtBin);
        TGraphAsymmErrors* TGDenPt = new TGraphAsymmErrors(nPtBin);

	double Num = 0;
	double Den = 0;
	double DenErrL = 0;
	double DenErrH = 0;
	double NumErrL = 0;
	double NumErrH = 0;
	double RatioPt = 0;

        for (Int_t i = 1; i < nPtBin+1; i++){
                Num = RecoEventsPt->GetBinContent(i);
                Den = GenEventsPt->GetBinContent(i);
                NumErrH = RecoEventsPt->GetBinErrorUp(i);
                NumErrL = RecoEventsPt->GetBinErrorLow(i);
                DenErrH = GenEventsPt->GetBinErrorUp(i);
                DenErrL = GenEventsPt->GetBinErrorLow(i);
                RatioPt = Num / Den;

                TGNumPt->SetPoint(i-1, ptBin[i-1], Num);
                TGNumPt->SetPointError(i-1, ptBinErr[i-1], ptBinErr[i-1], NumErrL, NumErrH);
                TGDenPt->SetPoint(i-1, ptBin[i-1], Den);
                TGDenPt->SetPointError(i-1, ptBinErr[i-1], ptBinErr[i-1], DenErrL, DenErrH);
	}

TCanvas *cCheck = new TCanvas("cCheck","cCheck",800,600);
cCheck->SetRightMargin(1);
cCheck->cd();

TGNumPt->SetName("TGNumPt");
TGNumPt->SetMarkerSize(1.0);
TGNumPt->SetMarkerColor(kRed);
TGNumPt->SetMarkerStyle(20);

TGNumPt->SetTitle("");
TGNumPt->GetYaxis()->SetTitle(Form("Number of events in %dS in %s",oniaMode, ispPb ? "pPb" : "PP"));
TGNumPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
TGNumPt->GetYaxis()->SetRangeUser(0,600000);
TGNumPt->GetXaxis()->SetRangeUser(0.0, 30.0);
TGNumPt->GetXaxis()->CenterTitle();
TGNumPt->GetYaxis()->CenterTitle();
TGNumPt->GetXaxis()->SetTitleOffset(1);
TGNumPt->GetYaxis()->SetTitleOffset(1);
TGNumPt->Draw("AP");

TGDenPt->SetName("TGDenPt");
TGDenPt->SetMarkerSize(1.0);
TGDenPt->SetMarkerColor(kBlue);
TGDenPt->SetMarkerStyle(20);
TGDenPt->Draw("Psame");

CMS_lumi(cCheck,iPeriod, iPos);
cCheck->Update();

cCheck->SaveAs(Form("eff_pPb9_29/NumDenPt_%dS_%s_9_29.png",oniaMode, ispPb ? "pPb" : "PP"));



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

c2->SaveAs(Form("eff_pPb9_29/EfficiencyPt_%dS_%s_9_29.png",oniaMode, ispPb ? "pPb" : "PP"));

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

c3->SaveAs(Form("eff_pPb9_29/EfficiencyRap_%dS_%s_9_29.png",oniaMode,ispPb ? "pPb" : "PP"));

/*
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

c4->SaveAs(Form("eff_pPb9_29/EfficiencyInt_%dS_%s_9_29.png",oniaMode, ispPb ? "pPb" : "PP"));
// */


TFile* MyFileEff;
MyFileEff = new TFile(Form("eff_pPb9_29/Eff_%s_%dS_9_29.root","pPb",oniaMode), "Recreate");
//hGenEventsD->Write();
//hRecoEventsD->Write();
//RecoEventsInt->Write();
RecoEventsPt->Write();
RecoEventsRap->Write();
//GenEventsInt->Write();
GenEventsPt->Write();
GenEventsRap->Write();
hCrossCheck->Write();
EffPt->Write();
EffRap->Write();
//EffInt->Write();
MyFileEff->Close();

// Writing out efficiencies
        for (Int_t i = 0; i < (nPtBin); i++){
        cout << "Pt" << EffPt->Eval(ptBin_arr[i]) << " , - " << EffPt->GetErrorYlow(i) << " , + " << EffPt->GetErrorYhigh(i) << endl;
	}
        for (Int_t i = 0; i < (nRapBin); i++){
        cout << "Rapidity" << EffRap->Eval(rapBin_arr[i]) << " , - " << EffRap->GetErrorYlow(i) << " , + " << EffRap->GetErrorYhigh(i) << endl;
        }

}    // end void



