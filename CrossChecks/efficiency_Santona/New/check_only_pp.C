#include "effCommon.h"
#include "tnp_weight.h"

void check_only_pp(
        int oniaMode = 3, //1 = 1S, 2 = 2S, 3 = 3S
        bool ispPb = 0 //true = pPb and false = pp
	){

	TChain *myTree_pp = new TChain("hionia/myTree");
        if ((oniaMode == 2) && !ispPb){
        myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups2SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }

        if ((oniaMode == 3) && !ispPb){
        myTree_pp->Add("/scratch_menkar/CMS_Trees/OniaTrees_2015_5TeV/pp_MC_Official/OniaTree_Ups3SMM_5p02TeV_TuneCUETP8M1_HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1.root");
        }

        Int_t           Gen_QQ_size;
//        Int_t           Gen_QQ_sign[45];   //[Gen_QQ_size]
        TClonesArray    *Gen_QQ_4mom;
        TClonesArray    *Gen_QQ_mupl_4mom;
        TClonesArray    *Gen_QQ_mumi_4mom;

        TBranch        *b_Gen_QQ_size;   //
//	TBranch        *b_Gen_QQ_sign; 
        TBranch        *b_Gen_QQ_4mom;   //!
        TBranch        *b_Gen_QQ_mupl_4mom;   //!
        TBranch        *b_Gen_QQ_mumi_4mom;   //!

//        Gen_QQ_size = 0;
        Gen_QQ_4mom = 0;
        Gen_QQ_mupl_4mom = 0;
        Gen_QQ_mumi_4mom = 0;

/*
        myTree_pPb->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
        myTree_pPb->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
        myTree_pPb->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
        myTree_pPb->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
        myTree_pPb->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
// */
        myTree_pp->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
        myTree_pp->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
        myTree_pp->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
        myTree_pp->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);

//        myTree_pPb->SetBranchStatus("*", 0);
        myTree_pp->SetBranchStatus("*", 0);


/*      myTree_pPb->SetBranchStatus("Reco_QQ_size", 1);
        myTree_pPb->SetBranchStatus("Reco_QQ_sign", 1);
        myTree_pPb->SetBranchStatus("Reco_QQ_4mom", 1);
        myTree_pPb->SetBranchStatus("Reco_QQ_mupl_4mom", 1);
        myTree_pPb->SetBranchStatus("Reco_QQ_mumi_4mom", 1);
// */

        myTree_pp->SetBranchStatus("Gen_QQ_size", 1);
        myTree_pp->SetBranchStatus("Gen_QQ_4mom", 1);
        myTree_pp->SetBranchStatus("Gen_QQ_mupl_4mom", 1);
        myTree_pp->SetBranchStatus("Gen_QQ_mumi_4mom", 1);


        Long64_t nentries_pp = myTree_pp->GetEntries();
        cout << nentries_pp << endl;

        for (Long64_t jentry = 0; jentry < nentries_pp; jentry++){
//                myTree_pPb->GetEntry(jentry);
                if(jentry%100000 == 0){
                        cout<<"--Processing Event: "<<jentry<<endl;
                }


//	cout << " Reco QQ size = " << Reco_QQ_size << endl;
	cout << " Gen QQ size = " << Gen_QQ_size << endl;
	}

}


