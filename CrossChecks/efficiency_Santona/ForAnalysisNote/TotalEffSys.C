#include "effCommon.h"

double RError(double A, double eA, double B, double eB);
double PError(double A, double eA, double B, double eB);

int iPeriod = 5;
int iPos = 33;

void TotalEffSys(int oniaMode = 1){
	setTDRStyle();


const int nPtBins1s  = 6;
const int nPtBins2s  = 3;
const int nPtBins3s  = 2;


        float           SystematicErrorCent[nCenBin];
        float           SystematicErrorPt[nPtBin];
        float           SystematicErrorRap[nRapBin];
        float           SystematicErrorInt[1];

        float           EffCentValues[nCenBin];
        float           xerrorCent[nCenBin];
        float           EffRapValues[nRapBin];
        float           xerrorRap[nRapBin];
        float           EffPtValues[nPtBin];
        float           xerrorPt[nPtBin];
        float           EffIntValues[1];
        float           xerrorInt[1];

	TGraphAsymmErrors* gEffPt_pPb;
	TGraphAsymmErrors* gEffPt_noReweight_pPb;
	TGraphAsymmErrors* gEffPt_tnpUp_pPb;
	TGraphAsymmErrors* gEffPt_tnpDown_pPb;

        TGraphAsymmErrors* gEffRap_pPb;
        TGraphAsymmErrors* gEffRap_noReweight_pPb;
        TGraphAsymmErrors* gEffRap_tnpUp_pPb;
        TGraphAsymmErrors* gEffRap_tnpDown_pPb;

        TGraphAsymmErrors* gEffPtRpA_pPb;
        TGraphAsymmErrors* gEffPtRpA_noReweight_pPb;
        TGraphAsymmErrors* gEffPtRpA_tnpUp_pPb;
        TGraphAsymmErrors* gEffPtRpA_tnpDown_pPb;


        TGraphAsymmErrors* gEffPt_pp;
        TGraphAsymmErrors* gEffPt_noReweight_pp;
        TGraphAsymmErrors* gEffPt_tnpUp_pp;
        TGraphAsymmErrors* gEffPt_tnpDown_pp;

        TGraphAsymmErrors* gEffRap_pp;
        TGraphAsymmErrors* gEffRap_noReweight_pp;
        TGraphAsymmErrors* gEffRap_tnpUp_pp;
        TGraphAsymmErrors* gEffRap_tnpDown_pp;;


	// Open pp files
        TFile* fEffNom_pp = new TFile(Form("./Nominal/eff_pp11_20_NewRpABin/Eff_pp_%dS_11_20_NewRpABin.root",oniaMode), "Open");
	fEffNom_pp->GetObject("EffPt", gEffPt_pp);
        fEffNom_pp->GetObject("EffRap", gEffRap_pp);

        TFile* fEffnoReweight_pp = new TFile(Form("./NoPtReweights/eff_pp11_21_NoPtReweight_NewRpABin/Eff_pp_%dS_11_20_NoPtReweight_NewRpABin.root",oniaMode), "Open");
        fEffnoReweight_pp->GetObject("EffPt", gEffPt_noReweight_pp);
        fEffnoReweight_pp->GetObject("EffRap", gEffRap_noReweight_pp);

        TFile* fEfftnpUp_pp = new TFile(Form("./TnPSysUP/eff_pp11_21_TnPup_NewRpABin/Eff_pp_%dS_11_21_TnPup_NewRpABin.root",oniaMode), "Open");
        fEfftnpUp_pp->GetObject("EffPt", gEffPt_tnpUp_pp);
        fEfftnpUp_pp->GetObject("EffRap", gEffRap_tnpUp_pp);

        TFile* fEfftnpDown_pp = new TFile(Form("./TnPSysDOWN/eff_pp11_21_TnPdown_NewRpABin/Eff_pp_%dS_11_21_TnPdown_NewRpABin.root",oniaMode), "Open");
        fEfftnpDown_pp->GetObject("EffPt", gEffPt_tnpDown_pp);
        fEfftnpDown_pp->GetObject("EffRap", gEffRap_tnpDown_pp);

	// Open pPb files
        TFile* fEffNom_pPb = new TFile(Form("./Nominal/eff_pPb11_20_NewRpABin/Eff_pPb_%dS_11_20_NewRpABin.root",oniaMode), "Open");
        fEffNom_pPb->GetObject("EffPt", gEffPt_pPb);
        fEffNom_pPb->GetObject("EffRap", gEffRap_pPb);
        fEffNom_pPb->GetObject("EffPtRpA", gEffPtRpA_pPb);

        TFile* fEffnoReweight_pPb = new TFile(Form("./NoPtReweights/eff_pPb11_21_NoPtReweight_NewRpABin/Eff_pPb_%dS_11_20_NoPtReweight_NewRpABin.root",oniaMode), "Open");
        fEffnoReweight_pPb->GetObject("EffPt", gEffPt_noReweight_pPb);
        fEffnoReweight_pPb->GetObject("EffRap", gEffRap_noReweight_pPb);
        fEffnoReweight_pPb->GetObject("EffPtRpA", gEffPtRpA_noReweight_pPb);

        TFile* fEfftnpUp_pPb = new TFile(Form("./TnPSysUP/eff_pPb11_21_TnPup_NewRpABin/Eff_pPb_%dS_11_21_TnPup_NewRpABin.root",oniaMode), "Open");
        fEfftnpUp_pPb->GetObject("EffPt", gEffPt_tnpUp_pPb);
        fEfftnpUp_pPb->GetObject("EffRap", gEffRap_tnpUp_pPb);
        fEfftnpUp_pPb->GetObject("EffPtRpA", gEffPtRpA_tnpUp_pPb);

        TFile* fEfftnpDown_pPb = new TFile(Form("./TnPSysDOWN/eff_pPb11_21_TnPdown_NewRpABin/Eff_pPb_%dS_11_21_TnPdown_NewRpABin.root",oniaMode), "Open");
        fEfftnpDown_pPb->GetObject("EffPt", gEffPt_tnpDown_pPb);
        fEfftnpDown_pPb->GetObject("EffRap", gEffRap_tnpDown_pPb);
        fEfftnpDown_pPb->GetObject("EffPtRpA", gEffPtRpA_tnpDown_pPb);


// Nominal
        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
	EffDen = 0;
	EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = gEffPtRpA_pPb->Eval(ptBin[i]);
                EffDen = gEffPt_pp->Eval(ptBin[i]);
                EffNumErrH = gEffPtRpA_pPb->GetErrorYhigh(i);
                EffNumErrL = gEffPtRpA_pPb->GetErrorYlow(i);
                EffDenErrH = gEffPt_pp->GetErrorYhigh(i);
                EffDenErrL = gEffPt_pp->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                gEffRat_Pt_Nom->SetPoint(i, ptBin[i], EffRatio);
                gEffRat_Pt_Nom->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);
        }

// No Pt Reweight
        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
        EffDen = 0;
        EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = gEffPtRpA_noReweight_pPb->Eval(ptBin[i]);
                EffDen = gEffPt_noReweight_pp->Eval(ptBin[i]);
                EffNumErrH = gEffPtRpA_noReweight_pPb->GetErrorYhigh(i);
                EffNumErrL = gEffPtRpA_noReweight_pPb->GetErrorYlow(i);
                EffDenErrH = gEffPt_noReweight_pp->GetErrorYhigh(i);
                EffDenErrL = gEffPt_noReweight_pp->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                gEffRat_Pt_NoReweight->SetPoint(i, ptBin[i], EffRatio);
                gEffRat_Pt_NoReweight->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);
        }

// TnP Up
        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
        EffDen = 0;
        EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = gEffPtRpA_tnpUp_pPb->Eval(ptBin[i]);
                EffDen = gEffPt_tnpUp_pp->Eval(ptBin[i]);
                EffNumErrH = gEffPtRpA_tnpUp_pPb->GetErrorYhigh(i);
                EffNumErrL = gEffPtRpA_tnpUp_pPb->GetErrorYlow(i);
                EffDenErrH = gEffPt_tnpUp_pp->GetErrorYhigh(i);
                EffDenErrL = gEffPt_tnpUp_pp->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                gEffRat_Pt_tnpUp->SetPoint(i, ptBin[i], EffRatio);
                gEffRat_Pt_tnpUp->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);
        }

// TnP Down
        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
        EffDen = 0;
        EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = gEffPtRpA_tnpDown_pPb->Eval(ptBin[i]);
                EffDen = gEffPt_tnpDown_pp->Eval(ptBin[i]);
                EffNumErrH = gEffPtRpA_tnpDown_pPb->GetErrorYhigh(i);
                EffNumErrL = gEffPtRpA_tnpDown_pPb->GetErrorYlow(i);
                EffDenErrH = gEffPt_tnpDown_pp->GetErrorYhigh(i);
                EffDenErrL = gEffPt_tnpDown_pp->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                gEffRat_Pt_tnpDown->SetPoint(i, ptBin[i], EffRatio);
                gEffRat_Pt_tnpDown->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);
        }











        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
	EffDen = 0;
	EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nRapBin; i++){
                EffNum = hEffNumRap->Eval(rapBin[i]);
                EffDen = hEffDenRap->Eval(rapBin[i]);
                EffNumErrH = hEffNumRap->GetErrorYhigh(i);
                EffNumErrL = hEffNumRap->GetErrorYlow(i);
                EffDenErrH = hEffDenRap->GetErrorYhigh(i);
                EffDenErrL = hEffDenRap->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL); 	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffRap2S1S->SetPoint(i, rapBin[i], EffRatio);
                hEffRap2S1S->SetPointError(i, rapBinErr[i], rapBinErr[i], EffRatioErrL, EffRatioErrH);


        }

        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
        EffDen = 0;
        EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nRapBin; i++){
                EffNum = hEffNumRapNoReweight->Eval(rapBin[i]);
                EffDen = hEffDenRapNoReweight->Eval(rapBin[i]);
                EffNumErrH = hEffNumRapNoReweight->GetErrorYhigh(i);
                EffNumErrL = hEffNumRapNoReweight->GetErrorYlow(i);
                EffDenErrH = hEffDenRapNoReweight->GetErrorYhigh(i);
                EffDenErrL = hEffDenRapNoReweight->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                hEffRap2S1SNoReweight->SetPoint(i, rapBin[i], EffRatio);
                hEffRap2S1SNoReweight->SetPointError(i, rapBinErr[i], rapBinErr[i], EffRatioErrL, EffRatioErrH);


        }



        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
	EffDen = 0;
	EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nIntBin; i++){
                EffNum = hEffNumInt->Eval(intBin[i]);
                EffDen = hEffDenInt->Eval(intBin[i]);
                EffNumErrH = hEffNumInt->GetErrorYhigh(i);
                EffNumErrL = hEffNumInt->GetErrorYlow(i);
                EffDenErrH = hEffDenInt->GetErrorYhigh(i);
                EffDenErrL = hEffDenInt->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffInt2S1S->SetPoint(i, intBin[i], EffRatio);
                hEffInt2S1S->SetPointError(i, intBinErr[i], intBinErr[i], EffRatioErrL, EffRatioErrH);


        }

        EffRatio = 0;
        EffNum = 0;
        EffNumErrH = 0;
        EffDenErrH = 0;
        EffNumErrL = 0;
        EffDenErrL = 0;
        EffDen = 0;
        EffRatioErrH = 0;
        EffRatioErrL = 0;

        for (Int_t i = 0; i < nIntBin; i++){
                EffNum = hEffNumIntNoReweight->Eval(intBin[i]);
                EffDen = hEffDenIntNoReweight->Eval(intBin[i]);
                EffNumErrH = hEffNumIntNoReweight->GetErrorYhigh(i);
                EffNumErrL = hEffNumIntNoReweight->GetErrorYlow(i);
                EffDenErrH = hEffDenIntNoReweight->GetErrorYhigh(i);
                EffDenErrL = hEffDenIntNoReweight->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                hEffInt2S1SNoReweight->SetPoint(i, intBin[i], EffRatio);
                hEffInt2S1SNoReweight->SetPointError(i, intBinErr[i], intBinErr[i], EffRatioErrL, EffRatioErrH);


        }



	TFile* OutFile;
        OutFile = new TFile(Form("EffSingleRatio_%s.root",isPbPb ? "PbPb": "PP"), "Recreate");
        if(isPbPb){hEffCen2S1S->Write();}
        gEffRat_Pt_->Write();
        hEffRap2S1S->Write();
        hEffInt2S1S->Write();
        
        if(isPbPb){hEffCen2S1SNoReweight->Write();}
        gEffRat_Pt_NoReweight->Write();
        hEffRap2S1SNoReweight->Write();
        hEffInt2S1SNoReweight->Write();

        OutFile->Close();

	if(!isPbPb){
	iPeriod = 6;
	}


for (Int_t i = 0; i < (nCenBin); i++){
        EffCentValues[i] = hEffCen2S1S->Eval(CenBin[i]);
        xerrorCent[i] = 1.2;
}

        for (Int_t i = 0; i < (nCenBin); i++){
        SystematicErrorCent[i] = TMath::Abs(EffCentValues[i] - hEffCen2S1SNoReweight->Eval(CenBin[i]));
        cout << "Systematic Error in cent bin " << i+1 << " is " <<  SystematicErrorCent[i] << endl;
        }

TGraphErrors *EffCentSys = new TGraphErrors(nCenBin, CenBin, EffCentValues, xerrorCent, SystematicErrorCent);


	if(isPbPb){
        TCanvas *c1 = new TCanvas("c1","c1",800,600);
        c1->SetRightMargin(1);
	c1->cd();

	//adding a line
	TLine* line1 = new TLine(0,1,100,1);
        line1->SetLineStyle(kDashed);

        hEffCen2S1S->SetMarkerSize(2.0);
        hEffCen2S1S->SetMarkerColor(kRed);
        hEffCen2S1S->SetMarkerStyle(20);
	hEffCen2S1S->SetMarkerSize(2.0);
        hEffCen2S1S->SetTitle("");
        hEffCen2S1S->GetXaxis()->SetTitle("Centrality");
        hEffCen2S1S->GetXaxis()->CenterTitle();
	hEffCen2S1S->GetYaxis()->CenterTitle();
//	hEffCen2S1S->GetXaxis()->SetTitleOffset(1.5);
//	hEffCen2S1S->GetYaxis()->SetTitleOffset(1.8);
	hEffCen2S1S->GetYaxis()->SetTitle(Form("Eff^{#varUpsilon(2S)/#varUpsilon(1S)}_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffCen2S1S->GetYaxis()->SetRangeUser(0.95, 1.05);
	hEffCen2S1S->GetXaxis()->SetRangeUser(0.0, 100);
        hEffCen2S1S->Draw("AP");
	line1->Draw("");
EffCentSys->SetFillColor(2);
EffCentSys->SetFillStyle(3001);
EffCentSys->Draw("2");
	CMS_lumi(c1,iPeriod, iPos);
	c1->Update();

        c1->SaveAs(Form("SingleRatioEff_Cent_%s.png",isPbPb ? "PbPb" : "PP"));	
        }






        TCanvas *c2 = new TCanvas("c2","c2",800,600);
        c2->SetRightMargin(1);
        c2->cd();

        for (Int_t i = 0; i < (nPtBin); i++){
        SystematicErrorPt[i] = TMath::Abs(gEffRat_Pt_->Eval(ptBin[i]) - gEffRat_Pt_NoReweight->Eval(ptBin[i]));
        cout << "Systematic Error in pt bin " << i+1 << " is " <<  SystematicErrorPt[i] << endl;
        }

for (Int_t i = 0; i < (nPtBin); i++){
        EffPtValues[i] = gEffRat_Pt_->Eval(ptBin[i]);
        xerrorPt[i] = 0.45;
}

TGraphErrors *EffPtSys = new TGraphErrors(nPtBin, ptBin, EffPtValues, xerrorPt, SystematicErrorPt);


	TLine* line2 = new TLine(0,1,30,1);
        line2->SetLineStyle(kDashed);

        gEffRat_Pt_->SetMarkerSize(2.0);
        gEffRat_Pt_->SetMarkerColor(kRed);
        gEffRat_Pt_->SetMarkerStyle(20);
        gEffRat_Pt_->SetMarkerSize(2.0);
        gEffRat_Pt_->SetTitle("");
        gEffRat_Pt_->GetXaxis()->SetTitle("p_{T}");
        gEffRat_Pt_->GetXaxis()->CenterTitle();
	gEffRat_Pt_->GetYaxis()->CenterTitle();
//        gEffRat_Pt_->GetXaxis()->SetTitleOffset(1.5);
//        gEffRat_Pt_->GetYaxis()->SetTitleOffset(1.8);
	gEffRat_Pt_->GetYaxis()->SetTitle(Form("Eff^{#varUpsilon(2S)/#varUpsilon(1S)}_{%s}",isPbPb ? "PbPb" : "PP"));
	gEffRat_Pt_->GetYaxis()->SetRangeUser(0.95, 1.05);
	gEffRat_Pt_->GetXaxis()->SetRangeUser(0.0, 30);
	gEffRat_Pt_->Draw("AP");
	line2->Draw("");
EffPtSys->SetFillColor(2);
EffPtSys->SetFillStyle(3001);
EffPtSys->Draw("2");
        CMS_lumi(c2,iPeriod, iPos);
        c2->Update();

        c2->SaveAs(Form("SingleRatioEff_Pt_%s.png",isPbPb ? "PbPb" : "PP"));	

        TCanvas *c3 = new TCanvas("c3","c3",800,600);
        c3->SetRightMargin(1);
        c3->cd();

        for (Int_t i = 0; i < (nRapBin); i++){
        SystematicErrorRap[i] = TMath::Abs(hEffRap2S1S->Eval(rapBin[i]) - hEffRap2S1SNoReweight->Eval(rapBin[i]));
        cout << "Systematic Error in rap bin " << i+1 << " is " <<  SystematicErrorRap[i] << endl;
        }


for (Int_t i = 0; i < (nRapBin); i++){
        EffRapValues[i] = hEffRap2S1S->Eval(rapBin[i]);
        xerrorRap[i] = 0.03;
}

TGraphErrors *EffRapSys = new TGraphErrors(nRapBin, rapBin, EffRapValues, xerrorRap, SystematicErrorRap);


	TLine* line3 = new TLine(0,1,2.4,1);
        line3->SetLineStyle(kDashed);

        hEffRap2S1S->SetMarkerSize(2.0);
        hEffRap2S1S->SetMarkerColor(kRed);
        hEffRap2S1S->SetMarkerStyle(20);
        hEffRap2S1S->SetMarkerSize(2.0);
        hEffRap2S1S->SetTitle("");
        hEffRap2S1S->GetXaxis()->SetTitle("|y|");
        hEffRap2S1S->GetXaxis()->CenterTitle();
	hEffRap2S1S->GetYaxis()->CenterTitle();
//        hEffRap2S1S->GetXaxis()->SetTitleOffset(1.5);
//        hEffRap2S1S->GetYaxis()->SetTitleOffset(1.8);
	hEffRap2S1S->GetYaxis()->SetTitle(Form("Eff^{#varUpsilon(2S)/#varUpsilon(1S)}_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffRap2S1S->GetYaxis()->SetRangeUser(0.95, 1.05);
	hEffRap2S1S->GetXaxis()->SetRangeUser(0.0, 2.4);
	hEffRap2S1S->Draw("AP");
	line3->Draw("");
EffRapSys->SetFillColor(2);
EffRapSys->SetFillStyle(3001);
EffRapSys->Draw("2");
        CMS_lumi(c3,iPeriod, iPos);
        c3->Update();

        c3->SaveAs(Form("SingleRatioEff_Rap_%s.png",isPbPb ? "PbPb" : "PP"));	


        TCanvas *c4 = new TCanvas("c4","c4",800,600);
        c4->SetRightMargin(1);
        c4->cd();


        for (Int_t i = 0; i < 1; i++){
        SystematicErrorInt[i] = TMath::Abs(hEffInt2S1S->Eval(intBin[i]) - hEffInt2S1SNoReweight->Eval(intBin[i]));
        cout << "Systematic Error in Int bin " << i+1 << " is " <<  SystematicErrorInt[i] << endl;
        }

for (Int_t i = 0; i < 1; i++){
        EffIntValues[i] = hEffInt2S1S->Eval(intBin[i]);
        xerrorInt[i] = 1.2;
}

TGraphErrors *EffIntSys = new TGraphErrors(1, intBin, EffIntValues, xerrorInt, SystematicErrorInt);


	TLine* line4 = new TLine(0,1,100,1);
        line4->SetLineStyle(kDashed);

        hEffInt2S1S->SetMarkerSize(2.0);
        hEffInt2S1S->SetMarkerColor(kRed);
        hEffInt2S1S->SetMarkerStyle(20);
        hEffInt2S1S->SetMarkerSize(2.0);
        hEffInt2S1S->SetTitle("");
        hEffInt2S1S->GetXaxis()->SetTitle("Integrated");
        hEffInt2S1S->GetXaxis()->CenterTitle();
	hEffInt2S1S->GetYaxis()->CenterTitle();
//        hEffInt2S1S->GetXaxis()->SetTitleOffset(1.5);
//        hEffInt2S1S->GetYaxis()->SetTitleOffset(1.8);
	hEffInt2S1S->GetYaxis()->SetTitle(Form("Eff^{#varUpsilon(2S)/#varUpsilon(1S)}_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffInt2S1S->GetYaxis()->SetRangeUser(0.95, 1.05);
	hEffInt2S1S->GetXaxis()->SetRangeUser(0.0, 100);
	hEffInt2S1S->Draw("AP");
	line4->Draw("");
EffIntSys->SetFillColor(2);
EffIntSys->SetFillStyle(3001);
EffIntSys->Draw("2");
        CMS_lumi(c4,iPeriod, iPos);
        c4->Update();

        c4->SaveAs(Form("SingleRatioEff_Int_%s.png",isPbPb ? "PbPb" : "PP"));	


        if(isPbPb){
        cout <<"doing Centrality"<<endl;
        for (Int_t i = 0; i < (nCenBin); i++){
        cout << hEffCen2S1S->Eval(CenBin[i]) << " , - " << hEffCen2S1S->GetErrorYlow(i) << " , + " << hEffCen2S1S->GetErrorYhigh(i) << endl;
        }
	}
        cout <<"doing Int"<<endl;
        for (Int_t i = 0; i < (1); i++){
        cout << hEffInt2S1S->Eval(intBin[i]) << " , - " << hEffInt2S1S->GetErrorYlow(i) << " , + " << hEffInt2S1S->GetErrorYhigh(i) << endl;
        }
        cout <<"doing PT"<<endl;
        for (Int_t i = 0; i < (nPtBin); i++){
        cout << gEffRat_Pt_->Eval(ptBin[i]) << " , - " << gEffRat_Pt_->GetErrorYlow(i) << " , + " << gEffRat_Pt_->GetErrorYhigh(i) << endl;
        }
        cout <<"doing Rap"<<endl;
        for (Int_t i = 0; i < (nRapBin); i++){
        cout << hEffRap2S1S->Eval(rapBin[i]) << " , - " << hEffRap2S1S->GetErrorYlow(i) << " , + " << hEffRap2S1S->GetErrorYhigh(i) << endl;
        }

        cout << "over" << endl;
	
	fEff1S->Close(); 
	fEff2S->Close(); 
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

// Relative Error
double RelError(double A, double B){
	double f = abs(A-B)/A;
	return f;
}
