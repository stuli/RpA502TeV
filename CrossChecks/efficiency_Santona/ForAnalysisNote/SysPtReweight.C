#include "effCommon.h"



double RError(double A, double eA, double B, double eB);
double PError(double A, double eA, double B, double eB);

const int  nCenBin = 9;
const int  nPtBin = 3;
const int  nRapBin = 2;
const int  nIntBin = 1;

int iPeriod = 5;
int iPos = 33;

void SystematicssingleRatioEff(bool isPbPb = false){
//        gROOT->Macro("logon.C+");
	setTDRStyle();

        float          CenBin[nCenBin] = {2.5,7.5,15,25,35,45,55,65,85};
        float          CenBinErr[nCenBin] = {2.5,2.5,5,5,5,5,5,5,15};

        float          ptBin[nPtBin] = {2.5,8.5,21};
        float          ptBinErr[nPtBin] = {2.5,3.5,9};

        float          rapBin[nRapBin] = {0.6,1.8};
        float          rapBinErr[nRapBin] = {0.6,0.6};

        float          intBin[1] = {50};
        float          intBinErr[1] = {50};

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

	TGraphAsymmErrors* hEffNumCent;
	TGraphAsymmErrors* hEffNumInt;
	TGraphAsymmErrors* hEffNumPt;
	TGraphAsymmErrors* hEffNumRap;

        TGraphAsymmErrors* hEffDenCent;
        TGraphAsymmErrors* hEffDenInt;
        TGraphAsymmErrors* hEffDenPt;
        TGraphAsymmErrors* hEffDenRap;

        TGraphAsymmErrors* hEffCen2S1S = new TGraphAsymmErrors(nCenBin);
	hEffCen2S1S->SetName("hEffCenSingle");
        TGraphAsymmErrors* hEffPt2S1S = new TGraphAsymmErrors(nPtBin);
	hEffPt2S1S->SetName("hEffPtSingle");
        TGraphAsymmErrors* hEffRap2S1S = new TGraphAsymmErrors(nRapBin);
	hEffRap2S1S->SetName("hEffRapSingle");
        TGraphAsymmErrors* hEffInt2S1S = new TGraphAsymmErrors(nIntBin);
	hEffInt2S1S->SetName("hEffIntSingle");

        TGraphAsymmErrors* hEffNumCentNoReweight;
        TGraphAsymmErrors* hEffNumIntNoReweight;
        TGraphAsymmErrors* hEffNumPtNoReweight;
        TGraphAsymmErrors* hEffNumRapNoReweight;

        TGraphAsymmErrors* hEffDenCentNoReweight;
        TGraphAsymmErrors* hEffDenIntNoReweight;
        TGraphAsymmErrors* hEffDenPtNoReweight;
        TGraphAsymmErrors* hEffDenRapNoReweight;

        TGraphAsymmErrors* hEffCen2S1SNoReweight = new TGraphAsymmErrors(nCenBin);
        hEffCen2S1SNoReweight->SetName("hEffCenSingleNoReweight");
        TGraphAsymmErrors* hEffPt2S1SNoReweight = new TGraphAsymmErrors(nPtBin);
        hEffPt2S1SNoReweight->SetName("hEffPtSingleNoReweight");
        TGraphAsymmErrors* hEffRap2S1SNoReweight = new TGraphAsymmErrors(nRapBin);
        hEffRap2S1SNoReweight->SetName("hEffRapSingleNoReweight");
        TGraphAsymmErrors* hEffInt2S1SNoReweight = new TGraphAsymmErrors(nIntBin);
        hEffInt2S1SNoReweight->SetName("hEffIntSingleNoReweight");


//--- pp 2S

        TFile* fEff2S = new TFile(Form("Eff_%s_2S.root",isPbPb ? "PbPb" : "PP"), "Open");
	fEff2S->GetObject("EffCent", hEffNumCent);
	fEff2S->GetObject("EffInt", hEffNumInt);
	fEff2S->GetObject("EffPt", hEffNumPt);
	fEff2S->GetObject("EffRap", hEffNumRap);

        fEff2S->GetObject("EffCentNoReweight", hEffNumCentNoReweight);
        fEff2S->GetObject("EffIntNoReweight", hEffNumIntNoReweight);
        fEff2S->GetObject("EffPtNoReweight", hEffNumPtNoReweight);
        fEff2S->GetObject("EffRapNoReweight", hEffNumRapNoReweight);

//--- pp 1S

	TFile* fEff1S = new TFile(Form("Eff_%s_1S.root",isPbPb ? "PbPb" : "PP"), "Open");
        fEff1S->GetObject("EffCent", hEffDenCent);
        fEff1S->GetObject("EffInt", hEffDenInt);
        fEff1S->GetObject("EffPt", hEffDenPt);
        fEff1S->GetObject("EffRap", hEffDenRap);

        fEff1S->GetObject("EffCentNoReweight", hEffDenCentNoReweight);
        fEff1S->GetObject("EffIntNoReweight", hEffDenIntNoReweight);
        fEff1S->GetObject("EffPtNoReweight", hEffDenPtNoReweight);
        fEff1S->GetObject("EffRapNoReweight", hEffDenRapNoReweight);


//--- Single Ratio 2S/1S Calculation

        double EffRatio;
        double EffNum;
        double EffNumErrH;
        double EffDenErrH;
        double EffNumErrL;
        double EffDenErrL;
        double EffDen;
	double EffRatioErrH;
        double EffRatioErrL;


	if(isPbPb){
        for (Int_t i = 0; i < nCenBin; i++){
                EffNum = hEffNumCent->Eval(CenBin[i]);
                EffDen = hEffDenCent->Eval(CenBin[i]);
                EffNumErrH = hEffNumCent->GetErrorYhigh(i);
                EffNumErrL = hEffNumCent->GetErrorYlow(i);
                EffDenErrH = hEffDenCent->GetErrorYhigh(i);
                EffDenErrL = hEffDenCent->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffCen2S1S->SetPoint(i, CenBin[i], EffRatio);
                hEffCen2S1S->SetPointError(i, CenBinErr[i], CenBinErr[i], EffRatioErrL, EffRatioErrH);

		
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

        for (Int_t i = 0; i < nCenBin; i++){
                EffNum = hEffNumCentNoReweight->Eval(CenBin[i]);
                EffDen = hEffDenCentNoReweight->Eval(CenBin[i]);
                EffNumErrH = hEffNumCentNoReweight->GetErrorYhigh(i);
                EffNumErrL = hEffNumCentNoReweight->GetErrorYlow(i);
                EffDenErrH = hEffDenCentNoReweight->GetErrorYhigh(i);
                EffDenErrL = hEffDenCentNoReweight->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                hEffCen2S1SNoReweight->SetPoint(i, CenBin[i], EffRatio);
                hEffCen2S1SNoReweight->SetPointError(i, CenBinErr[i], CenBinErr[i], EffRatioErrL, EffRatioErrH);


        }
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


        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = hEffNumPt->Eval(ptBin[i]);
                EffDen = hEffDenPt->Eval(ptBin[i]);
                EffNumErrH = hEffNumPt->GetErrorYhigh(i);
                EffNumErrL = hEffNumPt->GetErrorYlow(i);
                EffDenErrH = hEffDenPt->GetErrorYhigh(i);
                EffDenErrL = hEffDenPt->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffPt2S1S->SetPoint(i, ptBin[i], EffRatio);
                hEffPt2S1S->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);


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


        for (Int_t i = 0; i < nPtBin; i++){
                EffNum = hEffNumPtNoReweight->Eval(ptBin[i]);
                EffDen = hEffDenPtNoReweight->Eval(ptBin[i]);
                EffNumErrH = hEffNumPtNoReweight->GetErrorYhigh(i);
                EffNumErrL = hEffNumPtNoReweight->GetErrorYlow(i);
                EffDenErrH = hEffDenPtNoReweight->GetErrorYhigh(i);
                EffDenErrL = hEffDenPtNoReweight->GetErrorYlow(i);
                EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH);

                hEffPt2S1SNoReweight->SetPoint(i, ptBin[i], EffRatio);
                hEffPt2S1SNoReweight->SetPointError(i, ptBinErr[i], ptBinErr[i], EffRatioErrL, EffRatioErrH);


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
        hEffPt2S1S->Write();
        hEffRap2S1S->Write();
        hEffInt2S1S->Write();
        
        if(isPbPb){hEffCen2S1SNoReweight->Write();}
        hEffPt2S1SNoReweight->Write();
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
        SystematicErrorPt[i] = TMath::Abs(hEffPt2S1S->Eval(ptBin[i]) - hEffPt2S1SNoReweight->Eval(ptBin[i]));
        cout << "Systematic Error in pt bin " << i+1 << " is " <<  SystematicErrorPt[i] << endl;
        }

for (Int_t i = 0; i < (nPtBin); i++){
        EffPtValues[i] = hEffPt2S1S->Eval(ptBin[i]);
        xerrorPt[i] = 0.45;
}

TGraphErrors *EffPtSys = new TGraphErrors(nPtBin, ptBin, EffPtValues, xerrorPt, SystematicErrorPt);


	TLine* line2 = new TLine(0,1,30,1);
        line2->SetLineStyle(kDashed);

        hEffPt2S1S->SetMarkerSize(2.0);
        hEffPt2S1S->SetMarkerColor(kRed);
        hEffPt2S1S->SetMarkerStyle(20);
        hEffPt2S1S->SetMarkerSize(2.0);
        hEffPt2S1S->SetTitle("");
        hEffPt2S1S->GetXaxis()->SetTitle("p_{T}");
        hEffPt2S1S->GetXaxis()->CenterTitle();
	hEffPt2S1S->GetYaxis()->CenterTitle();
//        hEffPt2S1S->GetXaxis()->SetTitleOffset(1.5);
//        hEffPt2S1S->GetYaxis()->SetTitleOffset(1.8);
	hEffPt2S1S->GetYaxis()->SetTitle(Form("Eff^{#varUpsilon(2S)/#varUpsilon(1S)}_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffPt2S1S->GetYaxis()->SetRangeUser(0.95, 1.05);
	hEffPt2S1S->GetXaxis()->SetRangeUser(0.0, 30);
	hEffPt2S1S->Draw("AP");
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
        cout << hEffPt2S1S->Eval(ptBin[i]) << " , - " << hEffPt2S1S->GetErrorYlow(i) << " , + " << hEffPt2S1S->GetErrorYhigh(i) << endl;
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


