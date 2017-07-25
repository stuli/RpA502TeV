//Making code a littler cleaner
#include "effCommon.h"



double RError(double A, double eA, double B, double eB);
double PError(double A, double eA, double B, double eB);

const int  nCenBin = 9;
const int  nPtBin = 3;
const int  nRapBin = 2;
const int  nIntBin = 1;


void singleRatioEff(bool isPbPb = true){
        gROOT->Macro("logon.C+");

        double          CenBin[nCenBin] = {5,15,30,50,70,90,110,130,170};
        double          CenBinErr[nCenBin] = {5,5,10,10,10,10,10,10,30};

        double          ptBin[nPtBin] = {2.5,8.5,21};
        double          ptBinErr[nPtBin] = {2.5,3.5,9};

        double          rapBin[nRapBin] = {0.6,1.8};
        double          rapBinErr[nRapBin] = {0.6,0.6};

        double          intBin[1] = {100};
        double          intBinErr[1] = {100};

//Declaring histograms
/*
	TH1D  *RecoEvents;
        TH1D  *GenEvents;
        TH1D  *RecoEventsInt;
        TH1D  *GenEventsInt;
        TH1D  *RecoEventsPt;
        TH1D  *GenEventsPt;
        TH1D  *RecoEventsRap;
        TH1D  *GenEventsRap;
// */

	TH1D* hNumCent;
	TH1D* hNumInt;
	TH1D* hNumPt;
	TH1D* hNumRap;

        TH1D* hDenCent;
        TH1D* hDenInt;
        TH1D* hDenPt;
        TH1D* hDenRap;

        TH1D* hEffCen = new TH1D(nCenBin);
        TH1D* hEffPt = new TH1D(nPtBin);
        TH1D* hEffRap = new TH1D(nRapBin);
        TH1D* hEffInt = new TH1D(nIntBin);


//--- pp 2S

        TFile* fEff2S = new TFile(Form("Eff_%s_2S.root",isPbPb ? "PbPb" : "PP"), "Open");
	fEff2S->GetObject("RecoEvents", hNumCent);
        fEff2S->GetObject("GenEvents", hDenCent);
	fEff2S->GetObject("RecoEventsInt", hNumInt);
        fEff2S->GetObject("GenEventsInt", hDenInt);
	fEff2S->GetObject("RecoEventsPt", hNumPt);
         fEff2S->GetObject("GenEventsPt", hDenPt);
	fEff2S->GetObject("RecoEventsRap", hNumRap);
         fEff2S->GetObject("GenEventsRap", hDenRap);


//--- pp 1S

	TFile* fEff1S = new TFile(Form("Eff_%s_1S.root",isPbPb ? "PbPb" : "PP"), "Open");
        fEff1S->GetObject("EffCent", hDenCent);
        fEff1S->GetObject("EffInt", hDenInt);
        fEff1S->GetObject("EffPt", hDenPt);
        fEff1S->GetObject("EffRap", hDenRap);

// Getting Efficiency 

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
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  //Calculation for the combined efficiency	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffCen2S1S->SetPoint(i, CenBin[i], EffRatio);
                hEffCen2S1S->SetPointError(i, CenBinErr[i], CenBinErr[i], EffRatioErrL, EffRatioErrH);

		
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
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  //Calculation for the combined efficiency	
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

        for (Int_t i = 0; i < nRapBin; i++){
                EffNum = hEffNumRap->Eval(rapBin[i]);
                EffDen = hEffDenRap->Eval(rapBin[i]);
                EffNumErrH = hEffNumRap->GetErrorYhigh(i);
                EffNumErrL = hEffNumRap->GetErrorYlow(i);
                EffDenErrH = hEffDenRap->GetErrorYhigh(i);
                EffDenErrL = hEffDenRap->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  //Calculation for the combined efficiency	
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

        for (Int_t i = 0; i < nIntBin; i++){
                EffNum = hEffNumInt->Eval(intBin[i]);
                EffDen = hEffDenInt->Eval(intBin[i]);
                EffNumErrH = hEffNumInt->GetErrorYhigh(i);
                EffNumErrL = hEffNumInt->GetErrorYlow(i);
                EffDenErrH = hEffDenInt->GetErrorYhigh(i);
                EffDenErrL = hEffDenInt->GetErrorYlow(i);
		EffRatio = EffNum / EffDen;
                EffRatioErrH = RError(EffNum, EffNumErrH, EffDen, EffDenErrL);  //Calculation for the combined efficiency	
                EffRatioErrL = RError(EffNum, EffNumErrL, EffDen, EffDenErrH); 

                hEffInt2S1S->SetPoint(i, intBin[i], EffRatio);
                hEffInt2S1S->SetPointError(i, intBinErr[i], intBinErr[i], EffRatioErrL, EffRatioErrH);


        }



	TFile* OutFile;
        OutFile = new TFile(Form("EffSingleRatio_%s.root",isPbPb ? "PbPb": "PP"), "Recreate");
        if(isPbPb){hEffCen2S1S->Write();}
        hEffPt2S1S->Write();
        hEffRap2S1S->Write();
        hEffInt2S1S->Write();
        
        OutFile->Close();

	if(isPbPb){
        TCanvas* c1 = new TCanvas("c1", "Canvas with results1", 600, 600);
        c1->cd();

	//adding a line
	TLine* line1 = new TLine(0,1,200,1);
        line1->SetLineStyle(kDashed);

        //hEffCen2S1S->SetMarkerSize(2.0);
        hEffCen2S1S->SetMarkerColor(kRed);
        hEffCen2S1S->SetMarkerStyle(20);
        hEffCen2S1S->SetTitle("");
        hEffCen2S1S->GetXaxis()->SetTitle("Centrality");
        hEffCen2S1S->GetXaxis()->CenterTitle();
	hEffCen2S1S->GetYaxis()->CenterTitle();
	hEffCen2S1S->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(2S)/#varUpsilon(1S)]_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffCen2S1S->GetYaxis()->SetRangeUser(0.5, 1.5);
	hEffCen2S1S->GetXaxis()->SetRangeUser(0.0, 200);
        hEffCen2S1S->Draw("AP");
	line1->Draw("");

        c1->SaveAs(Form("SingleRatioEff_Cent_%s.png",isPbPb ? "PbPb" : "PP"));	
        }
        TCanvas* c2 = new TCanvas("c2", "Canvas with results1", 600, 600);
        c2->cd();

	TLine* line2 = new TLine(0,1,30,1);
        line2->SetLineStyle(kDashed);

        //hEffPt2S1S->SetMarkerSize(1.0);
        hEffPt2S1S->SetMarkerColor(kRed);
        hEffPt2S1S->SetMarkerStyle(20);
        hEffPt2S1S->SetTitle("");
        hEffPt2S1S->GetXaxis()->SetTitle("p_{T}");
        hEffPt2S1S->GetXaxis()->CenterTitle();
	hEffPt2S1S->GetYaxis()->CenterTitle();
	hEffPt2S1S->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(2S)/#varUpsilon(1S)]_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffPt2S1S->GetYaxis()->SetRangeUser(0.5, 1.5);
	hEffPt2S1S->GetXaxis()->SetRangeUser(0.0, 30);
	hEffPt2S1S->Draw("AP");
	line2->Draw("");

        c2->SaveAs(Form("SingleRatioEff_Pt_%s.png",isPbPb ? "PbPb" : "PP"));	

	TCanvas* c3 = new TCanvas("c3", "Canvas with results1", 600, 600);
        c3->cd();

	TLine* line3 = new TLine(0,1,2.4,1);
        line3->SetLineStyle(kDashed);

        //hEffRap2S1S->SetMarkerSize(2.0);
        hEffRap2S1S->SetMarkerColor(kRed);
        hEffRap2S1S->SetMarkerStyle(20);
        hEffRap2S1S->SetTitle("");
        hEffRap2S1S->GetXaxis()->SetTitle("|y|");
        hEffRap2S1S->GetXaxis()->CenterTitle();
	hEffRap2S1S->GetYaxis()->CenterTitle();
	hEffRap2S1S->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(2S)/#varUpsilon(1S)]_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffRap2S1S->GetYaxis()->SetRangeUser(0.5, 1.5);
	hEffRap2S1S->GetXaxis()->SetRangeUser(0.0, 2.4);
	hEffRap2S1S->Draw("AP");
	line3->Draw("");

        c3->SaveAs(Form("SingleRatioEff_Rap_%s.png",isPbPb ? "PbPb" : "PP"));	


	TCanvas* c4 = new TCanvas("c4", "Canvas with results1", 600, 600);
        c4->cd();

	TLine* line4 = new TLine(0,1,200,1);
        line4->SetLineStyle(kDashed);

        //hEffInt2S1S->SetMarkerSize(2.0);
        hEffInt2S1S->SetMarkerColor(kRed);
        hEffInt2S1S->SetMarkerStyle(20);
        hEffInt2S1S->SetTitle("");
        hEffInt2S1S->GetXaxis()->SetTitle("Integrated");
        hEffInt2S1S->GetXaxis()->CenterTitle();
	hEffInt2S1S->GetYaxis()->CenterTitle();
	hEffInt2S1S->GetYaxis()->SetTitle(Form("Efficiency[#varUpsilon(2S)/#varUpsilon(1S)]_{%s}",isPbPb ? "PbPb" : "PP"));
	hEffInt2S1S->GetYaxis()->SetRangeUser(0.5, 1.5);
	hEffInt2S1S->GetXaxis()->SetRangeUser(0.0, 200);
	hEffInt2S1S->Draw("AP");
	line4->Draw("");


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


