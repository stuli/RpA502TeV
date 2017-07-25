#include "effCommon.h"

double RError(double A, double eA, double B, double eB);


int  nPtBin;
int  nRapBin;
int nCentBin;
int nIntBin;

void Comparison(
        int oniaMode = 2, //1 = 1S, 2 = 2S, 3 = 3S
        bool isPbPb = true //true = PbPb and false = pp
        ){

  setTDRStyle();

// Define Cent, Int, Pt, Rap Arrays

const int nPtBins1s  = 5;
const int nPtBins2s  = 3; 
const int nPtBins3s  = 3;

const int nYBins1S  = 6;
const int nYBins2S  = 2;
const int nYBins3S  = 2;

const int nCentBins1s2s = 9;
const int nCentBins3s = 4;

        float           IntBin[1] = { 100 };
        float           IntBinEdges[2] = { 0, 100 };
	float		IntBinEdgesYongsun[2] = {0, 200};

if (oniaMode == 1){
nPtBin = nPtBins1s;
cout << "This is fine" << endl;
float PtBin[nPtBins1s] = {1.25,3.75,6.5,11.5,22.5};
        float          PtBinEdges[nPtBins1s + 1] ={0,2.5,5,8,15,30};
nRapBin = nYBins1S;
float          RapBin[nYBins1S] = {0.2,0.6,1.0,1.4,1.8,2.2};
float          RapBinEdges[nYBins1S + 1] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
nCentBin = nCentBins1s2s;
        float           CentBin[nCentBins1s2s] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65, 85 };
        float           CentBinEdges[nCentBins1s2s + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 100 };
	float		CentBinYongsun[nCentBins1s2s] = { 5, 15, 30, 50, 70, 90, 110, 130, 170 };
	float		CentBinEdgesYongsun[nCentBins1s2s+1] = {0,10,20,40,60,80,100,120,140,200};
}
if (oniaMode == 2){
nPtBin = nPtBins2s;
float PtBin[nPtBins2s] = {2.5,10,22.5};
        float          PtBinEdges[nPtBins2s + 1] ={0,5,15,30};
nRapBin = nYBins2S;
float          RapBin[nYBins2S] = { 0.6, 1.8 };
float          RapBinEdges[nYBins2S + 1] = { 0, 1.2, 2.4 };
nCentBin = nCentBins1s2s;
        float           CentBin[nCentBins1s2s] = { 2.5, 7.5, 15, 25, 35, 45, 55, 65, 85 };
        float           CentBinEdges[nCentBins1s2s + 1] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 100 };
        float           CentBinYongsun[nCentBins1s2s] = { 5, 15, 30, 50, 70, 90, 110, 130, 170 };
        float           CentBinEdgesYongsun[nCentBins1s2s+1] = {0,10,20,40,60,80,100,120,140,200};
}
if (oniaMode == 3){
nPtBin = nPtBins3s;
float PtBin[nPtBins3s] = {2.5,10,22.5};
        float          PtBinEdges[nPtBins3s + 1] ={0,5,15,30};
nRapBin = nYBins3S;
float          RapBin[nYBins3S] = { 0.6, 1.8 };
float          RapBinEdges[nYBins3S + 1] = { 0, 1.2, 2.4 };
nCentBin  = nCentBins3s;
        float           CentBin[nCentBins3s] = { 5, 20, 40, 75};
	float CentBinEdges[nCentBins3s+1] = {0,10,30,50,100};
        float           CentBinYongsun[nCentBins3s] = { 10, 40, 80, 150 };
        float           CentBinEdgesYongsun[nCentBins3s+1] = {0,20,60,100,200};
}


// I will get my efficiencies directly from the TGraphs and turn them into histograms. Then I will get Yongsun's efficiencies as histograms and find their raito.

// Mine
        TGraphAsymmErrors* hEffNumCent;
        TGraphAsymmErrors* hEffNumInt;
        TGraphAsymmErrors* hEffNumPt;
        TGraphAsymmErrors* hEffNumRap;

	if(isPbPb){
        TH1D* h1EffNumCent = new TH1D("h1_EffNumCent","",nCentBin,CentBinEdges);
	}
        TH1D* h1EffNumInt = new TH1D("h1_EffNumInt","",nIntBin,IntBinEdges);
        TH1D* h1EffNumPt = new TH1D("h1_EffNumPt","",nPtBin,PtBinEdges);
        TH1D* h1EffNumRap = new TH1D("h1_EffNumRap","",nRapBin,RapBinEdges);

// Yongsun's
        TGraphAsymmErrors* hEffDenCent;
        TGraphAsymmErrors* hEffDenInt;
        TGraphAsymmErrors* hEffDenPt;
        TGraphAsymmErrors* hEffDenRap;

	if(isPbPb){
        TH1D* h1EffDenCentYongsun = new TH1D("h1_EffDenCentYongsun","",nCentBin,CentBinEdgesYongsun);
	TH1D* h1EffDenCent = new TH1D("h1_EffDenCent","",nCentBin,CentBinEdges);
	}
        TH1D* h1EffDenInt = new TH1D("h1_EffDenInt","",nIntBin,IntBinEdges);
        TH1D* h1EffDenIntYongsun = new TH1D("h1_EffDenIntYongsun","",nIntBin,IntBinEdgesYongsun);
        TH1D* h1EffDenPt = new TH1D("h1_EffDenPt","",nPtBin,PtBinEdges);
        TH1D* h1EffDenRap = new TH1D("h1_EffDenRap","",nRapBin,RapBinEdges);


// Declaring Ratio histograms
        if (isPbPb){
        TH1D* hRatioCent = new TH1D("",";Centrality;Crosscheck Ratio",nCentBin, CentBinEdges);
	}
        TH1D* hRatioPt = new TH1D("",";p_{T} [GeV/c];Crosscheck Ratio",nPtBin, PtBinEdges);
        TH1D* hRatioRap = new TH1D("",";|y|;Crosscheck Ratio",nRapBin, RapBinEdges);
        TH1D* hRatioInt = new TH1D("",";Integrated;Crosscheck Ratio",nIntBin, IntBinEdges);

// Getting my TGraphs

        TFile* fEffSantona = new TFile(Form("Eff_%s_%sS.root",isPbPb ? "PbPb" : "PP", oniaMode), "Open");
	if (isPbPb){
        fEffSantona->GetObject("EffCent", hEffNumCent);
	}
        fEffSantona->GetObject("EffInt", hEffNumInt);
        fEffSantona->GetObject("EffPt", hEffNumPt);
        fEffSantona->GetObject("EffRap", hEffNumRap);

// Turning my TGraphs into Histograms

  if (isPbPb){

 //Cent
  for (Int_t i = 0; i < nCentBin ; i++){
    double xCent;
    double EffNumCent;
    hEffNumCent->GetPoint(i,xCent,EffNumCent);

    double EffNumErrH = hEffNumCent->GetErrorYhigh(i);
    double EffNumErrL = hEffNumCent->GetErrorYlow(i);

    cout << "x / y position of EffNumCent " << xCent << "\t" << EffNumCent << endl;
    if (xCent>CentBinEdges[nCentBin]) break;

    int bin = h1EffNumCent->FindBin(xCent);
    double error = EffNumErrH>EffNumErrL ? EffNumErrH : EffNumErrL;
    h1EffNumCent->SetBinContent(bin, EffNumCent);
    h1EffNumCent->SetBinError(bin, error);
    error = EffDenErrH>EffDenErrL ? EffDenErrH : EffDenErrL;

    cout << "bin:\t" << bin << " ";
    cout << "NumCent hist:\t" << h1EffNumCent->GetBinContent(bin)<<"\t";
    cout << endl;
  }
  }

 //Int
  for (Int_t i = 0; i < nIntBin ; i++){
    double xInt;
    double EffNumInt;
    hEffNumInt->GetPoint(i,xInt,EffNumInt);

    double EffNumErrH = hEffNumInt->GetErrorYhigh(i);
    double EffNumErrL = hEffNumInt->GetErrorYlow(i);

    cout << "x / y position of EffNumInt " << xInt << "\t" << EffNumInt << endl;
    if (xInt>IntBinEdges[nIntBin]) break;

    int bin = h1EffNumInt->FindBin(xInt);
    double error = EffNumErrH>EffNumErrL ? EffNumErrH : EffNumErrL;
    h1EffNumInt->SetBinContent(bin, EffNumInt);
    h1EffNumInt->SetBinError(bin, error);
    error = EffDenErrH>EffDenErrL ? EffDenErrH : EffDenErrL;

    cout << "bin:\t" << bin << " ";
    cout << "NumInt hist:\t" << h1EffNumInt->GetBinContent(bin)<<"\t";
    cout << endl;
  }

 //Pt
  for (Int_t i = 0; i < nPtBin ; i++){
    double xPt;
    double EffNumPt;
    hEffNumPt->GetPoint(i,xPt,EffNumPt);

    double EffNumErrH = hEffNumPt->GetErrorYhigh(i);
    double EffNumErrL = hEffNumPt->GetErrorYlow(i);

    cout << "x / y position of EffNumPt " << xPt << "\t" << EffNumPt << endl;
    if (xPt>PtBinEdges[nPtBin]) break;

    int bin = h1EffNumPt->FindBin(xPt);
    double error = EffNumErrH>EffNumErrL ? EffNumErrH : EffNumErrL;
    h1EffNumPt->SetBinContent(bin, EffNumPt);
    h1EffNumPt->SetBinError(bin, error);
    error = EffDenErrH>EffDenErrL ? EffDenErrH : EffDenErrL;

    cout << "bin:\t" << bin << " ";
    cout << "NumPt hist:\t" << h1EffNumPt->GetBinContent(bin)<<"\t";
    cout << endl;
  }

 //Rap
   for (Int_t i = 0; i < nRapBin ; i++){
    double xRap;
    double EffNumRap;
    hEffNumRap->GetPoint(i,xRap,EffNumRap);

    double EffNumErrH = hEffNumRap->GetErrorYhigh(i);
    double EffNumErrL = hEffNumRap->GetErrorYlow(i);

    cout << "x / y position of EffNumRap " << xRap << "\t" << EffNumRap << endl;
    if (xRap>RapBinEdges[nRapBin]) break;

    int bin = h1EffNumRap->FindBin(xRap);
    double error = EffNumErrH>EffNumErrL ? EffNumErrH : EffNumErrL;
    h1EffNumRap->SetBinContent(bin, EffNumRap);
    h1EffNumRap->SetBinError(bin, error);
    error = EffDenErrH>EffDenErrL ? EffDenErrH : EffDenErrL;

    cout << "bin:\t" << bin << " ";
    cout << "NumRap hist:\t" << h1EffNumRap->GetBinContent(bin)<<"\t";
    cout << endl;
  }


// Getting Yongsun's Histograms

        TFile* fEffYongsun = new TFile(Form("efficiency_ups%ss_MCDATA.root", oniaMode), "Open");
    if (!isPbPb) {
        fEffYongsun->GetObject("hptEffPP", h1EffDenPt);
        fEffYongsun->GetObject("hrapEffPP", h1EffDenRap);
        fEffYongsun->GetObject("hcentEffPP", h1EffDenIntYongsun);
    }
    else {
	fEffYongsun->GetObject("hptEffAA", h1EffDenPt);
        fEffYongsun->GetObject("hrapEffAA", h1EffDenRap);
        fEffYongsun->GetObject("hcentEffAA", h1EffDenCentYongsun);
        fEffYongsun->GetObject("hcentEffAA_int", h1EffDenIntYongsun);
    }

 // Need to change Yongsun's centrality binning to mine. Will do this by reassigning the bins in h1EffDen_

    if (isPbPb) {
	for (Int_t i = 0; i < nCentBin ; i++){
	    	double EffTmp;
		EffTmp = h1EffDenCentYongsun->GetBinContent(i);
		h1EffDenCent->SetBinContent(i,EffTmp);
	}
    }

        for (Int_t i = 0; i < nIntBin ; i++){
                double EffTmpInt;
                EffTmpInt = h1EffDenIntYongsun->GetBinContent(i);
                h1EffDenInt->SetBinContent(i,EffTmpInt);
        }


// Getting Ratio Histograms

        if (isPbPb){
        hRatioCent->Divide(h1EffNumCent,h1EffDenCent);
        } 
        hRatioInt->Divide(h1EffNumInt,h1EffDenInt);
        hRatioPt->Divide(h1EffNumPt,h1EffDenPt);
        hRatioRap->Divide(h1EffNumRap,h1EffDenRap);


// Plotting

  TCanvas *canv2 = new TCanvas("canv2","canv2",800,800);
  TCanvas *canv3 = new TCanvas("canv3","canv3",800,800);
  TCanvas *canv4 = new TCanvas("canv4","canv4",800,800);

  if(isPbPb){
  TCanvas *canv1 = new TCanvas("canv","canv",800,800);
  canv1->cd();
  hRatioCent->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioCent->GetXaxis()->SetRangeUser(CentBinEdges[0], CentBinEdges[nCentBin]);
  hRatioCent->Draw("AP");
  canv1->Update;
  canv1->SaveAs(Form("EfficiencyCrossCheckCent_%dS_%s.png",oniaMode,isPbPb ? "PbPb" : "PP"));
  }

  canv2->cd();
  hRatioInt->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioInt->GetXaxis()->SetRangeUser(IntBinEdges[0], IntBinEdges[nIntBin]);
  hRatioInt->Draw("AP");
  canv2->Update;
  canv2->SaveAs(Form("EfficiencyCrossCheckInt_%dS_%s.png",oniaMode,isPbPb ? "PbPb" : "PP"));

  canv3->cd();
  hRatioPt->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioPt->GetXaxis()->SetRangeUser(PtBinEdges[0], PtBinEdges[nPtBin]);
  hRatioPt->Draw("AP");
  canv3->Update;
  canv3->SaveAs(Form("EfficiencyCrossCheckPt_%dS_%s.png",oniaMode,isPbPb ? "PbPb" : "PP"));

  canv4->cd();
  hRatioRap->GetYaxis()->SetRangeUser(0.8,1.2);
  hRatioRap->GetXaxis()->SetRangeUser(RapBinEdges[0], RapBinEdges[nRapBin]);
  hRatioRap->Draw("AP");
  canv4->Update;
  canv4->SaveAs(Form("EfficiencyCrossCheckRap_%dS_%s.png",oniaMode,isPbPb ? "PbPb" : "PP"));


//Close Tfiles
  fEffSantona->Close();
  fEffYongsun->Close();

}

