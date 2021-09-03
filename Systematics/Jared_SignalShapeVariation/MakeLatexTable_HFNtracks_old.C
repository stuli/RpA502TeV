#include "TString.h"
#include "TNtuple.h"

void MakeLatexTable_HFNtracks(int whichUpsilon=1,int binmode=0) {
  //0=hfmode, 1=ntmode

  TString hfntracksbins = "_hfbins.root";
  if (binmode==1) hfntracksbins = "_ntracksbins.root";

  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is",whichUpsilon);
  filename = filename + hfntracksbins;
  cout << "%" << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntuplename = Form("ntuple%is;1",whichUpsilon);
  TNtuple* ntuple1s = (TNtuple*)inFile->Get(ntuplename);

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[5] = {0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[3] = {0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[2] = {0.0,1.93};
  }
  float hfbins[5] = {0,20,30,40,120};
  int ntracksbins[6] = {0,20,30,40,120,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntracksbins = sizeof(ntracksbins)/sizeof(float)-1;

  int binlowlim, binuplim;
  float binlowlimhf, binuplimhf, binlowlimy, binuplimy, pPbFErr, pPbBErr, RFBErr;
  TString strRFBErr = "RFBErr";
  TString binvar, tableName, GeVString;
  if (binmode==0) {
    const int numbins=numhfbins;
    binvar="\\sum E_T^{HF}";
    tableName = "HF";
    GeVString = "\\GeV";
  }
  else {
    const int numbins=numntracksbins;
    binvar="\\textrm{Ntracks}";
    tableName = "Ntracks";
    GeVString = "";
  }

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iS_%s}",whichUpsilon,tableName.Data()) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|cc|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("\\multicolumn{2}{|c|}{Bin}  & \\multicolumn{2}{l|}{%iS Yield Dev.($\\%$) } & $R_{FB}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << " & & Forward & Backward &  \\\\" << endl;
for (int ihf = 0; ihf<numbins; ihf++) {
  cout << "\\hline" << endl;
  if (binmode==0) {
    binlowlimhf = hfbins[ihf];
    binuplimhf = hfbins[ihf+1];
  }
  else {
    binlowlimhf = ntracksbins[ihf];
    binuplimhf = ntracksbins[ihf+1];
  }
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1s->GetEntry(ihf*numybins+iy);
    TLeaf *binlowlimLeaf = ntuple1s->GetLeaf("binlowlim");
    binlowlimy = (float)binlowlimLeaf->GetValue();
    TLeaf *binuplimLeaf = ntuple1s->GetLeaf("binuplim");
    binuplimy = (float)binuplimLeaf->GetValue();
    TLeaf *pPbFErrLeaf = ntuple1s->GetLeaf(Form("pPb%isFerr",whichUpsilon));
    pPbFErr = (float)pPbFErrLeaf->GetValue();
    TLeaf *pPbBErrLeaf = ntuple1s->GetLeaf(Form("pPb%isBerr",whichUpsilon));
    pPbBErr = (float)pPbBErrLeaf->GetValue();
    TLeaf *RFBErrLeaf = ntuple1s->GetLeaf("RFBErr");
    RFBErr = (float)RFBErrLeaf->GetValue();
    cout << Form("$%.0f<%s<%.0f %s$ & $%.2f<y<%.2f$ & %.2f & %.2f & %.2f \\\\",binlowlimhf,binvar.Data(),binuplimhf,GeVString.Data(),binlowlimy,binuplimy,pPbFErr,pPbBErr,RFBErr) << endl;
  }
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << Form("\\caption{Systematic uncertainties of forward and backward $\\Upsilon$(%iS) yields and $R_{FB}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << "\\end{table}" << endl << endl;

}
