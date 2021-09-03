#include "TString.h"
#include "TNtuple.h"

void MakeLatexTable_HFNtracks(int whichUpsilon=1,int binmode=0) {
  //0=hfmode, 1=ntmode

  TString hfntracksbins = "_hfbins.root";
  if (binmode==1) hfntracksbins = "_ntracksbins.root";

  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is",whichUpsilon);
  filename = filename + hfntracksbins;
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[3] = {-1.93,0.0,1.93};
  }
  float hfbins[5] = {0,20,30,40,120};
  int ntracksbins[6] = {0,20,30,40,120,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;
  const int numntracksbins = sizeof(ntracksbins)/sizeof(float)-1;

  int binlowlim, binuplim;
  float binlowlimhf, binuplimhf, binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iS}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("Bin  & \\multicolumn{1}{l}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << "  & pp & pPb &  \\\\" << endl;
for (int ihf = 0; ihf<numhfbins; ihf++) {
  cout << "\\hline" << endl;
  binlowlimhf = hfbins[ihf];
  binuplimhf = hfbins[ihf+1];
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1sy->GetEntry(ihf*numybins+iy);
    TLeaf *binlowlimLeaf = ntuple1sy->GetLeaf("binlowlim");
    binlowlimy = (float)binlowlimLeaf->GetValue();
    TLeaf *binuplimLeaf = ntuple1sy->GetLeaf("binuplim");
    binuplimy = (float)binuplimLeaf->GetValue();
    TLeaf *pp1sErrLeaf = ntuple1sy->GetLeaf("pp1sErr");
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf("pPb1sErr");
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1sy->GetLeaf("RpPbErr");
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    cout << Form("$%.2f<hf<%.2f,%.2f<y<%.2f$ & %.2f & %.2f & %.2f \\\\",binlowlimhf,binuplimhf,binlowlimy,binuplimy,pp1sErr,pPb1sErr,RpPbErr) << endl;
  }
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << Form("\\caption{Systematic uncertainties of %iS yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << "\\end{table}" << endl << endl;

}
