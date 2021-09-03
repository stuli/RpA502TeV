void MakeLatexTable_pt0to6to30(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is_pt0to6to30.root",whichUpsilon);
  cout << "%" << filename << endl;
  TFile *inFile = new TFile(filename);

  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);
  const int numybins = (int)ntuple1sy->GetEntries()/2;

  int binlowlim, binuplim;
  float binlowlimpt, binuplimpt, binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iS_pt0to6to30}",whichUpsilon) << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|cc|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("\\multicolumn{2}{|c|}{Bin}  & \\multicolumn{2}{l|}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << " & & pp & pPb &  \\\\" << endl;
for (int ipt = 0; ipt<2; ipt++) {
  cout << "\\hline" << endl;
  if (ipt==0) {
    binlowlimpt = 0;
    binuplimpt = 6;
  }
  else {
    binlowlimpt = 6;
    binuplimpt = 30;
  }
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1sy->GetEntry(iy+ipt*numybins);
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
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & %.2f & %.2f & %.2f \\\\",binlowlimpt,binuplimpt,binlowlimy,binuplimy,pp1sErr,pPb1sErr,RpPbErr) << endl;
  }
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

}
