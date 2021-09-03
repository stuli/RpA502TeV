void MakeLatexTable_y193to000to193(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%is_y193to000to193.root",whichUpsilon);
  cout << "%" << filename << endl;
  TFile *inFile = new TFile(filename);

  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries()/2;

  int binlowlim, binuplim;
  float binlowlimpt, binuplimpt, binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iS_y193to0to193}",whichUpsilon) << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|cc|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("\\multicolumn{2}{|c|}{Bin}  & \\multicolumn{2}{l|}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << " & & pp & pPb &  \\\\" << endl;
for (int iy = 0; iy<2; iy++) {
  cout << "\\hline" << endl;
  if (iy==0) {
    binlowlimy = -1.93;
    binuplimy = 0.0;
  }
  else {
    binlowlimy = 0.0;
    binuplimy = 1.93;
  }
  for (int ipt = 0; ipt<numptbins; ipt++) {
    ntuple1spt->GetEntry(ipt+iy*numptbins);
    TLeaf *binlowlimLeaf = ntuple1spt->GetLeaf("binlowlim");
    binlowlimpt = (float)binlowlimLeaf->GetValue();
    TLeaf *binuplimLeaf = ntuple1spt->GetLeaf("binuplim");
    binuplimpt = (float)binuplimLeaf->GetValue();
    TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf("pp1sErr");
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf("pPb1sErr");
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf("RpPbErr");
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & %.2f & %.2f & %.2f \\\\",binlowlimpt,binuplimpt,binlowlimy,binuplimy,pp1sErr,pPb1sErr,RpPbErr) << endl;
  }
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

}
