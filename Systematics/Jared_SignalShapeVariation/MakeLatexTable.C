void MakeLatexTable(int whichUpsilon=1) {
  
  TString filename = Form("SystematicErrorFromData%is.root",whichUpsilon);
  cout << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries();
  TString ntupleyname = Form("ntuple%isy;1",whichUpsilon);
  TNtuple* ntuple1sy = (TNtuple*)inFile->Get(ntupleyname);
  const int numybins = (int)ntuple1sy->GetEntries();

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iSData}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("Bin  & \\multicolumn{2}{l}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << "  & pp & pPb &  \\\\" << endl;
cout << "\\hline" << endl;
  ntuple1spt->GetEntry(0);
  TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf("pp1sErr");
  pp1sErr = (float)pp1sErrLeaf->GetValue();
  TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf("pPb1sErr");
  pPb1sErr = (float)pPb1sErrLeaf->GetValue();
  TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf("RpPbErr");
  RpPbErr = (float)RpPbErrLeaf->GetValue();
cout << Form("$\\pt$, y integrated & %.2f & %.2f & %.2f  \\\\",pp1sErr,pPb1sErr,RpPbErr) << endl;
cout << "\\hline" << endl;
for (int ipt = 1; ipt<numptbins; ipt++) {
  ntuple1spt->GetEntry(ipt);
  TLeaf *binlowlimLeaf = ntuple1spt->GetLeaf("binlowlim");
  binlowlim = (int)binlowlimLeaf->GetValue();
  TLeaf *binuplimLeaf = ntuple1spt->GetLeaf("binuplim");
  binuplim = (int)binuplimLeaf->GetValue();
  TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf("pp1sErr");
  pp1sErr = (float)pp1sErrLeaf->GetValue();
  TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf("pPb1sErr");
  pPb1sErr = (float)pPb1sErrLeaf->GetValue();
  TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf("RpPbErr");
  RpPbErr = (float)RpPbErrLeaf->GetValue();
  if (binlowlim>0)  cout << Form("$%i<\\pt<%i \\GeVc$ & %.2f & %.2f & %.2f \\\\",binlowlim,binuplim,pp1sErr,pPb1sErr,RpPbErr) << endl;
  else  cout << Form("$\\pt<%i \\GeVc$ & %.2f & %.2f & %.2f \\\\",binuplim,pp1sErr,pPb1sErr,RpPbErr) << endl;
}
cout << "\\hline" << endl;
for (int iy = 0; iy<numybins; iy++) {
  ntuple1sy->GetEntry(iy);
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
  cout << Form("$%.2f<\\y<%.2f$ & %.2f & %.2f & %.2f \\\\",binlowlimy,binuplimy,pp1sErr,pPb1sErr,RpPbErr) << endl;
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << Form("\\caption{Percent deviations of %iS yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian}",whichUpsilon) << endl;
cout << "\\end{table}" << endl << endl;

}
