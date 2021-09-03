void MakeLatexTable_y287to193(int whichUpsilon=1) {
  
  TString filename = Form("ErrorEstimates/SystematicErrorSignal%isXSBins.root",whichUpsilon);
  cout << "%" << filename << endl;
  TFile *inFile = new TFile(filename);
  TString ntupleptname = Form("ntuple%ispt;1",whichUpsilon);
  TNtuple* ntuple1spt = (TNtuple*)inFile->Get(ntupleptname);
  const int numptbins = (int)ntuple1spt->GetEntries();

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\label{sys:signalPDFChange%iSXSBins}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("Bin  & pPb %iS Yield Dev.($\\%$) \\\\",whichUpsilon) << endl;
cout << "\\hline" << endl;
  ntuple1spt->GetEntry(0);
  TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf("pPb1sErr");
  pPb1sErr = (float)pPb1sErrLeaf->GetValue();
cout << Form("$\\pt$ integrated & %.2f \\\\",pPb1sErr) << endl;
cout << "\\hline" << endl;
for (int ipt = 1; ipt<numptbins; ipt++) {
  ntuple1spt->GetEntry(ipt);
  TLeaf *binlowlimLeaf = ntuple1spt->GetLeaf("binlowlim");
  binlowlim = (int)binlowlimLeaf->GetValue();
  TLeaf *binuplimLeaf = ntuple1spt->GetLeaf("binuplim");
  binuplim = (int)binuplimLeaf->GetValue();
  TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf("pPb1sErr");
  pPb1sErr = (float)pPb1sErrLeaf->GetValue();
  if (binlowlim>0)  cout << Form("$%i<\\pt<%i \\GeVc$ & %.2f \\\\",binlowlim,binuplim,pPb1sErr) << endl;
  else  cout << Form("$\\pt<%i \\GeVc$ & %.2f\\\\",binuplim,pPb1sErr) << endl;
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yield due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << "\\end{table}" << endl << endl;

}
