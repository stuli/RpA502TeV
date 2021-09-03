void NewMakeLatexTable(int whichUpsilon=1, TString dir="ErrorEstimates", TString whichSyst="") {
  
  //TString filename = Form("ErrorEstimates/SysSig%isCombined.root",whichUpsilon);
  //TString filename = Form("ErrorEstimates/SysSig%is_ParamFixingOnly.root",whichUpsilon);
  TString filename = Form("%s/SysSig%is%s.root",dir.Data(),whichUpsilon,whichSyst.Data());
  //TString filename = "ErrorEstimates/SysSig1s_alpha.root";
  cout << "%" << filename << endl;
  TFile inFile(filename);
  TH1D* hintSysSigPP = (TH1D*)inFile.Get("hintSysSigPP;1");
  TH1D* hintSysSigPA = (TH1D*)inFile.Get("hintSysSigPA;1");
  TH1D* hintSysSigRpA = (TH1D*)inFile.Get("hintSysSigRpA;1");
  TH1D* hptSysSigPP = (TH1D*)inFile.Get("hptSysSigPP;1");
  TH1D* hptSysSigPA = (TH1D*)inFile.Get("hptSysSigPA;1");
  TH1D* hptSysSigRpA = (TH1D*)inFile.Get("hptSysSigRpA;1");
  TH1D* hySysSigPP = (TH1D*)inFile.Get("hySysSigPP;1");
  TH1D* hySysSigPA = (TH1D*)inFile.Get("hySysSigPA;1");
  TH1D* hySysSigRpA = (TH1D*)inFile.Get("hySysSigRpA;1");

  const int numptbins = hptSysSigPA->GetSize()-1;
  const int numybins = hySysSigPA->GetSize()-2;
  int binlowlim, binuplim;
  double binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << Form("\\label{sys:signalPDFChange%iS}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("Bin  & \\multicolumn{2}{|l|}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << "  & pp & pPb &  \\\\" << endl;
cout << "\\hline" << endl;
  pp1sErr = hintSysSigPP->GetBinContent(1)*100;
  pPb1sErr = hintSysSigPA->GetBinContent(1)*100;
  RpPbErr = hintSysSigRpA->GetBinContent(1)*100;
cout << Form("$\\pt$, $y$ integrated & ");
cout.precision(2);
cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
cout << "\\hline" << endl;
for (int ipt = 1; ipt<numptbins; ipt++) {
  binlowlim = (int)hptSysSigRpA->GetBinLowEdge(ipt);
  binuplim = (int)hptSysSigRpA->GetBinLowEdge(ipt+1);
  pp1sErr = hptSysSigPP->GetBinContent(ipt)*100;
  pPb1sErr = hptSysSigPA->GetBinContent(ipt)*100;
  RpPbErr = hptSysSigRpA->GetBinContent(ipt)*100;
  if (binlowlim>0) {
    cout << Form("$%i<\\pt<%i \\GeVc$ & ",binlowlim,binuplim);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
  else {
    cout << Form("$\\pt<%i \\GeVc$ & ",binuplim);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
}
cout << "\\hline" << endl;
for (int iy = 1; iy<numybins+1; iy++) {
  binlowlimy = (float)hySysSigRpA->GetBinLowEdge(iy);
  binuplimy = (float)hySysSigRpA->GetBinLowEdge(iy+1);
  pp1sErr = hySysSigPP->GetBinContent(iy)*100;
  pPb1sErr = hySysSigPA->GetBinContent(iy)*100;
  RpPbErr = hySysSigRpA->GetBinContent(iy)*100;
  cout << Form("$%.2f<y<%.2f$ & ",binlowlimy,binuplimy);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

  inFile.Close();
}
