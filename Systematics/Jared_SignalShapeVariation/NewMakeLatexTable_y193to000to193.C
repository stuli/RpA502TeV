void NewMakeLatexTable_y193to000to193(int whichUpsilon=1, TString dir="ErrorEstimates", TString whichSyst="") {

  //TString filename = Form("ErrorEstimates/SysSig%isCombined.root",whichUpsilon);
  //TString filename = Form("ErrorEstimates/SysSig%is_ParamFixingOnly.root",whichUpsilon);
  TString filename = Form("%s/SysSig%is%s.root",dir.Data(),whichUpsilon,whichSyst.Data());
  cout << "%" << filename << endl;
  TFile inFile(filename);
  TH1D* hptSysSigPPBackwardY = (TH1D*)inFile.Get("hptSysSigPPBackwardY;1");
  TH1D* hptSysSigPABackwardY = (TH1D*)inFile.Get("hptSysSigPABackwardY;1");
  TH1D* hptSysSigRpABackwardY = (TH1D*)inFile.Get("hptSysSigRpABackwardY;1");
  TH1D* hptSysSigPPForwardY = (TH1D*)inFile.Get("hptSysSigPPForwardY;1");
  TH1D* hptSysSigPAForwardY = (TH1D*)inFile.Get("hptSysSigPAForwardY;1");
  TH1D* hptSysSigRpAForwardY = (TH1D*)inFile.Get("hptSysSigRpAForwardY;1");

  const int numptbins = hptSysSigPABackwardY->GetSize()-2;
  int binlowlim, binuplim;
  float binlowlimpt, binuplimpt, binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << Form("\\label{sys:signalPDFChange%iS_y193to0to193}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|cc|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("\\multicolumn{2}{|c|}{Bin}  & \\multicolumn{2}{l|}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << " & & pp & pPb &  \\\\" << endl;
cout.precision(2);
  cout << "\\hline" << endl;
  binlowlimy = -1.93;
  binuplimy = 0.0;
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    binlowlimpt = (int)hptSysSigPABackwardY->GetBinLowEdge(ipt);
    binuplimpt = (int)hptSysSigPABackwardY->GetBinLowEdge(ipt+1);
    pp1sErr = hptSysSigPPBackwardY->GetBinContent(ipt)*100;
    pPb1sErr = hptSysSigPABackwardY->GetBinContent(ipt)*100;
    RpPbErr = hptSysSigRpABackwardY->GetBinContent(ipt)*100;
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & ",binlowlimpt,binuplimpt,binlowlimy,binuplimy);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
  binlowlimy = 0.0;
  binuplimy = 1.93;
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    binlowlimpt = (int)hptSysSigPAForwardY->GetBinLowEdge(ipt);
    binuplimpt = (int)hptSysSigPAForwardY->GetBinLowEdge(ipt+1);
    pp1sErr = hptSysSigPPForwardY->GetBinContent(ipt)*100;
    pPb1sErr = hptSysSigPAForwardY->GetBinContent(ipt)*100;
    RpPbErr = hptSysSigRpAForwardY->GetBinContent(ipt)*100;
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & ",binlowlimpt,binuplimpt,binlowlimy,binuplimy);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

  inFile.Close();
}
