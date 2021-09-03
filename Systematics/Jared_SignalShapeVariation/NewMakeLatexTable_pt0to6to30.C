void NewMakeLatexTable_pt0to6to30(int whichUpsilon=1, TString dir="ErrorEstimates", TString whichSyst="") {
  
  //TString filename = Form("ErrorEstimates/SysSig%isCombined.root",whichUpsilon);
  //TString filename = Form("ErrorEstimates/SysSig%is_ParamFixingOnly.root",whichUpsilon);
  TString filename = Form("%s/SysSig%is%s.root",dir.Data(),whichUpsilon,whichSyst.Data());
  cout << "%" << filename << endl;
  TFile inFile(filename);
  TH1D* hySysSigPPLowPt = (TH1D*)inFile.Get("hySysSigPPLowPt_in3Sbins;1");
  TH1D* hySysSigPALowPt = (TH1D*)inFile.Get("hySysSigPALowPt_in3Sbins;1");
  TH1D* hySysSigRpALowPt = (TH1D*)inFile.Get("hySysSigRpALowPt_in3Sbins;1");
  TH1D* hySysSigPPHighPt = (TH1D*)inFile.Get("hySysSigPPHighPt_in3Sbins;1");
  TH1D* hySysSigPAHighPt = (TH1D*)inFile.Get("hySysSigPAHighPt_in3Sbins;1");
  TH1D* hySysSigRpAHighPt = (TH1D*)inFile.Get("hySysSigRpAHighPt_in3Sbins;1");

  const int numybins = hySysSigPALowPt->GetSize()-2;
  int binlowlim, binuplim;
  float binlowlimpt, binuplimpt, binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\caption{Systematic uncertainties of $\\Upsilon$(%iS) yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << Form("\\label{sys:signalPDFChange%iS_pt0to6to30}",whichUpsilon) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|cc|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("\\multicolumn{2}{|c|}{Bin}  & \\multicolumn{2}{l|}{%iS Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << " & & pp & pPb &  \\\\" << endl;
cout.precision(2);
  cout << "\\hline" << endl;
  binlowlimpt = 0;
  binuplimpt = 6;
  for (int iy = 1; iy<numybins+1; iy++) {
    binlowlimy = (float)hySysSigPALowPt->GetBinLowEdge(iy);
    binuplimy = (float)hySysSigPALowPt->GetBinLowEdge(iy+1);
    pp1sErr = hySysSigPPLowPt->GetBinContent(iy)*100;
    pPb1sErr = hySysSigPALowPt->GetBinContent(iy)*100;
    RpPbErr = hySysSigRpALowPt->GetBinContent(iy)*100;
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & ",binlowlimpt,binuplimpt,binlowlimy,binuplimy);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
  cout << "\\hline" << endl;
  binlowlimpt = 6;
  binuplimpt = 30;
  for (int iy = 1; iy<numybins+1; iy++) {
    binlowlimy = (float)hySysSigPAHighPt->GetBinLowEdge(iy);
    binuplimy = (float)hySysSigPAHighPt->GetBinLowEdge(iy+1);
    pp1sErr = hySysSigPPHighPt->GetBinContent(iy)*100;
    pPb1sErr = hySysSigPAHighPt->GetBinContent(iy)*100;
    RpPbErr = hySysSigRpAHighPt->GetBinContent(iy)*100;
    cout << Form("$%.0f<\\pt<%.0f \\GeVc$ & $%.2f<y<%.2f$ & ",binlowlimpt,binuplimpt,binlowlimy,binuplimy);
    cout << pp1sErr << " & " << pPb1sErr << " & " << RpPbErr << " \\\\" << endl;
  }
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

  inFile.Close();
}
