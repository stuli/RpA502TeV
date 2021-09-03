void NewMakeLatexTable_HFNtracks(int whichUpsilon=1,int binmode=0, TString dir="ErrorEstimates", TString whichSyst="") {
  //0=hfmode, 1=ntmode

  TString hfntracksbins = "_hfbins.root";
  if (binmode==1) hfntracksbins = "_ntracksbins.root";

  //TString filename = Form("ErrorEstimates/SysSig%isCombined.root",whichUpsilon);
  //TString filename = Form("ErrorEstimates/SysSig%is_ParamFixingOnly.root",whichUpsilon);
  TString filename = Form("%s/SysSig%is%s.root",dir.Data(),whichUpsilon,whichSyst.Data());
  cout << "%" << filename << endl;
  TFile inFile(filename);
  TString histoname = "hHFSysSigRFB000to193;1";
  if (binmode==1) histoname = "hNtracksSysSigRFB000to193;1";
  TH1D* hHFSysSigRFB000to193 = (TH1D*)inFile.Get(histoname);

  //choose a set of bins
  float ybinsCM[3] = {-1.93,0.0,1.93};
  float hfbins3[3] = {0,12,120};
  int ntbins3[3] = {0,40,400};
  float hfbins1[5] = {0,12,19,27,120};
  int ntbins1[5] = {0,40,62,88,400};

  float* hfbinsptr;
  int* ntbinsptr;
  int numhfbinstemp, numntbinstemp;
  if (whichUpsilon==3) {
    hfbinsptr = &hfbins3[0];
    ntbinsptr = &ntbins3[0];
    numhfbinstemp = sizeof(hfbins3)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins3)/sizeof(float)-1;
  }
  else {
    hfbinsptr = &hfbins1[0];
    ntbinsptr = &ntbins1[0];
    numhfbinstemp = sizeof(hfbins1)/sizeof(float)-1;
    numntbinstemp = sizeof(ntbins1)/sizeof(float)-1;
  }

  const int numybins = 2;
  const int numhfbins = numhfbinstemp;
  const int numntbins = numntbinstemp;

  int binlowlim, binuplim, numbins;
  float binlowlimhf, binuplimhf, binlowlimy, binuplimy, pPbFErr, pPbBErr, RFBErr;
  TString strRFBErr = "RFBErr";
  TString binvar, tableName, GeVString;
  if (binmode==0) {
    numbins=numhfbins;
    binvar="\\sum E_T^{HF}";
    tableName = "HF";
    GeVString = "\\GeV";
  }
  else {
    numbins=numntbins;
    binvar="\\textrm{Ntracks}";
    tableName = "Ntracks";
    GeVString = "";
  }

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << Form("\\caption{Systematic uncertainties of forward and backward $\\Upsilon$(%iS) yields and $R_{FB}$ due to signal PDF change to Crystal Ball plus Gaussian, estimated from pseudo-experiments.}",whichUpsilon) << endl;
cout << Form("\\label{sys:signalPDFChange%iS_%s}",whichUpsilon,tableName.Data()) << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|c|}" << endl;
cout << "\\hline" << endl;
cout << Form("Bin  &  $R_{FB}$ Dev.($\\%$)\\\\",whichUpsilon) << endl;
cout << "\\hline" << endl;
cout.precision(2);
for (int ihf = 1; ihf<numbins+1; ihf++) {
  if (binmode==0) {
    binlowlimhf = *(hfbinsptr+ihf-1);
    binuplimhf = *(hfbinsptr+ihf);
  }
  else {
    binlowlimhf = *(ntbinsptr+ihf-1);
    binuplimhf = *(ntbinsptr+ihf);
  }
  binlowlimy = (float)hHFSysSigRFB000to193->GetBinLowEdge(ihf);
  binuplimy = (float)hHFSysSigRFB000to193->GetBinLowEdge(ihf+1);
  //pp1sErr = hHFSysSigRFB000to193->GetBinContent(ihf)*100;
  //pPbBErr = hHFSysSigRFB000to193->GetBinContent(ihf)*100;
  RFBErr = hHFSysSigRFB000to193->GetBinContent(ihf)*100;
  cout << Form("$%.0f<%s<%.0f %s$ & ",binlowlimhf,binvar.Data(),binuplimhf,GeVString.Data());
  cout << RFBErr << " \\\\" << endl;
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\end{table}" << endl << endl;

  inFile.Close();
}
