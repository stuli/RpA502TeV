

void PlotErrors(int whichUpsilon=2) {

  gStyle->SetOptStat(0);
  TString fixedParams[4] = {"alpha","n","x","f"};
  TFile* file[4];
  TH1D* hptSysSigPP[4];
  TCanvas* c1[26];
  TLegend* leg[26];
  TString HistoNames[26] = {"hintSysSigPP;1","hintSysSigPA;1","hintSysSigRpA;1","hptSysSigPP;1", "hptSysSigPA;1", "hptSysSigRpA;1", "hySysSigPP;1", "hySysSigPA;1", "hySysSigRpA;1", "hptSysSigPPBackwardY;1", "hptSysSigPABackwardY;1", "hptSysSigRpABackwardY;1", "hptSysSigPPForwardY;1", "hptSysSigPAForwardY;1", "hptSysSigRpAForwardY;1", "hySysSigPPLowPt;1", "hySysSigPALowPt;1", "hySysSigRpALowPt;1", "hySysSigPPHighPt;1", "hySysSigPAHighPt;1", "hySysSigRpAHighPt;1", "hptSysSigPA_y287to193;1", "hrapSysSigCross;1", "hHFSysSigRFB000to193;1", "hNtracksSysSigRFB000to193;1", "hSysSigRFBIntActivity;1"};

for (int j=0; j<26; j++) {
  cout << "here1" << endl;
  c1[j] = new TCanvas(Form("c1[%i]",j),"c1",4,40,400,400);
  leg[j] = new TLegend(0.7,0.7,0.9,0.9);
  cout << "here2" << endl;
  for (int i=0; i<4; i++) {
    cout << "reading file " << i << endl;
    cout << Form("SysSig%is_%s.root",whichUpsilon,fixedParams[i].Data()) << endl;
    file[i] = TFile::Open(Form("SysSig%is_%s.root",whichUpsilon,fixedParams[i].Data()),"READ");
    hptSysSigPP[i] = (TH1D*)file[i]->Get(HistoNames[j]);
    hptSysSigPP[i]->Draw("same");
    hptSysSigPP[i]->SetLineColor(i+1);
    hptSysSigPP[i]->SetMinimum(0);
    leg[j]->AddEntry(hptSysSigPP[i],fixedParams[i].Data(),"l");
  }
  leg[j]->Draw("same");
}


  c1[0]->Print(Form("MChistograms%is.pdf(",whichUpsilon),"pdf");
  for (int k=1; k<25; k++) {
    c1[k]->Print(Form("MChistograms%is.pdf",whichUpsilon),"pdf");
  }
  c1[25]->Print(Form("MChistograms%is.pdf)",whichUpsilon),"pdf");

}
