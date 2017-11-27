void ChangeNtuplestoHistos(int whichUpsilon=2) {
  
  TString filename = Form("SystematicErrorSignal%is.root",whichUpsilon);
  TString outfilename = Form("HistoSystematicErrorSignal%is.root",whichUpsilon);
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
  TString strpp1sErr = Form("pp%isErr",1);
  TString strpPb1sErr = Form("pPb%isErr",1);
  TString strRpPbErr = "RpPbErr";

  //Choose a set of bins
  if (whichUpsilon==1) {
    float ptbins[7] = {0,2,4,6,9,12,30};
    float ybins[9] = {-1.93,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ptbins[4] = {0,4,9,30};
    float ybins[5] = {-1.93,-0.8,0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ptbins[3] = {0,6,30};
    float ybins[3] = {-1.93,0.0,1.93};
  }

  const int numptbins = sizeof(ptbins)/sizeof(float)-1;
  const int numybins = sizeof(ybins)/sizeof(float)-1;

  //Declare histograms
  TH1F* hSignalErrptPP = new TH1F("hSignalErrptPP","hSignalErrptPP",numptbins,ptbins);
  TH1F* hSignalErryPP = new TH1F("hSignalErryPP","hSignalErryPP",numybins,ybins);
  TH1F* hSignalErrptPA = new TH1F("hSignalErrptPA","hSignalErrptPA",numptbins,ptbins);
  TH1F* hSignalErryPA = new TH1F("hSignalErryPA","hSignalErryPA",numybins,ybins);
  TH1F* hSignalErrptRpA = new TH1F("hSignalErrptRpA","hSignalErrptRpA",numptbins,ptbins);
  TH1F* hSignalErryRpA = new TH1F("hSignalErryRpA","hSignalErryRpA",numybins,ybins);

//Fill the histograms
  for (int ipt = 1; ipt<numptbins+1; ipt++) {
    ntuple1spt->GetEntry(ipt);
    TLeaf *pp1sErrLeaf = ntuple1spt->GetLeaf(strpp1sErr);
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *pPb1sErrLeaf = ntuple1spt->GetLeaf(strpPb1sErr);
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1spt->GetLeaf(strRpPbErr);
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    hSignalErrptPP->SetBinContent(ipt, TMath::Abs(pp1sErr)/100);
    hSignalErrptPA->SetBinContent(ipt, TMath::Abs(pPb1sErr)/100);
    hSignalErrptRpA->SetBinContent(ipt, TMath::Abs(RpPbErr)/100);
  }
  for (int iy = 0; iy<numybins; iy++) {
    ntuple1sy->GetEntry(iy);
    TLeaf *pp1sErrLeaf = ntuple1sy->GetLeaf(strpp1sErr);
    pp1sErr = (float)pp1sErrLeaf->GetValue();
    TLeaf *pPb1sErrLeaf = ntuple1sy->GetLeaf(strpPb1sErr);
    pPb1sErr = (float)pPb1sErrLeaf->GetValue();
    TLeaf *RpPbErrLeaf = ntuple1sy->GetLeaf(strRpPbErr);
    RpPbErr = (float)RpPbErrLeaf->GetValue();
    hSignalErryPP->SetBinContent(iy+1, TMath::Abs(pp1sErr)/100);
    hSignalErryPA->SetBinContent(iy+1, TMath::Abs(pPb1sErr)/100);
    hSignalErryRpA->SetBinContent(iy+1, TMath::Abs(RpPbErr)/100);
  }

  //save histograms
  TFile outFile(outfilename, "RECREATE");
  hSignalErrptPP->Write();
  hSignalErrptPA->Write();
  hSignalErrptRpA->Write();
  hSignalErryPP->Write();
  hSignalErryPA->Write();
  hSignalErryRpA->Write();
  outFile.Close();
}
