void SaveDataErrorsInFile() {

//1S
  int numptbins = 6;
  int numybins = 8;
  TFile outFile1("SystematicErrorFromData1s.root", "RECREATE");
  TNtuple* ntuple1spt = new TNtuple("ntuple1spt","Error estimates from data in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numptbins);
  TNtuple* ntuple1sy = new TNtuple("ntuple1sy","Error estimates from data in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr",numybins);

//                    BIN        PP 1S ERR   PA 1S ERR     RpA ERR
 ntuple1spt->Fill(0.00,30.00,       0.768,    0.888542,    0.119623);
  ntuple1spt->Fill(0.00,2.00,  -0.0372428,  -0.0145724,    0.022683);
  ntuple1spt->Fill(2.00,4.00,     1.45666,   -0.324359,    -1.75545);
  ntuple1spt->Fill(4.00,6.00,    0.721426,     1.75972,     1.03085);
  ntuple1spt->Fill(6.00,9.00,    0.365488,     2.97596,     2.60097);
 ntuple1spt->Fill(9.00,12.00,     0.10038,  -0.0292843,   -0.129534);
ntuple1spt->Fill(12.00,30.00,     1.04093,  -0.0820584,    -1.11142);
ntuple1sy->Fill(-1.93,-1.20,    -2.47019,     2.84724,     5.45211);
ntuple1sy->Fill(-1.20,-0.80,      1.7064,    0.108814,    -1.57078);
ntuple1sy->Fill(-0.80,-0.40,    0.183581,    0.440513,    0.256462);
 ntuple1sy->Fill(-0.40,0.00,   0.0282652,   0.0907903,   0.0625089);
  ntuple1sy->Fill(0.00,0.40,   0.0282652,    0.167802,    0.139494);
  ntuple1sy->Fill(0.40,0.80,    0.183581,    0.177548, -0.00602048);
  ntuple1sy->Fill(0.80,1.20,      1.7064,     1.67092,  -0.0348787);
  ntuple1sy->Fill(1.20,1.93,    -2.47019,     5.97424,     8.65831);

  int binlowlim, binuplim;
  float binlowlimy, binuplimy, pp1sErr, pPb1sErr, RpPbErr;

//Print out a table formatted for latex
cout << endl;
cout << "\\begin{table}[hbtp]" << endl;
cout << "\\label{sys:signalPDFChange1SData}" << endl;
cout << "\\centering" << endl;
cout << "\\begin{tabular}{|c|cc|c|}" << endl;
cout << "\\hline" << endl;
cout << "Bin  & \\multicolumn{2}{l}{1S Yield Dev.($\\%$) } & $R_{pA}$ Dev.($\\%$)\\\\" << endl;
cout << "  & pp & pPb &  \\\\" << endl;
cout << "\\hline" << endl;
cout << "$\\pt$, y integrated & 0.008 & 2.36 & 2.35  \\\\" << endl;
cout << "\\hline" << endl;
for (int ipt = 0; ipt<numptbins; ipt++) {
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
  cout << Form("$%i<\\pt<%i \\GeVc$ & %.2f & %.2f & %.2f \\\\",binlowlim,binuplim,pp1sErr,pPb1sErr,RpPbErr) << endl;
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
  cout << Form("$%.2f<\\y<%.2f \\GeVc$ & %.2f & %.2f & %.2f \\\\",binlowlimy,binuplimy,pp1sErr,pPb1sErr,RpPbErr) << endl;
}
cout << "\\hline" << endl;
cout << "\\end{tabular}" << endl;
cout << "\\caption{Percent deviations of 1S yields and $R_{pA}$ due to signal PDF change to Crystal Ball plus Gaussian}" << endl;
cout << "\\end{table}" << endl << endl;

  ntuple1spt->Write();
  ntuple1sy->Write();
  outFile1.Close();

//2S
  TFile outFile2("SystematicErrorFromData2s.root", "RECREATE");
  TNtuple* ntuple2spt = new TNtuple("ntuple2spt","Error estimates from data in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr:ppR21Err:pPbR21Err:DR21Err",3);
  TNtuple* ntuple2sy = new TNtuple("ntuple2sy","Error estimates from data in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr:ppR21Err:pPbR21Err:DR21Err",4);

//    BIN        PP 2S ERR   PA 2S ERR     RpA ERR  PP R21 ERR  PA R21 ERR      DR ERR
 ntuple2spt->Fill(0.00,30.00,    0.477228,     3.24137,     2.75101,    -0.28856,      2.3321,     2.62824);
  ntuple2spt->Fill(0.00,4.00,  0.00294522,    -2.40074,    -2.40361,    0.216678,     -3.3027,    -3.51177);
  ntuple2spt->Fill(4.00,9.00,    0.598469,    0.532294,  -0.0657878,   0.0994367,    0.236349,    0.136779);
 ntuple2spt->Fill(9.00,30.00,     3.72692,   -0.978194,    -4.53606,     6.73083,    -2.48544,    -8.63506);
ntuple2sy->Fill(-1.93,-0.80,    -1.36283,    0.795034,     2.18769,    -3.18302,    0.230041,     3.52527);
 ntuple2sy->Fill(-0.80,0.00,   -0.799227,   -0.438195,     0.36394,    -0.36434,    -1.88021,    -1.52141);
  ntuple2sy->Fill(0.00,0.80,   -0.799227,  -0.0904878,    0.714451,    -0.36434,  -0.0826176,    0.282749);
  ntuple2sy->Fill(0.80,1.93,    -1.36283,   -0.140009,     1.23972,    -3.18302,    -0.23134,     3.04872);

  ntuple2spt->Write();
  ntuple2sy->Write();
  outFile2.Close();

//3S
  TFile outFile3("SystematicErrorFromData3s.root", "RECREATE");
  TNtuple* ntuple3spt = new TNtuple("ntuple3spt","Error estimates from data in pt bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr:ppR21Err:pPbR21Err:DR21Err",2);
  TNtuple* ntuple3sy = new TNtuple("ntuple3sy","Error estimates from data in y bins","binlowlim:binuplim:pp1sErr:pPb1sErr:RpPbErr:ppR21Err:pPbR21Err:DR21Err",2);

//    BIN        PP 3S ERR   PA 3S ERR     RpA ERR  PP R31 ERR  PA R31 ERR      DR ERR
 ntuple3spt->Fill(0.00,30.00,    -1.64717,     4.36482,     6.11268,    -2.39676,     3.44567,     5.98589);
  ntuple3spt->Fill(0.00,6.00,     1.97508,    -10.4746,    -12.2086,     3.31599,    -11.3878,    -14.2319);
 ntuple3spt->Fill(6.00,30.00,   -0.235065,   0.0126428,    0.248294,   -0.463168,    0.010551,    0.475928);
 ntuple3sy->Fill(-1.93,0.00,    -1.64717,  -0.0337579,     1.64043,    -2.39676,  -0.0374606,     2.41723);
  ntuple3sy->Fill(0.00,1.93,    -1.64717,     1.33428,     3.03138,    -2.39676,     1.43605,     3.92692);

  ntuple3spt->Write();
  ntuple3sy->Write();
  outFile3.Close();

}
