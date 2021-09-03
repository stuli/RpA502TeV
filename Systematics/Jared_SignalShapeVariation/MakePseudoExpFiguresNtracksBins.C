#include "../../HeaderFiles/cutsAndBin.h"

void MakePseudoExpFiguresNtracksBins(int whichUpsilon=1) {

  //choose a set of bins
  if (whichUpsilon==1) {
    float ybins[5] = {0.0,0.4,0.8,1.2,1.93};
  }
  else if (whichUpsilon==2) {
    float ybins[3] = {0.0,0.8,1.93};
  }
  else if (whichUpsilon==3) {
    float ybins[2] = {0.0,1.93};
  }
  float ntbins[5] = {0,40,65,90,400};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numntbins = sizeof(ntbins)/sizeof(float)-1;

  float ptLow = 0;
  float ptHigh = 30;
  float yLow, yHigh;
  float hfLow = 0;
  float hfHigh = 400;
  int ntLow, ntHigh;
  TString caption;

  for (int ihf = 0; ihf<numntbins; ihf++) {

    ntLow = ntbins[ihf];
    ntHigh = ntbins[ihf+1];

    for (int i = 0; i<numybins; i++) {

      yLow = ybins[i];
      yHigh = ybins[i+1];

      caption = Form("	\\caption{Percent differences in $\\Upsilon$(%iS) yields between the nominal and alternative signal PDF fits to pseudodata generated in the bin: Ntracks $\\in$ [%i,%i], y $\\in$ [%.2f,%.2f].}",whichUpsilon,ntLow,ntHigh,yLow,yHigh);

      cout << "\\begin{figure}[hbtp]" << endl;
      cout << Form("	\\includegraphics[width=1.0\\textwidth]{{figures/systematics_SignalPDF/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0_hfsum%.2f-%.2f_ntracks%i-%i}.pdf}", whichUpsilon, ptLow, ptHigh, yLow, yHigh, hfLow, hfHigh, ntLow, ntHigh) << endl;
      cout << caption.Data() << endl;
      //cout << Form("	\\label{fig:PseudoExp%iSnt%ito%iy%.3ito%.3i}", whichUpsilon, ntLow, ntHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
      cout << "\\end{figure}" << endl;
      cout << "\%" << endl;
    }
  }
}
