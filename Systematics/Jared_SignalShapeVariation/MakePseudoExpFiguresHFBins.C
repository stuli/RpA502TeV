#include "../../HeaderFiles/cutsAndBin.h"

void MakePseudoExpFiguresHFBins(int whichUpsilon=1) {

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
  float hfbins[5] = {0,15,22,30,120};

  const int numybins = sizeof(ybins)/sizeof(float)-1;
  const int numhfbins = sizeof(hfbins)/sizeof(float)-1;

  float ptLow = 0;
  float ptHigh = 30;
  float yLow, yHigh, hfLow, hfHigh;
  int ntLow = 0;
  int ntHigh = 400;
  TString caption;

  for (int ihf = 0; ihf<numhfbins; ihf++) {

    hfLow = hfbins[ihf];
    hfHigh = hfbins[ihf+1];

    for (int i = 0; i<numybins; i++) {

      yLow = ybins[i];
      yHigh = ybins[i+1];

      caption = Form("	\\caption{Percent differences in $\\Upsilon$(%iS) yields between the nominal and alternative signal PDF fits to pseudodata generated in the bin: $E_T^{HF|\eta|>4}\\in$ [%.0f,%.0f], y $\\in$ [%.2f,%.2f].}",whichUpsilon,hfLow,hfHigh,yLow,yHigh);

      cout << "\\begin{figure}[hbtp]" << endl;
      cout << Form("	\\includegraphics[width=1.0\\textwidth]{{figures/systematics_SignalPDF/PseudoExpPlots%iS_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0_hfsum%.2f-%.2f_ntracks%i-%i}.pdf}", whichUpsilon, ptLow, ptHigh, yLow, yHigh, hfLow, hfHigh, ntLow, ntHigh) << endl;
      cout << caption.Data() << endl;
      //cout << Form("	\\label{fig:PseudoExp%iShf%.0fto%.0fy%.3ito%.3i}", whichUpsilon, hfLow, hfHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
      cout << "\\end{figure}" << endl;
      cout << "\%" << endl;
    }
  }
}
