#include "../../HeaderFiles/cutsAndBin.h"

void MakeNominalFiguresNtracksBins(int whichUpsilon=1) {

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
  float yLow, yHigh, yLowB, yHighB;
  int ntLow, ntHigh;
  float hfLow = 0;
  float hfHigh = 400;
  TString caption;

  for (int ihf = 0; ihf<numntbins; ihf++) {

    ntLow = ntbins[ihf];
    ntHigh = ntbins[ihf+1];

    for (int i = 0; i<numybins; i++) {

      yLow = ybins[i];
      yHigh = ybins[i+1];
      yLowB = -yHigh;
      if (yLow>0.0) yHighB = -yLow;
      else yHighB = 0.0;

      caption = Form("	\\caption{Double Crystal Ball Fit to pPb in backward (left) and forward (right) rapdity ranges for Ntracks $\\in$ [%i,%i], $|y|$ $\\in$ [%.2f,%.2f].}",ntLow,ntHigh,yLow,yHigh);

      cout << "\\begin{figure}[hbtp]" << endl;
      cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0_hfsum%.2f-%.2f_ntracks%i-%i}.pdf}", ptLow, ptHigh, yLowB, yHighB, hfLow, hfHigh, ntLow, ntHigh) << endl;
      cout << Form("	\\includegraphics[width=0.45\\textwidth]{{figures/NominalFits/nomfitresults_upsilon_PA_DATA_pt%.1f-%.1f_y%.2f-%.2f_muPt4.0_hfsum%.2f-%.2f_ntracks%i-%i}.pdf}", ptLow, ptHigh, yLow, yHigh, hfLow, hfHigh, ntLow, ntHigh) << endl;
      cout << caption.Data() << endl;
      //cout << Form("	\\label{fig:NomFitnt%ito%iy%.3ito%.3i}", ntLow, ntHigh, (int)(TMath::Abs(yLow)*100+0.5f), (int)(TMath::Abs(yHigh)*100+0.5f)) << endl;
      cout << "\\end{figure}" << endl;
      cout << "\%" << endl;
    }
  }
}
