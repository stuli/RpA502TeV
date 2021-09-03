#include "PlotFittedParams.C"
#include "PlotFittedParams_pt0to6to30.C"
#include "PlotFittedParams_y193to0to193.C"
#include "PlotFittedParams_HF.C"
#include "PlotFittedParams_Ntracks.C"

void PlotFittedParamsAll() {

  //Regular pt and rapidity bins.
  //PlotFittedParams(1);
  //PlotFittedParams(2);
  //PlotFittedParams(3);
  //PlotFittedParams(1,kPPDATA);
  //PlotFittedParams(2,kPPDATA);
  //PlotFittedParams(3,kPPDATA);

  //Pt 0-6-30 bins
  //PlotFittedParams_pt0to6to30(1);
  //PlotFittedParams_pt0to6to30(2);
  //PlotFittedParams_pt0to6to30(3);
  //PlotFittedParams_pt0to6to30(1,kPPDATA);
  //PlotFittedParams_pt0to6to30(2,kPPDATA);
  //PlotFittedParams_pt0to6to30(3,kPPDATA);

  //y-193-0.0-1.93 bins
  //PlotFittedParams_y193to0to193(1);
  //PlotFittedParams_y193to0to193(2);
  //PlotFittedParams_y193to0to193(3);
  //PlotFittedParams_y193to0to193(1,kPPDATA);
  //PlotFittedParams_y193to0to193(2,kPPDATA);
  //PlotFittedParams_y193to0to193(3,kPPDATA);

  //HF bins
  //PlotFittedParams_HF(1,0);
  //PlotFittedParams_HF(1,1);
  //PlotFittedParams_HF(1,2);
  //PlotFittedParams_HF(1,3);
  //PlotFittedParams_HF(2,0);
  //PlotFittedParams_HF(2,1);
  //PlotFittedParams_HF(3,0);

  //Ntracks bins
  //PlotFittedParams_Ntracks(1,0);
  //PlotFittedParams_Ntracks(1,1);
  //PlotFittedParams_Ntracks(1,2);
  //PlotFittedParams_Ntracks(1,3);
  //PlotFittedParams_Ntracks(2,0);
  PlotFittedParams_Ntracks(2,1);
  PlotFittedParams_Ntracks(3,0);

}
