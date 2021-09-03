#include "PlotFittedParamsSinglePlots.C"
#include "PlotFittedParamsXSSinglePlots.C"
#include "PlotFittedParams_pt0to6to30SinglePlots.C"
#include "PlotFittedParams_y193to0to193SinglePlots.C"
#include "PlotFittedParams_HFSinglePlots.C"
#include "PlotFittedParams_NtracksSinglePlots.C"
#include "PlotFittedParams_intBinSinglePlots.C"

void PlotFittedParamsAllSinglePlots() {

  //Integrated bin
  //PlotFittedParams_intBinSinglePlots();
  //PlotFittedParams_intBinSinglePlots(kPPDATA);

  //Regular pt and rapidity bins.
  //PlotFittedParamsSinglePlots(1);
  //PlotFittedParamsSinglePlots(2);
  //PlotFittedParamsSinglePlots(3);
  //PlotFittedParamsSinglePlots(1,kPPDATA);
  //PlotFittedParamsSinglePlots(2,kPPDATA);
  //PlotFittedParamsSinglePlots(3,kPPDATA);

  //pt bins in -2.87<y<1.93
  //PlotFittedParamsXSSinglePlots(1);
  //PlotFittedParamsXSSinglePlots(2);
  //PlotFittedParamsXSSinglePlots(3);

  //Pt 0-6-30 bins
  //PlotFittedParams_pt0to6to30SinglePlots(1);
  //PlotFittedParams_pt0to6to30SinglePlots(2);
  //PlotFittedParams_pt0to6to30SinglePlots(3);
  //PlotFittedParams_pt0to6to30SinglePlots(1,kPPDATA);
  //PlotFittedParams_pt0to6to30SinglePlots(2,kPPDATA);
  //PlotFittedParams_pt0to6to30SinglePlots(3,kPPDATA);

  //y-193-0.0-1.93 bins
  //PlotFittedParams_y193to0to193SinglePlots(1);
  //PlotFittedParams_y193to0to193SinglePlots(2);
  //PlotFittedParams_y193to0to193SinglePlots(3);
  //These ones are just the same as the regular pt bins:
  //PlotFittedParams_y193to0to193SinglePlots(1,kPPDATA);
  //PlotFittedParams_y193to0to193SinglePlots(2,kPPDATA);
  //PlotFittedParams_y193to0to193SinglePlots(3,kPPDATA);

  //HF bins
  //PlotFittedParams_HFSinglePlots(1);
  //PlotFittedParams_HFSinglePlots(2);
  //PlotFittedParams_HFSinglePlots(3);

  //Ntracks bins
  //PlotFittedParams_NtracksSinglePlots(1);
  //PlotFittedParams_NtracksSinglePlots(2);
  //PlotFittedParams_NtracksSinglePlots(3);

  gROOT->ProcessLine(".q");
}
