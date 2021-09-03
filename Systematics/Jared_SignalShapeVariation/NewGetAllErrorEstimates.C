#include "NewGetErrorEstimates.C"
#include "NewGetErrorEstimatespt0to6to30.C"
#include "NewGetErrorEstimatespt0to6to30_in3Sbins.C"
#include "NewGetErrorEstimatesy193to000to193.C"
#include "NewGetErrorEstimates_HFNtracks.C"
#include "NewGetErrorEstimatesXSBins.C"
#include "NewJaebeomStyle_ChangeNtuplestoHistos.C"
#include "NewMakeAllLatexTables.C"

void NewGetAllErrorEstimates() {

  NewGetErrorEstimates(1);
  NewGetErrorEstimates(2);
  NewGetErrorEstimates(3);
  NewGetErrorEstimatespt0to6to30(1);
  NewGetErrorEstimatespt0to6to30(2);
  NewGetErrorEstimatespt0to6to30(3);
  NewGetErrorEstimatespt0to6to30_in3Sbins(1);
  NewGetErrorEstimatespt0to6to30_in3Sbins(2);
  NewGetErrorEstimatespt0to6to30_in3Sbins(3);
  NewGetErrorEstimatesy193to000to193(1);
  NewGetErrorEstimatesy193to000to193(2);
  NewGetErrorEstimatesy193to000to193(3);
  NewGetErrorEstimates_HFNtracks(1,0);
  NewGetErrorEstimates_HFNtracks(2,0);
  NewGetErrorEstimates_HFNtracks(3,0);
  NewGetErrorEstimates_HFNtracks(1,1);
  NewGetErrorEstimates_HFNtracks(2,1);
  NewGetErrorEstimates_HFNtracks(3,1);
  NewGetErrorEstimatesXSBins(1);
  NewGetErrorEstimatesXSBins(2);
  NewGetErrorEstimatesXSBins(3);

  NewJaebeomStyle_ChangeNtuplestoHistos(1);
  NewJaebeomStyle_ChangeNtuplestoHistos(2);
  NewJaebeomStyle_ChangeNtuplestoHistos(3);

  NewMakeAllLatexTables();

}
