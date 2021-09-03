#include "GetErrorEstimates.C"
#include "GetErrorEstimatespt0to6to30.C"
#include "GetErrorEstimatesy193to000to193.C"
#include "GetErrorEstimates_HFNtracks.C"
#include "GetErrorEstimates_HFNtracksIntBin.C"
#include "GetErrorEstimatesXSBins.C"
#include "ChangeNtuplestoHistos.C"
#include "ChangeNtuplestoHistos_y287to193.C"
#include "ChangeNtuplestoHistos_y240to193.C"
#include "ChangeNtuplestoHistospt0to6to30.C"
#include "ChangeNtuplestoHistosy193to000to193.C"
#include "ChangeNtuplestoHistos_HFNtracks.C"
#include "ChangeNtuplestoHistosXSBins.C"
#include "MakeAllLatexTables.C"

void GetAllErrorEstimates() {

  GetErrorEstimates(1);
  GetErrorEstimates(2);
  GetErrorEstimates(3);
  GetErrorEstimatespt0to6to30(1);
  GetErrorEstimatespt0to6to30(2);
  GetErrorEstimatespt0to6to30(3);
  GetErrorEstimatesy193to000to193(1);
  GetErrorEstimatesy193to000to193(2);
  GetErrorEstimatesy193to000to193(3);
  GetErrorEstimates_HFNtracks(1,0);
  GetErrorEstimates_HFNtracks(2,0);
  GetErrorEstimates_HFNtracks(3,0);
  GetErrorEstimates_HFNtracks(1,1);
  GetErrorEstimates_HFNtracks(2,1);
  GetErrorEstimates_HFNtracks(3,1);
  GetErrorEstimates_HFNtracksIntBin(1);
  GetErrorEstimates_HFNtracksIntBin(2);
  GetErrorEstimates_HFNtracksIntBin(3);
  GetErrorEstimatesXSBins(1);
  GetErrorEstimatesXSBins(2);
  GetErrorEstimatesXSBins(3);

  ChangeNtuplestoHistos(1);
  ChangeNtuplestoHistos(2);
  ChangeNtuplestoHistos(3);
  ChangeNtuplestoHistos_y287to193(1);
  ChangeNtuplestoHistos_y287to193(2);
  ChangeNtuplestoHistos_y287to193(3);
  ChangeNtuplestoHistos_y240to193(1);
  ChangeNtuplestoHistos_y240to193(2);
  ChangeNtuplestoHistos_y240to193(3);
  ChangeNtuplestoHistospt0to6to30(1);
  ChangeNtuplestoHistospt0to6to30(2);
  ChangeNtuplestoHistospt0to6to30(3);
  ChangeNtuplestoHistosy193to000to193(1);
  ChangeNtuplestoHistosy193to000to193(2);
  ChangeNtuplestoHistosy193to000to193(3);
  ChangeNtuplestoHistos_HFNtracks(1,0);
  ChangeNtuplestoHistos_HFNtracks(2,0);
  ChangeNtuplestoHistos_HFNtracks(3,0);
  ChangeNtuplestoHistos_HFNtracks(1,1);
  ChangeNtuplestoHistos_HFNtracks(2,1);
  ChangeNtuplestoHistos_HFNtracks(3,1);
  ChangeNtuplestoHistosXSBins(1);
  ChangeNtuplestoHistosXSBins(2);
  ChangeNtuplestoHistosXSBins(3);

  MakeAllLatexTables();

}
