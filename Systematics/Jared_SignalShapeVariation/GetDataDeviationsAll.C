#include "GetDataDeviations.C"
#include "GetDataDeviationspt0to6to30.C"
#include "GetDataDeviationsy193to000to193.C"
#include "GetDataDeviations_HFNtracks.C"
#include "GetDataDeviationsXSBins.C"
#include "NewJaebeomStyle_fromData.C"

void GetDataDeviationsAll() {
/*
  GetDataDeviations(1);
  GetDataDeviations(2);
  GetDataDeviations(3);

  GetDataDeviationspt0to6to30(1);
  GetDataDeviationspt0to6to30(2);
  GetDataDeviationspt0to6to30(3);

  GetDataDeviationsy193to000to193(1);
  GetDataDeviationsy193to000to193(2);
  GetDataDeviationsy193to000to193(3);
*/
  GetDataDeviationsXSBins(1);
  GetDataDeviationsXSBins(2);
  GetDataDeviationsXSBins(3);

  GetDataDeviations_HFNtracks(1,0);
  GetDataDeviations_HFNtracks(2,0);
  GetDataDeviations_HFNtracks(3,0);

  GetDataDeviations_HFNtracks(1,1);
  GetDataDeviations_HFNtracks(2,1);
  GetDataDeviations_HFNtracks(3,1);

  NewJaebeomStyle_fromData(1);
  NewJaebeomStyle_fromData(2);
  NewJaebeomStyle_fromData(3);
}
