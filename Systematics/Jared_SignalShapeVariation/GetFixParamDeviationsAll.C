#include "GetFixParamDeviations.C"
#include "GetFixParamDeviationspt0to6to30.C"
#include "GetFixParamDeviationspt0to6to30_in3Sbins.C"
#include "GetFixParamDeviationsy193to000to193.C"
#include "GetFixParamDeviations_HFNtracks.C"
#include "GetFixParamDeviationsXSBins.C"
#include "GetFixParamDeviations_Combine.C"
#include "GetFixParamDeviationspt0to6to30_Combine.C"
#include "GetFixParamDeviationspt0to6to30_in3Sbins_Combine.C"
#include "GetFixParamDeviationsy193to000to193_Combine.C"
#include "GetFixParamDeviations_HFNtracks_Combine.C"
#include "GetFixParamDeviationsXSBins_Combine.C"
#include "NewJaebeomStyle_Combined.C"
#include "NewJaebeomStyle_OneParam.C"
#include "NewJaebeomStyle_ChangeNtuplestoHistos.C"
#include "NewJaebeomStyle_JustFixedParams.C"
//#include "NewMakeAllLatexTables.C"

void GetFixParamDeviationsAll() {

  //Get deviations in yields due to changing each parameter
  GetFixParamDeviations(1,"alpha");
  GetFixParamDeviations(2,"alpha");
  GetFixParamDeviations(3,"alpha");
  GetFixParamDeviations(1,"n");
  GetFixParamDeviations(2,"n");
  GetFixParamDeviations(3,"n");
  GetFixParamDeviations(1,"x");
  GetFixParamDeviations(2,"x");
  GetFixParamDeviations(3,"x");
  GetFixParamDeviations(1,"f");
  GetFixParamDeviations(2,"f");
  GetFixParamDeviations(3,"f");

  GetFixParamDeviationspt0to6to30(1,"alpha");
  GetFixParamDeviationspt0to6to30(2,"alpha");
  GetFixParamDeviationspt0to6to30(3,"alpha");
  GetFixParamDeviationspt0to6to30(1,"n");
  GetFixParamDeviationspt0to6to30(2,"n");
  GetFixParamDeviationspt0to6to30(3,"n");
  GetFixParamDeviationspt0to6to30(1,"x");
  GetFixParamDeviationspt0to6to30(2,"x");
  GetFixParamDeviationspt0to6to30(3,"x");
  GetFixParamDeviationspt0to6to30(1,"f");
  GetFixParamDeviationspt0to6to30(2,"f");
  GetFixParamDeviationspt0to6to30(3,"f");

  GetFixParamDeviationsy193to000to193(1,"alpha");
  GetFixParamDeviationsy193to000to193(2,"alpha");
  GetFixParamDeviationsy193to000to193(3,"alpha");
  GetFixParamDeviationsy193to000to193(1,"n");
  GetFixParamDeviationsy193to000to193(2,"n");
  GetFixParamDeviationsy193to000to193(3,"n");
  GetFixParamDeviationsy193to000to193(1,"x");
  GetFixParamDeviationsy193to000to193(2,"x");
  GetFixParamDeviationsy193to000to193(3,"x");
  GetFixParamDeviationsy193to000to193(1,"f");
  GetFixParamDeviationsy193to000to193(2,"f");
  GetFixParamDeviationsy193to000to193(3,"f");

  GetFixParamDeviationsXSBins(1,"alpha");
  GetFixParamDeviationsXSBins(2,"alpha");
  GetFixParamDeviationsXSBins(3,"alpha");
  GetFixParamDeviationsXSBins(1,"n");
  GetFixParamDeviationsXSBins(2,"n");
  GetFixParamDeviationsXSBins(3,"n");
  GetFixParamDeviationsXSBins(1,"x");
  GetFixParamDeviationsXSBins(2,"x");
  GetFixParamDeviationsXSBins(3,"x");
  GetFixParamDeviationsXSBins(1,"f");
  GetFixParamDeviationsXSBins(2,"f");
  GetFixParamDeviationsXSBins(3,"f");

  GetFixParamDeviations_HFNtracks(1,0,"alpha");
  GetFixParamDeviations_HFNtracks(2,0,"alpha");
  GetFixParamDeviations_HFNtracks(3,0,"alpha");
  GetFixParamDeviations_HFNtracks(1,0,"n");
  GetFixParamDeviations_HFNtracks(2,0,"n");
  GetFixParamDeviations_HFNtracks(3,0,"n");
  GetFixParamDeviations_HFNtracks(1,0,"x");
  GetFixParamDeviations_HFNtracks(2,0,"x");
  GetFixParamDeviations_HFNtracks(3,0,"x");
  GetFixParamDeviations_HFNtracks(1,0,"f");
  GetFixParamDeviations_HFNtracks(2,0,"f");
  GetFixParamDeviations_HFNtracks(3,0,"f");

  GetFixParamDeviations_HFNtracks(1,1,"alpha");
  GetFixParamDeviations_HFNtracks(2,1,"alpha");
  GetFixParamDeviations_HFNtracks(3,1,"alpha");
  GetFixParamDeviations_HFNtracks(1,1,"n");
  GetFixParamDeviations_HFNtracks(2,1,"n");
  GetFixParamDeviations_HFNtracks(3,1,"n");
  GetFixParamDeviations_HFNtracks(1,1,"x");
  GetFixParamDeviations_HFNtracks(2,1,"x");
  GetFixParamDeviations_HFNtracks(3,1,"x");
  GetFixParamDeviations_HFNtracks(1,1,"f");
  GetFixParamDeviations_HFNtracks(2,1,"f");
  GetFixParamDeviations_HFNtracks(3,1,"f");

  GetFixParamDeviationspt0to6to30_in3Sbins(1,"alpha");
  GetFixParamDeviationspt0to6to30_in3Sbins(2,"alpha");
  GetFixParamDeviationspt0to6to30_in3Sbins(3,"alpha");
  GetFixParamDeviationspt0to6to30_in3Sbins(1,"n");
  GetFixParamDeviationspt0to6to30_in3Sbins(2,"n");
  GetFixParamDeviationspt0to6to30_in3Sbins(3,"n");
  GetFixParamDeviationspt0to6to30_in3Sbins(1,"x");
  GetFixParamDeviationspt0to6to30_in3Sbins(2,"x");
  GetFixParamDeviationspt0to6to30_in3Sbins(3,"x");
  GetFixParamDeviationspt0to6to30_in3Sbins(1,"f");
  GetFixParamDeviationspt0to6to30_in3Sbins(2,"f");
  GetFixParamDeviationspt0to6to30_in3Sbins(3,"f");

  //Combine the results of all the parameters to get a systematic uncertainty
  GetFixParamDeviations_Combine(1);
  GetFixParamDeviations_Combine(2);
  GetFixParamDeviations_Combine(3);

  GetFixParamDeviationspt0to6to30_Combine(1);
  GetFixParamDeviationspt0to6to30_Combine(2);
  GetFixParamDeviationspt0to6to30_Combine(3);

  GetFixParamDeviationsy193to000to193_Combine(1);
  GetFixParamDeviationsy193to000to193_Combine(2);
  GetFixParamDeviationsy193to000to193_Combine(3);

  GetFixParamDeviationsXSBins_Combine(1);
  GetFixParamDeviationsXSBins_Combine(2);
  GetFixParamDeviationsXSBins_Combine(3);

  GetFixParamDeviations_HFNtracks_Combine(1,0);
  GetFixParamDeviations_HFNtracks_Combine(2,0);
  GetFixParamDeviations_HFNtracks_Combine(3,0);

  GetFixParamDeviations_HFNtracks_Combine(1,1);
  GetFixParamDeviations_HFNtracks_Combine(2,1);
  GetFixParamDeviations_HFNtracks_Combine(3,1);

  GetFixParamDeviationspt0to6to30_in3Sbins_Combine(1);
  GetFixParamDeviationspt0to6to30_in3Sbins_Combine(2);
  GetFixParamDeviationspt0to6to30_in3Sbins_Combine(3);

  //Combine the fixed param uncertainty with the signal pdf uncertainty in a root file
  NewJaebeomStyle_Combined(1);
  NewJaebeomStyle_Combined(2);
  NewJaebeomStyle_Combined(3);

  //Make root files containing just the fixed param uncertainty
  NewJaebeomStyle_JustFixedParams(1);
  NewJaebeomStyle_JustFixedParams(2);
  NewJaebeomStyle_JustFixedParams(3);

  //Make root files containing just the signal pdf uncertainty
  NewJaebeomStyle_ChangeNtuplestoHistos(1);
  NewJaebeomStyle_ChangeNtuplestoHistos(2);
  NewJaebeomStyle_ChangeNtuplestoHistos(3);

  //Make root files containing just the uncertainty due to each parameter individually
  NewJaebeomStyle_OneParam(1,"alpha");
  NewJaebeomStyle_OneParam(2,"alpha");
  NewJaebeomStyle_OneParam(3,"alpha");
  NewJaebeomStyle_OneParam(1,"n");
  NewJaebeomStyle_OneParam(2,"n");
  NewJaebeomStyle_OneParam(3,"n");
  NewJaebeomStyle_OneParam(1,"x");
  NewJaebeomStyle_OneParam(2,"x");
  NewJaebeomStyle_OneParam(3,"x");
  NewJaebeomStyle_OneParam(1,"f");
  NewJaebeomStyle_OneParam(2,"f");
  NewJaebeomStyle_OneParam(3,"f");

}
