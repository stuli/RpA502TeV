#include "NewGetFixParamDeviations.C"
#include "NewGetFixParamDeviationspt0to6to30.C"
/*#include "NewGetFixParamDeviationspt0to6to30_in3Sbins.C"
#include "NewGetFixParamDeviationsy193to000to193.C"
#include "NewGetFixParamDeviations_HFNtracks.C"
#include "NewGetFixParamDeviationsXSBins.C"
#include "GetFixParamDeviations_Combine.C"
#include "GetFixParamDeviationspt0to6to30_Combine.C"
#include "GetFixParamDeviationspt0to6to30_in3Sbins_Combine.C"
#include "GetFixParamDeviationsy193to000to193_Combine.C"
#include "GetFixParamDeviations_HFNtracks_Combine.C"
#include "GetFixParamDeviationsXSBins_Combine.C"
#include "NewJaebeomStyle_Combined.C"
#include "NewJaebeomStyle_OneParam.C"
#include "NewJaebeomStyle_ChangeNtuplestoHistos.C"
#include "NewJaebeomStyle_JustFixedParams.C"*/
//#include "NewMakeAllLatexTables.C"

void NewGetFixParamDeviationsAll() {

  //Get deviations in yields due to changing each parameter
  /*NewGetFixParamDeviations(1,"alpha");
  NewGetFixParamDeviations(2,"alpha");
  NewGetFixParamDeviations(3,"alpha");
  NewGetFixParamDeviations(1,"n");
  NewGetFixParamDeviations(2,"n");
  NewGetFixParamDeviations(3,"n");
  NewGetFixParamDeviations(1,"x");
  NewGetFixParamDeviations(2,"x");
  NewGetFixParamDeviations(3,"x");
  NewGetFixParamDeviations(1,"f");
  NewGetFixParamDeviations(2,"f");
  NewGetFixParamDeviations(3,"f");*/

  NewGetFixParamDeviationspt0to6to30(1,"alpha");
  NewGetFixParamDeviationspt0to6to30(2,"alpha");
  NewGetFixParamDeviationspt0to6to30(3,"alpha");
  NewGetFixParamDeviationspt0to6to30(1,"n");
  NewGetFixParamDeviationspt0to6to30(2,"n");
  NewGetFixParamDeviationspt0to6to30(3,"n");
  NewGetFixParamDeviationspt0to6to30(1,"x");
  NewGetFixParamDeviationspt0to6to30(2,"x");
  NewGetFixParamDeviationspt0to6to30(3,"x");
  NewGetFixParamDeviationspt0to6to30(1,"f");
  NewGetFixParamDeviationspt0to6to30(2,"f");
  NewGetFixParamDeviationspt0to6to30(3,"f");

 /* NewGetFixParamDeviationspt0to6to30_in3Sbins(1,"alpha");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(2,"alpha");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(3,"alpha");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(1,"n");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(2,"n");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(3,"n");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(1,"x");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(2,"x");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(3,"x");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(1,"f");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(2,"f");
  NewGetFixParamDeviationspt0to6to30_in3Sbins(3,"f");

  NewGetFixParamDeviationsy193to000to193(1,"alpha");
  NewGetFixParamDeviationsy193to000to193(2,"alpha");
  NewGetFixParamDeviationsy193to000to193(3,"alpha");
  NewGetFixParamDeviationsy193to000to193(1,"n");
  NewGetFixParamDeviationsy193to000to193(2,"n");
  NewGetFixParamDeviationsy193to000to193(3,"n");
  NewGetFixParamDeviationsy193to000to193(1,"x");
  NewGetFixParamDeviationsy193to000to193(2,"x");
  NewGetFixParamDeviationsy193to000to193(3,"x");
  NewGetFixParamDeviationsy193to000to193(1,"f");
  NewGetFixParamDeviationsy193to000to193(2,"f");
  NewGetFixParamDeviationsy193to000to193(3,"f");

  NewGetFixParamDeviationsXSBins(1,"alpha");
  NewGetFixParamDeviationsXSBins(2,"alpha");
  NewGetFixParamDeviationsXSBins(3,"alpha");
  NewGetFixParamDeviationsXSBins(1,"n");
  NewGetFixParamDeviationsXSBins(2,"n");
  NewGetFixParamDeviationsXSBins(3,"n");
  NewGetFixParamDeviationsXSBins(1,"x");
  NewGetFixParamDeviationsXSBins(2,"x");
  NewGetFixParamDeviationsXSBins(3,"x");
  NewGetFixParamDeviationsXSBins(1,"f");
  NewGetFixParamDeviationsXSBins(2,"f");
  NewGetFixParamDeviationsXSBins(3,"f");

  NewGetFixParamDeviations_HFNtracks(1,0,"alpha");
  NewGetFixParamDeviations_HFNtracks(2,0,"alpha");
  NewGetFixParamDeviations_HFNtracks(3,0,"alpha");
  NewGetFixParamDeviations_HFNtracks(1,0,"n");
  NewGetFixParamDeviations_HFNtracks(2,0,"n");
  NewGetFixParamDeviations_HFNtracks(3,0,"n");
  NewGetFixParamDeviations_HFNtracks(1,0,"x");
  NewGetFixParamDeviations_HFNtracks(2,0,"x");
  NewGetFixParamDeviations_HFNtracks(3,0,"x");
  NewGetFixParamDeviations_HFNtracks(1,0,"f");
  NewGetFixParamDeviations_HFNtracks(2,0,"f");
  NewGetFixParamDeviations_HFNtracks(3,0,"f");

  NewGetFixParamDeviations_HFNtracks(1,1,"alpha");
  NewGetFixParamDeviations_HFNtracks(2,1,"alpha");
  NewGetFixParamDeviations_HFNtracks(3,1,"alpha");
  NewGetFixParamDeviations_HFNtracks(1,1,"n");
  NewGetFixParamDeviations_HFNtracks(2,1,"n");
  NewGetFixParamDeviations_HFNtracks(3,1,"n");
  NewGetFixParamDeviations_HFNtracks(1,1,"x");
  NewGetFixParamDeviations_HFNtracks(2,1,"x");
  NewGetFixParamDeviations_HFNtracks(3,1,"x");
  NewGetFixParamDeviations_HFNtracks(1,1,"f");
  NewGetFixParamDeviations_HFNtracks(2,1,"f");
  NewGetFixParamDeviations_HFNtracks(3,1,"f");

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
*/
}
