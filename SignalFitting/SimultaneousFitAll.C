#include "SimultaneousFits2ptAndFull.C"
#include "SimultaneousFits3ptConstrained.C"
#include "SimultaneousFits2ptConstrained.C"
#include "SimultaneousFits2yConstrained.C"
#include "SimultaneousFits4hfConstrained.C"
#include "SimultaneousFits4ntracksConstrained.C"
#include "SimultaneousChangeToNormal2ptAndFull.C"
#include "SimultaneousChangeToNormal3pt.C"
#include "SimultaneousChangeToNormal2y.C"
#include "SimultaneousChangeToNormal2pt.C"
#include "SimultaneousChangeToNormal4hf.C"
#include "SimultaneousChangeToNormal4ntracks.C"


SimultaneousFitAll( 
       int collId = kPPDATA,  
       bool whichModel=0   // Nominal = 0. Alternative = 1.
			)  {

  float yLow = 0.0;
  float yHigh = 1.93;
  int cLow=0;
  int cHigh=200;
  float muPtCut=4.0;

//P T   B I N S
  //integrated + 3S pt bins
  /*SimultaneousFits2ptAndFull(collId, 0, 30, yLow, yHigh, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2ptAndFull(collId, 0, 30, yLow, yHigh, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2ptAndFull(collId, 0, 30, yLow, yHigh, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2ptAndFull(collId, 0, 30, yLow, yHigh, "_C", cLow, cHigh, muPtCut, whichModel);

  //2S pt bins
  SimultaneousFits3ptConstrained(collId, 0, 30, yLow, yHigh, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 30, yLow, yHigh, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 30, yLow, yHigh, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 30, yLow, yHigh, "_C", cLow, cHigh, muPtCut, whichModel);

  //1S pt bins
  SimultaneousFits3ptConstrained(collId, 0, 6, yLow, yHigh, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 6, yLow, yHigh, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 6, yLow, yHigh, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 0, 6, yLow, yHigh, "_C", cLow, cHigh, muPtCut, whichModel);

  SimultaneousFits3ptConstrained(collId, 6, 30, yLow, yHigh, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 6, 30, yLow, yHigh, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 6, 30, yLow, yHigh, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal3pt(collId, 6, 30, yLow, yHigh, "_C", cLow, cHigh, muPtCut, whichModel);

  //3S pt bins in forward and backward rapidity
  SimultaneousFits2yConstrained(collId, 0, 6, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 6, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 6, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 6, 30, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 6, 30, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 6, 30, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

  //2S pt bins in forward and backward rapidity
  SimultaneousFits2yConstrained(collId, 0, 4, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 4, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 4, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 4, 9, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 4, 9, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 4, 9, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 9, 30, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 9, 30, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 9, 30, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

  //1S pt bins in forward and backward rapidity
  SimultaneousFits2yConstrained(collId, 0, 2, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 2, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 2, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 2, 4, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 2, 4, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 2, 4, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 4, 6, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 4, 6, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 4, 6, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 6, 9, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 6, 9, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 6, 9, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 9, 12, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 9, 12, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 9, 12, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 12, 30, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 12, 30, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 12, 30, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

//Y   B I N S
  //3S bins
  SimultaneousFits2yConstrained(collId, 0, 30, -1.93, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  //2S bins
  SimultaneousFits2yConstrained(collId, 0, 30, -1.93, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 0, 30, 0.0, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.0, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.0, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  //1S bins
  /*SimultaneousFits2yConstrained(collId, 0, 30, -1.93, -0.8, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, -0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -1.93, -0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 0, 30, -0.8, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -0.8, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, -0.8, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 0, 30, 0.0, 0.8, cLow, cHigh, muPtCut, whichModel);*/
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.0, 0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.0, 0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2yConstrained(collId, 0, 30, 0.8, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.8, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2y(collId, 0, 30, 0.8, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

  //3S y bins in low and high pt (already covered)

  //2S y bins in low and high pt
  /*SimultaneousFits2ptConstrained(collId, 0, 30, -1.93, -0.8, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.93, -0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.93, -0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, -0.8, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.8, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.8, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);*/
  SimultaneousFits2ptConstrained(collId, 0, 30, 0.0, 0.8, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.0, 0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.0, 0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, 0.8, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.8, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.8, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

  //1S y bins in low and high pt
  /*SimultaneousFits2ptConstrained(collId, 0, 30, -1.93, -1.2, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.93, -1.2, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.93, -1.2, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, -1.2, -0.8, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.2, -0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -1.2, -0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, -0.8, -0.4, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.8, -0.4, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.8, -0.4, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, -0.4, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.4, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, -0.4, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);*/
  SimultaneousFits2ptConstrained(collId, 0, 30, 0.0, 0.4, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.0, 0.4, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.0, 0.4, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, 0.4, 0.8, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.4, 0.8, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.4, 0.8, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, 0.8, 1.2, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.8, 1.2, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 0.8, 1.2, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits2ptConstrained(collId, 0, 30, 1.2, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 1.2, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal2pt(collId, 0, 30, 1.2, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);

//H F   B I N S
/*  SimultaneousFits4hfConstrained(collId, 0, 30, -1.93, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, -1.93, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, -1.93, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, -1.93, 0.0, "_C", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, -1.93, 0.0, "_D", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits4hfConstrained(collId, 0, 30, 0.0, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, 0.0, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, 0.0, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, 0.0, 1.93, "_C", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4hf(collId, 0, 30, 0.0, 1.93, "_D", cLow, cHigh, muPtCut, whichModel);

//N T R A C K S   B I N S
  SimultaneousFits4ntracksConstrained(collId, 0, 30, -1.93, 0.0, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, -1.93, 0.0, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, -1.93, 0.0, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, -1.93, 0.0, "_C", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, -1.93, 0.0, "_D", cLow, cHigh, muPtCut, whichModel);
  SimultaneousFits4ntracksConstrained(collId, 0, 30, 0.0, 1.93, cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, 0.0, 1.93, "_A", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, 0.0, 1.93, "_B", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, 0.0, 1.93, "_C", cLow, cHigh, muPtCut, whichModel);
  SimultaneousChangeToNormal4ntracks(collId, 0, 30, 0.0, 1.93, "_D", cLow, cHigh, muPtCut, whichModel);*/

}
