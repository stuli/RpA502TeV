#include "GetChisquare.C"
#include "PlotChiSquare.C"

void PlotAllChisquare() {

  int collId=kPADATA;

  /*GetChisquare(collId,1,2);
  GetChisquare(collId,1,3);
  GetChisquare(collId,1,4);
  GetChisquare(collId,2,2);
  GetChisquare(collId,2,3);
  GetChisquare(collId,2,4);
  GetChisquare(collId,3,2);
  GetChisquare(collId,3,3);
  GetChisquare(collId,3,4);
*/
  for (int whichUpsilon=1; whichUpsilon<4; whichUpsilon++) {
    PlotChiSquare(collId,whichUpsilon);
  }

}
