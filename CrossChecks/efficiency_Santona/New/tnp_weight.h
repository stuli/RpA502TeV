#ifndef tnp_weight_h

#define tnp_weight_h



#include "TMath.h"



// IN THIS FILE YOU WILL FIND:

// ---------------------------

//

// - MuIdTrg:

//   * idx = 1..100: toy variations, stat. only

//   * idx = -1: syst variation, "new_MAX", +1 sigma

//   * idx = -2: syst variation, "new_MAX", -1 sigma

// - STA:

//   * central values are identical in pp and pbpb

//   * idx = 1..100, pp: toy variations, stat. only

//   * idx = 1..100, pbpb: toy variations from pp, stat increased by 4% everywhere in each toy

//   * idx = -1: syst variation, "new_MAX", +1 sigma

//   * idx = -2: syst variation, "new_MAX", -1 sigma



/////////////////////////////////////////////////////////////

//               M U I D + T R G    P b P b                //

/////////////////////////////////////////////////////////////

double tnp_weight_muidtrg_pbpb(double x, double eta, int idx)

{

   // denominator (from MC)

   double den=1;

   if (fabs(eta)<0.9) den = 0.9724*TMath::Erf((x-0.4114)/3.3775);

   else if (fabs(eta)<1.6) den = 0.9502*TMath::Erf((x-1.3857)/2.0757);

   else if (fabs(eta)<2.1) den = 0.8971*TMath::Erf((x-1.0984)/2.3510);

   else den = 0.7763*TMath::Erf((x-0.8419)/1.6742);



   // numerator (from data)

   double num=1;



   if (fabs(eta)<0.9)

   {

      if (idx==0) num = 0.9646*TMath::Erf((x-0.1260)/3.5155);

      else if (idx == -1   ) num =   0.9689*TMath::Erf((x-0.2792)/3.3611);

      else if (idx == -2   ) num =   0.9603*TMath::Erf((x-0.0062)/3.6371);

      else if (idx == 1   ) num = 0.9580*TMath::Erf((x-0.0001)/3.3987);

      else if (idx == 2   ) num = 0.9700*TMath::Erf((x-0.0001)/3.8239);

      else if (idx == 3   ) num = 0.9584*TMath::Erf((x-0.0010)/3.5980);

      else if (idx == 4   ) num = 0.9586*TMath::Erf((x-0.0001)/3.6194);

      else if (idx == 5   ) num = 0.9677*TMath::Erf((x-0.4618)/3.3114);

      else if (idx == 6   ) num = 0.9532*TMath::Erf((x-0.0032)/3.5537);

      else if (idx == 7   ) num = 0.9552*TMath::Erf((x-0.4666)/3.2115);

      else if (idx == 8   ) num = 0.9699*TMath::Erf((x-0.4640)/3.1621);

      else if (idx == 9   ) num = 0.9632*TMath::Erf((x-0.0075)/3.5914);

      else if (idx == 10  ) num = 0.9661*TMath::Erf((x-0.0000)/3.6240);

      else if (idx == 11  ) num = 0.9655*TMath::Erf((x-0.0555)/3.5126);

      else if (idx == 12  ) num = 0.9512*TMath::Erf((x-0.0018)/3.5403);

      else if (idx == 13  ) num = 0.9637*TMath::Erf((x-0.4283)/3.3204);

      else if (idx == 14  ) num = 0.9714*TMath::Erf((x-0.0028)/3.5921);

      else if (idx == 15  ) num = 0.9606*TMath::Erf((x-0.0079)/3.5818);

      else if (idx == 16  ) num = 0.9762*TMath::Erf((x-0.0001)/3.7796);

      else if (idx == 17  ) num = 0.9726*TMath::Erf((x-0.2160)/3.5547);

      else if (idx == 18  ) num = 0.9620*TMath::Erf((x-0.0001)/3.3431);

      else if (idx == 19  ) num = 0.9628*TMath::Erf((x-0.0025)/3.6435);

      else if (idx == 20  ) num = 0.9636*TMath::Erf((x-0.0343)/3.5076);

      else if (idx == 21  ) num = 0.9579*TMath::Erf((x-0.0021)/3.5104);

      else if (idx == 22  ) num = 0.9659*TMath::Erf((x-0.3535)/3.4405);

      else if (idx == 23  ) num = 0.9636*TMath::Erf((x-0.0000)/3.5128);

      else if (idx == 24  ) num = 0.9679*TMath::Erf((x-0.0001)/3.5950);

      else if (idx == 25  ) num = 0.9781*TMath::Erf((x-0.1624)/3.5441);

      else if (idx == 26  ) num = 0.9661*TMath::Erf((x-0.1073)/3.4927);

      else if (idx == 27  ) num = 0.9668*TMath::Erf((x-0.0008)/3.7302);

      else if (idx == 28  ) num = 0.9689*TMath::Erf((x-0.0148)/3.6462);

      else if (idx == 29  ) num = 0.9546*TMath::Erf((x-0.0001)/3.5153);

      else if (idx == 30  ) num = 0.9697*TMath::Erf((x-0.4878)/3.3219);

      else if (idx == 31  ) num = 0.9665*TMath::Erf((x-0.0047)/3.5859);

      else if (idx == 32  ) num = 0.9577*TMath::Erf((x-0.0023)/3.5728);

      else if (idx == 33  ) num = 0.9637*TMath::Erf((x-1.1384)/2.6215);

      else if (idx == 34  ) num = 0.9681*TMath::Erf((x-0.0000)/3.6458);

      else if (idx == 35  ) num = 0.9586*TMath::Erf((x-0.6170)/2.9608);

      else if (idx == 36  ) num = 0.9592*TMath::Erf((x-0.0356)/3.5520);

      else if (idx == 37  ) num = 0.9728*TMath::Erf((x-0.4801)/3.3046);

      else if (idx == 38  ) num = 0.9687*TMath::Erf((x-0.0048)/3.5840);

      else if (idx == 39  ) num = 0.9665*TMath::Erf((x-0.0031)/3.7074);

      else if (idx == 40  ) num = 0.9692*TMath::Erf((x-0.4004)/3.3315);

      else if (idx == 41  ) num = 0.9621*TMath::Erf((x-0.0010)/3.6909);

      else if (idx == 42  ) num = 0.9640*TMath::Erf((x-0.1493)/3.5289);

      else if (idx == 43  ) num = 0.9550*TMath::Erf((x-1.1878)/2.5480);

      else if (idx == 44  ) num = 0.9621*TMath::Erf((x-0.1472)/3.4247);

      else if (idx == 45  ) num = 0.9575*TMath::Erf((x-1.1371)/2.5032);

      else if (idx == 46  ) num = 0.9625*TMath::Erf((x-0.0021)/3.5871);

      else if (idx == 47  ) num = 0.9600*TMath::Erf((x-0.0000)/3.6078);

      else if (idx == 48  ) num = 0.9666*TMath::Erf((x-0.0697)/3.4470);

      else if (idx == 49  ) num = 0.9559*TMath::Erf((x-0.0080)/3.6261);

      else if (idx == 50  ) num = 0.9640*TMath::Erf((x-0.0007)/3.6154);

      else if (idx == 51  ) num = 0.9648*TMath::Erf((x-0.0000)/3.6528);

      else if (idx == 52  ) num = 0.9631*TMath::Erf((x-0.1540)/3.3775);

      else if (idx == 53  ) num = 0.9524*TMath::Erf((x-1.1822)/2.5716);

      else if (idx == 54  ) num = 0.9669*TMath::Erf((x-0.7193)/2.9984);

      else if (idx == 55  ) num = 0.9611*TMath::Erf((x-0.0020)/3.5976);

      else if (idx == 56  ) num = 0.9628*TMath::Erf((x-0.0049)/3.5352);

      else if (idx == 57  ) num = 0.9722*TMath::Erf((x-0.0000)/3.7260);

      else if (idx == 58  ) num = 0.9605*TMath::Erf((x-0.1179)/3.5161);

      else if (idx == 59  ) num = 0.9661*TMath::Erf((x-0.0001)/3.6233);

      else if (idx == 60  ) num = 0.9631*TMath::Erf((x-0.8113)/2.8734);

      else if (idx == 61  ) num = 0.9688*TMath::Erf((x-0.3675)/3.3700);

      else if (idx == 62  ) num = 0.9590*TMath::Erf((x-1.0540)/2.6185);

      else if (idx == 63  ) num = 0.9699*TMath::Erf((x-0.5861)/3.2081);

      else if (idx == 64  ) num = 0.9648*TMath::Erf((x-0.0587)/3.5212);

      else if (idx == 65  ) num = 0.9589*TMath::Erf((x-1.4212)/2.3152);

      else if (idx == 66  ) num = 0.9567*TMath::Erf((x-1.2436)/2.5294);

      else if (idx == 67  ) num = 0.9669*TMath::Erf((x-0.5060)/3.2238);

      else if (idx == 68  ) num = 0.9653*TMath::Erf((x-0.7611)/2.9642);

      else if (idx == 69  ) num = 0.9499*TMath::Erf((x-1.5697)/2.1802);

      else if (idx == 70  ) num = 0.9536*TMath::Erf((x-0.0007)/3.5021);

      else if (idx == 71  ) num = 0.9579*TMath::Erf((x-0.0001)/3.4084);

      else if (idx == 72  ) num = 0.9646*TMath::Erf((x-0.2299)/3.4096);

      else if (idx == 73  ) num = 0.9643*TMath::Erf((x-0.7476)/2.9606);

      else if (idx == 74  ) num = 0.9678*TMath::Erf((x-0.0091)/3.7134);

      else if (idx == 75  ) num = 0.9707*TMath::Erf((x-0.0002)/3.6366);

      else if (idx == 76  ) num = 0.9623*TMath::Erf((x-0.0000)/3.5325);

      else if (idx == 77  ) num = 0.9579*TMath::Erf((x-1.3155)/2.3733);

      else if (idx == 78  ) num = 0.9621*TMath::Erf((x-0.0132)/3.6240);

      else if (idx == 79  ) num = 0.9656*TMath::Erf((x-0.0631)/3.5508);

      else if (idx == 80  ) num = 0.9594*TMath::Erf((x-0.0977)/3.5039);

      else if (idx == 81  ) num = 0.9608*TMath::Erf((x-0.1607)/3.4106);

      else if (idx == 82  ) num = 0.9648*TMath::Erf((x-0.0007)/3.5032);

      else if (idx == 83  ) num = 0.9609*TMath::Erf((x-0.0015)/3.5921);

      else if (idx == 84  ) num = 0.9590*TMath::Erf((x-0.0270)/3.5173);

      else if (idx == 85  ) num = 0.9472*TMath::Erf((x-1.5578)/2.0177);

      else if (idx == 86  ) num = 0.9620*TMath::Erf((x-0.3578)/3.2582);

      else if (idx == 87  ) num = 0.9652*TMath::Erf((x-0.0517)/3.6790);

      else if (idx == 88  ) num = 0.9622*TMath::Erf((x-0.0043)/3.6722);

      else if (idx == 89  ) num = 0.9581*TMath::Erf((x-0.0000)/3.4619);

      else if (idx == 90  ) num = 0.9659*TMath::Erf((x-0.0148)/3.5563);

      else if (idx == 91  ) num = 0.9607*TMath::Erf((x-1.4799)/2.3074);

      else if (idx == 92  ) num = 0.9635*TMath::Erf((x-0.3932)/3.2461);

      else if (idx == 93  ) num = 0.9648*TMath::Erf((x-0.9900)/2.7390);

      else if (idx == 94  ) num = 0.9574*TMath::Erf((x-0.0075)/3.6679);

      else if (idx == 95  ) num = 0.9603*TMath::Erf((x-0.5632)/3.1398);

      else if (idx == 96  ) num = 0.9586*TMath::Erf((x-0.4193)/3.2136);

      else if (idx == 97  ) num = 0.9695*TMath::Erf((x-0.0000)/3.8195);

      else if (idx == 98  ) num = 0.9589*TMath::Erf((x-0.0020)/3.6145);

      else if (idx == 99  ) num = 0.9628*TMath::Erf((x-0.2956)/3.2855);

      else if (idx == 100 ) num = 0.9673*TMath::Erf((x-0.4306)/3.3085);

   }

   else if (fabs(eta)<1.6)

   {

      if (idx==0) num = 0.9725*TMath::Erf((x-1.0054)/2.3187);

      else if (idx == -1   ) num =   0.9747*TMath::Erf((x-0.8107)/2.4694);

      else if (idx == -2   ) num =   0.9698*TMath::Erf((x-1.1745)/2.1832);

      else if (idx==1   ) num = 0.9664*TMath::Erf((x-1.3309)/1.9451);

      else if (idx==2   ) num = 0.9632*TMath::Erf((x-1.3290)/1.9513);

      else if (idx==3   ) num = 0.9678*TMath::Erf((x-1.0148)/2.3376);

      else if (idx==4   ) num = 0.9848*TMath::Erf((x-0.9444)/2.3932);

      else if (idx==5   ) num = 0.9646*TMath::Erf((x-1.3502)/1.8603);

      else if (idx==6   ) num = 0.9765*TMath::Erf((x-0.9802)/2.2845);

      else if (idx==7   ) num = 0.9862*TMath::Erf((x-0.5519)/2.9480);

      else if (idx==8   ) num = 0.9588*TMath::Erf((x-1.3311)/1.9312);

      else if (idx==9   ) num = 0.9557*TMath::Erf((x-1.3483)/1.8008);

      else if (idx==10  ) num = 0.9711*TMath::Erf((x-1.3174)/1.9481);

      else if (idx==11  ) num = 0.9834*TMath::Erf((x-0.7039)/2.7280);

      else if (idx==12  ) num = 0.9727*TMath::Erf((x-1.0827)/2.2481);

      else if (idx==13  ) num = 0.9650*TMath::Erf((x-1.3450)/1.8623);

      else if (idx==14  ) num = 0.9626*TMath::Erf((x-1.0382)/2.2124);

      else if (idx==15  ) num = 0.9799*TMath::Erf((x-0.4489)/2.8839);

      else if (idx==16  ) num = 0.9967*TMath::Erf((x-0.0000)/3.5322);

      else if (idx==17  ) num = 0.9694*TMath::Erf((x-1.1527)/2.1037);

      else if (idx==18  ) num = 0.9679*TMath::Erf((x-0.8587)/2.4470);

      else if (idx==19  ) num = 0.9581*TMath::Erf((x-1.2331)/1.9368);

      else if (idx==20  ) num = 0.9735*TMath::Erf((x-0.6093)/2.7565);

      else if (idx==21  ) num = 0.9707*TMath::Erf((x-0.9561)/2.3765);

      else if (idx==22  ) num = 0.9710*TMath::Erf((x-1.4346)/1.9103);

      else if (idx==23  ) num = 0.9739*TMath::Erf((x-0.7883)/2.5295);

      else if (idx==24  ) num = 0.9584*TMath::Erf((x-1.0871)/2.2088);

      else if (idx==25  ) num = 0.9683*TMath::Erf((x-1.1529)/2.1280);

      else if (idx==26  ) num = 0.9696*TMath::Erf((x-1.3294)/2.0530);

      else if (idx==27  ) num = 0.9631*TMath::Erf((x-1.3926)/1.9192);

      else if (idx==28  ) num = 0.9743*TMath::Erf((x-0.8310)/2.4924);

      else if (idx==29  ) num = 0.9724*TMath::Erf((x-0.9429)/2.3947);

      else if (idx==30  ) num = 0.9780*TMath::Erf((x-0.6907)/2.6474);

      else if (idx==31  ) num = 0.9750*TMath::Erf((x-1.2125)/2.1596);

      else if (idx==32  ) num = 0.9875*TMath::Erf((x-0.4066)/3.0293);

      else if (idx==33  ) num = 0.9435*TMath::Erf((x-1.7034)/1.3893);

      else if (idx==34  ) num = 0.9647*TMath::Erf((x-1.2425)/1.9835);

      else if (idx==35  ) num = 0.9643*TMath::Erf((x-0.5486)/2.7589);

      else if (idx==36  ) num = 0.9696*TMath::Erf((x-1.2034)/2.0922);

      else if (idx==37  ) num = 0.9647*TMath::Erf((x-1.1622)/2.1077);

      else if (idx==38  ) num = 0.9638*TMath::Erf((x-1.1493)/2.1684);

      else if (idx==39  ) num = 0.9608*TMath::Erf((x-1.5260)/1.7004);

      else if (idx==40  ) num = 0.9784*TMath::Erf((x-0.9245)/2.4133);

      else if (idx==41  ) num = 0.9694*TMath::Erf((x-0.6601)/2.6798);

      else if (idx==42  ) num = 0.9694*TMath::Erf((x-1.0725)/2.2150);

      else if (idx==43  ) num = 0.9709*TMath::Erf((x-1.1948)/2.0999);

      else if (idx==44  ) num = 0.9994*TMath::Erf((x-0.7046)/2.8191);

      else if (idx==45  ) num = 0.9717*TMath::Erf((x-1.1572)/2.1618);

      else if (idx==46  ) num = 0.9668*TMath::Erf((x-1.2549)/2.0258);

      else if (idx==47  ) num = 0.9881*TMath::Erf((x-0.7124)/2.6168);

      else if (idx==48  ) num = 0.9762*TMath::Erf((x-0.7719)/2.5990);

      else if (idx==49  ) num = 0.9681*TMath::Erf((x-1.3814)/1.8730);

      else if (idx==50  ) num = 0.9741*TMath::Erf((x-0.7377)/2.6037);

      else if (idx==51  ) num = 0.9698*TMath::Erf((x-1.0154)/2.3301);

      else if (idx==52  ) num = 0.9739*TMath::Erf((x-0.6024)/2.6989);

      else if (idx==53  ) num = 0.9654*TMath::Erf((x-1.4714)/1.7744);

      else if (idx==54  ) num = 0.9788*TMath::Erf((x-0.5921)/2.8095);

      else if (idx==55  ) num = 0.9674*TMath::Erf((x-1.3385)/2.0008);

      else if (idx==56  ) num = 0.9753*TMath::Erf((x-0.7505)/2.6561);

      else if (idx==57  ) num = 0.9915*TMath::Erf((x-0.5006)/3.0452);

      else if (idx==58  ) num = 0.9789*TMath::Erf((x-0.3582)/3.0704);

      else if (idx==59  ) num = 0.9813*TMath::Erf((x-0.1746)/3.1458);

      else if (idx==60  ) num = 0.9937*TMath::Erf((x-0.5639)/2.8380);

      else if (idx==61  ) num = 0.9752*TMath::Erf((x-0.6244)/2.7183);

      else if (idx==62  ) num = 0.9933*TMath::Erf((x-0.4229)/3.0276);

      else if (idx==63  ) num = 0.9676*TMath::Erf((x-0.9515)/2.2616);

      else if (idx==64  ) num = 0.9711*TMath::Erf((x-1.1906)/2.1208);

      else if (idx==65  ) num = 0.9772*TMath::Erf((x-1.1955)/2.1415);

      else if (idx==66  ) num = 0.9557*TMath::Erf((x-1.5576)/1.5638);

      else if (idx==67  ) num = 0.9995*TMath::Erf((x-0.5510)/2.9726);

      else if (idx==68  ) num = 0.9754*TMath::Erf((x-1.0703)/2.3058);

      else if (idx==69  ) num = 0.9755*TMath::Erf((x-1.1692)/2.1517);

      else if (idx==70  ) num = 0.9659*TMath::Erf((x-1.4382)/1.8701);

      else if (idx==71  ) num = 0.9804*TMath::Erf((x-0.8016)/2.6120);

      else if (idx==72  ) num = 0.9775*TMath::Erf((x-0.6922)/2.7078);

      else if (idx==73  ) num = 0.9713*TMath::Erf((x-0.8208)/2.5738);

      else if (idx==74  ) num = 0.9683*TMath::Erf((x-1.0708)/2.2729);

      else if (idx==75  ) num = 0.9691*TMath::Erf((x-0.6992)/2.5727);

      else if (idx==76  ) num = 0.9764*TMath::Erf((x-0.1513)/3.1129);

      else if (idx==77  ) num = 0.9704*TMath::Erf((x-1.0191)/2.3289);

      else if (idx==78  ) num = 0.9716*TMath::Erf((x-1.0818)/2.1600);

      else if (idx==79  ) num = 0.9777*TMath::Erf((x-0.6776)/2.7371);

      else if (idx==80  ) num = 0.9881*TMath::Erf((x-0.3183)/3.1203);

      else if (idx==81  ) num = 0.9974*TMath::Erf((x-0.7883)/2.7769);

      else if (idx==82  ) num = 0.9673*TMath::Erf((x-1.5485)/1.6782);

      else if (idx==83  ) num = 0.9926*TMath::Erf((x-0.8995)/2.4962);

      else if (idx==84  ) num = 0.9683*TMath::Erf((x-0.8082)/2.5190);

      else if (idx==85  ) num = 0.9807*TMath::Erf((x-1.0941)/2.2703);

      else if (idx==86  ) num = 0.9655*TMath::Erf((x-0.9080)/2.4409);

      else if (idx==87  ) num = 0.9713*TMath::Erf((x-0.8185)/2.5144);

      else if (idx==88  ) num = 0.9629*TMath::Erf((x-1.0726)/2.2500);

      else if (idx==89  ) num = 0.9828*TMath::Erf((x-0.3977)/2.9414);

      else if (idx==90  ) num = 0.9747*TMath::Erf((x-1.0780)/2.2479);

      else if (idx==91  ) num = 0.9648*TMath::Erf((x-1.2186)/2.0635);

      else if (idx==92  ) num = 0.9729*TMath::Erf((x-1.0200)/2.3461);

      else if (idx==93  ) num = 0.9751*TMath::Erf((x-1.1198)/2.1788);

      else if (idx==94  ) num = 0.9618*TMath::Erf((x-1.2555)/1.9968);

      else if (idx==95  ) num = 0.9930*TMath::Erf((x-0.5569)/2.8619);

      else if (idx==96  ) num = 0.9479*TMath::Erf((x-1.1822)/1.9777);

      else if (idx==97  ) num = 0.9653*TMath::Erf((x-1.4351)/1.8350);

      else if (idx==98  ) num = 0.9877*TMath::Erf((x-1.1694)/2.3072);

      else if (idx==99  ) num = 0.9845*TMath::Erf((x-1.2104)/2.1558);

      else if (idx==100 ) num = 0.9621*TMath::Erf((x-1.2094)/1.9538);

   }

   else if (fabs(eta)<2.1)

   {

      if (idx==0) num = 0.9194*TMath::Erf((x-0.9733)/2.1374);

      else if (idx == -1   ) num =   0.9231*TMath::Erf((x-0.8453)/2.2609);

      else if (idx == -2   ) num =   0.9147*TMath::Erf((x-1.0894)/2.0185);

      else if (idx == 1   ) num = 0.9136*TMath::Erf((x-0.8890)/2.2189);

      else if (idx == 2   ) num = 0.9143*TMath::Erf((x-1.1614)/1.8663);

      else if (idx == 3   ) num = 0.9219*TMath::Erf((x-1.0305)/2.0959);

      else if (idx == 4   ) num = 0.9389*TMath::Erf((x-0.8146)/2.5125);

      else if (idx == 5   ) num = 0.9282*TMath::Erf((x-0.8912)/2.2689);

      else if (idx == 6   ) num = 0.9089*TMath::Erf((x-1.1956)/1.7182);

      else if (idx == 7   ) num = 0.9060*TMath::Erf((x-1.5305)/1.3611);

      else if (idx == 8   ) num = 0.9330*TMath::Erf((x-0.7772)/2.4646);

      else if (idx == 9   ) num = 0.9373*TMath::Erf((x-0.8982)/2.3104);

      else if (idx == 10  ) num = 0.9334*TMath::Erf((x-0.8014)/2.3380);

      else if (idx == 11  ) num = 0.9171*TMath::Erf((x-1.2998)/1.7800);

      else if (idx == 12  ) num = 0.9205*TMath::Erf((x-1.2009)/1.8539);

      else if (idx == 13  ) num = 0.9227*TMath::Erf((x-1.0014)/2.1622);

      else if (idx == 14  ) num = 0.9106*TMath::Erf((x-1.1427)/1.8989);

      else if (idx == 15  ) num = 0.9198*TMath::Erf((x-1.0855)/1.9159);

      else if (idx == 16  ) num = 0.9030*TMath::Erf((x-1.1239)/1.8559);

      else if (idx == 17  ) num = 0.9125*TMath::Erf((x-1.0541)/1.9829);

      else if (idx == 18  ) num = 0.9332*TMath::Erf((x-0.8008)/2.4238);

      else if (idx == 19  ) num = 0.9224*TMath::Erf((x-0.9452)/2.1581);

      else if (idx == 20  ) num = 0.9066*TMath::Erf((x-1.1043)/1.9367);

      else if (idx == 21  ) num = 0.9079*TMath::Erf((x-1.2229)/1.7043);

      else if (idx == 22  ) num = 0.9118*TMath::Erf((x-1.0026)/2.0035);

      else if (idx == 23  ) num = 0.9319*TMath::Erf((x-0.8512)/2.4270);

      else if (idx == 24  ) num = 0.9116*TMath::Erf((x-1.2237)/1.7587);

      else if (idx == 25  ) num = 0.9123*TMath::Erf((x-0.9056)/2.2180);

      else if (idx == 26  ) num = 0.9135*TMath::Erf((x-0.8956)/2.2129);

      else if (idx == 27  ) num = 0.9100*TMath::Erf((x-1.2214)/1.7898);

      else if (idx == 28  ) num = 0.9136*TMath::Erf((x-0.6790)/2.5150);

      else if (idx == 29  ) num = 0.9196*TMath::Erf((x-1.1365)/1.8568);

      else if (idx == 30  ) num = 0.9027*TMath::Erf((x-0.9514)/2.1689);

      else if (idx == 31  ) num = 0.9181*TMath::Erf((x-1.3820)/1.5986);

      else if (idx == 32  ) num = 0.9074*TMath::Erf((x-1.0125)/2.1852);

      else if (idx == 33  ) num = 0.9063*TMath::Erf((x-1.1293)/2.0074);

      else if (idx == 34  ) num = 0.9219*TMath::Erf((x-0.8299)/2.3550);

      else if (idx == 35  ) num = 0.9139*TMath::Erf((x-1.0402)/1.9759);

      else if (idx == 36  ) num = 0.9251*TMath::Erf((x-0.9864)/2.1405);

      else if (idx == 37  ) num = 0.9180*TMath::Erf((x-1.3252)/1.6162);

      else if (idx == 38  ) num = 0.9115*TMath::Erf((x-0.8794)/2.2442);

      else if (idx == 39  ) num = 0.9300*TMath::Erf((x-0.9588)/2.2023);

      else if (idx == 40  ) num = 0.9358*TMath::Erf((x-0.5255)/2.8061);

      else if (idx == 41  ) num = 0.9179*TMath::Erf((x-0.9750)/2.0591);

      else if (idx == 42  ) num = 0.9224*TMath::Erf((x-1.1167)/2.0221);

      else if (idx == 43  ) num = 0.9262*TMath::Erf((x-1.3012)/1.7483);

      else if (idx == 44  ) num = 0.9186*TMath::Erf((x-0.9977)/2.1139);

      else if (idx == 45  ) num = 0.9117*TMath::Erf((x-0.8147)/2.2814);

      else if (idx == 46  ) num = 0.9203*TMath::Erf((x-0.8888)/2.2076);

      else if (idx == 47  ) num = 0.9361*TMath::Erf((x-0.8648)/2.2962);

      else if (idx == 48  ) num = 0.9273*TMath::Erf((x-1.1748)/1.8930);

      else if (idx == 49  ) num = 0.9213*TMath::Erf((x-0.7159)/2.2602);

      else if (idx == 50  ) num = 0.9110*TMath::Erf((x-0.8012)/2.2979);

      else if (idx == 51  ) num = 0.9277*TMath::Erf((x-1.0984)/2.0348);

      else if (idx == 52  ) num = 0.9281*TMath::Erf((x-1.0283)/2.1330);

      else if (idx == 53  ) num = 0.9197*TMath::Erf((x-0.9655)/2.2792);

      else if (idx == 54  ) num = 0.9193*TMath::Erf((x-0.9371)/2.1170);

      else if (idx == 55  ) num = 0.9282*TMath::Erf((x-0.6647)/2.5108);

      else if (idx == 56  ) num = 0.9269*TMath::Erf((x-0.9371)/2.2086);

      else if (idx == 57  ) num = 0.9390*TMath::Erf((x-1.1537)/2.0261);

      else if (idx == 58  ) num = 0.8994*TMath::Erf((x-1.0848)/1.9149);

      else if (idx == 59  ) num = 0.8955*TMath::Erf((x-1.0163)/1.9935);

      else if (idx == 60  ) num = 0.9070*TMath::Erf((x-0.9321)/2.1018);

      else if (idx == 61  ) num = 0.9081*TMath::Erf((x-0.8588)/2.2016);

      else if (idx == 62  ) num = 0.9061*TMath::Erf((x-1.3879)/1.5605);

      else if (idx == 63  ) num = 0.9117*TMath::Erf((x-1.3374)/1.6472);

      else if (idx == 64  ) num = 0.9175*TMath::Erf((x-1.1181)/1.9847);

      else if (idx == 65  ) num = 0.9306*TMath::Erf((x-1.2060)/1.8942);

      else if (idx == 66  ) num = 0.9108*TMath::Erf((x-0.8332)/2.2913);

      else if (idx == 67  ) num = 0.9112*TMath::Erf((x-1.0143)/2.0736);

      else if (idx == 68  ) num = 0.9223*TMath::Erf((x-0.8414)/2.2866);

      else if (idx == 69  ) num = 0.9252*TMath::Erf((x-1.0402)/2.1944);

      else if (idx == 70  ) num = 0.9255*TMath::Erf((x-0.9230)/2.1973);

      else if (idx == 71  ) num = 0.9353*TMath::Erf((x-0.8712)/2.4505);

      else if (idx == 72  ) num = 0.9339*TMath::Erf((x-0.9489)/2.2496);

      else if (idx == 73  ) num = 0.9057*TMath::Erf((x-1.1914)/1.8161);

      else if (idx == 74  ) num = 0.9176*TMath::Erf((x-0.9172)/2.1535);

      else if (idx == 75  ) num = 0.9309*TMath::Erf((x-0.8589)/2.4201);

      else if (idx == 76  ) num = 0.9212*TMath::Erf((x-0.9165)/2.1658);

      else if (idx == 77  ) num = 0.9241*TMath::Erf((x-1.0771)/2.0346);

      else if (idx == 78  ) num = 0.9169*TMath::Erf((x-0.8655)/2.2281);

      else if (idx == 79  ) num = 0.9130*TMath::Erf((x-1.1189)/1.9109);

      else if (idx == 80  ) num = 0.9394*TMath::Erf((x-0.8063)/2.5163);

      else if (idx == 81  ) num = 0.9217*TMath::Erf((x-1.0160)/2.0265);

      else if (idx == 82  ) num = 0.9114*TMath::Erf((x-1.3426)/1.6690);

      else if (idx == 83  ) num = 0.9273*TMath::Erf((x-0.6706)/2.5357);

      else if (idx == 84  ) num = 0.8885*TMath::Erf((x-1.0660)/1.8581);

      else if (idx == 85  ) num = 0.9254*TMath::Erf((x-1.0400)/1.9650);

      else if (idx == 86  ) num = 0.9264*TMath::Erf((x-1.0145)/2.0810);

      else if (idx == 87  ) num = 0.9179*TMath::Erf((x-0.8640)/2.2218);

      else if (idx == 88  ) num = 0.9270*TMath::Erf((x-0.9804)/2.1203);

      else if (idx == 89  ) num = 0.9163*TMath::Erf((x-1.0244)/2.0294);

      else if (idx == 90  ) num = 0.9051*TMath::Erf((x-0.9415)/2.1941);

      else if (idx == 91  ) num = 0.9134*TMath::Erf((x-0.5733)/2.6140);

      else if (idx == 92  ) num = 0.9277*TMath::Erf((x-0.8292)/2.3615);

      else if (idx == 93  ) num = 0.9091*TMath::Erf((x-0.9803)/2.0822);

      else if (idx == 94  ) num = 0.9193*TMath::Erf((x-1.0026)/2.0792);

      else if (idx == 95  ) num = 0.9217*TMath::Erf((x-1.2577)/1.8122);

      else if (idx == 96  ) num = 0.9014*TMath::Erf((x-1.2271)/1.6965);

      else if (idx == 97  ) num = 0.9216*TMath::Erf((x-0.9389)/2.1746);

      else if (idx == 98  ) num = 0.9150*TMath::Erf((x-0.9944)/2.1404);

      else if (idx == 99  ) num = 0.9097*TMath::Erf((x-0.8439)/2.1981);

      else if (idx == 100 ) num = 0.9104*TMath::Erf((x-0.8171)/2.4001);

   }

   else

   {

      if (idx==0) num = 0.8079*TMath::Erf((x-0.9421)/0.8577);

      else if (idx == -1   ) num =   0.8211*TMath::Erf((x-0.9282)/0.8999);

      else if (idx == -2   ) num =   0.7941*TMath::Erf((x-0.9906)/0.7813);

      else if (idx == 1   ) num = 0.8056*TMath::Erf((x-1.6103)/0.2200);

      else if (idx == 2   ) num = 0.8363*TMath::Erf((x-0.0000)/1.9077);

      else if (idx == 3   ) num = 0.8078*TMath::Erf((x-1.1152)/0.8519);

      else if (idx == 4   ) num = 0.7736*TMath::Erf((x-0.6190)/0.4166);

      else if (idx == 5   ) num = 0.8128*TMath::Erf((x-1.4446)/0.3582);

      else if (idx == 6   ) num = 0.8259*TMath::Erf((x-1.1372)/0.9449);

      else if (idx == 7   ) num = 0.8460*TMath::Erf((x-0.6870)/0.3865);

      else if (idx == 8   ) num = 0.8077*TMath::Erf((x-0.5292)/1.0800);

      else if (idx == 9   ) num = 0.8057*TMath::Erf((x-0.9524)/0.9679);

      else if (idx == 10  ) num = 0.8470*TMath::Erf((x-0.0012)/1.7617);

      else if (idx == 11  ) num = 0.8019*TMath::Erf((x-1.1453)/0.7969);

      else if (idx == 12  ) num = 0.8197*TMath::Erf((x-0.9877)/0.8917);

      else if (idx == 13  ) num = 0.7789*TMath::Erf((x-0.0004)/1.6429);

      else if (idx == 14  ) num = 0.7852*TMath::Erf((x-0.6597)/0.3558);

      else if (idx == 15  ) num = 0.8214*TMath::Erf((x-0.0000)/1.8619);

      else if (idx == 16  ) num = 0.7848*TMath::Erf((x-1.5317)/0.3120);

      else if (idx == 17  ) num = 0.8037*TMath::Erf((x-1.5570)/0.2359);

      else if (idx == 18  ) num = 0.7886*TMath::Erf((x-1.1624)/0.5437);

      else if (idx == 19  ) num = 0.7824*TMath::Erf((x-0.5396)/1.1359);

      else if (idx == 20  ) num = 0.7860*TMath::Erf((x-1.4176)/0.3419);

      else if (idx == 21  ) num = 0.8193*TMath::Erf((x-1.1025)/0.7039);

      else if (idx == 22  ) num = 0.7929*TMath::Erf((x-1.6091)/0.2073);

      else if (idx == 23  ) num = 0.7997*TMath::Erf((x-0.0001)/1.4745);

      else if (idx == 24  ) num = 0.7877*TMath::Erf((x-0.0000)/1.6610);

      else if (idx == 25  ) num = 0.7838*TMath::Erf((x-0.6438)/0.3553);

      else if (idx == 26  ) num = 0.8296*TMath::Erf((x-0.6897)/0.3656);

      else if (idx == 27  ) num = 0.8230*TMath::Erf((x-1.0366)/0.9573);

      else if (idx == 28  ) num = 0.8098*TMath::Erf((x-0.8357)/0.9431);

      else if (idx == 29  ) num = 0.7889*TMath::Erf((x-1.2595)/0.5656);

      else if (idx == 30  ) num = 0.8530*TMath::Erf((x-0.3186)/1.9275);

      else if (idx == 31  ) num = 0.8247*TMath::Erf((x-1.5087)/0.2973);

      else if (idx == 32  ) num = 0.7911*TMath::Erf((x-1.5843)/0.2086);

      else if (idx == 33  ) num = 0.7873*TMath::Erf((x-0.5589)/1.2147);

      else if (idx == 34  ) num = 0.7971*TMath::Erf((x-1.4626)/0.2584);

      else if (idx == 35  ) num = 0.8005*TMath::Erf((x-0.0000)/1.6135);

      else if (idx == 36  ) num = 0.7791*TMath::Erf((x-0.4788)/0.2388);

      else if (idx == 37  ) num = 0.7930*TMath::Erf((x-0.0045)/1.5717);

      else if (idx == 38  ) num = 0.7977*TMath::Erf((x-0.9392)/0.8600);

      else if (idx == 39  ) num = 0.7937*TMath::Erf((x-1.0338)/0.7507);

      else if (idx == 40  ) num = 0.7907*TMath::Erf((x-1.1912)/0.8897);

      else if (idx == 41  ) num = 0.7962*TMath::Erf((x-0.0000)/1.5652);

      else if (idx == 42  ) num = 0.8132*TMath::Erf((x-0.7345)/0.5072);

      else if (idx == 43  ) num = 0.8098*TMath::Erf((x-1.4911)/0.2409);

      else if (idx == 44  ) num = 0.8182*TMath::Erf((x-0.6608)/0.3465);

      else if (idx == 45  ) num = 0.8074*TMath::Erf((x-0.0000)/1.5236);

      else if (idx == 46  ) num = 0.8176*TMath::Erf((x-1.5892)/0.2123);

      else if (idx == 47  ) num = 0.7789*TMath::Erf((x-0.7727)/0.9241);

      else if (idx == 48  ) num = 0.8146*TMath::Erf((x-0.0000)/1.9448);

      else if (idx == 49  ) num = 0.8159*TMath::Erf((x-0.8485)/0.9537);

      else if (idx == 50  ) num = 0.8444*TMath::Erf((x-0.0000)/2.3424);

      else if (idx == 51  ) num = 0.8261*TMath::Erf((x-1.2050)/0.6735);

      else if (idx == 52  ) num = 0.8154*TMath::Erf((x-0.1403)/1.3527);

      else if (idx == 53  ) num = 0.8402*TMath::Erf((x-0.0002)/2.2443);

      else if (idx == 54  ) num = 0.8308*TMath::Erf((x-0.0000)/1.7752);

      else if (idx == 55  ) num = 0.7957*TMath::Erf((x-1.5440)/0.2214);

      else if (idx == 56  ) num = 0.7966*TMath::Erf((x-0.0031)/1.2847);

      else if (idx == 57  ) num = 0.8248*TMath::Erf((x-0.0000)/1.7527);

      else if (idx == 58  ) num = 0.8137*TMath::Erf((x-1.1530)/0.5978);

      else if (idx == 59  ) num = 0.8026*TMath::Erf((x-1.1938)/0.6154);

      else if (idx == 60  ) num = 0.7613*TMath::Erf((x-1.4313)/0.4128);

      else if (idx == 61  ) num = 0.8146*TMath::Erf((x-1.6101)/0.1989);

      else if (idx == 62  ) num = 0.8238*TMath::Erf((x-1.2892)/0.7766);

      else if (idx == 63  ) num = 0.7815*TMath::Erf((x-1.1236)/0.3757);

      else if (idx == 64  ) num = 0.8347*TMath::Erf((x-0.0000)/2.3934);

      else if (idx == 65  ) num = 0.8223*TMath::Erf((x-0.9849)/0.8888);

      else if (idx == 66  ) num = 0.7940*TMath::Erf((x-0.9704)/0.8937);

      else if (idx == 67  ) num = 0.8240*TMath::Erf((x-0.0119)/1.3856);

      else if (idx == 68  ) num = 0.7962*TMath::Erf((x-0.0000)/1.4659);

      else if (idx == 69  ) num = 0.8461*TMath::Erf((x-0.3805)/1.4052);

      else if (idx == 70  ) num = 0.8013*TMath::Erf((x-1.1735)/0.5411);

      else if (idx == 71  ) num = 0.7967*TMath::Erf((x-0.0000)/1.0902);

      else if (idx == 72  ) num = 0.8220*TMath::Erf((x-1.5990)/0.2022);

      else if (idx == 73  ) num = 0.7906*TMath::Erf((x-1.0907)/0.6622);

      else if (idx == 74  ) num = 0.8162*TMath::Erf((x-1.0903)/0.7047);

      else if (idx == 75  ) num = 0.8398*TMath::Erf((x-0.0000)/1.8480);

      else if (idx == 76  ) num = 0.8208*TMath::Erf((x-0.7020)/1.1439);

      else if (idx == 77  ) num = 0.8119*TMath::Erf((x-1.1593)/0.7509);

      else if (idx == 78  ) num = 0.8124*TMath::Erf((x-1.2988)/0.6738);

      else if (idx == 79  ) num = 0.8094*TMath::Erf((x-0.9686)/0.8786);

      else if (idx == 80  ) num = 0.8533*TMath::Erf((x-1.0974)/1.0263);

      else if (idx == 81  ) num = 0.7965*TMath::Erf((x-1.4786)/0.5360);

      else if (idx == 82  ) num = 0.7923*TMath::Erf((x-0.6557)/0.3657);

      else if (idx == 83  ) num = 0.7796*TMath::Erf((x-0.4551)/0.1290);

      else if (idx == 84  ) num = 0.8040*TMath::Erf((x-0.6722)/0.3595);

      else if (idx == 85  ) num = 0.8045*TMath::Erf((x-1.4161)/0.5137);

      else if (idx == 86  ) num = 0.8101*TMath::Erf((x-1.5392)/0.2451);

      else if (idx == 87  ) num = 0.8055*TMath::Erf((x-0.7864)/0.6117);

      else if (idx == 88  ) num = 0.8407*TMath::Erf((x-0.0004)/1.8438);

      else if (idx == 89  ) num = 0.7927*TMath::Erf((x-1.4325)/0.3236);

      else if (idx == 90  ) num = 0.8128*TMath::Erf((x-0.9951)/1.0570);

      else if (idx == 91  ) num = 0.8446*TMath::Erf((x-1.6056)/0.2033);

      else if (idx == 92  ) num = 0.7835*TMath::Erf((x-0.6747)/0.9460);

      else if (idx == 93  ) num = 0.7984*TMath::Erf((x-0.8991)/1.0029);

      else if (idx == 94  ) num = 0.8150*TMath::Erf((x-0.0000)/1.4728);

      else if (idx == 95  ) num = 0.7871*TMath::Erf((x-0.9038)/0.9169);

      else if (idx == 96  ) num = 0.8038*TMath::Erf((x-1.3358)/0.3702);

      else if (idx == 97  ) num = 0.8265*TMath::Erf((x-0.3027)/1.7874);

      else if (idx == 98  ) num = 0.8565*TMath::Erf((x-0.0001)/2.4075);

      else if (idx == 99  ) num = 0.8269*TMath::Erf((x-0.6933)/1.0097);

      else if (idx == 100 ) num = 0.8245*TMath::Erf((x-1.1228)/0.8057);

   }



   // return

   return num/den;

}



/////////////////////////////////////////////////////////////

//               M U I D + T R G    P P                    //

/////////////////////////////////////////////////////////////

double tnp_weight_muidtrg_pp(double x, double eta, int idx)

{

   // denominator (from MC)

   double den=1;

   if (fabs(eta)<0.9) den = 0.9547*TMath::Erf((x-1.7776)/2.0497);

   else if (fabs(eta)<1.6) den = 0.9150*TMath::Erf((x-1.8502)/1.6651);

   else if (fabs(eta)<2.1) den = 0.8721*TMath::Erf((x-1.1449)/2.5504);

   else den = 0.6137*TMath::Erf((x-1.0202)/1.0729);



   // numerator (from data)

   double num=1;

   if (fabs(eta)<0.9)

   {

      if (idx==0) num = 0.9452*TMath::Erf((x-1.9895)/1.6646);

      else if (idx == -1   ) num =   0.9464*TMath::Erf((x-1.9423)/1.6791);

      else if (idx == -2   ) num =   0.9444*TMath::Erf((x-2.0398)/1.6463);

      else if (idx == 1   ) num = 0.9483*TMath::Erf((x-1.3362)/2.1673);

      else if (idx == 2   ) num = 0.9477*TMath::Erf((x-1.6542)/2.0202);

      else if (idx == 3   ) num = 0.9423*TMath::Erf((x-1.8122)/1.8067);

      else if (idx == 4   ) num = 0.9428*TMath::Erf((x-1.7459)/1.8705);

      else if (idx == 5   ) num = 0.9470*TMath::Erf((x-2.0749)/1.6280);

      else if (idx == 6   ) num = 0.9375*TMath::Erf((x-2.0581)/1.5661);

      else if (idx == 7   ) num = 0.9374*TMath::Erf((x-2.3451)/1.3376);

      else if (idx == 8   ) num = 0.9503*TMath::Erf((x-2.1375)/1.5256);

      else if (idx == 9   ) num = 0.9459*TMath::Erf((x-1.7659)/1.8494);

      else if (idx == 10  ) num = 0.9457*TMath::Erf((x-1.9827)/1.6586);

      else if (idx == 11  ) num = 0.9441*TMath::Erf((x-2.2392)/1.4087);

      else if (idx == 12  ) num = 0.9377*TMath::Erf((x-1.9663)/1.6504);

      else if (idx == 13  ) num = 0.9450*TMath::Erf((x-2.0095)/1.6808);

      else if (idx == 14  ) num = 0.9495*TMath::Erf((x-1.9847)/1.6555);

      else if (idx == 15  ) num = 0.9431*TMath::Erf((x-1.9126)/1.7116);

      else if (idx == 16  ) num = 0.9546*TMath::Erf((x-1.5625)/2.1135);

      else if (idx == 17  ) num = 0.9527*TMath::Erf((x-1.6592)/2.0238);

      else if (idx == 18  ) num = 0.9518*TMath::Erf((x-1.3956)/2.1133);

      else if (idx == 19  ) num = 0.9487*TMath::Erf((x-1.4774)/2.1432);

      else if (idx == 20  ) num = 0.9458*TMath::Erf((x-1.9738)/1.6501);

      else if (idx == 21  ) num = 0.9457*TMath::Erf((x-1.5396)/2.0235);

      else if (idx == 22  ) num = 0.9427*TMath::Erf((x-2.2034)/1.4990);

      else if (idx == 23  ) num = 0.9471*TMath::Erf((x-1.8383)/1.7690);

      else if (idx == 24  ) num = 0.9503*TMath::Erf((x-1.6669)/1.9466);

      else if (idx == 25  ) num = 0.9531*TMath::Erf((x-1.9542)/1.7182);

      else if (idx == 26  ) num = 0.9464*TMath::Erf((x-2.0405)/1.6082);

      else if (idx == 27  ) num = 0.9493*TMath::Erf((x-1.4188)/2.2079);

      else if (idx == 28  ) num = 0.9442*TMath::Erf((x-2.2467)/1.4232);

      else if (idx == 29  ) num = 0.9413*TMath::Erf((x-1.7924)/1.7962);

      else if (idx == 30  ) num = 0.9491*TMath::Erf((x-1.9612)/1.7506);

      else if (idx == 31  ) num = 0.9476*TMath::Erf((x-1.8944)/1.7391);

      else if (idx == 32  ) num = 0.9442*TMath::Erf((x-1.5842)/2.0101);

      else if (idx == 33  ) num = 0.9492*TMath::Erf((x-2.1584)/1.5558);

      else if (idx == 34  ) num = 0.9515*TMath::Erf((x-1.4319)/2.1722);

      else if (idx == 35  ) num = 0.9431*TMath::Erf((x-2.3194)/1.3397);

      else if (idx == 36  ) num = 0.9404*TMath::Erf((x-2.1599)/1.4844);

      else if (idx == 37  ) num = 0.9527*TMath::Erf((x-1.8584)/1.8439);

      else if (idx == 38  ) num = 0.9498*TMath::Erf((x-1.8268)/1.8022);

      else if (idx == 39  ) num = 0.9454*TMath::Erf((x-1.9100)/1.7528);

      else if (idx == 40  ) num = 0.9509*TMath::Erf((x-1.8071)/1.8729);

      else if (idx == 41  ) num = 0.9459*TMath::Erf((x-1.5921)/2.0436);

      else if (idx == 42  ) num = 0.9459*TMath::Erf((x-1.8862)/1.7767);

      else if (idx == 43  ) num = 0.9421*TMath::Erf((x-2.3702)/1.3372);

      else if (idx == 44  ) num = 0.9438*TMath::Erf((x-2.1081)/1.5347);

      else if (idx == 45  ) num = 0.9460*TMath::Erf((x-2.1597)/1.5050);

      else if (idx == 46  ) num = 0.9446*TMath::Erf((x-1.9137)/1.7179);

      else if (idx == 47  ) num = 0.9476*TMath::Erf((x-1.3248)/2.2502);

      else if (idx == 48  ) num = 0.9477*TMath::Erf((x-2.0664)/1.5615);

      else if (idx == 49  ) num = 0.9401*TMath::Erf((x-1.8205)/1.8024);

      else if (idx == 50  ) num = 0.9464*TMath::Erf((x-1.7524)/1.8710);

      else if (idx == 51  ) num = 0.9447*TMath::Erf((x-1.9361)/1.7098);

      else if (idx == 52  ) num = 0.9448*TMath::Erf((x-2.1506)/1.4861);

      else if (idx == 53  ) num = 0.9412*TMath::Erf((x-2.2622)/1.4482);

      else if (idx == 54  ) num = 0.9504*TMath::Erf((x-2.0114)/1.6768);

      else if (idx == 55  ) num = 0.9444*TMath::Erf((x-1.7792)/1.8362);

      else if (idx == 56  ) num = 0.9448*TMath::Erf((x-1.9716)/1.6499);

      else if (idx == 57  ) num = 0.9491*TMath::Erf((x-1.8497)/1.8145);

      else if (idx == 58  ) num = 0.9433*TMath::Erf((x-1.9355)/1.7135);

      else if (idx == 59  ) num = 0.9519*TMath::Erf((x-1.2819)/2.3067);

      else if (idx == 60  ) num = 0.9478*TMath::Erf((x-2.1621)/1.5231);

      else if (idx == 61  ) num = 0.9465*TMath::Erf((x-2.1703)/1.5236);

      else if (idx == 62  ) num = 0.9466*TMath::Erf((x-2.1639)/1.5152);

      else if (idx == 63  ) num = 0.9491*TMath::Erf((x-2.0715)/1.6405);

      else if (idx == 64  ) num = 0.9433*TMath::Erf((x-2.2480)/1.4026);

      else if (idx == 65  ) num = 0.9478*TMath::Erf((x-2.2398)/1.4695);

      else if (idx == 66  ) num = 0.9447*TMath::Erf((x-2.2079)/1.5108);

      else if (idx == 67  ) num = 0.9502*TMath::Erf((x-1.8212)/1.8593);

      else if (idx == 68  ) num = 0.9500*TMath::Erf((x-1.9555)/1.7335);

      else if (idx == 69  ) num = 0.9418*TMath::Erf((x-2.4049)/1.3157);

      else if (idx == 70  ) num = 0.9416*TMath::Erf((x-1.6954)/1.8805);

      else if (idx == 71  ) num = 0.9466*TMath::Erf((x-1.5939)/1.9504);

      else if (idx == 72  ) num = 0.9440*TMath::Erf((x-2.1999)/1.4639);

      else if (idx == 73  ) num = 0.9482*TMath::Erf((x-2.0565)/1.6300);

      else if (idx == 74  ) num = 0.9484*TMath::Erf((x-1.6706)/1.9870);

      else if (idx == 75  ) num = 0.9489*TMath::Erf((x-1.9007)/1.7411);

      else if (idx == 76  ) num = 0.9449*TMath::Erf((x-1.9233)/1.6891);

      else if (idx == 77  ) num = 0.9460*TMath::Erf((x-2.3657)/1.3302);

      else if (idx == 78  ) num = 0.9428*TMath::Erf((x-2.0271)/1.6259);

      else if (idx == 79  ) num = 0.9503*TMath::Erf((x-1.6080)/2.0177);

      else if (idx == 80  ) num = 0.9418*TMath::Erf((x-2.0468)/1.5975);

      else if (idx == 81  ) num = 0.9438*TMath::Erf((x-2.0536)/1.5859);

      else if (idx == 82  ) num = 0.9469*TMath::Erf((x-1.9686)/1.6451);

      else if (idx == 83  ) num = 0.9468*TMath::Erf((x-1.5151)/2.0741);

      else if (idx == 84  ) num = 0.9431*TMath::Erf((x-1.9347)/1.6883);

      else if (idx == 85  ) num = 0.9415*TMath::Erf((x-2.4854)/1.1732);

      else if (idx == 86  ) num = 0.9437*TMath::Erf((x-2.2034)/1.4558);

      else if (idx == 87  ) num = 0.9398*TMath::Erf((x-2.3238)/1.3617);

      else if (idx == 88  ) num = 0.9417*TMath::Erf((x-2.0500)/1.6110);

      else if (idx == 89  ) num = 0.9476*TMath::Erf((x-1.4043)/2.1322);

      else if (idx == 90  ) num = 0.9457*TMath::Erf((x-2.0425)/1.5920);

      else if (idx == 91  ) num = 0.9467*TMath::Erf((x-2.4370)/1.2966);

      else if (idx == 92  ) num = 0.9460*TMath::Erf((x-2.0574)/1.6029);

      else if (idx == 93  ) num = 0.9485*TMath::Erf((x-2.2023)/1.4971);

      else if (idx == 94  ) num = 0.9412*TMath::Erf((x-1.7708)/1.8630);

      else if (idx == 95  ) num = 0.9435*TMath::Erf((x-2.1435)/1.5407);

      else if (idx == 96  ) num = 0.9413*TMath::Erf((x-2.2512)/1.4161);

      else if (idx == 97  ) num = 0.9493*TMath::Erf((x-1.5259)/2.1535);

      else if (idx == 98  ) num = 0.9382*TMath::Erf((x-2.2689)/1.3880);

      else if (idx == 99  ) num = 0.9466*TMath::Erf((x-1.9826)/1.6576);

      else if (idx == 100 ) num = 0.9471*TMath::Erf((x-2.0518)/1.6388);

   }

   else if (fabs(eta)<1.6)

   {

      if (idx==0) num = 0.9296*TMath::Erf((x-1.8082)/1.4939);

      else if (idx == -1   ) num =   0.9321*TMath::Erf((x-1.7773)/1.5122);

      else if (idx == -2   ) num =   0.9273*TMath::Erf((x-1.8360)/1.4806);

      else if (idx == 1   ) num = 0.9264*TMath::Erf((x-1.8875)/1.3760);

      else if (idx == 2   ) num = 0.9235*TMath::Erf((x-1.8870)/1.3730);

      else if (idx == 3   ) num = 0.9278*TMath::Erf((x-1.7989)/1.5252);

      else if (idx == 4   ) num = 0.9378*TMath::Erf((x-1.7873)/1.5336);

      else if (idx == 5   ) num = 0.9253*TMath::Erf((x-1.8995)/1.3282);

      else if (idx == 6   ) num = 0.9337*TMath::Erf((x-1.8036)/1.4840);

      else if (idx == 7   ) num = 0.9407*TMath::Erf((x-1.6706)/1.7630);

      else if (idx == 8   ) num = 0.9215*TMath::Erf((x-1.8879)/1.3677);

      else if (idx == 9   ) num = 0.9176*TMath::Erf((x-1.9161)/1.2630);

      else if (idx == 10  ) num = 0.9300*TMath::Erf((x-1.8800)/1.3852);

      else if (idx == 11  ) num = 0.9365*TMath::Erf((x-1.7313)/1.6384);

      else if (idx == 12  ) num = 0.9313*TMath::Erf((x-1.8144)/1.4981);

      else if (idx == 13  ) num = 0.9253*TMath::Erf((x-1.8993)/1.3256);

      else if (idx == 14  ) num = 0.9229*TMath::Erf((x-1.8319)/1.4252);

      else if (idx == 15  ) num = 0.9368*TMath::Erf((x-1.6867)/1.6717);

      else if (idx == 16  ) num = 0.9333*TMath::Erf((x-1.6338)/1.7441);

      else if (idx == 17  ) num = 0.9283*TMath::Erf((x-1.8472)/1.4162);

      else if (idx == 18  ) num = 0.9276*TMath::Erf((x-1.7757)/1.5324);

      else if (idx == 19  ) num = 0.9200*TMath::Erf((x-1.8830)/1.3214);

      else if (idx == 20  ) num = 0.9312*TMath::Erf((x-1.7169)/1.6361);

      else if (idx == 21  ) num = 0.9275*TMath::Erf((x-1.8050)/1.4923);

      else if (idx == 22  ) num = 0.9242*TMath::Erf((x-1.9137)/1.3586);

      else if (idx == 23  ) num = 0.9266*TMath::Erf((x-1.7902)/1.4916);

      else if (idx == 24  ) num = 0.9229*TMath::Erf((x-1.8184)/1.4797);

      else if (idx == 25  ) num = 0.9286*TMath::Erf((x-1.8398)/1.4427);

      else if (idx == 26  ) num = 0.9280*TMath::Erf((x-1.8693)/1.4530);

      else if (idx == 27  ) num = 0.9223*TMath::Erf((x-1.9033)/1.3640);

      else if (idx == 28  ) num = 0.9313*TMath::Erf((x-1.7731)/1.5407);

      else if (idx == 29  ) num = 0.9276*TMath::Erf((x-1.8049)/1.4932);

      else if (idx == 30  ) num = 0.9341*TMath::Erf((x-1.7404)/1.5928);

      else if (idx == 31  ) num = 0.9339*TMath::Erf((x-1.8282)/1.5132);

      else if (idx == 32  ) num = 0.9384*TMath::Erf((x-1.6769)/1.7083);

      else if (idx == 33  ) num = 0.9127*TMath::Erf((x-1.9940)/1.1590);

      else if (idx == 34  ) num = 0.9255*TMath::Erf((x-1.8690)/1.3769);

      else if (idx == 35  ) num = 0.9210*TMath::Erf((x-1.7517)/1.5339);

      else if (idx == 36  ) num = 0.9285*TMath::Erf((x-1.8517)/1.4297);

      else if (idx == 37  ) num = 0.9260*TMath::Erf((x-1.8451)/1.4299);

      else if (idx == 38  ) num = 0.9240*TMath::Erf((x-1.8427)/1.4443);

      else if (idx == 39  ) num = 0.9244*TMath::Erf((x-1.9306)/1.3111);

      else if (idx == 40  ) num = 0.9290*TMath::Erf((x-1.8034)/1.4858);

      else if (idx == 41  ) num = 0.9271*TMath::Erf((x-1.7430)/1.5798);

      else if (idx == 42  ) num = 0.9297*TMath::Erf((x-1.8176)/1.4784);

      else if (idx == 43  ) num = 0.9298*TMath::Erf((x-1.8480)/1.4378);

      else if (idx == 44  ) num = 0.9389*TMath::Erf((x-1.7440)/1.6258);

      else if (idx == 45  ) num = 0.9294*TMath::Erf((x-1.8417)/1.4497);

      else if (idx == 46  ) num = 0.9249*TMath::Erf((x-1.8741)/1.3843);

      else if (idx == 47  ) num = 0.9308*TMath::Erf((x-1.7932)/1.4701);

      else if (idx == 48  ) num = 0.9326*TMath::Erf((x-1.7477)/1.5957);

      else if (idx == 49  ) num = 0.9282*TMath::Erf((x-1.8963)/1.3579);

      else if (idx == 50  ) num = 0.9305*TMath::Erf((x-1.7537)/1.5697);

      else if (idx == 51  ) num = 0.9296*TMath::Erf((x-1.7946)/1.5325);

      else if (idx == 52  ) num = 0.9276*TMath::Erf((x-1.7512)/1.5393);

      else if (idx == 53  ) num = 0.9267*TMath::Erf((x-1.9182)/1.3305);

      else if (idx == 54  ) num = 0.9363*TMath::Erf((x-1.6945)/1.6922);

      else if (idx == 55  ) num = 0.9283*TMath::Erf((x-1.8731)/1.4350);

      else if (idx == 56  ) num = 0.9322*TMath::Erf((x-1.7391)/1.6233);

      else if (idx == 57  ) num = 0.9301*TMath::Erf((x-1.7322)/1.6207);

      else if (idx == 58  ) num = 0.9331*TMath::Erf((x-1.6716)/1.7101);

      else if (idx == 59  ) num = 0.9349*TMath::Erf((x-1.6648)/1.6767);

      else if (idx == 60  ) num = 0.9346*TMath::Erf((x-1.7493)/1.5603);

      else if (idx == 61  ) num = 0.9306*TMath::Erf((x-1.7353)/1.5899);

      else if (idx == 62  ) num = 0.9371*TMath::Erf((x-1.7022)/1.6562);

      else if (idx == 63  ) num = 0.9261*TMath::Erf((x-1.8193)/1.4257);

      else if (idx == 64  ) num = 0.9305*TMath::Erf((x-1.8407)/1.4610);

      else if (idx == 65  ) num = 0.9296*TMath::Erf((x-1.8597)/1.4205);

      else if (idx == 66  ) num = 0.9202*TMath::Erf((x-1.9556)/1.2193);

      else if (idx == 67  ) num = 0.9459*TMath::Erf((x-1.6792)/1.7514);

      else if (idx == 68  ) num = 0.9302*TMath::Erf((x-1.8159)/1.5023);

      else if (idx == 69  ) num = 0.9281*TMath::Erf((x-1.8570)/1.4148);

      else if (idx == 70  ) num = 0.9272*TMath::Erf((x-1.9003)/1.3859);

      else if (idx == 71  ) num = 0.9323*TMath::Erf((x-1.7591)/1.5872);

      else if (idx == 72  ) num = 0.9330*TMath::Erf((x-1.7306)/1.6292);

      else if (idx == 73  ) num = 0.9311*TMath::Erf((x-1.7467)/1.6193);

      else if (idx == 74  ) num = 0.9289*TMath::Erf((x-1.8080)/1.5168);

      else if (idx == 75  ) num = 0.9263*TMath::Erf((x-1.7666)/1.5153);

      else if (idx == 76  ) num = 0.9324*TMath::Erf((x-1.6672)/1.6576);

      else if (idx == 77  ) num = 0.9292*TMath::Erf((x-1.8013)/1.5186);

      else if (idx == 78  ) num = 0.9300*TMath::Erf((x-1.8332)/1.4296);

      else if (idx == 79  ) num = 0.9324*TMath::Erf((x-1.7319)/1.6263);

      else if (idx == 80  ) num = 0.9422*TMath::Erf((x-1.6389)/1.7814);

      else if (idx == 81  ) num = 0.9373*TMath::Erf((x-1.7510)/1.6393);

      else if (idx == 82  ) num = 0.9211*TMath::Erf((x-1.9584)/1.2399);

      else if (idx == 83  ) num = 0.9369*TMath::Erf((x-1.7917)/1.5262);

      else if (idx == 84  ) num = 0.9268*TMath::Erf((x-1.7681)/1.5443);

      else if (idx == 85  ) num = 0.9285*TMath::Erf((x-1.8426)/1.4420);

      else if (idx == 86  ) num = 0.9246*TMath::Erf((x-1.7895)/1.5218);

      else if (idx == 87  ) num = 0.9288*TMath::Erf((x-1.7705)/1.5423);

      else if (idx == 88  ) num = 0.9227*TMath::Erf((x-1.8294)/1.4586);

      else if (idx == 89  ) num = 0.9372*TMath::Erf((x-1.6828)/1.6709);

      else if (idx == 90  ) num = 0.9326*TMath::Erf((x-1.8144)/1.4983);

      else if (idx == 91  ) num = 0.9227*TMath::Erf((x-1.8622)/1.3988);

      else if (idx == 92  ) num = 0.9313*TMath::Erf((x-1.7989)/1.5346);

      else if (idx == 93  ) num = 0.9328*TMath::Erf((x-1.8280)/1.4666);

      else if (idx == 94  ) num = 0.9235*TMath::Erf((x-1.8678)/1.3891);

      else if (idx == 95  ) num = 0.9352*TMath::Erf((x-1.7395)/1.5874);

      else if (idx == 96  ) num = 0.9129*TMath::Erf((x-1.8773)/1.3181);

      else if (idx == 97  ) num = 0.9267*TMath::Erf((x-1.9054)/1.3577);

      else if (idx == 98  ) num = 0.9320*TMath::Erf((x-1.8461)/1.4877);

      else if (idx == 99  ) num = 0.9329*TMath::Erf((x-1.8560)/1.4395);

      else if (idx == 100 ) num = 0.9220*TMath::Erf((x-1.8830)/1.3148);

   }

   else if (fabs(eta)<2.1)

   {

      if (idx==0) num = 0.8808*TMath::Erf((x-0.8275)/2.6569);

      else if (idx == -1   ) num =   0.8880*TMath::Erf((x-0.8487)/2.6109);

      else if (idx == -2   ) num =   0.8741*TMath::Erf((x-0.7909)/2.7254);

      else if (idx == 1   ) num = 0.8770*TMath::Erf((x-0.7880)/2.6943);

      else if (idx == 2   ) num = 0.8796*TMath::Erf((x-0.8924)/2.5517);

      else if (idx == 3   ) num = 0.8837*TMath::Erf((x-0.8413)/2.6633);

      else if (idx == 4   ) num = 0.8946*TMath::Erf((x-0.7660)/2.8588);

      else if (idx == 5   ) num = 0.8839*TMath::Erf((x-0.8219)/2.6615);

      else if (idx == 6   ) num = 0.8824*TMath::Erf((x-0.8453)/2.5985);

      else if (idx == 7   ) num = 0.8801*TMath::Erf((x-1.0358)/2.3720);

      else if (idx == 8   ) num = 0.8926*TMath::Erf((x-0.7354)/2.8534);

      else if (idx == 9   ) num = 0.8916*TMath::Erf((x-0.8156)/2.7134);

      else if (idx == 10  ) num = 0.8916*TMath::Erf((x-0.7433)/2.7871);

      else if (idx == 11  ) num = 0.8811*TMath::Erf((x-0.9594)/2.5088);

      else if (idx == 12  ) num = 0.8855*TMath::Erf((x-0.9109)/2.5561);

      else if (idx == 13  ) num = 0.8807*TMath::Erf((x-0.8590)/2.6326);

      else if (idx == 14  ) num = 0.8778*TMath::Erf((x-0.8646)/2.6030);

      else if (idx == 15  ) num = 0.8845*TMath::Erf((x-0.8554)/2.5860);

      else if (idx == 16  ) num = 0.8645*TMath::Erf((x-0.9157)/2.4368);

      else if (idx == 17  ) num = 0.8769*TMath::Erf((x-0.8501)/2.5900);

      else if (idx == 18  ) num = 0.8918*TMath::Erf((x-0.7382)/2.8471);

      else if (idx == 19  ) num = 0.8820*TMath::Erf((x-0.8223)/2.6482);

      else if (idx == 20  ) num = 0.8726*TMath::Erf((x-0.8663)/2.5738);

      else if (idx == 21  ) num = 0.8773*TMath::Erf((x-0.8907)/2.5166);

      else if (idx == 22  ) num = 0.8736*TMath::Erf((x-0.8550)/2.5404);

      else if (idx == 23  ) num = 0.8911*TMath::Erf((x-0.7634)/2.8452);

      else if (idx == 24  ) num = 0.8766*TMath::Erf((x-0.9200)/2.4898);

      else if (idx == 25  ) num = 0.8732*TMath::Erf((x-0.8114)/2.6580);

      else if (idx == 26  ) num = 0.8797*TMath::Erf((x-0.7584)/2.7608);

      else if (idx == 27  ) num = 0.8771*TMath::Erf((x-0.9050)/2.5382);

      else if (idx == 28  ) num = 0.8784*TMath::Erf((x-0.6728)/2.8952);

      else if (idx == 29  ) num = 0.8802*TMath::Erf((x-0.9025)/2.4973);

      else if (idx == 30  ) num = 0.8697*TMath::Erf((x-0.8121)/2.6681);

      else if (idx == 31  ) num = 0.8782*TMath::Erf((x-1.0297)/2.3387);

      else if (idx == 32  ) num = 0.8742*TMath::Erf((x-0.8084)/2.7484);

      else if (idx == 33  ) num = 0.8726*TMath::Erf((x-0.8769)/2.6209);

      else if (idx == 34  ) num = 0.8819*TMath::Erf((x-0.7743)/2.7490);

      else if (idx == 35  ) num = 0.8788*TMath::Erf((x-0.8461)/2.5855);

      else if (idx == 36  ) num = 0.8871*TMath::Erf((x-0.8165)/2.6974);

      else if (idx == 37  ) num = 0.8891*TMath::Erf((x-0.9215)/2.5337);

      else if (idx == 38  ) num = 0.8754*TMath::Erf((x-0.7864)/2.7009);

      else if (idx == 39  ) num = 0.8866*TMath::Erf((x-0.8428)/2.6524);

      else if (idx == 40  ) num = 0.8923*TMath::Erf((x-0.6470)/2.9881);

      else if (idx == 41  ) num = 0.8818*TMath::Erf((x-0.8133)/2.6415);

      else if (idx == 42  ) num = 0.8884*TMath::Erf((x-0.8576)/2.6823);

      else if (idx == 43  ) num = 0.8869*TMath::Erf((x-0.9677)/2.4840);

      else if (idx == 44  ) num = 0.8844*TMath::Erf((x-0.8030)/2.7200);

      else if (idx == 45  ) num = 0.8740*TMath::Erf((x-0.7566)/2.7149);

      else if (idx == 46  ) num = 0.8800*TMath::Erf((x-0.8041)/2.6571);

      else if (idx == 47  ) num = 0.8940*TMath::Erf((x-0.7826)/2.7504);

      else if (idx == 48  ) num = 0.8867*TMath::Erf((x-0.9197)/2.5293);

      else if (idx == 49  ) num = 0.8798*TMath::Erf((x-0.7474)/2.6376);

      else if (idx == 50  ) num = 0.8747*TMath::Erf((x-0.7553)/2.7202);

      else if (idx == 51  ) num = 0.8894*TMath::Erf((x-0.8685)/2.6462);

      else if (idx == 52  ) num = 0.8881*TMath::Erf((x-0.8480)/2.6765);

      else if (idx == 53  ) num = 0.8910*TMath::Erf((x-0.7353)/2.9356);

      else if (idx == 54  ) num = 0.8825*TMath::Erf((x-0.7976)/2.6684);

      else if (idx == 55  ) num = 0.8855*TMath::Erf((x-0.7173)/2.8026);

      else if (idx == 56  ) num = 0.8837*TMath::Erf((x-0.8344)/2.6504);

      else if (idx == 57  ) num = 0.8965*TMath::Erf((x-0.9054)/2.6317);

      else if (idx == 58  ) num = 0.8707*TMath::Erf((x-0.8363)/2.6046);

      else if (idx == 59  ) num = 0.8629*TMath::Erf((x-0.8368)/2.5692);

      else if (idx == 60  ) num = 0.8717*TMath::Erf((x-0.7978)/2.6381);

      else if (idx == 61  ) num = 0.8706*TMath::Erf((x-0.7922)/2.6360);

      else if (idx == 62  ) num = 0.8760*TMath::Erf((x-0.9813)/2.4250);

      else if (idx == 63  ) num = 0.8789*TMath::Erf((x-0.9652)/2.4536);

      else if (idx == 64  ) num = 0.8811*TMath::Erf((x-0.8680)/2.6235);

      else if (idx == 65  ) num = 0.8876*TMath::Erf((x-0.9485)/2.5030);

      else if (idx == 66  ) num = 0.8757*TMath::Erf((x-0.7507)/2.7548);

      else if (idx == 67  ) num = 0.8752*TMath::Erf((x-0.8480)/2.6100);

      else if (idx == 68  ) num = 0.8821*TMath::Erf((x-0.7743)/2.7201);

      else if (idx == 69  ) num = 0.8891*TMath::Erf((x-0.8295)/2.7646);

      else if (idx == 70  ) num = 0.8876*TMath::Erf((x-0.7989)/2.7085);

      else if (idx == 71  ) num = 0.8927*TMath::Erf((x-0.7799)/2.8487);

      else if (idx == 72  ) num = 0.8916*TMath::Erf((x-0.8202)/2.7260);

      else if (idx == 73  ) num = 0.8752*TMath::Erf((x-0.8830)/2.5678);

      else if (idx == 74  ) num = 0.8778*TMath::Erf((x-0.8098)/2.6368);

      else if (idx == 75  ) num = 0.8916*TMath::Erf((x-0.7567)/2.8626);

      else if (idx == 76  ) num = 0.8779*TMath::Erf((x-0.8353)/2.5928);

      else if (idx == 77  ) num = 0.8840*TMath::Erf((x-0.8769)/2.5987);

      else if (idx == 78  ) num = 0.8794*TMath::Erf((x-0.7735)/2.7067);

      else if (idx == 79  ) num = 0.8762*TMath::Erf((x-0.8888)/2.5356);

      else if (idx == 80  ) num = 0.9004*TMath::Erf((x-0.7071)/2.9781);

      else if (idx == 81  ) num = 0.8831*TMath::Erf((x-0.8390)/2.6131);

      else if (idx == 82  ) num = 0.8810*TMath::Erf((x-0.9536)/2.5064);

      else if (idx == 83  ) num = 0.8849*TMath::Erf((x-0.7082)/2.8367);

      else if (idx == 84  ) num = 0.8638*TMath::Erf((x-0.8066)/2.6004);

      else if (idx == 85  ) num = 0.8862*TMath::Erf((x-0.8462)/2.5902);

      else if (idx == 86  ) num = 0.8824*TMath::Erf((x-0.8834)/2.5490);

      else if (idx == 87  ) num = 0.8816*TMath::Erf((x-0.7643)/2.7252);

      else if (idx == 88  ) num = 0.8826*TMath::Erf((x-0.8591)/2.5872);

      else if (idx == 89  ) num = 0.8795*TMath::Erf((x-0.8398)/2.6132);

      else if (idx == 90  ) num = 0.8769*TMath::Erf((x-0.7443)/2.8197);

      else if (idx == 91  ) num = 0.8752*TMath::Erf((x-0.6779)/2.8430);

      else if (idx == 92  ) num = 0.8881*TMath::Erf((x-0.7564)/2.7943);

      else if (idx == 93  ) num = 0.8733*TMath::Erf((x-0.8334)/2.6063);

      else if (idx == 94  ) num = 0.8827*TMath::Erf((x-0.8274)/2.6512);

      else if (idx == 95  ) num = 0.8845*TMath::Erf((x-0.9412)/2.5259);

      else if (idx == 96  ) num = 0.8725*TMath::Erf((x-0.8932)/2.5069);

      else if (idx == 97  ) num = 0.8879*TMath::Erf((x-0.7633)/2.7854);

      else if (idx == 98  ) num = 0.8820*TMath::Erf((x-0.7927)/2.7485);

      else if (idx == 99  ) num = 0.8809*TMath::Erf((x-0.7087)/2.8046);

      else if (idx == 100 ) num = 0.8775*TMath::Erf((x-0.7242)/2.8572);

   }

   else

   {

      if (idx==0) num = 0.7180*TMath::Erf((x-0.8578)/0.8700);

      else if (idx == -1   ) num =   0.7237*TMath::Erf((x-1.0005)/0.7219);

      else if (idx == -2   ) num =   0.7120*TMath::Erf((x-0.7996)/0.9479);

      else if (idx == 1   ) num = 0.7166*TMath::Erf((x-1.2516)/0.5359);

      else if (idx == 2   ) num = 0.7318*TMath::Erf((x-0.0027)/1.7106);

      else if (idx == 3   ) num = 0.7214*TMath::Erf((x-0.7102)/1.0657);

      else if (idx == 4   ) num = 0.7039*TMath::Erf((x-0.8217)/0.8051);

      else if (idx == 5   ) num = 0.7211*TMath::Erf((x-1.0845)/0.6735);

      else if (idx == 6   ) num = 0.7278*TMath::Erf((x-0.7072)/1.1182);

      else if (idx == 7   ) num = 0.7470*TMath::Erf((x-0.0013)/1.6604);

      else if (idx == 8   ) num = 0.7173*TMath::Erf((x-0.8430)/0.8568);

      else if (idx == 9   ) num = 0.7170*TMath::Erf((x-0.7897)/0.9600);

      else if (idx == 10  ) num = 0.7344*TMath::Erf((x-0.4984)/1.2313);

      else if (idx == 11  ) num = 0.7145*TMath::Erf((x-0.8770)/0.8858);

      else if (idx == 12  ) num = 0.7216*TMath::Erf((x-0.8724)/0.8814);

      else if (idx == 13  ) num = 0.7024*TMath::Erf((x-0.6774)/0.9737);

      else if (idx == 14  ) num = 0.7078*TMath::Erf((x-0.6645)/0.9321);

      else if (idx == 15  ) num = 0.7246*TMath::Erf((x-0.0143)/1.6680);

      else if (idx == 16  ) num = 0.7053*TMath::Erf((x-1.2393)/0.5302);

      else if (idx == 17  ) num = 0.7155*TMath::Erf((x-1.2468)/0.5147);

      else if (idx == 18  ) num = 0.7086*TMath::Erf((x-1.0463)/0.6646);

      else if (idx == 19  ) num = 0.7055*TMath::Erf((x-0.8309)/0.8505);

      else if (idx == 20  ) num = 0.7057*TMath::Erf((x-1.1795)/0.5500);

      else if (idx == 21  ) num = 0.7241*TMath::Erf((x-0.9261)/0.8229);

      else if (idx == 22  ) num = 0.7086*TMath::Erf((x-1.4130)/0.3678);

      else if (idx == 23  ) num = 0.7126*TMath::Erf((x-0.5529)/1.0944);

      else if (idx == 24  ) num = 0.7086*TMath::Erf((x-0.3100)/1.3014);

      else if (idx == 25  ) num = 0.7014*TMath::Erf((x-1.3142)/0.3928);

      else if (idx == 26  ) num = 0.7292*TMath::Erf((x-0.4158)/1.1805);

      else if (idx == 27  ) num = 0.7267*TMath::Erf((x-0.7178)/1.0722);

      else if (idx == 28  ) num = 0.7166*TMath::Erf((x-0.8528)/0.8656);

      else if (idx == 29  ) num = 0.7091*TMath::Erf((x-1.0373)/0.7004);

      else if (idx == 30  ) num = 0.7331*TMath::Erf((x-0.3617)/1.4460);

      else if (idx == 31  ) num = 0.7264*TMath::Erf((x-1.1119)/0.6594);

      else if (idx == 32  ) num = 0.7082*TMath::Erf((x-1.5147)/0.2648);

      else if (idx == 33  ) num = 0.7074*TMath::Erf((x-0.7820)/0.9116);

      else if (idx == 34  ) num = 0.7107*TMath::Erf((x-1.3653)/0.3824);

      else if (idx == 35  ) num = 0.7141*TMath::Erf((x-0.7533)/0.9374);

      else if (idx == 36  ) num = 0.7079*TMath::Erf((x-0.3588)/1.2020);

      else if (idx == 37  ) num = 0.7100*TMath::Erf((x-0.7398)/0.9350);

      else if (idx == 38  ) num = 0.7145*TMath::Erf((x-0.8533)/0.8685);

      else if (idx == 39  ) num = 0.7099*TMath::Erf((x-0.9747)/0.7461);

      else if (idx == 40  ) num = 0.7075*TMath::Erf((x-0.8300)/0.9575);

      else if (idx == 41  ) num = 0.7119*TMath::Erf((x-0.4227)/1.2104);

      else if (idx == 42  ) num = 0.7296*TMath::Erf((x-0.0441)/1.5493);

      else if (idx == 43  ) num = 0.7158*TMath::Erf((x-1.4735)/0.2933);

      else if (idx == 44  ) num = 0.7205*TMath::Erf((x-0.7557)/0.8847);

      else if (idx == 45  ) num = 0.7197*TMath::Erf((x-0.3019)/1.3366);

      else if (idx == 46  ) num = 0.7181*TMath::Erf((x-1.4132)/0.3676);

      else if (idx == 47  ) num = 0.7027*TMath::Erf((x-0.9052)/0.7791);

      else if (idx == 48  ) num = 0.7160*TMath::Erf((x-0.5359)/1.1692);

      else if (idx == 49  ) num = 0.7210*TMath::Erf((x-0.8620)/0.8733);

      else if (idx == 50  ) num = 0.7266*TMath::Erf((x-0.0249)/1.7357);

      else if (idx == 51  ) num = 0.7245*TMath::Erf((x-1.0032)/0.7698);

      else if (idx == 52  ) num = 0.7220*TMath::Erf((x-0.7383)/0.9563);

      else if (idx == 53  ) num = 0.7239*TMath::Erf((x-0.4208)/1.3443);

      else if (idx == 54  ) num = 0.7300*TMath::Erf((x-0.1245)/1.5572);

      else if (idx == 55  ) num = 0.7099*TMath::Erf((x-1.3518)/0.4030);

      else if (idx == 56  ) num = 0.7127*TMath::Erf((x-0.8248)/0.8412);

      else if (idx == 57  ) num = 0.7245*TMath::Erf((x-0.5167)/1.1879);

      else if (idx == 58  ) num = 0.7209*TMath::Erf((x-0.9908)/0.7445);

      else if (idx == 59  ) num = 0.7166*TMath::Erf((x-0.9852)/0.7568);

      else if (idx == 60  ) num = 0.6962*TMath::Erf((x-1.1946)/0.5513);

      else if (idx == 61  ) num = 0.7198*TMath::Erf((x-1.4611)/0.3291);

      else if (idx == 62  ) num = 0.7255*TMath::Erf((x-0.8597)/0.9723);

      else if (idx == 63  ) num = 0.7052*TMath::Erf((x-1.1680)/0.5282);

      else if (idx == 64  ) num = 0.7187*TMath::Erf((x-0.2292)/1.5254);

      else if (idx == 65  ) num = 0.7240*TMath::Erf((x-0.8342)/0.9212);

      else if (idx == 66  ) num = 0.7119*TMath::Erf((x-0.8580)/0.8737);

      else if (idx == 67  ) num = 0.7247*TMath::Erf((x-0.7476)/0.9426);

      else if (idx == 68  ) num = 0.7120*TMath::Erf((x-0.6907)/0.9676);

      else if (idx == 69  ) num = 0.7366*TMath::Erf((x-0.6924)/1.0656);

      else if (idx == 70  ) num = 0.7145*TMath::Erf((x-0.9956)/0.7190);

      else if (idx == 71  ) num = 0.7154*TMath::Erf((x-0.5779)/1.0455);

      else if (idx == 72  ) num = 0.7234*TMath::Erf((x-1.3667)/0.4170);

      else if (idx == 73  ) num = 0.7078*TMath::Erf((x-1.0463)/0.6724);

      else if (idx == 74  ) num = 0.7196*TMath::Erf((x-1.0048)/0.7375);

      else if (idx == 75  ) num = 0.7302*TMath::Erf((x-0.4342)/1.3019);

      else if (idx == 76  ) num = 0.7233*TMath::Erf((x-0.7556)/0.9829);

      else if (idx == 77  ) num = 0.7194*TMath::Erf((x-0.8984)/0.8662);

      else if (idx == 78  ) num = 0.7196*TMath::Erf((x-0.9647)/0.8319);

      else if (idx == 79  ) num = 0.7177*TMath::Erf((x-0.8629)/0.8744);

      else if (idx == 80  ) num = 0.7390*TMath::Erf((x-0.6601)/1.1948);

      else if (idx == 81  ) num = 0.7124*TMath::Erf((x-1.0556)/0.7668);

      else if (idx == 82  ) num = 0.7085*TMath::Erf((x-1.1835)/0.5101);

      else if (idx == 83  ) num = 0.6993*TMath::Erf((x-0.7678)/0.7849);

      else if (idx == 84  ) num = 0.7196*TMath::Erf((x-0.3082)/1.2741);

      else if (idx == 85  ) num = 0.7185*TMath::Erf((x-0.8952)/0.8965);

      else if (idx == 86  ) num = 0.7186*TMath::Erf((x-1.1725)/0.5826);

      else if (idx == 87  ) num = 0.7163*TMath::Erf((x-0.9281)/0.7490);

      else if (idx == 88  ) num = 0.7296*TMath::Erf((x-0.6708)/1.0681);

      else if (idx == 89  ) num = 0.7110*TMath::Erf((x-1.0808)/0.6446);

      else if (idx == 90  ) num = 0.7182*TMath::Erf((x-0.7241)/1.0576);

      else if (idx == 91  ) num = 0.7340*TMath::Erf((x-1.2615)/0.5320);

      else if (idx == 92  ) num = 0.7063*TMath::Erf((x-0.9237)/0.7618);

      else if (idx == 93  ) num = 0.7128*TMath::Erf((x-0.7948)/0.9402);

      else if (idx == 94  ) num = 0.7241*TMath::Erf((x-0.3625)/1.2925);

      else if (idx == 95  ) num = 0.7076*TMath::Erf((x-0.8464)/0.8645);

      else if (idx == 96  ) num = 0.7150*TMath::Erf((x-1.0450)/0.6694);

      else if (idx == 97  ) num = 0.7211*TMath::Erf((x-0.4410)/1.3040);

      else if (idx == 98  ) num = 0.7290*TMath::Erf((x-0.3473)/1.4477);

      else if (idx == 99  ) num = 0.7269*TMath::Erf((x-0.8277)/0.9027);

      else if (idx == 100 ) num = 0.7260*TMath::Erf((x-0.8638)/0.9165);

   }



   // return

   return num/den;

}



/////////////////////////////////////////////////////////////

//                     STA    P b P b                      //

/////////////////////////////////////////////////////////////

double tnp_weight_sta_pbpb(double x, double eta, int idx)

{

   // denominator (from MC)

   double den=1;

   if (fabs(eta)<1.6) den = 0.9911*TMath::Erf((x-1.4336)/2.8548);

   else den = den = 0.9132*TMath::Erf((x-0.8045)/1.8366);

   

   // numerator (from data)

   double num=1;

   if (fabs(eta)<1.6)

   {

      if (idx==0) num = 0.9891*TMath::Erf((x-1.4814)/2.5014);

      else if (idx == -1  ) num = 0.9996*TMath::Erf((x-1.4885)/2.4666);

      else if (idx == -2  ) num = 0.9764*TMath::Erf((x-1.4746)/2.5279);

      else if (idx == 1   ) num = 0.9948*TMath::Erf((x-1.0027)/3.0825);

      else if (idx == 2   ) num = 0.9658*TMath::Erf((x-1.8116)/1.9534);

      else if (idx == 3   ) num = 0.9601*TMath::Erf((x-1.5043)/2.5557);

      else if (idx == 4   ) num = 0.9755*TMath::Erf((x-1.5490)/2.4079);

      else if (idx == 5   ) num = 0.9580*TMath::Erf((x-1.4780)/2.3290);

      else if (idx == 6   ) num = 0.9459*TMath::Erf((x-1.5368)/2.1954);

      else if (idx == 7   ) num = 0.9847*TMath::Erf((x-1.6312)/2.2940);

      else if (idx == 8   ) num = 0.9696*TMath::Erf((x-1.6106)/2.3524);

      else if (idx == 9   ) num = 0.9579*TMath::Erf((x-1.6494)/2.1508);

      else if (idx == 10  ) num = 0.9611*TMath::Erf((x-1.3995)/2.3846);

      else if (idx == 11  ) num = 0.9553*TMath::Erf((x-1.5769)/2.4887);

      else if (idx == 12  ) num = 0.9927*TMath::Erf((x-1.2934)/2.8614);

      else if (idx == 13  ) num = 0.9689*TMath::Erf((x-1.6733)/2.1308);

      else if (idx == 14  ) num = 0.9677*TMath::Erf((x-1.5914)/2.3259);

      else if (idx == 15  ) num = 0.9498*TMath::Erf((x-1.6213)/2.2265);

      else if (idx == 16  ) num = 0.9981*TMath::Erf((x-1.5474)/2.6294);

      else if (idx == 17  ) num = 0.9566*TMath::Erf((x-1.6260)/2.2666);

      else if (idx == 18  ) num = 0.9637*TMath::Erf((x-1.7857)/1.9085);

      else if (idx == 19  ) num = 0.9891*TMath::Erf((x-1.1183)/3.0663);

      else if (idx == 20  ) num = 0.9625*TMath::Erf((x-1.7861)/1.9473);

      else if (idx == 21  ) num = 0.9816*TMath::Erf((x-1.4461)/2.6316);

      else if (idx == 22  ) num = 0.9533*TMath::Erf((x-1.6919)/1.9712);

      else if (idx == 23  ) num = 0.9781*TMath::Erf((x-1.6831)/2.1241);

      else if (idx == 24  ) num = 0.9391*TMath::Erf((x-1.6315)/2.0061);

      else if (idx == 25  ) num = 0.9685*TMath::Erf((x-1.6208)/2.2084);

      else if (idx == 26  ) num = 0.9378*TMath::Erf((x-1.6176)/1.8588);

      else if (idx == 27  ) num = 0.9851*TMath::Erf((x-1.5616)/2.5608);

      else if (idx == 28  ) num = 0.9715*TMath::Erf((x-1.4177)/2.4401);

      else if (idx == 29  ) num = 0.9565*TMath::Erf((x-1.4087)/2.6610);

      else if (idx == 30  ) num = 0.9866*TMath::Erf((x-1.1504)/3.0086);

      else if (idx == 31  ) num = 0.9743*TMath::Erf((x-1.1424)/2.8948);

      else if (idx == 32  ) num = 0.9435*TMath::Erf((x-1.5286)/2.2686);

      else if (idx == 33  ) num = 0.9657*TMath::Erf((x-1.5931)/2.2059);

      else if (idx == 34  ) num = 0.9857*TMath::Erf((x-1.2805)/2.9793);

      else if (idx == 35  ) num = 0.9663*TMath::Erf((x-1.3391)/2.6491);

      else if (idx == 36  ) num = 0.9759*TMath::Erf((x-1.4775)/2.5575);

      else if (idx == 37  ) num = 0.9835*TMath::Erf((x-1.4982)/2.6455);

      else if (idx == 38  ) num = 0.9946*TMath::Erf((x-1.3276)/2.7413);

      else if (idx == 39  ) num = 0.9515*TMath::Erf((x-1.6369)/2.1667);

      else if (idx == 40  ) num = 0.9508*TMath::Erf((x-1.4388)/2.3989);

      else if (idx == 41  ) num = 0.9888*TMath::Erf((x-1.4619)/2.4613);

      else if (idx == 42  ) num = 0.9731*TMath::Erf((x-1.2972)/2.7154);

      else if (idx == 43  ) num = 0.9985*TMath::Erf((x-1.4898)/2.6072);

      else if (idx == 44  ) num = 0.9320*TMath::Erf((x-1.4094)/2.3230);

      else if (idx == 45  ) num = 0.9511*TMath::Erf((x-1.3477)/2.4102);

      else if (idx == 46  ) num = 0.9864*TMath::Erf((x-1.4347)/2.6130);

      else if (idx == 47  ) num = 0.9957*TMath::Erf((x-1.2196)/2.8369);

      else if (idx == 48  ) num = 0.9479*TMath::Erf((x-1.6938)/2.1571);

      else if (idx == 49  ) num = 0.9630*TMath::Erf((x-1.4457)/2.5090);

      else if (idx == 50  ) num = 0.9838*TMath::Erf((x-1.5105)/2.4755);

      else if (idx == 51  ) num = 0.9864*TMath::Erf((x-1.4726)/2.5555);

      else if (idx == 52  ) num = 0.9537*TMath::Erf((x-1.7316)/1.9708);

      else if (idx == 53  ) num = 0.9763*TMath::Erf((x-1.2341)/2.6906);

      else if (idx == 54  ) num = 0.9688*TMath::Erf((x-1.6741)/2.1406);

      else if (idx == 55  ) num = 0.9708*TMath::Erf((x-1.6162)/2.2292);

      else if (idx == 56  ) num = 0.9565*TMath::Erf((x-1.2409)/2.7221);

      else if (idx == 57  ) num = 0.9670*TMath::Erf((x-1.3479)/2.6662);

      else if (idx == 58  ) num = 0.9569*TMath::Erf((x-1.8988)/1.8177);

      else if (idx == 59  ) num = 0.9848*TMath::Erf((x-1.0690)/2.9668);

      else if (idx == 60  ) num = 0.9516*TMath::Erf((x-1.5590)/2.4386);

      else if (idx == 61  ) num = 0.9798*TMath::Erf((x-1.0860)/3.0184);

      else if (idx == 62  ) num = 0.9580*TMath::Erf((x-1.5088)/2.5848);

      else if (idx == 63  ) num = 0.9396*TMath::Erf((x-1.6843)/2.1154);

      else if (idx == 64  ) num = 0.9742*TMath::Erf((x-0.9365)/3.1109);

      else if (idx == 65  ) num = 0.9481*TMath::Erf((x-1.6703)/2.0737);

      else if (idx == 66  ) num = 0.9642*TMath::Erf((x-1.6919)/2.1213);

      else if (idx == 67  ) num = 0.9828*TMath::Erf((x-1.6155)/2.4307);

      else if (idx == 68  ) num = 0.9851*TMath::Erf((x-1.5165)/2.5912);

      else if (idx == 69  ) num = 0.9627*TMath::Erf((x-1.3952)/2.6783);

      else if (idx == 70  ) num = 0.9556*TMath::Erf((x-1.5954)/2.0434);

      else if (idx == 71  ) num = 0.9558*TMath::Erf((x-1.7714)/1.9476);

      else if (idx == 72  ) num = 0.9417*TMath::Erf((x-1.3330)/2.4873);

      else if (idx == 73  ) num = 0.9641*TMath::Erf((x-1.4129)/2.4777);

      else if (idx == 74  ) num = 0.9753*TMath::Erf((x-1.1892)/2.8455);

      else if (idx == 75  ) num = 0.9616*TMath::Erf((x-1.4862)/2.3338);

      else if (idx == 76  ) num = 0.9850*TMath::Erf((x-1.2442)/2.9777);

      else if (idx == 77  ) num = 0.9609*TMath::Erf((x-1.6084)/2.2445);

      else if (idx == 78  ) num = 0.9847*TMath::Erf((x-1.4876)/2.6138);

      else if (idx == 79  ) num = 0.9382*TMath::Erf((x-1.4622)/2.3917);

      else if (idx == 80  ) num = 0.9697*TMath::Erf((x-1.5235)/2.2411);

      else if (idx == 81  ) num = 0.9709*TMath::Erf((x-1.7719)/1.9731);

      else if (idx == 82  ) num = 0.9544*TMath::Erf((x-1.8381)/1.7623);

      else if (idx == 83  ) num = 0.9970*TMath::Erf((x-1.2574)/2.7917);

      else if (idx == 84  ) num = 0.9995*TMath::Erf((x-1.3149)/2.8497);

      else if (idx == 85  ) num = 0.9650*TMath::Erf((x-1.3187)/2.4415);

      else if (idx == 86  ) num = 0.9700*TMath::Erf((x-1.6363)/2.1275);

      else if (idx == 87  ) num = 0.9342*TMath::Erf((x-1.6477)/2.1985);

      else if (idx == 88  ) num = 0.9425*TMath::Erf((x-1.6008)/2.1196);

      else if (idx == 89  ) num = 0.9784*TMath::Erf((x-1.2755)/2.7493);

      else if (idx == 90  ) num = 0.9691*TMath::Erf((x-1.5775)/2.3070);

      else if (idx == 91  ) num = 0.9697*TMath::Erf((x-1.7054)/2.1071);

      else if (idx == 92  ) num = 0.9575*TMath::Erf((x-1.3237)/2.6191);

      else if (idx == 93  ) num = 0.9691*TMath::Erf((x-1.4130)/2.5344);

      else if (idx == 94  ) num = 0.9677*TMath::Erf((x-1.4688)/2.4914);

      else if (idx == 95  ) num = 0.9685*TMath::Erf((x-1.6630)/1.9997);

      else if (idx == 96  ) num = 0.9705*TMath::Erf((x-1.5190)/2.3908);

      else if (idx == 97  ) num = 0.9537*TMath::Erf((x-1.5912)/2.0281);

      else if (idx == 98  ) num = 0.9409*TMath::Erf((x-1.3793)/2.3750);

      else if (idx == 99  ) num = 0.9441*TMath::Erf((x-1.5825)/2.1456);

      else if (idx == 100 ) num = 0.9794*TMath::Erf((x-1.3770)/2.5317);

   }

   else

   {

      if (idx==0) num = 0.8956*TMath::Erf((x-0.5162)/1.7646);

      else if (idx == -1  ) num = 0.9153*TMath::Erf((x-0.5049)/1.7401);

      else if (idx == -2  ) num = 0.8762*TMath::Erf((x-0.5250)/1.7942);

      else if (idx == 1   ) num = 0.8422*TMath::Erf((x-0.2660)/1.6859);

      else if (idx == 2   ) num = 0.9189*TMath::Erf((x-0.8973)/1.4428);

      else if (idx == 3   ) num = 0.8999*TMath::Erf((x-0.5813)/1.7805);

      else if (idx == 4   ) num = 0.9285*TMath::Erf((x-0.5231)/1.9030);

      else if (idx == 5   ) num = 0.9033*TMath::Erf((x-0.5506)/1.7935);

      else if (idx == 6   ) num = 0.9028*TMath::Erf((x-0.8598)/1.3702);

      else if (idx == 7   ) num = 0.9259*TMath::Erf((x-0.5527)/2.0627);

      else if (idx == 8   ) num = 0.8697*TMath::Erf((x-0.6804)/1.4575);

      else if (idx == 9   ) num = 0.9116*TMath::Erf((x-0.9879)/1.3085);

      else if (idx == 10  ) num = 0.9324*TMath::Erf((x-0.0001)/2.2847);

      else if (idx == 11  ) num = 0.8728*TMath::Erf((x-0.0000)/2.0122);

      else if (idx == 12  ) num = 0.8753*TMath::Erf((x-0.4843)/1.7428);

      else if (idx == 13  ) num = 0.8826*TMath::Erf((x-0.6228)/1.6058);

      else if (idx == 14  ) num = 0.8695*TMath::Erf((x-0.3728)/1.7623);

      else if (idx == 15  ) num = 0.8912*TMath::Erf((x-0.0001)/2.3161);

      else if (idx == 16  ) num = 0.9216*TMath::Erf((x-0.6802)/1.6184);

      else if (idx == 17  ) num = 0.8704*TMath::Erf((x-0.0003)/2.1381);

      else if (idx == 18  ) num = 0.9062*TMath::Erf((x-0.7328)/1.4465);

      else if (idx == 19  ) num = 0.9094*TMath::Erf((x-0.4722)/1.7021);

      else if (idx == 20  ) num = 0.8952*TMath::Erf((x-0.9183)/1.5092);

      else if (idx == 21  ) num = 0.8595*TMath::Erf((x-0.1234)/1.8473);

      else if (idx == 22  ) num = 0.8858*TMath::Erf((x-0.8695)/1.3805);

      else if (idx == 23  ) num = 0.8362*TMath::Erf((x-0.6018)/1.5297);

      else if (idx == 24  ) num = 0.9436*TMath::Erf((x-0.0430)/2.4950);

      else if (idx == 25  ) num = 0.9408*TMath::Erf((x-0.9158)/1.4521);

      else if (idx == 26  ) num = 0.8901*TMath::Erf((x-0.6716)/1.8363);

      else if (idx == 27  ) num = 0.9394*TMath::Erf((x-0.7987)/1.7748);

      else if (idx == 28  ) num = 0.8981*TMath::Erf((x-0.9823)/1.3582);

      else if (idx == 29  ) num = 0.8957*TMath::Erf((x-0.3658)/1.8645);

      else if (idx == 30  ) num = 0.8753*TMath::Erf((x-0.3396)/1.7724);

      else if (idx == 31  ) num = 0.8746*TMath::Erf((x-0.7573)/1.3503);

      else if (idx == 32  ) num = 0.9015*TMath::Erf((x-0.1244)/2.1035);

      else if (idx == 33  ) num = 0.8943*TMath::Erf((x-0.8617)/1.2818);

      else if (idx == 34  ) num = 0.8754*TMath::Erf((x-0.8873)/1.2571);

      else if (idx == 35  ) num = 0.9094*TMath::Erf((x-0.0000)/2.2913);

      else if (idx == 36  ) num = 0.8748*TMath::Erf((x-0.3532)/1.8459);

      else if (idx == 37  ) num = 0.8762*TMath::Erf((x-0.1559)/1.7802);

      else if (idx == 38  ) num = 0.8949*TMath::Erf((x-0.8620)/1.3595);

      else if (idx == 39  ) num = 0.9234*TMath::Erf((x-0.3123)/2.0381);

      else if (idx == 40  ) num = 0.9314*TMath::Erf((x-0.3319)/2.1803);

      else if (idx == 41  ) num = 0.8983*TMath::Erf((x-0.7128)/1.5919);

      else if (idx == 42  ) num = 0.9096*TMath::Erf((x-0.7037)/1.7531);

      else if (idx == 43  ) num = 0.8963*TMath::Erf((x-0.0000)/2.1340);

      else if (idx == 44  ) num = 0.8698*TMath::Erf((x-0.6676)/1.5945);

      else if (idx == 45  ) num = 0.9066*TMath::Erf((x-0.5411)/1.9654);

      else if (idx == 46  ) num = 0.9148*TMath::Erf((x-0.2244)/2.3106);

      else if (idx == 47  ) num = 0.9024*TMath::Erf((x-0.5429)/1.6936);

      else if (idx == 48  ) num = 0.9351*TMath::Erf((x-0.7964)/1.5699);

      else if (idx == 49  ) num = 0.8367*TMath::Erf((x-0.0000)/2.0225);

      else if (idx == 50  ) num = 0.8829*TMath::Erf((x-0.6488)/1.5750);

      else if (idx == 51  ) num = 0.9118*TMath::Erf((x-0.3236)/1.9814);

      else if (idx == 52  ) num = 0.9122*TMath::Erf((x-0.9721)/1.3455);

      else if (idx == 53  ) num = 0.9488*TMath::Erf((x-0.0001)/2.6158);

      else if (idx == 54  ) num = 0.9088*TMath::Erf((x-0.0000)/1.9977);

      else if (idx == 55  ) num = 0.9224*TMath::Erf((x-0.6847)/1.8960);

      else if (idx == 56  ) num = 0.8904*TMath::Erf((x-0.8242)/1.4509);

      else if (idx == 57  ) num = 0.8992*TMath::Erf((x-0.7761)/1.4855);

      else if (idx == 58  ) num = 0.8634*TMath::Erf((x-0.6169)/1.6313);

      else if (idx == 59  ) num = 0.8947*TMath::Erf((x-0.0010)/2.2935);

      else if (idx == 60  ) num = 0.8724*TMath::Erf((x-0.0008)/2.1255);

      else if (idx == 61  ) num = 0.8777*TMath::Erf((x-1.1245)/1.1490);

      else if (idx == 62  ) num = 0.8892*TMath::Erf((x-0.1847)/2.0135);

      else if (idx == 63  ) num = 0.9294*TMath::Erf((x-0.7816)/1.5228);

      else if (idx == 64  ) num = 0.9270*TMath::Erf((x-0.0001)/2.2717);

      else if (idx == 65  ) num = 0.8819*TMath::Erf((x-0.1885)/1.8378);

      else if (idx == 66  ) num = 0.9206*TMath::Erf((x-0.2796)/1.8883);

      else if (idx == 67  ) num = 0.9468*TMath::Erf((x-0.9246)/1.6092);

      else if (idx == 68  ) num = 0.8844*TMath::Erf((x-0.0000)/2.2454);

      else if (idx == 69  ) num = 0.8823*TMath::Erf((x-0.7625)/1.5269);

      else if (idx == 70  ) num = 0.8882*TMath::Erf((x-0.7528)/1.3078);

      else if (idx == 71  ) num = 0.8822*TMath::Erf((x-0.0041)/2.3017);

      else if (idx == 72  ) num = 0.9142*TMath::Erf((x-0.8162)/1.4932);

      else if (idx == 73  ) num = 0.8767*TMath::Erf((x-0.3554)/1.7210);

      else if (idx == 74  ) num = 0.8689*TMath::Erf((x-0.0000)/2.1210);

      else if (idx == 75  ) num = 0.9012*TMath::Erf((x-0.0000)/1.9986);

      else if (idx == 76  ) num = 0.9120*TMath::Erf((x-0.0038)/2.4376);

      else if (idx == 77  ) num = 0.9272*TMath::Erf((x-0.1844)/2.0625);

      else if (idx == 78  ) num = 0.8869*TMath::Erf((x-1.0120)/1.1948);

      else if (idx == 79  ) num = 0.9220*TMath::Erf((x-0.5897)/1.8112);

      else if (idx == 80  ) num = 0.9181*TMath::Erf((x-0.0000)/2.3025);

      else if (idx == 81  ) num = 0.9095*TMath::Erf((x-1.0638)/1.2207);

      else if (idx == 82  ) num = 0.9628*TMath::Erf((x-0.0001)/2.8342);

      else if (idx == 83  ) num = 0.9052*TMath::Erf((x-0.7234)/1.6366);

      else if (idx == 84  ) num = 0.8753*TMath::Erf((x-0.7512)/1.4707);

      else if (idx == 85  ) num = 0.9292*TMath::Erf((x-0.3041)/2.2660);

      else if (idx == 86  ) num = 0.8963*TMath::Erf((x-1.1185)/1.1034);

      else if (idx == 87  ) num = 0.9214*TMath::Erf((x-0.0000)/2.5208);

      else if (idx == 88  ) num = 0.9333*TMath::Erf((x-1.1350)/1.2039);

      else if (idx == 89  ) num = 0.9076*TMath::Erf((x-0.1215)/2.3650);

      else if (idx == 90  ) num = 0.8832*TMath::Erf((x-1.1867)/0.9271);

      else if (idx == 91  ) num = 0.8747*TMath::Erf((x-0.5156)/1.5641);

      else if (idx == 92  ) num = 0.8997*TMath::Erf((x-0.8770)/1.5595);

      else if (idx == 93  ) num = 0.8804*TMath::Erf((x-1.0645)/1.3025);

      else if (idx == 94  ) num = 0.9037*TMath::Erf((x-0.4540)/1.8448);

      else if (idx == 95  ) num = 0.8897*TMath::Erf((x-0.0002)/2.2927);

      else if (idx == 96  ) num = 0.8879*TMath::Erf((x-0.3107)/1.5989);

      else if (idx == 97  ) num = 0.8923*TMath::Erf((x-0.9246)/1.3586);

      else if (idx == 98  ) num = 0.8734*TMath::Erf((x-0.5715)/1.7266);

      else if (idx == 99  ) num = 0.8638*TMath::Erf((x-0.8665)/1.3054);

      else if (idx == 100 ) num = 0.8544*TMath::Erf((x-0.8422)/1.1597);

   }



   // return

   return num/den;

}



/////////////////////////////////////////////////////////////

//                     STA    P P                          //

/////////////////////////////////////////////////////////////

double tnp_weight_sta_pp(double x, double eta, int idx)

{

   // denominator (from MC)

   double den=1;

   if (fabs(eta)<1.6) den = 0.9911*TMath::Erf((x-1.4336)/2.8548);

   else den = den = 0.9132*TMath::Erf((x-0.8045)/1.8366);

   

   // numerator (from data)

   double num=1;

   if (fabs(eta)<1.6)

   {

      if (idx==0) num = 0.9891*TMath::Erf((x-1.4814)/2.5014);

      else if (idx == -1  ) num = 0.9996*TMath::Erf((x-1.4885)/2.4666);

      else if (idx == -2  ) num = 0.9764*TMath::Erf((x-1.4746)/2.5279);

      else if (idx == 1   ) num = 0.9997*TMath::Erf((x-1.2305)/2.8470);

      else if (idx == 2   ) num = 0.9808*TMath::Erf((x-1.6531)/2.2358);

      else if (idx == 3   ) num = 0.9841*TMath::Erf((x-1.5142)/2.5031);

      else if (idx == 4   ) num = 0.9885*TMath::Erf((x-1.5256)/2.4544);

      else if (idx == 5   ) num = 0.9825*TMath::Erf((x-1.4605)/2.4849);

      else if (idx == 6   ) num = 0.9792*TMath::Erf((x-1.4838)/2.4373);

      else if (idx == 7   ) num = 0.9958*TMath::Erf((x-1.5510)/2.4435);

      else if (idx == 8   ) num = 0.9863*TMath::Erf((x-1.5458)/2.4393);

      else if (idx == 9   ) num = 0.9798*TMath::Erf((x-1.5607)/2.3604);

      else if (idx == 10  ) num = 0.9781*TMath::Erf((x-1.4435)/2.4518);

      else if (idx == 11  ) num = 0.9760*TMath::Erf((x-1.5528)/2.4327);

      else if (idx == 12  ) num = 0.9988*TMath::Erf((x-1.4017)/2.6770);

      else if (idx == 13  ) num = 0.9878*TMath::Erf((x-1.5813)/2.3393);

      else if (idx == 14  ) num = 0.9823*TMath::Erf((x-1.5285)/2.4336);

      else if (idx == 15  ) num = 0.9714*TMath::Erf((x-1.5666)/2.3340);

      else if (idx == 16  ) num = 1.0000*TMath::Erf((x-1.5067)/2.5766);

      else if (idx == 17  ) num = 0.9794*TMath::Erf((x-1.5532)/2.4038);

      else if (idx == 18  ) num = 0.9831*TMath::Erf((x-1.6256)/2.2598);

      else if (idx == 19  ) num = 0.9972*TMath::Erf((x-1.3103)/2.7837);

      else if (idx == 20  ) num = 0.9777*TMath::Erf((x-1.6366)/2.2356);

      else if (idx == 21  ) num = 0.9921*TMath::Erf((x-1.4602)/2.5783);

      else if (idx == 22  ) num = 0.9787*TMath::Erf((x-1.5735)/2.2953);

      else if (idx == 23  ) num = 0.9908*TMath::Erf((x-1.5606)/2.3858);

      else if (idx == 24  ) num = 0.9686*TMath::Erf((x-1.5499)/2.2851);

      else if (idx == 25  ) num = 0.9839*TMath::Erf((x-1.5491)/2.3801);

      else if (idx == 26  ) num = 0.9704*TMath::Erf((x-1.5020)/2.3000);

      else if (idx == 27  ) num = 0.9962*TMath::Erf((x-1.5253)/2.5331);

      else if (idx == 28  ) num = 0.9877*TMath::Erf((x-1.4312)/2.5322);

      else if (idx == 29  ) num = 0.9809*TMath::Erf((x-1.4750)/2.5360);

      else if (idx == 30  ) num = 0.9964*TMath::Erf((x-1.3416)/2.7330);

      else if (idx == 31  ) num = 0.9882*TMath::Erf((x-1.3116)/2.7141);

      else if (idx == 32  ) num = 0.9734*TMath::Erf((x-1.4982)/2.4234);

      else if (idx == 33  ) num = 0.9837*TMath::Erf((x-1.5326)/2.3887);

      else if (idx == 34  ) num = 0.9927*TMath::Erf((x-1.3951)/2.6987);

      else if (idx == 35  ) num = 0.9864*TMath::Erf((x-1.4228)/2.5776);

      else if (idx == 36  ) num = 0.9866*TMath::Erf((x-1.4885)/2.5036);

      else if (idx == 37  ) num = 0.9932*TMath::Erf((x-1.4908)/2.5735);

      else if (idx == 38  ) num = 0.9996*TMath::Erf((x-1.3992)/2.6534);

      else if (idx == 39  ) num = 0.9794*TMath::Erf((x-1.5634)/2.3570);

      else if (idx == 40  ) num = 0.9768*TMath::Erf((x-1.4562)/2.4804);

      else if (idx == 41  ) num = 0.9984*TMath::Erf((x-1.4645)/2.5347);

      else if (idx == 42  ) num = 0.9864*TMath::Erf((x-1.3932)/2.6099);

      else if (idx == 43  ) num = 0.9996*TMath::Erf((x-1.4791)/2.5683);

      else if (idx == 44  ) num = 0.9684*TMath::Erf((x-1.4392)/2.4518);

      else if (idx == 45  ) num = 0.9763*TMath::Erf((x-1.4217)/2.4716);

      else if (idx == 46  ) num = 0.9948*TMath::Erf((x-1.4627)/2.5648);

      else if (idx == 47  ) num = 0.9998*TMath::Erf((x-1.3571)/2.6758);

      else if (idx == 48  ) num = 0.9752*TMath::Erf((x-1.5850)/2.3455);

      else if (idx == 49  ) num = 0.9818*TMath::Erf((x-1.4742)/2.4910);

      else if (idx == 50  ) num = 0.9945*TMath::Erf((x-1.4668)/2.5607);

      else if (idx == 51  ) num = 0.9940*TMath::Erf((x-1.4846)/2.5226);

      else if (idx == 52  ) num = 0.9768*TMath::Erf((x-1.6064)/2.2607);

      else if (idx == 53  ) num = 0.9916*TMath::Erf((x-1.3411)/2.6676);

      else if (idx == 54  ) num = 0.9842*TMath::Erf((x-1.5723)/2.3520);

      else if (idx == 55  ) num = 0.9890*TMath::Erf((x-1.5482)/2.3999);

      else if (idx == 56  ) num = 0.9759*TMath::Erf((x-1.4152)/2.5205);

      else if (idx == 57  ) num = 0.9845*TMath::Erf((x-1.4162)/2.5878);

      else if (idx == 58  ) num = 0.9759*TMath::Erf((x-1.7060)/2.1569);

      else if (idx == 59  ) num = 0.9952*TMath::Erf((x-1.2909)/2.7459);

      else if (idx == 60  ) num = 0.9740*TMath::Erf((x-1.5456)/2.4154);

      else if (idx == 61  ) num = 0.9921*TMath::Erf((x-1.2999)/2.7547);

      else if (idx == 62  ) num = 0.9790*TMath::Erf((x-1.5088)/2.5043);

      else if (idx == 63  ) num = 0.9705*TMath::Erf((x-1.5901)/2.3125);

      else if (idx == 64  ) num = 0.9907*TMath::Erf((x-1.2406)/2.7966);

      else if (idx == 65  ) num = 0.9776*TMath::Erf((x-1.5751)/2.3212);

      else if (idx == 66  ) num = 0.9871*TMath::Erf((x-1.5730)/2.3779);

      else if (idx == 67  ) num = 0.9925*TMath::Erf((x-1.5635)/2.4470);

      else if (idx == 68  ) num = 0.9951*TMath::Erf((x-1.4928)/2.5688);

      else if (idx == 69  ) num = 0.9850*TMath::Erf((x-1.4419)/2.5940);

      else if (idx == 70  ) num = 0.9827*TMath::Erf((x-1.4959)/2.4007);

      else if (idx == 71  ) num = 0.9793*TMath::Erf((x-1.6164)/2.2787);

      else if (idx == 72  ) num = 0.9694*TMath::Erf((x-1.4133)/2.4945);

      else if (idx == 73  ) num = 0.9874*TMath::Erf((x-1.4365)/2.5503);

      else if (idx == 74  ) num = 0.9917*TMath::Erf((x-1.3423)/2.6975);

      else if (idx == 75  ) num = 0.9807*TMath::Erf((x-1.4708)/2.4584);

      else if (idx == 76  ) num = 0.9938*TMath::Erf((x-1.4043)/2.6676);

      else if (idx == 77  ) num = 0.9818*TMath::Erf((x-1.5339)/2.4156);

      else if (idx == 78  ) num = 0.9967*TMath::Erf((x-1.4764)/2.5805);

      else if (idx == 79  ) num = 0.9674*TMath::Erf((x-1.4919)/2.4139);

      else if (idx == 80  ) num = 0.9867*TMath::Erf((x-1.4662)/2.4670);

      else if (idx == 81  ) num = 0.9877*TMath::Erf((x-1.6161)/2.3018);

      else if (idx == 82  ) num = 0.9770*TMath::Erf((x-1.6503)/2.1851);

      else if (idx == 83  ) num = 0.9997*TMath::Erf((x-1.3749)/2.6550);

      else if (idx == 84  ) num = 1.0000*TMath::Erf((x-1.4090)/2.6545);

      else if (idx == 85  ) num = 0.9845*TMath::Erf((x-1.3795)/2.5452);

      else if (idx == 86  ) num = 0.9855*TMath::Erf((x-1.5484)/2.3563);

      else if (idx == 87  ) num = 0.9614*TMath::Erf((x-1.5913)/2.2771);

      else if (idx == 88  ) num = 0.9702*TMath::Erf((x-1.5391)/2.3327);

      else if (idx == 89  ) num = 0.9880*TMath::Erf((x-1.3903)/2.6074);

      else if (idx == 90  ) num = 0.9873*TMath::Erf((x-1.5271)/2.4336);

      else if (idx == 91  ) num = 0.9873*TMath::Erf((x-1.5820)/2.3645);

      else if (idx == 92  ) num = 0.9783*TMath::Erf((x-1.4245)/2.5306);

      else if (idx == 93  ) num = 0.9831*TMath::Erf((x-1.4624)/2.4908);

      else if (idx == 94  ) num = 0.9872*TMath::Erf((x-1.4791)/2.5078);

      else if (idx == 95  ) num = 0.9885*TMath::Erf((x-1.5323)/2.3743);

      else if (idx == 96  ) num = 0.9855*TMath::Erf((x-1.5063)/2.4558);

      else if (idx == 97  ) num = 0.9748*TMath::Erf((x-1.5128)/2.3255);

      else if (idx == 98  ) num = 0.9771*TMath::Erf((x-1.4117)/2.5177);

      else if (idx == 99  ) num = 0.9759*TMath::Erf((x-1.5152)/2.3957);

      else if (idx == 100 ) num = 0.9944*TMath::Erf((x-1.4125)/2.5834);

   }

   else

   {

      if (idx==0) num = 0.8956*TMath::Erf((x-0.5162)/1.7646);

      else if (idx == -1  ) num = 0.9153*TMath::Erf((x-0.5049)/1.7401);

      else if (idx == -2  ) num = 0.8762*TMath::Erf((x-0.5250)/1.7942);

      else if (idx == 1   ) num = 0.8598*TMath::Erf((x-0.3794)/1.6775);

      else if (idx == 2   ) num = 0.9071*TMath::Erf((x-0.7765)/1.5351);

      else if (idx == 3   ) num = 0.9070*TMath::Erf((x-0.5316)/1.8221);

      else if (idx == 4   ) num = 0.9148*TMath::Erf((x-0.5456)/1.8124);

      else if (idx == 5   ) num = 0.8999*TMath::Erf((x-0.5368)/1.7840);

      else if (idx == 6   ) num = 0.9028*TMath::Erf((x-0.7388)/1.5306);

      else if (idx == 7   ) num = 0.9112*TMath::Erf((x-0.5761)/1.8778);

      else if (idx == 8   ) num = 0.8828*TMath::Erf((x-0.6303)/1.5635);

      else if (idx == 9   ) num = 0.9021*TMath::Erf((x-0.8599)/1.4212);

      else if (idx == 10  ) num = 0.9214*TMath::Erf((x-0.0992)/2.1977);

      else if (idx == 11  ) num = 0.8966*TMath::Erf((x-0.0000)/2.1930);

      else if (idx == 12  ) num = 0.8840*TMath::Erf((x-0.5171)/1.7189);

      else if (idx == 13  ) num = 0.8872*TMath::Erf((x-0.5657)/1.6991);

      else if (idx == 14  ) num = 0.8801*TMath::Erf((x-0.4565)/1.7340);

      else if (idx == 15  ) num = 0.8986*TMath::Erf((x-0.0001)/2.3079);

      else if (idx == 16  ) num = 0.9161*TMath::Erf((x-0.5990)/1.7086);

      else if (idx == 17  ) num = 0.8810*TMath::Erf((x-0.1226)/2.0724);

      else if (idx == 18  ) num = 0.9073*TMath::Erf((x-0.5643)/1.7015);

      else if (idx == 19  ) num = 0.9056*TMath::Erf((x-0.4883)/1.7268);

      else if (idx == 20  ) num = 0.9011*TMath::Erf((x-0.7915)/1.6093);

      else if (idx == 21  ) num = 0.8770*TMath::Erf((x-0.2307)/1.8953);

      else if (idx == 22  ) num = 0.8914*TMath::Erf((x-0.6992)/1.5914);

      else if (idx == 23  ) num = 0.8560*TMath::Erf((x-0.6063)/1.5601);

      else if (idx == 24  ) num = 0.9262*TMath::Erf((x-0.2751)/2.1548);

      else if (idx == 25  ) num = 0.9307*TMath::Erf((x-0.7603)/1.6270);

      else if (idx == 26  ) num = 0.8949*TMath::Erf((x-0.6696)/1.7413);

      else if (idx == 27  ) num = 0.9239*TMath::Erf((x-0.7273)/1.7403);

      else if (idx == 28  ) num = 0.8983*TMath::Erf((x-0.8582)/1.4675);

      else if (idx == 29  ) num = 0.8970*TMath::Erf((x-0.3959)/1.8485);

      else if (idx == 30  ) num = 0.8809*TMath::Erf((x-0.4183)/1.7454);

      else if (idx == 31  ) num = 0.8817*TMath::Erf((x-0.6713)/1.5019);

      else if (idx == 32  ) num = 0.8976*TMath::Erf((x-0.2735)/1.9515);

      else if (idx == 33  ) num = 0.8952*TMath::Erf((x-0.7057)/1.5066);

      else if (idx == 34  ) num = 0.8842*TMath::Erf((x-0.7611)/1.4477);

      else if (idx == 35  ) num = 0.9058*TMath::Erf((x-0.0001)/2.3016);

      else if (idx == 36  ) num = 0.8841*TMath::Erf((x-0.4137)/1.8044);

      else if (idx == 37  ) num = 0.8802*TMath::Erf((x-0.2550)/1.7980);

      else if (idx == 38  ) num = 0.8932*TMath::Erf((x-0.7425)/1.5000);

      else if (idx == 39  ) num = 0.9204*TMath::Erf((x-0.2807)/2.0567);

      else if (idx == 40  ) num = 0.9197*TMath::Erf((x-0.4595)/1.9504);

      else if (idx == 41  ) num = 0.8961*TMath::Erf((x-0.6696)/1.6281);

      else if (idx == 42  ) num = 0.9087*TMath::Erf((x-0.6564)/1.7282);

      else if (idx == 43  ) num = 0.9029*TMath::Erf((x-0.0140)/2.2072);

      else if (idx == 44  ) num = 0.8831*TMath::Erf((x-0.6211)/1.6572);

      else if (idx == 45  ) num = 0.9050*TMath::Erf((x-0.5887)/1.8157);

      else if (idx == 46  ) num = 0.9094*TMath::Erf((x-0.3107)/2.1052);

      else if (idx == 47  ) num = 0.9037*TMath::Erf((x-0.4793)/1.8075);

      else if (idx == 48  ) num = 0.9198*TMath::Erf((x-0.7007)/1.6440);

      else if (idx == 49  ) num = 0.8651*TMath::Erf((x-0.0000)/2.1206);

      else if (idx == 50  ) num = 0.8829*TMath::Erf((x-0.6340)/1.5785);

      else if (idx == 51  ) num = 0.9117*TMath::Erf((x-0.3544)/1.9652);

      else if (idx == 52  ) num = 0.9085*TMath::Erf((x-0.8290)/1.4986);

      else if (idx == 53  ) num = 0.9339*TMath::Erf((x-0.0001)/2.5367);

      else if (idx == 54  ) num = 0.9055*TMath::Erf((x-0.0051)/2.1052);

      else if (idx == 55  ) num = 0.9139*TMath::Erf((x-0.6377)/1.8184);

      else if (idx == 56  ) num = 0.8956*TMath::Erf((x-0.7133)/1.5832);

      else if (idx == 57  ) num = 0.8944*TMath::Erf((x-0.6918)/1.5760);

      else if (idx == 58  ) num = 0.8740*TMath::Erf((x-0.6102)/1.6201);

      else if (idx == 59  ) num = 0.8934*TMath::Erf((x-0.2343)/2.0419);

      else if (idx == 60  ) num = 0.8862*TMath::Erf((x-0.1511)/2.0487);

      else if (idx == 61  ) num = 0.8825*TMath::Erf((x-0.9751)/1.2993);

      else if (idx == 62  ) num = 0.8925*TMath::Erf((x-0.2182)/2.0271);

      else if (idx == 63  ) num = 0.9230*TMath::Erf((x-0.6998)/1.6171);

      else if (idx == 64  ) num = 0.9294*TMath::Erf((x-0.0001)/2.3928);

      else if (idx == 65  ) num = 0.8908*TMath::Erf((x-0.2356)/1.9359);

      else if (idx == 66  ) num = 0.9198*TMath::Erf((x-0.2721)/1.9908);

      else if (idx == 67  ) num = 0.9297*TMath::Erf((x-0.8119)/1.6569);

      else if (idx == 68  ) num = 0.8875*TMath::Erf((x-0.0197)/2.2234);

      else if (idx == 69  ) num = 0.8853*TMath::Erf((x-0.7150)/1.5560);

      else if (idx == 70  ) num = 0.8926*TMath::Erf((x-0.6122)/1.5496);

      else if (idx == 71  ) num = 0.8866*TMath::Erf((x-0.1591)/2.1265);

      else if (idx == 72  ) num = 0.9078*TMath::Erf((x-0.7423)/1.5587);

      else if (idx == 73  ) num = 0.8853*TMath::Erf((x-0.3719)/1.7744);

      else if (idx == 74  ) num = 0.8848*TMath::Erf((x-0.0002)/2.2106);

      else if (idx == 75  ) num = 0.9000*TMath::Erf((x-0.0000)/2.1305);

      else if (idx == 76  ) num = 0.9025*TMath::Erf((x-0.2675)/2.0987);

      else if (idx == 77  ) num = 0.9209*TMath::Erf((x-0.2320)/2.0596);

      else if (idx == 78  ) num = 0.8915*TMath::Erf((x-0.8565)/1.3909);

      else if (idx == 79  ) num = 0.9107*TMath::Erf((x-0.5860)/1.7424);

      else if (idx == 80  ) num = 0.9137*TMath::Erf((x-0.0004)/2.3112);

      else if (idx == 81  ) num = 0.9028*TMath::Erf((x-0.9218)/1.3760);

      else if (idx == 82  ) num = 0.9411*TMath::Erf((x-0.0002)/2.6224);

      else if (idx == 83  ) num = 0.9021*TMath::Erf((x-0.6618)/1.6607);

      else if (idx == 84  ) num = 0.8834*TMath::Erf((x-0.6852)/1.5732);

      else if (idx == 85  ) num = 0.9159*TMath::Erf((x-0.3801)/2.0638);

      else if (idx == 86  ) num = 0.8949*TMath::Erf((x-0.9764)/1.2627);

      else if (idx == 87  ) num = 0.9155*TMath::Erf((x-0.0003)/2.4445);

      else if (idx == 88  ) num = 0.9227*TMath::Erf((x-0.9742)/1.3795);

      else if (idx == 89  ) num = 0.9004*TMath::Erf((x-0.4428)/1.9254);

      else if (idx == 90  ) num = 0.8901*TMath::Erf((x-0.9824)/1.2106);

      else if (idx == 91  ) num = 0.8849*TMath::Erf((x-0.5116)/1.6491);

      else if (idx == 92  ) num = 0.8991*TMath::Erf((x-0.7910)/1.5912);

      else if (idx == 93  ) num = 0.8842*TMath::Erf((x-0.9408)/1.3832);

      else if (idx == 94  ) num = 0.9000*TMath::Erf((x-0.5123)/1.7608);

      else if (idx == 95  ) num = 0.8953*TMath::Erf((x-0.1440)/2.1545);

      else if (idx == 96  ) num = 0.8901*TMath::Erf((x-0.3508)/1.7106);

      else if (idx == 97  ) num = 0.8954*TMath::Erf((x-0.8025)/1.5071);

      else if (idx == 98  ) num = 0.8822*TMath::Erf((x-0.5702)/1.6965);

      else if (idx == 99  ) num = 0.8738*TMath::Erf((x-0.7460)/1.4603);

      else if (idx == 100 ) num = 0.8709*TMath::Erf((x-0.7599)/1.3440);

   }



   // return

   return num/den;

}



#endif //#ifndef tnp_weight_h