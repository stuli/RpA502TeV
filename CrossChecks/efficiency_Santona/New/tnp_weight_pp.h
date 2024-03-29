#ifndef tnp_weight_h
#define tnp_weight_h

#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++
//
// - Trigger: (tnp_weight_trg_*)
//   * idx = 0:  nominal
//   * idx = 1..100: toy variations, stat. only
//   * idx = -1: syst variation, "new_MAX", +1 sigma
//   * idx = -2: syst variation, "new_MAX", -1 sigma
//   * idx = -10: binned
// - MuID, STA: (tnp_weight_muid_*, tnp_weight_sta_pp)
//   * same organisation. BUT these SF are for systematics only.
// - Inner tracking: 
//   * idx = 0: nominal correction (pt, eta independent)
//   * idx = -1: +1 sigma
//   * idx = -2: -1 sigma
//
// x is pT
//
// Corrections are needed for trg and trk only. For muid and sta are needed for systematics only. 


// THE INDIVIDUAL SFs
// ++++++++++++++++++
double tnp_weight_trg_pp(double x, double eta, int idx=0);
double tnp_weight_muid_pp(double x, double eta, int idx=0);
double tnp_weight_sta_pp(double x, double eta, int idx=0);
double tnp_weight_trk_pp(int idx);


///////////////////////////////////////////////////
//               T R G    P P                    //
///////////////////////////////////////////////////
double tnp_weight_trg_pp(double x, double eta, int idx)
{
   // binned
   if (idx == -10) {
      if (fabs(eta) < 1.2) {
         // 0 < |eta| < 1.2
         if (x<4) return 1.05133;
         else if (x<4.5) return 1.01493;
         else if (x<5) return 1.00637;
         else if (x<5.5) return 1.00866;
         else if (x<6) return 1.00403;
         else if (x<7) return 0.999526;
         else if (x<8) return 0.993882;
         else if (x<10.5) return 0.999195;
         else if (x<14) return 0.994448;
         else if (x<18) return 0.997814;
         else return 0.99535;
      } else if (fabs(eta)<1.8) {
         // 1.2 < |eta| < 1.8
         if (x<3) return 1.20399;
         else if (x<3.5) return 1.08496;
         else if (x<4) return 1.0579;
         else if (x<4.5) return 1.02943;
         else if (x<5) return 1.02136;
         else if (x<6) return 1.01855;
         else if (x<7) return 1.01709;
         else if (x<9) return 1.00907;
         else if (x<14) return 1.00748;
         else if (x<18) return 1.02484;
         else return 1.02539;
      } else if (fabs(eta)<2.1) {
         // 1.8 < |eta| < 2.1
         if (x<2.5) return 1.23506;
         else if (x<3) return 1.07046;
         else if (x<3.5) return 1.05211;
         else if (x<4) return 1.04089;
         else if (x<4.5) return 1.04115;
         else if (x<5) return 1.03067;
         else if (x<6) return 1.02442;
         else if (x<7) return 1.025;
         else if (x<9) return 1.02442;
         else if (x<13) return 0.997836;
         else return 0.944194;
      } else {
         // 2.1 < |eta| < 2.4
         if (x<2.5) return 1.21185;
         else if (x<3) return 1.10746;
         else if (x<3.5) return 1.07255;
         else if (x<4) return 1.03562;
         else if (x<4.5) return 1.03939;
         else if (x<5) return 1.02434;
         else if (x<6.5) return 1.02967;
         else if (x<8.5) return 1.02523;
         else if (x<12) return 1.03297;
         else return 0.990749;
      }
   }

   // denominator (from MC)
   double den=1;
   if (fabs(eta)<1.2) den = (0.51387*(TMath::Erf((x-2.95780)/1.01717)*TMath::Exp(-0.00001*x)))+0.46186;
   else if (fabs(eta)<1.8) den = (0.87963*(TMath::Erf((x-2.33663)/1.22198)*TMath::Exp(-0.00013*x)))+0.06241;
   else if (fabs(eta)<2.1) den = (0.53171*(TMath::Erf((x-1.80661)/1.67849)*TMath::Exp(-0.00000*x)))+0.36279;
   else den = (0.38994*(TMath::Erf((x-1.55863)/1.70815)*TMath::Exp(-0.00003*x)))+0.48718;

   // numerator (from data)
   double num=1;
   if (fabs(eta)<1.2)
   {
      if (idx==0) num = (0.52580*(TMath::Erf((x-2.68070)/1.16799)*TMath::Exp(-0.00234*x)))+0.46172;
      else if (idx == -1   ) num = (0.52587*(TMath::Erf((x-2.68165)/1.16925)*TMath::Exp(-0.00229*x)))+0.46179;
      else if (idx == -2   ) num = (0.52570*(TMath::Erf((x-2.67968)/1.16728)*TMath::Exp(-0.00239*x)))+0.46163;
      else if (idx == 1   ) num = (0.52452*(TMath::Erf((x-2.35489)/1.45504)*TMath::Exp(-0.00197*x)))+0.46344;
      else if (idx == 2   ) num = (0.52934*(TMath::Erf((x-2.77208)/1.09665)*TMath::Exp(-0.00408*x)))+0.46478;
      else if (idx == 3   ) num = (0.52466*(TMath::Erf((x-2.68561)/1.17908)*TMath::Exp(-0.00204*x)))+0.46070;
      else if (idx == 4   ) num = (0.52281*(TMath::Erf((x-2.75691)/1.07074)*TMath::Exp(-0.00159*x)))+0.46007;
      else if (idx == 5   ) num = (0.52514*(TMath::Erf((x-2.67519)/1.16231)*TMath::Exp(-0.00247*x)))+0.46113;
      else if (idx == 6   ) num = (0.52817*(TMath::Erf((x-2.59227)/1.29525)*TMath::Exp(-0.00295*x)))+0.46265;
      else if (idx == 7   ) num = (0.52828*(TMath::Erf((x-2.69555)/1.17103)*TMath::Exp(-0.00370*x)))+0.46376;
      else if (idx == 8   ) num = (0.52621*(TMath::Erf((x-2.68257)/1.16886)*TMath::Exp(-0.00226*x)))+0.46210;
      else if (idx == 9   ) num = (0.52731*(TMath::Erf((x-2.59959)/1.28505)*TMath::Exp(-0.00239*x)))+0.46214;
      else if (idx == 10  ) num = (0.52414*(TMath::Erf((x-2.67112)/1.16549)*TMath::Exp(-0.00216*x)))+0.46028;
      else if (idx == 11  ) num = (0.52552*(TMath::Erf((x-2.60357)/1.24309)*TMath::Exp(-0.00193*x)))+0.46096;
      else if (idx == 12  ) num = (0.52633*(TMath::Erf((x-2.52728)/1.31946)*TMath::Exp(-0.00227*x)))+0.46254;
      else if (idx == 13  ) num = (0.52310*(TMath::Erf((x-2.66254)/1.16432)*TMath::Exp(-0.00128*x)))+0.45950;
      else if (idx == 14  ) num = (0.52642*(TMath::Erf((x-2.68076)/1.16826)*TMath::Exp(-0.00263*x)))+0.46224;
      else if (idx == 15  ) num = (0.52526*(TMath::Erf((x-2.75226)/1.09164)*TMath::Exp(-0.00233*x)))+0.46198;
      else if (idx == 16  ) num = (0.52919*(TMath::Erf((x-2.76695)/1.09741)*TMath::Exp(-0.00407*x)))+0.46459;
      else if (idx == 17  ) num = (0.52728*(TMath::Erf((x-2.68522)/1.16877)*TMath::Exp(-0.00304*x)))+0.46296;
      else if (idx == 18  ) num = (0.52206*(TMath::Erf((x-2.81136)/1.04090)*TMath::Exp(-0.00110*x)))+0.45924;
      else if (idx == 19  ) num = (0.52457*(TMath::Erf((x-2.66344)/1.16210)*TMath::Exp(-0.00193*x)))+0.46075;
      else if (idx == 20  ) num = (0.52589*(TMath::Erf((x-2.67629)/1.16444)*TMath::Exp(-0.00235*x)))+0.46182;
      else if (idx == 21  ) num = (0.52703*(TMath::Erf((x-2.58564)/1.28918)*TMath::Exp(-0.00224*x)))+0.46212;
      else if (idx == 22  ) num = (0.52829*(TMath::Erf((x-2.61746)/1.23702)*TMath::Exp(-0.00361*x)))+0.46300;
      else if (idx == 23  ) num = (0.52490*(TMath::Erf((x-2.67810)/1.16528)*TMath::Exp(-0.00194*x)))+0.46095;
      else if (idx == 24  ) num = (0.52293*(TMath::Erf((x-2.66434)/1.18072)*TMath::Exp(-0.00111*x)))+0.45940;
      else if (idx == 25  ) num = (0.52724*(TMath::Erf((x-2.84878)/1.00287)*TMath::Exp(-0.00332*x)))+0.46312;
      else if (idx == 26  ) num = (0.52744*(TMath::Erf((x-2.70500)/1.17644)*TMath::Exp(-0.00298*x)))+0.46300;
      else if (idx == 27  ) num = (0.52612*(TMath::Erf((x-2.68385)/1.17128)*TMath::Exp(-0.00226*x)))+0.46201;
      else if (idx == 28  ) num = (0.52708*(TMath::Erf((x-2.44817)/1.41574)*TMath::Exp(-0.00302*x)))+0.46477;
      else if (idx == 29  ) num = (0.52813*(TMath::Erf((x-2.60511)/1.25966)*TMath::Exp(-0.00328*x)))+0.46273;
      else if (idx == 30  ) num = (0.52656*(TMath::Erf((x-2.45324)/1.41476)*TMath::Exp(-0.00260*x)))+0.46426;
      else if (idx == 31  ) num = (0.52394*(TMath::Erf((x-2.54061)/1.29309)*TMath::Exp(-0.00132*x)))+0.45918;
      else if (idx == 32  ) num = (0.52373*(TMath::Erf((x-2.67145)/1.17424)*TMath::Exp(-0.00145*x)))+0.46001;
      else if (idx == 33  ) num = (0.52636*(TMath::Erf((x-2.76497)/1.08028)*TMath::Exp(-0.00272*x)))+0.46277;
      else if (idx == 34  ) num = (0.52462*(TMath::Erf((x-2.67079)/1.16798)*TMath::Exp(-0.00160*x)))+0.46079;
      else if (idx == 35  ) num = (0.52286*(TMath::Erf((x-2.81377)/1.00901)*TMath::Exp(-0.00171*x)))+0.46024;
      else if (idx == 36  ) num = (0.52523*(TMath::Erf((x-2.66794)/1.16193)*TMath::Exp(-0.00199*x)))+0.46131;
      else if (idx == 37  ) num = (0.52312*(TMath::Erf((x-2.74102)/1.09489)*TMath::Exp(-0.00139*x)))+0.46015;
      else if (idx == 38  ) num = (0.52600*(TMath::Erf((x-2.68107)/1.16878)*TMath::Exp(-0.00232*x)))+0.46190;
      else if (idx == 39  ) num = (0.52232*(TMath::Erf((x-2.65428)/1.15834)*TMath::Exp(-0.00099*x)))+0.45888;
      else if (idx == 40  ) num = (0.52700*(TMath::Erf((x-2.68694)/1.17244)*TMath::Exp(-0.00310*x)))+0.46271;
      else if (idx == 41  ) num = (0.52538*(TMath::Erf((x-2.82397)/1.00989)*TMath::Exp(-0.00267*x)))+0.46217;
      else if (idx == 42  ) num = (0.52622*(TMath::Erf((x-2.69059)/1.17828)*TMath::Exp(-0.00223*x)))+0.46208;
      else if (idx == 43  ) num = (0.52754*(TMath::Erf((x-2.77303)/1.08091)*TMath::Exp(-0.00316*x)))+0.46369;
      else if (idx == 44  ) num = (0.52447*(TMath::Erf((x-2.68015)/1.18223)*TMath::Exp(-0.00190*x)))+0.46059;
      else if (idx == 45  ) num = (0.52586*(TMath::Erf((x-2.67962)/1.16711)*TMath::Exp(-0.00233*x)))+0.46178;
      else if (idx == 46  ) num = (0.52551*(TMath::Erf((x-2.67197)/1.15913)*TMath::Exp(-0.00248*x)))+0.46147;
      else if (idx == 47  ) num = (0.52278*(TMath::Erf((x-2.83594)/0.99002)*TMath::Exp(-0.00115*x)))+0.45945;
      else if (idx == 48  ) num = (0.52288*(TMath::Erf((x-2.73721)/1.10714)*TMath::Exp(-0.00090*x)))+0.45985;
      else if (idx == 49  ) num = (0.52140*(TMath::Erf((x-2.64424)/1.16427)*TMath::Exp(-0.00055*x)))+0.45820;
      else if (idx == 50  ) num = (0.52635*(TMath::Erf((x-2.75543)/1.10865)*TMath::Exp(-0.00262*x)))+0.46234;
      else if (idx == 51  ) num = (0.52330*(TMath::Erf((x-2.86218)/0.95880)*TMath::Exp(-0.00166*x)))+0.45946;
      else if (idx == 52  ) num = (0.52373*(TMath::Erf((x-2.86971)/0.97688)*TMath::Exp(-0.00127*x)))+0.45826;
      else if (idx == 53  ) num = (0.52422*(TMath::Erf((x-2.74165)/1.09627)*TMath::Exp(-0.00186*x)))+0.46085;
      else if (idx == 54  ) num = (0.52561*(TMath::Erf((x-2.68303)/1.17037)*TMath::Exp(-0.00238*x)))+0.46153;
      else if (idx == 55  ) num = (0.52823*(TMath::Erf((x-2.57200)/1.30657)*TMath::Exp(-0.00304*x)))+0.46310;
      else if (idx == 56  ) num = (0.52560*(TMath::Erf((x-2.47978)/1.36119)*TMath::Exp(-0.00226*x)))+0.46240;
      else if (idx == 57  ) num = (0.52256*(TMath::Erf((x-2.65770)/1.16469)*TMath::Exp(-0.00145*x)))+0.45901;
      else if (idx == 58  ) num = (0.52987*(TMath::Erf((x-2.57489)/1.30466)*TMath::Exp(-0.00352*x)))+0.46417;
      else if (idx == 59  ) num = (0.52598*(TMath::Erf((x-2.69172)/1.17819)*TMath::Exp(-0.00227*x)))+0.46184;
      else if (idx == 60  ) num = (0.52207*(TMath::Erf((x-2.79244)/1.01993)*TMath::Exp(-0.00114*x)))+0.45945;
      else if (idx == 61  ) num = (0.52589*(TMath::Erf((x-2.69265)/1.17934)*TMath::Exp(-0.00230*x)))+0.46176;
      else if (idx == 62  ) num = (0.52564*(TMath::Erf((x-2.69262)/1.17797)*TMath::Exp(-0.00198*x)))+0.46156;
      else if (idx == 63  ) num = (0.52588*(TMath::Erf((x-2.49548)/1.35761)*TMath::Exp(-0.00200*x)))+0.46150;
      else if (idx == 64  ) num = (0.52445*(TMath::Erf((x-2.67730)/1.16713)*TMath::Exp(-0.00197*x)))+0.46054;
      else if (idx == 65  ) num = (0.52580*(TMath::Erf((x-2.44552)/1.40849)*TMath::Exp(-0.00258*x)))+0.46364;
      else if (idx == 66  ) num = (0.52354*(TMath::Erf((x-2.87352)/0.94759)*TMath::Exp(-0.00132*x)))+0.45788;
      else if (idx == 67  ) num = (0.52817*(TMath::Erf((x-2.69552)/1.17103)*TMath::Exp(-0.00320*x)))+0.46370;
      else if (idx == 68  ) num = (0.52298*(TMath::Erf((x-2.78556)/1.07754)*TMath::Exp(-0.00140*x)))+0.46030;
      else if (idx == 69  ) num = (0.52216*(TMath::Erf((x-2.83751)/0.98565)*TMath::Exp(-0.00126*x)))+0.45871;
      else if (idx == 70  ) num = (0.52433*(TMath::Erf((x-2.67610)/1.17537)*TMath::Exp(-0.00189*x)))+0.46048;
      else if (idx == 71  ) num = (0.52754*(TMath::Erf((x-2.41440)/1.46757)*TMath::Exp(-0.00329*x)))+0.46541;
      else if (idx == 72  ) num = (0.52696*(TMath::Erf((x-2.59618)/1.26238)*TMath::Exp(-0.00291*x)))+0.46181;
      else if (idx == 73  ) num = (0.52646*(TMath::Erf((x-2.48416)/1.35982)*TMath::Exp(-0.00286*x)))+0.46385;
      else if (idx == 74  ) num = (0.52460*(TMath::Erf((x-2.67009)/1.15880)*TMath::Exp(-0.00206*x)))+0.46070;
      else if (idx == 75  ) num = (0.52595*(TMath::Erf((x-2.58353)/1.26549)*TMath::Exp(-0.00224*x)))+0.46110;
      else if (idx == 76  ) num = (0.52562*(TMath::Erf((x-2.50314)/1.34485)*TMath::Exp(-0.00173*x)))+0.46138;
      else if (idx == 77  ) num = (0.52198*(TMath::Erf((x-2.86934)/0.94978)*TMath::Exp(-0.00118*x)))+0.45928;
      else if (idx == 78  ) num = (0.52220*(TMath::Erf((x-2.65921)/1.17117)*TMath::Exp(-0.00077*x)))+0.45880;
      else if (idx == 79  ) num = (0.52313*(TMath::Erf((x-2.65808)/1.15260)*TMath::Exp(-0.00159*x)))+0.45949;
      else if (idx == 80  ) num = (0.52668*(TMath::Erf((x-2.59023)/1.26609)*TMath::Exp(-0.00244*x)))+0.46165;
      else if (idx == 81  ) num = (0.52635*(TMath::Erf((x-2.67755)/1.16111)*TMath::Exp(-0.00264*x)))+0.46218;
      else if (idx == 82  ) num = (0.52602*(TMath::Erf((x-2.68928)/1.17250)*TMath::Exp(-0.00255*x)))+0.46185;
      else if (idx == 83  ) num = (0.52582*(TMath::Erf((x-2.68361)/1.17183)*TMath::Exp(-0.00230*x)))+0.46173;
      else if (idx == 84  ) num = (0.52293*(TMath::Erf((x-2.66843)/1.16907)*TMath::Exp(-0.00140*x)))+0.45930;
      else if (idx == 85  ) num = (0.52434*(TMath::Erf((x-2.67209)/1.17375)*TMath::Exp(-0.00194*x)))+0.46050;
      else if (idx == 86  ) num = (0.52803*(TMath::Erf((x-2.59703)/1.28077)*TMath::Exp(-0.00253*x)))+0.46249;
      else if (idx == 87  ) num = (0.52411*(TMath::Erf((x-2.77073)/1.06247)*TMath::Exp(-0.00173*x)))+0.46099;
      else if (idx == 88  ) num = (0.52422*(TMath::Erf((x-2.51623)/1.31767)*TMath::Exp(-0.00088*x)))+0.46003;
      else if (idx == 89  ) num = (0.53135*(TMath::Erf((x-2.36660)/1.53977)*TMath::Exp(-0.00494*x)))+0.46920;
      else if (idx == 90  ) num = (0.52736*(TMath::Erf((x-2.59354)/1.24906)*TMath::Exp(-0.00274*x)))+0.46224;
      else if (idx == 91  ) num = (0.52696*(TMath::Erf((x-2.60887)/1.23753)*TMath::Exp(-0.00236*x)))+0.46194;
      else if (idx == 92  ) num = (0.52292*(TMath::Erf((x-2.82867)/1.00100)*TMath::Exp(-0.00210*x)))+0.46016;
      else if (idx == 93  ) num = (0.52697*(TMath::Erf((x-2.54916)/1.29555)*TMath::Exp(-0.00282*x)))+0.46144;
      else if (idx == 94  ) num = (0.52955*(TMath::Erf((x-2.52152)/1.33797)*TMath::Exp(-0.00342*x)))+0.46417;
      else if (idx == 95  ) num = (0.52509*(TMath::Erf((x-2.77296)/1.07610)*TMath::Exp(-0.00201*x)))+0.46190;
      else if (idx == 96  ) num = (0.52587*(TMath::Erf((x-2.54922)/1.28166)*TMath::Exp(-0.00232*x)))+0.46101;
      else if (idx == 97  ) num = (0.52282*(TMath::Erf((x-2.66902)/1.16301)*TMath::Exp(-0.00108*x)))+0.45924;
      else if (idx == 98  ) num = (0.52473*(TMath::Erf((x-2.67440)/1.17207)*TMath::Exp(-0.00173*x)))+0.46087;
      else if (idx == 99  ) num = (0.52958*(TMath::Erf((x-2.55928)/1.29032)*TMath::Exp(-0.00399*x)))+0.46380;
      else if (idx == 100 ) num = (0.52348*(TMath::Erf((x-2.89211)/0.94363)*TMath::Exp(-0.00106*x)))+0.45791;
   }
   else if (fabs(eta)<1.8)
   {
      if (idx==0) num = (0.91504*(TMath::Erf((x-2.23999)/1.16962)*TMath::Exp(-0.00000*x)))+0.04202;
      else if (idx == -1   ) num = (0.91559*(TMath::Erf((x-2.24020)/1.16984)*TMath::Exp(-0.00000*x)))+0.04248;
      else if (idx == -2   ) num = (0.91450*(TMath::Erf((x-2.23958)/1.16980)*TMath::Exp(-0.00000*x)))+0.04156;
      else if (idx == 1   ) num = (0.91624*(TMath::Erf((x-2.21327)/1.23381)*TMath::Exp(-0.00069*x)))+0.04514;
      else if (idx == 2   ) num = (0.91450*(TMath::Erf((x-2.23961)/1.16855)*TMath::Exp(-0.00000*x)))+0.04154;
      else if (idx == 3   ) num = (0.91504*(TMath::Erf((x-2.22873)/1.18396)*TMath::Exp(-0.00000*x)))+0.04261;
      else if (idx == 4   ) num = (0.91564*(TMath::Erf((x-2.24467)/1.17760)*TMath::Exp(-0.00000*x)))+0.04234;
      else if (idx == 5   ) num = (0.91569*(TMath::Erf((x-2.23605)/1.16241)*TMath::Exp(-0.00055*x)))+0.04270;
      else if (idx == 6   ) num = (0.91614*(TMath::Erf((x-2.26412)/1.13703)*TMath::Exp(-0.00000*x)))+0.04080;
      else if (idx == 7   ) num = (0.91508*(TMath::Erf((x-2.24315)/1.17493)*TMath::Exp(-0.00000*x)))+0.04191;
      else if (idx == 8   ) num = (0.91477*(TMath::Erf((x-2.25059)/1.15967)*TMath::Exp(-0.00000*x)))+0.04084;
      else if (idx == 9   ) num = (0.91606*(TMath::Erf((x-2.24308)/1.17389)*TMath::Exp(-0.00000*x)))+0.04278;
      else if (idx == 10  ) num = (0.91281*(TMath::Erf((x-2.22255)/1.17429)*TMath::Exp(-0.00008*x)))+0.04093;
      else if (idx == 11  ) num = (0.91419*(TMath::Erf((x-2.26934)/1.10729)*TMath::Exp(-0.00001*x)))+0.03862;
      else if (idx == 12  ) num = (0.91426*(TMath::Erf((x-2.22107)/1.19868)*TMath::Exp(-0.00000*x)))+0.04264;
      else if (idx == 13  ) num = (0.91481*(TMath::Erf((x-2.22626)/1.18360)*TMath::Exp(-0.00000*x)))+0.04253;
      else if (idx == 14  ) num = (0.91491*(TMath::Erf((x-2.23372)/1.19243)*TMath::Exp(-0.00000*x)))+0.04222;
      else if (idx == 15  ) num = (0.91538*(TMath::Erf((x-2.24101)/1.17106)*TMath::Exp(-0.00000*x)))+0.04227;
      else if (idx == 16  ) num = (0.91761*(TMath::Erf((x-2.24399)/1.16619)*TMath::Exp(-0.00115*x)))+0.04394;
      else if (idx == 17  ) num = (0.91343*(TMath::Erf((x-2.22052)/1.19040)*TMath::Exp(-0.00000*x)))+0.04144;
      else if (idx == 18  ) num = (0.91516*(TMath::Erf((x-2.24204)/1.17186)*TMath::Exp(-0.00019*x)))+0.04202;
      else if (idx == 19  ) num = (0.94768*(TMath::Erf((x-2.18250)/1.22848)*TMath::Exp(-0.00000*x)))+0.01365;
      else if (idx == 20  ) num = (0.91407*(TMath::Erf((x-2.24078)/1.16954)*TMath::Exp(-0.00017*x)))+0.04099;
      else if (idx == 21  ) num = (0.91554*(TMath::Erf((x-2.25665)/1.14497)*TMath::Exp(-0.00000*x)))+0.04086;
      else if (idx == 22  ) num = (0.91499*(TMath::Erf((x-2.22687)/1.18002)*TMath::Exp(-0.00000*x)))+0.04254;
      else if (idx == 23  ) num = (0.91553*(TMath::Erf((x-2.26896)/1.11767)*TMath::Exp(-0.00028*x)))+0.03978;
      else if (idx == 24  ) num = (0.91773*(TMath::Erf((x-2.26930)/1.13960)*TMath::Exp(-0.00000*x)))+0.04181;
      else if (idx == 25  ) num = (0.91575*(TMath::Erf((x-2.22574)/1.20824)*TMath::Exp(-0.00000*x)))+0.04364;
      else if (idx == 26  ) num = (0.91746*(TMath::Erf((x-2.24844)/1.17975)*TMath::Exp(-0.00000*x)))+0.04385;
      else if (idx == 27  ) num = (0.91822*(TMath::Erf((x-2.25975)/1.17822)*TMath::Exp(-0.00061*x)))+0.04363;
      else if (idx == 28  ) num = (0.91996*(TMath::Erf((x-2.23135)/1.20719)*TMath::Exp(-0.00148*x)))+0.04720;
      else if (idx == 29  ) num = (0.91598*(TMath::Erf((x-2.27156)/1.12497)*TMath::Exp(-0.00003*x)))+0.04016;
      else if (idx == 30  ) num = (0.91693*(TMath::Erf((x-2.20717)/1.21986)*TMath::Exp(-0.00000*x)))+0.04609;
      else if (idx == 31  ) num = (0.91593*(TMath::Erf((x-2.26413)/1.13096)*TMath::Exp(-0.00000*x)))+0.04060;
      else if (idx == 32  ) num = (0.91661*(TMath::Erf((x-2.22189)/1.20778)*TMath::Exp(-0.00000*x)))+0.04468;
      else if (idx == 33  ) num = (0.91499*(TMath::Erf((x-2.23510)/1.16195)*TMath::Exp(-0.00020*x)))+0.04217;
      else if (idx == 34  ) num = (0.91585*(TMath::Erf((x-2.26657)/1.13471)*TMath::Exp(-0.00000*x)))+0.04051;
      else if (idx == 35  ) num = (0.91645*(TMath::Erf((x-2.23531)/1.20068)*TMath::Exp(-0.00000*x)))+0.04363;
      else if (idx == 36  ) num = (0.91473*(TMath::Erf((x-2.24554)/1.15959)*TMath::Exp(-0.00000*x)))+0.04112;
      else if (idx == 37  ) num = (0.93865*(TMath::Erf((x-2.23319)/1.16246)*TMath::Exp(-0.00042*x)))+0.02170;
      else if (idx == 38  ) num = (0.91606*(TMath::Erf((x-2.24886)/1.14822)*TMath::Exp(-0.00000*x)))+0.04169;
      else if (idx == 39  ) num = (0.91572*(TMath::Erf((x-2.24370)/1.17367)*TMath::Exp(-0.00000*x)))+0.04243;
      else if (idx == 40  ) num = (0.91710*(TMath::Erf((x-2.25494)/1.15265)*TMath::Exp(-0.00016*x)))+0.04255;
      else if (idx == 41  ) num = (0.91560*(TMath::Erf((x-2.27996)/1.13411)*TMath::Exp(-0.00007*x)))+0.03964;
      else if (idx == 42  ) num = (0.91635*(TMath::Erf((x-2.24468)/1.17510)*TMath::Exp(-0.00000*x)))+0.04298;
      else if (idx == 43  ) num = (0.91613*(TMath::Erf((x-2.24046)/1.17170)*TMath::Exp(-0.00029*x)))+0.04286;
      else if (idx == 44  ) num = (0.91553*(TMath::Erf((x-2.21496)/1.23001)*TMath::Exp(-0.00000*x)))+0.04429;
      else if (idx == 45  ) num = (0.91385*(TMath::Erf((x-2.20187)/1.22045)*TMath::Exp(-0.00000*x)))+0.04343;
      else if (idx == 46  ) num = (0.91595*(TMath::Erf((x-2.25898)/1.13996)*TMath::Exp(-0.00000*x)))+0.04117;
      else if (idx == 47  ) num = (0.91479*(TMath::Erf((x-2.20434)/1.23089)*TMath::Exp(-0.00000*x)))+0.04429;
      else if (idx == 48  ) num = (0.93296*(TMath::Erf((x-2.20352)/1.19159)*TMath::Exp(-0.00000*x)))+0.02307;
      else if (idx == 49  ) num = (0.91515*(TMath::Erf((x-2.23952)/1.16767)*TMath::Exp(-0.00000*x)))+0.04210;
      else if (idx == 50  ) num = (0.91576*(TMath::Erf((x-2.23112)/1.19488)*TMath::Exp(-0.00000*x)))+0.04320;
      else if (idx == 51  ) num = (0.91760*(TMath::Erf((x-2.30165)/1.07416)*TMath::Exp(-0.00000*x)))+0.03914;
      else if (idx == 52  ) num = (0.91412*(TMath::Erf((x-2.23643)/1.16489)*TMath::Exp(-0.00000*x)))+0.04133;
      else if (idx == 53  ) num = (0.91592*(TMath::Erf((x-2.23210)/1.19744)*TMath::Exp(-0.00000*x)))+0.04334;
      else if (idx == 54  ) num = (0.91455*(TMath::Erf((x-2.22019)/1.19593)*TMath::Exp(-0.00000*x)))+0.04291;
      else if (idx == 55  ) num = (0.91551*(TMath::Erf((x-2.22002)/1.18231)*TMath::Exp(-0.00000*x)))+0.04358;
      else if (idx == 56  ) num = (0.91670*(TMath::Erf((x-2.24853)/1.17987)*TMath::Exp(-0.00000*x)))+0.04311;
      else if (idx == 57  ) num = (0.91605*(TMath::Erf((x-2.25072)/1.16032)*TMath::Exp(-0.00000*x)))+0.04195;
      else if (idx == 58  ) num = (0.91527*(TMath::Erf((x-2.23331)/1.19945)*TMath::Exp(-0.00000*x)))+0.04261;
      else if (idx == 59  ) num = (0.91595*(TMath::Erf((x-2.28381)/1.10885)*TMath::Exp(-0.00017*x)))+0.03923;
      else if (idx == 60  ) num = (0.91509*(TMath::Erf((x-2.23106)/1.18060)*TMath::Exp(-0.00000*x)))+0.04250;
      else if (idx == 61  ) num = (0.91412*(TMath::Erf((x-2.22950)/1.15380)*TMath::Exp(-0.00020*x)))+0.04158;
      else if (idx == 62  ) num = (0.91689*(TMath::Erf((x-2.23513)/1.18675)*TMath::Exp(-0.00000*x)))+0.04401;
      else if (idx == 63  ) num = (0.91306*(TMath::Erf((x-2.19505)/1.23155)*TMath::Exp(-0.00000*x)))+0.04328;
      else if (idx == 64  ) num = (0.91614*(TMath::Erf((x-2.25895)/1.13011)*TMath::Exp(-0.00000*x)))+0.04094;
      else if (idx == 65  ) num = (0.91615*(TMath::Erf((x-2.24869)/1.14070)*TMath::Exp(-0.00000*x)))+0.04183;
      else if (idx == 66  ) num = (0.91686*(TMath::Erf((x-2.21650)/1.22437)*TMath::Exp(-0.00000*x)))+0.04562;
      else if (idx == 67  ) num = (0.91474*(TMath::Erf((x-2.24662)/1.15118)*TMath::Exp(-0.00000*x)))+0.04093;
      else if (idx == 68  ) num = (0.91697*(TMath::Erf((x-2.24931)/1.18305)*TMath::Exp(-0.00000*x)))+0.04335;
      else if (idx == 69  ) num = (0.91395*(TMath::Erf((x-2.19152)/1.23844)*TMath::Exp(-0.00000*x)))+0.04445;
      else if (idx == 70  ) num = (0.91637*(TMath::Erf((x-2.24735)/1.17976)*TMath::Exp(-0.00000*x)))+0.04287;
      else if (idx == 71  ) num = (0.91499*(TMath::Erf((x-2.24305)/1.17595)*TMath::Exp(-0.00000*x)))+0.04186;
      else if (idx == 72  ) num = (0.91715*(TMath::Erf((x-2.23241)/1.20396)*TMath::Exp(-0.00000*x)))+0.04455;
      else if (idx == 73  ) num = (0.91610*(TMath::Erf((x-2.22280)/1.18545)*TMath::Exp(-0.00036*x)))+0.04403;
      else if (idx == 74  ) num = (0.91522*(TMath::Erf((x-2.24363)/1.17444)*TMath::Exp(-0.00000*x)))+0.04200;
      else if (idx == 75  ) num = (0.91583*(TMath::Erf((x-2.24799)/1.16094)*TMath::Exp(-0.00007*x)))+0.04200;
      else if (idx == 76  ) num = (0.91508*(TMath::Erf((x-2.24902)/1.15314)*TMath::Exp(-0.00000*x)))+0.04101;
      else if (idx == 77  ) num = (0.91531*(TMath::Erf((x-2.21356)/1.21929)*TMath::Exp(-0.00000*x)))+0.04423;
      else if (idx == 78  ) num = (0.91829*(TMath::Erf((x-2.26671)/1.13958)*TMath::Exp(-0.00064*x)))+0.04262;
      else if (idx == 79  ) num = (0.91452*(TMath::Erf((x-2.23449)/1.16208)*TMath::Exp(-0.00000*x)))+0.04178;
      else if (idx == 80  ) num = (0.91883*(TMath::Erf((x-2.26213)/1.14875)*TMath::Exp(-0.00000*x)))+0.04345;
      else if (idx == 81  ) num = (0.91560*(TMath::Erf((x-2.25953)/1.13956)*TMath::Exp(-0.00000*x)))+0.04061;
      else if (idx == 82  ) num = (0.91408*(TMath::Erf((x-2.23898)/1.16638)*TMath::Exp(-0.00003*x)))+0.04108;
      else if (idx == 83  ) num = (0.91687*(TMath::Erf((x-2.24909)/1.18231)*TMath::Exp(-0.00000*x)))+0.04327;
      else if (idx == 84  ) num = (0.91638*(TMath::Erf((x-2.26237)/1.13921)*TMath::Exp(-0.00000*x)))+0.04122;
      else if (idx == 85  ) num = (0.91490*(TMath::Erf((x-2.23999)/1.17185)*TMath::Exp(-0.00000*x)))+0.04189;
      else if (idx == 86  ) num = (0.91427*(TMath::Erf((x-2.22135)/1.19577)*TMath::Exp(-0.00000*x)))+0.04235;
      else if (idx == 87  ) num = (0.91765*(TMath::Erf((x-2.20885)/1.22987)*TMath::Exp(-0.00135*x)))+0.04667;
      else if (idx == 88  ) num = (0.91493*(TMath::Erf((x-2.23762)/1.16496)*TMath::Exp(-0.00000*x)))+0.04199;
      else if (idx == 89  ) num = (0.91354*(TMath::Erf((x-2.24043)/1.16732)*TMath::Exp(-0.00000*x)))+0.04061;
      else if (idx == 90  ) num = (0.91669*(TMath::Erf((x-2.24022)/1.17043)*TMath::Exp(-0.00000*x)))+0.04354;
      else if (idx == 91  ) num = (0.91815*(TMath::Erf((x-2.24499)/1.17417)*TMath::Exp(-0.00119*x)))+0.04451;
      else if (idx == 92  ) num = (0.91556*(TMath::Erf((x-2.22880)/1.19251)*TMath::Exp(-0.00000*x)))+0.04313;
      else if (idx == 93  ) num = (0.91590*(TMath::Erf((x-2.25187)/1.15998)*TMath::Exp(-0.00000*x)))+0.04179;
      else if (idx == 94  ) num = (0.91475*(TMath::Erf((x-2.22202)/1.19631)*TMath::Exp(-0.00000*x)))+0.04285;
      else if (idx == 95  ) num = (0.91545*(TMath::Erf((x-2.23933)/1.16607)*TMath::Exp(-0.00000*x)))+0.04236;
      else if (idx == 96  ) num = (0.91638*(TMath::Erf((x-2.27566)/1.12569)*TMath::Exp(-0.00000*x)))+0.04032;
      else if (idx == 97  ) num = (0.91827*(TMath::Erf((x-2.25203)/1.18323)*TMath::Exp(-0.00000*x)))+0.04448;
      else if (idx == 98  ) num = (0.91873*(TMath::Erf((x-2.22044)/1.21240)*TMath::Exp(-0.00159*x)))+0.04667;
      else if (idx == 99  ) num = (0.91738*(TMath::Erf((x-2.26136)/1.14473)*TMath::Exp(-0.00000*x)))+0.04214;
      else if (idx == 100 ) num = (0.91590*(TMath::Erf((x-2.24495)/1.17456)*TMath::Exp(-0.00000*x)))+0.04252;
   }
   else if (fabs(eta)<2.1)
   {
      if (idx==0) num = (0.34844*(TMath::Erf((x-2.42840)/1.83954)*TMath::Exp(-0.04675*x)))+0.66739;
      else if (idx == -1   ) num = (0.41766*(TMath::Erf((x-2.10696)/2.12556)*TMath::Exp(-0.03957*x)))+0.60680;
      else if (idx == -2   ) num = (0.31867*(TMath::Erf((x-2.56424)/1.67845)*TMath::Exp(-0.04915*x)))+0.68905;
      else if (idx == 1   ) num = (0.98606*(TMath::Erf((x-0.34625)/3.52546)*TMath::Exp(-0.02754*x)))+0.12205;
      else if (idx == 2   ) num = (0.92956*(TMath::Erf((x-0.39451)/3.29197)*TMath::Exp(-0.02271*x)))+0.13735;
      else if (idx == 3   ) num = (0.55061*(TMath::Erf((x-1.59642)/2.42317)*TMath::Exp(-0.03070*x)))+0.47850;
      else if (idx == 4   ) num = (0.66956*(TMath::Erf((x-1.16638)/2.75609)*TMath::Exp(-0.02383*x)))+0.36856;
      else if (idx == 5   ) num = (0.19726*(TMath::Erf((x-3.12221)/0.95883)*TMath::Exp(-0.05524*x)))+0.77702;
      else if (idx == 6   ) num = (0.17669*(TMath::Erf((x-3.04011)/0.92640)*TMath::Exp(-0.02881*x)))+0.76831;
      else if (idx == 7   ) num = (0.36331*(TMath::Erf((x-2.42567)/2.00163)*TMath::Exp(-0.05129*x)))+0.67213;
      else if (idx == 8   ) num = (0.20325*(TMath::Erf((x-2.90687)/1.11376)*TMath::Exp(-0.02529*x)))+0.74331;
      else if (idx == 9   ) num = (0.35753*(TMath::Erf((x-2.40106)/2.01298)*TMath::Exp(-0.05180*x)))+0.67100;
      else if (idx == 10  ) num = (0.22486*(TMath::Erf((x-2.98020)/1.19159)*TMath::Exp(-0.05309*x)))+0.75860;
      else if (idx == 11  ) num = (0.19997*(TMath::Erf((x-3.05150)/1.06847)*TMath::Exp(-0.05300*x)))+0.77138;
      else if (idx == 12  ) num = (0.17315*(TMath::Erf((x-3.08381)/0.93228)*TMath::Exp(-0.03695*x)))+0.77731;
      else if (idx == 13  ) num = (0.82763*(TMath::Erf((x-0.75201)/2.84778)*TMath::Exp(-0.01819*x)))+0.19915;
      else if (idx == 14  ) num = (0.56221*(TMath::Erf((x-1.53363)/2.32000)*TMath::Exp(-0.02122*x)))+0.44494;
      else if (idx == 15  ) num = (0.51934*(TMath::Erf((x-1.63827)/2.54485)*TMath::Exp(-0.03045*x)))+0.51000;
      else if (idx == 16  ) num = (0.56042*(TMath::Erf((x-1.77027)/2.44976)*TMath::Exp(-0.04624*x)))+0.51466;
      else if (idx == 17  ) num = (0.34032*(TMath::Erf((x-2.43542)/1.68754)*TMath::Exp(-0.04503*x)))+0.66524;
      else if (idx == 18  ) num = (0.92468*(TMath::Erf((x-0.52108)/3.09878)*TMath::Exp(-0.02195*x)))+0.13162;
      else if (idx == 19  ) num = (0.97589*(TMath::Erf((x-0.42845)/3.29236)*TMath::Exp(-0.02397*x)))+0.10661;
      else if (idx == 20  ) num = (0.52171*(TMath::Erf((x-1.66719)/2.20460)*TMath::Exp(-0.02439*x)))+0.47960;
      else if (idx == 21  ) num = (0.95714*(TMath::Erf((x-0.50765)/2.93934)*TMath::Exp(-0.01645*x)))+0.07665;
      else if (idx == 22  ) num = (0.52812*(TMath::Erf((x-1.71999)/2.29539)*TMath::Exp(-0.03288*x)))+0.50165;
      else if (idx == 23  ) num = (0.21806*(TMath::Erf((x-3.05126)/1.07145)*TMath::Exp(-0.06760*x)))+0.77147;
      else if (idx == 24  ) num = (0.20033*(TMath::Erf((x-3.11978)/0.99375)*TMath::Exp(-0.06188*x)))+0.77702;
      else if (idx == 25  ) num = (0.22251*(TMath::Erf((x-3.05029)/1.13496)*TMath::Exp(-0.05607*x)))+0.76696;
      else if (idx == 26  ) num = (0.35061*(TMath::Erf((x-2.43754)/1.77577)*TMath::Exp(-0.04655*x)))+0.66730;
      else if (idx == 27  ) num = (0.99507*(TMath::Erf((x-0.36481)/3.67458)*TMath::Exp(-0.03004*x)))+0.13955;
      else if (idx == 28  ) num = (0.78858*(TMath::Erf((x-0.88956)/2.84612)*TMath::Exp(-0.02149*x)))+0.25020;
      else if (idx == 29  ) num = (0.45086*(TMath::Erf((x-1.96659)/2.11859)*TMath::Exp(-0.03585*x)))+0.56900;
      else if (idx == 30  ) num = (0.97895*(TMath::Erf((x-0.43611)/3.36280)*TMath::Exp(-0.02556*x)))+0.11776;
      else if (idx == 31  ) num = (0.22547*(TMath::Erf((x-2.91153)/1.28905)*TMath::Exp(-0.03692*x)))+0.74581;
      else if (idx == 32  ) num = (0.91683*(TMath::Erf((x-0.56920)/2.96710)*TMath::Exp(-0.01798*x)))+0.11662;
      else if (idx == 33  ) num = (0.28717*(TMath::Erf((x-2.89872)/1.54706)*TMath::Exp(-0.07190*x)))+0.74097;
      else if (idx == 34  ) num = (0.21731*(TMath::Erf((x-2.90506)/1.30542)*TMath::Exp(-0.03655*x)))+0.74660;
      else if (idx == 35  ) num = (0.22306*(TMath::Erf((x-2.89206)/1.21125)*TMath::Exp(-0.04111*x)))+0.74501;
      else if (idx == 36  ) num = (0.35652*(TMath::Erf((x-2.36220)/1.83709)*TMath::Exp(-0.04330*x)))+0.65707;
      else if (idx == 37  ) num = (0.30338*(TMath::Erf((x-2.62022)/1.60007)*TMath::Exp(-0.04676*x)))+0.69976;
      else if (idx == 38  ) num = (0.80758*(TMath::Erf((x-0.88769)/3.15926)*TMath::Exp(-0.03204*x)))+0.29159;
      else if (idx == 39  ) num = (0.93133*(TMath::Erf((x-0.45163)/3.08495)*TMath::Exp(-0.01711*x)))+0.10271;
      else if (idx == 40  ) num = (0.28332*(TMath::Erf((x-2.57355)/1.59773)*TMath::Exp(-0.03197*x)))+0.69258;
      else if (idx == 41  ) num = (0.23537*(TMath::Erf((x-2.93664)/1.25978)*TMath::Exp(-0.05776*x)))+0.75510;
      else if (idx == 42  ) num = (0.33319*(TMath::Erf((x-2.40753)/1.82646)*TMath::Exp(-0.03711*x)))+0.66759;
      else if (idx == 43  ) num = (0.18912*(TMath::Erf((x-3.08258)/1.07116)*TMath::Exp(-0.03446*x)))+0.76874;
      else if (idx == 44  ) num = (0.93223*(TMath::Erf((x-0.41033)/3.10451)*TMath::Exp(-0.01441*x)))+0.09671;
      else if (idx == 45  ) num = (0.50117*(TMath::Erf((x-1.67507)/2.22412)*TMath::Exp(-0.02517*x)))+0.50153;
      else if (idx == 46  ) num = (0.38227*(TMath::Erf((x-2.46335)/1.90357)*TMath::Exp(-0.06190*x)))+0.67122;
      else if (idx == 47  ) num = (0.22082*(TMath::Erf((x-3.03582)/1.17328)*TMath::Exp(-0.05979*x)))+0.76753;
      else if (idx == 48  ) num = (0.88416*(TMath::Erf((x-0.60615)/3.08067)*TMath::Exp(-0.02096*x)))+0.16958;
      else if (idx == 49  ) num = (0.23618*(TMath::Erf((x-2.85780)/1.29683)*TMath::Exp(-0.04258*x)))+0.74345;
      else if (idx == 50  ) num = (0.29170*(TMath::Erf((x-2.63632)/1.57637)*TMath::Exp(-0.04471*x)))+0.70373;
      else if (idx == 51  ) num = (0.35723*(TMath::Erf((x-2.41712)/1.84671)*TMath::Exp(-0.05244*x)))+0.66849;
      else if (idx == 52  ) num = (0.85447*(TMath::Erf((x-0.74043)/2.88116)*TMath::Exp(-0.02185*x)))+0.18992;
      else if (idx == 53  ) num = (0.31080*(TMath::Erf((x-2.73857)/1.63979)*TMath::Exp(-0.05951*x)))+0.71785;
      else if (idx == 54  ) num = (0.19118*(TMath::Erf((x-3.05885)/0.97559)*TMath::Exp(-0.04578*x)))+0.77077;
      else if (idx == 55  ) num = (0.23273*(TMath::Erf((x-3.03149)/1.15138)*TMath::Exp(-0.07406*x)))+0.76988;
      else if (idx == 56  ) num = (0.32840*(TMath::Erf((x-2.39547)/1.89654)*TMath::Exp(-0.03672*x)))+0.66717;
      else if (idx == 57  ) num = (0.94312*(TMath::Erf((x-0.37747)/3.25109)*TMath::Exp(-0.01917*x)))+0.10701;
      else if (idx == 58  ) num = (0.98626*(TMath::Erf((x-0.39342)/3.35103)*TMath::Exp(-0.02390*x)))+0.10410;
      else if (idx == 59  ) num = (0.28311*(TMath::Erf((x-2.82334)/1.50227)*TMath::Exp(-0.06294*x)))+0.73293;
      else if (idx == 60  ) num = (0.48544*(TMath::Erf((x-1.84079)/2.22102)*TMath::Exp(-0.03223*x)))+0.53660;
      else if (idx == 61  ) num = (0.35954*(TMath::Erf((x-2.45024)/1.89194)*TMath::Exp(-0.05225*x)))+0.66844;
      else if (idx == 62  ) num = (0.20433*(TMath::Erf((x-3.05319)/1.12290)*TMath::Exp(-0.04776*x)))+0.77020;
      else if (idx == 63  ) num = (0.33934*(TMath::Erf((x-2.44293)/1.78934)*TMath::Exp(-0.04174*x)))+0.66498;
      else if (idx == 64  ) num = (0.73553*(TMath::Erf((x-0.94356)/2.86190)*TMath::Exp(-0.02107*x)))+0.29593;
      else if (idx == 65  ) num = (0.19263*(TMath::Erf((x-3.04487)/1.02329)*TMath::Exp(-0.04121*x)))+0.76970;
      else if (idx == 66  ) num = (0.31499*(TMath::Erf((x-2.51611)/1.78366)*TMath::Exp(-0.04105*x)))+0.68159;
      else if (idx == 67  ) num = (0.97314*(TMath::Erf((x-0.29400)/3.40047)*TMath::Exp(-0.02084*x)))+0.09566;
      else if (idx == 68  ) num = (0.54872*(TMath::Erf((x-1.63796)/2.41435)*TMath::Exp(-0.03636*x)))+0.48874;
      else if (idx == 69  ) num = (0.84812*(TMath::Erf((x-0.72677)/2.97132)*TMath::Exp(-0.02049*x)))+0.19740;
      else if (idx == 70  ) num = (0.97740*(TMath::Erf((x-0.24212)/3.33105)*TMath::Exp(-0.01759*x)))+0.06617;
      else if (idx == 71  ) num = (0.44042*(TMath::Erf((x-2.16681)/2.16466)*TMath::Exp(-0.04903*x)))+0.61069;
      else if (idx == 72  ) num = (0.29106*(TMath::Erf((x-2.59438)/1.59128)*TMath::Exp(-0.03653*x)))+0.69700;
      else if (idx == 73  ) num = (0.21354*(TMath::Erf((x-2.98755)/1.21430)*TMath::Exp(-0.04496*x)))+0.75912;
      else if (idx == 74  ) num = (0.19416*(TMath::Erf((x-3.10478)/0.93802)*TMath::Exp(-0.05252*x)))+0.77691;
      else if (idx == 75  ) num = (0.98515*(TMath::Erf((x-0.32343)/3.14781)*TMath::Exp(-0.01562*x)))+0.05130;
      else if (idx == 76  ) num = (0.19066*(TMath::Erf((x-2.98663)/1.04503)*TMath::Exp(-0.03299*x)))+0.76024;
      else if (idx == 77  ) num = (0.60121*(TMath::Erf((x-1.44503)/2.48964)*TMath::Exp(-0.02813*x)))+0.43059;
      else if (idx == 78  ) num = (0.95018*(TMath::Erf((x-0.24431)/3.28054)*TMath::Exp(-0.01554*x)))+0.08055;
      else if (idx == 79  ) num = (0.29152*(TMath::Erf((x-2.71664)/1.58926)*TMath::Exp(-0.04936*x)))+0.71420;
      else if (idx == 80  ) num = (0.83645*(TMath::Erf((x-0.61438)/3.10626)*TMath::Exp(-0.01951*x)))+0.19894;
      else if (idx == 81  ) num = (0.21382*(TMath::Erf((x-2.93329)/1.10851)*TMath::Exp(-0.04674*x)))+0.75659;
      else if (idx == 82  ) num = (0.91802*(TMath::Erf((x-0.38594)/3.38512)*TMath::Exp(-0.01954*x)))+0.14618;
      else if (idx == 83  ) num = (0.64480*(TMath::Erf((x-1.35441)/2.80624)*TMath::Exp(-0.03460*x)))+0.42445;
      else if (idx == 84  ) num = (0.29208*(TMath::Erf((x-2.77496)/1.45841)*TMath::Exp(-0.06376*x)))+0.72629;
      else if (idx == 85  ) num = (0.50329*(TMath::Erf((x-1.71188)/2.23493)*TMath::Exp(-0.02810*x)))+0.50209;
      else if (idx == 86  ) num = (0.36565*(TMath::Erf((x-2.45407)/1.81635)*TMath::Exp(-0.05470*x)))+0.66867;
      else if (idx == 87  ) num = (0.60620*(TMath::Erf((x-1.34506)/2.63652)*TMath::Exp(-0.02745*x)))+0.42722;
      else if (idx == 88  ) num = (0.39672*(TMath::Erf((x-2.10028)/2.00563)*TMath::Exp(-0.03004*x)))+0.60122;
      else if (idx == 89  ) num = (0.32743*(TMath::Erf((x-2.40928)/1.80779)*TMath::Exp(-0.03853*x)))+0.66401;
      else if (idx == 90  ) num = (0.34622*(TMath::Erf((x-2.42513)/1.78578)*TMath::Exp(-0.04874*x)))+0.66571;
      else if (idx == 91  ) num = (0.18033*(TMath::Erf((x-3.06503)/0.91306)*TMath::Exp(-0.03577*x)))+0.77136;
      else if (idx == 92  ) num = (0.58441*(TMath::Erf((x-1.39171)/2.45221)*TMath::Exp(-0.02100*x)))+0.41650;
      else if (idx == 93  ) num = (0.49394*(TMath::Erf((x-1.76932)/2.24226)*TMath::Exp(-0.02970*x)))+0.52318;
      else if (idx == 94  ) num = (0.95013*(TMath::Erf((x-0.50591)/3.08801)*TMath::Exp(-0.02079*x)))+0.10818;
      else if (idx == 95  ) num = (0.95441*(TMath::Erf((x-0.38257)/3.18652)*TMath::Exp(-0.01665*x)))+0.08774;
      else if (idx == 96  ) num = (0.46372*(TMath::Erf((x-1.91637)/2.29474)*TMath::Exp(-0.03974*x)))+0.56942;
      else if (idx == 97  ) num = (0.21306*(TMath::Erf((x-3.04485)/1.08285)*TMath::Exp(-0.05445*x)))+0.76924;
      else if (idx == 98  ) num = (0.15324*(TMath::Erf((x-3.01940)/0.74429)*TMath::Exp(-0.01371*x)))+0.76954;
      else if (idx == 99  ) num = (0.21822*(TMath::Erf((x-3.05144)/1.01727)*TMath::Exp(-0.06135*x)))+0.76794;
      else if (idx == 100 ) num = (0.21049*(TMath::Erf((x-2.91620)/1.24286)*TMath::Exp(-0.02410*x)))+0.74188;
   }
   else
   {
      if (idx==0) num = (0.46966*(TMath::Erf((x-0.49035)/2.40643)*TMath::Exp(-0.00715*x)))+0.45423;
      else if (idx == -1   ) num = (0.46984*(TMath::Erf((x-0.48605)/2.39891)*TMath::Exp(-0.00716*x)))+0.45439;
      else if (idx == -2   ) num = (0.46958*(TMath::Erf((x-0.49357)/2.41355)*TMath::Exp(-0.00710*x)))+0.45415;
      else if (idx == 1   ) num = (0.46600*(TMath::Erf((x-0.67911)/2.09606)*TMath::Exp(-0.00517*x)))+0.44736;
      else if (idx == 2   ) num = (0.47692*(TMath::Erf((x-0.11926)/3.11360)*TMath::Exp(-0.00875*x)))+0.46590;
      else if (idx == 3   ) num = (0.47115*(TMath::Erf((x-0.18845)/2.81839)*TMath::Exp(-0.00399*x)))+0.44898;
      else if (idx == 4   ) num = (0.45346*(TMath::Erf((x-0.65391)/1.94557)*TMath::Exp(-0.00006*x)))+0.43924;
      else if (idx == 5   ) num = (0.46830*(TMath::Erf((x-0.45314)/2.39493)*TMath::Exp(-0.00645*x)))+0.45352;
      else if (idx == 6   ) num = (0.48089*(TMath::Erf((x-0.46824)/2.19316)*TMath::Exp(-0.00001*x)))+0.41390;
      else if (idx == 7   ) num = (0.46104*(TMath::Erf((x-0.39623)/2.40150)*TMath::Exp(-0.00000*x)))+0.44890;
      else if (idx == 8   ) num = (0.47551*(TMath::Erf((x-0.28373)/2.80729)*TMath::Exp(-0.00847*x)))+0.46073;
      else if (idx == 9   ) num = (0.63735*(TMath::Erf((x-0.02429)/2.67383)*TMath::Exp(-0.01024*x)))+0.30407;
      else if (idx == 10  ) num = (0.46123*(TMath::Erf((x-0.80634)/1.84176)*TMath::Exp(-0.00184*x)))+0.44265;
      else if (idx == 11  ) num = (0.45932*(TMath::Erf((x-0.54875)/2.19503)*TMath::Exp(-0.00000*x)))+0.44414;
      else if (idx == 12  ) num = (0.60191*(TMath::Erf((x-0.01151)/2.76483)*TMath::Exp(-0.01259*x)))+0.34727;
      else if (idx == 13  ) num = (0.47341*(TMath::Erf((x-0.32608)/2.62191)*TMath::Exp(-0.00668*x)))+0.45387;
      else if (idx == 14  ) num = (0.46790*(TMath::Erf((x-0.57957)/2.28587)*TMath::Exp(-0.00826*x)))+0.45173;
      else if (idx == 15  ) num = (0.52837*(TMath::Erf((x-0.00926)/3.39781)*TMath::Exp(-0.02144*x)))+0.45092;
      else if (idx == 16  ) num = (0.46895*(TMath::Erf((x-0.42796)/2.46395)*TMath::Exp(-0.00477*x)))+0.45491;
      else if (idx == 17  ) num = (0.46201*(TMath::Erf((x-0.55432)/2.33964)*TMath::Exp(-0.00001*x)))+0.44784;
      else if (idx == 18  ) num = (0.49652*(TMath::Erf((x-0.19666)/3.16119)*TMath::Exp(-0.02232*x)))+0.47798;
      else if (idx == 19  ) num = (0.46152*(TMath::Erf((x-0.40297)/2.45333)*TMath::Exp(-0.00000*x)))+0.44650;
      else if (idx == 20  ) num = (0.48696*(TMath::Erf((x-0.20566)/3.10981)*TMath::Exp(-0.01637*x)))+0.47157;
      else if (idx == 21  ) num = (0.49084*(TMath::Erf((x-0.17707)/3.10697)*TMath::Exp(-0.02098*x)))+0.47431;
      else if (idx == 22  ) num = (0.64005*(TMath::Erf((x-0.00424)/6.34049)*TMath::Exp(-0.06996*x)))+0.58062;
      else if (idx == 23  ) num = (0.47072*(TMath::Erf((x-0.60809)/2.20996)*TMath::Exp(-0.00861*x)))+0.44910;
      else if (idx == 24  ) num = (0.46043*(TMath::Erf((x-0.77055)/1.90434)*TMath::Exp(-0.00190*x)))+0.44335;
      else if (idx == 25  ) num = (0.45431*(TMath::Erf((x-0.84448)/1.76636)*TMath::Exp(-0.00007*x)))+0.43736;
      else if (idx == 26  ) num = (0.50013*(TMath::Erf((x-0.61451)/2.53455)*TMath::Exp(-0.02409*x)))+0.47504;
      else if (idx == 27  ) num = (0.47302*(TMath::Erf((x-0.52479)/2.42150)*TMath::Exp(-0.00736*x)))+0.45643;
      else if (idx == 28  ) num = (0.46546*(TMath::Erf((x-0.50500)/2.30282)*TMath::Exp(-0.00467*x)))+0.44461;
      else if (idx == 29  ) num = (0.47186*(TMath::Erf((x-0.51065)/2.37613)*TMath::Exp(-0.00882*x)))+0.45530;
      else if (idx == 30  ) num = (0.47089*(TMath::Erf((x-0.35575)/2.58180)*TMath::Exp(-0.00833*x)))+0.45550;
      else if (idx == 31  ) num = (0.48483*(TMath::Erf((x-0.45881)/2.48702)*TMath::Exp(-0.01341*x)))+0.45923;
      else if (idx == 32  ) num = (0.50439*(TMath::Erf((x-0.26054)/3.20335)*TMath::Exp(-0.02606*x)))+0.48375;
      else if (idx == 33  ) num = (0.47287*(TMath::Erf((x-0.00347)/3.67745)*TMath::Exp(-0.02178*x)))+0.50667;
      else if (idx == 34  ) num = (0.47189*(TMath::Erf((x-0.05281)/3.05598)*TMath::Exp(-0.01064*x)))+0.46401;
      else if (idx == 35  ) num = (0.48844*(TMath::Erf((x-0.17639)/3.15514)*TMath::Exp(-0.01550*x)))+0.47146;
      else if (idx == 36  ) num = (0.48026*(TMath::Erf((x-0.32777)/2.66725)*TMath::Exp(-0.01229*x)))+0.45485;
      else if (idx == 37  ) num = (0.46781*(TMath::Erf((x-0.00017)/3.51845)*TMath::Exp(-0.02025*x)))+0.50037;
      else if (idx == 38  ) num = (0.46386*(TMath::Erf((x-0.81109)/1.90202)*TMath::Exp(-0.00446*x)))+0.44498;
      else if (idx == 39  ) num = (0.48310*(TMath::Erf((x-0.04417)/3.28801)*TMath::Exp(-0.01630*x)))+0.47334;
      else if (idx == 40  ) num = (0.46622*(TMath::Erf((x-0.52946)/2.27725)*TMath::Exp(-0.00193*x)))+0.44539;
      else if (idx == 41  ) num = (0.48692*(TMath::Erf((x-0.01393)/3.14005)*TMath::Exp(-0.00909*x)))+0.45273;
      else if (idx == 42  ) num = (0.45904*(TMath::Erf((x-0.68323)/2.07581)*TMath::Exp(-0.00000*x)))+0.44290;
      else if (idx == 43  ) num = (0.45722*(TMath::Erf((x-0.52994)/2.21623)*TMath::Exp(-0.00000*x)))+0.44433;
      else if (idx == 44  ) num = (0.46025*(TMath::Erf((x-0.49986)/2.24493)*TMath::Exp(-0.00000*x)))+0.44574;
      else if (idx == 45  ) num = (0.47353*(TMath::Erf((x-0.50368)/2.35622)*TMath::Exp(-0.00958*x)))+0.45658;
      else if (idx == 46  ) num = (0.51162*(TMath::Erf((x-0.29791)/3.29705)*TMath::Exp(-0.02664*x)))+0.48929;
      else if (idx == 47  ) num = (0.46059*(TMath::Erf((x-0.45416)/2.32819)*TMath::Exp(-0.00177*x)))+0.44501;
      else if (idx == 48  ) num = (0.46328*(TMath::Erf((x-0.61101)/2.15933)*TMath::Exp(-0.00000*x)))+0.44737;
      else if (idx == 49  ) num = (0.47179*(TMath::Erf((x-0.13700)/2.83970)*TMath::Exp(-0.01033*x)))+0.45983;
      else if (idx == 50  ) num = (0.46290*(TMath::Erf((x-0.00564)/3.90898)*TMath::Exp(-0.03018*x)))+0.53395;
      else if (idx == 51  ) num = (0.47020*(TMath::Erf((x-0.35082)/2.60623)*TMath::Exp(-0.00634*x)))+0.45336;
      else if (idx == 52  ) num = (0.48449*(TMath::Erf((x-0.08916)/3.15882)*TMath::Exp(-0.01672*x)))+0.47330;
      else if (idx == 53  ) num = (0.45951*(TMath::Erf((x-0.74966)/1.89145)*TMath::Exp(-0.00024*x)))+0.44015;
      else if (idx == 54  ) num = (0.46645*(TMath::Erf((x-0.83215)/1.86917)*TMath::Exp(-0.00564*x)))+0.44633;
      else if (idx == 55  ) num = (0.46416*(TMath::Erf((x-0.42221)/2.41527)*TMath::Exp(-0.00597*x)))+0.45067;
      else if (idx == 56  ) num = (0.45854*(TMath::Erf((x-0.45275)/2.23926)*TMath::Exp(-0.00087*x)))+0.43791;
      else if (idx == 57  ) num = (0.50448*(TMath::Erf((x-0.38715)/2.28237)*TMath::Exp(-0.00017*x)))+0.39662;
      else if (idx == 58  ) num = (0.46983*(TMath::Erf((x-0.69735)/2.07846)*TMath::Exp(-0.00744*x)))+0.45068;
      else if (idx == 59  ) num = (0.45840*(TMath::Erf((x-0.65511)/2.01268)*TMath::Exp(-0.00261*x)))+0.44209;
      else if (idx == 60  ) num = (0.49470*(TMath::Erf((x-0.16296)/3.18671)*TMath::Exp(-0.02064*x)))+0.47676;
      else if (idx == 61  ) num = (0.45982*(TMath::Erf((x-0.85657)/1.78932)*TMath::Exp(-0.00049*x)))+0.43992;
      else if (idx == 62  ) num = (0.50030*(TMath::Erf((x-0.11141)/3.42504)*TMath::Exp(-0.02278*x)))+0.48362;
      else if (idx == 63  ) num = (0.47046*(TMath::Erf((x-0.06463)/3.15032)*TMath::Exp(-0.00769*x)))+0.46261;
      else if (idx == 64  ) num = (0.47209*(TMath::Erf((x-0.52422)/2.37018)*TMath::Exp(-0.00902*x)))+0.45529;
      else if (idx == 65  ) num = (0.46755*(TMath::Erf((x-0.28964)/2.69359)*TMath::Exp(-0.00414*x)))+0.45495;
      else if (idx == 66  ) num = (0.45490*(TMath::Erf((x-0.72072)/1.93714)*TMath::Exp(-0.00003*x)))+0.44058;
      else if (idx == 67  ) num = (0.50886*(TMath::Erf((x-0.26254)/2.37796)*TMath::Exp(-0.00023*x)))+0.38832;
      else if (idx == 68  ) num = (0.46238*(TMath::Erf((x-0.32866)/2.50795)*TMath::Exp(-0.00000*x)))+0.44648;
      else if (idx == 69  ) num = (0.46478*(TMath::Erf((x-0.85936)/1.84542)*TMath::Exp(-0.00507*x)))+0.44428;
      else if (idx == 70  ) num = (0.46208*(TMath::Erf((x-0.61722)/2.17406)*TMath::Exp(-0.00000*x)))+0.44608;
      else if (idx == 71  ) num = (0.48371*(TMath::Erf((x-0.29986)/2.85392)*TMath::Exp(-0.01443*x)))+0.46607;
      else if (idx == 72  ) num = (0.47121*(TMath::Erf((x-0.37454)/2.59811)*TMath::Exp(-0.00615*x)))+0.45302;
      else if (idx == 73  ) num = (0.47199*(TMath::Erf((x-0.38083)/2.66019)*TMath::Exp(-0.00722*x)))+0.45770;
      else if (idx == 74  ) num = (0.91844*(TMath::Erf((x-0.00227)/9.87151)*TMath::Exp(-0.09257*x)))+0.60562;
      else if (idx == 75  ) num = (0.47380*(TMath::Erf((x-0.54123)/2.38738)*TMath::Exp(-0.00891*x)))+0.45655;
      else if (idx == 76  ) num = (0.45930*(TMath::Erf((x-0.22934)/2.72916)*TMath::Exp(-0.00000*x)))+0.45137;
      else if (idx == 77  ) num = (0.46061*(TMath::Erf((x-0.82184)/1.87764)*TMath::Exp(-0.00013*x)))+0.44277;
      else if (idx == 78  ) num = (0.45559*(TMath::Erf((x-0.37117)/2.43872)*TMath::Exp(-0.00000*x)))+0.44788;
      else if (idx == 79  ) num = (0.50413*(TMath::Erf((x-0.00011)/4.68558)*TMath::Exp(-0.03772*x)))+0.55203;
      else if (idx == 80  ) num = (0.45922*(TMath::Erf((x-0.79219)/1.84311)*TMath::Exp(-0.00072*x)))+0.44096;
      else if (idx == 81  ) num = (0.45807*(TMath::Erf((x-0.24367)/2.57270)*TMath::Exp(-0.00019*x)))+0.45012;
      else if (idx == 82  ) num = (0.45663*(TMath::Erf((x-0.66609)/2.00535)*TMath::Exp(-0.00001*x)))+0.44188;
      else if (idx == 83  ) num = (0.47970*(TMath::Erf((x-0.38146)/2.59412)*TMath::Exp(-0.00778*x)))+0.45202;
      else if (idx == 84  ) num = (0.47141*(TMath::Erf((x-0.87187)/1.85984)*TMath::Exp(-0.00889*x)))+0.45001;
      else if (idx == 85  ) num = (0.46090*(TMath::Erf((x-0.80677)/1.93369)*TMath::Exp(-0.00259*x)))+0.44323;
      else if (idx == 86  ) num = (0.45641*(TMath::Erf((x-0.80082)/1.75842)*TMath::Exp(-0.00013*x)))+0.43797;
      else if (idx == 87  ) num = (0.45564*(TMath::Erf((x-1.02489)/1.52097)*TMath::Exp(-0.00000*x)))+0.43741;
      else if (idx == 88  ) num = (0.47429*(TMath::Erf((x-0.43210)/2.52407)*TMath::Exp(-0.00802*x)))+0.45352;
      else if (idx == 89  ) num = (0.46163*(TMath::Erf((x-0.63181)/2.03682)*TMath::Exp(-0.00243*x)))+0.44447;
      else if (idx == 90  ) num = (0.47158*(TMath::Erf((x-0.15015)/2.92293)*TMath::Exp(-0.00796*x)))+0.46129;
      else if (idx == 91  ) num = (0.46155*(TMath::Erf((x-0.53068)/2.24367)*TMath::Exp(-0.00000*x)))+0.44687;
      else if (idx == 92  ) num = (0.50277*(TMath::Erf((x-0.12519)/3.44181)*TMath::Exp(-0.02250*x)))+0.48563;
      else if (idx == 93  ) num = (0.46909*(TMath::Erf((x-0.56454)/2.26183)*TMath::Exp(-0.00379*x)))+0.44491;
      else if (idx == 94  ) num = (0.49849*(TMath::Erf((x-0.00730)/3.66396)*TMath::Exp(-0.02361*x)))+0.49015;
      else if (idx == 95  ) num = (0.45773*(TMath::Erf((x-0.67700)/2.01777)*TMath::Exp(-0.00013*x)))+0.44259;
      else if (idx == 96  ) num = (0.46384*(TMath::Erf((x-0.42174)/2.42757)*TMath::Exp(-0.00590*x)))+0.45047;
      else if (idx == 97  ) num = (0.47086*(TMath::Erf((x-0.61911)/2.15541)*TMath::Exp(-0.00351*x)))+0.44176;
      else if (idx == 98  ) num = (0.46237*(TMath::Erf((x-0.20028)/2.79741)*TMath::Exp(-0.00000*x)))+0.45492;
      else if (idx == 99  ) num = (0.46102*(TMath::Erf((x-0.70485)/2.01826)*TMath::Exp(-0.00002*x)))+0.44133;
      else if (idx == 100 ) num = (0.52455*(TMath::Erf((x-0.00249)/3.35289)*TMath::Exp(-0.02111*x)))+0.45048;
   }

   // return
   return num/den;
}

///////////////////////////////////////////////////
//                 M U I D    P P                //
///////////////////////////////////////////////////

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!  ONLY FOR SYSTEMATICS! DO NOT APPLY FOR THE NOMINAL CORRECTION!!! !!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double tnp_weight_muid_pp(double x, double eta, int idx) {
   // binned
   if (idx == -10) {
      if (fabs(eta) < 1.2) {
         // 0 < |eta| < 1.2
         if (x<4) return 0.998994;
         else if (x<4.5) return 1.00322;
         else if (x<5) return 1.00128;
         else if (x<5.5) return 0.995761;
         else if (x<6) return 0.996457;
         else if (x<7) return 0.999095;
         else if (x<8) return 0.997853;
         else if (x<10.5) return 0.996943;
         else if (x<14) return 0.998211;
         else if (x<18) return 1.00255;
         else return 1.00194;
      } else if (fabs(eta)<1.8) {
         // 1.2 < |eta| < 1.8
         if (x<3) return 1.00148;
         else if (x<3.5) return 0.998129;
         else if (x<4) return 1.00036;
         else if (x<4.5) return 0.999244;
         else if (x<5) return 0.99598;
         else if (x<6) return 0.999206;
         else if (x<7) return 0.998106;
         else if (x<9) return 1.00012;
         else if (x<14) return 1.00015;
         else if (x<18) return 0.985758;
         else return 0.984038;
      } else if (fabs(eta)<2.1) {
         // 1.8 < |eta| < 2.1
         if (x<2.5) return 1.00626;
         else if (x<3) return 1.00149;
         else if (x<3.5) return 1.00119;
         else if (x<4) return 0.999079;
         else if (x<4.5) return 1.00136;
         else if (x<5) return 0.998673;
         else if (x<6) return 1.00312;
         else if (x<7) return 0.998298;
         else if (x<9) return 0.997686;
         else if (x<13) return 1.00114;
         else return 0.982895;
      } else {
         // 2.1 < |eta| < 2.4
         if (x<2.5) return 1.00188;
         else if (x<3) return 0.99958;
         else if (x<3.5) return 1.00014;
         else if (x<4) return 0.998119;
         else if (x<4.5) return 1.00171;
         else if (x<5) return 0.998421;
         else if (x<6.5) return 1.00079;
         else if (x<8.5) return 0.996052;
         else if (x<12) return 0.999585;
         else return 0.99234;
      }
   }

   // denominator (from MC)
   double den=1;
   if (fabs(eta)<1.2) den = (0.98383*TMath::Erf((x+0.47673)/2.94587))+0.00653;
   else if (fabs(eta)<1.8) den = (0.02768*TMath::Erf((x+8.66696)/28.28792))+0.97612;
   else if (fabs(eta)<2.1) den = (0.74373*TMath::Erf((x-0.30826)/1.07049))+0.25156;
   else den = (0.82009*TMath::Erf((x-1.22789)/0.48924))+0.17143;

   // numerator (from data)
   double num=1;
   if (fabs(eta)<1.2)
   {
      if (idx==0) num = (0.97976*TMath::Erf((x+0.74721)/3.10382))+0.00985;
      else if (idx == -1   ) num = (0.97990*TMath::Erf((x+0.75177)/3.10073))+0.01000;
      else if (idx == -2   ) num = (0.97962*TMath::Erf((x+0.74218)/3.10736))+0.00972;
      else if (idx == 1   ) num = (0.98034*TMath::Erf((x+2.80697)/4.45288))+0.01191;
      else if (idx == 2   ) num = (0.98037*TMath::Erf((x-0.86725)/1.99017))+0.00654;
      else if (idx == 3   ) num = (0.97941*TMath::Erf((x+0.72288)/3.11992))+0.00952;
      else if (idx == 4   ) num = (0.97946*TMath::Erf((x+0.76939)/3.09098))+0.00957;
      else if (idx == 5   ) num = (0.97895*TMath::Erf((x+0.78240)/3.08157))+0.00911;
      else if (idx == 6   ) num = (0.97986*TMath::Erf((x+0.71392)/3.12731))+0.00994;
      else if (idx == 7   ) num = (0.97921*TMath::Erf((x+0.75601)/3.09874))+0.00935;
      else if (idx == 8   ) num = (0.98012*TMath::Erf((x+0.74874)/3.09937))+0.01023;
      else if (idx == 9   ) num = (0.98034*TMath::Erf((x+0.69031)/3.13848))+0.01042;
      else if (idx == 10  ) num = (0.97896*TMath::Erf((x+1.91862)/3.86757))+0.01000;
      else if (idx == 11  ) num = (0.98008*TMath::Erf((x+2.11350)/4.04508))+0.01169;
      else if (idx == 12  ) num = (0.98008*TMath::Erf((x+0.76590)/3.09158))+0.01020;
      else if (idx == 13  ) num = (0.97983*TMath::Erf((x+0.75360)/3.09995))+0.00993;
      else if (idx == 14  ) num = (0.97997*TMath::Erf((x+1.31859)/3.49726))+0.01014;
      else if (idx == 15  ) num = (0.97986*TMath::Erf((x+0.75641)/3.09748))+0.00996;
      else if (idx == 16  ) num = (0.97973*TMath::Erf((x-0.80809)/2.02720))+0.00711;
      else if (idx == 17  ) num = (0.97953*TMath::Erf((x+0.76185)/3.09692))+0.00963;
      else if (idx == 18  ) num = (0.97988*TMath::Erf((x+0.71155)/3.12622))+0.00996;
      else if (idx == 19  ) num = (0.97964*TMath::Erf((x+0.78575)/3.08041))+0.00974;
      else if (idx == 20  ) num = (0.97972*TMath::Erf((x+0.76339)/3.09446))+0.00981;
      else if (idx == 21  ) num = (0.98036*TMath::Erf((x+0.71760)/3.12110))+0.01047;
      else if (idx == 22  ) num = (0.97869*TMath::Erf((x+0.78971)/3.07753))+0.00890;
      else if (idx == 23  ) num = (0.97986*TMath::Erf((x+0.75234)/3.10015))+0.00996;
      else if (idx == 24  ) num = (0.97989*TMath::Erf((x+2.01694)/4.03072))+0.01154;
      else if (idx == 25  ) num = (0.97931*TMath::Erf((x-0.45786)/2.26456))+0.00805;
      else if (idx == 26  ) num = (0.97977*TMath::Erf((x+0.72151)/3.11724))+0.00986;
      else if (idx == 27  ) num = (0.98017*TMath::Erf((x+0.74084)/3.10719))+0.01027;
      else if (idx == 28  ) num = (0.98025*TMath::Erf((x+1.36375)/3.52191))+0.01052;
      else if (idx == 29  ) num = (0.97905*TMath::Erf((x+0.76549)/3.09429))+0.00920;
      else if (idx == 30  ) num = (0.98034*TMath::Erf((x+1.99489)/3.99493))+0.01156;
      else if (idx == 31  ) num = (0.97913*TMath::Erf((x+3.06124)/4.68785))+0.01271;
      else if (idx == 32  ) num = (0.98027*TMath::Erf((x+1.86330)/3.90489))+0.01118;
      else if (idx == 33  ) num = (0.97967*TMath::Erf((x+0.77052)/3.08783))+0.00977;
      else if (idx == 34  ) num = (0.98031*TMath::Erf((x+0.74244)/3.10719))+0.01043;
      else if (idx == 35  ) num = (0.97936*TMath::Erf((x+0.78025)/3.08291))+0.00947;
      else if (idx == 36  ) num = (0.97990*TMath::Erf((x+0.77674)/3.08624))+0.01000;
      else if (idx == 37  ) num = (0.97994*TMath::Erf((x+0.75412)/3.09892))+0.01004;
      else if (idx == 38  ) num = (0.97988*TMath::Erf((x+0.75525)/3.09815))+0.00998;
      else if (idx == 39  ) num = (0.97811*TMath::Erf((x+4.17632)/5.41568))+0.01444;
      else if (idx == 40  ) num = (0.97906*TMath::Erf((x+0.12271)/2.66119))+0.00879;
      else if (idx == 41  ) num = (0.97916*TMath::Erf((x+0.78958)/3.07611))+0.00929;
      else if (idx == 42  ) num = (0.98034*TMath::Erf((x+0.71172)/3.12549))+0.01044;
      else if (idx == 43  ) num = (0.97977*TMath::Erf((x-0.36353)/2.33829))+0.00897;
      else if (idx == 44  ) num = (0.97964*TMath::Erf((x+0.71058)/3.12869))+0.00973;
      else if (idx == 45  ) num = (0.97968*TMath::Erf((x+0.02509)/2.59241))+0.00930;
      else if (idx == 46  ) num = (0.97920*TMath::Erf((x+0.78531)/3.07976))+0.00933;
      else if (idx == 47  ) num = (0.98028*TMath::Erf((x+1.90520)/3.90684))+0.01145;
      else if (idx == 48  ) num = (0.98069*TMath::Erf((x+0.71365)/3.12361))+0.01085;
      else if (idx == 49  ) num = (0.97786*TMath::Erf((x+4.31008)/5.48301))+0.01492;
      else if (idx == 50  ) num = (0.97979*TMath::Erf((x+0.73236)/3.11137))+0.00988;
      else if (idx == 51  ) num = (0.97935*TMath::Erf((x+0.77039)/3.09245))+0.00947;
      else if (idx == 52  ) num = (0.97989*TMath::Erf((x+0.71901)/3.11909))+0.00997;
      else if (idx == 53  ) num = (0.97971*TMath::Erf((x+0.75937)/3.09576))+0.00981;
      else if (idx == 54  ) num = (0.97951*TMath::Erf((x+0.75256)/3.09918))+0.00961;
      else if (idx == 55  ) num = (0.97970*TMath::Erf((x+0.73916)/3.10877))+0.00979;
      else if (idx == 56  ) num = (0.97976*TMath::Erf((x+2.37276)/4.15895))+0.01121;
      else if (idx == 57  ) num = (0.97909*TMath::Erf((x+2.19973)/4.06462))+0.01044;
      else if (idx == 58  ) num = (0.97988*TMath::Erf((x+0.75481)/3.09853))+0.00998;
      else if (idx == 59  ) num = (0.98004*TMath::Erf((x+0.71848)/3.11927))+0.01012;
      else if (idx == 60  ) num = (0.97975*TMath::Erf((x+0.77204)/3.08983))+0.00985;
      else if (idx == 61  ) num = (0.98004*TMath::Erf((x+0.70843)/3.13003))+0.01012;
      else if (idx == 62  ) num = (0.98029*TMath::Erf((x+0.70596)/3.13035))+0.01039;
      else if (idx == 63  ) num = (0.98008*TMath::Erf((x+2.33435)/4.23779))+0.01202;
      else if (idx == 64  ) num = (0.97962*TMath::Erf((x+0.74263)/3.10726))+0.00972;
      else if (idx == 65  ) num = (0.97961*TMath::Erf((x+2.28882)/4.15209))+0.01125;
      else if (idx == 66  ) num = (0.97965*TMath::Erf((x+2.27300)/4.13958))+0.01122;
      else if (idx == 67  ) num = (0.97990*TMath::Erf((x+0.75184)/3.10034))+0.01000;
      else if (idx == 68  ) num = (0.97994*TMath::Erf((x+0.71192)/3.12476))+0.01002;
      else if (idx == 69  ) num = (0.97933*TMath::Erf((x+0.75599)/3.09987))+0.00945;
      else if (idx == 70  ) num = (0.97956*TMath::Erf((x+0.73436)/3.11271))+0.00966;
      else if (idx == 71  ) num = (0.97969*TMath::Erf((x+2.01570)/3.98558))+0.01076;
      else if (idx == 72  ) num = (0.97907*TMath::Erf((x+0.76096)/3.09840))+0.00922;
      else if (idx == 73  ) num = (0.97966*TMath::Erf((x+2.91346)/4.46914))+0.01097;
      else if (idx == 74  ) num = (0.97925*TMath::Erf((x+0.77710)/3.08710))+0.00937;
      else if (idx == 75  ) num = (0.97960*TMath::Erf((x+0.75692)/3.09891))+0.00970;
      else if (idx == 76  ) num = (0.98061*TMath::Erf((x+1.59162)/3.69649))+0.01127;
      else if (idx == 77  ) num = (0.98012*TMath::Erf((x+2.53578)/4.30562))+0.01131;
      else if (idx == 78  ) num = (0.98000*TMath::Erf((x+0.74152)/3.11007))+0.01010;
      else if (idx == 79  ) num = (0.97939*TMath::Erf((x+2.87409)/4.46236))+0.01075;
      else if (idx == 80  ) num = (0.97975*TMath::Erf((x+0.75080)/3.10202))+0.00984;
      else if (idx == 81  ) num = (0.97947*TMath::Erf((x+0.77254)/3.09114))+0.00958;
      else if (idx == 82  ) num = (0.97957*TMath::Erf((x+0.73845)/3.10918))+0.00967;
      else if (idx == 83  ) num = (0.97972*TMath::Erf((x+0.74545)/3.10497))+0.00982;
      else if (idx == 84  ) num = (0.97959*TMath::Erf((x+0.74064)/3.10878))+0.00969;
      else if (idx == 85  ) num = (0.97973*TMath::Erf((x+1.81417)/3.82785))+0.01037;
      else if (idx == 86  ) num = (0.98055*TMath::Erf((x+0.71181)/3.12616))+0.01067;
      else if (idx == 87  ) num = (0.98013*TMath::Erf((x+0.75509)/3.10030))+0.01025;
      else if (idx == 88  ) num = (0.98087*TMath::Erf((x+2.52163)/4.29365))+0.01251;
      else if (idx == 89  ) num = (0.97949*TMath::Erf((x+0.76139)/3.09589))+0.00960;
      else if (idx == 90  ) num = (0.97947*TMath::Erf((x+0.78473)/3.08177))+0.00958;
      else if (idx == 91  ) num = (0.97995*TMath::Erf((x+0.75814)/3.09683))+0.01005;
      else if (idx == 92  ) num = (0.97874*TMath::Erf((x+0.77514)/3.08771))+0.00894;
      else if (idx == 93  ) num = (0.97895*TMath::Erf((x+0.78670)/3.08197))+0.00911;
      else if (idx == 94  ) num = (0.97997*TMath::Erf((x+1.37165)/3.50437))+0.01017;
      else if (idx == 95  ) num = (0.98031*TMath::Erf((x+0.73502)/3.11004))+0.01043;
      else if (idx == 96  ) num = (0.97925*TMath::Erf((x+0.81140)/3.06554))+0.00937;
      else if (idx == 97  ) num = (0.98006*TMath::Erf((x+0.74433)/3.10588))+0.01016;
      else if (idx == 98  ) num = (0.98017*TMath::Erf((x+0.74246)/3.10629))+0.01028;
      else if (idx == 99  ) num = (0.97846*TMath::Erf((x+0.83701)/3.04796))+0.00869;
      else if (idx == 100 ) num = (0.98013*TMath::Erf((x+0.72817)/3.11401))+0.01022;
   }
   else if (fabs(eta)<1.8)
   {
      if (idx==0) num = (0.65913*TMath::Erf((x+7.42715)/3.84453))+0.32953;
      else if (idx == -1   ) num = (0.65951*TMath::Erf((x+7.73063)/3.66662))+0.32992;
      else if (idx == -2   ) num = (0.65873*TMath::Erf((x+7.75539)/3.58881))+0.32913;
      else if (idx == 1   ) num = (0.65830*TMath::Erf((x+5.65712)/3.54639))+0.32871;
      else if (idx == 2   ) num = (0.65888*TMath::Erf((x+7.75487)/3.56581))+0.32929;
      else if (idx == 3   ) num = (0.65918*TMath::Erf((x+7.95689)/3.55897))+0.32959;
      else if (idx == 4   ) num = (0.65936*TMath::Erf((x+7.25007)/5.05766))+0.32977;
      else if (idx == 5   ) num = (0.65852*TMath::Erf((x+8.06675)/3.51869))+0.32893;
      else if (idx == 6   ) num = (0.65908*TMath::Erf((x+7.97942)/3.56818))+0.32949;
      else if (idx == 7   ) num = (0.65916*TMath::Erf((x+7.37090)/4.91615))+0.32957;
      else if (idx == 8   ) num = (0.65883*TMath::Erf((x+5.52751)/3.62770))+0.32924;
      else if (idx == 9   ) num = (0.65953*TMath::Erf((x+7.67868)/2.55408))+0.32993;
      else if (idx == 10  ) num = (0.65839*TMath::Erf((x+8.61149)/3.11054))+0.32880;
      else if (idx == 11  ) num = (0.65820*TMath::Erf((x+8.81841)/2.99068))+0.32861;
      else if (idx == 12  ) num = (0.65895*TMath::Erf((x+7.79233)/3.26253))+0.32936;
      else if (idx == 13  ) num = (0.65922*TMath::Erf((x+7.82067)/3.64682))+0.32963;
      else if (idx == 14  ) num = (0.65868*TMath::Erf((x+5.36040)/3.72951))+0.32909;
      else if (idx == 15  ) num = (0.00003*TMath::Erf((x+9.34815)/9.44508))+0.98903;
      else if (idx == 16  ) num = (0.65800*TMath::Erf((x+8.54428)/3.06931))+0.32841;
      else if (idx == 17  ) num = (0.65848*TMath::Erf((x+8.12618)/3.26218))+0.32889;
      else if (idx == 18  ) num = (0.65843*TMath::Erf((x+6.55875)/3.52166))+0.32884;
      else if (idx == 19  ) num = (0.00004*TMath::Erf((x+9.45640)/15.46847))+0.99057;
      else if (idx == 20  ) num = (0.65839*TMath::Erf((x+7.47845)/3.72144))+0.32880;
      else if (idx == 21  ) num = (0.65906*TMath::Erf((x+7.92821)/3.58792))+0.32947;
      else if (idx == 22  ) num = (0.65923*TMath::Erf((x+7.82306)/3.64996))+0.32964;
      else if (idx == 23  ) num = (0.65831*TMath::Erf((x+8.76448)/3.04681))+0.32872;
      else if (idx == 24  ) num = (0.66058*TMath::Erf((x+9.49218)/6.66294))+0.33105;
      else if (idx == 25  ) num = (0.66006*TMath::Erf((x+9.55025)/6.56137))+0.33048;
      else if (idx == 26  ) num = (0.66004*TMath::Erf((x+7.29759)/4.57431))+0.33045;
      else if (idx == 27  ) num = (0.65918*TMath::Erf((x+6.71529)/5.12029))+0.32961;
      else if (idx == 28  ) num = (0.65823*TMath::Erf((x+6.82733)/3.49450))+0.32864;
      else if (idx == 29  ) num = (0.65901*TMath::Erf((x+7.60584)/3.77243))+0.32942;
      else if (idx == 30  ) num = (0.00016*TMath::Erf((x+9.45101)/20.58607))+0.99126;
      else if (idx == 31  ) num = (0.65920*TMath::Erf((x+7.77921)/3.66701))+0.32961;
      else if (idx == 32  ) num = (0.00209*TMath::Erf((x+9.45910)/19.94435))+0.98856;
      else if (idx == 33  ) num = (0.65887*TMath::Erf((x+7.76943)/3.67871))+0.32928;
      else if (idx == 34  ) num = (0.65894*TMath::Erf((x+7.22870)/3.62746))+0.32935;
      else if (idx == 35  ) num = (0.66039*TMath::Erf((x+9.92039)/6.98058))+0.33087;
      else if (idx == 36  ) num = (0.65867*TMath::Erf((x+8.26705)/3.31191))+0.32908;
      else if (idx == 37  ) num = (0.65937*TMath::Erf((x+7.72104)/3.71122))+0.32978;
      else if (idx == 38  ) num = (0.65964*TMath::Erf((x+7.72885)/3.69714))+0.33005;
      else if (idx == 39  ) num = (0.65927*TMath::Erf((x+7.67557)/3.70011))+0.32967;
      else if (idx == 40  ) num = (0.65951*TMath::Erf((x+7.70722)/3.71684))+0.32992;
      else if (idx == 41  ) num = (0.65845*TMath::Erf((x+4.39991)/3.64432))+0.32887;
      else if (idx == 42  ) num = (0.65961*TMath::Erf((x+6.79096)/4.44302))+0.33001;
      else if (idx == 43  ) num = (0.65919*TMath::Erf((x+7.82513)/3.63965))+0.32960;
      else if (idx == 44  ) num = (0.66011*TMath::Erf((x+8.70932)/6.21011))+0.33056;
      else if (idx == 45  ) num = (0.65898*TMath::Erf((x+7.83436)/3.61512))+0.32939;
      else if (idx == 46  ) num = (0.65932*TMath::Erf((x+7.71203)/3.71036))+0.32972;
      else if (idx == 47  ) num = (0.00019*TMath::Erf((x+9.58115)/8.79196))+0.98884;
      else if (idx == 48  ) num = (0.65886*TMath::Erf((x+7.77063)/3.67122))+0.32926;
      else if (idx == 49  ) num = (0.65931*TMath::Erf((x+7.73612)/3.69291))+0.32972;
      else if (idx == 50  ) num = (0.65958*TMath::Erf((x+7.79937)/5.06362))+0.32999;
      else if (idx == 51  ) num = (0.65966*TMath::Erf((x+5.36707)/4.01020))+0.33007;
      else if (idx == 52  ) num = (0.65893*TMath::Erf((x+7.77995)/3.65931))+0.32934;
      else if (idx == 53  ) num = (0.65981*TMath::Erf((x+8.00261)/5.47730))+0.33022;
      else if (idx == 54  ) num = (0.65906*TMath::Erf((x+7.99536)/3.52879))+0.32947;
      else if (idx == 55  ) num = (0.66004*TMath::Erf((x+7.73399)/3.69697))+0.33045;
      else if (idx == 56  ) num = (0.66012*TMath::Erf((x+5.80936)/4.53863))+0.33053;
      else if (idx == 57  ) num = (0.65970*TMath::Erf((x+8.72419)/5.54925))+0.33009;
      else if (idx == 58  ) num = (0.65955*TMath::Erf((x+8.73281)/6.08808))+0.32996;
      else if (idx == 59  ) num = (0.65857*TMath::Erf((x+5.37007)/3.68700))+0.32898;
      else if (idx == 60  ) num = (0.65922*TMath::Erf((x+7.78856)/3.66736))+0.32963;
      else if (idx == 61  ) num = (0.65870*TMath::Erf((x+8.32684)/3.37184))+0.32911;
      else if (idx == 62  ) num = (0.02169*TMath::Erf((x+9.39993)/8.81451))+0.96950;
      else if (idx == 63  ) num = (0.65861*TMath::Erf((x+8.41608)/2.83607))+0.32902;
      else if (idx == 64  ) num = (0.65935*TMath::Erf((x+7.73712)/3.69013))+0.32976;
      else if (idx == 65  ) num = (0.65964*TMath::Erf((x+7.73015)/3.70227))+0.33005;
      else if (idx == 66  ) num = (0.65986*TMath::Erf((x+7.71948)/4.41047))+0.33027;
      else if (idx == 67  ) num = (0.65872*TMath::Erf((x+8.51285)/3.01665))+0.32913;
      else if (idx == 68  ) num = (0.66047*TMath::Erf((x+8.41592)/6.16326))+0.33094;
      else if (idx == 69  ) num = (0.65922*TMath::Erf((x+7.82428)/3.62519))+0.32963;
      else if (idx == 70  ) num = (0.66047*TMath::Erf((x+9.99355)/7.10794))+0.33099;
      else if (idx == 71  ) num = (0.65885*TMath::Erf((x+5.47069)/3.64802))+0.32926;
      else if (idx == 72  ) num = (0.00282*TMath::Erf((x-5.38298)/1.24327))+0.99004;
      else if (idx == 73  ) num = (0.65944*TMath::Erf((x+7.74039)/3.69609))+0.32985;
      else if (idx == 74  ) num = (0.65915*TMath::Erf((x+5.52818)/3.83173))+0.32956;
      else if (idx == 75  ) num = (0.65915*TMath::Erf((x+7.91105)/3.60900))+0.32956;
      else if (idx == 76  ) num = (0.65894*TMath::Erf((x+7.76803)/3.64816))+0.32934;
      else if (idx == 77  ) num = (0.65947*TMath::Erf((x+7.54041)/2.14953))+0.32988;
      else if (idx == 78  ) num = (0.65886*TMath::Erf((x+8.02060)/3.25923))+0.32927;
      else if (idx == 79  ) num = (0.65937*TMath::Erf((x+7.74931)/3.69068))+0.32977;
      else if (idx == 80  ) num = (0.61448*TMath::Erf((x+5.72793)/4.61157))+0.37840;
      else if (idx == 81  ) num = (0.65921*TMath::Erf((x+7.73428)/3.68511))+0.32962;
      else if (idx == 82  ) num = (0.65869*TMath::Erf((x+7.74361)/3.66920))+0.32910;
      else if (idx == 83  ) num = (0.66062*TMath::Erf((x+9.27937)/6.73552))+0.33108;
      else if (idx == 84  ) num = (0.65939*TMath::Erf((x+7.71195)/3.70838))+0.32979;
      else if (idx == 85  ) num = (0.65894*TMath::Erf((x+7.65272)/2.69727))+0.32935;
      else if (idx == 86  ) num = (0.65895*TMath::Erf((x+7.79596)/3.60989))+0.32936;
      else if (idx == 87  ) num = (0.65778*TMath::Erf((x+8.10325)/2.85733))+0.32819;
      else if (idx == 88  ) num = (0.65920*TMath::Erf((x+7.88649)/3.60984))+0.32961;
      else if (idx == 89  ) num = (0.65821*TMath::Erf((x+7.07368)/3.23135))+0.32862;
      else if (idx == 90  ) num = (0.00952*TMath::Erf((x+1.92845)/6.04223))+0.98200;
      else if (idx == 91  ) num = (0.65783*TMath::Erf((x+6.57855)/3.32438))+0.32824;
      else if (idx == 92  ) num = (0.65936*TMath::Erf((x+9.93265)/3.72713))+0.32977;
      else if (idx == 93  ) num = (0.65944*TMath::Erf((x+6.55330)/4.03829))+0.32977;
      else if (idx == 94  ) num = (0.65922*TMath::Erf((x+7.77667)/3.64848))+0.32963;
      else if (idx == 95  ) num = (0.65951*TMath::Erf((x+7.73059)/3.69493))+0.32992;
      else if (idx == 96  ) num = (0.65931*TMath::Erf((x+5.14790)/3.79332))+0.32972;
      else if (idx == 97  ) num = (0.66136*TMath::Erf((x+9.60062)/7.11838))+0.33221;
      else if (idx == 98  ) num = (0.65776*TMath::Erf((x+8.44103)/2.98451))+0.32817;
      else if (idx == 99  ) num = (0.00004*TMath::Erf((x+9.37311)/17.50347))+0.98965;
      else if (idx == 100 ) num = (0.65958*TMath::Erf((x+9.05568)/5.76041))+0.32999;
   }
   else if (fabs(eta)<2.1)
   {
      if (idx==0) num = (0.99392*TMath::Erf((x-0.14813)/0.77299))+0.00001;
      else if (idx == -1   ) num = (0.99383*TMath::Erf((x-0.06281)/0.73936))+0.00090;
      else if (idx == -2   ) num = (0.99282*TMath::Erf((x-0.01734)/0.71980))+0.00000;
      else if (idx == 1   ) num = (0.99003*TMath::Erf((x-0.01641)/0.71921))+0.00000;
      else if (idx == 2   ) num = (0.99214*TMath::Erf((x-0.00807)/0.71582))+0.00000;
      else if (idx == 3   ) num = (0.99314*TMath::Erf((x-0.01454)/0.71798))+0.00000;
      else if (idx == 4   ) num = (0.99392*TMath::Erf((x+0.15381)/0.66026))+0.00077;
      else if (idx == 5   ) num = (0.99371*TMath::Erf((x+0.00480)/0.71015))+0.00000;
      else if (idx == 6   ) num = (0.99395*TMath::Erf((x-0.00429)/0.71574))+0.00081;
      else if (idx == 7   ) num = (0.99395*TMath::Erf((x+0.10132)/0.67246))+0.00051;
      else if (idx == 8   ) num = (0.99386*TMath::Erf((x-0.69994)/0.75965))+0.00058;
      else if (idx == 9   ) num = (0.99276*TMath::Erf((x-0.00843)/0.71620))+0.00000;
      else if (idx == 10  ) num = (0.99390*TMath::Erf((x+0.03019)/0.69672))+0.00028;
      else if (idx == 11  ) num = (0.99394*TMath::Erf((x+0.12484)/0.66375))+0.00011;
      else if (idx == 12  ) num = (0.99393*TMath::Erf((x+0.19755)/0.63575))+0.00083;
      else if (idx == 13  ) num = (0.99370*TMath::Erf((x-0.36375)/0.76862))+0.00000;
      else if (idx == 14  ) num = (0.99402*TMath::Erf((x-0.67955)/0.78997))+0.00175;
      else if (idx == 15  ) num = (0.99395*TMath::Erf((x+0.19719)/0.63756))+0.00029;
      else if (idx == 16  ) num = (0.99083*TMath::Erf((x-0.19115)/0.66635))+0.00000;
      else if (idx == 17  ) num = (0.99323*TMath::Erf((x-0.45554)/0.76272))+0.00000;
      else if (idx == 18  ) num = (0.99099*TMath::Erf((x+0.13363)/0.62059))+0.00000;
      else if (idx == 19  ) num = (0.69301*TMath::Erf((x+4.85094)/2.40209))+0.29792;
      else if (idx == 20  ) num = (0.99371*TMath::Erf((x-0.62167)/0.78530))+0.00000;
      else if (idx == 21  ) num = (0.99398*TMath::Erf((x-0.56243)/0.76409))+0.00000;
      else if (idx == 22  ) num = (0.99280*TMath::Erf((x-0.43228)/0.74336))+0.00000;
      else if (idx == 23  ) num = (0.99316*TMath::Erf((x-0.00869)/0.71589))+0.00000;
      else if (idx == 24  ) num = (0.99301*TMath::Erf((x-0.00872)/0.71585))+0.00000;
      else if (idx == 25  ) num = (0.99396*TMath::Erf((x-0.48122)/0.79575))+0.00056;
      else if (idx == 26  ) num = (0.99390*TMath::Erf((x-0.00388)/0.71461))+0.00044;
      else if (idx == 27  ) num = (0.99086*TMath::Erf((x-0.01110)/0.71665))+0.00000;
      else if (idx == 28  ) num = (0.78863*TMath::Erf((x+7.34670)/2.21655))+0.20383;
      else if (idx == 29  ) num = (0.99358*TMath::Erf((x-0.01252)/0.71735))+0.00000;
      else if (idx == 30  ) num = (0.71215*TMath::Erf((x+0.73969)/0.43620))+0.27887;
      else if (idx == 31  ) num = (0.99397*TMath::Erf((x-0.07640)/0.74441))+0.00072;
      else if (idx == 32  ) num = (0.99142*TMath::Erf((x-0.51018)/0.75512))+0.00000;
      else if (idx == 33  ) num = (0.99252*TMath::Erf((x-0.01322)/0.71731))+0.00000;
      else if (idx == 34  ) num = (0.76191*TMath::Erf((x+5.49038)/0.01681))+0.23212;
      else if (idx == 35  ) num = (0.99364*TMath::Erf((x-0.02966)/0.72385))+0.00000;
      else if (idx == 36  ) num = (0.99402*TMath::Erf((x-0.06566)/0.73988))+0.00066;
      else if (idx == 37  ) num = (0.99396*TMath::Erf((x-0.06064)/0.73641))+0.00001;
      else if (idx == 38  ) num = (0.99234*TMath::Erf((x-0.00872)/0.71613))+0.00000;
      else if (idx == 39  ) num = (0.70397*TMath::Erf((x+5.11462)/1.44750))+0.28846;
      else if (idx == 40  ) num = (0.99396*TMath::Erf((x-0.08404)/0.74394))+0.00032;
      else if (idx == 41  ) num = (0.99397*TMath::Erf((x+0.09112)/0.67567))+0.00023;
      else if (idx == 42  ) num = (0.99396*TMath::Erf((x+0.04633)/0.69537))+0.00161;
      else if (idx == 43  ) num = (0.99393*TMath::Erf((x+0.19912)/0.63639))+0.00082;
      else if (idx == 44  ) num = (0.99399*TMath::Erf((x-0.08128)/0.74426))+0.00100;
      else if (idx == 45  ) num = (0.99403*TMath::Erf((x-0.06157)/0.73821))+0.00067;
      else if (idx == 46  ) num = (0.99293*TMath::Erf((x-0.02771)/0.72327))+0.00000;
      else if (idx == 47  ) num = (0.99398*TMath::Erf((x+0.02654)/0.70198))+0.00003;
      else if (idx == 48  ) num = (0.99298*TMath::Erf((x-0.01088)/0.71678))+0.00000;
      else if (idx == 49  ) num = (0.99396*TMath::Erf((x-0.45066)/0.77140))+0.00142;
      else if (idx == 50  ) num = (0.99369*TMath::Erf((x-0.37351)/0.76883))+0.00000;
      else if (idx == 51  ) num = (0.99303*TMath::Erf((x-0.01076)/0.71730))+0.00000;
      else if (idx == 52  ) num = (0.99191*TMath::Erf((x-0.04425)/0.72952))+0.00000;
      else if (idx == 53  ) num = (0.99393*TMath::Erf((x-0.51304)/0.78600))+0.00027;
      else if (idx == 54  ) num = (0.99383*TMath::Erf((x+0.07826)/0.68150))+0.00001;
      else if (idx == 55  ) num = (0.99358*TMath::Erf((x-0.00072)/0.71275))+0.00000;
      else if (idx == 56  ) num = (0.99394*TMath::Erf((x+0.19600)/0.63666))+0.00030;
      else if (idx == 57  ) num = (0.99260*TMath::Erf((x-0.00860)/0.71583))+0.00000;
      else if (idx == 58  ) num = (0.70668*TMath::Erf((x+4.63980)/2.11225))+0.28541;
      else if (idx == 59  ) num = (0.99317*TMath::Erf((x-0.01400)/0.71779))+0.00000;
      else if (idx == 60  ) num = (0.99388*TMath::Erf((x-0.07753)/0.74401))+0.00000;
      else if (idx == 61  ) num = (0.99264*TMath::Erf((x-0.01654)/0.71871))+0.00000;
      else if (idx == 62  ) num = (0.99400*TMath::Erf((x-0.00623)/0.71490))+0.00098;
      else if (idx == 63  ) num = (0.99341*TMath::Erf((x-0.56366)/0.76660))+0.00000;
      else if (idx == 64  ) num = (0.99323*TMath::Erf((x-0.01278)/0.71801))+0.00000;
      else if (idx == 65  ) num = (0.99393*TMath::Erf((x+0.20077)/0.63620))+0.00087;
      else if (idx == 66  ) num = (0.99279*TMath::Erf((x-0.01178)/0.71702))+0.00000;
      else if (idx == 67  ) num = (0.99229*TMath::Erf((x-0.00973)/0.71664))+0.00000;
      else if (idx == 68  ) num = (0.99058*TMath::Erf((x-0.33125)/0.72509))+0.00000;
      else if (idx == 69  ) num = (0.99275*TMath::Erf((x-0.28144)/0.71377))+0.00000;
      else if (idx == 70  ) num = (0.99190*TMath::Erf((x-0.01149)/0.71739))+0.00000;
      else if (idx == 71  ) num = (0.99295*TMath::Erf((x-0.01491)/0.71790))+0.00000;
      else if (idx == 72  ) num = (0.99397*TMath::Erf((x-0.46573)/0.81423))+0.00127;
      else if (idx == 73  ) num = (0.99395*TMath::Erf((x+0.10113)/0.67215))+0.00058;
      else if (idx == 74  ) num = (0.99394*TMath::Erf((x+0.18992)/0.63532))+0.00032;
      else if (idx == 75  ) num = (0.99397*TMath::Erf((x+0.11443)/0.66646))+0.00012;
      else if (idx == 76  ) num = (0.99392*TMath::Erf((x+0.08396)/0.67716))+0.00023;
      else if (idx == 77  ) num = (0.99329*TMath::Erf((x-0.08224)/0.74593))+0.00000;
      else if (idx == 78  ) num = (0.99321*TMath::Erf((x-0.00790)/0.71561))+0.00000;
      else if (idx == 79  ) num = (0.99393*TMath::Erf((x+0.45528)/1.20694))+0.00003;
      else if (idx == 80  ) num = (0.99205*TMath::Erf((x-0.01206)/0.71720))+0.00000;
      else if (idx == 81  ) num = (0.99393*TMath::Erf((x+0.20424)/0.63683))+0.00086;
      else if (idx == 82  ) num = (0.99392*TMath::Erf((x+0.16977)/0.64610))+0.00095;
      else if (idx == 83  ) num = (0.99221*TMath::Erf((x-0.00925)/0.71609))+0.00000;
      else if (idx == 84  ) num = (0.99297*TMath::Erf((x-0.49743)/0.75112))+0.00000;
      else if (idx == 85  ) num = (0.99165*TMath::Erf((x-0.42585)/0.71699))+0.00000;
      else if (idx == 86  ) num = (0.99323*TMath::Erf((x-0.42550)/0.75024))+0.00000;
      else if (idx == 87  ) num = (0.99375*TMath::Erf((x+0.02742)/0.70175))+0.00000;
      else if (idx == 88  ) num = (0.99405*TMath::Erf((x-0.46837)/0.78630))+0.00072;
      else if (idx == 89  ) num = (0.99282*TMath::Erf((x-0.02588)/0.72339))+0.00000;
      else if (idx == 90  ) num = (0.99263*TMath::Erf((x-0.01165)/0.71692))+0.00000;
      else if (idx == 91  ) num = (0.99402*TMath::Erf((x-0.04233)/0.73004))+0.00102;
      else if (idx == 92  ) num = (0.99247*TMath::Erf((x-0.01259)/0.67362))+0.00000;
      else if (idx == 93  ) num = (0.99394*TMath::Erf((x-0.00465)/0.71699))+0.00063;
      else if (idx == 94  ) num = (0.76450*TMath::Erf((x+6.85423)/3.28641))+0.22779;
      else if (idx == 95  ) num = (0.99391*TMath::Erf((x+0.06211)/0.68637))+0.00021;
      else if (idx == 96  ) num = (0.99304*TMath::Erf((x-0.00966)/0.71680))+0.00000;
      else if (idx == 97  ) num = (0.99394*TMath::Erf((x-0.00376)/0.71546))+0.00063;
      else if (idx == 98  ) num = (0.99386*TMath::Erf((x-0.00429)/0.71520))+0.00091;
      else if (idx == 99  ) num = (0.99322*TMath::Erf((x-0.61439)/0.76578))+0.00000;
      else if (idx == 100 ) num = (0.99400*TMath::Erf((x-0.58399)/0.78876))+0.00077;
   }
   else
   {
      if (idx==0) num = (0.99005*TMath::Erf((x-0.34628)/0.89806))+0.00000;
      else if (idx == -1   ) num = (0.99023*TMath::Erf((x-0.36059)/0.90574))+0.00122;
      else if (idx == -2   ) num = (0.98868*TMath::Erf((x-0.30155)/0.87534))+0.00000;
      else if (idx == 1   ) num = (0.98991*TMath::Erf((x-0.40243)/0.91905))+0.00000;
      else if (idx == 2   ) num = (0.99016*TMath::Erf((x-0.37861)/0.91389))+0.00092;
      else if (idx == 3   ) num = (0.98966*TMath::Erf((x-0.36754)/0.90978))+0.00000;
      else if (idx == 4   ) num = (0.98875*TMath::Erf((x-0.15875)/0.79508))+0.00000;
      else if (idx == 5   ) num = (0.98990*TMath::Erf((x-0.25715)/0.85026))+0.00051;
      else if (idx == 6   ) num = (0.98819*TMath::Erf((x-0.34147)/0.89823))+0.00000;
      else if (idx == 7   ) num = (0.98998*TMath::Erf((x-0.30322)/0.87565))+0.00265;
      else if (idx == 8   ) num = (0.99013*TMath::Erf((x-0.37021)/0.91050))+0.00055;
      else if (idx == 9   ) num = (0.98836*TMath::Erf((x-0.23087)/0.83610))+0.00000;
      else if (idx == 10  ) num = (0.98997*TMath::Erf((x-0.31582)/0.88252))+0.00052;
      else if (idx == 11  ) num = (0.99018*TMath::Erf((x-0.35910)/0.90468))+0.00067;
      else if (idx == 12  ) num = (0.98907*TMath::Erf((x-0.07519)/0.74473))+0.00000;
      else if (idx == 13  ) num = (0.99006*TMath::Erf((x-0.32196)/0.88502))+0.00089;
      else if (idx == 14  ) num = (0.98756*TMath::Erf((x-0.33858)/0.89723))+0.00000;
      else if (idx == 15  ) num = (0.98714*TMath::Erf((x-0.31840)/0.88711))+0.00000;
      else if (idx == 16  ) num = (0.98990*TMath::Erf((x-0.28029)/0.86329))+0.00223;
      else if (idx == 17  ) num = (0.99042*TMath::Erf((x-0.43028)/0.94519))+0.00137;
      else if (idx == 18  ) num = (0.98775*TMath::Erf((x-0.04865)/0.73572))+0.00000;
      else if (idx == 19  ) num = (0.99003*TMath::Erf((x-0.35982)/0.90551))+0.00100;
      else if (idx == 20  ) num = (0.98798*TMath::Erf((x-0.31120)/0.88467))+0.00002;
      else if (idx == 21  ) num = (0.98687*TMath::Erf((x-0.14734)/0.79117))+0.00000;
      else if (idx == 22  ) num = (0.98825*TMath::Erf((x-0.26109)/0.85279))+0.00000;
      else if (idx == 23  ) num = (0.98780*TMath::Erf((x-0.30613)/0.87955))+0.00000;
      else if (idx == 24  ) num = (0.99008*TMath::Erf((x-0.35728)/0.90364))+0.00035;
      else if (idx == 25  ) num = (0.98842*TMath::Erf((x-0.36682)/0.90919))+0.00000;
      else if (idx == 26  ) num = (0.98840*TMath::Erf((x-0.29196)/0.86951))+0.00000;
      else if (idx == 27  ) num = (0.99019*TMath::Erf((x-0.37071)/0.91151))+0.00104;
      else if (idx == 28  ) num = (0.98810*TMath::Erf((x-0.36737)/0.90567))+0.00000;
      else if (idx == 29  ) num = (0.98954*TMath::Erf((x-0.30877)/0.87814))+0.00000;
      else if (idx == 30  ) num = (0.98885*TMath::Erf((x-0.20481)/0.82198))+0.00000;
      else if (idx == 31  ) num = (0.98982*TMath::Erf((x-0.20604)/0.82252))+0.00030;
      else if (idx == 32  ) num = (0.98751*TMath::Erf((x-0.27035)/0.85858))+0.00000;
      else if (idx == 33  ) num = (0.98801*TMath::Erf((x-0.08930)/0.75352))+0.00000;
      else if (idx == 34  ) num = (0.98894*TMath::Erf((x-0.08298)/0.74829))+0.00000;
      else if (idx == 35  ) num = (0.98887*TMath::Erf((x-0.35147)/0.90320))+0.00000;
      else if (idx == 36  ) num = (0.98727*TMath::Erf((x-0.25415)/0.84950))+0.00000;
      else if (idx == 37  ) num = (0.98745*TMath::Erf((x-0.11079)/0.76611))+0.00000;
      else if (idx == 38  ) num = (0.98979*TMath::Erf((x-0.40399)/0.92019))+0.00000;
      else if (idx == 39  ) num = (0.98765*TMath::Erf((x-0.17171)/0.80282))+0.00000;
      else if (idx == 40  ) num = (0.99007*TMath::Erf((x-0.35502)/0.90323))+0.00124;
      else if (idx == 41  ) num = (0.98985*TMath::Erf((x-0.27721)/0.86121))+0.00002;
      else if (idx == 42  ) num = (0.64946*TMath::Erf((x-0.54467)/0.91622))+0.34098;
      else if (idx == 43  ) num = (0.98987*TMath::Erf((x-0.26704)/0.85600))+0.00044;
      else if (idx == 44  ) num = (0.98999*TMath::Erf((x-0.32088)/0.88473))+0.00186;
      else if (idx == 45  ) num = (0.98980*TMath::Erf((x-0.19189)/0.81458))+0.00020;
      else if (idx == 46  ) num = (0.98951*TMath::Erf((x-0.33046)/0.89084))+0.00000;
      else if (idx == 47  ) num = (0.98942*TMath::Erf((x-0.29028)/0.86895))+0.00000;
      else if (idx == 48  ) num = (0.99002*TMath::Erf((x-0.36099)/0.90564))+0.00331;
      else if (idx == 49  ) num = (0.98842*TMath::Erf((x-0.16422)/0.79819))+0.00000;
      else if (idx == 50  ) num = (0.98629*TMath::Erf((x-0.10145)/0.76389))+0.00000;
      else if (idx == 51  ) num = (0.98972*TMath::Erf((x-0.24760)/0.84479))+0.00001;
      else if (idx == 52  ) num = (0.98844*TMath::Erf((x-0.16410)/0.79766))+0.00000;
      else if (idx == 53  ) num = (0.99004*TMath::Erf((x-0.34517)/0.89750))+0.00077;
      else if (idx == 54  ) num = (0.98994*TMath::Erf((x-0.36053)/0.90490))+0.00000;
      else if (idx == 55  ) num = (0.98845*TMath::Erf((x-0.22460)/0.83252))+0.00000;
      else if (idx == 56  ) num = (0.98812*TMath::Erf((x-0.24517)/0.84399))+0.00000;
      else if (idx == 57  ) num = (0.98985*TMath::Erf((x-0.37598)/0.91174))+0.00000;
      else if (idx == 58  ) num = (0.98986*TMath::Erf((x-0.29370)/0.87024))+0.00003;
      else if (idx == 59  ) num = (0.98836*TMath::Erf((x-0.25562)/0.84984))+0.00000;
      else if (idx == 60  ) num = (0.98815*TMath::Erf((x-0.21971)/0.82957))+0.00000;
      else if (idx == 61  ) num = (0.99020*TMath::Erf((x-0.39317)/0.91817))+0.00050;
      else if (idx == 62  ) num = (0.98842*TMath::Erf((x-0.21847)/0.82904))+0.00000;
      else if (idx == 63  ) num = (0.98899*TMath::Erf((x-0.34358)/0.89912))+0.00000;
      else if (idx == 64  ) num = (0.98992*TMath::Erf((x-0.36324)/0.90558))+0.00000;
      else if (idx == 65  ) num = (0.99002*TMath::Erf((x-0.28116)/0.86321))+0.00099;
      else if (idx == 66  ) num = (0.98903*TMath::Erf((x-0.31134)/0.88125))+0.00000;
      else if (idx == 67  ) num = (0.98817*TMath::Erf((x-0.27425)/0.86092))+0.00000;
      else if (idx == 68  ) num = (0.99002*TMath::Erf((x-0.31759)/0.88362))+0.00141;
      else if (idx == 69  ) num = (0.98893*TMath::Erf((x-0.40290)/0.91558))+0.00000;
      else if (idx == 70  ) num = (0.99029*TMath::Erf((x-0.45137)/0.92178))+0.00213;
      else if (idx == 71  ) num = (0.98887*TMath::Erf((x-0.32863)/0.88928))+0.00000;
      else if (idx == 72  ) num = (0.98999*TMath::Erf((x-0.34876)/0.89979))+0.00000;
      else if (idx == 73  ) num = (0.99001*TMath::Erf((x-0.35885)/0.90491))+0.00000;
      else if (idx == 74  ) num = (0.98930*TMath::Erf((x+0.03258)/0.69978))+0.00003;
      else if (idx == 75  ) num = (0.99003*TMath::Erf((x-0.35238)/0.90118))+0.00000;
      else if (idx == 76  ) num = (0.99026*TMath::Erf((x-0.40797)/0.92386))+0.00080;
      else if (idx == 77  ) num = (0.99034*TMath::Erf((x-0.44258)/0.92515))+0.00141;
      else if (idx == 78  ) num = (0.98993*TMath::Erf((x-0.29084)/0.86901))+0.00018;
      else if (idx == 79  ) num = (0.98948*TMath::Erf((x-0.28391)/0.86584))+0.00000;
      else if (idx == 80  ) num = (0.98997*TMath::Erf((x-0.30888)/0.87873))+0.00016;
      else if (idx == 81  ) num = (0.98987*TMath::Erf((x-0.25678)/0.84999))+0.00153;
      else if (idx == 82  ) num = (0.98994*TMath::Erf((x-0.29598)/0.87183))+0.00015;
      else if (idx == 83  ) num = (0.99018*TMath::Erf((x-0.36271)/0.90649))+0.00066;
      else if (idx == 84  ) num = (0.98926*TMath::Erf((x-0.33985)/0.89515))+0.00000;
      else if (idx == 85  ) num = (0.98921*TMath::Erf((x-0.42443)/0.91967))+0.00000;
      else if (idx == 86  ) num = (0.98985*TMath::Erf((x-0.24274)/0.84210))+0.00009;
      else if (idx == 87  ) num = (0.98991*TMath::Erf((x-0.38700)/0.91432))+0.00000;
      else if (idx == 88  ) num = (0.98962*TMath::Erf((x-0.33612)/0.89378))+0.00000;
      else if (idx == 89  ) num = (0.98981*TMath::Erf((x-0.19458)/0.81551))+0.00007;
      else if (idx == 90  ) num = (0.98984*TMath::Erf((x-0.24411)/0.84265))+0.00030;
      else if (idx == 91  ) num = (0.99004*TMath::Erf((x-0.35316)/0.90167))+0.00255;
      else if (idx == 92  ) num = (0.98971*TMath::Erf((x-0.31044)/0.87978))+0.00000;
      else if (idx == 93  ) num = (0.99032*TMath::Erf((x-0.40482)/0.92386))+0.00001;
      else if (idx == 94  ) num = (0.98777*TMath::Erf((x-0.02262)/0.72190))+0.00000;
      else if (idx == 95  ) num = (0.99007*TMath::Erf((x-0.33958)/0.89492))+0.00040;
      else if (idx == 96  ) num = (0.98799*TMath::Erf((x-0.23252)/0.83726))+0.00000;
      else if (idx == 97  ) num = (0.99019*TMath::Erf((x-0.38996)/0.91732))+0.00066;
      else if (idx == 98  ) num = (0.98957*TMath::Erf((x-0.33942)/0.89736))+0.00292;
      else if (idx == 99  ) num = (0.99019*TMath::Erf((x-0.38857)/0.91848))+0.00066;
      else if (idx == 100 ) num = (0.98634*TMath::Erf((x-0.23997)/0.84378))+0.00000;
   }

   // return
   return num/den;
}

///////////////////////////////////////////////////
//                   S T A    P P                //
///////////////////////////////////////////////////

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!  ONLY FOR SYSTEMATICS! DO NOT APPLY FOR THE NOMINAL CORRECTION!!! !!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double tnp_weight_sta_pp(double x, double eta, int idx) {
   // binned
   if (idx == -10) {
      if (fabs(eta) < 1.2) {
         // 0 < |eta| < 1.2
         if (x<4) return 1.03808;
         else if (x<4.5) return 1.00984;
         else if (x<5) return 1.01512;
         else if (x<5.5) return 1.00469;
         else if (x<6) return 1.00776;
         else if (x<7) return 0.997825;
         else if (x<8) return 0.996708;
         else if (x<9) return 1.00285;
         else if (x<14) return 0.998507;
         else return 1.00917;
      } else if (fabs(eta)<1.8) {
         // 1.2 < |eta| < 1.8
         if (x<3) return 1.04355;
         else if (x<3.5) return 1.02278;
         else if (x<4) return 1.01099;
         else if (x<4.5) return 0.999804;
         else if (x<5) return 1.00963;
         else if (x<6) return 1.01372;
         else if (x<7) return 1.00494;
         else if (x<8) return 1.02321;
         else if (x<9) return 1.00911;
         else if (x<14) return 1.01161;
         else return 1.03605;
      } else if (fabs(eta)<2.1) {
         // 1.8 < |eta| < 2.1
         if (x<2.5) return 1.03342;
         else if (x<3) return 1.00154;
         else if (x<3.5) return 1.00246;
         else if (x<4) return 1.01347;
         else if (x<4.5) return 0.997808;
         else if (x<5) return 1.01391;
         else if (x<6) return 1.01296;
         else if (x<7) return 1.01524;
         else if (x<10) return 0.992407;
         else return 1.0035;
      } else {
         // 2.1 < |eta| < 2.4
         if (x<2.5) return 1.02921;
         else if (x<3) return 1.00856;
         else if (x<3.5) return 0.989795;
         else if (x<4) return 0.96595;
         else if (x<4.5) return 0.974834;
         else if (x<5) return 0.975238;
         else if (x<6) return 1.01298;
         else if (x<7) return 1.01219;
         else if (x<10) return 1.00975;
         else return 1.00634;
      }
   }

   // denominator (from MC)
   double den=1;
   if (fabs(eta)<1.2) den = (0.58746*TMath::Erf((x-3.17585)/0.96599))+0.39714;
   else if (fabs(eta)<1.8) den = (0.51632*TMath::Erf((x-0.56984)/1.78604))+0.45932;
   else if (fabs(eta)<2.1) den = (0.51122*TMath::Erf((x-0.01449)/1.43367))+0.47344;
   else den = (0.43706*TMath::Erf((x-0.50936)/1.62531))+0.54785;

   // numerator (from data)
   double num=1;
   if (fabs(eta)<1.2)
   {
      if (idx==0) num = (0.52806*TMath::Erf((x-3.14214)/1.01709))+0.45965;
      else if (idx == -1   ) num = (0.53145*TMath::Erf((x-3.13538)/1.00320))+0.46288;
      else if (idx == -2   ) num = (0.52425*TMath::Erf((x-3.14802)/1.03003))+0.45617;
      else if (idx == 1   ) num = (0.52732*TMath::Erf((x-3.00821)/1.18646))+0.46382;
      else if (idx == 2   ) num = (0.52542*TMath::Erf((x-3.10948)/1.07637))+0.45856;
      else if (idx == 3   ) num = (0.52817*TMath::Erf((x-3.18153)/0.96716))+0.45842;
      else if (idx == 4   ) num = (0.52530*TMath::Erf((x-3.15326)/0.97512))+0.45642;
      else if (idx == 5   ) num = (0.52946*TMath::Erf((x-3.15344)/1.02933))+0.46068;
      else if (idx == 6   ) num = (0.52722*TMath::Erf((x-3.08052)/1.09655))+0.46092;
      else if (idx == 7   ) num = (0.52754*TMath::Erf((x-3.20101)/0.89540))+0.45682;
      else if (idx == 8   ) num = (0.52511*TMath::Erf((x-3.08851)/1.07923))+0.45888;
      else if (idx == 9   ) num = (0.52922*TMath::Erf((x-3.17724)/0.98635))+0.45965;
      else if (idx == 10  ) num = (0.52734*TMath::Erf((x-3.09522)/1.10036))+0.46086;
      else if (idx == 11  ) num = (0.52746*TMath::Erf((x-3.00920)/1.25762))+0.46420;
      else if (idx == 12  ) num = (0.52841*TMath::Erf((x-3.24713)/0.83499))+0.45586;
      else if (idx == 13  ) num = (0.52637*TMath::Erf((x-3.05892)/1.11269))+0.46097;
      else if (idx == 14  ) num = (0.52859*TMath::Erf((x-3.20595)/0.93288))+0.45771;
      else if (idx == 15  ) num = (0.52808*TMath::Erf((x-3.19212)/0.94549))+0.45802;
      else if (idx == 16  ) num = (0.52798*TMath::Erf((x-3.13891)/1.01023))+0.45964;
      else if (idx == 17  ) num = (0.52796*TMath::Erf((x-3.16654)/0.98850))+0.45883;
      else if (idx == 18  ) num = (0.52837*TMath::Erf((x-3.14532)/1.02015))+0.45986;
      else if (idx == 19  ) num = (0.52768*TMath::Erf((x-3.13564)/1.00829))+0.45944;
      else if (idx == 20  ) num = (0.52767*TMath::Erf((x-3.10299)/1.05943))+0.46069;
      else if (idx == 21  ) num = (0.52563*TMath::Erf((x-3.05365)/1.10282))+0.46033;
      else if (idx == 22  ) num = (0.52688*TMath::Erf((x-3.17277)/0.95607))+0.45747;
      else if (idx == 23  ) num = (0.52522*TMath::Erf((x-3.04635)/1.14833))+0.46042;
      else if (idx == 24  ) num = (0.52730*TMath::Erf((x-3.13277)/1.00565))+0.45917;
      else if (idx == 25  ) num = (0.52607*TMath::Erf((x-3.05410)/1.16424))+0.46121;
      else if (idx == 26  ) num = (0.52596*TMath::Erf((x-3.13463)/1.00934))+0.45793;
      else if (idx == 27  ) num = (0.52539*TMath::Erf((x-3.13123)/1.00658))+0.45751;
      else if (idx == 28  ) num = (0.52789*TMath::Erf((x-3.13850)/1.01511))+0.45960;
      else if (idx == 29  ) num = (0.52795*TMath::Erf((x-3.03997)/1.19275))+0.46352;
      else if (idx == 30  ) num = (0.52651*TMath::Erf((x-3.12471)/0.99341))+0.45862;
      else if (idx == 31  ) num = (0.52820*TMath::Erf((x-3.14198)/1.01983))+0.45980;
      else if (idx == 32  ) num = (0.52668*TMath::Erf((x-3.24001)/0.85045))+0.45463;
      else if (idx == 33  ) num = (0.52550*TMath::Erf((x-3.09583)/1.07317))+0.45911;
      else if (idx == 34  ) num = (0.52525*TMath::Erf((x-3.14903)/0.97112))+0.45815;
      else if (idx == 35  ) num = (0.52760*TMath::Erf((x-3.22200)/0.90168))+0.45642;
      else if (idx == 36  ) num = (0.52516*TMath::Erf((x-3.16276)/0.98789))+0.45819;
      else if (idx == 37  ) num = (0.52894*TMath::Erf((x-3.19286)/0.95460))+0.45869;
      else if (idx == 38  ) num = (0.52388*TMath::Erf((x-3.04652)/1.13780))+0.45902;
      else if (idx == 39  ) num = (0.52798*TMath::Erf((x-3.20118)/0.96850))+0.45780;
      else if (idx == 40  ) num = (0.52629*TMath::Erf((x-3.20059)/0.94776))+0.45624;
      else if (idx == 41  ) num = (0.52615*TMath::Erf((x-3.09321)/1.05043))+0.45945;
      else if (idx == 42  ) num = (0.52609*TMath::Erf((x-3.22029)/0.90937))+0.45522;
      else if (idx == 43  ) num = (0.52577*TMath::Erf((x-3.07245)/1.06466))+0.45993;
      else if (idx == 44  ) num = (0.53007*TMath::Erf((x-3.16965)/0.98177))+0.46062;
      else if (idx == 45  ) num = (0.52822*TMath::Erf((x-3.17170)/0.99844))+0.45890;
      else if (idx == 46  ) num = (0.52831*TMath::Erf((x-3.13452)/1.00817))+0.46004;
      else if (idx == 47  ) num = (0.52871*TMath::Erf((x-3.18797)/0.93995))+0.45864;
      else if (idx == 48  ) num = (0.52668*TMath::Erf((x-3.14675)/1.02375))+0.45834;
      else if (idx == 49  ) num = (0.52721*TMath::Erf((x-3.13721)/1.00894))+0.45898;
      else if (idx == 50  ) num = (0.52590*TMath::Erf((x-3.17402)/0.97309))+0.45674;
      else if (idx == 51  ) num = (0.52704*TMath::Erf((x-3.13816)/1.01042))+0.45880;
      else if (idx == 52  ) num = (0.52669*TMath::Erf((x-3.20236)/0.90052))+0.45611;
      else if (idx == 53  ) num = (0.52376*TMath::Erf((x-3.16352)/0.95051))+0.45498;
      else if (idx == 54  ) num = (0.52664*TMath::Erf((x-3.26836)/0.81717))+0.45377;
      else if (idx == 55  ) num = (0.52925*TMath::Erf((x-3.24998)/0.87026))+0.45681;
      else if (idx == 56  ) num = (0.52776*TMath::Erf((x-3.05991)/1.14238))+0.46248;
      else if (idx == 57  ) num = (0.52692*TMath::Erf((x-3.14825)/0.97047))+0.45978;
      else if (idx == 58  ) num = (0.52726*TMath::Erf((x-3.18855)/0.93893))+0.45722;
      else if (idx == 59  ) num = (0.52595*TMath::Erf((x-3.19613)/0.91852))+0.45584;
      else if (idx == 60  ) num = (0.52909*TMath::Erf((x-3.23347)/0.85368))+0.45705;
      else if (idx == 61  ) num = (0.52868*TMath::Erf((x-3.15242)/1.03012))+0.46002;
      else if (idx == 62  ) num = (0.52593*TMath::Erf((x-3.13785)/1.01431))+0.45784;
      else if (idx == 63  ) num = (0.52665*TMath::Erf((x-3.13810)/1.01365))+0.45847;
      else if (idx == 64  ) num = (0.52870*TMath::Erf((x-3.16098)/0.98932))+0.46115;
      else if (idx == 65  ) num = (0.52733*TMath::Erf((x-3.17061)/0.95244))+0.45792;
      else if (idx == 66  ) num = (0.52801*TMath::Erf((x-3.08834)/1.07034))+0.46150;
      else if (idx == 67  ) num = (0.52831*TMath::Erf((x-3.11586)/1.04241))+0.46074;
      else if (idx == 68  ) num = (0.52886*TMath::Erf((x-3.22580)/0.91974))+0.45742;
      else if (idx == 69  ) num = (0.52688*TMath::Erf((x-3.13549)/1.00923))+0.45873;
      else if (idx == 70  ) num = (0.52812*TMath::Erf((x-3.17774)/0.95976))+0.45849;
      else if (idx == 71  ) num = (0.52871*TMath::Erf((x-3.19362)/0.92308))+0.45833;
      else if (idx == 72  ) num = (0.52795*TMath::Erf((x-3.09509)/1.04845))+0.46120;
      else if (idx == 73  ) num = (0.52290*TMath::Erf((x-3.08736)/1.03130))+0.45672;
      else if (idx == 74  ) num = (0.52663*TMath::Erf((x-3.13495)/1.00769))+0.45852;
      else if (idx == 75  ) num = (0.52552*TMath::Erf((x-3.16908)/0.95582))+0.45637;
      else if (idx == 76  ) num = (0.52698*TMath::Erf((x-3.14258)/1.01763))+0.45867;
      else if (idx == 77  ) num = (0.52785*TMath::Erf((x-3.15813)/0.98634))+0.45904;
      else if (idx == 78  ) num = (0.52667*TMath::Erf((x-3.21879)/0.89072))+0.45550;
      else if (idx == 79  ) num = (0.52815*TMath::Erf((x-3.23146)/0.89567))+0.45654;
      else if (idx == 80  ) num = (0.52943*TMath::Erf((x-3.19792)/0.94807))+0.45896;
      else if (idx == 81  ) num = (0.52302*TMath::Erf((x-3.16785)/0.97467))+0.45429;
      else if (idx == 82  ) num = (0.52715*TMath::Erf((x-3.14094)/1.01698))+0.45886;
      else if (idx == 83  ) num = (0.52619*TMath::Erf((x-3.06379)/1.12714))+0.46080;
      else if (idx == 84  ) num = (0.52914*TMath::Erf((x-3.19483)/0.95476))+0.45893;
      else if (idx == 85  ) num = (0.52834*TMath::Erf((x-3.14959)/1.02914))+0.45977;
      else if (idx == 86  ) num = (0.52556*TMath::Erf((x-3.16307)/0.94253))+0.45654;
      else if (idx == 87  ) num = (0.52595*TMath::Erf((x-3.13316)/1.00810))+0.45796;
      else if (idx == 88  ) num = (0.52899*TMath::Erf((x-3.13670)/1.00812))+0.46060;
      else if (idx == 89  ) num = (0.52841*TMath::Erf((x-3.18265)/0.95994))+0.45845;
      else if (idx == 90  ) num = (0.52808*TMath::Erf((x-3.15169)/1.03009))+0.45950;
      else if (idx == 91  ) num = (0.52716*TMath::Erf((x-3.23915)/0.85238))+0.45525;
      else if (idx == 92  ) num = (0.52749*TMath::Erf((x-3.18529)/0.96399))+0.45773;
      else if (idx == 93  ) num = (0.52412*TMath::Erf((x-3.12691)/0.99839))+0.45644;
      else if (idx == 94  ) num = (0.52580*TMath::Erf((x-3.04481)/1.16104))+0.46127;
      else if (idx == 95  ) num = (0.52865*TMath::Erf((x-3.09725)/1.08272))+0.46189;
      else if (idx == 96  ) num = (0.52657*TMath::Erf((x-3.12943)/0.99887))+0.45857;
      else if (idx == 97  ) num = (0.52813*TMath::Erf((x-3.18301)/0.96180))+0.45830;
      else if (idx == 98  ) num = (0.52524*TMath::Erf((x-3.22204)/0.85851))+0.45403;
      else if (idx == 99  ) num = (0.52420*TMath::Erf((x-3.05893)/1.12943))+0.45923;
      else if (idx == 100 ) num = (0.52633*TMath::Erf((x-3.16584)/0.97034))+0.45726;
   }
   else if (fabs(eta)<1.8)
   {
      if (idx==0) num = (0.01753*TMath::Erf((x-5.24739)/9.24400))+0.98259;
      else if (idx == -1   ) num = (0.00345*TMath::Erf((x-8.43141)/9.74621))+0.99431;
      else if (idx == -2   ) num = (0.20268*TMath::Erf((x-0.06569)/3.15975))+0.78340;
      else if (idx == 1   ) num = (0.02123*TMath::Erf((x-2.21280)/5.60816))+0.96816;
      else if (idx == 2   ) num = (0.07298*TMath::Erf((x-0.74511)/3.86904))+0.91953;
      else if (idx == 3   ) num = (0.00008*TMath::Erf((x-0.23320)/3.20845))+0.98116;
      else if (idx == 4   ) num = (0.01284*TMath::Erf((x-2.85194)/2.91393))+0.97533;
      else if (idx == 5   ) num = (0.01024*TMath::Erf((x-4.74096)/1.95632))+0.98025;
      else if (idx == 6   ) num = (0.07849*TMath::Erf((x-1.43892)/1.68762))+0.90076;
      else if (idx == 7   ) num = (0.00536*TMath::Erf((x-6.47107)/0.00603))+0.98578;
      else if (idx == 8   ) num = (0.01119*TMath::Erf((x-9.56506)/9.99459))+0.98739;
      else if (idx == 9   ) num = (0.00297*TMath::Erf((x-9.99732)/9.97056))+0.98424;
      else if (idx == 10  ) num = (0.00092*TMath::Erf((x-3.40668)/3.90063))+0.97684;
      else if (idx == 11  ) num = (0.17787*TMath::Erf((x-2.42570)/0.39688))+0.80253;
      else if (idx == 12  ) num = (0.00000*TMath::Erf((x-5.50635)/9.91138))+0.98182;
      else if (idx == 13  ) num = (0.00006*TMath::Erf((x-2.16945)/8.90821))+0.98502;
      else if (idx == 14  ) num = (0.01142*TMath::Erf((x-4.09115)/1.38966))+0.97578;
      else if (idx == 15  ) num = (0.00631*TMath::Erf((x-5.04864)/0.60484))+0.97723;
      else if (idx == 16  ) num = (0.03447*TMath::Erf((x-1.94682)/2.89737))+0.95615;
      else if (idx == 17  ) num = (0.00072*TMath::Erf((x-8.44164)/9.99718))+0.97435;
      else if (idx == 18  ) num = (0.00216*TMath::Erf((x-6.34180)/9.99970))+0.98002;
      else if (idx == 19  ) num = (0.01930*TMath::Erf((x-2.68888)/2.22092))+0.96763;
      else if (idx == 20  ) num = (0.07824*TMath::Erf((x-2.35046)/0.90841))+0.90876;
      else if (idx == 21  ) num = (0.07431*TMath::Erf((x-0.04659)/4.09791))+0.91602;
      else if (idx == 22  ) num = (0.00470*TMath::Erf((x-9.77349)/0.67485))+0.98877;
      else if (idx == 23  ) num = (0.01892*TMath::Erf((x-0.19310)/1.23145))+0.95943;
      else if (idx == 24  ) num = (0.04034*TMath::Erf((x-1.16846)/5.20853))+0.95131;
      else if (idx == 25  ) num = (0.00057*TMath::Erf((x-0.45482)/0.32997))+0.97534;
      else if (idx == 26  ) num = (0.00685*TMath::Erf((x-2.93533)/9.99594))+0.98459;
      else if (idx == 27  ) num = (0.00000*TMath::Erf((x-8.45186)/5.09233))+0.98098;
      else if (idx == 28  ) num = (0.01979*TMath::Erf((x-6.97256)/9.64327))+0.97954;
      else if (idx == 29  ) num = (0.03926*TMath::Erf((x-2.71520)/0.23755))+0.94710;
      else if (idx == 30  ) num = (0.02550*TMath::Erf((x-2.72328)/1.51738))+0.95471;
      else if (idx == 31  ) num = (0.01655*TMath::Erf((x-9.99304)/9.96405))+0.98043;
      else if (idx == 32  ) num = (0.00700*TMath::Erf((x-9.88581)/9.99487))+0.98603;
      else if (idx == 33  ) num = (0.00957*TMath::Erf((x-6.08627)/9.98522))+0.98657;
      else if (idx == 34  ) num = (0.00460*TMath::Erf((x-7.19233)/9.99964))+0.97524;
      else if (idx == 35  ) num = (0.00000*TMath::Erf((x-6.85833)/9.96900))+0.97882;
      else if (idx == 36  ) num = (0.01833*TMath::Erf((x-0.27358)/8.75652))+0.97590;
      else if (idx == 37  ) num = (0.01063*TMath::Erf((x-5.87262)/9.95885))+0.98109;
      else if (idx == 38  ) num = (0.92605*TMath::Erf((x-2.19442)/0.39761))+0.05928;
      else if (idx == 39  ) num = (0.01128*TMath::Erf((x-4.49926)/9.97057))+0.98323;
      else if (idx == 40  ) num = (0.01215*TMath::Erf((x-4.36385)/0.03160))+0.97481;
      else if (idx == 41  ) num = (0.00000*TMath::Erf((x-4.62980)/9.93975))+0.98127;
      else if (idx == 42  ) num = (0.01316*TMath::Erf((x-3.26365)/0.15704))+0.96935;
      else if (idx == 43  ) num = (0.00521*TMath::Erf((x-0.03623)/0.38553))+0.97682;
      else if (idx == 44  ) num = (0.02152*TMath::Erf((x-3.00242)/0.49726))+0.96400;
      else if (idx == 45  ) num = (0.03872*TMath::Erf((x-0.00691)/6.46206))+0.95523;
      else if (idx == 46  ) num = (0.00750*TMath::Erf((x-6.75129)/9.93606))+0.97827;
      else if (idx == 47  ) num = (0.00376*TMath::Erf((x-3.78416)/8.99843))+0.98098;
      else if (idx == 48  ) num = (0.04502*TMath::Erf((x-1.01609)/4.21667))+0.93773;
      else if (idx == 49  ) num = (0.02707*TMath::Erf((x-0.00118)/9.44463))+0.96760;
      else if (idx == 50  ) num = (0.00962*TMath::Erf((x-6.02643)/9.99919))+0.97619;
      else if (idx == 51  ) num = (0.42449*TMath::Erf((x-0.87898)/1.40779))+0.55945;
      else if (idx == 52  ) num = (0.01550*TMath::Erf((x-4.46179)/9.78694))+0.98374;
      else if (idx == 53  ) num = (0.00685*TMath::Erf((x-7.94336)/9.97777))+0.98141;
      else if (idx == 54  ) num = (0.00687*TMath::Erf((x-4.96863)/9.99857))+0.98158;
      else if (idx == 55  ) num = (0.02018*TMath::Erf((x-6.43505)/6.01923))+0.97942;
      else if (idx == 56  ) num = (0.00000*TMath::Erf((x-0.36676)/2.37384))+0.98214;
      else if (idx == 57  ) num = (0.00000*TMath::Erf((x-9.92242)/9.99939))+0.97309;
      else if (idx == 58  ) num = (0.00668*TMath::Erf((x-9.97428)/9.99977))+0.98347;
      else if (idx == 59  ) num = (0.00988*TMath::Erf((x-9.29244)/0.75776))+0.98906;
      else if (idx == 60  ) num = (0.01449*TMath::Erf((x-5.37556)/0.00593))+0.97489;
      else if (idx == 61  ) num = (0.01614*TMath::Erf((x-7.36244)/9.78174))+0.97844;
      else if (idx == 62  ) num = (0.01818*TMath::Erf((x-2.98224)/0.08243))+0.96564;
      else if (idx == 63  ) num = (0.00020*TMath::Erf((x-9.81342)/0.50523))+0.97746;
      else if (idx == 64  ) num = (0.00756*TMath::Erf((x-4.62848)/0.16948))+0.97370;
      else if (idx == 65  ) num = (0.00000*TMath::Erf((x-1.79164)/8.86880))+0.97914;
      else if (idx == 66  ) num = (0.00042*TMath::Erf((x-6.74663)/9.99658))+0.97908;
      else if (idx == 67  ) num = (0.00312*TMath::Erf((x-3.76840)/4.51264))+0.97982;
      else if (idx == 68  ) num = (0.01730*TMath::Erf((x-4.94032)/7.12698))+0.98331;
      else if (idx == 69  ) num = (0.02059*TMath::Erf((x-0.01048)/9.78680))+0.96953;
      else if (idx == 70  ) num = (0.10724*TMath::Erf((x-0.99565)/0.58495))+0.87941;
      else if (idx == 71  ) num = (0.00515*TMath::Erf((x-9.37961)/3.01867))+0.99175;
      else if (idx == 72  ) num = (0.06026*TMath::Erf((x-0.00048)/5.28534))+0.93200;
      else if (idx == 73  ) num = (0.00668*TMath::Erf((x-8.59713)/0.04927))+0.98436;
      else if (idx == 74  ) num = (0.00603*TMath::Erf((x-4.97847)/9.99977))+0.98153;
      else if (idx == 75  ) num = (0.00000*TMath::Erf((x-9.04356)/9.99461))+0.97762;
      else if (idx == 76  ) num = (0.00445*TMath::Erf((x-8.68461)/0.02187))+0.98634;
      else if (idx == 77  ) num = (0.03153*TMath::Erf((x-2.43133)/1.50034))+0.95343;
      else if (idx == 78  ) num = (0.01121*TMath::Erf((x-6.67375)/2.47009))+0.98573;
      else if (idx == 79  ) num = (0.80029*TMath::Erf((x-2.30555)/0.36271))+0.18062;
      else if (idx == 80  ) num = (0.17088*TMath::Erf((x-0.97699)/1.86668))+0.81245;
      else if (idx == 81  ) num = (0.00000*TMath::Erf((x-7.50633)/9.97689))+0.97833;
      else if (idx == 82  ) num = (0.00185*TMath::Erf((x-4.69943)/9.99992))+0.98482;
      else if (idx == 83  ) num = (0.00169*TMath::Erf((x-9.99764)/9.94047))+0.98101;
      else if (idx == 84  ) num = (0.00275*TMath::Erf((x-9.97517)/9.96644))+0.98804;
      else if (idx == 85  ) num = (0.04852*TMath::Erf((x-2.10701)/1.16846))+0.93910;
      else if (idx == 86  ) num = (0.00102*TMath::Erf((x-3.51616)/9.99134))+0.98391;
      else if (idx == 87  ) num = (0.00000*TMath::Erf((x-7.45945)/9.98449))+0.97997;
      else if (idx == 88  ) num = (0.00305*TMath::Erf((x-9.98948)/9.99131))+0.97964;
      else if (idx == 89  ) num = (0.00936*TMath::Erf((x-4.53171)/0.04084))+0.98181;
      else if (idx == 90  ) num = (0.05330*TMath::Erf((x-0.42369)/4.81723))+0.93541;
      else if (idx == 91  ) num = (0.00010*TMath::Erf((x-0.49921)/3.17135))+0.98043;
      else if (idx == 92  ) num = (0.00413*TMath::Erf((x-5.06056)/9.99450))+0.98216;
      else if (idx == 93  ) num = (0.10634*TMath::Erf((x-0.06331)/0.37948))+0.87904;
      else if (idx == 94  ) num = (0.21357*TMath::Erf((x-1.95245)/0.83414))+0.76931;
      else if (idx == 95  ) num = (0.00806*TMath::Erf((x-3.05679)/9.96035))+0.98511;
      else if (idx == 96  ) num = (0.00000*TMath::Erf((x-2.06480)/9.91022))+0.97691;
      else if (idx == 97  ) num = (0.00237*TMath::Erf((x-9.99984)/9.99880))+0.98051;
      else if (idx == 98  ) num = (0.00904*TMath::Erf((x-7.15537)/4.97318))+0.98257;
      else if (idx == 99  ) num = (0.08340*TMath::Erf((x-0.84331)/2.48669))+0.90425;
      else if (idx == 100 ) num = (0.03027*TMath::Erf((x-1.60165)/2.79258))+0.95561;
   }
   else if (fabs(eta)<2.1)
   {
      if (idx==0) num = (0.52738*TMath::Erf((x-0.00120)/0.58543))+0.46394;
      else if (idx == -1   ) num = (0.53055*TMath::Erf((x-0.00032)/0.57456))+0.46714;
      else if (idx == -2   ) num = (0.51981*TMath::Erf((x-1.37234)/0.48277))+0.45680;
      else if (idx == 1   ) num = (0.51877*TMath::Erf((x-0.95937)/0.46548))+0.45555;
      else if (idx == 2   ) num = (0.52019*TMath::Erf((x-0.00031)/0.57439))+0.45689;
      else if (idx == 3   ) num = (0.51956*TMath::Erf((x-0.00031)/0.57437))+0.45629;
      else if (idx == 4   ) num = (0.52039*TMath::Erf((x-0.00031)/0.57439))+0.45708;
      else if (idx == 5   ) num = (0.52473*TMath::Erf((x-0.00032)/0.57445))+0.46131;
      else if (idx == 6   ) num = (0.52185*TMath::Erf((x-0.00031)/0.57441))+0.45849;
      else if (idx == 7   ) num = (0.52206*TMath::Erf((x-0.98475)/0.48490))+0.45870;
      else if (idx == 8   ) num = (0.52635*TMath::Erf((x-1.22640)/0.57798))+0.46291;
      else if (idx == 9   ) num = (0.52265*TMath::Erf((x-0.00031)/0.57438))+0.45927;
      else if (idx == 10  ) num = (0.51934*TMath::Erf((x-1.11083)/0.42686))+0.45613;
      else if (idx == 11  ) num = (0.52094*TMath::Erf((x-0.97008)/0.49055))+0.45762;
      else if (idx == 12  ) num = (0.52065*TMath::Erf((x-0.96974)/0.46913))+0.45735;
      else if (idx == 13  ) num = (0.51928*TMath::Erf((x-0.00031)/0.57436))+0.45602;
      else if (idx == 14  ) num = (0.52297*TMath::Erf((x-0.00031)/0.57441))+0.45958;
      else if (idx == 15  ) num = (0.52560*TMath::Erf((x-0.00032)/0.57448))+0.46217;
      else if (idx == 16  ) num = (0.52544*TMath::Erf((x-0.00032)/0.57445))+0.46201;
      else if (idx == 17  ) num = (0.52444*TMath::Erf((x-0.00032)/0.57444))+0.46102;
      else if (idx == 18  ) num = (0.51991*TMath::Erf((x-1.14726)/0.32956))+0.45663;
      else if (idx == 19  ) num = (0.52441*TMath::Erf((x-0.00032)/0.57447))+0.46099;
      else if (idx == 20  ) num = (0.52119*TMath::Erf((x-0.00031)/0.57439))+0.45785;
      else if (idx == 21  ) num = (0.52641*TMath::Erf((x-0.00032)/0.57448))+0.46297;
      else if (idx == 22  ) num = (0.52124*TMath::Erf((x-0.00031)/0.57439))+0.45791;
      else if (idx == 23  ) num = (0.52331*TMath::Erf((x-0.00031)/0.57415))+0.45991;
      else if (idx == 24  ) num = (0.52498*TMath::Erf((x-0.00032)/0.57438))+0.46155;
      else if (idx == 25  ) num = (0.52146*TMath::Erf((x-0.00031)/0.57440))+0.45811;
      else if (idx == 26  ) num = (0.52360*TMath::Erf((x-0.00032)/0.57410))+0.46020;
      else if (idx == 27  ) num = (0.51932*TMath::Erf((x-0.93423)/0.45709))+0.45606;
      else if (idx == 28  ) num = (0.52175*TMath::Erf((x-0.98176)/0.47628))+0.45840;
      else if (idx == 29  ) num = (0.51952*TMath::Erf((x-0.93517)/0.46966))+0.45625;
      else if (idx == 30  ) num = (0.52309*TMath::Erf((x-0.00031)/0.57444))+0.45970;
      else if (idx == 31  ) num = (0.52467*TMath::Erf((x-0.00032)/0.57447))+0.46125;
      else if (idx == 32  ) num = (0.52538*TMath::Erf((x-0.00032)/0.57448))+0.46194;
      else if (idx == 33  ) num = (0.51991*TMath::Erf((x-0.00031)/0.57437))+0.45663;
      else if (idx == 34  ) num = (0.52042*TMath::Erf((x-1.11789)/0.43008))+0.45717;
      else if (idx == 35  ) num = (0.51948*TMath::Erf((x-0.00031)/0.57437))+0.45622;
      else if (idx == 36  ) num = (0.51865*TMath::Erf((x-0.97828)/0.46262))+0.45543;
      else if (idx == 37  ) num = (0.52495*TMath::Erf((x-0.00032)/0.57441))+0.46152;
      else if (idx == 38  ) num = (0.51404*TMath::Erf((x-1.06130)/0.37719))+0.45115;
      else if (idx == 39  ) num = (0.52158*TMath::Erf((x-0.00031)/0.57438))+0.45824;
      else if (idx == 40  ) num = (0.52467*TMath::Erf((x-0.00032)/0.57447))+0.46124;
      else if (idx == 41  ) num = (0.52213*TMath::Erf((x-0.00031)/0.57423))+0.45876;
      else if (idx == 42  ) num = (0.52410*TMath::Erf((x-0.00032)/0.57425))+0.46069;
      else if (idx == 43  ) num = (0.52116*TMath::Erf((x-0.00031)/0.57441))+0.45783;
      else if (idx == 44  ) num = (0.52132*TMath::Erf((x-0.00031)/0.57437))+0.45798;
      else if (idx == 45  ) num = (0.52111*TMath::Erf((x-0.00031)/0.57440))+0.45777;
      else if (idx == 46  ) num = (0.52325*TMath::Erf((x-0.53684)/0.33949))+0.45985;
      else if (idx == 47  ) num = (0.51999*TMath::Erf((x-0.00031)/0.57436))+0.45670;
      else if (idx == 48  ) num = (0.52562*TMath::Erf((x-0.85114)/0.50692))+0.46218;
      else if (idx == 49  ) num = (0.52296*TMath::Erf((x-0.00031)/0.57443))+0.45957;
      else if (idx == 50  ) num = (0.52566*TMath::Erf((x-0.00032)/0.57448))+0.46222;
      else if (idx == 51  ) num = (0.52431*TMath::Erf((x-1.02324)/0.50634))+0.46090;
      else if (idx == 52  ) num = (0.52118*TMath::Erf((x-0.00031)/0.57426))+0.45785;
      else if (idx == 53  ) num = (0.52279*TMath::Erf((x-0.73273)/0.45814))+0.45941;
      else if (idx == 54  ) num = (0.52295*TMath::Erf((x-0.79158)/0.04821))+0.45956;
      else if (idx == 55  ) num = (0.51912*TMath::Erf((x-0.00031)/0.57428))+0.45587;
      else if (idx == 56  ) num = (0.51974*TMath::Erf((x-1.10957)/0.40618))+0.45653;
      else if (idx == 57  ) num = (0.52663*TMath::Erf((x-1.13158)/0.57474))+0.46319;
      else if (idx == 58  ) num = (0.51990*TMath::Erf((x-0.00031)/0.57433))+0.45661;
      else if (idx == 59  ) num = (0.52183*TMath::Erf((x-0.73312)/0.38747))+0.45847;
      else if (idx == 60  ) num = (0.52221*TMath::Erf((x-0.00031)/0.57442))+0.45885;
      else if (idx == 61  ) num = (0.52339*TMath::Erf((x-1.00870)/0.49868))+0.46000;
      else if (idx == 62  ) num = (0.52165*TMath::Erf((x-0.00031)/0.57441))+0.45830;
      else if (idx == 63  ) num = (0.52299*TMath::Erf((x-0.00031)/0.57439))+0.45960;
      else if (idx == 64  ) num = (0.52159*TMath::Erf((x-0.00031)/0.57440))+0.45824;
      else if (idx == 65  ) num = (0.52271*TMath::Erf((x-0.00031)/0.57441))+0.45933;
      else if (idx == 66  ) num = (0.52019*TMath::Erf((x-0.35377)/0.30739))+0.45689;
      else if (idx == 67  ) num = (0.52043*TMath::Erf((x-1.36689)/0.47622))+0.45740;
      else if (idx == 68  ) num = (0.52432*TMath::Erf((x-0.00032)/0.57443))+0.46091;
      else if (idx == 69  ) num = (0.52160*TMath::Erf((x-0.00031)/0.57406))+0.45825;
      else if (idx == 70  ) num = (0.52214*TMath::Erf((x-0.73918)/0.43745))+0.45877;
      else if (idx == 71  ) num = (0.52642*TMath::Erf((x-0.00032)/0.57450))+0.46297;
      else if (idx == 72  ) num = (0.52588*TMath::Erf((x-0.00032)/0.57448))+0.46244;
      else if (idx == 73  ) num = (0.52047*TMath::Erf((x-0.00031)/0.57437))+0.45716;
      else if (idx == 74  ) num = (0.52033*TMath::Erf((x-0.00031)/0.57454))+0.45703;
      else if (idx == 75  ) num = (0.52265*TMath::Erf((x-0.99006)/0.48302))+0.45928;
      else if (idx == 76  ) num = (0.52443*TMath::Erf((x-1.06422)/0.54060))+0.46101;
      else if (idx == 77  ) num = (0.51846*TMath::Erf((x-0.00031)/0.57425))+0.45525;
      else if (idx == 78  ) num = (0.52100*TMath::Erf((x-0.00031)/0.57435))+0.45768;
      else if (idx == 79  ) num = (0.52339*TMath::Erf((x-1.03258)/0.52512))+0.46000;
      else if (idx == 80  ) num = (0.52042*TMath::Erf((x-0.00031)/0.57437))+0.45711;
      else if (idx == 81  ) num = (0.52567*TMath::Erf((x-0.00032)/0.57447))+0.46224;
      else if (idx == 82  ) num = (0.52411*TMath::Erf((x-0.81900)/0.49001))+0.46069;
      else if (idx == 83  ) num = (0.52404*TMath::Erf((x-0.99263)/0.48609))+0.46063;
      else if (idx == 84  ) num = (0.52064*TMath::Erf((x-0.96857)/0.47125))+0.45734;
      else if (idx == 85  ) num = (0.52259*TMath::Erf((x-0.00031)/0.57443))+0.45921;
      else if (idx == 86  ) num = (0.52390*TMath::Erf((x-1.00840)/0.49325))+0.46050;
      else if (idx == 87  ) num = (0.52390*TMath::Erf((x-0.00032)/0.57440))+0.46049;
      else if (idx == 88  ) num = (0.52006*TMath::Erf((x-0.80364)/0.49253))+0.45677;
      else if (idx == 89  ) num = (0.52737*TMath::Erf((x-1.12955)/0.57368))+0.46392;
      else if (idx == 90  ) num = (0.52483*TMath::Erf((x-0.00032)/0.57447))+0.46140;
      else if (idx == 91  ) num = (0.52471*TMath::Erf((x-0.60374)/0.45598))+0.46129;
      else if (idx == 92  ) num = (0.52402*TMath::Erf((x-0.00032)/0.57443))+0.46061;
      else if (idx == 93  ) num = (0.00003*TMath::Erf((x-9.82161)/3.05417))+0.99437;
      else if (idx == 94  ) num = (0.52142*TMath::Erf((x-0.00031)/0.57439))+0.45807;
      else if (idx == 95  ) num = (0.52203*TMath::Erf((x-0.87043)/0.49791))+0.45867;
      else if (idx == 96  ) num = (0.52302*TMath::Erf((x-0.66262)/0.45686))+0.45963;
      else if (idx == 97  ) num = (0.51874*TMath::Erf((x-0.98373)/0.46248))+0.45552;
      else if (idx == 98  ) num = (0.52210*TMath::Erf((x-0.00031)/0.57442))+0.45873;
      else if (idx == 99  ) num = (0.52014*TMath::Erf((x-0.00031)/0.57432))+0.45685;
      else if (idx == 100 ) num = (0.52173*TMath::Erf((x-0.00031)/0.57442))+0.45837;
   }
   else
   {
      if (idx==0) num = (0.52196*TMath::Erf((x-0.00031)/1.71500))+0.45823;
      else if (idx == -1   ) num = (0.52518*TMath::Erf((x-0.00001)/0.94216))+0.46239;
      else if (idx == -2   ) num = (0.51944*TMath::Erf((x-0.97798)/1.06998))+0.44553;
      else if (idx == 1   ) num = (0.65158*TMath::Erf((x-1.11230)/0.73330))+0.30682;
      else if (idx == 2   ) num = (0.13140*TMath::Erf((x-0.01972)/4.79355))+0.86149;
      else if (idx == 3   ) num = (0.51388*TMath::Erf((x-0.00011)/1.23801))+0.45091;
      else if (idx == 4   ) num = (0.41467*TMath::Erf((x-0.06535)/1.73932))+0.55160;
      else if (idx == 5   ) num = (0.51476*TMath::Erf((x-0.00009)/1.29729))+0.45177;
      else if (idx == 6   ) num = (0.02811*TMath::Erf((x-4.40424)/0.46762))+0.94925;
      else if (idx == 7   ) num = (0.51314*TMath::Erf((x-0.00269)/1.56981))+0.45013;
      else if (idx == 8   ) num = (0.51326*TMath::Erf((x-1.54979)/0.41108))+0.44631;
      else if (idx == 9   ) num = (0.51787*TMath::Erf((x-0.04374)/1.64907))+0.45438;
      else if (idx == 10  ) num = (0.51226*TMath::Erf((x-0.00015)/1.17491))+0.44942;
      else if (idx == 11  ) num = (0.93320*TMath::Erf((x-1.40057)/0.49704))+0.04443;
      else if (idx == 12  ) num = (0.52308*TMath::Erf((x-1.60571)/0.41712))+0.45428;
      else if (idx == 13  ) num = (0.51946*TMath::Erf((x-0.00009)/1.49651))+0.45606;
      else if (idx == 14  ) num = (0.52254*TMath::Erf((x-1.60456)/0.39333))+0.44965;
      else if (idx == 15  ) num = (0.51573*TMath::Erf((x-0.00113)/1.71759))+0.45398;
      else if (idx == 16  ) num = (0.51085*TMath::Erf((x-0.00057)/1.38865))+0.44795;
      else if (idx == 17  ) num = (0.51781*TMath::Erf((x-0.00010)/1.51935))+0.45439;
      else if (idx == 18  ) num = (0.51316*TMath::Erf((x-0.00020)/1.31865))+0.45014;
      else if (idx == 19  ) num = (0.51294*TMath::Erf((x-0.00348)/1.44121))+0.44979;
      else if (idx == 20  ) num = (0.51286*TMath::Erf((x-0.00071)/1.47696))+0.44976;
      else if (idx == 21  ) num = (0.51637*TMath::Erf((x-1.66219)/0.36424))+0.44717;
      else if (idx == 22  ) num = (0.51243*TMath::Erf((x-0.00016)/0.77164))+0.45004;
      else if (idx == 23  ) num = (0.52195*TMath::Erf((x-0.00000)/1.36475))+0.45862;
      else if (idx == 24  ) num = (0.04780*TMath::Erf((x-9.92774)/8.47034))+0.98426;
      else if (idx == 25  ) num = (0.51829*TMath::Erf((x-0.00006)/1.43923))+0.45493;
      else if (idx == 26  ) num = (0.51769*TMath::Erf((x-0.00113)/1.67750))+0.45426;
      else if (idx == 27  ) num = (0.51315*TMath::Erf((x-0.45679)/1.51073))+0.44848;
      else if (idx == 28  ) num = (0.51210*TMath::Erf((x-0.00060)/1.46649))+0.44912;
      else if (idx == 29  ) num = (0.51580*TMath::Erf((x-0.00000)/1.51162))+0.45254;
      else if (idx == 30  ) num = (0.51478*TMath::Erf((x-0.00088)/1.47608))+0.45150;
      else if (idx == 31  ) num = (0.02112*TMath::Erf((x-6.47299)/2.27533))+0.98159;
      else if (idx == 32  ) num = (0.51785*TMath::Erf((x-0.00152)/1.65385))+0.45434;
      else if (idx == 33  ) num = (0.51807*TMath::Erf((x-0.00011)/1.50878))+0.45474;
      else if (idx == 34  ) num = (0.03416*TMath::Erf((x-9.62949)/9.99901))+0.97408;
      else if (idx == 35  ) num = (0.03485*TMath::Erf((x-9.90403)/9.51530))+0.98199;
      else if (idx == 36  ) num = (0.58230*TMath::Erf((x-1.56927)/0.41627))+0.39202;
      else if (idx == 37  ) num = (0.51206*TMath::Erf((x-0.49001)/1.27927))+0.44740;
      else if (idx == 38  ) num = (0.03530*TMath::Erf((x-4.95109)/0.75042))+0.95920;
      else if (idx == 39  ) num = (0.51949*TMath::Erf((x-0.00006)/1.52838))+0.45606;
      else if (idx == 40  ) num = (0.51755*TMath::Erf((x-0.00018)/1.58100))+0.45412;
      else if (idx == 41  ) num = (0.51759*TMath::Erf((x-0.00000)/1.85016))+0.45602;
      else if (idx == 42  ) num = (0.51894*TMath::Erf((x-0.00024)/1.59899))+0.45545;
      else if (idx == 43  ) num = (0.51885*TMath::Erf((x-0.00039)/1.26987))+0.45576;
      else if (idx == 44  ) num = (0.07124*TMath::Erf((x-1.87602)/9.01167))+0.92698;
      else if (idx == 45  ) num = (0.51681*TMath::Erf((x-0.00007)/1.34320))+0.45362;
      else if (idx == 46  ) num = (0.51484*TMath::Erf((x-0.00069)/1.47100))+0.45163;
      else if (idx == 47  ) num = (0.51389*TMath::Erf((x-0.00046)/0.83783))+0.45132;
      else if (idx == 48  ) num = (0.52610*TMath::Erf((x-1.66186)/0.39669))+0.45612;
      else if (idx == 49  ) num = (0.51383*TMath::Erf((x-0.00011)/0.75089))+0.45133;
      else if (idx == 50  ) num = (0.52412*TMath::Erf((x-0.00002)/1.59718))+0.46056;
      else if (idx == 51  ) num = (0.51872*TMath::Erf((x-0.00006)/1.52323))+0.45530;
      else if (idx == 52  ) num = (0.01853*TMath::Erf((x-5.48158)/0.34381))+0.97339;
      else if (idx == 53  ) num = (0.51475*TMath::Erf((x-0.00026)/1.40703))+0.45158;
      else if (idx == 54  ) num = (0.13880*TMath::Erf((x-0.55029)/4.03962))+0.84579;
      else if (idx == 55  ) num = (0.51762*TMath::Erf((x-0.54838)/1.38486))+0.45186;
      else if (idx == 56  ) num = (0.51130*TMath::Erf((x-0.00014)/0.70692))+0.44876;
      else if (idx == 57  ) num = (0.51474*TMath::Erf((x-0.06578)/1.69583))+0.45173;
      else if (idx == 58  ) num = (0.51878*TMath::Erf((x-0.00018)/1.39072))+0.45549;
      else if (idx == 59  ) num = (0.54114*TMath::Erf((x-1.53529)/0.45255))+0.43299;
      else if (idx == 60  ) num = (0.52454*TMath::Erf((x-1.63488)/0.41427))+0.45343;
      else if (idx == 61  ) num = (0.08065*TMath::Erf((x-1.81474)/4.52511))+0.91779;
      else if (idx == 62  ) num = (0.51655*TMath::Erf((x-0.00072)/1.56834))+0.45312;
      else if (idx == 63  ) num = (0.51719*TMath::Erf((x-0.00010)/1.46958))+0.45385;
      else if (idx == 64  ) num = (0.51246*TMath::Erf((x-0.05341)/0.57703))+0.45008;
      else if (idx == 65  ) num = (0.51646*TMath::Erf((x-0.00015)/1.47931))+0.45308;
      else if (idx == 66  ) num = (0.51812*TMath::Erf((x-0.00010)/1.36804))+0.45491;
      else if (idx == 67  ) num = (0.52001*TMath::Erf((x-0.92896)/1.00961))+0.45417;
      else if (idx == 68  ) num = (0.51292*TMath::Erf((x-0.00019)/0.85289))+0.45023;
      else if (idx == 69  ) num = (0.50982*TMath::Erf((x-0.00053)/1.18656))+0.44705;
      else if (idx == 70  ) num = (0.51808*TMath::Erf((x-0.00014)/1.57221))+0.45467;
      else if (idx == 71  ) num = (0.51328*TMath::Erf((x-0.00013)/1.28801))+0.45029;
      else if (idx == 72  ) num = (0.51483*TMath::Erf((x-0.00015)/1.27967))+0.45178;
      else if (idx == 73  ) num = (0.51681*TMath::Erf((x-1.64161)/0.43335))+0.44791;
      else if (idx == 74  ) num = (0.51502*TMath::Erf((x-0.00010)/1.34015))+0.45187;
      else if (idx == 75  ) num = (0.52054*TMath::Erf((x-0.00004)/1.51463))+0.45706;
      else if (idx == 76  ) num = (0.51746*TMath::Erf((x-1.44169)/0.61086))+0.44972;
      else if (idx == 77  ) num = (0.51008*TMath::Erf((x-0.04970)/0.54049))+0.44757;
      else if (idx == 78  ) num = (0.51388*TMath::Erf((x-0.00017)/1.01973))+0.45104;
      else if (idx == 79  ) num = (0.51576*TMath::Erf((x-0.00023)/1.47230))+0.45251;
      else if (idx == 80  ) num = (0.51453*TMath::Erf((x-0.00003)/1.61643))+0.45187;
      else if (idx == 81  ) num = (0.89680*TMath::Erf((x-0.21475)/1.31634))+0.06023;
      else if (idx == 82  ) num = (0.51753*TMath::Erf((x-0.00010)/1.51397))+0.45416;
      else if (idx == 83  ) num = (0.51639*TMath::Erf((x-0.01292)/1.54253))+0.45304;
      else if (idx == 84  ) num = (0.52343*TMath::Erf((x-0.42147)/1.44618))+0.45521;
      else if (idx == 85  ) num = (0.51805*TMath::Erf((x-0.00009)/1.57594))+0.45466;
      else if (idx == 86  ) num = (0.20201*TMath::Erf((x-0.01147)/2.93439))+0.78637;
      else if (idx == 87  ) num = (0.50736*TMath::Erf((x-0.00005)/0.75234))+0.44502;
      else if (idx == 88  ) num = (0.52209*TMath::Erf((x-0.00000)/1.68505))+0.45842;
      else if (idx == 89  ) num = (0.51186*TMath::Erf((x-0.02469)/1.73065))+0.45003;
      else if (idx == 90  ) num = (0.51634*TMath::Erf((x-0.00022)/1.72547))+0.45366;
      else if (idx == 91  ) num = (0.04311*TMath::Erf((x-9.90480)/9.88915))+0.99340;
      else if (idx == 92  ) num = (0.51561*TMath::Erf((x-0.00015)/0.89445))+0.45294;
      else if (idx == 93  ) num = (0.51986*TMath::Erf((x-0.00000)/0.79247))+0.45734;
      else if (idx == 94  ) num = (0.52155*TMath::Erf((x-0.00000)/1.37730))+0.45821;
      else if (idx == 95  ) num = (0.52384*TMath::Erf((x-0.10106)/1.84166))+0.46014;
      else if (idx == 96  ) num = (0.50630*TMath::Erf((x-0.00547)/1.27866))+0.44390;
      else if (idx == 97  ) num = (0.50809*TMath::Erf((x-0.00005)/0.98670))+0.44559;
      else if (idx == 98  ) num = (0.51113*TMath::Erf((x-0.01912)/0.54316))+0.44877;
      else if (idx == 99  ) num = (0.51183*TMath::Erf((x-0.03253)/0.53965))+0.44923;
      else if (idx == 100 ) num = (0.52161*TMath::Erf((x-0.00217)/1.70550))+0.45789;
   }

   // return
   return num/den;
}


///////////////////////////////////////////////////
//                 T R K    P P                  //
///////////////////////////////////////////////////

double tnp_weight_trk_pp(int idx) {
   if (idx==-1) return 0.998;
   if (idx==-2) return 0.992;
   return 0.995;
}

#endif //#ifndef tnp_weight_h
