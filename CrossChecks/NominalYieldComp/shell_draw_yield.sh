#!/bin/bash

#int draw_yield_rap_comp(TString szAA = "PP", int states =1, int DrawOpt = 2, bool isDrawPub = true )
# DrawOpt = 0 : draw all fits
# DrawOpt = 1 : draw jaebeom's fit
# DrawOpt = 2 : draw jared's fit
# isDrawPub = 1 : with 16-008 result


for (( i=0 ; i<=2 ; i++ ))
do
  for (( j=1 ; j<=3 ; j++ ))
  do
    for ((k=0 ; k<=1 ; k++ ))
    do
      root -l -q -b 'draw_yield_rap_comp.C+("PP",'$j','$i','$k')'
    done
    root -l -q -b 'draw_yield_rap_comp.C+("PA",'$j','$i',0)'
    root -l -q -b 'draw_yield_pt_comp.C+("PP",'$j','$i')'
    root -l -q -b 'draw_yield_pt_comp.C+("PA",'$j','$i')'
  done
done



rm *.pcm *.so *_C.d
