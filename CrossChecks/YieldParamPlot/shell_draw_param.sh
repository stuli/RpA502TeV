#!/bin/bash

#int draw_bkgparam_rap_comp(TString szAA = "PP", int states =1, int DrawOpt = 2 )
# DrawOpt = 0 : draw all fits
# DrawOpt = 1 : draw jaebeom's fit
# DrawOpt = 2 : draw jared's fit


for (( i=0 ; i<=2 ; i++ ))
do
  for (( j=1 ; j<=3 ; j++ ))
  do
    root -l -q -b 'draw_bkgparam_rap_comp.C+("PA",'$j','$i')'
    root -l -q -b 'draw_bkgparam_rap_comp.C+("PP",'$j','$i')'
    root -l -q -b 'draw_sigparam_rap_comp.C+("PA",'$j','$i')'
    root -l -q -b 'draw_sigparam_rap_comp.C+("PP",'$j','$i')'
    root -l -q -b 'draw_bkgparam_pt_comp.C+("PA",'$j','$i')'
    root -l -q -b 'draw_bkgparam_pt_comp.C+("PP",'$j','$i')'
    root -l -q -b 'draw_sigparam_pt_comp.C+("PA",'$j','$i')'
    root -l -q -b 'draw_sigparam_pt_comp.C+("PP",'$j','$i')'
  done
done



rm *.pcm *.so *_C.d
