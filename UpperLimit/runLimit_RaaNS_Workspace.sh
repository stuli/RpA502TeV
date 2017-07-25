#!/bin/bash

#outputDir='results'
poiname='raa3'

mkdir -p results/limit

root -l -b <<EOF
.L Limit_RaaNS_Workspace.C
.q
EOF

echo entering to the loop...

#for pt in '0,2' '2,4' '4,6' '6,8' '8,12' '12,16' '16,30'  '0,8' '8,30' '0,30'
for pt in '4.0,6.0'
do
  #for rap in '0,1.2' '1.2,2.4' '0,2.4'
  for rap in '0,1.2'
  do
    #echo  $outputDir 'Limit_RaaNS_Workspace.C('$pt','$rap','0,200')'
    #root -l -b -q $outputDir 'Limit_RaaNS_Workspace.C++('$pt','$rap','0,200')'
    echo  'Limit_RaaNS_Workspace.C('$pt','$rap','0,200')'
    root -l -b -q  'Limit_RaaNS_Workspace.C++('$pt','$rap','0,200')'
  done
done
mv c1.pdf results/limit/c1_Limit_raa3_pt$pt-y$rap-cent0-200.pdf
mv c1.root results/limit/c1_Limit_raa3_pt$pt-y$rap-cent0-200.root
mv Asym*.root results/limit/
