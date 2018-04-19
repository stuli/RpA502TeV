#void getEfficiencyUpsilon(int state = 1, bool useDataWeight=true, bool useTnpWeigt=false, int tnpIdx=0) { 





#void getEfficiencyUpsilon(int state = 1, bool useDataWeight=true,
#                         int trgIdx=0,  int trkIdx=0, int muIdx=0,  int staIdx=0 //  id = -100 means no correction on this



root -l -b <<EOF
.L getEfficiencyUpsilon.C++
.q
EOF

# NOMINAL
for state in 1 2 3
do
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',0,0,0,-100,-100)'  #TNP Nominal
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-100,-100)'
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',0,-100,-100,-100,-100)'  #Pure MC
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,-100,-100,-100,-100)'
done

# SYSTEMATICS
#    int trgIdx, trkIdx, muIdx, staIdx
for state in 1 2 3
do
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,-1,0,-100,-100)'   
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,-2,0,-100,-100)'   
    
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,-1,-100,-100)'   
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,-2,-100,-100)'   

    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-1,-100)'   
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-2,-100)'   

    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-100,-1)'   
    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-100,-2)'   

    root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,-10,0,-100,-100)'   

    for (( i=1 ; i<=100 ; i++ ))
    do
        root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,'$i',0,-100,-100)'   # stat for trigger
        root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,'$i',-100)'     # stat for muid
        root -l -b -q 'getEfficiencyUpsilon.C+('$state',1,0,0,-100,'$i')'     # stat for muid
    done


done    




