#!/bin/bash
#$1 = OniaMode (Y state): {1,2,3} or 0 for all, $2 = ispPb (Do you want to run on pPb?) {1,0}, $3 = Identifying name (eg. day_month)
if [[ $2=1 ]]; then
	coll_type='pPb'
else
	coll_type='pp'
fi

mkdir eff_${coll_type}$3

if [[ $1=0 ]]; then
	for i in 1 2 3; do
		cp dimuEff_RpA_ana.C dimuEff_oniaMode$i_${coll_type}_$3.C
		sed -i "" "s/VVV/$i/g;" dimuEff_oniaMode$i_${coll_type}_$3.C
		sed -i "" "s/WWW/$2/g;" dimuEff_oniaMode$i_${coll_type}_$3.C
		sed -i "" "s/XXX/${coll_type}/g;" dimuEff_oniaMode$i_${coll_type}_$3.C
		sed -i "" "s/RpA_ana/oniaMode$i_${coll_type}_$3/g;" dimuEff_oniaMode$i_${coll_type}_$3.C
		sed -i "" "s/TAG/$3/g;" dimuEff_oniaMode$i_${coll_type}_$3.C
		cp dimuEff_oniaMode$i_${coll_type}_$3.C eff_${coll_type}$3/dimuEff_oniaMode$i_${coll_type}_$3.C
		root -l -q -b ./dimuEff_oniaMode$i_${coll_type}_$3.C
		rm dimuEff_oniaMode$i_${coll_type}_$3.C

	done
else
	cp dimuEff_RpA_ana.C dimuEff_oniaMode$1_${coll_type}_$3.C
	sed -i "" "s/VVV/$1/g;" dimuEff_oniaMode$1_${coll_type}_$3.C
	sed -i "" "s/WWW/$2/g;" dimuEff_oniaMode$1_${coll_type}_$3.C
	sed -i "" "s/XXX/${coll_type}/g;" dimuEff_oniaMode$1_${coll_type}_$3.C
	sed -i "" "s/RpA_ana/oniaMode$1_${coll_type}_$3/g;" dimuEff_oniaMode$1_${coll_type}_$3.C
	sed -i "" "s/TAG/$3/g;" dimuEff_oniaMode$1_${coll_type}_$3.C
	cp dimuEff_oniaMode$1_${coll_type}_$3.C eff_${coll_type}$3/dimuEff_oniaMode$1_${coll_type}_$3.C
	root -l -q -b ./dimuEff_oniaMode$1_${coll_type}_$3.C
	#This line is optional ^^^ You can run later. But if you do, don't remove it below!
	rm dimuEff_oniaMode$1_${coll_type}_$3.C

fi 
