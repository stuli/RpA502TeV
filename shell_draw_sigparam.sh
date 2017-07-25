#int draw_sigparam_pt_noerr(TString szAA = "PP", int states =1)
### pt_noerr ###
root -l -b -q 'draw_sigparam_pt_noerr.C+("PP",1)'
root -l -b -q 'draw_sigparam_pt_noerr.C+("PP",2)'
root -l -b -q 'draw_sigparam_pt_noerr.C+("PP",3)'
root -l -b -q 'draw_sigparam_pt_noerr.C+("AA",1)'
root -l -b -q 'draw_sigparam_pt_noerr.C+("AA",2)'
root -l -b -q 'draw_sigparam_pt_noerr.C+("AA",3)'
### rap_noerr ###
root -l -b -q 'draw_sigparam_rap_noerr.C+("PP",1)'
root -l -b -q 'draw_sigparam_rap_noerr.C+("PP",2)'
root -l -b -q 'draw_sigparam_rap_noerr.C+("PP",3)'
root -l -b -q 'draw_sigparam_rap_noerr.C+("AA",1)'
root -l -b -q 'draw_sigparam_rap_noerr.C+("AA",2)'
root -l -b -q 'draw_sigparam_rap_noerr.C+("AA",3)'

#int draw_sigparam_pt(TString szAA = "PP", int states =1)
### pt ###
#root -l -b -q 'draw_sigparam_pt.C+("PP",1)'
#root -l -b -q 'draw_sigparam_pt.C+("PP",2)'
#root -l -b -q 'draw_sigparam_pt.C+("PP",3)'
#root -l -b -q 'draw_sigparam_pt.C+("AA",1)'
#root -l -b -q 'draw_sigparam_pt.C+("AA",2)'
#root -l -b -q 'draw_sigparam_pt.C+("AA",3)'
### rap ### No Input file - Don't use!
##root -l -b -q 'draw_sigparam_rap.C+("PP",1)'
##root -l -b -q 'draw_sigparam_rap.C+("PP",2)'
##root -l -b -q 'draw_sigparam_rap.C+("PP",3)'
##root -l -b -q 'draw_sigparam_rap.C+("AA",1)'
##root -l -b -q 'draw_sigparam_rap.C+("AA",2)'
##root -l -b -q 'draw_sigparam_rap.C+("AA",3)'

