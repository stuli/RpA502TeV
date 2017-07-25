#int draw_yield_pt(TString szAA = "PP", int states =1)
### pt ###
root -l -b -q 'draw_yield_pt.C+("PP",1)'
root -l -b -q 'draw_yield_pt.C+("PP",2)'
root -l -b -q 'draw_yield_pt.C+("PP",3)'
root -l -b -q 'draw_yield_pt.C+("AA",1)'
root -l -b -q 'draw_yield_pt.C+("AA",2)'
root -l -b -q 'draw_yield_pt.C+("AA",3)'
### rap ###
root -l -b -q 'draw_yield_rap.C+("PP",1)'
root -l -b -q 'draw_yield_rap.C+("PP",2)'
root -l -b -q 'draw_yield_rap.C+("PP",3)'
root -l -b -q 'draw_yield_rap.C+("AA",1)'
root -l -b -q 'draw_yield_rap.C+("AA",2)'
root -l -b -q 'draw_yield_rap.C+("AA",3)'
### cent ###
#root -l -b -q 'draw_yield_cent.C+("AA",1)'
#root -l -b -q 'draw_yield_cent.C+("AA",2)' ## same as 1S
#root -l -b -q 'draw_yield_cent.C+("AA",3)'

