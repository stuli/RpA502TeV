#int draw_yield_pt(TString szAA = "PP", int states =1)
### pt ###
root -l -b -q 'draw_yield_pt.C+("PA",1)'
root -l -b -q 'draw_yield_pt.C+("PA",2)'
root -l -b -q 'draw_yield_pt.C+("PA",3)'
### rap ###
root -l -b -q 'draw_yield_rap.C+("PA",1)'
root -l -b -q 'draw_yield_rap.C+("PA",2)'
root -l -b -q 'draw_yield_rap.C+("PA",3)'
### cent ###
#root -l -b -q 'draw_yield_cent.C+("AA",1)'
#root -l -b -q 'draw_yield_cent.C+("AA",2)' ## same as 1S
#root -l -b -q 'draw_yield_cent.C+("AA",3)'

rm *.so *.pcm *.d
