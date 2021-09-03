#!/bin/bash

whichRound=R3a

# PA DATA
#regular y bins

./DoFit.sh kPADATA 0.0 30.0 -2.87 -1.93 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -1.93 -1.2 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -1.2 -0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.8 -0.4 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.4 0.0 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $whichRound

./DoFit.sh kPADATA 0.0 30.0 -1.93 -0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.8 0.0 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $whichRound

whichRound=R3b

# PA DATA
#regular y bins
<<COMMENT
./DoFit.sh kPADATA 0.0 30.0 -2.87 -1.93 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -1.93 -1.2 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -1.2 -0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.8 -0.4 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.4 0.0 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.0 0.4 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.4 0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.8 1.2 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 1.2 1.93 0 200 4.0 0 120 0 400 0 $whichRound

./DoFit.sh kPADATA 0.0 30.0 -1.93 -0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 -0.8 0.0 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.0 0.8 0 200 4.0 0 120 0 400 0 $whichRound
./DoFit.sh kPADATA 0.0 30.0 0.8 1.93 0 200 4.0 0 120 0 400 0 $whichRound
COMMENT
