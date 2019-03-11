#!/bin/bash

#cp grids_bak/* grids/
#../repopt/repopt  rep1.in  | tee rep1.out

for i in {1..1}; do
  ../repopt/repopt  rep2.in  | tee rep2.out_$i
done

#cp grids_bak/* grids/
#../../gaserepfit/gaserepfit/repfit_v5_ref/gasrepfit repo2.in >& repo2.out_old
#
 cp grids_bak/* grids/
 ../repopt/repopt  rep2.in  >& rep2.out_tnew

#rm pop.final.dat  pop.initial.dat  score.dat
#diff rep1.out rep1.out_ref
