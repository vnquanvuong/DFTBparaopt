#!/bin/bash

export PATH=/home/qv0/Dropbox/projects/dftbopt/utils:$PATH

path=`pwd`

echo " h_h_ "; SplineAnsch -a  1.4:2.72  $path/H-H.skf
echo " c_h_ "; SplineAnsch -a  1.5:3.54  $path/C-H.skf
echo " h_n_ "; SplineAnsch -a  1.5:3.40  $path/H-N.skf
echo " h_o_ "; SplineAnsch -a  1.5:3.29  $path/H-O.skf
echo " c_c_ "; SplineAnsch -a  1.5:4.37  $path/C-C.skf
echo " c_n_ "; SplineAnsch -a  1.5:4.22  $path/C-N.skf
echo " c_o_ "; SplineAnsch -a  1.5:4.11  $path/C-O.skf
echo " n_n_ "; SplineAnsch -a  1.5:4.08  $path/N-N.skf
echo " n_o_ "; SplineAnsch -a  1.5:3.97  $path/N-O.skf
echo " n_o_ "; SplineAnsch -a  1.5:3.97  $path/N-O.skf
echo " o_o_ "; SplineAnsch -a  1.5:3.86  $path/O-O.skf


