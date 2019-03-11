#!/bin/bash

export PATH=/home/qv0/Dropbox/projects/dftbopt/utils:$PATH

rep2XabSpl ../rep2.out
for file in *.4abSpl; do
  name=${file%.4abSpl}
  a1=${name:0:2}
  a2=${name:2:2}
  if [ $a1 != $a2 ]; then
    cp $file ${a2}${a1}.4abSpl
  fi
done

for file in *.4abSpl; do
  name=${file%.4abSpl}
  a1=${name:0:2}
  a2=${name:2:2}
  e1=${a1//_}
  e1=`echo $e1 | tr a-z A-Z`
  e2=${a2//_}
  e2=`echo $e2 | tr a-z A-Z`
  echo $e1-$e2
  xabSpl2spl $file ../eskf/$e1-$e2.skf 1
  mv ${a1}${a2}.spl $e1-$e2.skf
done
rm *.4abSpl

