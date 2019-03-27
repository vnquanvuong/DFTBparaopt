#/bin/bash

rm -rf 3.3.7.tar.gz eigen3.37 galib247 galib247.tgz DFTBparaopt_on.rc

cd repopt
make clean 
cd ../erepopt
make clean 
cd ../utils 
make clean 
cd ../examples
make clean 
cd ../doc/manual/
make clean 
cd ../../

