#/bin/bash

CXX="icpc"
#CXX="g++"

if [ ! -f galib247.tgz ]; then
  wget http://lancet.mit.edu/ga/dist/galib247.tgz
fi
rm -rf galib247
tar -xf galib247.tgz
find galib247 -type f -exec perl -pi -w -e "s/float/double/g" {} \;
sed -i'' -e 's|$(INSTALL) $(LIB) $(LIB_DEST_DIR)|$(MKDIR) $(LIB_DEST_DIR)\n\t$(MKDIR) $(HDR_DEST_DIR)\n\t$(INSTALL) $(LIB) $(LIB_DEST_DIR)|' galib247/ga/makefile
sed -i'' -e 's|TMPDIR=/var/tmp|TMPDIR=.|' galib247/makefile
sed -i'' -e 's|DESTDIR=/usr/local|DESTDIR=../|' galib247/makevars 
sed -i'' -e 's|CXXFLAGS    = -g -Wall|CXXFLAGS    = -g -fpermissive -Wall|' galib247/makevars 
if [ $CXX == "g++" ]; then
  echo " " 
elif [ $CXX == "icpc" ]; then
  sed -i'' -e 's|g++|icpc -fPIC|' galib247/makevars 
else
  echo "the compiler $CXX is not supported. Please select g++ or icpc"
  exit
fi
#cd galib247
#make clean
#make install | tee make.log
#cd ../
#
#if [ ! -f 3.3.7.tar.gz ]; then
#  wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
#fi
#rm -rf eigen-eigen-323c052e1731
#tar -xf 3.3.7.tar.gz
#rm -rf eigen3.37
#mv eigen-eigen-323c052e1731 eigen3.37
#
#cd repopt
#if [ $CXX == "g++" ]; then
#  cp makefile.gnu makefile
#elif [ $CXX == "icpc" ]; then
#  cp makefile.intel makefile
#else
#  echo "the compiler $CXX is not supported. Please select g++ or icpc"
#  exit
#fi
#make clean 
#make all | tee make.log
#cd ../
#
#cd erepopt
#if [ $CXX == "g++" ]; then
#  cp makefile.gnu makefile
#elif [ $CXX == "icpc" ]; then
#  cp makefile.intel makefile
#else
#  echo "the compiler $CXX is not supported. Please select g++ or icpc"
#  exit
#fi
#make clean 
#make all | tee make.log
#cd ../
#
#cd utils 
#if [ $CXX == "g++" ]; then
#  cp makefile.gnu makefile
#elif [ $CXX == "icpc" ]; then
#  cp makefile.intel makefile
#else
#  echo "the compiler $CXX is not supported. Please select g++ or icpc"
#  exit
#fi
#make clean 
#make all | tee make.log
#cd ../
