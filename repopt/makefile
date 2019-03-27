CXX=icpc
CXXFLAGS= -std=c++11 -O3 -fpermissive -fPIC #-DEIGEN_USE_MKL_ALL
INCLUDES= -I../eigen3.37 -I../galib247/include
LDFLAGS= -L../galib247/lib
LIBS= -lga #-mkl

all: repopt

repopt: repopt.o allequations.o molecule.o tools.o input.o output.o ga.o allequations.hpp auxiliary.hpp ga.hpp tools.hpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o repopt repopt.o allequations.o molecule.o tools.o input.o output.o ga.o $(LDFLAGS) $(LIBS) 

repopt.o: repopt.cpp allequations.hpp auxiliary.hpp tools.hpp ga.hpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c repopt.cpp  

allequations.o: allequations.cpp allequations.hpp auxiliary.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c allequations.cpp  

molecule.o: molecule.cpp allequations.hpp auxiliary.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c molecule.cpp  

input.o: input.cpp allequations.hpp auxiliary.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c input.cpp  

output.o: output.cpp allequations.hpp auxiliary.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c output.cpp  

tools.o: tools.cpp allequations.hpp tools.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c tools.cpp  

ga.o:  ga.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ga.cpp 

clean:
	rm -rf repopt *.o *.log *.d

