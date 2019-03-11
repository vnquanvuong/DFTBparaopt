CXX=g++
CXXFLAGS= -std=c++11 -O3  
INCLUDES= -I../eigen3.37
LDFLAGS=
LIBS=

all: ord2abSpl

ord2abSpl: ord2abSpl.o 
	$(CXX) -o ord2abSpl ord2abSpl.o $(LDFLAGS) $(LIBS) 

ord2abSpl.o: ord2abSpl.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ord2abSpl.cpp  

clean:
	rm -rf ord2abSpl ord2abSpl.o

