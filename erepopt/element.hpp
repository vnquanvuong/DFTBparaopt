#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <string>
#include <vector>
#include <cstring>
#include <cmath>

using namespace std;

class Radius {
public:
  double r, minr, maxr, delta;
  int precision;
};

class Element {
public:
  string name;
  int optype, lmax;
  vector<Radius> radius; 
};

class Atomparameters {
public:
  int power;
  string name,shellresolved;
  vector<string> atomconfig; 
};

#endif /* ELEMENT_HPP */


