#ifndef DDH_HPP
#define DDH_HPP

#include <string>
#include<vector>

using namespace std;

struct sdatapoint {
  double    min,value,max,delta;
  int       precision;
  string    name,type;
};

struct sddh {
  bool td3;
  vector<sdatapoint> d3;
  bool tdamph;
  sdatapoint damph;
  bool thubbardderivs;
  vector<sdatapoint> hubbardderivs;
  bool thubbards;
  vector<sdatapoint> hubbards;
  bool tvorbes;
  vector<sdatapoint> vorbes;
};

#endif /* DDH_HPP */

