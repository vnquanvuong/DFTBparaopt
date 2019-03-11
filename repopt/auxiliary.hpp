#ifndef AUXILIARY_HPP_INCLUDED
#define AUXILIARY_HPP_INCLUDED

#include <string>
#include <vector>
#include <chrono>

using namespace std;

struct sdatapoint {
  string  name;
  int     precision;
  double  min,value,max;
};

struct sddh {
  bool                td3;
  vector<sdatapoint>  d3;
  bool                tdamph;
  sdatapoint          damph;
  bool                thubbardderivs;
  vector<sdatapoint>  hubbardderivs;
};

class sdivision {
public:
  double p[2],coeff[10];
};

class spot{
public:
  bool    read;
  int     aa; 
  string  potname;
  string  gridname;
  int     ordspl;
  int     nknots; 
  vector<double>  vr;
  vector<sdivision> division; 
  double  minr,minRbond,max_step,expA,expB,expC;
  int     smooth_order;
};

class Timer {
private:
  using clock_t = chrono::high_resolution_clock;
  using second_t = chrono::duration<double, ratio<1> >;
  chrono::time_point<clock_t> m_beg;
 
public:
  Timer() : m_beg(clock_t::now()) {}
  void reset() {
    m_beg = clock_t::now();
  }
  double elapsed() const {
    return chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
  }
};

#endif /* AUXILIARY_HPP_INCLUDED */

