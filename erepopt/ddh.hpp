#ifndef DDH_HPP
#define DDH_HPP

#include <string>
#include <vector>
#include <chrono>

using namespace std;

struct sdatapoint {
  string  name,type;
  int     precision;
  double  min,value,max,delta;
};

struct sddh {
  bool                td3;
  vector<sdatapoint>  d3;
  bool                tdamph;
  sdatapoint          damph;
  bool                thubbardderivs;
  bool                thirdorderfull;
  vector<sdatapoint>  hubbardderivs;
  bool                tdamphver2;
  vector<sdatapoint>  damphver2;
  bool                thubbards;
  vector<sdatapoint>  hubbards;
  bool                tvorbes;
  vector<sdatapoint>  vorbes;
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


#endif /* DDH_HPP */

