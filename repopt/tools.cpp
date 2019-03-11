#include <string>
#include <cmath>
#include <algorithm>
#include "allequations.hpp"

using namespace std;

string itoa(int value, int base) {

  string buf;
  if(base < 2 || base > 16) return buf;

  enum { kMaxDigits = 35 };
  buf.reserve( kMaxDigits );

  int quotient = value;

  do{
    buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
    quotient /= base;
  } while ( quotient );

  if( value < 0) buf += '-';

  reverse( buf.begin(), buf.end() );
  return buf;
}

string formattime(string timestr){
  string tflag;
  if( timestr.substr(0,2) != "00" ) {
    tflag = " hours";
    return timestr + tflag;
  }
  if( timestr.substr(3,2) != "00" ) {
    tflag = " min";
    return timestr.substr(3) + tflag;
  }
  if( timestr.substr(6,1) != "0" ) {
    tflag = " sec";
    return timestr.substr(6) + tflag;
  }
  tflag = " sec";
  return timestr.substr(7) + tflag;
}

void remove_comment(char *cstr){
  int i,ipos;
  string sstr;
  sstr=cstr;
  string::size_type pos = sstr.find('#');
  if (pos!= string::npos){
    ipos=pos;
    for(i=ipos;i<=sstr.length();i++){
      cstr[i]='\0';
    }
  }
}

double absv(double v[3]) {
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void normv(double v[3]) {
  int i;
  double r;
  r=absv(v);
  for(i=0;i<3;i++){
    v[i]=v[i]/r;
  }
}

double dotv(double v1[3],double v2[3]) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void crossv(double v1[3],double v2[3],double v[3]) {
  v[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}
int randomi(int min, int max){
//srand(time(0));
  return min + rand() % (max+1-min);
}

double randomd(double min, double max){
//srand(time(0));
  double random = (double(rand())) / double(RAND_MAX);
  double range = max - min;
  return (random*range) + min;
}
double polynomial1st(double r0, double r, double s0, double s1){
  r=r-r0;
  return s0 + s1*r;
}

double polynomial2nd(double r0, double r, double s0, double s1, double s2){
  double r2;
  r=r-r0;
  r2=r*r;
  return s0 + s1*r + s2*r2;
}

double polynomial3rd(double r0, double r, double s0, double s1, double s2, double s3){
  double r2, r3;
  r=r-r0;
  r2=r*r;
  r3=r2*r;
  return s0 + s1*r + s2*r2 + s3*r3;
}

double polynomial4th(double r0, double r, double s0, double s1, double s2, double s3, double s4){
  double r2, r3, r4;
  r=r-r0;
  r2=r*r;
  r3=r2*r;
  r4=r2*r2;
  return s0 + s1*r + s2*r2 + s3*r3 + s4*r4;
}

double polynomial5th(double r0, double r, double s0, double s1, double s2, double s3, double s4, double s5){
  double r2, r3, r4, r5;
  r=r-r0;
  r2=r*r;
  r3=r2*r;
  r4=r2*r2;
  r5=r4*r;
  return s0 + s1*r + s2*r2 + s3*r3 + s4*r4 + s5*r5;
}

