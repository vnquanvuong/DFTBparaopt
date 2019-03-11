#ifndef TOOLS_HPP_INCLUDED
#define TOOLS_HPP_INCLUDED

#include <string>

using namespace std;

string itoa(int value, int base);

string formattime(string timestr);

void remove_comment(char *cstr);

void   normv(double v[3]);
double absv(double v[3]);
double dotv(double v1[3],double v2[3]);
void   crossv(double v1[3],double v2[3],double v[3]);
int    randomi(int min, int max);
double randomd(double min, double max);
double polynomial1st(double r0, double r, double s0, double s1);
double polynomial2nd(double r0, double r, double s0, double s1, double s2);
double polynomial3rd(double r0, double r, double s0, double s1, double s2, double s3);
double polynomial4th(double r0, double r, double s0, double s1, double s2, double s3, double s4);
double polynomial5th(double r0, double r, double s0, double s1, double s2, double s3, double s4, double s5);

#endif // TOOLS_HPP_INCLUDED
