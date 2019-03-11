#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <map>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double AA_Bohr  = 0.5291772085936;
const double kcal_H   = 627.509474;

extern bool runtest;

class Molecule {
private:
  vector<string> atomname;
  vector<string> elemname;
  vector<string> elemnameext;
  vector<int>    atomindex;
  MatrixXd coord;
  void writegen(const string& tmp_dir);
  void writedftb(const string& tmp_dir, const string& skf_dir);
public:
  
  string    name;
  double    weight,start,end,gauss,sgrid,error,b0;
  vector<double> dosrefs;
  vector<double> dosrefd;
  vector<double> evref;
  vector<double> occref;
  int       icharge, nev, nelem, natom, nelectron, npelectron; 

  void init(const string& xyzfilename, const int icharge_, const double start_, const double end_, const double gauss_, const string& dosfilename, const double weight_);
  double score(int id, int score_type, const string& dftbversion, bool getresidual);
  void getev(int id, int &pos, double * evarray, const string& dftbversion);
};

#endif /* MOLECULE_HPP */


