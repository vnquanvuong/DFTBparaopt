#ifndef EREPFIT_H
#define EREPFIT_H

#include <vector>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "molecule.hpp"
#include "element.hpp"

class Erepobj {
protected:
  int score_type,ngrid; 
  double dgrid;
  vector<Molecule>     vmol;
  void remove_comment(char *cstr);
  void ifnotfile(const string filename);
  void readinp(const string inputfile);
  void writeskgen(const string& tmp_dir, const GAGenome& g );
public:
  int  fit_type,skdefversion; 
  int  nelem, ncompound;
  bool lc; 
  double omega;
  string skgen,onecent,twocent,libdir,libadddir,xcfunctional,superposition;
//bool repexist[100][100],addskf,skfexist[100][100];
  bool addskf,addrep;
  string dftbversion,repadddir,outfilename,fgrid,frep_in,repopt;
  double restotS,rescS,resevS; // MSE of residuals (mean signed error)
  double restotU,rescU,resevU; // MUE of residuals (mean unsigned error)
  double restot2,resc2,resev2;  // RMS of residuals (root mean square)
  double restot4,resc4,resev4;  // RMS2 of residuals (root mean square)
  double restot8,resc8,resev8;  // RMS4 of residuals (root mean square)
  double evscaling;
  Erepobj();
  Erepobj(const string inputfile);
  void prepare(const string inputfile);
  double score(int id);
  void get_residual(const GAGenome& g);
  void get_residual();
  bool checkskf(const string filename, int type);
  void makeskf(const GAGenome& g );
  void writeout();
  vector<Element> velem;
  vector<double> atomicenergy;
  vector<string> atomicname;
  vector<Atomparameters> atomparameters;
  vector<string> onecenterparameters;
  vector<string> twocenterparameters;
  ofstream fout;

};

#endif
