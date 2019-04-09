#ifndef ALLEQUATIONS_HPP_INCLUDED
#define ALLEQUATIONS_HPP_INCLUDED

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <map>
#include <Eigen/LU>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double AA_Bohr  = 0.5291772085936;
const double kcal_H   = 627.509474;

extern bool runtest;

class Molecule {
protected:

private:
  MatrixXd fel;
  MatrixXd fref;

public:
  string    name, dftbin, abbr, finp;
  double    ebind, eel, teel, error, ebindmel;
  int       nelem;
  int       natom,ncontraineda;
  double    eweight, fweight;   // weight for energy and force condition (if =0 no equation is set up)
  bool      pbc;
  int       pbctype;

  vector<string> atomname;
  vector<string> elemname;
  vector<string> elemnameext;
  vector<int>    atomindex;
  vector<int>    contraineda;

  double latvec[3][3];
  MatrixXd coord;
  MatrixXd dist;
  MatrixXd frefmel;

  void init(const string& xyzfilename, const double ediss_, const double eweight_, const double fweight_,
	          const string& dftbin_,     const string abbr_ , const string& dftbversion, const string& finp_);

private:
  void writegen(const string& filename);
};

class Potential {
protected:

private:
  
public:
  struct MolDist {
    public:
      string molname;
      double dist;
  };

  string          potname;
  string          gridname;
  int             ordspl;
  int             nknots;     // number of knots, cutoff included
  vector<double>  vr;
  double          minRbond;
  vector<MolDist> vmoldist;
  double          expA,expB,expC;
  string          tsmooth;
  double          psmooth;
  int             derivsmooth,neighsmooth;
  int             smooth_order;
  int             aa;
  double          max_expand,max_step;

  void init(const string potname_, const string gridname_, const int ordspl_, 
            const string tsmooth_, const double psmooth_, const int derivsmooth_, const int neighsmooth_) {
    using namespace std;
    
    ifstream potfile;
    double rvalue;
    string str;

    potname  = potname_;
    for(int i=0; i<potname.size(); i++) {
      potname[i] = tolower(potname[i]);
    }

    gridname = gridname_;
    ordspl   = ordspl_;
    tsmooth  = tsmooth_;
    psmooth  = psmooth_;
    derivsmooth = derivsmooth_;
    neighsmooth = neighsmooth_;

    if (gridname.substr(0,8)!="autogrid"){
      potfile.open(gridname.c_str());
      for (nknots=0; potfile >> rvalue && potfile.good();nknots++ ){
        // besser mit push_back arbeiten!
        vr.resize(nknots+1);
        vr[nknots] = rvalue/AA_Bohr;
      }
      potfile.close();
    }
   
    if (tsmooth!="no" && neighsmooth > nknots-1){
      cerr << endl << "ERROR: " << endl << endl;
      cerr << "in inputfile: neighbour for smoothing must be in the range of (2) - (nknots-1). "
           << "Erroneous input for potential \"" << potname << "\"." << endl << "exit repopt" << endl << endl;
      exit(1);
    }
      
  }

  void insert_moldist(const string molname, const double dist){
    double eps=1e-4;
    MolDist dummymoldist;
    dummymoldist.molname = molname;
    dummymoldist.dist    = dist;
    if (dist-eps > vmoldist[vmoldist.size()-1].dist) vmoldist.push_back(dummymoldist) ;
    else {
      for (int imd=0;imd<vmoldist.size();imd++){
        if (dist+eps < vmoldist[imd].dist){
          vmoldist.insert( vmoldist.begin()+imd, dummymoldist ) ;
          break;
        }
        else if (dist-eps < vmoldist[imd].dist && molname == vmoldist[imd].molname) { 
          break; 
        }
        else if (dist-eps < vmoldist[imd].dist && molname != vmoldist[imd].molname) { 
          vmoldist.insert( vmoldist.begin()+imd, dummymoldist ) ;
          break; 
        }
      }
    }
  }
  
  void init_vmoldist(){
    vmoldist.resize(1);
    vmoldist[0].dist = -1;
  }

  void remove_dummy_vmoldist(){
    vmoldist.erase( vmoldist.begin() );
  }

  void get_autogrid(){

    if (gridname.substr(8,2)=="A-"){
      get_autogridA();
    } else if (gridname.substr(8,2)=="B-"){
      get_autogridB();
    } else if (gridname.substr(8,2)=="C-"){
      get_autogridC();
    } else{
      cerr << endl << "ERROR: " << endl << endl;
      cerr << gridname <<  ": This type of autogrid is currently not implemented." << endl << endl;
      exit(1);
    }
  }
 
  void get_autogridA(){
    // gridpoints at every first+x*0.1 is set until cutoff 
 
    double eps    = 1e-6;
    double first  = atof( gridname.substr(10,3).data() );
    double cutoff = atof( gridname.substr(14,30).data() );
    nknots = floor((cutoff-first)*10+eps)+1;
    vr.resize(nknots);
    for (int i=0;i<nknots;i++) vr[i] = first+i*0.1;
    if (first+(nknots-1)*0.1+eps < cutoff) {
      vr.push_back(cutoff);
      nknots++;
    }
  }

  void get_autogridB(){
    // create grid with 5 knots for every interatomic distance,
    // 2 knots lower and 3 knots greater than the interatomic distances on a fixed mesh of multiple 0.1
    // skip multiples

    double cutoff = atof( gridname.substr(10,30).data() );
    vr.resize(1);
    vr[0]=cutoff;
    for (int idist=0;idist<vmoldist.size();idist++){
      for (int iknot=-2;iknot<=2;iknot++){
        double dist = ceil( vmoldist[idist].dist * 10 ) * 0.1 + 0.1*iknot;
        if (dist+0.2 >= cutoff) break;
        for (int ir=0;ir<vr.size();ir++){
          if (dist == vr[ir]) break;
          if (dist  < vr[ir]) {
            vr.insert(vr.begin()+ir,dist);
            break;
          }
        }
      }
    }
    nknots = vr.size();
  }

  void get_autogridC(){
    // create grid with 3 knots for every interatomic distance,
    // 2 knots lower and 3 knots greater than the interatomic distances on a fixed mesh of multiple 0.1
    // skip multiples

    double cutoff = atof( gridname.substr(10,30).data() );
    vr.resize(1);
    vr[0]=cutoff;
    for (int idist=0;idist<vmoldist.size();idist++){
      for (int iknot=-1;iknot<=1;iknot++){
        double dist = ceil( vmoldist[idist].dist * 10 ) * 0.1 + 0.1*iknot;
        if (dist+0.2 >= cutoff) break;
        for (int ir=0;ir<vr.size();ir++){
          if (dist == vr[ir]) break;
          if (dist  < vr[ir]) {
            vr.insert(vr.begin()+ir,dist);
            break;
          }
        }
      }
    }
    nknots = vr.size();
  }
};


class Reaction 
{
public:
  vector<double>  vcoeff;
  vector<string>  vabbr;
  vector<int>     vmol;   // index of mol in order to find infos in Reamol
  double          reaE,treaE, error;
  double          reaweight;
  string          reastr;

  void init(const string& reastr_, const vector<Molecule>& vreamol){
    istringstream   ss(reastr_);
    string    str;
    int     nreamol,ind;

    reastr = reastr_;

    // use chemical notation for stochiometric factor!  3 c2h2 - 1 c6h6 -> rea_energy
    // => coeff(c2h2)=-3    coeff(c6h6)=+1
    for (nreamol=0; ss >> str && str != "->" ; nreamol++){
      vcoeff.push_back(-atoi(str.c_str()));
      ss >> str;
      vabbr.push_back(str);

      for (ind=0; ind < vreamol.size(); ind++) 
        if (vabbr[nreamol] == vreamol[ind].abbr ) break;
      if ( ind == vreamol.size() ) {
        cout  << endl << "ERROR: " << endl << endl 
              << "Molecule \"" << vabbr[nreamol] << "\" used in the reaction section needs a definition!" 
              << endl << "exit repopt" << endl << endl;
        exit(1);
      }
      vmol.push_back(ind);
    }
    nreamol++;
    ss >> reaE >> reaweight;
    reaE = reaE/kcal_H;
  }
};

struct Colind {
  public:
  int k,c,p,e;      // knot, column, potential, error
};

class Allequations {
////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                            //
// expression to be minimized:                                                                //
//                                                                                            //
// | eqmat * vector(unknowns) - vector(references) | = 0     (after weighting eqmat and vref) //
//                                                                                            //
//                                                                                            //
//            spline equations            ( here references are 0)                            //
//            energy equations      ( here references are reference dissoziation energies)    //
// eqmat  =   force  equations          ( here references are "-F_el" )                       //
//            ...                                                                             //
//                                                                                            //
//                                                                                            //
// unknowns = spline coefficients, atomic energies                                            //
//                                                                                            //
//    spline:  f_i(x) = a_i^0 + a_i^1 * (x-x_i) + a_i^2 * (x-x_i)^2 + ...                     //
//                                                                                            //
//    spline coefficients: a_1^0(potential 1), a_1^1(potential 1),... for all potentials      //
//                                                                                            //
//    atomic energies: e.g. E_Carbon, E_Hydrogen,...                                          //
//                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////
protected:
  class ExclPot {
    public:
    string potname;
    int    ndist;     // number of distances which could not be calculated
    double mindist;   // minimum of distance b/w the two atoms of the potential
    string molname;   // name where mindist appears
  };
  class AddEq {       // additional equation
    public:
    string potname;   // e.g. c_c_
    string comment;   // to know what this additional equation is for
    int    deriv;     // which derivative 0st,1st,2nd or 3rd (0,1,2,3)
    double dist,val;  // distance and value of deriv
    double weight;    // weight of condition
  };

  //int nrows, ncols, nspleq, neeq, nfeq, naddeq, nreaeq, nsmootheq, nallequations, nelem, natom;
  int natfit;   // number of elements where eatom is to be fitted
  int nzl;      // number of zero-lines in case of not enough equations
  int nfree;    // number of spline parameters that are free for fitting (one per interval)
  string dftbversion;
  vector<Molecule>    vmol;     // sorted by the way they appear in the inputfile
  vector<AddEq>       vadd;     // sorted by the way they appear in the inputfile
  vector<string>      velem;    // sorted by the way they appear
  vector<double>      veatom;   // first entries sorted as in velem (energy of atoms)
  vector<ExclPot>     vexclpot; // in calculation excluded potentials because no input
  vector<Molecule>    vreamol;  // molecules of reactions
  vector<Reaction>    vrea;     // vector of reactions
  
  MatrixXd eqmat;    // rows: as vector<double> vref; cols: as in vector<double> vunknown
  VectorXd vref;     // does not have to be sorted, just equivalent as eqmat
  VectorXd vweight;  // weighing vector
  VectorXd vres;     // residual = vref-eqmat*vunknown
  VectorXd sigma;    // vector with singular values in descending order       
  VectorXd sigma_e_f;// vector with singular values of reduced equation system

  int svderror;
  double cond1, cond2;          // condition number of sigma before and after truncation
  int trunc;                    // number of singular values which are kept, veqmat.size2()-trunc are number of truncated singular values

private:
  void   readinp(const string inputfile);
  void   sort_inputdist();
  void   get_autogrids();
  void   sizeeqsys();
  void   get_minRbond();

  int   include_splineeq(const int ieq_);
  int   include_energyeq(const int ieq_);
  int   include_forceeq (const int ieq_);
  int   include_addeq   (const int ieq_);
  int   include_reaeq   (const int ieq_);
  int   include_smootheq(const int ieq_);
  //void   include_frequencyeq();
  void   weighteqsys();
  void   svd_all();
  void   svd_fulfill_spleq();
  void   get_residual();

  double myfac(const int a, const int b);
  void   add_energy_lhs(const int ieq, const Molecule& mol,const double coeff);
  Colind findcolind(const string  atomname1, const string atomname2, const double dist,
                          const string& molname,   const int at1,          const int at2     );
  string findpotname(const string at1, const string at2) const;
  void   excludepot(const string potname, const double dist, const string& molname);
  void   ifnotfile(const string filename) const;
  void   addextef(const int ieq,const Molecule& mol,const int at1, 
          const int at2,const string ef,const double coeff);

public:
  Allequations();
  Allequations(const string inputfile);
  vector<Potential>   vpot;     // sorted by the way they appear in the inputfile
  VectorXd            vunknown; // first spline coeff, then atomic energies, latter, as sortid in vector<string> velem
  double restotS,reseS,resfS,resaddS,resreaS,ressmoothS ; // MSE of residuals (mean signed error)
  double restotU,reseU,resfU,resaddU,resreaU,ressmoothU ; // MUE of residuals (mean unsigned error)
  double restot2,rese2,resf2,resadd2,resrea2,ressmooth2;  // RMS of residuals (root mean square)
  double restot4,rese4,resf4,resadd4,resrea4,ressmooth4;  // RMS2 of residuals (root mean square)
  double restot8,rese8,resf8,resadd8,resrea8,ressmooth8;  // RMS4 of residuals (root mean square)
  int    nrows, ncols, nspleq, neeq, nfeq, naddeq, nreaeq, nsmootheq, nallequations, nelem, natom;
  void calculate(const string inputfile);
  void prepare(const string inputfile);
  void rereadmol(const string inputfile);
  void reset();
  double score();
  void writeout() const;
};


#endif /* ALLEQUATIONS_HPP_INCLUDED */

