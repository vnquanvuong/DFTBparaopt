#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/bindings/lapack/lapack.hpp>
//#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class OrigRepr {
public:
  int           nspl;
  double         cutoff;
  double         A,B,C;   
  vector<double> splcoeff3;
  vector<double> splcoeff5;

  OrigRepr(){
    nspl    = 0;
    cutoff  = 0;
    A       = 0;
    B       = 0;
    C       = 0;
    splcoeff5.resize(6);
    for (int i=0;i<6;i++){
      splcoeff5[i] = 0;
    }
  }

};

void ifnotfile(const string filename);
OrigRepr calc_orig_mod1(const vector<double> &vr, const vector<double> &va, const double c1, const double cn);
OrigRepr calc_orig_mod2(const vector<double> &vr, const vector<double> &v2nd);
void writeorigabSpline (const vector<double> &vr, const vector<double> &va, const double c1, const double cn, const int mod);
int ilmsfit=2, idecompose=5;
int main(int argc, char** argv){
// writes repulsiv potential to standardout in .abSpline-format (from original dftb 1998)
//
// argv[1]: needs inputfile with splines of nth order in the following format (similar to .abSpline-format)
// argv[2]: also needs order of splines
// argv[3]: modus: 1 is take values and 2nd deriv at first and last point, 2 is second deriv at given r
// argv[4]: optional: grid (knot vecotr) for new spline in the following format (more information than grid is not needed)
//
// ...
// Spline
// NumberOfSplines CutoffDistance
// "something"
// r1 r2 a1 b1 c1 d1 e1 ...
// r2 r3 a2 b2 c2 d2 e2 ...
// ...
// rn CutoffDitance an bn cn dn en ...

  if ( argc < 4 ){
    cerr << "usage: " << argv[0] << "  xspl-file  OrderOfSpline  modus [optional: grid for new spline as .spl-filename or similar]" << endl;
    cerr << "where modus=1 is fitting splines taking function values and second deriv at first and last point and " << endl;
    cerr << "      modus=2 is fitting splines taking only second derivatives and assuming 3rd and 4th deriv be zero at cutoff" << endl;
    cerr << "      grid: default is grid as in xspl-file" << endl;
    cerr << "      NOTE: the grid has to contain more than 3 spline-segments! If within the xspl-file there are defined less than 3," << endl;
    cerr << "            a grid with more than 3 spline-segments is required as argument!" << endl;
    exit(1);
  }

  string   str = argv[1];
  int      ordspl = atoi (argv[2]);
  int     mod = atoi (argv[3]);
  int     nspl;
  double   cutoff;

  // check if modus is correct
  if ( mod != 1 && mod != 2 ) {
    cerr << "Modus " << mod << " is not available." << endl;
    exit(1);
  }
  
  // read xspl-file
  ifnotfile(str);
  ifstream fin(str.c_str());
  while ( getline(fin,str) ){
    if ( str == "Spline" ) break;
    cout << str << endl;
  }

  fin >> nspl >> cutoff;
  getline(fin,str);
  getline(fin,str);
  
  vector<double> vr(nspl+1);
  vector<double> v2ch((ordspl+1)*nspl);
  
  for (int i=0; i < nspl; i++){  
    fin >> vr[i] >> str;
    for (int j=0; j <= ordspl; j++){
      fin >> v2ch[i*(ordspl+1)+j];
    }
    getline(fin,str);
  }

  fin.close();
  vr[nspl] = cutoff;

  // determine parameters to ereprepr::calc_orig
  vector<double> vr_neu;
  vector<double> va;
  int           nspl_neu=0;
  double         c1=0,cn=0,cutoff_neu=0;
  if (argc >= 5){
    str = argv[4];
    ifnotfile(str);
    
    ifstream fin(str.c_str());
    while ( getline(fin,str) ) { if (str == "Spline") break; }
    fin >> nspl_neu >> cutoff_neu;

    if ( ! getline(fin,str) ) {
      cerr << argv[4] << ": not in correct format, there might be one or more spaces after \"Spline\" -> delete these!" << endl;
      exit(1);
    }
    if ( cutoff_neu != cutoff ) {
       cerr << endl 
            << "WARNING: cutoff of grid from file " << argv[4] << " is not identical to cutoff of file " << argv[1] 
            << endl << endl;}

    getline(fin,str);

    // read vr_hlp
    vector<double> vr_hlp(nspl_neu+1);
    for (int i=0; i < nspl_neu; i++){
      fin >> vr_hlp[i];
      getline(fin,str);
    }
    fin.close();
    vr_hlp[nspl_neu] = cutoff_neu;
    
    // create vr_neu 
    int ii; 
    for (ii=0; ii < nspl_neu; ii++){  
      if ( vr_hlp[ii] > vr[0] ) break;
    }
    vr_neu.resize(nspl_neu + 1 - ii+1); 
    vr_neu[0] = vr[0];
    for (int j=1; j < vr_neu.size(); j++){
      vr_neu[j] = vr_hlp[j-1+ii];
    }
    nspl_neu = nspl_neu - ii+1;

    // calculate c1, va, cn
    c1 = v2ch[2];
    va.resize(nspl_neu);
    for (int i=0; i < nspl_neu; i++){
      va[i] = 0;
      for (ii=0; ii<vr.size()-1;ii++){
        if      ( ii == vr.size()-2 && vr[ii] <= vr_neu[i] ) { break;} 
        else if ( vr[ii] == vr_neu[i] ) break;
        else if ( vr[ii] >  vr_neu[i] ) { ii--; break;}
      }
      double delta = vr_neu[i] - vr[ii];
      for (int j=0; j <= ordspl; j++){
        if ( vr_neu[i] >= cutoff ) {va[i] = 0; break;} 
        if (mod == 1) {
          va[i] += v2ch[ii*(ordspl+1)+j] * pow(delta,j);
        }
        if (mod == 2) {
          if ( delta == 0 ) {
            va[i] =  2*v2ch[ii*(ordspl+1)+2];
          }
          else  va[i] += j*(j-1)*v2ch[ii*(ordspl+1)+j] * pow(delta,j-2);
        }
      }
      if ( i == nspl_neu-1 ){  
        for (int j=2; j<=ordspl; j++){
          if ( vr_neu[i] >= cutoff ) { cn = 0; break;} 
          cn += j*(j-1)* v2ch[v2ch.size()-ordspl-1+j] * pow(delta,j-2) / 2;
        }
      }     
    }
  }
  else {
    vr_neu.resize(vr.size());
    vr_neu = vr;
          va.resize(vr.size()-1);
          for (int i=0; i<va.size(); i++){
      if (mod == 1 ) {
        va[i] = v2ch[i*(ordspl+1)];
      }
      if (mod == 2 ) {
        va[i] = v2ch[i*(ordspl+1)+2];
      }
    }
    c1 = v2ch[2];
    cn = v2ch[v2ch.size()-ordspl+1];
        }

  writeorigabSpline(vr_neu,va,c1,cn,mod);

  return 0;
}

void ifnotfile(const string filename) {
  ifstream file;
  file.open(filename.c_str());
  if ( !file ) {
    cerr << filename << ": file not found" << endl;
    exit(1);
  }
  file.close();
}

OrigRepr calc_orig_mod1(const vector<double> &vr, const vector<double> &va, const double c1, const double cn){   
// calculates repulsive potential represented in the original way (dftb 1998)
// 
// vr    : knot vector inclusive cutoff
// va[i] : value of curve to fit at vr[i]
// c1    : second derivative devided by two at the first point of the curve
// cn    : second derivative at vr[vr.size()-2] (one before the cutoff point)
//
// e(x)   = exp(-Ax+B) + C
// S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
// L(x)   = a + b*(x-x_n) + c*(x-x_n)^2 + d*(x-x_n)^3 + e*(x-x_n)^4 + f*(x-x_n)^5
//
// h_i = x_i+1 - x_i
//
// S''_i(x_i+1) = S''_i+1(x_i+1) 
//    ==> d_i = (c_i+1 - c_i) / 3h_i
//
// S_i(x_i+1) = S_i+1(x_i+1) 
//    ==> b_i = (a_i+1 - a_i) / h_i  -  (c_i+1 + 2*c_i)*h_i / 3
// 
// S'_i(x_i+1) = S'_i+1(x_i+1) 
//    ==>
//    ( 2(h_2 + h_1)    h_2         )(c_2  )   ( 3*( (a_3 - a_2)/h_2 - (a_2 - a_1)/h_1) - c_1*h_1  )
//    (   h_2         2(h_3 + h_2)   h_3                      )(c_3  )   ( 3*( (a_4 - a_3)/h_3 - (a_3 - a_2)/h_2)            )
//    (     ..  ..  ..    )(..   )   ( ..                  )
//    (       ..  ..  ..  )(..   ) = ( ..                                                )
//    (       h_n-2   2(h_n-1 + h_n-2))(c_n-1)   (3*((a_n -a_n-1)/h_n-1 -(a_n-1 -a_n-2)/h_1) -c_n*h_n-1)
// with the right side known, all coeffiecients can be calculated
//
// L  (x_n+1) = 0 L  (x_n) = a
// L' (x_n+1) = 0       L' (x_n) = b  = b_n-1 + 2*c_n-1*h_n-1 + 3*d_n-1*(h_n-1)^2
// L''(x_n_1) = 0       L''(x_n) = 2c = 2cn
// calculate a,b,c first, then setup equationsystem for left three conditions and solve
//
// e  (x_1) = a_1
// e' (x_1) = b_1
// e''(x_1) = 2*c_1  
// three conditions for three unknowns A,B,C
// 

  //using namespace boost::numeric::ublas;
  //using namespace boost::numeric::bindings::lapack;
  //namespace bnu = boost::numeric::ublas;

  OrigRepr origrepr;

  origrepr.nspl   = vr.size()-1;
  origrepr.cutoff = vr[vr.size()-1];

  origrepr.splcoeff3.resize(4*(origrepr.nspl-1)); 
  for (int i=0; i < 4*(origrepr.nspl-1); i++){
    origrepr.splcoeff3[i] = 0;
  }

  int cmatsize = origrepr.nspl - 2;
  MatrixXd cmat = MatrixXd::Zero(cmatsize,cmatsize);
  MatrixXd cvec = MatrixXd::Zero(cmatsize,1);
  // for more efficiency change cmat to tridiagonal matrix and solve eqsystem...
  // see Numerical Recipies in Fortran 77 2nd ed, Press, Teukolsky,... p.43
  for (int i=0; i<cmatsize; i++){
    cvec(i,0) = 0;
    for(int j=0; j<cmatsize;j++){
      cmat(i,j) = 0;
    }
  }

  // setup cmat

  double delta1 = vr[1]-vr[0];  
  double delta2 = vr[2]-vr[1];
  double a1     = va[0];
  double a2     = va[1];
  double a3     = va[2];

  cmat(0,0) = 2*(delta1+delta2);
  cmat(0,1) = delta2;
  cvec(0,0) = 3*( (a3-a2)/delta2 - (a2-a1)/delta1 ) - c1*delta1;

  for (int i=1; i < cmatsize-1; i++){
    delta1 = delta2;
    delta2 = vr[i+2]-vr[i+1];
    a1     = a2;
    a2     = a3;
    a3     = va[i+2]; 

    cmat(i,i-1) = delta1;
    cmat(i,i  ) = 2*(delta2+delta1);
    cmat(i,i+1) = delta2;
    cvec(i,0)   = 3*( (a3-a2)/delta2 - (a2-a1)/delta1 ); 
  }

  delta1 = delta2;
  delta2 = vr[cmatsize+1]-vr[cmatsize];
  a1     = a2;
  a2     = a3;
  a3     = va[va.size()-1];

  cmat(cmatsize-1,cmatsize-2) = delta1;
  cmat(cmatsize-1,cmatsize-1) = 2*(delta1+delta2);
  cvec(cmatsize-1,0)          = 3*((a3-a2)/delta2 - (a2-a1)/delta1 ) - cn*delta2;
  
  // solve cmat * cvec_out = cvec_in
  VectorXi ipiv(cmatsize);
  //int error = getrf(cmat,ipiv);
  //error     = getrs(cmat,ipiv,cvec);

  MatrixXd cmatinvcvec;
  if(idecompose==1){
    cmatinvcvec = cmat.ldlt().solve(cvec);
  }else if(idecompose==2){
    cmatinvcvec = cmat.partialPivLu().solve(cvec);
  }else if(idecompose==3){
    cmatinvcvec = cmat.fullPivLu().solve(cvec);
  }else if(idecompose==4){
    cmatinvcvec = cmat.householderQr().solve(cvec);
  }else if(idecompose==5){
    cmatinvcvec = cmat.colPivHouseholderQr().solve(cvec);
  }else if(idecompose==6){
    cmatinvcvec = cmat.fullPivHouseholderQr().solve(cvec);
  }else if(idecompose==7){
    cmatinvcvec = cmat.completeOrthogonalDecomposition().solve(cvec);
  }

  // create splcoeff3
  for (int i=0; i < origrepr.nspl-1; i++){
    a1     = va[i];
    a2     = va[i+1];
    delta1 = vr[i+1]-vr[i];

    origrepr.splcoeff3[4*i  ] = a1;
    if ( i==0 ){
      origrepr.splcoeff3[1] = (a2-a1)/delta1 - (cmatinvcvec(0,0)+2*c1)*delta1/3;
      origrepr.splcoeff3[2] = c1;
      origrepr.splcoeff3[3] = (cmatinvcvec(0,0)-c1)/(delta1*3);
    }
    else if ( i == origrepr.nspl-2){
      origrepr.splcoeff3[4*i+1] = (a2-a1)/delta1 - (cn+2*cmatinvcvec(i-1,0))*delta1/3;
      origrepr.splcoeff3[4*i+2] = cmatinvcvec(i-1,0);
      origrepr.splcoeff3[4*i+3] = (cn-cmatinvcvec(i-1,0))/(delta1*3);
    }
    else {
      origrepr.splcoeff3[4*i+1] = (a2-a1)/delta1 - (cmatinvcvec(i,0)+2*cmatinvcvec(i-1,0))*delta1/3;
      origrepr.splcoeff3[4*i+2] = cmatinvcvec(i-1,0);
      origrepr.splcoeff3[4*i+3] = (cmatinvcvec(i,0)-cmatinvcvec(i-1,0))/(delta1*3);
    }
  }
  
  // now create splcoeff5
  double d1,d2,d3,d4,d5;
  d1 = vr[vr.size()-1] - vr[vr.size()-2];
  d2 = d1*d1;
  d3 = d2*d1;
  d4 = d2*d2;
  d5 = d4*d1;
  double bnm1,cnm1,dnm1;
  delta1 = vr[vr.size()-2] - vr[vr.size()-3];
  bnm1   = origrepr.splcoeff3[4*(origrepr.nspl-2)+1];    
  cnm1   = origrepr.splcoeff3[4*(origrepr.nspl-2)+2];    
  dnm1   = origrepr.splcoeff3[4*(origrepr.nspl-2)+3];    
  
  origrepr.splcoeff5[0] = va[va.size()-1];
  origrepr.splcoeff5[1] = bnm1 + 2*cnm1*delta1 + 3*dnm1*delta1*delta1;
  origrepr.splcoeff5[2] = cn;
  double a = origrepr.splcoeff5[0];
  double b = origrepr.splcoeff5[1];
  double c = origrepr.splcoeff5[2];
  origrepr.splcoeff5[3] = -10/d3 * a - 6/d2 * b - 3/d1 * c;
  origrepr.splcoeff5[4] =  15/d4 * a + 8/d3 * b + 3/d2 * c;
  origrepr.splcoeff5[5] = - 6/d5 * a - 3/d4 * b - 1/d3 * c;

  
  // now create A,B,C - exponential coefficients
  a = origrepr.splcoeff3[0];
  b = origrepr.splcoeff3[1];
  c = origrepr.splcoeff3[2];

  origrepr.A = -2*c/b;
  origrepr.B = log(-b/origrepr.A) + origrepr.A*vr[0];
  origrepr.C = a-exp(-origrepr.A*vr[0]+origrepr.B);

  return origrepr;
}


OrigRepr calc_orig_mod2(const vector<double> &vr, const vector<double> &v2nd){   
//------------------- 
//        dr      = cutoff - vr[n-1];                             // fifth order spline in the end y = ax5 + bx4 + cx3 + dx2 + ex + f
//        a       = - veneu[n-1] / ( 20 * pow(dr, 3) );
//        b       =   veneu[n-1] / ( 4  * dr * dr    );
//        c       = - veneu[n-1] / ( 2  * dr         );
//        d       =   veneu[n-1] / ( 2               );
//        vc[n-1] = - veneu[n-1] * dr / 4;                        // = e
//        vd[n-1] =   veneu[n-1] * dr * dr / 20;                  // = f
//
//
//        for ( int i=n-2; i>=0 ; i--) {                          // cubic splines before y = ax3 + bx2 + cx + d
//                dr    =  vr[i+1]    - vr[i];
//                va[i] = (veneu[i+1] - veneu[i]) / (6*dr);
//                vb[i] =  veneu[i  ]             /  2;
//                vc[i] =  vc   [i+1] - va[i] * 3 * dr*dr  - vb[i] * 2  * dr;
//                vd[i] =  vd   [i+1] - va[i] * pow(dr, 3) - vb[i] * dr * dr - vc[i] * dr;
//        }
//
//        A = -2*vb[0]/vc[0];                                     // exponential fit on the left margin
//        B = log(-vc[0]/A) + A*vr[0];
//        C = vd[0]-exp(-A*vr[0]+B);
//-------------------

  ///// initialization /////
  OrigRepr origrepr;

  int n          = vr.size() -1;
  origrepr.nspl   = vr.size()-1;
  origrepr.cutoff = vr[vr.size()-1];

  origrepr.splcoeff3.resize(4*(n-1)); 
  for (int i=0; i < 4*(n-1); i++){
    origrepr.splcoeff3[i] = 0;
  }

  ///// fifth order spline in the end  y = a + bx + cx2 + dx3 + ex4 + fx5 /////
  double dr = origrepr.cutoff - vr[n-1];  
  origrepr.splcoeff5[0] =   v2nd[n-1] *   dr * dr / 20   ;
  origrepr.splcoeff5[1] = - v2nd[n-1] *   dr      /  4   ;
  origrepr.splcoeff5[2] =   v2nd[n-1] /    2             ;
  origrepr.splcoeff5[3] = - v2nd[n-1] / (  2 * dr       );
  origrepr.splcoeff5[4] =   v2nd[n-1] / (  4 * dr * dr  );
  origrepr.splcoeff5[5] = - v2nd[n-1] / ( 20 * pow(dr,3));
  
  ///// cubic splines y = a + bx + cx2 + dx3 /////
  double a = origrepr.splcoeff5[0];
  double b = origrepr.splcoeff5[1];
  double c,d;
  for ( int i=n-2; i>=0 ; i-- ) {
    dr = vr[i+1] - vr[i];
    d  = (v2nd[i+1] - v2nd[i]) / (6*dr)    ;
    c  =  v2nd[i  ]            /  2        ;
    b  =  b - d*3*dr*dr   - c*2*dr         ;
    a  =  a - d*pow(dr,3) - c*dr*dr - b*dr ;

    origrepr.splcoeff3[4*i  ] = a;
    origrepr.splcoeff3[4*i+1] = b;
    origrepr.splcoeff3[4*i+2] = c;
    origrepr.splcoeff3[4*i+3] = d;
  }

  ///// exponential fit on left margin /////
  origrepr.A = -2*c/b;
              origrepr.B = log(-b/origrepr.A) + origrepr.A*vr[0];
              origrepr.C = a-exp(-origrepr.A*vr[0]+origrepr.B);
  
  return origrepr;
}

void writeorigabSpline (const vector<double> &vr, const vector<double> &va, const double c1, const double cn, const int mod){
  // writes repulsiv potential to standardout in .abSpline-format (from original dftb 1998)
  //
  // vr    : knot vector inclusive cutoff
  // va[i] : value of curve to fit at vr[i]
  // c1    : second derivative devided by two at the first point of the curve
  // cn    : second derivative at vr[vr.size()-2] (one before the cutoff point)

  int icoeff=0;

  OrigRepr origrepr;
  if ( mod == 1 ) origrepr = calc_orig_mod1(vr,va,c1,cn);
  if ( mod == 2 ) origrepr = calc_orig_mod2(vr,va);

  
  // first three lines
  cout << "Spline" << endl;
  cout << origrepr.nspl << "  " << origrepr.cutoff << endl;
  cout << scientific << setprecision(15) << origrepr.A << "  " << origrepr.B << "  " << origrepr.C << endl;
  
  // 3rd order splines
  for (int i=0; i < vr.size()-2; i++ ){
    cout << setiosflags(ios::showpoint) << resetiosflags(ios::scientific);
          cout << setprecision(6) << setw(8) << vr[i] << "  " << setw(8) << vr[i+1] << "   ";
    for (int j=0; j <= 3; j++){
                cout << setprecision(15) << scientific << setw(23) << origrepr.splcoeff3[icoeff] << " ";
          icoeff++;
    }
    cout << endl;
  }

  // 5th order spline (last)
  cout << setiosflags(ios::showpoint) << resetiosflags(ios::scientific);
        cout << setprecision(6) << setw(8) << vr[vr.size()-2] << "  " << setw(8) << vr[vr.size()-1] << "   ";
  for (int i=0; i <= 5; i++){
    cout << setprecision(15) << scientific << setw(23) << origrepr.splcoeff5[i] << " "; 
  }
  cout << endl;
}
