#include <iostream>     // std::cout
#include <iterator>     // std::ostream_iterator
#include <vector>       // std::vector
#include <algorithm>    // std::copy
#include <fstream>
#include <string>
#include "auxiliary.hpp"
#include "allequations.hpp"

using namespace std;
extern bool runtest;
extern sddh ddh;
extern void remove_comment(char *cstr);
extern int dftbout_new;
void Molecule::init(const string& xyzfilename, const double ediss_, const double eweight_, const double fweight_,
                    const string& dftbin_,     const string abbr_ , const string& dftbversion, const string& finp_) {
  char ctline[512];
  double dtmp;
  string stmp;
  fstream infile,outfile;
  name    = xyzfilename;
  ebind   = -ediss_/kcal_H;
  eweight = eweight_;
  fweight = fweight_;
  dftbin  = dftbin_;
  abbr    = abbr_;
  finp    = finp_;
  nelem   = 0;
  natom   = 0;

  map<string,int> index;
  index.clear();

  ifstream fin(xyzfilename.c_str());
  string dummy;

  fin >> natom;
  
  coord.resize(natom,3);
  atomname.resize(natom);
  atomindex.resize(natom);

  getline(fin,dummy);
  getline(fin,dummy);
  
  for (int iat=0; iat<natom ; iat++) {
    double x, y, z;
    fin >> dummy >> x >> y >> z ;
    
    coord(iat,0) = x/AA_Bohr ;
    coord(iat,1) = y/AA_Bohr ;
    coord(iat,2) = z/AA_Bohr ;

    for(int i=0; i<dummy.size(); i++) {
      dummy[i] = tolower(dummy[i]);
    }
    if ( dummy.length() == 1 )  dummy += "_";
    toupper(dummy[0]);   // or is there really a need for that??
    atomname[iat] = dummy;
    index.insert(make_pair(dummy,0));
  }

  nelem = index.size() ;

  int i=0;
  for(map<string,int>::iterator p=index.begin();p!=index.end();++p) {
    (*p).second = i;
    dummy = (*p).first;
    elemname.push_back(dummy);
    dummy[0] = toupper(dummy[0]);
    if (dummy.substr(1,1)=="_") dummy.erase(1,1);
    elemnameext.push_back(dummy);
    ++i;
  }

  for (int iat=0; iat<natom ; iat++) {
    atomindex[iat] = index[atomname[iat]];
  }

  fin.close();

  // size force matrices and initialize to 0
  fel.resize(natom,3);
  fref.resize(natom,3);
  frefmel.resize(natom,3);
  for (int i=0;i<natom;i++) {
    for (int j=0;j<3;j++) {
      fel(i,j)=0;
      fref(i,j)=0;
      frefmel(i,j)=0;
    }
  }

  // read reference force fref
  if (finp != "0") {
    fin.open(finp.c_str());
    for (int iat=0; iat<natom ; iat++) {
      fin >> fref(iat,0) >> fref(iat,1) >> fref(iat,2);
    }
    fin.close();
  }

  //calculate eel and ebind - eel = ebindmel and get forces fel
  //string stilltochange  = "rm -rf dftb.tmp/; mkdir -p dftb.tmp; auto_xyz2gen " + name + " > dftb.tmp/in.gen";

  string stilltochange  = "rm -rf dftb.tmp/; mkdir -p dftb.tmp;";
  int iq = system(stilltochange.c_str()); 

  writegen("dftb.tmp/in.gen");

  //string call ="cd dftb.tmp/; ";
  //string call = "cp " + dftbin + " dftb.tmp/dftb_in.hsd; " + "cd dftb.tmp/; ";
  //if ( dftbin.substr(0,5) == "/home" )          // still not all possible mistakes considered... 
  //{
  //  call += dftbversion + " < " + dftbin + " > /dev/null; " ;//+ name + ".out; "; 
  //}else{
  //  call += dftbversion + " < ../" + dftbin + " > /dev/null; "; //../ + name + ".out; "; 
  //}
 
  infile.open(dftbin.c_str(),ios::in);
  if(!infile.is_open()){
    cout<<"unable to open "<<dftbin<<" file.\n";
    exit(1);
  }
  outfile.open("dftb.tmp/dftb_in.hsd",ios::out);
  if(!outfile.is_open()){
    cout<<"unable to open dftb.tmp/dftb_in.hsd file.\n";
    exit(1);
  }
  outfile<<std::fixed;
  while(infile.getline(ctline,512)){
    outfile<<ctline<<endl;
    remove_comment(ctline);
    stmp=ctline;
    if(stmp.find("Hamiltonian")!=string::npos){
      if(ddh.td3){
        outfile<<"  Dispersion = DftD3{\n"; 
        outfile<<"    Damping = BeckeJohnson{\n"; 
        for(i=0;i<ddh.d3.size();i++){ 
          if(ddh.d3[i].name=="a1" ||ddh.d3[i].name=="a2"){
            outfile.precision(ddh.d3[i].precision);
            outfile<<"      "<<ddh.d3[i].name<<" = "<<ddh.d3[i].value<<endl;
          }
        }
        outfile<<"    "<<"}\n"; 
        for(i=0;i<ddh.d3.size();i++){ 
          if(ddh.d3[i].name=="s6" ||ddh.d3[i].name=="s8"){
            outfile.precision(ddh.d3[i].precision);
            outfile<<"    "<<ddh.d3[i].name<<" = "<<ddh.d3[i].value<<endl;
          }
        }
        outfile<<"  "<<"}\n"; 
      }
      if(ddh.tdamph){
        outfile<<"  DampXH = Yes\n"; 
        outfile.precision(ddh.damph.precision);
        outfile<<"  DampXHExponent = "<<ddh.damph.value<<endl; 
      }
      if(ddh.thubbardderivs){
        if(ddh.thirdorderfull) outfile<<"  ThirdOrderFull = Yes\n"; 
        else  outfile<<"  ThirdOrder = Yes\n"; 
        outfile<<"  HubbardDerivs{\n"; 
        for(i=0;i<ddh.hubbardderivs.size();i++){ 
          outfile.precision(ddh.hubbardderivs[i].precision);
          outfile<<"    "<<ddh.hubbardderivs[i].name<<" = "<<ddh.hubbardderivs[i].value<<endl;
        }
        outfile<<"  "<<"}\n"; 
      }
      if(ddh.tdamphver2){
        outfile<<"  HCorrection = DampingVer2{\n"; 
        for(i=0;i<ddh.damphver2.size();i++){ 
          outfile.precision(ddh.damphver2[i].precision);
          outfile<<"    "<<ddh.damphver2[i].name<<" = "<<ddh.damphver2[i].value<<endl;
        }
        outfile<<"  "<<"}\n"; 
      }
    }
  }
  infile.close();
  outfile.close();

  string call = "cd dftb.tmp/; "+ dftbversion + " > log; cd ../ ";
  //cout << call << endl;
  iq = system(call.c_str());
        
  fin.open("dftb.tmp/detailed.out");
  if ( !fin ) {
    cerr << endl << "ERROR: " << endl << endl;
    cerr << "The system call: \"" << call << "\" did not produce the file: \"dftb.tmp/detailed.out\""
         << endl << "calculation of electronic energy for molecule \"" << name << "\" failed."
         << endl << "exit repopt" << endl << endl;
    exit(1);
  }
  string stemp;
  bool passscc,passe,passf;
  char cline[512],ctemp1[512],ctemp2[512],ctemp3[512];
  ctemp1[0]=ctemp2[0]=ctemp3[0]='\0';
  passscc=true;
  passe=passf=false;
  while(fin.getline(cline,512)){
    stemp=cline;
    if(!runtest){
      //if(stemp.find("Total Electronic energy:")!=string::npos){
      //  sscanf(cline,"%s%s%s%le",ctemp1,ctemp2,ctemp3,&eel);
      if(stemp.find("Total Mermin free energy:")!=string::npos){
        sscanf(cline,"%s%s%s%s%le",ctemp1,ctemp2,ctemp1,ctemp2,&eel);
        passe=true;
      }else if(stemp.find("Total Forces")!=string::npos){
        passf=true;
        for (int i=0; i<natom ; i++) {
          if(dftbout_new) fin >>  dtmp;
          fin >> fel(i,0) >> fel(i,1) >> fel(i,2);
        }
      }else if(stemp.find("SCC is NOT converged")!=string::npos){
        passscc=false;
      }
    }else{
      if(stemp.find("Total Mermin free energy:")!=string::npos){
        sscanf(cline,"%s%s%s%s%le",ctemp1,ctemp2,ctemp1,ctemp2,&eel);
        passe=true;
      }else if(stemp.find("Total Forces")!=string::npos){
        passf=true;
        for (int i=0; i<natom ; i++) {
          if(dftbout_new) fin >>  dtmp;
          fin >> fel(i,0) >> fel(i,1) >> fel(i,2);
        }
      }else if(stemp.find("SCC is NOT converged")!=string::npos){
        passscc=false;
      }
    }
  }
  fin.close();
  
  if ( eel == 0 ) {
    //cerr << endl << "ERROR: " << name << " could not be calculated by " << dftbversion << endl 
    //     << "exit repopt"    << endl << endl; exit(0);
    passscc=false;
  }
  if (!passe) {
    //cerr << endl << "ERROR: Energy of " << name << " could not be calculated by " << dftbversion << endl 
    //     << "exit repopt"    << endl << endl; exit(0);
    passscc=false;
  }
  if (!passf) {
    //cerr << endl << "ERROR: Force of " << name << " could not be calculated by " << dftbversion << endl 
    //     << "exit repopt"    << endl << endl; exit(0);
    passscc=false;
  }
  if(!passscc){
    eel=999999.9;
    for (int i=0; i<natom; i++) {
      fel(i,0)=999999.9;
      fel(i,1)=999999.9; 
      fel(i,2)=999999.9;
    }
  }
  
  //iq = system("rm -rf dftb.tmp/");

  ebindmel = ebind - eel;
  for (int i=0; i<natom ; i++) {
    for (int j=0;j<3;j++) {
      frefmel(i,j)=fref(i,j)-fel(i,j);
      //cerr << i << "  " << j << "    " << fref(i,j) << "   "  << frefmel(i,j) << endl;
    }
  }

  // setup distance matrix (upper triangular)
  dist.resize(natom,natom);
  for (int i=0; i<natom; i++) {
    for (int j=0; j<natom; j++) {
      dist(i,j) = 0;
    }
  }
  double hlp;
  for (int at1=0; at1 < natom-1; at1++) {
    for (int at2=at1+1; at2 < natom; at2++) {
      hlp = 0;
      for (int i=0; i<3; i++) {
        hlp += ( coord(at1,i)-coord(at2,i) ) * ( coord(at1,i)-coord(at2,i) );
      }
      dist(at1,at2) = sqrt(hlp);
    }
  }
}

void Molecule::writegen(const string& filename){
  ofstream fout(filename.c_str());
  fout << natom << " C" << endl;
  copy(elemnameext.begin(),elemnameext.end(),ostream_iterator<string>(fout, " "));
  fout << endl;
  for(int i=0;i<natom;++i) {
    fout << i+1 << " " << atomindex[i]+1 << " " 
         << coord(i,0)*AA_Bohr << " " << coord(i,1)*AA_Bohr << " " <<  coord(i,2)*AA_Bohr << endl;
  }
  fout.close();
}

