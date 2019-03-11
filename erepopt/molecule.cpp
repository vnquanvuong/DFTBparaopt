#include <iostream>     // std::cout
#include <iterator>     // std::ostream_iterator
#include <vector>       // std::vector
#include <algorithm>    // std::copy
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include <iterator> 
#include <map>
#include "tools.hpp"
#include "molecule.hpp"
#include "erepobj.hpp"

using namespace std;

extern vector<string> addHamiltonian;
extern Erepobj erepobj;
extern string scratch;

void Molecule::writegen(const string& tmp_dir) {
  string stemp = tmp_dir+"/"+"in.gen";
  ofstream fout(stemp.c_str());
  fout << natom << " C" << endl;
  copy(elemnameext.begin(),elemnameext.end(),ostream_iterator<string>(fout, " "));
  fout << endl;
  for(int i=0;i<natom;++i) {
    fout << i+1 << " " << atomindex[i]+1 << " " 
         << coord(i,0)*AA_Bohr << " " << coord(i,1)*AA_Bohr << " " <<  coord(i,2)*AA_Bohr << endl;
    //cout << i+1 << " " << atomindex[i]+1 << " " 
    //                 << coord(i,0) << " " << coord(i,1) << " " <<  coord(i,2) << endl;
  }
  fout.close();
}

void Molecule::writedftb(const string& tmp_dir, const string& skf_dir) {

  string stemp = tmp_dir+"/"+"dftb_in.hsd";
  ofstream fout(stemp.c_str());
  fout<<"Geometry = GenFormat {\n\
          <<< \"in.gen\"\n\
        }\n\
        Hamiltonian = DFTB {\n\
          Eigensolver = RelativelyRobust{}\n\
          Charge = "<<icharge<<"\n\
          SCC = Yes\n\
          SlaterKosterFiles = Type2Filenames{\n\
             Prefix = \"../"<<skf_dir<<"/\"\n\
             Suffix = \".skf\"\n\
             Separator = \"-\"\n\
             LowerCaseTypeName = No\n\
          }\n";
  for(int i=0; i<addHamiltonian.size();i++){
    fout<<addHamiltonian[i]<<endl;
  }
  fout<<"}\n\
        ParserOptions {\n\
          IgnoreUnprocessedNodes = Yes\n\
        }\n";

  fout.close();
}

void Molecule::init(const string& xyzfilename, const int icharge_, const double start_, const double end_, const double gauss_, const string& dosfilename, const double weight_) {

  int i,j;
  double x,y,z;
  char cline[512],ctemp[512];
  double d1,d2;

  string dummy;
  name    = xyzfilename;
  icharge = icharge_;
  weight  = weight_;
  start   = start_;
  end     = end_;  
  gauss   = gauss_;  

  nelem=natom=0; 
  nelectron=-icharge;

  map<string,int> index;

  ifstream fin(xyzfilename.c_str());
  fin >> natom;
  
  coord.resize(natom,3);
  atomname.resize(natom);
  atomindex.resize(natom);

  getline(fin,dummy);
  getline(fin,dummy);
  
  for (i=0; i<natom ; i++) {
    fin >> dummy >> x >> y >> z ;
    
    coord(i,0) = x/AA_Bohr ;
    coord(i,1) = y/AA_Bohr ;
    coord(i,2) = z/AA_Bohr ;

    //to_lower(dummy);
    for(int it=0; it<dummy.size(); it++) {
      dummy[it] = tolower(dummy[it]);
    }
    if ( dummy == "h" ) nelectron+=1;
    else if ( dummy == "c" ) nelectron+=4;
    else if ( dummy == "n" ) nelectron+=5;
    else if ( dummy == "o" ) nelectron+=6;
    else if ( dummy == "p" ) nelectron+=5;
    else if ( dummy == "s" ) nelectron+=6;

    atomname[i] = dummy;
    index.insert(make_pair(dummy,0));
  }

  dosrefs.resize(0);
  dosrefd.resize(0);
  evref.resize(0);
  occref.resize(0);

  nelem = index.size() ;

  i=0;
  for(map<string,int>::iterator p=index.begin();p!=index.end();++p) {
    (*p).second = i;
    dummy = (*p).first;
    elemname.push_back(dummy);
    dummy[0] = toupper(dummy[0]);
    if (dummy.substr(1,1)=="_") dummy.erase(1,1);
    elemnameext.push_back(dummy);
    ++i;
  }

  for (i=0; i<natom ; i++) {
    atomindex[i] = index[atomname[i]];
  }
  
  fin.close();

  // read reference
  fin.open(dosfilename.c_str());
  if(erepobj.fit_type==1){
    while(fin.getline(cline,512)){
      if(sscanf(cline,"%s",ctemp)<0) continue;
      if(ctemp[0]=='#' ) continue;
      if(ctemp[0]=='K' ) continue;
      sscanf(cline,"%le%le",&d1,&d2);
      evref.push_back(d1);
      occref.push_back(d2);
    }
    b0=0.0;
    for (i=0; i<natom ; i++) {
      for (j=0; j<erepobj.atomicname.size(); j++) {
        if(atomname[i] == erepobj.atomicname[j]);
//      cout<<atomname[i]<<erepobj.atomicname[j]<<erepobj.atomicenergy[j]<<endl;
          b0=b0+erepobj.atomicenergy[j];
          break;
      }
    }
    b0=b0/double(nelectron);
//  weight=2.0*weight/double(nelectron);
  }else if(erepobj.fit_type==2){ 
    while(fin.getline(cline,512)){
      if(sscanf(cline,"%s",ctemp)<0) continue;
      if(ctemp[0]=='#' ) continue;
      if(ctemp[0]=='K' ) continue;
      sscanf(cline,"%le%le",&d1,&d2);
      dosrefs.push_back(d1);
      dosrefd.push_back(d2);
    }
    sgrid=dosrefs[2]-dosrefs[1];
  }
  fin.close();
}

double Molecule::score(int id, int score_type, const string& dftbversion, bool getresidual){
  bool pass;
  char cline[512],ctemp[512];
  int i,iq,j,iref,inew;
  double score,d1,d2,dtmp,dmin,tweight;
  double reftmp,newtmp,delta; 
  string stemp,dftb_dir,skf_dir;
  
  vector<double> doss;
  vector<double> dosd;

  doss.resize(0);
  dosd.resize(0);

  dftb_dir=scratch+"/dftb.tmp_"+itoa(id,10);
  skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);
  string call  = "rm -rf "+dftb_dir+"/; mkdir -p "+dftb_dir+";";
  iq = system(call.c_str());	

  writegen(dftb_dir);
  writedftb(dftb_dir,skf_dir);

  ostringstream ss,ss2;
  ss.str(std::string());
  ss2.str(std::string());
  ss.precision(6); ss<<std::fixed<<gauss;
  ss2.precision(6); ss2<<std::fixed<<sgrid;
  call = "cd "+dftb_dir+"/; " + dftbversion + " > log; dp_dos -b "+ss.str()+" -g "+ss2.str()+" band.out dos.dat; " + "cd ../ ";
  iq = system(call.c_str());
  
  stemp=dftb_dir+"/dos.dat";
  ifstream fin(stemp.c_str());
  if ( !fin ) {
    cerr << endl << "ERROR: " << endl << endl;
    cerr << "The system call: \"" << call << "\" did not produce the file: \"dftb.tmp/dos.dat\""
         << endl << "calculation of electronic energy for molecule \"" << name << "\" failed."
         << endl << "exit erepobj" << endl << endl;
//  exit(1);
    return 999999999.9;
  }
  pass=false;
  while(fin.getline(cline,512)){
    stemp=cline;
    if(sscanf(cline,"%s",ctemp)<0) continue;
    if(ctemp[0]=='#' ) continue;
    pass=true;
    sscanf(cline,"%le%le",&d1,&d2);
    doss.push_back(d1);
    dosd.push_back(d2);
  }
  fin.close();

  
  call  = "rm -rf "+dftb_dir+"/;";
  iq = system(call.c_str());	


//for(iref=0; iref<dosrefs.size();iref++){
//  if(abs(dosrefs[iref])<0.000001) break;
//}

  dmin=999999999.9;
  for(i=0; i<dosrefs.size();i++){
    dtmp=abs(dosrefs[i]-end);
    if(dtmp<dmin){
      dmin=dtmp; 
      iref=i;
    }
  }

//for(inew=0; inew<doss.size();inew++){
//  if(abs(doss[inew])<0.000001) break;
//}

  dmin=999999999.9;
  for(i=0; i<doss.size();i++){
    dtmp=abs(doss[i]-end);
    if(dtmp<dmin){
      dmin=dtmp; 
      inew=i;
    }
  }

//cout<<iref<<"  "<<inew<<endl;
  score=0.0;
  i=iref; 
  j=inew;
//delta=dosrefs[i]-dosrefs[i-1];
//dtmp=doss[j]-doss[j-1];
//if (abs(dtmp-delta)>0.000001) {
//  cerr << endl << "ERROR: " << endl << endl;
//  cerr << "grid is not consistent\""
//       << endl << "calculation of electronic energy for molecule \"" << name << "\" failed."
//       << endl << "exit erepobj" << endl << endl;
//  exit(1);

//}
//dtmp=dosrefs[i];
//while((dtmp>start) && (dtmp<end) ){
//  if(i>(dosrefs.size()-1))reftmp=0.0;
//  else reftmp=dosrefd[i];
//  
//  if(j>(doss.size()-1))newtmp=0.0;
//  else newtmp=dosd[j];
//  
//  if(score_type==1) score+=abs(newtmp-reftmp);
//  else if(score_type==2) score+=pow((newtmp-reftmp),2.0);
//  else if(score_type==4) score+=pow((newtmp-reftmp),4.0);
//  i++;
//  j++;
//  dtmp=dtmp+delta;
//}
  i=iref-1; 
  j=inew-1;
  dtmp=dosrefs[i];
  while((dtmp>start) && (dtmp<end) ){
    if(i<0)reftmp=0.0;
    else reftmp=dosrefd[i];
    
    if(j<0)newtmp=0.0;
    else newtmp=dosd[j];
   
    if(abs(weight+1.0)<0.000001) tweight=(dtmp-start)/(end-start);
    else if(abs(weight+2.0)<0.000001) tweight=pow((dtmp-start)/(end-start),2.0);
    else if(abs(weight+4.0)<0.000001) tweight=pow((dtmp-start)/(end-start),4.0);
    else tweight=weight;

    if(score_type==1) score+=tweight*abs(newtmp-reftmp);
    else if(score_type==2) score+=tweight*pow((newtmp-reftmp),2.0);
    else if(score_type==4) score+=tweight*pow((newtmp-reftmp),4.0);
    i--;
    j--;
    dtmp=dtmp-sgrid;
  }

  if(getresidual){
    error=2.0*score/double(nelectron);
  }

  return 2.0*score/double(nelectron);
}

//double Molecule::getev(int id, , const string& dftbversion, bool getresidual){
void Molecule::getev(int id, int &pos, double * evarray, const string& dftbversion){

  bool pass;
  char cline[512],ctemp[512];
  int i,iq,j,iref,inew;
  double score,d1,d2,dtmp,dmin,tweight;
  double reftmp,newtmp,delta; 
  string stemp,dftb_dir,skf_dir;

  dftb_dir=scratch+"/dftb.tmp_"+itoa(id,10);
  skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);
  string call  = "rm -rf "+dftb_dir+"/; mkdir -p "+dftb_dir+";";
  iq = system(call.c_str());	

  writegen(dftb_dir);
  writedftb(dftb_dir,skf_dir);

  call = "cd "+dftb_dir+"/; " + dftbversion + " > log; cd ../ ";
  iq = system(call.c_str());
  
  stemp=dftb_dir+"/band.out";
  ifstream fin(stemp.c_str());
  if ( !fin ) {
//  cerr << endl << "ERROR: " << endl << endl;
//  cerr << "The system call: \"" << call << "\" did not produce the file: \"dftb.tmp/band.out\""
//       << endl << "calculation of electronic energy for molecule \"" << name << "\" failed."
//       << endl << "exit erepobj" << endl << endl;
    for(i=0; i<nelectron/2; i++){
      evarray[pos]=999999999.9;
      pos++;
    }
  }else{
    fin.getline(cline,512);
    for(i=0; i<nelectron/2; i++){
      fin.getline(cline,512);
      sscanf(cline,"%le%le",&d1,&d2);
      evarray[pos]=d1;
      pos++;
    }
    fin.close();
  }
  call  = "rm -rf "+dftb_dir+"/;";
  iq = system(call.c_str());	
}

