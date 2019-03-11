#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <unistd.h>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "tools.hpp"
#include "ddh.hpp"
#include "erepobj.hpp"
#include "mpi.h"

extern int power;
extern sddh ddh;
extern string scratch;
extern int cpu_number;
using namespace std;
const double approximaterezo=0.000001; 
bool Erepobj::checkskf(const string ifilename, int type){
  ifstream ifile;
  ofstream ofile;
  string filename;
  bool passlc,passspline,passgrid,calculating,pass; 
  char cline[5120],ctemp1[512];
  int nline;
  double otmp;
  pass=passlc=passspline=passgrid=calculating=false;
  filename=libdir+"/"+ifilename;
  ifile.open(filename.c_str());
  if ( !ifile ) {
    ifile.close();
  }else{
    nline=0;
    string stemp;
    while(ifile.getline(cline,5120)){
      nline++;
      stemp=cline;
      if(stemp.find("Spline")!=string::npos){
        passspline=true; 
        if(nline>ngrid) passgrid=true;
      }else if(stemp.find("LC")!=string::npos){
        passlc=true;  
        sscanf(cline,"%s %le",ctemp1,&otmp);
        if(abs(otmp-omega)>approximaterezo) passlc=false;
      }
    }
    ifile.close();
  }

  if(lc){
    if(passspline && passgrid && passlc) pass=true;
  }else{
    if(passspline && passgrid ) pass=true;
  }

  if(type==1){
    if(!pass){
      cerr << endl << "ERROR: " << endl << endl;
      cerr << libdir<<"/"<<ifilename << ": is not correct" << endl << "exit score" << endl ;
      exit(1);
    }else return true;
  }else{
    if(!pass) return false;  
    else return true;
  }
} 
void Erepobj::writeskgen(const string& tmp_dir, const GAGenome& g) {
  GA1DArrayGenome<double>& gnome = (GA1DArrayGenome<double>&)g;
  int i,j,nvars,idx,idx2,idx0;
  string sorbital;
  string stemp = tmp_dir+"/"+"skdef.hsd";
  ofstream fout(stemp.c_str());
  bool hp,cp,np,op;
  hp=cp=np=op=false;

  idx0=0; 
  for(i=0;i<velem.size();i++) idx0=idx0+velem[i].lmax+2;
  if(ddh.td3) idx0=idx0+ddh.d3.size();
  if(ddh.tdamph) idx0++;
  if(ddh.thubbardderivs) idx0=idx0+ddh.hubbardderivs.size();

  fout.precision(6);
  fout<<"SkdefVersion = 1\n\
        \n\
        Globals {\n\
          Superposition = density\n";
  if(lc){
    fout<<"          func = xchybrid {omega = "<<omega<<"}\n";
  }else{
    fout<<"          func = pbe\n";
  }
  fout<<"        }\n\
        \n\
        AtomParameters {\n\
        \n\
          $OccupationsNe {\n\
            1S = 1.0 1.0\n\
            2S = 1.0 1.0\n\
            2P = 3.0 3.0\n\
          }\n\
          \n\
          $OccupationsAr {\n\
            $OccupationsNe\n\
            3S = 1.0 1.0\n\
            3P = 3.0 3.0\n\
          }\n\
          \n\
          $OccupationsKr {\n\
            $OccupationsAr\n\
            3D = 5.0 5.0\n\
            4S = 1.0 1.0\n\
            4P = 3.0 3.0\n\
          }\n\
          \n\
          $OccupationsXe {\n\
            $OccupationsKr\n\
            4D = 5.0 5.0\n\
            5S = 1.0 1.0\n\
            5P = 3.0 3.0\n\
          }\n\
          \n\
          $OccupationsHg {\n\
            $OccupationsXe\n\
            4F = 7.0 7.0\n\
            5D = 5.0 5.0\n\
            6S = 1.0 1.0\n\
          }\n\
          \n\
          $OccupationsRn {\n\
            $OccupationsHg\n\
            6P = 3.0 3.0\n\
          }\n\
          \n";
  
        idx=0;
        for(i=0; i<velem.size();i++){
          nvars=velem[i].lmax+2;
          if(velem[i].name=="H"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              H {\n\
                AtomConfig {\n\
                  AtomicNumber = 1\n\
                  Mass = 1.008\n";
          if(velem[i].lmax==0){ 
            fout<<"Occupations {\n\
                    1S = 1.0 0.0\n\
                  }\n\
                  ValenceShells = 1s\n";
          }else if(velem[i].lmax==1){ 
            hp=true;
            fout<<"Occupations {\n\
                    1S = 1.0 0.0\n\
                    2S = 0.0 0.0\n\
                    2P = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 1s 2p\n";
          }
            fout<<"Relativistics = None\n\
                }\n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="H" || ddh.hubbards[j].name=="h"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="H" || ddh.vorbes[j].name=="h"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="C"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              C {\n\
                AtomConfig {\n\
                  AtomicNumber = 6\n\
                  Mass = 12.01\n";
          if(velem[i].lmax==1){ 
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 1.0 1.0\n\
                  }\n\
                  ValenceShells = 2s 2p\n";
          }else if(velem[i].lmax==2){ 
            cp=true;
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 1.0 1.0\n\
                    3S = 0.0 0.0\n\
                    3P = 0.0 0.0\n\
                    3D = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 2s 2p 3d\n";
          }
            fout<<"Relativistics = None\n\
                }\n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="C" || ddh.hubbards[j].name=="c"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="C" || ddh.vorbes[j].name=="c"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="N"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              N {\n\
                AtomConfig {\n\
                  AtomicNumber = 7\n\
                  Mass = 14.007\n";
          if(velem[i].lmax==1){ 
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 2.0 1.0\n\
                  }\n\
                  ValenceShells = 2s 2p\n";
          }else if(velem[i].lmax==2){ 
            np=true;
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 2.0 1.0\n\
                    3S = 0.0 0.0\n\
                    3P = 0.0 0.0\n\
                    3D = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 2s 2p 3d\n";
          }
            fout<<"Relativistics = None\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="N" || ddh.hubbards[j].name=="n"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="N" || ddh.vorbes[j].name=="n"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="O"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              O {\n\
                AtomConfig {\n\
                  AtomicNumber = 8\n\
                  Mass = 16.01\n";
          if(velem[i].lmax==1){ 
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 2.0 2.0\n\
                  }\n\
                  ValenceShells = 2s 2p\n";
          }else if(velem[i].lmax==2){ 
            op=true;
            fout<<"Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    2P = 2.0 2.0\n\
                    3S = 0.0 0.0\n\
                    3P = 0.0 0.0\n\
                    3D = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 2s 2p 3d\n";
          }
            fout<<"Relativistics = None\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="O" || ddh.hubbards[j].name=="o"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="O" || ddh.vorbes[j].name=="o"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="P"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              P {\n\
                AtomConfig {\n\
                  AtomicNumber = 15\n\
                  Mass = 30.974\n\
                  Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    3S = 1.0 1.0\n\
                    2P = 3.0 3.0\n\
                    3P = 2.0 1.0\n\
                    3D = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 3s 3p 3d\n\
                  Relativistics = None\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    else if(j==4) sorbital="F";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="P" || ddh.hubbards[j].name=="p"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="P" || ddh.vorbes[j].name=="p"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="S"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              S {\n\
                AtomConfig {\n\
                  AtomicNumber = 16\n\
                  Mass = 32.065\n\
                  Occupations {\n\
                    1S = 1.0 1.0\n\
                    2S = 1.0 1.0\n\
                    3S = 1.0 1.0\n\
                    2P = 3.0 3.0\n\
                    3P = 2.0 2.0\n\
                    3D = 0.0 0.0\n\
                  }\n\
                  ValenceShells = 3s 3p 3d\n\
                  Relativistics = None\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    else if(j==4) sorbital="F";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="S" || ddh.hubbards[j].name=="s"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="S" || ddh.vorbes[j].name=="s"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          } else if(velem[i].name=="Br"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              Br {\n\
                AtomConfig {\n\
                  AtomicNumber = 35\n\
                  Mass = 79.904\n\
                  Occupations {\n\
                    $OccupationsAr\n\
                    3D = 5.0 5.0\n\
                    4S = 1.0 1.0\n\
                    4P = 3.0 2.0\n\
                  }\n\
                  ValenceShells = 4s 4p 4d\n\
                  Relativistics = Zora\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    else if(j==4) sorbital="F";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="Br" || ddh.hubbards[j].name=="br"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="Br" || ddh.vorbes[j].name=="br"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          }
          else if(velem[i].name=="I"){ 
            fout.precision(velem[i].radius[0].precision);
            fout<<"\n\
              I {\n\
                AtomConfig {\n\
                  AtomicNumber = 53\n\
                  Mass = 126.90447\n\
                  Occupations {\n\
                    $OccupationsKr\n\
                    4D = 5.0 5.0\n\
                    5S = 1.0 1.0\n\
                    5P = 3.0 2.0\n\
                  }\n\
                  ValenceShells = 5s 5p 5d\n\
                  Relativistics = Zora\n\
                }\n\
                \n\
                DftbAtom {\n\
                  ShellResolved = No\n\
                  DensityCompression = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n\
                  WaveCompressions = SingleAtomCompressions {\n";
                  idx++;
                  for(j=1; j<nvars; j++){
                    fout.precision(velem[i].radius[j].precision);
                    if(j==1) sorbital="S";
                    else if(j==2) sorbital="P";
                    else if(j==3) sorbital="D";
                    else if(j==4) sorbital="F";
                    fout<<sorbital<<" = PowerCompression { Power = "<<power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
                    idx++;
                  }
            fout<<"}\n";
            idx2=idx0;
            if(ddh.thubbards){
              fout<<"CustomizedHubbards {\n";
              for(j=0;j<ddh.hubbards.size();j++){ 
                if(ddh.hubbards[j].name=="I" || ddh.hubbards[j].name=="i"){
                  fout.precision(ddh.hubbards[j].precision);
                  fout<<ddh.hubbards[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            if(ddh.tvorbes){
              fout<<"CustomizedOnsites {\n";
              for(j=0;j<ddh.vorbes.size();j++){ 
                if(ddh.vorbes[j].name=="I" || ddh.vorbes[j].name=="i"){
                  fout.precision(ddh.vorbes[j].precision);
                  fout<<ddh.vorbes[j].type<<" = "<< gnome.gene(idx2)<<endl;
                }
                idx2++;
              }
              fout<<"}\n";
            }
            fout<<"}\n\
              }\n";
          }
        }

       
  fout.precision(6);
  fout<<"}\n\
        \n\
        \n\
        OnecenterParameters {\n\
          \n\
          $StandardDeltaFilling {\n\
            DeltaFilling = 0.001\n\
          }\n\
          \n\
          H {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n";
if(!hp){
  fout<<"     Exponents {\n\
                S = 0.50 1.0 2.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n";
}else{
  fout<<"     Exponents {\n\
                S = 0.50 1.0 2.0\n\
                P = 0.50 1.0 2.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n";
}
 fout<<"      }\n\
            }\n\
          }\n\
          \n\
          C {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n";
if(!cp){
  fout<<"     Exponents {\n\
                S = 0.5 1.14 2.62 6.0\n\
                P = 0.5 1.14 2.62 6.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n";
}else{        
  fout<<"     Exponents {\n\
                 S = 0.5 1.14 2.62 6.0\n\
                 P = 0.5 1.14 2.62 6.0\n\
                 D = 0.5 1.14 2.62 6.0\n\
               }\n\
               MaxPowers {\n\
                 S = 3\n\
                 P = 3\n\
                 D = 3\n";
}         
  fout<<"     }\n\
            }\n\
          }\n\
          \n\
          N {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n";
if(!np){
  fout<<"     Exponents {\n\
                S = 0.5 1.21 2.9 7.0\n\
                P = 0.5 1.21 2.9 7.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n";
}else{   
  fout<<"     Exponents {\n\
                S = 0.5 1.21 2.9 7.0\n\
                P = 0.5 1.21 2.9 7.0\n\
                D = 0.5 1.21 2.9 7.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n";
}        
  fout<<"     }\n\
            }\n\
          }\n\
          \n\
          O {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n";
if(!op){
  fout<<"     Exponents {\n\
                S = 0.5 1.26 3.17 8.0\n\
                P = 0.5 1.26 3.17 8.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n";
}else{   
  fout<<"     Exponents {\n\
                S = 0.5 1.26 3.17 8.0\n\
                P = 0.5 1.26 3.17 8.0\n\
                D = 0.5 1.26 3.17 8.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n";
}        
  fout<<"     }\n\
            }\n\
          }\n\
          \n\
          P {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n\
              Exponents {\n\
                S = 0.5 1.17 2.74 6.41 15.0\n\
                P = 0.5 1.17 2.74 6.41 15.0\n\
                D = 0.5 1.17 2.74 6.41 15.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n\
              }\n\
            }\n\
          }\n\
          \n\
          S {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n\
              Exponents {\n\
                S = 0.5 1.19 2.83 6.73 16.0\n\
                P = 0.5 1.19 2.83 6.73 16.0\n\
                D = 0.5 1.19 2.83 6.73 16.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n\
              }\n\
            }\n\
          }\n\
          \n\
          Br {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n\
              Exponents = {\n\
                S = 0.5 1.17 2.74 6.4 35.0\n\
                P = 0.5 1.17 2.74 6.4 35.0\n\
                D = 0.5 1.17 2.74 6.4 35.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n\
              }\n\
            }\n\
          }\n\
          \n\
          I {\n\
            $StandardDeltaFilling\n\
            Calculator = SlaterAtom {\n\
              Exponents = {\n\
                S = 0.5 1.27 3.23 8.21 53.0\n\
                P = 0.5 1.27 3.23 8.21 53.0\n\
                D = 0.5 1.27 3.23 8.21 53.0\n\
              }\n\
              MaxPowers {\n\
                S = 3\n\
                P = 3\n\
                D = 3\n\
              }\n\
            }\n\
          }\n\
          \n\
        }\n\
        \n\
        \n\
        TwoCenterParameters {\n\
        \n\
          $EqGrid = EquidistantGrid {\n\
              GridStart = 0.4\n\
              GridSeparation = "<<dgrid<<"\n\
              Tolerance = 5e-9\n\
              MaxDistance = "<<ngrid<<"\n\
          }\n\
          \n\
          $SkTwocnt_300_150 = Sktwocnt {\n\
            IntegrationPoints = 200 50\n\
          }\n\
          \n\
          $SkTwocnt_400_200 = Sktwocnt {\n\
            IntegrationPoints = 400 200\n\
          }\n\
          \n\
          H-H { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          H-C { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          H-N { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          H-O { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          H-P { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          H-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          H-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          H-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          C-C { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          C-N { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          C-O { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          C-P { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          C-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          C-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          C-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          N-N { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          N-O { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          N-P { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          N-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          N-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          N-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          O-O { Grid = $EqGrid; Calculator = $SkTwocnt_300_150 }\n\
          O-P { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          O-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          O-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          O-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          P-P { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          P-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          P-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          P-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          S-S { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          S-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          S-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 } \n\
          Br-Br{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          Br-I{ Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          I-I { Grid = $EqGrid; Calculator = $SkTwocnt_400_200 }\n\
          \n\
        }\n\n";
} 

void Erepobj::makeskf(const GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int ie1,ie2,i,j,k,it,idx;
  bool pass,pass2;
  ostringstream ss;
  int id=0;
  double tmp,tmps,tmpp;
  string sfilenametmp;

//ie1=-1;
//for(i=0;i<velem.size();i++){
//  ie1=ie1+2;  
//  tmp=genome.gene(ie1);
//  if(velem[i].optype==2){
//    genome.gene(ie1+1,tmp);
//  }else if(velem[i].optype==3){
//    genome.gene(ie1+1,tmp);
//    genome.gene(ie1+2,tmp);
//  }
//  ie1=ie1+velem[i].lmax; 
//}
  ie1=-1;
  for(i=0;i<velem.size();i++){
    ie1=ie1+2;  
    tmps=genome.gene(ie1);
    tmpp=genome.gene(ie1+1);
    if(velem[i].optype==2){
      genome.gene(ie1+1,tmps);
    }else if(velem[i].optype==3){
      genome.gene(ie1+1,tmps);
      genome.gene(ie1+2,tmps);
    }else if(velem[i].optype==12){
      genome.gene(ie1,GAMin(tmps, tmpp));
      genome.gene(ie1+1,GAMax(tmps, tmpp));
    }else if(velem[i].optype==14){
      genome.gene(ie1,GAMin(tmps, tmpp));
      genome.gene(ie1+1,GAMax(tmps, tmpp));
      if(genome.gene(ie1+1)-genome.gene(ie1)<0.4) genome.gene(ie1,genome.gene(ie1+1)-0.4); 
    }else if(velem[i].optype==21){
      genome.gene(ie1,GAMax(tmps, tmpp));
      genome.gene(ie1+1,GAMin(tmps, tmpp));
    }else if(velem[i].optype==11){
      genome.gene(ie1-1,tmps);
    }else if(velem[i].optype==111){
      genome.gene(ie1-1,tmps);
      genome.gene(ie1+1,tmps);
    }
    ie1=ie1+velem[i].lmax; 
  }
 
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  string call,stemp,stemp2,part1,part2,part3; 
  string skf_dir=scratch+"/buildsks.tmp_"+itoa(id,10);
  call="rm -rf "+skf_dir+"; mkdir -p "+skf_dir+";";
	it=system(call.c_str());	
  writeskgen(skf_dir,genome);

  idx=0; 
  for(i=0;i<velem.size();i++) idx=idx+velem[i].lmax+2;
  if(ddh.td3) idx=idx+ddh.d3.size();
  if(ddh.tdamph) idx++;
  if(ddh.thubbardderivs) idx=idx+ddh.hubbardderivs.size();

  part3="";
  if(ddh.thubbards){
    part3+="_hub";
    for(i=0;i<ddh.hubbards.size();i++){ 
      ss.str(std::string());
      ss.precision(ddh.hubbards[i].precision); 
      ss<<std::fixed<<genome.gene(idx);
      part3+="_"+ddh.hubbards[i].name+ss.str();
      idx++;
    }
  }
  if(ddh.tvorbes){
    part3+="_vor";
    for(i=0;i<ddh.vorbes.size();i++){ 
      ss.str(std::string());
      ss.precision(ddh.vorbes[i].precision); 
      ss<<std::fixed<<genome.gene(idx);
      part3+="_"+ddh.vorbes[i].name+ss.str();
      idx++;
    }
  }

  ie1=0;
  for(i=0;i<velem.size();i++){
    part1=velem[i].name;  
    ss.str(std::string());
    ss.precision(velem[i].radius[0].precision); ss<<std::fixed<<genome.gene(ie1);
    part1+="_d"+ss.str();
    ie1++;  
    for(k=1;k<velem[i].lmax+2.;k++){
      ss.str(std::string());
      ss.precision(velem[i].radius[k].precision); ss<<std::fixed<<genome.gene(ie1);
      part1+="_w"+ss.str();  
      ie1++; 
    }
    ie2=0;
    for(j=0;j<velem.size();j++){
      part2=velem[j].name;  
      ss.str(std::string());
      ss.precision(velem[j].radius[0].precision); ss<<std::fixed<<genome.gene(ie2);
      part2+="_d"+ss.str();  
      ie2++;
      for(k=1;k<velem[j].lmax+2.;k++){
        ss.str(std::string());
        ss.precision(velem[j].radius[k].precision); ss<<std::fixed<<genome.gene(ie2);
        part2+="_w"+ss.str();  
        ie2++;
      }
      if(lc){
        ss.str(std::string());
        ss.precision(2);ss<<std::fixed<<omega;
        stemp="omega_"+ss.str()+"_"+part1+"-"+part2+part3+".skf";
        stemp2="omega_"+ss.str()+"_"+part2+"-"+part1+part3+".skf";
      }else{
        stemp=part1+"-"+part2+part3+".skf";
        stemp2=part2+"-"+part1+part3+".skf";
      }
      if( i!=j){
//      #pragma omp critical (checkskf02)
//      if(skfexist[i][j]){
        if(addskf) sfilenametmp=libadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
        if(addskf && !access(sfilenametmp.c_str(), F_OK )){
          pass=true;
          pass2=true;
        }else{
          pass=checkskf(stemp,0);
          pass2=checkskf(stemp2,0);
        }
        if((!pass)&&(!pass2)){
          call=" cd "+skf_dir+"; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";";
          call+= skgen +" -l error -o "+onecent+" -t "+twocent+" sktable "+velem[i].name+" "+velem[j].name+" >log;";
          if(lc){
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                  sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/RangeSep/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  sed -i '/LC/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  cat "+ repadddir+"/"+velem[j].name+"-"+velem[i].name+".skf >> "+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"RangeSep\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"$omega\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }else{
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                  sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/RangeSep/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  sed -i '/LC/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"Spline\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"1 9.99\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"RangeSep\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"$omega\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }
          }else{
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                  cat "+ repadddir+"/"+velem[j].name+"-"+velem[i].name+".skf >> "+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }else{
              call+="echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"Spline\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"1 9.99\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }
          }
          it=system(call.c_str());	
        }
      }else{
 //     #pragma omp critical (checkskf01)
//      if(skfexist[i][j]){
        if(addskf) sfilenametmp=libadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
        if(addskf && !access(sfilenametmp.c_str(), F_OK )){
          pass=true;
        }else{
          pass=checkskf(stemp,0);
        }
        if(!pass){
          call=" cd "+skf_dir+"; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";";
          call+= skgen +" -l error -o "+onecent+" -t "+twocent+" sktable "+velem[i].name+" "+velem[j].name+" >log;";
          if(lc){
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                   sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }else{
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                   sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }
          }else{
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }else{
              call+="echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }
          }
          //cout<<"####################\n"<<call<<endl;
          it=system(call.c_str());	
        }
      }
    }
  }

  call="rm -rf "+skf_dir+";";
	it=system(call.c_str());	

} 


          //IntegrationPoints = 200 50\n\
          //IntegrationPoints = 300 150\n\
            
            
