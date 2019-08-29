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
    //if(passspline && passgrid && passlc) pass=true;
    if(passspline && passlc) pass=true;
  }else{
    //if(passspline && passgrid ) pass=true;
    if(passspline) pass=true;
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
  int i,j,k,nvars,idx,idx2,idx0;
  string sorbital;
  string stemp = tmp_dir+"/"+"skdef.hsd";
  ofstream fout(stemp.c_str());
  bool hp,cp,np,op;
  bool pass; 
  hp=cp=np=op=false;

  idx0=0; 
  for(i=0;i<velem.size();i++) idx0=idx0+velem[i].lmax+2;
  if(ddh.td3) idx0=idx0+ddh.d3.size();
  if(ddh.tdamph) idx0++;
  if(ddh.thubbardderivs) idx0=idx0+ddh.hubbardderivs.size();
  if(ddh.tdamphver2) idx0=idx0+ddh.damphver2.size();

  fout.precision(6);

  fout<<"SkdefVersion = "<<skdefversion<<"\n\n";

  fout<<"Globals {\n";
  if(lc){
    fout<<"  func = xchybrid {omega = "<<omega<<"}\n";
  }else{
    fout<<"  XCFunctional ="<<xcfunctional<<"\n";
  }
  fout<<"  Superposition = "<<superposition<<"\n";
  fout<<"}"<<"\n\n";

  fout<<"AtomParameters {\n\n";
  fout<<"\
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
  }\n\n";
 
  idx=0;
  for(i=0; i<velem.size();i++){
    pass=false;
    for(j=0; j<atomparameters.size();j++){
      if(velem[i].name==atomparameters[j].name){pass=true; break;}
    }
    if(!pass){
      cerr <<"\nERROR:\nCan not find atomparameters for " << velem[i].name << endl ;
      exit(1);
    }
    fout<<"  "<<atomparameters[j].name<<" {\n";
    for(k=0; k<atomparameters[j].atomconfig.size();k++){
      fout<<atomparameters[j].atomconfig[k]<<"\n";
    }
    fout<<"    DftbAtom {\n";
    fout<<"      ShellResolved = "<<atomparameters[j].shellresolved<<"\n";
    fout.precision(velem[i].radius[0].precision);
    fout<<"      DensityCompression = PowerCompression { Power = "<<atomparameters[j].power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
    idx++;
    fout<<"      WaveCompressions = SingleAtomCompressions {\n";
    nvars=velem[i].lmax+2;
    for(k=1; k<nvars; k++){
      fout.precision(velem[i].radius[k].precision);
      if(k==1) sorbital="S"; else if(k==2) sorbital="P"; else if(k==3) sorbital="D";
      fout<<"        "<<sorbital<<" = PowerCompression { Power = "<<atomparameters[j].power<<"; Radius = "<<std::fixed<<gnome.gene(idx)<<"}\n";
      idx++;
    }
    fout<<"      }\n";
    idx2=idx0;
    if(ddh.thubbards){
      fout<<"      CustomizedHubbards {\n";
      for(k=0;k<ddh.hubbards.size();k++){ 
        if(velem[i].name==ddh.hubbards[k].name){
          fout.precision(ddh.hubbards[k].precision);
          fout<<"        "<<ddh.hubbards[k].type<<" = "<< gnome.gene(idx2)<<endl;
        }
        idx2++;
      }
      fout<<"      }\n";
    }
    if(ddh.tvorbes){
      fout<<"      CustomizedOnsites {\n";
      for(k=0;k<ddh.vorbes.size();k++){ 
        if(velem[i].name==ddh.vorbes[k].name){
          fout.precision(ddh.vorbes[k].precision);
          fout<<"        "<<ddh.vorbes[k].type<<" = "<< gnome.gene(idx2)<<endl;
        }
        idx2++;
      }
      fout<<"      }\n";
    }
    fout<<"    }\n";
    fout<<"  }\n";
  } 
  fout<<"}\n\n";
       
  fout<<"OnecenterParameters {\n\n";
  for(k=0;k<onecenterparameters.size();k++){
    fout<<onecenterparameters[k]<<"\n";
  }
  fout<<"}\n\n";

  fout<<"TwoCenterParameters {\n\n";
  for(k=0;k<twocenterparameters.size();k++){
    fout<<twocenterparameters[k]<<"\n";
  }
  fout<<"}\n\n";
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
  //  tmps=genome.gene(ie1);
  //  tmpp=genome.gene(ie1+1);
  //  if(velem[i].optype==2){
  //    genome.gene(ie1+1,tmps);
  //  }else if(velem[i].optype==3){
  //    genome.gene(ie1+1,tmps);
  //    genome.gene(ie1+2,tmps);
  //  }else if(velem[i].optype==12){
  //    genome.gene(ie1,GAMin(tmps, tmpp));
  //    genome.gene(ie1+1,GAMax(tmps, tmpp));
  //  }else if(velem[i].optype==14){
  //    genome.gene(ie1,GAMin(tmps, tmpp));
  //    genome.gene(ie1+1,GAMax(tmps, tmpp));
  //    if(genome.gene(ie1+1)-genome.gene(ie1)<0.4) genome.gene(ie1,genome.gene(ie1+1)-0.4); 
  //  }else if(velem[i].optype==21){
  //    genome.gene(ie1,GAMax(tmps, tmpp));
  //    genome.gene(ie1+1,GAMin(tmps, tmpp));
  //  }else if(velem[i].optype==11){
  //    genome.gene(ie1-1,tmps);
  //  }else if(velem[i].optype==115){
  //    genome.gene(ie1-1,1.5*tmps);
  //  }else if(velem[i].optype==111){
  //    genome.gene(ie1-1,tmps);
  //    genome.gene(ie1+1,tmps);
  //  }else if(velem[i].optype==1115){
  //    genome.gene(ie1-1,1.5*tmps);
  //    genome.gene(ie1+1,tmps);
  //  }
  //  ie1=ie1+velem[i].lmax; 
  //}
 
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
  if(ddh.tdamphver2) idx=idx+ddh.damphver2.size();

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

