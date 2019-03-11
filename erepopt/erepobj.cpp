#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include "ga/ga.h"
#include <ga/std_stream.h>
#include "erepobj.hpp"
#include "tools.hpp"
#include "ddh.hpp"

extern sddh ddh;
extern string scratch;
extern int cpu_number;

using namespace std;
extern bool runtest;

Erepobj::Erepobj(){
}

Erepobj::Erepobj(const string inputfile) {
  readinp(inputfile); 
}

void Erepobj::prepare(const string inputfile) {
  readinp(inputfile); 
}

double Erepobj::score(int id){
  int i,j,iq,ia;
  double a,score,tscore,tmscore;
  string stemp,skf_dir,call;
  char cline[512];
  skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);


  if (fit_type==1){
    int index,length;
    double  * ev;
    length=0;
    for(i=0;i<vmol.size();i++){
      length=length+ vmol[i].nelectron/2;
    }
    ev = new double [length];
    index=0;
    for(i=0;i<vmol.size();i++){
      vmol[i].getev(id,index,ev,dftbversion);
    }
    a=0.1;
    score=999999999.9;
    for(ia=0;ia<9000;ia++){
      index=0;
      tscore=0.0;
      a=a+0.0001;
      a=1.0;
      for(i=0;i<vmol.size();i++){
        tmscore=0.0;
        for(j=0;j<vmol[i].nelectron/2;j++){
          if(score_type==0) tmscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==1) tmscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==2) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
          else if(score_type==4) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
          index++;
        }
        tscore=tscore+abs(tmscore);
      }
      if(tscore < score) score=tscore;
    }
    delete  [] ev;
  }else if (fit_type==2){
    score=0.0;
    for(i=0;i<vmol.size();i++){
      score+=vmol[i].score(id,score_type,dftbversion,false);
    } 
  }else{
    score=999999999.9;
    call  = "cp -rf "+fgrid+" "+skf_dir+";";
//  call += "cp -rf "+fdftb_inp+" "+skf_dir+";";
//  call += "cp "+frep_in+" "+skf_dir+"/rep.in;";
//  call += "cd "+skf_dir+"/; "+ gasrepfit +" "+ frep_in +">rep.out; " + "cd ../ ";
    call += "cd "+skf_dir+"/; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";"+ gasrepfit +" rep.in >rep.out; " + "cd ../ ";

    iq = system(call.c_str());	

    stemp=skf_dir+"/rep.out";
    ifstream fin(stemp.c_str());
    if ( !fin ) {
      cerr << endl << "ERROR: " << endl << endl;
      cerr << "The system call: \"" << call << "\" did not produce the file: \""<<skf_dir<<"/rep.out\""
           << endl << "calculation of gasrepfit was failed."
           << endl << "exit erepobj" << endl << endl;
      exit(1);
    }
    while(fin.getline(cline,512)){
      stemp=cline;
      if(stemp.find("final_score:")!=string::npos){
        fin >> score;
        break;
      }
    }
    fin.close();
  }

  return score;
}

void Erepobj::get_residual(const GAGenome& g) {
  double tmp,tmps,tmpp;
  restotU=0; restotS=0; restot2=0; restot4=0;
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int ie1,ie2,i,j,k,ia;
  double tscore,tmscore;
  int id=0,idx;

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
  
//if (fit_type!=0){
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
 
    ostringstream ss;
//  id=omp_get_thread_num();
    id=0;
    string call,stemp,part1,part2,part3; 
    string skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);
//  call="rm -rf "+skf_dir+"; mkdir -p "+skf_dir+"; ";
    call=" mkdir -p "+skf_dir+"; ";

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
      velem[i].radius[0].r=genome.gene(ie1);
      part1+="_d"+ss.str();
      ie1++;  
      for(k=1;k<velem[i].lmax+2.;k++){
        ss.str(std::string());
        ss.precision(velem[i].radius[k].precision); ss<<std::fixed<<genome.gene(ie1);
        velem[i].radius[k].r=genome.gene(ie1);
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
        }else{
          stemp=part1+"-"+part2+part3+".skf";
        }
        if (fit_type!=0){
           checkskf(stemp,1);
           call+=" cp "+libdir+"/"+stemp+" "+ skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf; ";
        }
      }
    }
    if(addskf) call+=" cp -f "+libadddir+"/*.skf " + skf_dir+"/; ";

    i=system(call.c_str());	
//}
  
  //tscore=score(id); 

  //tscore=0.0;
  //for(i=0;i<vmol.size();i++){
  //  tscore+=vmol[i].score(id,score_type,dftbversion,true);
  //}

  if (fit_type==1){
    int index,length;
    double  * ev,a,af,fscore;
    length=0;
    for(i=0;i<vmol.size();i++){
      length=length+ vmol[i].nelectron/2;
    }
    ev = new double [length];
    index=0;
    for(i=0;i<vmol.size();i++){
      vmol[i].getev(id,index,ev,dftbversion);
    }
    //index=0;
    //for(i=0;i<vmol.size();i++){
    //  for(j=0;j<vmol[i].nelectron/2;j++){
    //    cout<<ev[index]<<endl;
    //    index++;
    //  }
    //}
    af=a=0.1;
    fscore=999999999.9;
    for(ia=0;ia<9000;ia++){
      index=0;
      tscore=0.0;
      a=a+0.0001;
      a=1.0;
      for(i=0;i<vmol.size();i++){
        tmscore=0.0;
        for(j=0;j<vmol[i].nelectron/2;j++){
          if(score_type==0) tmscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==1) tmscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==2) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
          else if(score_type==4) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
          index++;
        }
        tscore=tscore+abs(tmscore);
      }
      if(tscore < fscore) {fscore=tscore; af=a;}
    }
    a=af;
    index=0;
    for(i=0;i<vmol.size();i++){
      tscore=0.0;
      for(j=0;j<vmol[i].nelectron/2;j++){
        if(score_type==0) tscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
        else if(score_type==1) tscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
        else if(score_type==2) tscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
        else if(score_type==4) tscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
        index++;
      }
      vmol[i].error=tscore;
    }
    evscaling=af;
    delete  [] ev;
  }else if (fit_type==2){
    tscore=0.0;
    for(i=0;i<vmol.size();i++){
      tscore+=vmol[i].score(id,score_type,dftbversion,true);
    }
  }
}

void Erepobj::get_residual() {
  double tmp;
  restotU=0; restotS=0; restot2=0; restot4=0;

}

