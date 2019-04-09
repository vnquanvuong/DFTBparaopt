#include <iostream>
#include <vector>
#include <string>
#include <sstream> 
#include <omp.h> 
#include "erepobj.hpp"
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "tools.hpp"
#include "ga.hpp"
#include "ddh.hpp"
#include "mpi.h"

using namespace std;

double MyObjective(GAGenome& gnome);
double MyComparator(const GAGenome&, const GAGenome&);
void  MyInitializer(GAGenome&);
int   MyMutator(GAGenome&, double);
int   MyOnePointCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
int   MyTwoPointCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);

vector<string> addHamiltonian;
string popfinalfile="pop.final.dat";
string popinitialfile="pop.initial.dat";
string scratch="/tmp";
int cpu_number=1;
int power=2;
int ga_popsize=100, ga_ngen=300, ga_scoref=1, ga_flushf=1;
int preserved_num=1, seed=0;
double ga_pmut=0.2, ga_pcross=0.9;
double gtol=0.000001;
bool ga=false,readr=true,restart=false;
bool runtest=false,skfclean=true;
bool endgen=false,initialgen=true,firsteval=true;
fstream infile,outfile,restartfile;
sddh ddh;
Erepobj erepobj;
int main(int argc, char** argv){
  int i,j,k,idx,length,rank;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 if(rank==0){
  if ( argc < 2 ){
    cerr << "usage: erepobj inputfile" << endl;
    exit(1);
  }
 }
  ddh.td3=ddh.tdamph=ddh.thubbardderivs=ddh.thirdorderfull=ddh.tdamphver2=ddh.thubbards=ddh.tvorbes=false;
  
  erepobj.prepare(argv[1]);
 if(rank==0){
  erepobj.fout.open(erepobj.outfilename.c_str());
  erepobj.fout<<std::fixed;
 }
  if(ga && !runtest){
    GARandomSeed(seed);
    length=0;
    for(i=0;i<erepobj.velem.size();i++){ 
      length+=erepobj.velem[i].lmax+2;
    }
    if(ddh.td3) length+=ddh.d3.size(); 
    if(ddh.tdamph) length+=1;
    if(ddh.thubbardderivs) length+=ddh.hubbardderivs.size(); 
    if(ddh.tdamphver2) length+=ddh.damphver2.size(); 
    if(ddh.thubbards) length+=ddh.hubbards.size(); 
    if(ddh.tvorbes) length+=ddh.vorbes.size(); 

    GA1DArrayGenome<double> genome(length, MyObjective);
    genome.initializer(::MyInitializer);
    genome.mutator(::MyMutator);
    genome.comparator(::MyComparator);
    MyGASimpleGA ga(genome);
    ga.crossover(::MyTwoPointCrossover);
  //ga.maximize();
    ga.minimize();
  //GANoScaling scale;
    GASigmaTruncationScaling scale;
    ga.scaling(scale);
    GARouletteWheelSelector select;
    ga.selector(select);
    ga.populationSize(ga_popsize);
    ga.nGenerations(ga_ngen);
    ga.pMutation(ga_pmut);
    ga.pCrossover(ga_pcross);
    ga.scoreFilename("score.dat");  // name of file for scores
    ga.scoreFrequency(ga_scoref); // keep the scores of every 10th generation
    ga.flushFrequency(ga_flushf); // specify how often to write the score to disk
    ga.selectScores(GAStatistics::Minimum);
  //ga.parameters(argc, argv, gaTrue); // parse commands, complain if bogus args
    
 if(rank==0){
    erepobj.fout << "#initializing...\n"; erepobj.fout.flush();
 }

  //string call="rm -rf skftmp; mkdir -p skftmp;";
  //i=system(call.c_str());	
  
    if(restart==true) restartfile.open(popfinalfile.c_str(), ios::in);
    ga.initialize();
    if(restart==true) restartfile.close();
    initialgen=false;
    ga.initialstep();
   if(rank==0){
    outfile.open(popinitialfile.c_str(), (STD_IOS_OUT | STD_IOS_TRUNC));
    for(i=0; i<ga.population().size(); i++){
      genome = ga.population().individual(i);
      idx=0;
      outfile<<std::fixed;
      for(j=0;j<erepobj.velem.size();j++){
        for(k=0;k<erepobj.velem[j].lmax+2.;k++){
          outfile.precision(erepobj.velem[j].radius[k].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;  
        }
      }
      if(ddh.td3){
        for(j=0;j<ddh.d3.size();j++){ 
          outfile.precision(ddh.d3[j].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
      }
      if(ddh.tdamph){
        outfile.precision(ddh.damph.precision);
        outfile.width(5);
        outfile << genome.gene(idx) << "\t";
        idx++;
      }
      if(ddh.thubbardderivs){
        for(j=0;j<ddh.hubbardderivs.size();j++){ 
          outfile.precision(ddh.hubbardderivs[j].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
      }
      if(ddh.tdamphver2){
        for(j=0;j<ddh.damphver2.size();j++){ 
          outfile.precision(ddh.damphver2[j].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
      }
      if(ddh.thubbards){
        for(j=0;j<ddh.hubbards.size();j++){ 
          outfile.precision(ddh.hubbards[j].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
      }
      if(ddh.tvorbes){
        for(j=0;j<ddh.vorbes.size();j++){ 
          outfile.precision(ddh.vorbes[j].precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
      }
      outfile.precision(9);
      outfile.width(20);
      outfile << genome.score() << "\n";
    }
    outfile.close();

    erepobj.fout << "#evolving...\n"; erepobj.fout.flush();
  }
    while(!ga.done()) {
     if(rank==0){
    //cout<<"generation:  "<<ga.generation()<<" start\n"; cout.flush();
      erepobj.fout.precision(2);
      erepobj.fout<<double(ga.generation())/double(ga_ngen)<<"   "; erepobj.fout.flush();
     }
      ga.step();
     if(rank==0){
      ga.flushScores();
      erepobj.fout.precision(9);
//    erepobj.fout<<ga.population().individual(0).score()<<endl; cout.flush();
      erepobj.fout<<ga.population().best().score()<<endl; cout.flush();
 
      outfile.open(popfinalfile.c_str(), (STD_IOS_OUT | STD_IOS_TRUNC));
      for(i=0; i<ga.population().size(); i++){
        genome = ga.population().individual(i);
        idx=0;
        outfile<<std::fixed;
        for(j=0;j<erepobj.velem.size();j++){
          for(k=0;k<erepobj.velem[j].lmax+2.;k++){
            outfile.precision(erepobj.velem[j].radius[k].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;  
          }
        }
        if(ddh.td3){
          for(j=0;j<ddh.d3.size();j++){ 
            outfile.precision(ddh.d3[j].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;
          }
        }
        if(ddh.tdamph){
          outfile.precision(ddh.damph.precision);
          outfile.width(5);
          outfile << genome.gene(idx) << "\t";
          idx++;
        }
        if(ddh.thubbardderivs){
          for(j=0;j<ddh.hubbardderivs.size();j++){ 
            outfile.precision(ddh.hubbardderivs[j].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;
          }
        }
        if(ddh.tdamphver2){
          for(j=0;j<ddh.damphver2.size();j++){ 
            outfile.precision(ddh.damphver2[j].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;
          }
        }
        if(ddh.thubbards){
          for(j=0;j<ddh.hubbards.size();j++){ 
            outfile.precision(ddh.hubbards[j].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;
          }
        }
        if(ddh.tvorbes){
          for(j=0;j<ddh.vorbes.size();j++){ 
            outfile.precision(ddh.vorbes[j].precision);
            outfile.width(5);
            outfile << genome.gene(idx) << "\t";
            idx++;
          }
        }
        outfile.precision(9);
        outfile.width(20);
        outfile << genome.score() << "\n";
//      for(j=0;j<genome.length();j++) {
//        outfile << genome.gene(j) << "\t";
//      }
//      outfile << genome.score() << "\n";
      }
      outfile.close();
      if(skfclean){
        string call="rm -f "+erepobj.libdir+"/*";
        i=system(call.c_str());	
      }
     }
    }
   if(rank==0){
    ga.flushScores();
    erepobj.fout<<"#ga end!\n";
    genome = ga.population().best();
    endgen=true;
    erepobj.get_residual(genome); 
    idx=0;
    for(j=0;j<erepobj.velem.size();j++){
      for(k=0;k<erepobj.velem[j].lmax+2.;k++){
         erepobj.velem[j].radius[k].r=genome.gene(idx);
        idx++;  
      }
    }
    if(ddh.td3){
      for(j=0;j<ddh.d3.size();j++){ 
        ddh.d3[j].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.tdamph){
      ddh.damph.value=genome.gene(idx);
      idx++;
    }
    if(ddh.thubbardderivs){
      for(j=0;j<ddh.hubbardderivs.size();j++){ 
         ddh.hubbardderivs[j].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.tdamphver2){
      for(j=0;j<ddh.damphver2.size();j++){ 
         ddh.damphver2[j].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.thubbards){
      for(j=0;j<ddh.hubbards.size();j++){ 
         ddh.hubbards[j].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.tvorbes){
      for(j=0;j<ddh.vorbes.size();j++){ 
         ddh.vorbes[j].value=genome.gene(idx);
        idx++;
      }
    }
   }
  }

 if(rank==0){
  erepobj.writeout();

  erepobj.fout << "** erepobj normal termination **" << endl << endl;
  erepobj.fout.close();
 }
  MPI_Finalize();
  return 0;
}



double MyObjective(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int ie1,ie2,i,j,k,idx,nvars;
  double value,tscore,tmp,tmps,tmpp;
  int id=0;
  char ctline[512];

  idx=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    nvars=erepobj.velem[i].lmax+2;
    for(j=0; j<nvars; j++){
      value  = genome.gene(idx);
      value  = GAMax(erepobj.velem[i].radius[j].minr, value);
      value  = GAMin(erepobj.velem[i].radius[j].maxr, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.d3[i].min, value);
      value  = GAMin(ddh.d3[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tdamph){
    value  = genome.gene(idx);
    value  = GAMax(ddh.damph.min, value);
    value  = GAMin(ddh.damph.max, value);
    genome.gene(idx,value);
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.hubbardderivs[i].min, value);
      value  = GAMin(ddh.hubbardderivs[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.damphver2[i].min, value);
      value  = GAMin(ddh.damphver2[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.thubbards){
    for(i=0;i<ddh.hubbards.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.hubbards[i].min, value);
      value  = GAMin(ddh.hubbards[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tvorbes){
    for(i=0;i<ddh.vorbes.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.vorbes[i].min, value);
      value  = GAMin(ddh.vorbes[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }



  ie1=-1;
  for(i=0;i<erepobj.velem.size();i++){
    ie1=ie1+2;  
    tmps=genome.gene(ie1);
    tmpp=genome.gene(ie1+1);
    if(erepobj.velem[i].optype==2){
      genome.gene(ie1+1,tmps);
    }else if(erepobj.velem[i].optype==3){
      genome.gene(ie1+1,tmps);
      genome.gene(ie1+2,tmps);
    }else if(erepobj.velem[i].optype==12){
      genome.gene(ie1,GAMin(tmps, tmpp));
      genome.gene(ie1+1,GAMax(tmps, tmpp));
    }else if(erepobj.velem[i].optype==14){
      genome.gene(ie1,GAMin(tmps, tmpp));
      genome.gene(ie1+1,GAMax(tmps, tmpp));
      if(genome.gene(ie1+1)-genome.gene(ie1)<0.4) genome.gene(ie1,genome.gene(ie1+1)-0.4); 
    }else if(erepobj.velem[i].optype==21){
      genome.gene(ie1,GAMax(tmps, tmpp));
      genome.gene(ie1+1,GAMin(tmps, tmpp));
    }else if(erepobj.velem[i].optype==11){
      genome.gene(ie1-1,tmps);
    }else if(erepobj.velem[i].optype==111){
      genome.gene(ie1-1,tmps);
      genome.gene(ie1+1,tmps);
    }
    ie1=ie1+erepobj.velem[i].lmax; 
  }

  if(!initialgen){
    ostringstream ss;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    string call,stemp,part1,part2,part3; 
    string skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);

    call="rm -rf "+skf_dir+"; mkdir -p "+skf_dir+"; ";
//  call=" mkdir -p "+skf_dir+"; ";
     
    idx=0; 
    for(i=0;i<erepobj.velem.size();i++) idx=idx+erepobj.velem[i].lmax+2;
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
    for(i=0;i<erepobj.velem.size();i++){
      part1=erepobj.velem[i].name;  
      ss.str(std::string());
      ss.precision(erepobj.velem[i].radius[0].precision); ss<<std::fixed<<genome.gene(ie1);
      part1+="_d"+ss.str();
      ie1++;  
      for(k=1;k<erepobj.velem[i].lmax+2.;k++){
        ss.str(std::string());
        ss.precision(erepobj.velem[i].radius[k].precision); ss<<std::fixed<<genome.gene(ie1);
        part1+="_w"+ss.str();  
        ie1++; 
      }
      ie2=0;
      for(j=0;j<erepobj.velem.size();j++){
        part2=erepobj.velem[j].name;  
        ss.str(std::string());
        ss.precision(erepobj.velem[j].radius[0].precision); ss<<std::fixed<<genome.gene(ie2);
        part2+="_d"+ss.str();  
        ie2++;
        for(k=1;k<erepobj.velem[j].lmax+2.;k++){
          ss.str(std::string());
          ss.precision(erepobj.velem[j].radius[k].precision); ss<<std::fixed<<genome.gene(ie2);
          part2+="_w"+ss.str();  
          ie2++;
        }
        if(erepobj.lc){
          ss.str(std::string());
          ss.precision(2);ss<<std::fixed<<erepobj.omega;
          stemp="omega_"+ss.str()+"_"+part1+"-"+part2+part3+".skf";
        }else{
          stemp=part1+"-"+part2+part3+".skf";
        }
        erepobj.checkskf(stemp,1);
        call+=" cp "+erepobj.libdir+"/"+stemp+" "+ skf_dir+"/"+erepobj.velem[i].name+"-"+erepobj.velem[j].name+".skf; ";
      }
    }
    if(erepobj.addskf) call+=" cp -f "+erepobj.libadddir+"/*.skf " + skf_dir+"/; ";

    i=system(call.c_str());	

    if (erepobj.fit_type==0){
      fstream tinfile, toutfile;
      tinfile.open(erepobj.frep_in.c_str(),ios::in);
      if(!tinfile.is_open()){
        cout<<"unable to open "<<erepobj.frep_in<<" file.\n";
        exit(1);
      }
      call = skf_dir+"/rep.in";
      toutfile.open(call.c_str(), ios::out);
      if(!toutfile.is_open()){
        cout<<"unable to open "<<call<<" file.\n";
        exit(1);
      }
      toutfile<<std::fixed;
      if( ddh.td3 || ddh.tdamph || ddh.thubbardderivs || ddh.tdamphver2){
        idx=ie2; 
        if(ddh.td3){
          toutfile<<"\n$d3:\n"; 
          for(i=0;i<ddh.d3.size();i++){ 
            toutfile<<ddh.d3[i].name<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<ddh.d3[i].precision<<endl;
            idx++;
          }
          toutfile<<"$end\n"; 
        }
        if(ddh.tdamph){
          toutfile<<"\n$damph:\n"; 
          toutfile<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<ddh.damph.precision<<endl;
          idx++;
          toutfile<<"$end\n"; 
        }
        if(ddh.thubbardderivs){
          if(ddh.thirdorderfull) toutfile<<"\n$hubbardderivsfull:\n"; 
          else toutfile<<"\n$hubbardderivs:\n"; 
          for(i=0;i<ddh.hubbardderivs.size();i++){ 
            toutfile<<ddh.hubbardderivs[i].name<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<ddh.hubbardderivs[i].precision<<endl;
            idx++;
          }
          toutfile<<"$end\n"; 
        }
        if(ddh.tdamphver2){
          toutfile<<"\n$damphver2:\n"; 
          for(i=0;i<ddh.damphver2.size();i++){ 
            toutfile<<ddh.damphver2[i].name<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<genome.gene(idx)<<" "<<ddh.damphver2[i].precision<<endl;
            idx++;
          }
          toutfile<<"$end\n"; 
        }
      }
      while(tinfile.getline(ctline,512)){
        toutfile<<ctline<<endl;
      }
      tinfile.close();
      toutfile.close();
    }

    tscore=erepobj.score(id); 
   
  //call="rm -rf "+skf_dir+";";
  //i=system(call.c_str());	

    return tscore;
  }else{
////////////Critical////////////////////////////////
    //if(firsteval){
    //  firsteval=false;
    //  return 0.00000099;
    //} 
    return 0.000001;
////////////////////////////////////////////
  }
}

  //for(i=0;i<erepobj.velem.size();i++){
  //  part1=erepobj.velem[i].name;  
  //  ss.precision(erepobj.velem[i].radius[0].precision); ss<<std::fixed<<erepobj.velem[i].radius[0].r;
  //  part1+="_d"+ss.str(); 
  //  for(k=1;k<erepobj.velem[i].lmax+2.;k++){
  //    ss.precision(erepobj.velem[i].radius[k].precision); ss<<std::fixed<<erepobj.velem[i].radius[k].r;
  //    part1+="_w"+ss.str();  
  //  }
  //  for(j=0;j<erepobj.velem.size();j++){
  //    part2=erepobj.velem[j].name;  
  //    ss.precision(erepobj.velem[j].radius[0].precision); ss<<std::fixed<<erepobj.velem[j].radius[0].r;
  //    part2+="_d"+ss.str();  
  //    for(k=1;k<erepobj.velem[j].lmax+2.;k++){
  //      ss.precision(erepobj.velem[j].radius[k].precision); ss<<std::fixed<<erepobj.velem[j].radius[k].r;
  //      part2+="_w"+ss.str();  
  //    }
  //    stemp=erepobj.libdir+"/"+part1+"-"+part2+".skf";
  //    erepobj.checkskf(stemp);
  //    call+=" cp "+stemp+" "+ skf_dir+"; ";
  //  }
  //}

void MyInitializer(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  bool accept; 
  int ie1,i,j,k,idx,nknots,npb,nvars;
  double tmp,value;

  if(readr==true){
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        genome.gene(idx,erepobj.velem[i].radius[k].r);
        idx++;
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        genome.gene(idx,ddh.d3[i].value);
        idx++;
      }
    }
    if(ddh.tdamph){
      genome.gene(idx,ddh.damph.value);
      idx++;
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        genome.gene(idx,ddh.hubbardderivs[i].value);
        idx++;
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        genome.gene(idx,ddh.damphver2[i].value);
        idx++;
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        genome.gene(idx,ddh.hubbards[i].value);
        idx++;
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        genome.gene(idx,ddh.vorbes[i].value);
        idx++;
      }
    }
    readr=false;
  }else{
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        genome.gene(idx,GARandomFloat(erepobj.velem[i].radius[k].minr,erepobj.velem[i].radius[k].maxr));
        idx++;
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        genome.gene(idx,GARandomFloat(ddh.d3[i].min,ddh.d3[i].max));
        idx++;
      }
    }
    if(ddh.tdamph){
      genome.gene(idx,GARandomFloat(ddh.damph.min,ddh.damph.max));
      idx++;
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        genome.gene(idx,GARandomFloat(ddh.hubbardderivs[i].min,ddh.hubbardderivs[i].max));
        idx++;
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        genome.gene(idx,GARandomFloat(ddh.damphver2[i].min,ddh.damphver2[i].max));
        idx++;
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        genome.gene(idx,GARandomFloat(ddh.hubbards[i].min,ddh.hubbards[i].max));
        idx++;
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        genome.gene(idx,GARandomFloat(ddh.vorbes[i].min,ddh.vorbes[i].max));
        idx++;
      }
    }
  }
    
  if(restart==true){
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    if(ddh.tdamph){
      tmp=999999999.9;
      restartfile>>tmp;
      if(tmp!=999999999.9){
        genome.gene(idx,tmp);
        idx++;
      }
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
          idx++;
        }
      }
    }
    restartfile>>tmp;
  }


  idx=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    nvars=erepobj.velem[i].lmax+2;
    for(j=0; j<nvars; j++){
      value  = genome.gene(idx);
      value  = GAMax(erepobj.velem[i].radius[j].minr, value);
      value  = GAMin(erepobj.velem[i].radius[j].maxr, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.d3[i].min, value);
      value  = GAMin(ddh.d3[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tdamph){
    value  = genome.gene(idx);
    value  = GAMax(ddh.damph.min, value);
    value  = GAMin(ddh.damph.max, value);
    genome.gene(idx,value);
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.hubbardderivs[i].min, value);
      value  = GAMin(ddh.hubbardderivs[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.damphver2[i].min, value);
      value  = GAMin(ddh.damphver2[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.thubbards){
    for(i=0;i<ddh.hubbards.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.hubbards[i].min, value);
      value  = GAMin(ddh.hubbards[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }
  if(ddh.tvorbes){
    for(i=0;i<ddh.vorbes.size();i++){ 
      value  = genome.gene(idx);
      value  = GAMax(ddh.vorbes[i].min, value);
      value  = GAMin(ddh.vorbes[i].max, value);
      genome.gene(idx,value);
      idx++;
    }
  }


  //erepobj.makeskf(genome);
}

int MyMutator(GAGenome& g, double pmut){
  GA1DArrayGenome<double>& child = (GA1DArrayGenome<double>&)g;
  int i,j,k,idx,nvars,iMut;
  double value;

 // cout<<"in mutator start\n";
  if(pmut <= 0.0) return(0);
  idx=0;
  iMut=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    nvars=erepobj.velem[i].lmax+2;
    for(j=0; j<nvars; j++){
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(erepobj.velem[i].radius[j].delta);
        value += child.gene(idx);
        value  = GAMax(erepobj.velem[i].radius[j].minr, value);
        value  = GAMin(erepobj.velem[i].radius[j].maxr, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.d3[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.d3[i].min, value);
        value  = GAMin(ddh.d3[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tdamph){
    if(GAFlipCoin(pmut)){
      value  = GAGaussianFloat(ddh.damph.delta);
      value += child.gene(idx);
      value  = GAMax(ddh.damph.min, value);
      value  = GAMin(ddh.damph.max, value);
      child.gene(idx,value);
      iMut++;
    }
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.hubbardderivs[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.hubbardderivs[i].min, value);
        value  = GAMin(ddh.hubbardderivs[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.damphver2[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.damphver2[i].min, value);
        value  = GAMin(ddh.damphver2[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.thubbards){
    for(i=0;i<ddh.hubbards.size();i++){ 
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.hubbards[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.hubbards[i].min, value);
        value  = GAMin(ddh.hubbards[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tvorbes){
    for(i=0;i<ddh.vorbes.size();i++){ 
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.vorbes[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.vorbes[i].min, value);
        value  = GAMin(ddh.vorbes[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  return(iMut);
}


double MyComparator(const GAGenome& g1, const GAGenome& g2) {
  GA1DArrayGenome<double>& gnome1 = (GA1DArrayGenome<double>&)g1;
  GA1DArrayGenome<double>& gnome2 = (GA1DArrayGenome<double>&)g2;
  double diff=0.0; 
  int i;
  for(i=0;i<gnome1.length();i++){ 
    diff+=abs(gnome2.gene(i)-gnome1.gene(i));
  }
//return diff;
  return 1; 
}


int MyOnePointCrossover(const GAGenome& p1, const GAGenome& p2,
  GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<double> &mom=DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int nc=0;
  unsigned int momsite, momlen;
  unsigned int dadsite, dadlen;

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE &&
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() ||
      sis.length() != bro.length() ||
      sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite = dadsite = 3*GARandomInt(0, mom.length()/3);
      momlen = dadlen = mom.length() - momsite;
    }else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameBehavReqd);
      return nc;
    }else{
      momsite = 3*GARandomInt(0, mom.length()/3);
      dadsite = 3*GARandomInt(0, dad.length()/3);
      momlen = mom.length() - momsite;
      dadlen = dad.length() - dadsite;
      sis.resize(momsite+dadlen);
      bro.resize(dadsite+momlen);
    }

    sis.copy(mom, 0, 0, momsite);
    sis.copy(dad, momsite, dadsite, dadlen);
    bro.copy(dad, 0, 0, dadsite);
    bro.copy(mom, dadsite, momsite, momlen);
    nc = 2;
  }else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
    DYN_CAST(GA1DArrayGenome<double> &, *c1) :
    DYN_CAST(GA1DArrayGenome<double> &, *c2));

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite = dadsite = 3*GARandomInt(0, mom.length()/3);
      momlen = dadlen = mom.length() - momsite;
    }else{
      momsite = 3*GARandomInt(0, mom.length()/3);
      dadsite = 3*GARandomInt(0, dad.length()/3);
      momlen = mom.length() - momsite;
      dadlen = dad.length() - dadsite;
      sis.resize(momsite+dadlen);
    }

    if(GARandomBit()){
      sis.copy(mom, 0, 0, momsite);
      sis.copy(dad, momsite, dadsite, dadlen);
    }else{
      sis.copy(dad, 0, 0, dadsite);
      sis.copy(mom, dadsite, momsite, momlen);
    }
    nc = 1;
  }
  return nc;
}

int MyTwoPointCrossover(const GAGenome& p1, const GAGenome& p2,
  GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<double> &mom=DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int nc=0;
  unsigned int momsite[2], momlen[2];
  unsigned int dadsite[2], dadlen[2];

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE &&
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() ||
      sis.length() != bro.length() ||
      sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "two-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      //cout<<"mom.length()/3:   "<<mom.length()/3<<endl;
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0] = momlen[0];
      dadlen[1] = momlen[1];
    }else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      return nc;
    }else{
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = 3*GARandomInt(0, dad.length()/3);
      dadsite[1] = 3*GARandomInt(0, dad.length()/3);
      if(dadsite[0] > dadsite[1]) SWAP(dadsite[0], dadsite[1]);
      dadlen[0] = dadsite[1] - dadsite[0];
      dadlen[1] = dad.length() - dadsite[1];

      sis.resize(momsite[0]+dadlen[0]+momlen[1]);
      bro.resize(dadsite[0]+momlen[0]+dadlen[1]);
    }

    sis.copy(mom, 0, 0, momsite[0]);
    sis.copy(dad, momsite[0], dadsite[0], dadlen[0]);
    sis.copy(mom, momsite[0]+dadlen[0], momsite[1], momlen[1]);
    bro.copy(dad, 0, 0, dadsite[0]);
    bro.copy(mom, dadsite[0], momsite[0], momlen[0]);
    bro.copy(dad, dadsite[0]+momlen[0], dadsite[1], dadlen[1]);

    nc = 2;
   
  }else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
    DYN_CAST(GA1DArrayGenome<double> &, *c1) :
    DYN_CAST(GA1DArrayGenome<double> &, *c2));

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "two-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0] = momlen[0];
      dadlen[1] = momlen[1];
    }else{
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = 3*GARandomInt(0, dad.length()/3);
      dadsite[1] = 3*GARandomInt(0, dad.length()/3);
      if(dadsite[0] > dadsite[1]) SWAP(dadsite[0], dadsite[1]);
      dadlen[0] = dadsite[1] - dadsite[0];
      dadlen[1] = dad.length() - dadsite[1];

      sis.resize(momsite[0]+dadlen[0]+momlen[1]);
    }

    if(GARandomBit()){
      sis.copy(mom, 0, 0, momsite[0]);
      sis.copy(dad, momsite[0], dadsite[0], dadlen[0]);
      sis.copy(mom, momsite[0]+dadlen[0], momsite[1], momlen[1]);
    }else{
      sis.copy(dad, 0, 0, dadsite[0]);
      sis.copy(mom, dadsite[0], momsite[0], momlen[0]);
      sis.copy(dad, dadsite[0]+momlen[0], dadsite[1], dadlen[1]);
    }

    nc = 1;
  
  }

  return nc;
}


