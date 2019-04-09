#include <iostream>
#include <iomanip>  
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <ga/ga.h>
#include "ga.hpp"
#include "tools.hpp"
#include "auxiliary.hpp"
#include "allequations.hpp"


using namespace std;

bool   ga=true,read_spline=true,fsmooth_spline=true,aa_spline=true;
bool   fderivative1st=true,fderivative2nd=true,fderivative3rd=false;
bool   runtest=false,endgen=false,grid_update=true;
int    ga_popsize=3000,ga_ngen=1000,ga_scoref=1,ga_flushf=1;
int    score_type=2,preserved_num=300,destroy_num=30,popsizemin=2,seed=0;
double ga_pmut=0.2,ga_pcross=0.9,gauss_dev=0.025,gtol=0.000001,deltar=0.0001; 
int    ilmsfit=4,idecompose=6,nreplicate=1;

sddh         ddh;
vector<spot> svpot; 
Allequations allequations;

double MyObjective(GAGenome&);
void   MyInitializer(GAGenome&);
int    MyMutator(GAGenome&,double);

string  inputfilename;  
fstream infile, outfile;

int main(int argc, char** argv){
  int i,j,length,k,idx; 
  double tmp;
  Timer gtime;
  ddh.td3=ddh.tdamph=ddh.thubbardderivs=ddh.thirdorderfull=ddh.tdamphver2=false;
  if ( argc < 2 ){
    cerr << "usage: repopt inputfile" << endl;
    exit(1);
  }
  inputfilename = argv[1];
  allequations.prepare(inputfilename.c_str());

  if(ga && !runtest){
    GARandomSeed(seed);
    length=0;
    for(i=0;i<allequations.vpot.size();i++){ 
      svpot.resize(i+1);
      svpot[i].potname       =allequations.vpot[i].potname;
      svpot[i].gridname      =allequations.vpot[i].gridname;
      svpot[i].nknots        =allequations.vpot[i].nknots;
      svpot[i].ordspl        =allequations.vpot[i].ordspl;
      svpot[i].minRbond      =allequations.vpot[i].minRbond;
      svpot[i].minr          =allequations.vpot[i].minRbond-allequations.vpot[i].max_expand;
      svpot[i].max_step      =allequations.vpot[i].max_step;
      svpot[i].smooth_order  =allequations.vpot[i].smooth_order;
      svpot[i].aa            =allequations.vpot[i].aa;
      svpot[i].division.resize(svpot[i].nknots-1);
      length+=allequations.vpot[i].nknots;
      for(j=0;j<allequations.vpot[i].nknots;j++){ 
        svpot[i].vr.push_back(allequations.vpot[i].vr[j]);
      }
    }
    if(ddh.td3) length+=ddh.d3.size(); 
    if(ddh.tdamph) length+=1;
    if(ddh.thubbardderivs) length+=ddh.hubbardderivs.size(); 
    if(ddh.tdamphver2) length+=ddh.damphver2.size(); 

    GA1DArrayGenome<double> genome(length, MyObjective);
    genome.initializer(::MyInitializer);
    genome.mutator(::MyMutator);
    genome.comparator(::MyComparator);
    MyGASimpleGA ga(genome);
    ga.crossover(::MyTwoPointCrossover);
    ga.minimize();
    GASigmaTruncationScaling scale;
    ga.scaling(scale);
    GARouletteWheelSelector select;
    ga.selector(select);
    ga.populationSize(ga_popsize);
    ga.nGenerations(ga_ngen);
    ga.pMutation(ga_pmut);
    ga.pCrossover(ga_pcross);
    ga.scoreFilename("score.dat");  
    ga.scoreFrequency(ga_scoref);   
    ga.flushFrequency(ga_flushf);   
    ga.selectScores(GAStatistics::Minimum);
    
    cout << "#initializing...\n"; cout.flush();
   
    ga.initialize();

    outfile.open("pop.initial.dat", (STD_IOS_OUT | STD_IOS_TRUNC));
    for(i=0; i<ga.population().size(); i++){
      genome = ga.population().individual(i);
      for(j=0;j<genome.length();j++) {
        outfile<<fixed<<setprecision(6)<<setw(10)<< genome.gene(j) << " ";
      }
      outfile<<fixed<<setprecision(9)<<setw(20)<< genome.score() << "\n";
    }
    outfile.close();

    cout << "#evolving...\n"; cout.flush();
    while(!ga.done()) {
      cout<<fixed<<setprecision(6)<<left<<setw(10)<<double(ga.generation())/double(ga_ngen); 
      cout.flush();
      ga.step();
      cout<<fixed<<setprecision(9)<<right<<setw(20)<<ga.population().individual(0).score()<<endl; 
      cout.flush();
    }
    cout.precision(9);
    cout<<"final_score: \n"<<ga.population().individual(0).score()<<endl; cout.flush();
    //cout << "\n\n";
    ga.flushScores();
    cout<<"#ga end!\n";
    outfile.open("pop.final.dat", (STD_IOS_OUT | STD_IOS_TRUNC));
    for(i=0; i<ga.population().size(); i++){
      genome = ga.population().individual(i);
      for(j=0;j<genome.length();j++) {
        outfile<<fixed<<setprecision(6)<<setw(10)<< genome.gene(j) << " ";
      }
      outfile<<fixed<<setprecision(9)<<setw(20)<< genome.score() << "\n";
    }
    outfile.close();

    genome = ga.population().best();
    idx=0;
    for(i=0;i<allequations.vpot.size();i++){ 
      for(j=0;j<allequations.vpot[i].nknots;j++){ 
        allequations.vpot[i].vr[j]=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        ddh.d3[i].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.tdamph){
      ddh.damph.value=genome.gene(idx);
      idx++;
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        ddh.hubbardderivs[i].value=genome.gene(idx);
        idx++;
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        ddh.damphver2[i].value=genome.gene(idx);
        idx++;
      }
    }
    for(i=0;i<allequations.vpot.size();i++){ 
      for(j=1;j<allequations.vpot[i].nknots-1;j++){
        for(k=j+1;k<allequations.vpot[i].nknots-1;k++){
          if(allequations.vpot[i].vr[k]<allequations.vpot[i].vr[j]){
            tmp=allequations.vpot[i].vr[j];
            allequations.vpot[i].vr[j]=allequations.vpot[i].vr[k];
            allequations.vpot[i].vr[k]=tmp;
          }
        }
      }
    }
    idx=0;
    for(i=0;i<allequations.vpot.size();i++){ 
      for(j=0;j<allequations.vpot[i].nknots;j++){ 
        genome.gene(idx,allequations.vpot[i].vr[j]);
        idx++;
      }
    }
    if(grid_update){
      for(i=0;i<allequations.vpot.size();i++){ 
        outfile.open(allequations.vpot[i].gridname.c_str(),ios::out);
        if(!outfile.is_open()){
          cout<<"unable to open"<< allequations.vpot[i].gridname<< " file.\n";
          return 0;
        }else{
          for(j=0;j<allequations.vpot[i].nknots;j++){ 
            outfile<<fixed<<setprecision(6)<<allequations.vpot[i].vr[j]*AA_Bohr<<endl;
          }
          outfile.close();
        }
      }
    }
    endgen=true;
    MyObjective(genome); 
  }

  allequations.reset();
  allequations.score();
  if(runtest){
    cout<<"final_score: \n"<<sqrt(allequations.restot2/(allequations.nrows-allequations.nspleq))<<endl; cout.flush();
  }
  allequations.writeout();

  cout<<"\n\nrepopt  total run time:  " <<gtime.elapsed() << " seconds to run.\n";
  cout << "** repopt normal termination **" << endl << endl;

  return 0;
}

double MyObjective(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  const double zero = 0.000000000000001;
  bool   checknan,checkexponential,checkpotential,checkderivative1st,checkderivative2nd,checkderivative3rd,checkderivative4th;
  int    i,j,k,nknots,idx,icoeff,ipoint,npoint,ntmp,naa0,naa1,naa2,naa3,naa4,icoeff0;
  double maxderivative2nd,addscore,ascale;
  double potential0,derivative1st0,derivative2nd0,derivative3rd0,derivative4th0;
  double potentialt,derivative1stt,derivative2ndt,derivative3rdt,derivative4tht;
  double score,potential,derivative1st,derivative2nd,derivative3rd,derivative4th,coeff0,coeff1st,coeff2nd,coeff3rd,coeff4th,r,tmp;

  checknan=checkexponential=checkpotential=checkderivative1st=checkderivative2nd=checkderivative3rd=checkderivative4th=true;

  idx=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots;j++){ 
      allequations.vpot[i].vr[j]=genome.gene(idx);
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      ddh.d3[i].value=genome.gene(idx);
      idx++;
    }
  }
  if(ddh.tdamph){
    ddh.damph.value=genome.gene(idx);
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      ddh.hubbardderivs[i].value=genome.gene(idx);
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      ddh.damphver2[i].value=genome.gene(idx);
      idx++;
    }
  }

  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots-1;j++){
      for(k=j+1;k<allequations.vpot[i].nknots-1;k++){
        if(allequations.vpot[i].vr[k]<allequations.vpot[i].vr[j]){
          tmp=allequations.vpot[i].vr[j];
          allequations.vpot[i].vr[j]=allequations.vpot[i].vr[k];
          allequations.vpot[i].vr[k]=tmp;
        }
      }
    }
  }

  for(i=0;i<allequations.vpot.size();i++){ 
    nknots=allequations.vpot[i].nknots;
    allequations.vpot[i].vr[0] = GAMin(svpot[i].minRbond, allequations.vpot[i].vr[0]);
    allequations.vpot[i].vr[0] = GAMax(svpot[i].minr, allequations.vpot[i].vr[0]);

    for(j=1;j<nknots-1;j++){ 
      allequations.vpot[i].vr[j] = GAMin(svpot[i].vr[nknots-1]-svpot[i].max_step, allequations.vpot[i].vr[j]);
      allequations.vpot[i].vr[j] = GAMax(svpot[i].minRbond+svpot[i].max_step, allequations.vpot[i].vr[j]);
    }
  
    for(j=1;j<nknots-1;j++){ 
      for(k=j+1;k<nknots-1;k++){
        if(abs(allequations.vpot[i].vr[k]-allequations.vpot[i].vr[j])<svpot[i].max_step){
          allequations.vpot[i].vr[k]=GAMin(allequations.vpot[i].vr[j]+svpot[i].max_step,svpot[i].vr[nknots-1]-svpot[i].max_step);
        }
      }
    }
    for(j=nknots-1;j>0;j--){
      for(k=j-1;k>=0;k--){
        if(abs(allequations.vpot[i].vr[k]-allequations.vpot[i].vr[j])<svpot[i].max_step){
          allequations.vpot[i].vr[k]-=svpot[i].max_step;
        }
      }
    }

  }

  idx=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots;j++){ 
      genome.gene(idx,allequations.vpot[i].vr[j]);
      idx++;
    }
  }

  allequations.reset();
  allequations.score();

  if(score_type==0) score=allequations.restotS/(allequations.nrows-allequations.nspleq);
  else if(score_type==1) score=allequations.restotU/(allequations.nrows-allequations.nspleq);
  else if(score_type==2) score=pow(allequations.restot2/(allequations.nrows-allequations.nspleq),0.5);
  else if(score_type==4) score=pow(allequations.restot4/(allequations.nrows-allequations.nspleq),0.25);
  else if(score_type==8) score=pow(allequations.restot8/(allequations.nrows-allequations.nspleq),0.125);
  else score=pow(allequations.restot2/(allequations.nrows-allequations.nspleq),0.5);

  if(fsmooth_spline){ 
    icoeff=0;
    for(i=0; i<allequations.vpot.size(); i++){
      nknots=allequations.vpot[i].nknots;
      svpot[i].expA=allequations.vpot[i].expA;
      svpot[i].expB=allequations.vpot[i].expB;
      svpot[i].expC=allequations.vpot[i].expC;
      for(j=0; j<nknots-1; j++){
        svpot[i].division[j].p[0]=allequations.vpot[i].vr[j];
        svpot[i].division[j].p[1]=allequations.vpot[i].vr[j+1];
        for(k=0; k<=allequations.vpot[i].ordspl; k++){
          svpot[i].division[j].coeff[k]=allequations.vunknown[icoeff];
          icoeff++;
        }
      }
    }
  }

  ascale=0.0;
  icoeff0=0;
  addscore=0.0;
  maxderivative2nd=0.0; 
  for(i=0; i<allequations.vpot.size(); i++){
    naa0=naa1=naa2=naa3=naa4=svpot[i].aa;
    if(fsmooth_spline){ 
      maxderivative2nd=0.0; 
      coeff0=coeff1st=coeff2nd=coeff3rd=coeff4th=1.0;
      checknan=checkpotential=checkderivative1st=checkderivative2nd=checkderivative3rd=checkderivative4th=true;
      r=allequations.vpot[i].vr[0]+deltar;
      nknots=allequations.vpot[i].nknots;
      npoint=(allequations.vpot[i].vr[nknots-1]-allequations.vpot[i].vr[0])/deltar;

      potential0    =potentialt    =potential    = 1.0;
      derivative1st0=derivative1stt=derivative1st=-1.0;
      derivative2nd0=derivative2ndt=derivative2nd= 1.0;
      derivative3rd0=derivative3rdt=derivative3rd=-1.0;
      derivative4th0=derivative4tht=derivative4th= 1.0;

      for(ipoint=0;ipoint<npoint-5;ipoint++){
        r+=deltar;
        potential0    =potentialt;
        derivative1st0=derivative1stt;
        derivative2nd0=derivative2ndt;
        derivative3rd0=derivative3rdt;
        derivative4th0=derivative4tht;
        if(r<svpot[i].division[0].p[0]){
          potential=exp(-svpot[i].expA*r+svpot[i].expB)+svpot[i].expC;
          derivative1st=-svpot[i].expA*(potential-svpot[i].expC);
          derivative2nd=-svpot[i].expA*derivative1st;
          derivative3rd=-svpot[i].expA*derivative2nd;
          derivative4th=-svpot[i].expA*derivative3rd;
        }else{
          for(j=0;j<nknots-1;j++){
            if((r>=svpot[i].division[j].p[0])&&(r<svpot[i].division[j].p[1])) break;
          }
          potential=polynomial4th(svpot[i].division[j].p[0],r,svpot[i].division[j].coeff[0],svpot[i].division[j].coeff[1],svpot[i].division[j].coeff[2],svpot[i].division[j].coeff[3],svpot[i].division[j].coeff[4]);
          derivative1st=polynomial3rd(svpot[i].division[j].p[0],r,svpot[i].division[j].coeff[1],2.0*svpot[i].division[j].coeff[2],3.0*svpot[i].division[j].coeff[3],4.0*svpot[i].division[j].coeff[4]);
          derivative2nd=polynomial2nd(svpot[i].division[j].p[0],r,2.0*svpot[i].division[j].coeff[2],6.0*svpot[i].division[j].coeff[3],12.0*svpot[i].division[j].coeff[4]);
          derivative3rd=polynomial1st(svpot[i].division[j].p[0],r,6.0*svpot[i].division[j].coeff[3],24.0*svpot[i].division[j].coeff[4]);
          derivative4th=24.0*svpot[i].division[j].coeff[4];
        }

        if((abs(potential)    <zero) || (abs(potential0*potential)        <zero*zero)) continue;
        if((abs(derivative1st)<zero) || (abs(derivative1st0*derivative1st)<zero*zero)) continue;
        if((abs(derivative2nd)<zero) || (abs(derivative2nd0*derivative2nd)<zero*zero)) continue;
        if((abs(derivative3rd)<zero) || (abs(derivative3rd0*derivative3rd)<zero*zero)) continue;
        if((abs(derivative4th)<zero) || (abs(derivative4th0*derivative4th)<zero*zero)) continue;

        potentialt    =potential;
        derivative1stt=derivative1st;
        derivative2ndt=derivative2nd;
        derivative3rdt=derivative3rd;
        derivative4tht=derivative4th;

        if(abs(derivative2nd)>maxderivative2nd) maxderivative2nd = abs(derivative2nd); 

        if((naa0>0)&&(potential0*potential<0.0))         {coeff0  =-1.0*coeff0;   naa0=naa0-1;}
        if((naa1>0)&&(derivative1st0*derivative1st<0.0)) {coeff1st=-1.0*coeff1st; naa1=naa1-1;} 
        if((naa2>0)&&(derivative2nd0*derivative2nd<0.0)) {coeff2nd=-1.0*coeff2nd; naa3=naa2-1;} 
        if((naa3>0)&&(derivative3rd0*derivative3rd<0.0)) {coeff3rd=-1.0*coeff3rd; naa4=naa3-1;} 
        if((naa4>0)&&(derivative4th0*derivative4th<0.0)) {coeff4th=-1.0*coeff4th; naa4=naa4-1;} 

        potential=coeff0*potential;
        derivative1st=coeff1st*derivative1st;
        derivative2nd=coeff2nd*derivative2nd;
        derivative3rd=coeff3rd*derivative3rd;
        derivative4th=coeff4th*derivative4th;
        
        if(potential<0.0){ascale+=3.2; checkpotential=false;break;}
        if((svpot[i].smooth_order>=1) && (derivative1st>0.0)){ascale+=1.6; checkderivative1st=false;break;}
        if((svpot[i].smooth_order>=2) && (derivative2nd<0.0)){ascale+=0.8; checkderivative2nd=false;break;}
        if((svpot[i].smooth_order>=3) && (derivative3rd>0.0)){ascale+=0.4; checkderivative3rd=false;break;}
        if((svpot[i].smooth_order>=4) && (derivative4th<0.0)){ascale+=0.2; checkderivative4th=false;break;}
      }
    }
  
    if(allequations.vunknown[icoeff0+2] < 0.0) {
      checknan=false;
      ascale+=1000; 
      if(endgen) cout<<"#Warning NAN"<< allequations.vpot[i].potname<<":  "<<allequations.vunknown[icoeff0+2]<<"  "<<maxderivative2nd<<endl;
    }
    if(2.0*allequations.vunknown[icoeff0+2] < 0.999*maxderivative2nd) {
      checkexponential=false; 
      ascale+=100; 
      ascale += maxderivative2nd - 2.0*allequations.vunknown[icoeff0+2];
      if(endgen) cout<<"#Warning Exponent"<< allequations.vpot[i].potname<<":  "<<allequations.vunknown[icoeff0+2]<<"  "<<maxderivative2nd<<endl;
    }
    if(!checkpotential && endgen)     cout<<"#Warning PotE"<< allequations.vpot[i].potname<<endl;
    if(!checkderivative1st && endgen) cout<<"#Warning DerE1"<< allequations.vpot[i].potname<<endl;
    if(!checkderivative2nd && endgen) cout<<"#Warning DerE2"<< allequations.vpot[i].potname<<endl;
    if(!checkderivative3rd && endgen) cout<<"#Warning DerE3"<< allequations.vpot[i].potname<<endl;
    if(!checkderivative4th && endgen) cout<<"#Warning DerE4"<< allequations.vpot[i].potname<<endl;

    if(std::isnan(allequations.vpot[i].expA)||std::isnan(allequations.vpot[i].expB)||std::isnan(allequations.vpot[i].expC)){checknan=false;}
    if(std::isnan(allequations.vpot[i].expA)||std::isnan(allequations.vpot[i].expB)||std::isnan(allequations.vpot[i].expC)){checkexponential=false;}
    if(!checknan) ascale+=1000; 
    if(!checkexponential) ascale+=1000; 

    //if(std::isnan(allequations.vpot[i].expA)||std::isnan(allequations.vpot[i].expB)||std::isnan(allequations.vpot[i].expC)){checknan=false;}
    icoeff0 += (allequations.vpot[i].vr.size()-1)*(allequations.vpot[i].ordspl+1);
  }
  
  if(std::isnan(score)) score=999999999.9;

  return ((1.0+ascale)*score);
}

void MyInitializer(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  bool accept; 
  int i,j,k,idx,nknots,ntmp; 
  double tmp;
  if(read_spline){
    for(i=0;i<allequations.vpot.size();i++){ 
      nknots=svpot[i].nknots;
      for(j=0;j<allequations.vpot[i].nknots;j++){ 
        allequations.vpot[i].vr[j]=svpot[i].vr[j];
      }
      allequations.vpot[i].vr[0] = GAMax(svpot[i].minr, allequations.vpot[i].vr[0]);
      allequations.vpot[i].vr[0] = GAMin(svpot[i].minRbond, allequations.vpot[i].vr[0]);
      for(j=1;j<nknots-1;j++){ 
        allequations.vpot[i].vr[j] = GAMax(svpot[i].minRbond+svpot[i].max_step, allequations.vpot[i].vr[j]);
        allequations.vpot[i].vr[j] = GAMin(svpot[i].vr[nknots-1]-svpot[i].max_step, allequations.vpot[i].vr[j]);
      }
    }
    read_spline=false;
    //cout<<"read\n";
  }else{
    idx=0;
    for(i=0;i<allequations.vpot.size();i++){ 
      nknots=svpot[i].nknots;
      allequations.vpot[i].vr[0]=GARandomFloat(svpot[i].minr,svpot[i].minRbond);
      allequations.vpot[i].vr[nknots-1]=svpot[i].vr[nknots-1];
      for(j=1;j<nknots-1;j++){
        allequations.vpot[i].vr[j]=0.0;
      }
      for(j=1;j<nknots-1;j++){ 
        accept=false;
        while(!accept){
          accept=true;
          tmp=GARandomFloat(svpot[i].minRbond+svpot[i].max_step,svpot[i].vr[nknots-1]-svpot[i].max_step);
          for(k=1;k<nknots;k++){
            if(abs(allequations.vpot[i].vr[k]-tmp)<svpot[i].max_step){
              accept=false;
              break;
            }
          }
        }
        allequations.vpot[i].vr[j]=tmp;
      }
     
      for(j=1;j<nknots-1;j++){
        for(k=j+1;k<nknots-1;k++){
          if(allequations.vpot[i].vr[k]<allequations.vpot[i].vr[j]){
            tmp=allequations.vpot[i].vr[j];
            allequations.vpot[i].vr[j]=allequations.vpot[i].vr[k];
            allequations.vpot[i].vr[k]=tmp;
          }
        }
      }
    }
    if(ddh.td3){
      for(i=0; i<ddh.d3.size(); i++){
        tmp=GARandomFloat(ddh.d3[i].min,ddh.d3[i].max);
        ddh.d3[i].value = tmp;
      }
    }
    if(ddh.tdamph){
      tmp=GARandomFloat(ddh.damph.min,ddh.damph.max);
      ddh.damph.value = tmp;
    }
    if(ddh.thubbardderivs){
      for(i=0; i<ddh.hubbardderivs.size(); i++){
        tmp=GARandomFloat(ddh.hubbardderivs[i].min,ddh.hubbardderivs[i].max);
        ddh.hubbardderivs[i].value = tmp;
      }
    }
    if(ddh.tdamphver2){
      for(i=0; i<ddh.damphver2.size(); i++){
        tmp=GARandomFloat(ddh.damphver2[i].min,ddh.damphver2[i].max);
        ddh.damphver2[i].value = tmp;
      }
    }
  }
  idx=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots;j++){ 
      genome.gene(idx,allequations.vpot[i].vr[j]);
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
}

int MyMutator(GAGenome& g, double pmut){
  GA1DArrayGenome<double>& child = (GA1DArrayGenome<double>&)g;
  int i,j,k,idx,nknots,ntmp,iMut,tnMut;
  double fMut,value,tmp;

  //cout<<"in mutator start\n";
  if(pmut <= 0.0) return(0);
  //get data from gene
  idx=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots;j++){ 
      allequations.vpot[i].vr[j]=child.gene(idx);
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      ddh.d3[i].value=child.gene(idx);
      idx++;
    }
  }
  if(ddh.tdamph){
    ddh.damph.value=child.gene(idx);
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      ddh.hubbardderivs[i].value=child.gene(idx);
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      ddh.damphver2[i].value=child.gene(idx);
      idx++;
    }
  }
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots-1;j++){
      for(k=j+1;k<allequations.vpot[i].nknots-1;k++){
        if(allequations.vpot[i].vr[k]<allequations.vpot[i].vr[j]){
          tmp=allequations.vpot[i].vr[j];
          allequations.vpot[i].vr[j]=allequations.vpot[i].vr[k];
          allequations.vpot[i].vr[k]=tmp;
        }
      }
    }
  }

  //start to mutate data
  tnMut=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    nknots=svpot[i].nknots;
    fMut = pmut*(double)(nknots-1.0) ;

    if(fMut < 1.0){   // we have to do a flip test on each element
      iMut = 0;
      for(j=0; j<(nknots-1); j++){
        if(GAFlipCoin(pmut)){
          value  = GAGaussianFloat(svpot[i].max_step);
          value += allequations.vpot[i].vr[j];
          value = GAMax(svpot[i].minr, value);
          value = GAMin(svpot[i].vr[nknots-1]-svpot[i].max_step, value);
          allequations.vpot[i].vr[j] = value;
          iMut++;
        }
      }
    }else{   // only mutate the ones we need to
      iMut = (int)(fMut); 
      for(k=0; k<iMut; k++){
        j = GARandomInt(0,nknots-2);
        value  = GAGaussianFloat(svpot[i].max_step);
        value += allequations.vpot[i].vr[j];
        value = GAMax(svpot[i].minr, value);
        value = GAMin(svpot[i].vr[nknots-1]-svpot[i].max_step, value);
        allequations.vpot[i].vr[j] = value;
      }
    }
    tnMut+=iMut;
  }
  if(ddh.td3){
    ntmp=ddh.d3.size();
    fMut = pmut*(double)(ntmp) ;
    if(fMut < 1.0){
      iMut = 0;
      for(i=0; i<ntmp; i++){
        if(GAFlipCoin(pmut)){
          value  = GAGaussianFloat(gauss_dev);
          value += ddh.d3[i].value;
          value = GAMax(ddh.d3[i].min, value);
          value = GAMin(ddh.d3[i].max, value);
          ddh.d3[i].value = value;
          iMut++;
        }
      }
    }else{
      iMut = (int)(fMut); 
      for(k=0; k<iMut; k++){
        i = GARandomInt(0,ntmp-1);
        value  = GAGaussianFloat(gauss_dev);
        value += ddh.d3[i].value;
        value = GAMax(ddh.d3[i].min, value);
        value = GAMin(ddh.d3[i].max, value);
        ddh.d3[i].value = value;
      }
    }
    tnMut+=iMut;
  }
  if(ddh.tdamph){
    if(GAFlipCoin(pmut)){
      value  = GAGaussianFloat(gauss_dev);
      value += ddh.damph.value;
      value = GAMax(ddh.damph.min, value);
      value = GAMin(ddh.damph.max, value);
      ddh.damph.value = value;
      tnMut++;
    }
  }
  if(ddh.thubbardderivs){
    ntmp=ddh.hubbardderivs.size();
    fMut = pmut*(double)(ntmp) ;
    if(fMut < 1.0){
      iMut = 0;
      for(i=0; i<ntmp; i++){
        if(GAFlipCoin(pmut)){
          value  = GAGaussianFloat(gauss_dev);
          value += ddh.hubbardderivs[i].value;
          value = GAMax(ddh.hubbardderivs[i].min, value);
          value = GAMin(ddh.hubbardderivs[i].max, value);
          ddh.hubbardderivs[i].value = value;
          iMut++;
        }
      }
    }else{
      iMut = (int)(fMut); 
      for(k=0; k<iMut; k++){
        i = GARandomInt(0,ntmp-1);
        value  = GAGaussianFloat(gauss_dev);
        value += ddh.hubbardderivs[i].value;
        value = GAMax(ddh.hubbardderivs[i].min, value);
        value = GAMin(ddh.hubbardderivs[i].max, value);
        ddh.hubbardderivs[i].value = value;
      }
    }
    tnMut+=iMut;
  }
  if(ddh.tdamphver2){
    ntmp=ddh.damphver2.size();
    fMut = pmut*(double)(ntmp) ;
    if(fMut < 1.0){
      iMut = 0;
      for(i=0; i<ntmp; i++){
        if(GAFlipCoin(pmut)){
          value  = GAGaussianFloat(gauss_dev);
          value += ddh.damphver2[i].value;
          value = GAMax(ddh.damphver2[i].min, value);
          value = GAMin(ddh.damphver2[i].max, value);
          ddh.damphver2[i].value = value;
          iMut++;
        }
      }
    }else{
      iMut = (int)(fMut); 
      for(k=0; k<iMut; k++){
        i = GARandomInt(0,ntmp-1);
        value  = GAGaussianFloat(gauss_dev);
        value += ddh.damphver2[i].value;
        value = GAMax(ddh.damphver2[i].min, value);
        value = GAMin(ddh.damphver2[i].max, value);
        ddh.damphver2[i].value = value;
      }
    }
    tnMut+=iMut;
  }
 
  //convert from data to gene
  idx=0;
  for(i=0;i<allequations.vpot.size();i++){ 
    for(j=0;j<allequations.vpot[i].nknots;j++){ 
      child.gene(idx,allequations.vpot[i].vr[j]);
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      child.gene(idx,ddh.d3[i].value);
      idx++;
    }
  }
  if(ddh.tdamph){
    child.gene(idx,ddh.damph.value);
    idx++;
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      child.gene(idx,ddh.hubbardderivs[i].value);
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      child.gene(idx,ddh.damphver2[i].value);
      idx++;
    }
  }
 
  return(tnMut);
}


