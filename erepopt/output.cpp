#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include"erepobj.hpp"
#include"ddh.hpp"

//const double AA_Bohr = 0.5291772085936 ;

using namespace std;

extern Erepobj erepobj;
extern vector<string> addHamiltonian;
extern int ga_popsize, ga_ngen, ga_scoref, ga_flushf;
extern int preserved_num, seed;
extern double gtol;
extern double ga_pmut, ga_pcross;
extern bool ga,readr;
extern bool runtest;
extern sddh ddh;


void Erepobj::writeout(){

  fout<<std::fixed;
  restotU=0; restotS=0; restot2=0; restot4=0;
  resevU=0; resevS=0; resev2=0; resev4=0;
  rescU=0; rescS=0; resc2=0; resc4=0;


  fout << "** erepobj **" << endl << endl;
  fout << "interpreted input: " << endl << endl;
  fout << "dftb-version used: " << endl;
  fout << dftbversion << endl;
  fout << endl;
  fout << "genetic_algrithm:# " << endl;
  fout << "popsize"<<"  "<<ga_popsize<<endl;
  fout << "ngen"<<"  "<<ga_ngen<<endl;
  fout << "preserved_num"<<"  "<<preserved_num<<endl;
  fout << "scoref"<<"  "<<ga_scoref<<endl;
  fout << "flushf"<<"  "<<ga_flushf<<endl;
  fout << "score_type"<<"  "<<score_type<<endl;
  fout << "seed"<<"  "<<seed<<endl;
  fout << "pmut"<<"  "<<ga_pmut<<endl;
  fout << "pcross"<<"  "<<ga_pcross<<endl;
  fout << "tol"<<"  "<<gtol<<endl;
  fout << "ga"<<"  "<<ga<<endl;
  fout << "runtest"<<"  "<<runtest<<endl;
  fout << "readr"<<"  "<<readr<<endl;
  fout << endl;
 

  fout<<"OPTIMAL RADII:"<<endl;
  int i,k;
  fout<<"$element_types:\n"; 
  for(i=0;i<velem.size();i++){
    fout<<std::fixed<<std::left;
    fout.width(3); fout<<velem[i].name;
    fout<<std::fixed<<std::right;
    fout.width(3); fout<<velem[i].optype;
    fout.width(3); fout<<velem[i].lmax;  
    for(k=0;k<velem[i].lmax+2.;k++){
      fout.precision(velem[i].radius[k].precision); 
    //fout.width(9); fout<<velem[i].radius[k].minr;
      fout.width(9); fout<<velem[i].radius[k].r;
      fout.width(5); fout<<velem[i].radius[k].r;
      fout.width(5); fout<<velem[i].radius[k].r;
    //fout.width(9); fout<<velem[i].radius[k].maxr;
      fout.width(3); fout<<velem[i].radius[k].precision;
    }
    fout<<endl;
  }
  fout<<"$end\n"; 
    
  if(ddh.td3){
    fout<<"$d3:\n"; 
    for(int i=0;i<ddh.d3.size();i++){ 
      fout.precision(ddh.d3[i].precision) ;
    //fout<<ddh.d3[i].name<<" "<<ddh.d3[i].min<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].max<<" "<<ddh.d3[i].precision<<endl;
      fout<<ddh.d3[i].name<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].precision<<endl;
    }
    fout<<"$end\n"; 
  }
  if(ddh.tdamph){
    fout<<"$damph:\n"; 
    fout.precision(ddh.damph.precision);
    //fout<<ddh.damph.min<<" "<<ddh.damph.value<<" "<<ddh.damph.max<<" "<<ddh.damph.precision<<endl;
      fout<<ddh.damph.value<<" "<<ddh.damph.value<<" "<<ddh.damph.value<<" "<<ddh.damph.precision<<endl;
    fout<<"$end\n"; 
  }
  if(ddh.thubbardderivs){
    fout<<"$hubbardderivs:\n"; 
    for(int i=0;i<ddh.hubbardderivs.size();i++){ 
      fout.precision(ddh.hubbardderivs[i].precision);
    //fout<<ddh.hubbardderivs[i].name<<" "<<ddh.hubbardderivs[i].min<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].max<<" "<<ddh.hubbardderivs[i].precision<<endl;
      fout<<ddh.hubbardderivs[i].name<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].precision<<endl;
    }
    fout<<"$end\n"; 
  }
  if(ddh.tdamphver2){
    fout<<"$damphver2:\n"; 
    for(int i=0;i<ddh.damphver2.size();i++){ 
      fout.precision(ddh.damphver2[i].precision);
    //fout<<ddh.damphver2[i].name<<" "<<ddh.damphver2[i].min<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].max<<" "<<ddh.damphver2[i].precision<<endl;
      fout<<ddh.damphver2[i].name<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].precision<<endl;
    }
    fout<<"$end\n"; 
  }
  if(ddh.thubbards){
    fout<<"$hubbards:\n"; 
    for(int i=0;i<ddh.hubbards.size();i++){ 
      fout.precision(ddh.hubbards[i].precision);
    //fout<<ddh.hubbards[i].name<<" "<<ddh.hubbards[i].type<<" "<<ddh.hubbards[i].min<<" "<<ddh.hubbards[i].value<<" "<<ddh.hubbards[i].max<<" "<<ddh.hubbards[i].precision<<endl;
      fout<<ddh.hubbards[i].name<<" "<<ddh.hubbards[i].type<<" "<<ddh.hubbards[i].value<<" "<<ddh.hubbards[i].value<<" "<<ddh.hubbards[i].value<<" "<<ddh.hubbards[i].precision<<endl;
    }
    fout<<"$end\n"; 
  }
  if(ddh.tvorbes){
    fout<<"$vorbes:\n"; 
    for(int i=0;i<ddh.vorbes.size();i++){ 
      fout.precision(ddh.vorbes[i].precision);
    //fout<<ddh.vorbes[i].name<<" "<<ddh.vorbes[i].type<<" "<<ddh.vorbes[i].min<<" "<<ddh.vorbes[i].value<<" "<<ddh.vorbes[i].max<<" "<<ddh.vorbes[i].precision<<endl;
      fout<<ddh.vorbes[i].name<<" "<<ddh.vorbes[i].type<<" "<<ddh.vorbes[i].value<<" "<<ddh.vorbes[i].value<<" "<<ddh.vorbes[i].value<<" "<<ddh.vorbes[i].precision<<endl;
    }
    fout<<"$end\n"; 
  }

  double dnatom,dnev,delta,adelta,delta2,delta4; 
  dnatom=dnev=0.0;

  fout<<"evscaling:  "<<evscaling<<endl;
  for(int i=0;i<vmol.size();i++){
    fout<<std::fixed<<std::left;
    fout.width(80);
    fout<<vmol[i].name;
    fout<<std::fixed<<std::right;
    fout.width(20);
    fout.precision(6);
    fout<<vmol[i].b0<<endl;
  }
  fout << endl;

  fout<<"ERROR:"<<endl;
  for(int i=0;i<vmol.size();i++){
    fout<<std::fixed<<std::left;
    fout.width(80);
    fout<<vmol[i].name;
    fout<<std::fixed<<std::right;
    fout.width(20);
    fout.precision(6);
    fout<<vmol[i].error<<endl;
  }
  fout << endl;
//for(int i=0;i<vmol.size();i++){
//  dnev+=vmol[i].nev;
//  fout<<vmol[i].name<<endl;
//  for(int j=0;j<vmol[i].nev;j++){
//    delta=vmol[i].evnew[j]-vmol[i].evref[j];
//    fout<<vmol[i].evref[j]<<"  "<<vmol[i].evnew[j]<<"  "<<delta<<endl;
//    adelta=abs(delta);
//    delta2=delta*delta;
//    delta4=delta2*delta2;
//    resevS+=delta;
//    resevU+=adelta;
//    resev2+=delta2;
//    resev4+=delta4;
//  }
//}

//restotS=(rescS+resevS)/(dnatom+dnev);
//restotU=(rescU+resevU)/(dnatom+dnev);
//restot2=(resc2+resev2)/(dnatom+dnev);
//restot4=(resc4+resev4)/(dnatom+dnev);
//rescS/=dnatom; 
//rescU/=dnatom; 
//resc2/=dnatom; 
//resc4/=dnatom; 
//resevS/=dnev; 
//resevU/=dnev; 
//resev2/=dnev; 
//resev4/=dnev; 

//fout<<endl;
//fout<<endl;
//fout.precision(6);
//
//fout<<"charge: "<<rescS<<" "<<rescU<<" "<<pow(resc2,0.5)<<" "<<pow(resc4,0.25)<<endl;
//fout<<"ev:     "<<resevS<<" "<<resevU<<" "<<pow(resev2,0.5)<<" "<<pow(resev4,0.25)<<endl;
//fout<<"total:  "<<restotS<<" "<<restotU<<" "<<pow(restot2,0.5)<<" "<<pow(restot4,0.25)<<endl;

}


