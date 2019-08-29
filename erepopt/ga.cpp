#include <cstdio>
#include <new>
#include <cmath>
#include <limits>
#include <omp.h>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "erepobj.hpp"
#include "ddh.hpp"
#include "mpi.h"

using namespace std;

#define cout STD_COUT
#define cerr STD_CERR
#define endl STD_ENDL
#define ofstream STD_OFSTREAM
#if !defined(GALIB_USE_AUTO_INST)
#include <ga/GA1DArrayGenome.C>
GALIB_INSTANTIATION_PREFIX GA1DArrayGenome<double>;
#endif

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#include "ga.hpp"
extern int cpu_number;
extern int preserved_num;
extern Erepobj erepobj;
extern sddh ddh;
extern int idrefr;
extern double s1,s2,s3,s4,s5,s6,s7,s8,s9,sdenswf;
void MyGASimpleGA::initialstep(){
  int i,j,k,length,idx,impi, start, end, size, rank,irank;
  int ie1,ie2,nvars;
  double value,tmpsref,tmps,tmpp;
  double buffer;
  double * abuffer;
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

  abuffer = new double [length];
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  start = rank*(pop->size()/size);
  if(rank==(size-1)){
    end = pop->size();
  }else{
    end = start + (pop->size()/size);
  }

if(rank==0){
  for(k=0; k<pop->size(); k++){
    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(k);
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
      tmpsref=genome.gene(idrefr);
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
      }else if(erepobj.velem[i].optype==115){
        tmps=double(int(tmps*10.0))/10.0;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1151){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s1*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1152){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s2*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1153){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s3*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1154){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s4*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==111){
        genome.gene(ie1-1,tmps);
        genome.gene(ie1+1,tmps);
      }else if(erepobj.velem[i].optype==1115){
        tmps=double(int(tmps*10.0))/10.0;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11151){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s1*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11152){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s2*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11153){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s3*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11154){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s4*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }
      ie1=ie1+erepobj.velem[i].lmax; 
    }
  }
}
///////////////////////////////////////////////////////////// 
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=0; i<pop->size(); i++){
    if(rank == 0 && i>=end) {
      GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
      for(j=0;j<length;j++){
        abuffer[j]=genome.gene(j);
      }
      irank=i/(pop->size()/size);
      MPI_Send(abuffer, length, MPI_DOUBLE, irank, 11111+(i+1), MPI_COMM_WORLD);
    }else if(rank != 0 && i>=start && i<end){
      MPI_Recv(abuffer, length, MPI_DOUBLE, 0, 11111+(i+1), MPI_COMM_WORLD, &status);
      GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
      for(j=0;j<length;j++){
        genome.gene(j,abuffer[j]);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
///////////////////////////////////////////////////////////// 
  for(impi=start; impi< end; impi++){
    erepobj.makeskf(pop->individual(impi));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(impi=start; impi< end; impi++){
    pop->individual(impi).evaluate(gaTrue);
  }
///////////////////////////////////////////////////////////// 
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0; i<pop->size(); i++){
    if(rank != 0 && i>=start && i<end){
      buffer=pop->individual(i).score();
      MPI_Send(&buffer, 1, MPI_DOUBLE, 0, 22222+(i+1), MPI_COMM_WORLD);
//    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
//    for(j=0;j<length;j++){
//      abuffer[j]=genome.gene(j);
//    }
//    MPI_Send(abuffer, length, MPI_DOUBLE, 0, 333333*(i+1), MPI_COMM_WORLD);
    }else if(rank == 0 && i>=end) {
      irank=i/(pop->size()/size);
      MPI_Recv(&buffer, 1, MPI_DOUBLE, irank, 22222+(i+1), MPI_COMM_WORLD, &status);
      pop->individual(i).score(buffer);
//    MPI_Recv(abuffer, length, MPI_DOUBLE, irank, 333333*(i+1), MPI_COMM_WORLD, &status);
//    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
//    for(j=0;j<length;j++){
//      genome.gene(j,abuffer[j]);
//    }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
///////////////////////////////////////////////////////////// 
//pop->evaluate(gaTrue);  // get info about current pop for next time
  if(rank==0){
    pop->sort(gaTrue);
  }
  stats.update(*pop);   // update the statistics by one generation
  delete abuffer;
}

void MyGASimpleGA::step() {
 // cout<<"step -1\n";
  int mut, c1, c2;
  int i,j,k,length,idx,impi, start, end, size, rank,irank;
  int ie1,ie2,nvars;
  double value,tmpsref,tmps,tmpp;
  double buffer;
  double * abuffer;
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

  abuffer = new double [length];
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
  // cout<<"step -1\n";
if(rank==0){

  GAGenome *mom, *dad;          // tmp holders for selected genomes

  GAPopulation *tmppop;   // Swap the old population with the new pop.
  tmppop = oldPop;    // When we finish the ++ we want the newly 
  oldPop = pop;     // generated population to be current (for
  pop = tmppop;     // references to it from member functions).

// Generate the individuals in the temporary population from individuals in 
// the main population.

 // cout<<"step 0\n";
  for(i=0; i<pop->size()-1; i+=2){  // takes care of odd population
    mom = &(oldPop->select());  
    dad = &(oldPop->select());
    stats.numsel += 2;    // keep track of number of selections

    c1 = c2 = 0;
    if(GAFlipCoin(pCrossover())){
      stats.numcro += (*scross)(*mom, *dad,
        &pop->individual(i), &pop->individual(i+1));
      c1 = c2 = 1;
    }
    else{
      pop->individual( i ).copy(*mom);
      pop->individual(i+1).copy(*dad);
    }
    stats.nummut += (mut = pop->individual( i ).mutate(pMutation()));
    if(mut > 0) c1 = 1;
    stats.nummut += (mut = pop->individual(i+1).mutate(pMutation()));
    if(mut > 0) c2 = 1;

    stats.numeval += c1 + c2;
  }
 // cout<<"step 1\n";
  if(pop->size() % 2 != 0){ // do the remaining population member
    mom = &(oldPop->select());  
    dad = &(oldPop->select());
    stats.numsel += 2;    // keep track of number of selections

    c1 = 0;
    if(GAFlipCoin(pCrossover())){
      stats.numcro += (*scross)(*mom, *dad, &pop->individual(i), (GAGenome*)0);
      c1 = 1;
    }
    else{
      if(GARandomBit()) pop->individual( i ).copy(*mom);
      else pop->individual( i ).copy(*dad);
    }
    stats.nummut += (mut = pop->individual( i ).mutate(pMutation()));
    if(mut > 0) c1 = 1;

    stats.numeval += c1;
  }
 // cout<<"step 2\n";

  stats.numrep += pop->size();

  for(k=0; k<pop->size(); k++){
    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(k);
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
      tmpsref=genome.gene(idrefr);
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
      }else if(erepobj.velem[i].optype==115){
        tmps=double(int(tmps*10.0))/10.0;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1151){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s1*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1152){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s2*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1153){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s3*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==1154){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s4*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==111){
        genome.gene(ie1-1,tmps);
        genome.gene(ie1+1,tmps);
      }else if(erepobj.velem[i].optype==1115){
        tmps=double(int(tmps*10.0))/10.0;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11151){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s1*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11152){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s2*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11153){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s3*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }else if(erepobj.velem[i].optype==11154){
        tmpsref=double(int(tmpsref*10.0))/10.0;
        genome.gene(idrefr,tmpsref);
        tmps=s4*tmpsref;
        genome.gene(ie1-1,s5*tmps);
        genome.gene(ie1+1,tmps);
        genome.gene(ie1,tmps);
      }
      ie1=ie1+erepobj.velem[i].lmax; 
    }
  }
}

  start = rank*(pop->size()/size);
  if(rank==(size-1)){
    end = pop->size();
  }else{
    end = start + (pop->size()/size);
  }
///////////////////////////////////////////////////////////// 
  MPI_Barrier(MPI_COMM_WORLD);
  for(i=0; i<pop->size(); i++){
    if(rank == 0 && i>=end) {
      GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
      for(j=0;j<length;j++){
        abuffer[j]=genome.gene(j);
      }
      irank=i/(pop->size()/size);
      MPI_Send(abuffer, length, MPI_DOUBLE, irank, 11111+(i+1), MPI_COMM_WORLD);
    }else if(rank != 0 && i>=start && i<end){
      MPI_Recv(abuffer, length, MPI_DOUBLE, 0, 11111+(i+1), MPI_COMM_WORLD, &status);
      idx=0;
      GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
      for(j=0;j<length;j++){
        genome.gene(j,abuffer[j]);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
///////////////////////////////////////////////////////////// 
  for(impi=start; impi< end; impi++){
    erepobj.makeskf(pop->individual(impi));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for(impi=start; impi< end; impi++){
    pop->individual(impi).evaluate(gaTrue);
  }
///////////////////////////////////////////////////////////// 
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0; i<pop->size(); i++){
    if(rank != 0 && i>=start && i<end){
      buffer=pop->individual(i).score();
      MPI_Send(&buffer, 1, MPI_DOUBLE, 0, 22222+(i+1), MPI_COMM_WORLD);
//    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
//    for(j=0;j<length;j++){
//      abuffer[j]=genome.gene(j);
//    }
//    MPI_Send(abuffer, length, MPI_DOUBLE, 0, 333333*(i+1), MPI_COMM_WORLD);
    }else if(rank == 0 && i>=end) {
      irank=i/(pop->size()/size);
      MPI_Recv(&buffer, 1, MPI_DOUBLE, irank, 22222+(i+1), MPI_COMM_WORLD, &status);
      pop->individual(i).score(buffer);
//    MPI_Recv(abuffer, length, MPI_DOUBLE, irank, 333333*(i+1), MPI_COMM_WORLD, &status);
//    GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)pop->individual(i);
//    for(j=0;j<length;j++){
//      genome.gene(j,abuffer[j]);
//    }
//    pop->individual(i).copy(genome);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
///////////////////////////////////////////////////////////// 

//pop->evaluate(gaTrue);  // get info about current pop for next time
//cout<<"step 3\n";

// If we are supposed to be elitist, carry the best individual from the old
// population into the current population.  Be sure to check whether we are
// supposed to minimize or maximize.

  if(rank==0){
    pop->sort(gaTrue);
    if(minimaxi() == GAGeneticAlgorithm::MAXIMIZE) {
      //if(el && oldPop->best().score() > pop->best().score()) oldPop->replace(pop->replace(&(oldPop->best()), GAPopulation::WORST), GAPopulation::BEST);
      for(i=0;i<preserved_num;i++) {
        if(oldPop->best(i).score() > pop->worst(i).score()) {
          pop->worst(i).copy(oldPop->best(i));
          stats.numrep += 1;
        } else break;
      }
    } else {
      //if(el && oldPop->best().score() < pop->best().score()) oldPop->replace(pop->replace(&(oldPop->best()), GAPopulation::WORST), GAPopulation::BEST);
      for(i=0;i<preserved_num;i++) {
        if(oldPop->best(i).score() < pop->worst(i).score()) {
          pop->worst(i).copy(oldPop->best(i));
          stats.numrep += 1;
        } else break;
      }
    }
    pop->sort(gaTrue);
  }
  stats.update(*pop);   // update the statistics by one generation
  delete abuffer;
}


