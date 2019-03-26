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
void MyGASimpleGA::initialstep(){
  int i,j,k,length,idx,impi, start, end, size, rank,irank;
  double buffer;
  double * abuffer;
  length=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    length+=erepobj.velem[i].lmax+2;
  }
  if(ddh.td3) length+=ddh.d3.size(); 
  if(ddh.tdamph) length+=1;
  if(ddh.thubbardderivs) length+=ddh.hubbardderivs.size(); 
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
  double buffer;
  double * abuffer;
  length=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    length+=erepobj.velem[i].lmax+2;
  }
  if(ddh.td3) length+=ddh.d3.size(); 
  if(ddh.tdamph) length+=1;
  if(ddh.thubbardderivs) length+=ddh.hubbardderivs.size(); 
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


