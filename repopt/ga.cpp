#include <cstdio>
#include <cmath>
#include <limits>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "ga.hpp"

using namespace std;

#define cout STD_COUT
#define cerr STD_CERR
#define endl STD_ENDL
#define ofstream STD_OFSTREAM
#if !defined(GALIB_USE_AUTO_INST)
#include <ga/GA1DArrayGenome.C>
GALIB_INSTANTIATION_PREFIX GA1DArrayGenome<double>;
#endif

extern int preserved_num;
extern int destroy_num;
extern int popsizemin;

void
MyGASimpleGA::step() {
  int i, nsize, mut, c1, c2;
  GAGenome *mom, *dad;         // tmp holders for selected genomes

  GAPopulation *tmppop;        // Swap the old population with the new pop.
  oldPop->size(pop->size());   // get info about current pop for next time
  tmppop = oldPop;             // When we finish the ++ we want the newly 
  oldPop = pop;                // generated population to be current (for
  pop    = tmppop;             // references to it from member functions).

  // Generate the individuals in the temporary population from individuals in 
  // the main population.

  for(i=0; i<pop->size()-1; i+=2){  // takes care of odd population
    mom = &(oldPop->select());  
    dad = &(oldPop->select());
    stats.numsel += 2;              // keep track of number of selections
    c1 = c2 = 0;
    if(GAFlipCoin(pCrossover())){
      stats.numcro += (*scross)(*mom, *dad, &pop->individual(i), &pop->individual(i+1));
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
  if(pop->size() % 2 != 0){      // do the remaining population member
    mom = &(oldPop->select());  
    dad = &(oldPop->select());
    stats.numsel += 2;           // keep track of number of selections
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

  stats.numrep += pop->size();
  pop->evaluate(gaTrue);         // get info about current pop for next time

  pop->sort(gaTrue);             //
  if(pop->size()-destroy_num>popsizemin) {
    nsize=pop->size()-destroy_num;
  }else nsize=popsizemin;
  pop->size(nsize); 

  pop->sort(gaTrue);             //
  if(preserved_num>pop->size()) preserved_num=pop->size();
  if(minimaxi() == GAGeneticAlgorithm::MAXIMIZE) {
    for(i=0;i<preserved_num;i++) {
      if(oldPop->best(i).score() > pop->worst(i).score()) {
        pop->worst(i).copy(oldPop->best(i));
        stats.numrep += 1;
      } else break;
    }
  }
  else {
    for(i=0;i<preserved_num;i++) {
      if(oldPop->best(i).score() < pop->worst(i).score()) {
        pop->worst(i).copy(oldPop->best(i));
        stats.numrep += 1;
      } else break;
    }
  }
  pop->sort(gaTrue);    // get info about current pop for next time
  stats.update(*pop);   // update the statistics by one generation
}

double MyComparator(const GAGenome& g1, const GAGenome& g2) {
  GA1DArrayGenome<double>& gnome1 = (GA1DArrayGenome<double>&)g1;
  GA1DArrayGenome<double>& gnome2 = (GA1DArrayGenome<double>&)g2;
  double dist=0.0; 
  int i;
  for(i=0;i<gnome1.length();i++){ 
    dist+=abs(gnome2.gene(i)-gnome1.gene(i));
  }
  return dist;
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
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0]  = momsite[1] - momsite[0];
      momlen[1]  = mom.length() - momsite[1];
      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0]  = momlen[0];
      dadlen[1]  = momlen[1];
    }else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      return nc;
    }else{
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0]  = momsite[1] - momsite[0];
      momlen[1]  = mom.length() - momsite[1];
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
      momlen[0]  = momsite[1] - momsite[0];
      momlen[1]  = mom.length() - momsite[1];
      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0]  = momlen[0];
      dadlen[1]  = momlen[1];
    }else{
      momsite[0] = 3*GARandomInt(0, mom.length()/3);
      momsite[1] = 3*GARandomInt(0, mom.length()/3);
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0]  = momsite[1] - momsite[0];
      momlen[1]  = mom.length() - momsite[1];

      dadsite[0] = 3*GARandomInt(0, dad.length()/3);
      dadsite[1] = 3*GARandomInt(0, dad.length()/3);
      if(dadsite[0] > dadsite[1]) SWAP(dadsite[0], dadsite[1]);
      dadlen[0]  = dadsite[1] - dadsite[0];
      dadlen[1]  = dad.length() - dadsite[1];
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

