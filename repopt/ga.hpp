#ifndef GA_HPP_INCLUDED
#define GA_HPP_INCLUDED

class MyGASimpleGA : public GASimpleGA {
public:
  GADefineIdentity("MyGASimpleGA", 200);
  MyGASimpleGA(const GAGenome& g) : GASimpleGA(g) {}
  virtual ~MyGASimpleGA() {}
  virtual void step();
  MyGASimpleGA & operator++() { step(); return *this; }
};

int    MyTwoPointCrossover(const GAGenome&,const GAGenome&,GAGenome*,GAGenome*);
double MyComparator(const GAGenome&,const GAGenome&);

#endif /* GA_HPP_INCLUDED */
