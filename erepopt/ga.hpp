#ifndef GA_H_INCLUDED
#define GA_H_INCLUDED

class MyGASimpleGA : public GASimpleGA {
public:
  GADefineIdentity("MyGASimpleGA", 200);
  MyGASimpleGA(const GAGenome& g) : GASimpleGA(g) {}
  
  virtual ~MyGASimpleGA() {}
  virtual void step();
  void initialstep();
  MyGASimpleGA & operator++() { step(); return *this; }
};

#endif
