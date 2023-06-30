#ifndef LIFDELTA_H
#define LIFDELTA_H

#include <cstdlib>
#include <iostream>

class LIFDelta //Single LIF neuron with delta coupling
{
public:
  LIFDelta();
  LIFDelta(double gamma_, double th_,double Vinf_);
  double pseudo_spike_time(double t0_,double tE_);
  void reset();
  double Vvalue(double t0_,double tE_);
  
  double V0;
  
private:
  double gamma, th; // 1/time, Threshold
  double Vinf; //Potential for t->inf
  //double tsp; //Pseudospikezeit, koennte man sich merken und nur updaten, wenn ein Input gekommen ist. (Sparse implementierung.)
  
  
  friend std::ostream& operator<<(std::ostream& o , LIFDelta const & neuron);
};

#endif // CURRENBASEDLIFDELTA_H
