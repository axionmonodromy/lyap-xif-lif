
#ifndef SPIKE_H
#define SPIKE_H

#include <cstdlib>
#include <iostream>
#include <iomanip>

class Spike
{
public:
    double time;
    int sender;

    Spike();

    Spike(double time_, int sender_);
    
    Spike(double time_, int sender_, double dummy_);

    friend std::ostream& operator<<(std::ostream& o , Spike const & sp);
    friend int operator<(Spike const & s1 , Spike const & s2);
    friend int operator==(Spike const & s1 , Spike const & s2);
};


class DesSpike : public Spike
{

public:
    double eps;

    DesSpike();
    DesSpike(double time_, int sender_, double eps_=0.);


    friend std::ostream& operator<<(std::ostream& o , DesSpike const & sp);
};




#endif
