#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include "spike.h"


Spike::Spike()
{}

Spike::Spike(double time_, int sender_) : time(time_), sender(sender_)
{}

Spike::Spike(double time_, int sender_, double dummy_) : time(time_), sender(sender_)
{}

std::ostream& operator<<(std::ostream& o , Spike const & sp)
{
    o << std::setprecision(18) << sp.time << ',' << sp.sender;
    return o;
}

int operator<(Spike const & s1 , Spike const & s2) {
    return (s1.time < s2.time);
}
int operator==(Spike const & s1 , Spike const & s2) {
    return (s1.time == s2.time);
}



DesSpike::DesSpike()
{}

DesSpike::DesSpike(double time_, int sender_, double eps_)
        : Spike(time_, sender_), eps(eps_)
{}

std::ostream& operator<<(std::ostream& o , DesSpike const & dsp)
{
    o<<"Spike: time="<<dsp.time << ", sender=" << dsp.sender << ", tol.="<<dsp.eps;
    return o;
}
