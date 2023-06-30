#include "lifdelta.h"

#include <math.h>

using namespace std;

LIFDelta::LIFDelta()
{
}


LIFDelta::LIFDelta(double gamma_, double th_,double Vinf_)
        : gamma(gamma_),th(th_), Vinf(Vinf_)
{
 //   cout << "gamma=" << gamma << ", th=" << th << ", Vinf=" << Vinf << endl;
}


void LIFDelta::reset()
{

    V0=0;

}


double LIFDelta::pseudo_spike_time(double t0_,double tE_)
{

//    cout << "V0=" << V0 << ", th=" << th << endl;

    if (V0>=th) return(t0_);

//    if (Vinf>th) //(V0<th) sicher erfuellt wenn wir bis hierhin kommen
//    {
        double tsp;
        tsp=t0_+1/gamma*log((Vinf-V0)/(Vinf-th));
        return(tsp);
//    }

//    return(tE_+1);//(V0<th)&&(Vinf<=th) sicher erfuellt wenn wir bis hierhin kommen

}


double LIFDelta::Vvalue(double t0_,double tE_)
{

    double VE;
    
    VE=(V0-Vinf)*exp(-gamma*(tE_-t0_))+Vinf;

    return(VE);//(V0<th)&&(Vinf<=th) sicher erfuellt wenn wir bis hierhin kommen

}


ostream& operator<<(ostream& o , LIFDelta const & neuron)
{
    o<<"Current based LIF neuron: Vt0="<<neuron.V0 << ", gamma="<< neuron.gamma<< ", th="<< neuron.th<< ", Vinf="<< neuron.Vinf<< endl;
    return o;
}



