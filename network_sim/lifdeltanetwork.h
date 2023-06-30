#ifndef NETWORKLIFDELTA_H
#define NETWORKLIFDELTA_H

#include <vector>
#include "jmatrix.h"
#include "lifdelta.h"
#include "spiketrain.h"
#include <iomanip>

class LIFDeltaNetwork
{
 public:
    LIFDeltaNetwork(int AN_,int AExt_, JMatrix *J_, double tau_, vector<double> gamma_, double th_,vector<double> Vinf_, int RN_, bool tracesFlag_=false);
    void VariablesToInitialValues(vector<double> V0s_, SpikeTrain<Spike> *InputSpikes_);
    void Evol ( double tanf, double tend, vector<double> V_inf_, int RN_ );    
    SpikeTrain<Spike> *InputSpikes;
    SpikeTrain<Spike> RecurrentSpikes;
    vector< deque<double> > traces;
    void traces_output();
    double end_time;
private:
    int AN;//Number of recurrent neurons
    int AExt;//Number of external (feedforward) neurons, including Poisson and perturbation neurons
    vector<LIFDelta> Neurons;//Recurrent neurons
    double tau;//delay
    JMatrix *J;// Pointer to the coupling matrix, non-sparse implementation
    bool tracesFlag;
 
friend std::ostream& operator<<(std::ostream& o , LIFDeltaNetwork const &network);
};


 

#endif // NETWORKLIFDELTA_H
