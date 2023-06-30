#ifndef SPIKETRAIN_H
#define SPIKETRAIN_H

#include <deque>
#include <vector>
#include "spike.h"


using namespace std;

template <class SpikeType> class SpikeTrain : public vector< SpikeType >
{
public:
    SpikeTrain()
    {}
    virtual ~SpikeTrain()
    {}

    friend std::ostream& operator<<(std::ostream& o , SpikeTrain <Spike> spiketrain);
    friend std::ostream& operator<<(std::ostream& o , SpikeTrain <DesSpike> spiketrain);
};


class InputSpikeTrains : public vector< SpikeTrain<Spike> >
{
public:
  InputSpikeTrains(int AN,int AExt,int P,string & input_spikes_fn);
  void output();
};

#endif // SPIKETRAIN_H
