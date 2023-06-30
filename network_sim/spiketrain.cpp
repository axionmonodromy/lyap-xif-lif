#include "spiketrain.h"
#include "networkdynioutilities.h"

#include <iostream>


/*
std::ostream& operator<<(std::ostream& o , SpikeTrain<Spike> spiketrain)
{
  int i=1;
  for (SpikeTrain<Spike>::iterator currentspike=spiketrain.begin();currentspike!=spiketrain.end();currentspike++ )
                  {	 
                        o<< (*currentspike) << endl;
			i++;
                  }
cout << endl;
    return o;
}
*/

std::ostream& operator<<(std::ostream& o, SpikeTrain<Spike> spiketrain)
{
   for (SpikeTrain<Spike>::iterator currentspike=spiketrain.begin();currentspike!=spiketrain.end();currentspike++)
        {
            o<<std::setprecision(18) << (*(currentspike)).time - (*(currentspike-1)).time << ',' << (*currentspike).sender << endl;
        }
    cout << endl;
    return o; 

}


/*
std::vector<double> delta_t(SpikeTrain<Spike> spiketrain)
{
    std::vector<double> delta_t = vector<double>(spiketrain.size());
    int i=0;
    for(SpikeTrain<Spike>::iterator currentspike=spiketrain.begin();currentspike!=spiketrain.end();currentspike++){
        delta_t[i] =  *(currentspike+1).time - *currentspike.time ;
        i++;
    }
    
    return delta_t;
    
}
*/

std::ostream& operator<<(std::ostream& o , SpikeTrain<DesSpike> spiketrain)
{
  int i=1;
  for (SpikeTrain<DesSpike>::iterator currentspike=spiketrain.begin();currentspike!=spiketrain.end();currentspike++ )
                  {	 
                        o<< i << "." << (*currentspike) << ", ";
			i++;
                  }
cout << endl;
    return o;
}





InputSpikeTrains :: InputSpikeTrains(int AN,int AExt,int P,string & input_spikes_fn) : vector< SpikeTrain<Spike> >(P,SpikeTrain<Spike>())
{
  SpikeTrainsFromFile<Spike>(AN,AN+AExt,AN, P, input_spikes_fn,(*this));
}

void InputSpikeTrains :: output()
{
  for(int mu=0;mu<this->size();mu++)
  {
		std::cout << "Inputspiketrain " << mu+1 << ": " << endl;
                std::cout << (*this)[mu];
  }
}
