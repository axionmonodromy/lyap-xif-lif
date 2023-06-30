#ifndef NETWORKDYNIOUTILITIES_H
#define NETWORKDYNIOUTILITIES_H

#include <iostream>
#include <deque>
#include <stdexcept>
#include <string>
#include <fstream>

#include "spiketrain.h"
#include "lifdelta.h"
#include "lifdeltanetwork.h"


using namespace std;

template <class SpikeType> void SpikeTrainsFromFile ( int FirstNeuron, int LastNeuronP1, int AddToSender, int P, string & spiketrains_fn, vector< SpikeTrain<SpikeType> > & SpikeTrains )
{
    ifstream in_file ( spiketrains_fn.c_str(), std::ios::binary );
    if ( in_file.is_open() )
    {
        char int_data[sizeof ( int ) ],double_data[sizeof ( double ) ];
        in_file.read ( int_data, sizeof ( int ) );
        int NN = * ( reinterpret_cast<int *> ( int_data ) );
        in_file.read ( int_data, sizeof ( int ) );
        int PP = * ( reinterpret_cast<int *> ( int_data ) );
        if ( ( NN != LastNeuronP1-FirstNeuron ) || ( PP != P ) )
        {
            cout<<NN<<" "<<LastNeuronP1<<" "<<FirstNeuron<<", "<<PP<<" "<<P<<std::endl;
            throw std::runtime_error ( std::string ( "ERROR: Wrong parameters ((NN!=Number of neurons)||(PP!=number of patterns_))" ) );
        }
        SpikeTrains.resize ( PP );
        for ( int mu = 0; mu < PP; ++mu )
        {
            SpikeTrains[mu]=SpikeTrain<SpikeType>();
            in_file.read ( int_data, sizeof ( int ) );
            int N_spikes = * ( reinterpret_cast<int *> ( int_data ) );
            //std::cout<<"\n"<<mu<<" "<<N_spikes<<std::endl;
            int sender;
            double t,eps;
            for ( int s = 0; s < N_spikes; ++s )
            {
                in_file.read ( int_data, sizeof ( int ) );
                sender = ( * ( reinterpret_cast<int *> ( int_data ) ) ) +AddToSender;
                in_file.read ( double_data, sizeof ( double ) );
                t = * ( reinterpret_cast<double *> ( double_data ) );
                in_file.read ( double_data, sizeof ( double ) );
                eps = * ( reinterpret_cast<double *> ( double_data ) );
                //std::cout<<aff<<" "<<t<<std::endl;
                if ( ( sender >= LastNeuronP1 ) || ( sender < FirstNeuron ) )
                {
                    std::cout<<mu<<" "<<s<<" "<<N_spikes<<", "<< FirstNeuron << " < "<<sender<<" < "<<LastNeuronP1<< " ? "<<std::endl;
                    throw std::runtime_error ( std::string ( "ERROR: Too small/large input neuron index" ) );
                }
                SpikeTrains[mu].push_back ( SpikeType ( t,sender,eps ) );
            }
        }
    }
    else
        throw std::runtime_error ( std::string ( "ERROR: Failed to open input spikes file" ) );
}


class CurrentBasedLIFsFromFile : public vector<LIFDelta>
{
public:
    CurrentBasedLIFsFromFile ( int AN, string & neuron_params_fn );
};


template<class T>
class ParameterVector : public vector< T >
{
public :
    ParameterVector ( int N, string fn ) : vector< T > ( N )
    {
        char int_data[sizeof ( int ) ], T_data[sizeof ( T ) ];
        ifstream I_f ( fn.c_str(), ios::binary );
        if ( I_f.is_open() )
        {
            I_f.read ( int_data, sizeof ( int ) );
            int NN=* ( reinterpret_cast<int *> ( int_data ) );
            if ( NN != N )
                throw runtime_error ( string ( "Wrong N in parameter file" ) );
            for ( int i=0; i<N; ++i )
            {
                I_f.read ( T_data, sizeof ( T ) );
                ( *this ) [i]= * ( reinterpret_cast< T * > ( T_data ) );
            }
        }
        else
            throw std::runtime_error ( string ( "Could not open parameter file" ) );
        I_f.close();
    }
};

class MixedOutputUtililties
{

public :
  MixedOutputUtililties();
    void SuccessOutputFile(double err, double t, long needed_runs,string & learn_success_fn,string name_extension="_success");
    void SpikeTrainOutputBinaryFile(SpikeTrain<Spike> & spiketrain,string & output_fn,string name_extension="_spikes");
    void ThetaOutputBinaryFile(vector<double> thetas,string & output_fn,string name_extension="_spikes");
    void EpilepticOutputFile(bool epileptic,double end_time,string & output_fn,string name_extension="_epileptic");
};





#endif // NETWORKDYNIOUTILITIES_H
