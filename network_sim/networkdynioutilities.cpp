#include "networkdynioutilities.h"

#include <iterator>
#include "spike.h"

using namespace std;

CurrentBasedLIFsFromFile::CurrentBasedLIFsFromFile ( int AN, string & neuron_params_fn ) : vector<LIFDelta> ( AN )
{
  char int_data[sizeof ( int ) ], double_data[sizeof ( double ) ];
  ifstream neuron_param_file ( neuron_params_fn.c_str(), ios::binary );
  if ( neuron_param_file.is_open() )
    {
      neuron_param_file.read ( int_data, sizeof ( int ) );
      if ( * ( reinterpret_cast<int *> ( int_data ) ) != AN )
        throw runtime_error ( string ( "Wrong AN in neuron parameter file" ) );
      double taum, taus;
      for ( int i = 0; i < AN; ++i )
        {
          neuron_param_file.read ( double_data, sizeof ( double ) );
          taum = * ( reinterpret_cast<double *> ( double_data ) );
          neuron_param_file.read ( double_data, sizeof ( double ) );
          taus = * ( reinterpret_cast<double *> ( double_data ) );
     //     ( *this ) [i] = LIFDelta ( taum,gamma ); Neu!
        }
    }
  else
    throw std::runtime_error ( string ( "Could not open neuron parameter file" ) );
  neuron_param_file.close();
}


MixedOutputUtililties::MixedOutputUtililties()
{}

void MixedOutputUtililties::SuccessOutputFile ( double err, double t, long needed_runs,string & learn_success_fn,string name_extension )
{
  ofstream learn_success_file ( ( learn_success_fn+name_extension ).c_str() );
  if ( err==0 )
    {
      cout<< "Successful learning: "<< t << " runs needed." <<endl;
      learn_success_file<<"1 "<<t<<endl;
    }
  else
    {
      cout<< "Unsuccessful learning after "<< t << " runs." <<endl;
      learn_success_file<<"0 "<<t<<endl;
    }
  learn_success_file.close();
}





void MixedOutputUtililties::SpikeTrainOutputBinaryFile ( SpikeTrain<Spike> & spiketrain,string & output_fn,string name_extension )
{
  cout<<endl<<"Writing spike train to file..."<<flush;
  ofstream spikes_file;
   int P ( 1 );
  int n_spikes = spiketrain.size();
  spikes_file.open ( ( output_fn+name_extension ).c_str(),ios::binary );
  spikes_file.write ( reinterpret_cast<const char *> ( &P ),sizeof ( int ) );
  spikes_file.write ( reinterpret_cast<char *> ( &n_spikes ), sizeof ( int ) );
  for ( int sp = 0; sp < n_spikes; ++sp )
    {
      spikes_file.write ( reinterpret_cast<char *> ( & ( spiketrain[sp].sender ) ), sizeof ( int ) );
      spikes_file.write ( reinterpret_cast<char *> ( & ( spiketrain[sp].time ) ), sizeof ( double ) );
    }
  spikes_file.close();
}

void MixedOutputUtililties::ThetaOutputBinaryFile(vector<double> thetas,string & output_fn,string name_extension)
{
    cout<<endl<<"Writing spike thetas to file..."<<flush;
  int AN(thetas.size());
  ofstream thetas_file;
   thetas_file.open ( ( output_fn+name_extension ).c_str(),ios::binary );
      thetas_file.write(reinterpret_cast<char *> (&AN), sizeof (int));
   for (int i=0; i<AN; ++i)
        thetas_file.write(reinterpret_cast<char *>(&(thetas[i])),sizeof(double));
    thetas_file.close();
}

  void MixedOutputUtililties::EpilepticOutputFile ( bool epileptic,double end_time,string & output_fn,string name_extension )
  {
    ofstream epileptic_file ( ( output_fn+name_extension ).c_str() );
    epileptic_file<<epileptic<<" "<<end_time<<endl;
    epileptic_file.close();
  }



