#include <fstream>
#include "lifdeltanetwork.h"
#include <iomanip>

LIFDeltaNetwork::LIFDeltaNetwork(int AN_, int AExt_, JMatrix* J_, double tau_, std::vector< double > gamma_, double th_, vector<double> Vinf_, int RN_, bool tracesFlag_) : AN ( AN_ ),AExt ( AExt_ ),J(J_),tau(tau_),tracesFlag(tracesFlag_)
{
    Neurons = vector<LIFDelta> (AN);

    for ( int i=0; i<AN; i++ )
    {
      Neurons[i]=LIFDelta(gamma_[i],th_,Vinf_[i]);
// 	if (i<RN_)
// 	  Neurons[i]=LIFDelta(-1.7*gamma_*(1.0+rand()/((double)RAND_MAX)),th_,-Vinf_);
// 	else
// 	  Neurons[i]=LIFDelta(gamma_*(1.0+rand()/((double)RAND_MAX)),th_,Vinf_);
    }
}

void LIFDeltaNetwork::VariablesToInitialValues(vector<double> V0s_,SpikeTrain<Spike> *InputSpikes_)
{
    InputSpikes=InputSpikes_;
    RecurrentSpikes.clear();
    for ( int i=0; i<AN; i++ )
    {
        (Neurons[i]).V0=V0s_[i];
        deque<double> trace; // Create empty voltage trace
        traces.push_back(trace); // Add to traces
    }
    deque<double> trace; // Create empty time trace
    traces.push_back(trace); // Add to traces (N+1)th element
}


void LIFDeltaNetwork::Evol ( double tanf, double tend, vector<double> Vinf_, int RN_ )
{
    int no=0;
    SpikeTrain<Spike>::iterator NextArrivingInput=InputSpikes->begin();
    int NextArrivingRecurrentInt(0);
    double t ( tanf );
    double NextArrivingInput_time,NextArrivingRecurrent_time,NextSpike_time,NextEvent_time;
    int NextSpike_sender,ArrivingSpike_sender,NextEvent;
    std::vector<bool> rescue (AN - RN_);
    ofstream f;
    f.open("output2.txt");
    while ( t<tend )
    {
        if ( NextArrivingInput!=InputSpikes->end() ) NextArrivingInput_time=NextArrivingInput->time;
        else NextArrivingInput_time=tend+1;
		
        if ( NextArrivingRecurrentInt!=RecurrentSpikes.size() ) NextArrivingRecurrent_time=RecurrentSpikes[NextArrivingRecurrentInt].time;
        else NextArrivingRecurrent_time=tend+1;

        NextSpike_time=tend+1;
        NextSpike_sender=AN+1;
        for (int i=0; i<AN; i++)
        {
            if ((Neurons[i]).pseudo_spike_time(t,tend)<NextSpike_time)
            {
                NextSpike_time=(Neurons[i]).pseudo_spike_time(t,tend);
                NextSpike_sender=i;
            }
        }

 //       cout << NextArrivingInput_time <<", " << NextArrivingRecurrent_time <<", " << NextSpike_time <<", " << NextSpike_sender << endl;

        if (NextArrivingInput_time<min(NextArrivingRecurrent_time,NextSpike_time))
        {
            NextEvent=0; //Input spike arrival
            NextEvent_time=NextArrivingInput_time;
            ArrivingSpike_sender=NextArrivingInput->sender;
            NextArrivingInput++;
        } else {
            if (NextArrivingRecurrent_time<NextSpike_time)
            {
                NextEvent=1; //Recurrent spike arrival
                NextEvent_time=NextArrivingRecurrent_time;
                ArrivingSpike_sender=RecurrentSpikes[NextArrivingRecurrentInt].sender;
		NextArrivingRecurrentInt++;
            } else {
                NextEvent=2; //Spike sending (werden prioritaer behandelt)
                NextEvent_time=NextSpike_time;
            }
        }

        if (NextEvent_time>tend)
        {
            t=tend;
            break;
        } else {//Evolution bis zum Event
            for (int i=0; i<AN; i++)
            {
                (Neurons[i]).V0=(Neurons[i]).Vvalue(t,NextEvent_time);
            }
            t=NextEvent_time;
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////
        if ((NextEvent==0)||(NextEvent==1))//Spike receival
        {
            for (int i=0; i<AN; i++)
            {
		        if (i<RN_ || (Neurons[i]).V0>0.){
                    (Neurons[i]).V0=(Neurons[i]).V0+( *J )[i][ArrivingSpike_sender];
                    if(i >= RN_){
                    rescue[i-RN_] = 0;
                    }
             //       rescue[i]=0; DELET this line when going back to only rescuing convex neurons
                }
                else if((*J)[i][ArrivingSpike_sender] != 0.){
                   rescue[i-RN_] = 1; 
                        
                  }
            }
        }

        if (NextEvent==2)//Spike sending
        {
            Spike NewSpike(NextEvent_time+tau,NextSpike_sender);
            RecurrentSpikes.push_back(NewSpike);
            (Neurons[NextSpike_sender]).V0=0;
        }
	
//	vector<deque<double>::iterator> it;
//	it=vector<deque<double>::iterator>(AN+1);
//	for (int r=0;r<=AN;r++)
//	  it[r]=traces[r].begin();
	
        if(tracesFlag)
        {
            for (int i=0; i<AN; i++)
            {
                //(Neurons[i]).V0=(Neurons[i]).V0+( *J )[i][ArrivingSpike_sender];
                traces[i].push_back((Neurons[i]).V0);
//		g << traces[i][no] << ' ';
            }
            traces[AN].push_back(t);
//	    g << traces[AN][no] << endl;
//	    no++;
	    
        }


//	cout << t << ": ";
//        cout << RecurrentSpikes;
    if ((NextEvent==0)||(NextEvent==1))//Spike receival
    {
    for(int i=RN_; i>= RN_ && i<AN;i++){  // THIS IS FOR ONLY RESCUING CONVEX NEURONS 
        cout << rescue[i-RN_] <<" ";
    }
    cout << endl;
    }
  /* 
    if ((NextEvent==0)||(NextEvent==1))//Spike receival
    {
    for(int i=0; i<AN;i++){  // THIS IS FOR RESCUING ALL NEURONS 
        cout << rescue[i] <<" ";
    }
    cout << endl;
    }   
*/
    /*
    if (NextEvent_time < 410000){
        for(int i=0;i<AN;i++){ // the if clause here is hardcoded to fit with the networksim.m script
            f << std::setprecision(18) << Neurons[i].V0 <<" ";
        }
     f << t <<" ";
     f << no;
     f << endl;
    
     no++;
     }
     */
    }
    

    end_time=t;
    f.close();
}


void LIFDeltaNetwork::traces_output()
{
    ofstream g;
    
    g.open("Traces.dat");
    
    for (int i=0; i<=AN; i++)
    {
        for (deque<double>::iterator traceIt = traces[i].begin(); traceIt != traces[i].end(); traceIt++)
        {
            g << ',' << *traceIt;
        }
        g << endl;
    }
    
    g.close();
}


ostream& operator<< ( ostream& o , LIFDeltaNetwork const &network )
{
    JMatrix J;
    J=*network.J;
    o<< "Network "<< endl;
    o<< network.AN << " rec. neurons, "<< network.AExt << " ext. neurons." << endl;
    o<< "Rec. neurons: "<< endl;
    vector<LIFDelta> Neurons=network.Neurons;
    for ( vector<LIFDelta>::iterator Neuron=Neurons.begin(); Neuron!=Neurons.end(); Neuron++ )
    {
        o<< ( *Neuron );
    }

    o<<"Coupling matrix" << endl;
    for ( vector<vector<double> >::iterator postsynneuron=J.begin(); postsynneuron!=J.end(); postsynneuron++ )
    {
        for ( vector<double>::iterator matrixelement= ( *postsynneuron ).begin(); matrixelement!= ( *postsynneuron ).end(); matrixelement++ )
        {
            o<< ( *matrixelement ) << " ";
        }
        o << endl;
    }

    return o;
}

