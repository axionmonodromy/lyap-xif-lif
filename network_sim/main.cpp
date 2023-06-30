#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "lifdelta.h"
#include "jmatrix.h"
#include "lifdeltanetwork.h"
#include "spiketrain.h"
//#include "Jacobian.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
//#include "spike.h"


using namespace std;



int main(int argc, char **argv) {
  
    srand(time(NULL));
 //   cout << "Simulation LIF neurons, delta coupling." << endl;
    
 //   cout << argc << endl;
    for (int qaz=0;qaz<argc;qaz++)
      puts(argv[qaz]);
    
//    cout << atoi(argv[1]) << endl;
    
    int RN(atoi(argv[7]));
    
    int AN(atoi(argv[1]));
//    int AN(5);
    int AExt(atoi(argv[2]));
//    int AExt(2);
    //
    
     vector<double> gamma;//1/membrane time constant
     gamma=vector<double>(AN);
     
     ifstream gm;
     gm.open(argv[10]);
     
     int gmc=0;
     double gmd;
     
      while (gmc<AN)
      {
	gm>>gmd;
	gamma[gmc]=gmd;
	gmc++;
//	cout << gmd << '\t' << gmc << endl;
      }
      
      gm.close();
     
//     double gamma(-0.1);
//      cout << "gamma " << gamma <<endl;
    double th(atof(argv[4]));//threshold
//     double th(1.5);
//     cout << "th " << th <<endl;
//     double Vinf(-2.0);
     
     vector<double> Vinf;//Potential that is assumed for t->\infty without inputs
     Vinf=vector<double>(AN);
     
     ifstream vf;
     vf.open(argv[11]);
     
     int vfc=0;
     double vfd;
     
     while (vfc<AN)
     {
       vf>>vfd;
       Vinf[vfc]=vfd;
       vfc++;
  //     cout << vfd << '\t' <<vfc << endl;
     }
     
     vf.close();
     
//     cout << "Vinf " << Vinf <<endl;
    double tau=atof(argv[3]);//delay
//    double tau(1.2);
//    cout << "tau " << tau <<endl;
//    cout << "===============================================================================================" << endl;
//     cout << "Vinf=" << Vinf << ", th=" << th << endl;

    vector<double> V0s;//Initial potentials
    V0s=vector<double>(AN);
//    V0s[0]=0.871434;
//    V0s[1]=1.2818;
//    V0s[2]=1.06876;
//    V0s[3]=0.972767;
//    V0s[4]=0.764489;
    ifstream z;
    z.open(argv[8]);
    
    int d=0;
    char buffer[100];
    strcpy(buffer,"0");
    double c;
    
    while (d<AN)
    {
      z>>c;
      V0s[d]=c;
      d++;
 //     cout << c << '\t' << d << endl;
    }
 //   cout << endl << d << endl;

    
/*    for (int i=0;i<AN;i++)
    {
      V0s[i]=rand()/(double)RAND_MAX*th;
      
      cout << V0s[i] << endl;
    }
*/    
    SpikeTrain<Spike> InputSpikes;//Some small example external input spike train
    Spike SingleInputSpike;
    
    int No_Spikes=5000;
    
    double ta=atof(argv[5]);
    double tb=atof(argv[6]);
    
//     ofstream l;
//     
//     l.open("InputSpiketimes.dat");
    
    ifstream wsx;
    
    wsx.open("Inputs.dat");
    
    double gh;
    
    if (AExt>0)
    {
      while (!wsx.eof())
      {
//	cout << "qaz " << AN+rand()/double(RAND_MAX) << endl;
	SingleInputSpike.sender=AN+int(rand()/double(RAND_MAX)*AExt);
// 	SingleInputSpike.time=(rand()/double(RAND_MAX)+i)*(tb-ta)/(No_Spikes);
	wsx >> gh;
	SingleInputSpike.time=gh;
// 	l << SingleInputSpike;
	InputSpikes.push_back(SingleInputSpike);
      }
    }
    
    wsx.close();
    
  //  cout << "=============================================================" << endl;
  //  cout << InputSpikes;
    
     
//      l << "Arindam";
     
//      l.close();
    
//    JMatrix J=JMatrix(AN, AExt, -0.2, -0.1);//Set up an example coupling matrix
    JMatrix J=JMatrix(AN, AExt, argv[9]);
    
    J.output();
    
  //  cout<<"==================================================================================="<<endl;
    LIFDeltaNetwork network=LIFDeltaNetwork(AN,AExt,&J,tau,gamma,th,Vinf,RN);//Create the network
  //  cout<<"==================================================================================="<<endl;
    network.VariablesToInitialValues(V0s,&InputSpikes);//Initialize the network
  //  cout<<"==================================================================================="<<endl;
  //  cout << network;
  //  cout<<"==================================================================================="<<endl;
    
    
     network.Evol(atof(argv[5]),atof(argv[6]),Vinf,RN);//Simulate the network
     
//     double final[1000][2];
//     int g;
//     double temp;
     
/*     for (g=0; g<network.RecurrentSpikes.size(); g++)
     {
       final[g][0]=network.RecurrentSpikes[g].time;
       final[g][1]=network.RecurrentSpikes[g].sender;
     }
*/     
//      for (g=network.RecurrentSpikes.size();g>=0; g--)
//        for (h=0; h<g-1;h++)
// 	 if (final[h][0]>final[h+1][0])
// 	 {
// 	   temp=final[h][0];
// 	   final[h][0]=final[h+1][0];
// 	   final[h+1][0]=temp;
// 	   
// 	   temp=final[h][1];
// 	   final[h][1]=final[h+1][1];
// 	   final[h+1][1]=temp;
// 	 }
	 
/*	 
      cout << "Printing final" << endl; 
      for (g=0; g<network.RecurrentSpikes.size(); g++)
	cout << final[g][0] << '\t' << final[g][1] << endl;
*/      
/*      double refsender;
      double reftime;
      
      refsender=final[0][1];
      reftime=final[0][0];
      
      double modfinal[500][3];
      
      modfinal[0][0]=0.0;
      modfinal[0][1]=0.0;
      modfinal[0][2]=final[0][1];
      double y=0;
      
      for (g=1; g<network.RecurrentSpikes.size(); g++)
      {
	if (final[g][1]==refsender)
	{
	  reftime=final[g][0];
	  y++;
	}
	modfinal[g][0]=final[g][0]-reftime;
	modfinal[g][1]=y;
	modfinal[g][2]=final[g][1];
      }

      cout << "Printing modfinal" << endl; 
      for (g=0; g<network.RecurrentSpikes.size(); g++)
	cout << g+1 << "." << '\t' << modfinal[g][0] << '\t' << modfinal[g][1] << ' ' << modfinal[g][2] << endl;
      
      ofstream z;
      
      z.open("Firetimes.dat");
      for (g=0; g<network.RecurrentSpikes.size(); g++)
	z << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
      z.close();
*/      
  //   cout << network.RecurrentSpikes;
     ofstream f;
     
     f.open("Spiketimes.dat");
     
     f << network.RecurrentSpikes;
     
     f.close();
     
     
  //   cout << "Rec. Spikes " << network.RecurrentSpikes.size() << endl;
     
     
//     network.traces_output();
     
    ofstream l;
    
    l.open("InputSpiketimes.dat");
    
    l << InputSpikes;
    
    l.close();
     
/*     f.open("Spike0.dat");
     
     for (g=0;g<network.RecurrentSpikes.size(); g++)
       if (modfinal[g][2]==0)
	 f << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
       
     f.close();
    
    f.open("Spike1.dat");
     
     for (g=0;g<network.RecurrentSpikes.size(); g++)
       if (modfinal[g][2]==1)
	 f << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
       
     f.close();
    
     f.open("Spike2.dat");
     
     for (g=0;g<network.RecurrentSpikes.size(); g++)
       if (modfinal[g][2]==2)
	 f << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
       
     f.close();
     
     f.open("Spike3.dat");
     
     for (g=0;g<network.RecurrentSpikes.size(); g++)
       if (modfinal[g][2]==3)
	 f << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
       
     f.close();
     
     f.open("Spike4.dat");
     
     for (g=0;g<network.RecurrentSpikes.size(); g++)
       if (modfinal[g][2]==4)
	 f << modfinal[g][1] << '\t' << modfinal[g][0] << endl;
       
     f.close(); 
    
*/    
//    J.output();
/*
double t,tnext,delta_t,tf;
int sender;
double Vk;
double a[AN],sum[AN],L[AN];
double norm;

ofstream zxcv;
zxcv.open("Lyapunov.dat");

for (int i=0;i<AN;i++)
{
  a[i]=0.0;
  sum[i]=0.0;
}
// gsl_matrix *M=gsl_matrix_alloc(AN,AN);

double M[AN*AN];
double dV[AN*AN],dVt[AN*AN];

for (int i=0;i<AN;i++)
      for (int j=0;j<AN;j++)
      {
	if (i==j)
	{
	  dV[i*AN+j]=1.0;
	  dVt[i*AN+j]=0.0;
	}
	else
	{
	  dV[i*AN+j]=0.0;
	  dVt[i*AN+j]=0.0;
	}
      }
    

 for (int i=0;i<AN;i++)
   for (int j=0;j<AN;j++)
   {
     if (i==j)
       M[i*AN+j]=1.0;
     else
//      gsl_matrix_set(M,i,j,0.0);
       M[i*AN+j]=0.0;
   }
   
  ifstream q;
  
  q.open("Spiketimes.dat");

  t=0.0;
  tf=atof(argv[8])-100;
  char asdf[100];
  while(t<tf)
  {
//     cout << "=========================================" <<endl;
    q.getline(asdf, 100,',');
    tnext=atof(asdf);
    q.getline(asdf, 100);
    sender=atof(asdf);
    
    delta_t=tnext-t;
    cout << delta_t << '\t' <<tnext << '\t' << t << endl;
    for (int i=0;i<AN;i++)
      for (int j=0;j<AN;j++)
      {
	if (i==j)
	  M[i*AN+j]=1.0;
	else
	  M[i*AN+j]=0.0;
      }

//     cout << sender << endl;
    
    for (int i=0;i<RN;i++)
    {
      M[i*AN+i]=exp(gamma*delta_t);
//      cout << i << '\t' << exp(gamma*delta_t) << endl;
    }
    for (int i=RN;i<AN;i++)
    {
      M[i*AN+i]=exp(-gamma*delta_t);
//      cout << i << '\t' << exp(-gamma*delta_t) << endl;
    }
    
    if (sender<RN)
      Vk=(th+Vinf*(1-exp(gamma*delta_t)))/exp(gamma*delta_t);
    else
      Vk=(th-Vinf*(1-exp(-gamma*delta_t)))/exp(-gamma*delta_t);
    
    for (int i=0;i<RN;i++)
      M[i*AN+sender]=J[i][sender]/(Vk+Vinf);
    for (int i=RN;i<AN;i++)
      M[i*AN+sender]=J[i][sender]/(Vk-Vinf);
    
    if (sender<RN)
        M[sender*AN+sender]=exp(gamma*delta_t)-th/(Vk+Vinf);
    else
        M[sender*AN+sender]=exp(-gamma*delta_t)-th/(Vk-Vinf);

    
//     for (int i=0;i<AN;i++)
//     {
//       for (int j=0;j<AN;j++)
//  	cout<<M[i*AN+j]<<'\t';
//       cout << endl;
//     }
    
//     cout << M[0] << endl;
    
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,AN,AN,AN,1.0,M,AN,dV,AN,1.0,dVt,AN);
    
    for (int i=0;i<AN;i++)
      for (int j=0;j<AN;j++)
	dV[i*AN+j]=dVt[i*AN+j];
      
    for (int i=0;i<AN;i++)
      for (int j=0;j<AN;j++)
	dVt[i*AN+j]=0.0;
    
//     for (int i=0;i<AN;i++)
//     {
//       for (int j=0;j<AN;j++)
// 	cout<<dV[i*AN+j]<<'\t';
//       cout << endl;
//     }
    
    for (int i=0;i<AN;i++)
    {
      for (int j=0;j<i;j++)
      {
// 	cout << "i = " << i << '\t' <<"j = " << j << endl;
	c=0.0;
	for (int k=0;k<AN;k++)
	  c=c+dV[k*AN+j]*dV[k*AN+i];
// 	cout << c << endl;
	for (int k=0;k<AN;k++)
	{
// 	  cout << dV[k*AN+i] << '\t' << c*dV[k*AN+j] << endl;
	  dV[k*AN+i]=dV[k*AN+i]-c*dV[k*AN+j];
 	}
 	
      }
     
     norm=0.0;
     for (int j=0;j<AN;j++)
       norm=norm+dV[j*AN+i]*dV[j*AN+i];
     
     norm=sqrt(norm);
     a[i]=norm;
     
     for (int j=0;j<AN;j++)
       dV[j*AN+i]=dV[j*AN+i]/norm;
     
//      for (int j=0;j<AN;j++)
       sum[i]=sum[i]+log(a[i]);
     
    }

//     cout << endl;
// 
//     for (int i=0;i<AN;i++)
//       cout << a[i] << '\t';
//     cout << endl;
//     
//     cout << endl;
//     
//     for (int i=0;i<AN;i++)
//     {
//       for (int j=0;j<AN;j++)
// 	cout<<dV[i*AN+j]<<'\t';
//       cout << endl;
//     }
    
    for (int i=0;i<AN;i++)
      L[i]=sum[i]/t;

    zxcv << t << '\t';
    for (int i=0;i<AN;i++)
      zxcv << L[i] << '\t';
    zxcv << endl;
    
//     for (int i=0;i<AN;i++)
//       cout <<L[i]<<'\t';
//     
     cout << t << endl;
    
    t=tnext;
  }

  ofstream zaq;
  zaq.open("lfinal.dat");
  
  for (int i=0;i<AN;i++)
      zaq <<L[i]<<endl;
//   cout << endl;
  */
//  cout << "Arindam" <<endl;
 
//   cout << delta_t << endl;
   
    return 0;
}
