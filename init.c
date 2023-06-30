#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "utilities.h"
#include <time.h>







int main () {

srand( time(NULL));
    
int RN=60;					// The number of concave rise neurons
int IN=100;					    // The total number of neurons
int EX=0;					    	// The number of output neurons (currently must be set to 0)
int tau=0;					    	// The synaptic delay (for the current model, we work for tau=0)
double tf=3000000.0;					// Final time
double target = -10.0;					// The sum of all inputs to any one neuron
double t_trans=300000.;			          //transient time for simulation
int indegree = 50; 
double weight = 1.0;


double * gamma = malloc(IN*sizeof(double));
for(int i=0;i<IN;i++){
    if(i<RN){
    *(gamma+i) =  0.01685681377845419;
    }
    else {
    *(gamma+i) = -0.01; 
    }
}


double * vinf = malloc(IN*sizeof(double));

for(int i=0;i<IN;i++){
    if(i<RN){
        *(vinf+i) = 2.0; //-1. + 2*(double)(rand())/(double)(RAND_MAX);

    }
    else{
        *(vinf+i) = -2.0;// - 1. + 2*(double)(rand())/(double)(RAND_MAX);
    }
}
double th=1.0;						// The threshold potential
int ti=0;						// Inital time
    
double * V_init = malloc(IN*sizeof(double));
for(int i=0;i<IN;i++){
    *(V_init +i) = th*(double)(rand())/(double)(RAND_MAX);
}


double ** C_matrix;
C_matrix = malloc(IN*sizeof(*C_matrix));
for(int i=0; i<IN; i++){
    *(C_matrix+i) = malloc(IN*sizeof(double));
}

/*
for(int i=0;i<IN;i++){
    for(int j=0; j<IN; j++){
        if(i==j){
        *(*(C_matrix+j)+i) = 0.;    
        }
        else{
            *(*(C_matrix+j)+i) = -1.;
        }
    }
}


for(int i=0;i<IN;i++){
    for(int j=0;j<IN;j++){
        if( (double)(rand())/(double)(RAND_MAX) < 0.2 )
            *(*(C_matrix+j)+i) = 0.;
    }
}

*/
int * synapses = malloc((IN-1)*sizeof(double));
for(int i=0;i<IN-1;i++){
   if(i<indegree){
        synapses[i] = 1;
    }
    else{
        synapses[i]= 0;
    }
}

for(int i=0;i<IN;i++){
    shuffle(synapses,IN-1);
    for(int j=0;j<IN;j++){
        if(j==i) *(*(C_matrix+i)+j)= 0.; 
        else if(j<i){
        *(*(C_matrix+i)+j) = synapses[j]*weight;
        }
        else{
        *(*(C_matrix+i)+j) = synapses[j-1]*weight;
        }
    }
}




 
         


for(int i=0;i<IN;i++){
    double norm=0;
    for(int j=0; j<IN;j++){
        norm += *(*(C_matrix+i)+j);
    }
    for(int j=0; j<IN;j++){
        *(*(C_matrix+i)+j) = *(*(C_matrix+i)+j)*(target)/norm;
    }
}


    FILE *f;
    
    f = fopen("./Gamma.dat","w");
    for(int i=0;i<IN;i++){
        fprintf(f,"%14.12f\n",*(gamma+i));
    }

    fclose(f);
    
    f = fopen("./parameters.txt","w");
    fprintf(f,"%i %i %i %f %i %f %i",IN,EX,tau,th,ti,tf,RN);
    
    fclose(f);
    
    
    f = fopen("./InitVolt.txt","w");
    for(int i=0;i<IN;i++){
        fprintf(f,"%lf\n",*(V_init+i));
    }
    
    fclose(f);
    
    f = fopen("./JMatrix.txt","w");
    
    for(int i=0;i<IN;i++){
        for(int j=0;j<IN;j++){
            fprintf(f,"%lf ", *(*(C_matrix+i)+j));
        }
        fprintf(f,"\n");
    }
    
    fclose(f);
    
    f = fopen("./VInf.dat","w");
    for(int i=0;i<IN;i++){
        fprintf(f,"%lf\n", *(vinf+i));
    }
    
    fclose(f);
    
    
    f = fopen("./ttrans.txt","w");
    fprintf(f,"%lf\n",t_trans);
    fclose(f);
    

    
    
    return(0);
}
