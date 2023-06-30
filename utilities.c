#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spmatrix.h>
#include <math.h>
#include "utilities.h"
#include <omp.h>
#include <stdlib.h>

/*
struct single_spike_jacobian{
    int spiking_neuron;
    double *diagonal;
    double *off_diagonal;
};
*/
gsl_matrix * matrix_mult(gsl_matrix *A, gsl_matrix *B, int N) {
    gsl_matrix *C = gsl_matrix_alloc(N,N);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, A, B,0.0, C);
    return C;
}


double * back_sub (int N, double * R, double * c){
    double * vec = malloc(N*sizeof(double));
    int l = (N*(N+1))/2-1;
    int i;
    for(i=N-1; i > -1;i--){
        *(vec+i) = *(c+i);
        for(int j = N-1; j > i; j--){
            *(vec+i) -= *(R+l)*(*(vec+j));
            l--;
        }
        *(vec+i) = *(vec+i)/ *(R+l);
        l--;
    }
    free(c);
    return vec;
}


int sign(double x){
    return (x > 0) - (x < 0);
}

/*
gsl_matrix * Jacobian(int sender, double * gamma, double delta_t, int N, double V_th, double **C){
    
    gsl_matrix * J  = gsl_matrix_alloc(N,N);
    double V_inf    = 3.*sign(*(gamma+sender));
   // double V_infb   = 3.*sign(*(gamma+sender+1));
    double V_k      = (V_th - V_inf * (1.-exp(-*(gamma+sender)*delta_t)))/exp(-*(gamma+sender) * delta_t);
    double temp;
    double diag;
    gsl_matrix_set_zero(J);
    
    for(int i=0; i<N; i++){
        diag = exp(-*(gamma+i)*delta_t);
        gsl_matrix_set (J, i, i, diag);
    }
    
    for(int i=0; i<N; i++){
        temp = (*(gamma+i) / *(gamma+sender)) * (*(*(C+i)+sender)/(V_k-V_inf));
        gsl_matrix_set (J, i, sender, temp);
        
    }
    
    temp = exp(- *(gamma+sender) * delta_t) - V_th/(V_k-V_inf) ;
    gsl_matrix_set(J, sender, sender, temp);
            
    return J;
}
*/

gsl_matrix * Jacobian(int sender, double * gamma, double delta_t, int N, double V_th, double **C, double * V_inf){
    
    gsl_matrix * J  = gsl_matrix_alloc(N,N);
   // double V_inf    = 3.*sign(*(gamma+sender));
   // double V_infb   = 3.*sign(*(gamma+sender+1));
    double temp;
    double diag;
  //  gsl_matrix_set_zero(J);
    
    for(int i=0; i<N; i++){
        diag = exp(-*(gamma+i)*delta_t);
        gsl_matrix_set (J, i, i, diag);
    }
    
    for(int i=0; i<N; i++){
        temp = (*(gamma+i) / *(gamma+sender)) * (*(*(C+i)+sender) * exp(-*(gamma+sender)*delta_t)/(V_th-*(V_inf+sender)));
        gsl_matrix_set (J, i, sender, temp);
        
    }
    
    temp = exp(- *(gamma+sender) * delta_t) - exp(- *(gamma+sender) * delta_t)*V_th/(V_th-*(V_inf+sender));
    gsl_matrix_set(J, sender, sender, temp);
          
    return J;
}

struct single_spike_jacobian ss_jacobian(int sender, double *gamma, double delta_t, int N, double V_th, double **C, double * V_inf, int *rescue,int NC) {
    struct single_spike_jacobian  ssj;
    ssj.diagonal = malloc(N*sizeof(double));
    ssj.off_diagonal = malloc(N*sizeof(double));
    ssj.spiking_neuron = sender;
    for(int i=0;i<N;i++){
        if(i==sender){
            ssj.diagonal[sender] = exp(- *(gamma+sender) * delta_t) - exp(- *(gamma+sender) * delta_t)*V_th/(V_th-*(V_inf+sender));
        }
        else{
            ssj.diagonal[i] = exp(-gamma[i]*delta_t);
        }
    }
    
    for(int i=0;i<N;i++){
        if(i==sender){
            ssj.off_diagonal[sender] = 0.;
        }
        else{
            ssj.off_diagonal[i] =(*(gamma+i) / *(gamma+sender)) * (*(*(C+i)+sender) * exp(-*(gamma+sender)*delta_t)/(V_th-*(V_inf+sender)));
        }
        if(i >= N-NC){
         if(rescue[i-(N-NC)] == 1){
             ssj.off_diagonal[i] = 0.;
         }
        }
    }
    return ssj;
}


void ssj_mult (struct single_spike_jacobian ssj,int N, gsl_matrix * matrix) {
    gsl_matrix * A = gsl_matrix_alloc(N,N);
    gsl_matrix * B = gsl_matrix_alloc(N,N);
    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B); 
    int i,j,k,l;
    double temp1,temp2;
    
//    #pragma omp parallel for private(i,j,temp1) shared(A,matrix,ssj) collapse(2)
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            temp1 = gsl_matrix_get(matrix,i,j) * ssj.diagonal[i];
            gsl_matrix_set(A,i,j,temp1);
        }
    }   

    
//    #pragma omp parallel for private(l,k,temp2) shared(B,matrix,ssj) collapse(2)
    for(k=0;k<N;k++){
        for(l=0;l<N;l++){
            temp2 = gsl_matrix_get(matrix,ssj.spiking_neuron,l) * ssj.off_diagonal[k];
            gsl_matrix_set(B,k,l,temp2);
        }
    }
    

    gsl_matrix_add(A,B);
    gsl_matrix_memcpy(matrix,A);
    gsl_matrix_free(B);
    gsl_matrix_free(A);
}

void ssj_free (struct single_spike_jacobian ssj){
    free(ssj.diagonal);
    free(ssj.off_diagonal);
}

double norm(double *X, int N){
    double temp=0;
    for(int i=0;i<N;i++){
        temp += (*(X+i))*(*(X+i));
    }
    return sqrt(temp);
}


void printut(double *R, int N){
    int loop = N;
    int j=0;
    double x=0;
    while(loop>0){
        for(int i=0;i<N-loop;i++){
            printf("%f ",x);
        }
        for(int i=0;i<loop;i++){
            printf("%f ", *(R+j));
            j++;
        }
            printf("\n");
            loop--;
    }
}

void printmat(gsl_matrix * M, int N){
    for(int i=0; i<N;i++){
        for(int j=0;j<N;j++){
            printf("%f ",gsl_matrix_get(M,i,j));
        }
        printf("\n");
    }
}


void normalize(double *V, int N){
    double temp = norm(V,N);
    for(int i=0;i<N;i++){
        *(V+i) = *(V+i) /temp;
    }
}
    
    
void matrixcopy(gsl_matrix *A, double **B, int N){
    for(int i=0; i<N;i++){
        for(int j=0; j<N;j++){
      gsl_matrix_set(A,i,j,*(*(B+j)+i));
            
        }
   }
}

void matrix_normalize(gsl_matrix *A, int N){
    double *a = malloc(N*sizeof(double));
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            *(a+j) = gsl_matrix_get(A,j,i);
        }
        normalize(a,N);
        for(int j=0;j<N;j++){
            gsl_matrix_set(A,j,i,*(a+j));
        }
    }
}
            
    
    
    
 
void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}
