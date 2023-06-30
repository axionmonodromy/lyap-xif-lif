#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "utilities.h"




int main (){
    int N;                                                        //total number of neurons
    int NX;                                                      //number of convex neurons 
    double t_trans;                                           //transient time for simulation
    int r_trans=0;
    int bases=1;                                //number of gs bases we want to save for later computation of CLV
    double t_total;
    int l=0;
    int diag=0;
    int spikes=0;                                                 //total number of spikes during simulation
    double V_th;
    double ddummy;
    int idummy;
    
    FILE *f;
    
    f=fopen("parameters.txt","r");
    fscanf(f,"%i 0 0 %lf 0 %lf %i\n",&N,&V_th,&t_total,&NX);
    fclose(f);
    

    f=fopen("ttrans.txt","r");
    fscanf(f,"%lf\n",&t_trans);
    fclose(f);
    
    f = fopen("Spiketimes.dat","r");
    
    if((f = fopen("./Spiketimes.dat","r")) != NULL) {
        while( fscanf(f,"%lf,%i\n",&ddummy,&idummy) == 2){
            spikes++;
        }
        
    }
    else{
        printf("Error 404 Spiketimes\n");
        return 0;
    }
    rewind(f);
    
    double * delta_t = malloc(spikes*sizeof(double));
    double total_time = 0.;
    int *sender = malloc(spikes*sizeof(int));                           
    double *spiketime = malloc(spikes*sizeof(double));                  
    
    for(int i=0;i<spikes-2;i++){
        fscanf(f,"%lf,%i\n",delta_t+i,sender+i);
        total_time += *(delta_t+i);
        *(spiketime+i)= total_time;
        if(total_time > t_trans && r_trans == 0){
            r_trans = i;
                }
    }


    fclose(f);
   
   
    double *gamma;                                                      //inverse membrane constants
    gamma = malloc(N*sizeof(double));
    
    f=fopen("Gamma.dat", "r");
    for(int i = 0; i<N;i++){
        fscanf(f,"%lf\n",gamma+i);
    }
    fclose(f);
    
    
    double **R_storage;                                                 //storage for the upper triangular matrices
    gsl_matrix ** Q_storage;                                            //storage for the gram schmidt matrices
    double **C;                                                         //coupling matrix
    
    
    R_storage = malloc(spikes*sizeof(*R_storage));
    for(int i=0;i<spikes;i++){
       *(R_storage +i) = malloc((N*(N+1)/2)*sizeof(double));
    }
    
    
    C = malloc(N*sizeof(*C));
    for(int i=0;i<N;i++){
        *(C+i) = malloc(N*sizeof(double));
    }
    
    f = fopen("JMatrix.txt","r");
    for(int i = 0; i < N*N; i++){
        fscanf(f,"%lf ", *(C+i/N) + i%N);
    }
    
    fclose(f);
    
    double * V_inf = malloc(N*sizeof(double));
    
    f = fopen("VInf.dat","r");
    for(int i = 0; i<N; i++){
        fscanf(f,"%lf\n",V_inf+i);
    }
    fclose(f);


                                                                        //matrices to be used for the gsl QR decomposition
    gsl_matrix *Q = gsl_matrix_alloc(N,N);                              
    gsl_matrix *R = gsl_matrix_alloc(N,N);
    gsl_matrix_set_identity(Q);
    gsl_matrix *temp_matrix = gsl_matrix_alloc(N,N);
    gsl_vector *tau = gsl_vector_alloc(N);
     
    Q_storage = malloc(bases*sizeof(*Q));
    for(int i=0;i<bases;i++){
        *(Q_storage+i) = gsl_matrix_alloc(N,N);
    }
    
    double * lyapunov_spectrum = malloc(N * sizeof(double));
    for(int i=0;i<N;i++){
        *(lyapunov_spectrum+i) = 0;
    }

    int **rescue;
    rescue = malloc((spikes-2)*sizeof(*rescue));
    for(int i=0;i<spikes-2;i++){
        *(rescue+i) = malloc((N-NX)*sizeof(int));
     }

    f = fopen("output.txt","r");
    for(int i=0;i< (spikes-2);i++){
        for(int j=0;j<N-NX;j++){
            fscanf(f,"%i ", *(rescue+i)+j);
            }
        fscanf(f,"\n");
    }
    fclose(f);


 

/////////////// Algorithm for LE starts here ///////////////////////////           
        
    for(int r=0;r<spikes-5;r++){
       
        struct single_spike_jacobian SSJ = ss_jacobian (*(sender+r), gamma, *(delta_t+r), N, V_th, C,V_inf,*(rescue+r),N-NX);
        ssj_mult(SSJ,N,Q);
        ssj_free(SSJ);
        free(*(rescue+r)); 

        gsl_linalg_QR_decomp(Q,tau);                                                    //QR decomposition
        gsl_matrix_memcpy(temp_matrix,Q);                
        gsl_linalg_QR_unpack(temp_matrix,tau,Q,R);
       
        if(*(spiketime+r)>t_trans-1){
		int loop = 0;
                for(int i=0;i<N;i++){
                   for(int j=i; j<N;j++){
                       *(*(R_storage+l)+loop) = gsl_matrix_get(R,i,j);
		       loop++;
                   }

		}
                
                if(l<bases){
                    gsl_matrix_memcpy(*(Q_storage+l),Q);
                    }
                    
                l++;       
        }

   } 
      
      

    for(int i=0;i<N;i++){
        for(int j=1;j<l-1;j++){
             *(lyapunov_spectrum+i) += log( fabs( *( *(R_storage + j)+ diag ) ) );
        }
        
        *(lyapunov_spectrum+i) = *(lyapunov_spectrum+i) / ( total_time  - t_trans );
        diag += N-i;
    }
   
    
   

    for(int i=0;i<N;i++){
        printf("%f\n",*(lyapunov_spectrum+i));
    }
    
    
    
    /////////////////// compute the CLV ///////////////////////
   
   
    double ** BW;
    BW = malloc(N*sizeof(*BW));
    for (int i=0;i<N;i++){
        *(BW +i) = malloc(N*sizeof(double));
    }
    
    
    
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i>j){
                *(*(BW+j)+i) = 0.;
            }
            else {
                *(*(BW+j)+i) = 1.;
            }
        }
    }
    

    gsl_matrix ** UT;
    UT = malloc(bases*sizeof(*UT));
    for(int i = 0; i<bases;i++){
        *(UT+i) = gsl_matrix_alloc(N,N);
    }
    
    double temp = 0;
    int iter;
    for(int i=l-2; i>1;i--){
        for(iter=0; iter<N;iter++){
            temp = norm(*(BW+iter),N);
            for(int k=0;k<N;k++){
                *(*(BW+iter)+k) = *(*(BW+iter)+k)/temp;
            }
            *(BW+iter) = back_sub(N,*(R_storage+i),*(BW+iter));
        }
        free(*(R_storage+i));
        if(i<bases){
            matrixcopy(*(UT+i),BW,N);
        }
    }
   
   
   gsl_matrix ** CLV;
   CLV = malloc(bases*sizeof(*CLV));
   for(int i=0;i<N;i++){
       *(CLV+i) = gsl_matrix_alloc(N,N);
   }
   

   for(int i=0;i<bases-2;i++){
       *(CLV+i) = matrix_mult(*(Q_storage+i+1),*(UT+i+2),N); 
       matrix_normalize(*(CLV+i),N);
   }
   
  
    

   
   
    f=fopen("CLV.txt", "w");
    for(int k=0;k<bases-2;k++){
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                fprintf(f,"%.12f   ",gsl_matrix_get(*(CLV+k),i,j));
            }
        fprintf(f,"\n");   
        }
    }
    

    
   
return 0;    
    
}
