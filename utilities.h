#ifndef UTILITIES_H
#define UTILITIES_H


#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spmatrix.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>

struct single_spike_jacobian {
    int spiking_neuron;
    double *diagonal;
    double *off_diagonal;
};
gsl_matrix * matrix_mult(gsl_matrix *A, gsl_matrix *B, int N);
double * back_sub (int N, double * R, double * c);
int sign(double x);
gsl_matrix * Jacobian(int sender, double * gamma, double delta_t, int N, double V_th, double **C, double * V_inf);
double norm(double *X, int N);
void printut(double *R, int N);
void printmat(gsl_matrix * M, int N);
void normalize(double *V, int N);
void matrixcopy(gsl_matrix *A, double **B, int N);
void matrix_normalize(gsl_matrix *A, int N);
struct single_spike_jacobian ss_jacobian(int sender, double *gamma, double delta_t, int N, double V_th, double **C, double * V_inf, int *rescue,int NC);
void ssj_mult (struct single_spike_jacobian ssj, int N ,gsl_matrix * matrix);
void ssj_free(struct single_spike_jacobian ssj);
void shuffle(int * array, size_t n);






#endif
