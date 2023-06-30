#include "jmatrix.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>


using namespace std;

JMatrix::JMatrix()
{}


/*JMatrix::JMatrix(int N, int N_input, string fn) : vector< vector <double> >(N,vector<double>(N+N_input)) //First Rows, then columns, Coupling matrix from file
{
    char int_data[sizeof (int) ], double_data[sizeof (double) ];
    ifstream J_f(fn.c_str(), std::ios::binary);

    if (J_f.is_open()) {
        J_f.read(int_data, sizeof (int));
        int NN_input=*(reinterpret_cast<int *> (int_data));
        J_f.read(int_data, sizeof (int));
        int NN=*(reinterpret_cast<int *> (int_data));
        if ((NN != N)||(NN_input!=N_input)) throw runtime_error(string("Wrong N or N_input in J file"));
        for (int j=N; j<N+N_input; ++j)
            for (int i=0; i<N; ++i) {
                J_f.read(double_data, sizeof (double));
                (*this)[i][j]= *(reinterpret_cast<double *> (double_data));
            }
        for (int j=0; j<N; ++j)
            for (int i=0; i<N; ++i) {
                J_f.read(double_data, sizeof (double));
                (*this)[i][j]= *(reinterpret_cast<double *> (double_data));
            }
    } else
        throw std::runtime_error(string("Could not open J file"));
    J_f.close();
}
*/
JMatrix::JMatrix(int N, int N_input, double jl, double jr) : vector< vector <double> >(N,vector<double>(N+N_input)) //Random coupling matrix
{
    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N+N_input; ++j)
	{
            (*this)[i][j]=jl+rand()/(double)RAND_MAX*(jr-jl);
        }
    }
}

JMatrix::JMatrix(int N, int N_input, double jl) : vector< vector <double> >(N,vector<double>(N+N_input)) //Random coupling matrix
{
    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N+N_input; ++j)
	{
	  if (j==i+1)
            (*this)[i][j]=jl;
	  else
	    (*this)[i][j]=0.0;
        }
    }
}

JMatrix::JMatrix(int N, int N_input, char* file) : vector< vector <double> >(N,vector<double>(N+N_input)) //Random coupling matrix
{
    ifstream f;
    double b;
    f.open(file);
    for (int i=0; i<N; ++i)
    {
        for (int j=0; j<N+N_input; ++j)
	{
	  f>>b;
          (*this)[i][j]=b;
        }
    }
}


void JMatrix::output()
{
    int rows=this->size();
    int columns(0);
    columns=((*this)[0]).size();
    for (int i=0;i<rows;i++)
    {
        for (int j=0;j<columns;j++)
        {
   //         cout << setw(7) << (*this)[i][j] << "  ";
        }
   //     cout << endl;
    }
}
