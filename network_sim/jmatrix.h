#ifndef JMATRIX_H
#define JMATRIX_H

#include <vector>
#include <string>


using namespace std;

class JMatrix : public vector< vector<double> > { //Coupling matrix
public:
JMatrix();
//JMatrix(int N, int N_input, string fn);
JMatrix(int N, int N_input, double jl, double jr);
JMatrix(int N, int N_input, double jl);
JMatrix(int N, int N_input, char* file);
void output();
};


#endif // JMATRIX_H
