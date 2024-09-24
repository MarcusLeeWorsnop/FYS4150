#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main() {

    int N = 6;  

    mat A = zeros<mat>(N,N);

    for (int i = 0; i < N; i++) {
        A(i,i) = 2;  
        if (i < N - 1) {
            A(i+1,i) = -1;  
            A(i,i+1) = -1; 
        }
    }

    cout << "Matrix A: " << endl << A << endl;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);

    cout << "Numerical Eigenvalues and Eigenvectors from Armadillo:" << endl;
    for (int i = 0; i < N; i++) {
        cout << "Eigenvalue " << i+1 << ": " << eigval[i] << endl;
        cout << "Corresponding Eigenvector: " << endl;
        cout << eigvec.col(i).t() << endl;  
    }

    double pi = M_PI;

    cout << "Analytical Eigenvalues and Eigenvectors:" << endl;
    for (int j = 1; j <= N; j++) {  
        
        double lambda = 2.0 + 2.0 * cos(j * pi / (N + 1));
        cout << "Analytical Eigenvalue " << j << ": " << lambda << endl;

        vec v(N);
        for (int k = 1; k <= N; k++) {
            v(k-1) = sin(k * j * pi / (N + 1));  
        }
        
        v = normalise(v);

        cout << "Analytical Eigenvector: "<< endl;
        cout << v.t() << endl;
    }

    return 0;
}
