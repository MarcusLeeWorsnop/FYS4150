#include <iostream>
#include <armadillo>


using namespace std;
using namespace arma;

int main() {

	int N = 6;

	mat A = zeros<mat>(N,N);

	for (int i = 0; i < N; i++){
		A(i,i) = 2;
		if (i < N -1){
			A(i+1,i) = -1;
			A(i,i+1) = -1;
		}
	}

	cout<<"Matrix A: "<<endl<< A;

	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, A);

    for (int i = 0; i < N; i++) {
    cout << "Eigenvalue " << i+1 << ": " << eigval[i] << endl;
    cout << "Corresponding Eigenvector: " << endl;
    cout << eigvec.col(i) << endl;
}

    return 0;
}
