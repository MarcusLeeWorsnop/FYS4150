#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


double max_offdiag_symmetric(const mat& A, int& k, int &l){
	double max_value;
	int N = A.n_rows;

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
        	if (i != j){
        		double abs_value = abs(A(i,j));
				if (abs_value > max_value){
					max_value = abs_value;
					k=i;
					l=j;
				}

        	}
        }
    }

    return max_value;
}


int main(){

	mat A = {{1,0,0,0.5},
			 {0,1,-0.7,0},
			 {0,-0.7,1,0},
			 {0.5,0,0,1}};

	int k = 0, l = 0;

	double max_val = max_offdiag_symmetric(A,k,l);

	cout << "Largest off-diagonal element: " << max_val << endl;
	cout << "Indicies of largest off-diagonal element: (" <<k+1<<","<<l+1<<")"<<endl;
	return 0;
}


