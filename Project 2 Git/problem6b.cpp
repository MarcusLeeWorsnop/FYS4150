#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

double max_offdiag_symmetric(const mat& A, int& k, int& l) {
    double max_value = 0.0;  // Initialize to 0
    int N = A.n_rows;

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double abs_value = abs(A(i, j));
            if (abs_value > max_value) {
                max_value = abs_value;
                k = i;
                l = j;
            }
        }
    }
    return max_value;
}

void jacobi_rotate(mat& A, mat& R, int k, int l) {
    if (A(k, l) != 0) {
        double tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
        double t_theta;
        if (tau > 0) {
            t_theta = 1 / (tau + sqrt(1 + pow(tau, 2)));
        } else {
            t_theta = 1 / (tau - sqrt(1 + pow(tau, 2)));
        }

        double c_theta = 1 / sqrt(1 + pow(t_theta, 2));
        double s_theta = c_theta * t_theta;

        double a_kk = A(k, k);
        double a_ll = A(l, l);
        double a_kl = A(k, l);

        A(k, k) = a_kk * pow(c_theta, 2) - 2 * a_kl * c_theta * s_theta + a_ll * pow(s_theta, 2);
        A(l, l) = a_ll * pow(c_theta, 2) + 2 * a_kl * c_theta * s_theta + a_kk * pow(s_theta, 2);
        A(k, l) = 0;
        A(l, k) = 0;

        for (int i = 0; i < A.n_rows; i++) {
            if (i != k && i != l) {
                double a_ik = A(i, k);
                double a_il = A(i, l);
                A(i, k) = c_theta * a_ik - s_theta * a_il;
                A(k, i) = A(i, k);
                A(i, l) = c_theta * a_il + s_theta * a_ik;
                A(l, i) = A(i, l);
            }
        }

        for (int i = 0; i < R.n_rows; i++) {
            double r_ik = R(i, k);
            double r_il = R(i, l);
            R(i, k) = c_theta * r_ik - s_theta * r_il;
            R(i, l) = c_theta * r_il + s_theta * r_ik;
        }
    }
}

void jacobi_eigensolver(mat& A, vec& eigenvalues, mat& eigenvectors, double tol = 1e-8, int max_iterations = 100000000) {
    int N = A.n_rows;
    eigenvectors = eye<mat>(N, N);
    int k, l;
    double max_offdiag;
    int iterations = 0;

    do {
        max_offdiag = max_offdiag_symmetric(A, k, l);
        if (max_offdiag > tol) {
            jacobi_rotate(A, eigenvectors, k, l);
        }
        iterations++;
    } while (max_offdiag > tol && iterations < max_iterations);

    eigenvalues = A.diag();

    //cout << "Converged in " << iterations << " iterations.\n";
}


void columnextract(mat eigenvectors,int& column){
  for(int i = 0; i < eigenvectors.n_rows; i++){
    cout << eigenvectors(column * eigenvectors.n_rows + i) << endl;
  }
}


int main() {

    double n = 11; // Set this to 101 for N = 100;  setting n as int causes h to become 0
    double h = 1/n; // step size as (x_max - x_min)/n_steps
    int N = int(n) - 1;

    double a = -1/pow(h,2);
    double d = 2/pow(h,2);

    mat A = d * eye<mat>(N, N);
    for (int i = 0; i < N - 1; i++) {
        A(i, i+1) = A(i+1, i) = a;
    }

    //cout << "Original matrix A:\n" << A << endl;

    vec eigenvalues;
    mat eigenvectors;

    jacobi_eigensolver(A, eigenvalues, eigenvectors);

    eigenvectors = normalise(eigenvectors);

    int eigmin = (abs(eigenvalues)).index_min(); // lowest eigenvalue (lowest absolute val)

    vec eig2 = eigenvalues;
    eig2(eigmin) = (abs(eigenvalues)).max();
    int eigmin2 = (abs(eig2)).index_min();
    vec eig3 = eig2;
    eig3(eigmin2) = (abs(eigenvalues)).max();
    int eigmin3 = (abs(eig3)).index_min();


    cout << "Lowest eigenvalue (index: "<< eigmin <<"): " << eigenvalues(eigmin) << endl;
    cout << "Second lowest (index: "<< eigmin2 <<"): " << eig2(eigmin2) << endl;
    cout << "Third lowest (index: "<< eigmin3 <<"): " << eig3(eigmin3) << endl;

    //cout << eigenvalues << endl;

    //columnextract(eigenvectors,eigmin);
    //cout << endl;
    //columnextract(eigenvectors,eigmin2);
    //cout << endl;
    //columnextract(eigenvectors,eigmin3);

    cout << "\nEigenvalues:\n" << eigenvalues << endl;
    cout << "\nEigenvectors (columns of R):\n" << eigenvectors << endl;



    return 0;
}
