#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void printMatrix(vector<vector<double>>& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void printVector(vector<double>& vec, int n) {
    for (int i = 0; i < n; i++) {
        cout << vec[i] << "\t";
    }
    cout << endl;
}

void swap_rows(vector<vector<double>>& A, int i, int j) {
    vector<double> temp = A[i];
    A[i] = A[j];
    A[j] = temp;
}

void swap_el(vector<int>& P, int i, int j) {
    int temp = P[i];
    P[i] = P[j];
    P[j] = temp;
}

int findPivot(vector<vector<double>>& matrix, int col, int n) {
    int pivotIdx = col;
    double maxVal = fabs(matrix[col][col]);
    for (int i = col + 1; i < n; i++) {
        if (fabs(matrix[i][col]) > maxVal) {
            maxVal = fabs(matrix[i][col]);
            pivotIdx = i;
        }
    }
    return pivotIdx;
}

void luDecomposition(vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P, int n, int& swaps) {

    for (int i = 0; i < n; i++) {
        P[i] = i;
        L[i][i] = 1.0; 
    }

    for (int k = 0; k < n; k++) {
        int pivot = findPivot(A, k, n);
        if (pivot != k) {
            swaps++;
            swap_rows(A, k, pivot);
            swap_el(P, k, pivot);
            for (int i = 0; i < k; i++) {
                double temp = L[k][i];
                L[k][i] = L[pivot][i];
                L[pivot][i] = temp;
            }
        }

        U[k][k] = A[k][k];
        for (int i = k + 1; i < n; i++) {
            L[i][k] = A[i][k] / U[k][k];
            U[k][i] = A[k][i];
        }

        for (int i = k + 1; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
}

vector<double> solveSystem(vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P, vector<double>& b, int n) {
    vector<double> y(n, 0);
    vector<double> x(n, 0);

    for (int i = 0; i < n; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

double determinant(vector<vector<double>>& U, vector<int>& P, int n, int swaps) {
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        det *= U[i][i];
    }
    return det * pow(-1, swaps); 
}

vector<vector<double>> inverseMatrix(vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P, int n) {
    vector<vector<double>> inv(n, vector<double>(n, 0));
    vector<double> e(n, 0);


    for (int col = 0; col < n; col++) {
        e.assign(n, 0);
        e[col] = 1.0;
        vector<double> x = solveSystem(L, U, P, e, n);
        for (int i = 0; i < n; i++) {
            inv[i][col] = x[i];
        }
    }
    return inv;
}

int main() {
    int n;
    cout << "Enter the size of the matrix (n): ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n, 0));
    vector<double> b(n, 0);

    cout << "Enter the elements of matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "A[" << i << "][" << j << "]: ";
            cin >> A[i][j];
        }
    }

    cout << "Enter the elements of vector b:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "b[" << i << "]: ";
        cin >> b[i];
    }

    int swaps = 0;
    //vector<vector<double>> A = {{-1, -8, 0, 5}, {6, -6, 2, 4}, {9, -5, -6, 4}, {-5, 0, -9, 1}};
    //vector<double> b = {-60, -10, 65, 18}; 

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    vector<int> P(n);

    cout << "Matrix A:" << endl;
    printMatrix(A, n);

    luDecomposition(A, L, U, P, n,swaps);

    cout << "Matrix L:" << endl;
    printMatrix(L, n);
    cout << "Matrix U:" << endl;
    printMatrix(U, n);
    vector<double> x = solveSystem(L, U, P, b, n);
    cout << "x:" << endl;
    printVector(x, n);

    double det = determinant(U, P, n, swaps);
    cout << "Determinant: " << det << endl;

    vector<vector<double>> inv = inverseMatrix(L, U, P, n);
    cout << "Inverse Matrix:" << endl;
    printMatrix(inv, n);

    return 0;
}