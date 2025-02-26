#include <iostream>
#include <vector>

using namespace std;

void extractDiagonals(const vector<vector<double>>& matrix, vector<double>& a, vector<double>& b, vector<double>& c) {
    size_t n = matrix.size();
    a.resize(n);
    b.resize(n);
    c.resize(n);

    for (size_t i = 0; i < n; ++i) {
        b[i] = matrix[i][i];
        if (i > 0) a[i] = matrix[i][i - 1];
        if (i < n - 1) c[i] = matrix[i][i + 1];
    }
}

void algorithm(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d, vector<double>& x) {
    size_t n = d.size();
    vector<double> P(n, 0.0); 
    vector<double> Q(n, 0.0); 

    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    for (size_t i = 1; i < n; ++i) {
        double denominator = b[i] + a[i] * P[i - 1];
        P[i] = -c[i] / denominator;
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denominator;
    }

    x[n - 1] = Q[n - 1];

    for (int i = n - 2; i >= 0; --i) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }
}


int main() {
    size_t n;
    cout << "Enter the size of the matrix (n): ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n, 0.0));
    cout << "Enter the elements of matrix A:" << endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            cout << "A[" << i << "][" << j << "]: ";
            cin >> matrix[i][j];
        }
    }

    vector<double> d(n, 0.0);
    cout << "Enter the elements of vector b:" << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << "b[" << i << "]: ";
        cin >> d[i];
    }

    vector<double> a, b, c;
    extractDiagonals(matrix, a, b, c);

    vector<double> x(n, 0.0);
    algorithm(a, b, c, d, x);

    cout << "Решение системы:" << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}