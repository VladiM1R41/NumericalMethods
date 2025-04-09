#include <bits/stdc++.h>
using namespace std;
using Complex = complex<double>;
using Matrix = vector<vector<Complex>>;


Matrix householder_transformation(const vector<Complex>& a) {
    int n = a.size();
    vector<Complex> v = a;
    double norm_a = 0.0;
    for (const auto& val : a) {
        norm_a += norm(val);
    }
    norm_a = sqrt(norm_a);
    v[0] += (a[0].real() >= 0 ? 1 : -1) * norm_a;
    double norm_v = 0.0;
    for (const auto& val : v) {
        norm_v += norm(val);
    }
    norm_v = sqrt(norm_v);
    for (auto& val : v) {
        val /= norm_v;
    }
    Matrix H(n, vector<Complex>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        H[i][i] = 1.0;
        for (int j = 0; j < n; ++j) {
            H[i][j] -= 2.0 * v[i] * conj(v[j]);
        }
    }
    return H;
}

pair<Matrix, Matrix> qr_decomposition(const Matrix& A, int n) {
    Matrix Q(n, vector<Complex>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        Q[i][i] = 1.0;
    }
    Matrix R = A;
    for (int i = 0; i < n - 1; ++i) {
        vector<Complex> column(R.size() - i);
        for (int j = i; j < R.size(); ++j) {
            column[j - i] = R[j][i];
        }
        Matrix H = householder_transformation(column);
        Matrix temp(R.size() - i, vector<Complex>(R.size() - i, 0.0));
        for (int k = i; k < n; ++k) {
            for (int l = i; l < n; ++l) {
                Complex sum = 0.0;
                for (int m = i; m < n; ++m) {
                    sum += H[k - i][m - i] * R[m][l];
                }
                temp[k - i][l - i] = sum;
            }
        }
        
        for (int k = i; k < n; ++k) {
            for (int l = i; l < n; ++l) {
                R[k][l] = temp[k - i][l - i];
            }
        }
        
        temp = Matrix(n, vector<Complex>(n - i, 0.0));
        for (int k = 0; k < n; ++k) {
            for (int l = i; l < n; ++l) {
                Complex sum = 0.0;
                for (int m = i; m < n; ++m) {
                    sum += Q[k][m] * conj(H[m - i][l - i]);
                }
                temp[k][l - i] = sum;
            }
        }
        for (int k = 0; k < n; ++k) {
            for (int l = i; l < n; ++l) {
                Q[k][l] = temp[k][l - i];
            }
        }
    }
    return {Q, R};
}

pair<Complex, Complex> complex_solve(Complex a11, Complex a12, Complex a21, Complex a22, double eps) {
    Complex a = 1.0;
    Complex b = -a11 - a22;
    Complex c = a11 * a22 - a12 * a21;
    Complex d = b * b - 4.0 * a * c;
    if (d.real() > eps) {
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    }
    Complex d_c = Complex(0, sqrt(-d.real()));
    Complex x1 = (-b + d_c) / (2.0 * a);
    Complex x2 = (-b - d_c) / (2.0 * a);
    
    return {x1, x2};
}

pair<vector<Complex>, int> qr_eigenvalues(const Matrix& A, double e, int n) {
    Matrix Ak = A;
    vector<Complex> eig_vals(n, 0.0);
    int iter = 0;
    while (true) {
        iter++;
        auto [Q, R] = qr_decomposition(Ak, n);
        
        Matrix temp(n, vector<Complex>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    temp[i][j] += R[i][k] * Q[k][j];
                }
            }
        }
        Ak = temp;
        
        bool conv = true;
        int i = 0;
        
        while (i < n) {
            if (i < n - 1 && abs(Ak[i + 1][i]) > e) {
                auto [eig1, eig2] = complex_solve(Ak[i][i], Ak[i][i + 1], Ak[i + 1][i], Ak[i + 1][i + 1], e);
                
                if (!isnan(eig1.real()) && !isnan(eig2.real())) {
                    eig_vals[i] = eig1;
                    eig_vals[i + 1] = eig2;
                    i += 1;
                } else {
                    conv = false;
                }
            } else {
                eig_vals[i] = Ak[i][i];
            }
            i += 1;
        }
        if (conv) {
            break;
        }
    }
    
    return {eig_vals, iter};
}

void printMatrix(const Matrix& m) {
    for (const auto& row : m) {
        for (const auto& val : row) {
            cout << fixed << setprecision(2) << val << " ";
        }
        cout << endl;
    }
}

int main() {
    double e;
    int n;
    
    cout << "Enter precision (e): ";
    cin >> e;
    cout << "Enter matrix size (n): ";
    cin >> n;
    
    Matrix A(n, vector<Complex>(n));
    cout << "Enter matrix elements (row-wise):" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double val;
            cin >> val;
            A[i][j] = val;
        }
    }
    
    auto [Q, R] = qr_decomposition(A, n);
    auto [h, iter] = qr_eigenvalues(A, e, n);
    
    cout << "\nEigenvalues: ";
    for (const auto& val : h) {
        cout << val << " ";
    }
    cout << endl;
    
    cout << "\nQ:" << endl;
    printMatrix(Q);
    
    cout << "\nR:" << endl;
    printMatrix(R);
    
    Matrix QR(n, vector<Complex>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                QR[i][j] += Q[i][k] * R[k][j];
            }
        }
    }
    
    cout << "\nQ * R:" << endl;
    printMatrix(QR);
    
    cout << "\nNumber of iterations: " << iter << endl;
    
    return 0;
}   