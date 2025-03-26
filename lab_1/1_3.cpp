#include <bits/stdc++.h>
using namespace std;

const double INF = 1e9;

double c_vec_norm(const vector<double>& vec) {
    double norm = 0;
    for (double val : vec) {
        norm = max(norm, abs(val));
    }
    return norm;
}

double c_mat_norm(const vector<vector<double>>& mat) {
    double norm = 0;
    for (const auto& row : mat) {
        double row_sum = 0;
        for (double val : row) {
            row_sum += abs(val);
        }
        norm = max(norm, row_sum);
    }
    return norm;
}

pair<vector<double>, int> simple_iterations(const vector<vector<double>>& A, 
                                         const vector<double>& b, 
                                         double epsilon, int n) {
    vector<vector<double>> alpha(n, vector<double>(n, 0));
    vector<double> beta(n);
    vector<double> x_prev(n);
    
    for (int i = 0; i < n; ++i) {
        beta[i] = b[i] / A[i][i];
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
    }
    
    x_prev = beta;
    int iterations = 0;
    double error = epsilon + 1; 
    
    while (error > epsilon) {
        vector<double> x_current(n);
        
        for (int i = 0; i < n; ++i) {
            x_current[i] = beta[i];
            for (int j = 0; j < n; ++j) {
                x_current[i] += alpha[i][j] * x_prev[j];
            }
        }
        
        vector<double> diff(n);
        for (int i = 0; i < n; ++i) {
            diff[i] = x_current[i] - x_prev[i];
        }
        
        error = (c_mat_norm(alpha) / (1 - c_mat_norm(alpha))) * c_vec_norm(diff);
        x_prev = x_current;
        iterations++;
    }
    
    return make_pair(x_prev, iterations);
}

pair<vector<double>, int> seidel_method(const vector<vector<double>>& A, 
                                        const vector<double>& b, 
                                        double epsilon, int n) {
    vector<vector<double>> alpha(n, vector<double>(n, 0));
    vector<double> beta(n);
    vector<double> x_prev(n);

    for (int i = 0; i < n; ++i) {
        beta[i] = b[i] / A[i][i];
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
    }

    x_prev = beta;
    int iterations = 0;
    double error = epsilon + 1;

    while (error > epsilon) {
        vector<double> x_current = x_prev;

        x_current[0] = beta[0];
        for (int j = 0; j < n; ++j) {
            x_current[0] += alpha[0][j] * x_prev[j];
        }

        for (int i = 1; i < n; ++i) {
            x_current[i] = beta[i];
            for (int j = 0; j < i; ++j) {
                x_current[i] += alpha[i][j] * x_current[j];
            }
            for (int j = i; j < n; ++j) {
                x_current[i] += alpha[i][j] * x_prev[j];
            }
        }

        vector<double> diff(n);
        for (int i = 0; i < n; ++i) {
            diff[i] = x_current[i] - x_prev[i];
        }

        error = (c_mat_norm(alpha) / (1 - c_mat_norm(alpha))) * c_vec_norm(diff);
        x_prev = x_current;
        iterations++;
    }

    return make_pair(x_prev, iterations);
}

int main() {
    int n = 4;
    vector<vector<double>> A = {{14, -4, -2, 3}, {-3, 23, -6, -9}, {-7, -8, 21, -5}, {-2, -2, 8, 18}};
    vector<double> b = {38, -195, -27, 142};
    double epsilon = 1e-6;

    auto [x_simple, iter_simple] = simple_iterations(A, b, epsilon, n);
    auto [x_seidel, iter_seidel] = seidel_method(A, b, epsilon, n);
    
    cout << "Method of simple iterations:" << endl;
    cout << "Solution: ";
    for (double val : x_simple) {
        cout << val << " ";
    }
    cout << "\nNumber of iterations: " << iter_simple << endl;
    
    cout << "\nSeidel method:" << endl;
    cout << "Solution: ";
    for (double val : x_seidel) {
        cout << val << " ";
    }
    cout << "\nNumber of iterations: " << iter_seidel << endl;
    
    return 0;
}