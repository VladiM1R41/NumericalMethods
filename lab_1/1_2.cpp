#include <iostream>
#include <vector>


void extractDiagonals(const std::vector<std::vector<double>>& matrix, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
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

void algorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& x) {
    size_t n = d.size();
    std::vector<double> P(n, 0.0); 
    std::vector<double> Q(n, 0.0); 

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
    std::cout << "Enter the size of the matrix (n): ";
    std::cin >> n;

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));
    std::cout << "Enter the elements of matrix A:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << "A[" << i << "][" << j << "]: ";
            std::cin >> matrix[i][j];
        }
    }

    std::vector<double> d(n, 0.0);
    std::cout << "Enter the elements of vector b:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        std::cout << "b[" << i << "]: ";
        std::cin >> d[i];
    }

    std::vector<double> a, b, c;
    extractDiagonals(matrix, a, b, c);

    std::vector<double> x(n, 0.0);
    algorithm(a, b, c, d, x);

    std::cout << "Решение системы:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}