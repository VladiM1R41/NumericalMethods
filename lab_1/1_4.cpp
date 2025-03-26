#include <bits/stdc++.h>

using namespace std;

double EPS = 0.0001;

double sumUpDiagonal(const vector<vector<double>>& matrix) {
    double sum = 0.0;
    int n = matrix.size();
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            sum += matrix[i][j] * matrix[i][j];
    return sqrt(sum);
}

pair<int, int> maxNonDiagonalElement(const vector<vector<double>>& A) {
    double maxElement = -__DBL_MAX__;
    pair<int, int> maxIndex = make_pair(0, 0);
    for (int i = 0; i < A.size() - 1; i++)
        for (int j = i + 1; j < A.size(); j++)
            if (abs(A[i][j]) > maxElement) {
                maxElement = abs(A[i][j]);
                maxIndex = make_pair(i, j);
            }
    return maxIndex;
}

double Get_Phi(int max_i, int max_j, const vector<vector<double>>& A) {
    if(A[max_i][max_i] == A[max_j][max_j])
        return M_PI / 4;
    else
        return 0.5 * atan(2 * A[max_i][max_j] / (A[max_i][max_i] - A[max_j][max_j]));
}

vector<vector<double>> Transpose_Matrix(const vector<vector<double>>& A) {
    vector<vector<double>> result(A[0].size(), vector<double>(A.size(), 0));
    for (int i = 0; i < A.size(); i++)
        for (int j = 0; j < A[i].size(); j++)
            result[j][i] = A[i][j];
    return result;
}

vector<vector<double>> Matrix_Multiplication(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int n = A.size();
    int m = B[0].size();
    int p = B.size();

    vector<vector<double>> C(n, vector<double>(m, 0.0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int k = 0; k < p; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

vector<vector<double>> Initialize_U(const vector<vector<double>> A)
{
    int im = maxNonDiagonalElement(A).first;
    int jm = maxNonDiagonalElement(A).second;
    double phi = Get_Phi(im, jm, A);

    vector<vector<double>> U(A.size(), vector<double> (A.size(), 0));
    for (int i = 0; i < U.size(); i++)
        U[i][i] = 1;
    U[im][jm] = -sin(phi);
    U[jm][im] = sin(phi);
    U[im][im] = U[jm][jm] = cos(phi);
    return U;
}

void Checking_Results(vector<vector<double>> A, vector<double> lambda, vector<vector<double>> V)
{
    cout << "Proverca resultatov" << endl;
    int n = lambda.size();

    for(int num = 0; num < n; num++) {
        cout << endl << "nomer sobstvennogo znacheniya:" << num << endl;

        vector<double> V_j (n, 0);
        for(int i = 0; i < n; i++)
            V_j[i] = V[i][num];

        vector<double> AV_j (n, 0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                AV_j[i] += A[i][j] * V_j[j];

        vector<double> lambdaV_j (n, 0);
        for (int i = 0; i < n; i++)
            lambdaV_j[i] += lambda[num] * V_j[i];
        
        for(int i = 0; i < n; i++)
            cout << "\t" << AV_j[i] - lambdaV_j[i] <<endl;
    }
}

void Jacobi_Eigenvalue(vector<vector<double>> A)
{
    vector<vector<double>> A_k = A;
    vector<vector<double>> V(A.size(), vector<double> (A.size(), 0));
    for(int i = 0; i < V.size(); i++)
        V[i][i] = 1;
    while(sumUpDiagonal(A_k) > EPS)
    {
        vector<vector<double>> U = Initialize_U(A_k);
        vector<vector<double>> U_t = Transpose_Matrix(U);
        A_k = Matrix_Multiplication(Matrix_Multiplication(U_t, A_k), U);
        V = Matrix_Multiplication(V, U);
    }
    
    vector<double> lambda (A_k.size());
    for(int i = 0; i < A_k.size(); i++)
        lambda[i] = A_k[i][i];
    
    cout << "Sobstvennie znacheniya:" << endl;
    for(int i = 0; i < lambda.size(); i++)
        cout << "\t" << lambda[i] << endl;
    cout << "Sobstvennie vectora:" << endl;
    for(int j = 0; j < V.size(); j++){
        cout << j << ":" << endl;
        for(int i = 0; i < V.size(); i++)
            cout << "\t" << V[i][j] << endl;
    }

    Checking_Results(A, lambda, V);
}

int main()
{
    vector<vector<double>> A = 
    {
        {7, 3, -1},
        {3, -7, -8},
        {-1, -8, -2}
    };
    Jacobi_Eigenvalue(A);
    return 0;
}