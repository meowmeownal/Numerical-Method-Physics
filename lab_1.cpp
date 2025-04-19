#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <Eigen/Eigenvalues>

#define N2 200

void siatka(std::vector<double>& xk, std::vector<double>& yk, double a, int n) {
    double delta_x = a * 2.0 / (n - 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int k = i * n + j;
            xk[k] = -a + delta_x * i;
            yk[k] = -a + delta_x * j;
        }
    }
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& vec)
 {
    for (const auto& row : vec)
     {
        for (size_t i = 0; i < row.size(); i++)
         {
            os << row[i];
            if (i < row.size() - 1) os << ", "; 
        }
        os << "\n"; 
    }
    return os;
}

double gauss(double x, double y, double xk, double yk, double alfa_x, double alfa_y) {
    return (1.0 / std::pow(alfa_x * M_PI, 0.25)) * exp( -(x - xk) * (x - xk) / (alfa_x * 2) ) *
           (1.0 / std::pow(alfa_y * M_PI, 0.25)) * exp( -(y - yk) * (y - yk) / (alfa_y * 2) );
}

double overlap(double xk, double yk, double xl, double yl, double alpha_x, double alpha_y) {
    return exp(-((xk - xl) * (xk - xl)) / (4 * alpha_x) - ((yk - yl) * (yk - yl)) / (4 * alpha_y));
}

double K(double xk, double yk, double xl, double yl, double alpha_x, double alpha_y, double m) {
    return -1.0 / (2.0 * m) * ( ((xk - xl) * (xk - xl) - 2.0 * alpha_x) / (4.0*alpha_x*alpha_x) +
                                ((yk - yl) * (yk - yl) - 2.0 * alpha_y) / (4.0 * alpha_y * alpha_y) );
}

double V(double xk, double yk, double xl, double yl, double m, double wx, double wy) {
    return m / 2.0* (wx * wx * ((xk + xl) * (xk + xl) + 2.0 / (wx * m)) / 4.0 +
                    wy * wy * ((yk + yl) * (yk + yl) + 2.0 / (wy * m)) / 4.0);
}

int main() {
    int n = 9;
    const double Eh = 27.211;
    const double wx = 0.08 / Eh;
    const double wy = 0.2 / Eh;
    const double m = 0.24;
    const double delta_x = 2.0 / 0.0529;
    const double delta_x2 = 1.0 / 0.0529;
    const double a = delta_x * (n - 1) / 2.0;
    const double a2 = delta_x2 * (n - 1) / 2.0;
    
    std::vector<double> xk1(n * n, 0), yk1(n * n, 0);
    std::vector<double> xk2(N2 * N2, 0), yk2(N2 * N2, 0);

    
    siatka(xk1, yk1, a, n);
    siatka(xk2, yk2, a, N2);
    
    std::vector<int> k_val{0, 8, 9};
    size_t t = 0;
    
    for (int k = 0; k < 10; k++) {
        if (t < k_val.size() && k == k_val[t]) 
        {
            std::ofstream file1("k" + std::to_string(k) + ".csv");
            for (int i = 0; i < N2; i++) 
            {
                for (int j = 0; j < N2; j++) 
                {
                    file1 << gauss(xk2[i*N2 + j], yk2[i*N2 + j], xk1[k], yk1[k], 1.0 / (wx * m), 1.0 / (wy * m));
                    if (j < N2 - 1) file1 << ",";
                }
                file1 << "\n";
            }
            t++;
        }
    }
    
    Eigen::MatrixXd H_matrix(n * n, n * n);
    Eigen::MatrixXd S_matrix(n * n, n * n);
    Eigen::MatrixXd K_matrix(n * n, n * n);
    Eigen::MatrixXd V_matrix(n * n, n * n);
    
    H_matrix.setZero();
    S_matrix.setZero();
    K_matrix.setZero();
    V_matrix.setZero();

    std::vector<double> xk3(n * n, 0), yk3(n * n, 0);
    std::vector<double> xk4(N2 * N2, 0), yk4(N2 * N2, 0);
    siatka(xk3, yk3, a2, n);
    siatka(xk4, yk4, a2, N2);
    
    for (int k = 0; k < n * n; k++) 
    {
        for (int l = 0; l < n * n; l++) 
        {
            S_matrix(k, l) = overlap(xk3[k], yk3[k], xk3[l], yk3[l], 1.0 / (wx * m), 1.0 / (wy * m));
            K_matrix(k, l) = K(xk3[k], yk3[k], xk3[l], yk3[l], 1.0 / (wx * m), 1.0 / (wy * m), m) * S_matrix(k, l);
            V_matrix(k, l) = V(xk3[k], yk3[k], xk3[l], yk3[l], m, wx, wy) * S_matrix(k, l);
            H_matrix(k, l) = K_matrix(k, l) + V_matrix(k, l);
        }
    }
    
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_matrix, S_matrix);
    Eigen::VectorXd wartosci_wlasne = solver.eigenvalues();
    Eigen::MatrixXd wektory_wlasne = solver.eigenvectors();
    Eigen::MatrixXd psi_squared = Eigen::MatrixXd::Zero(N2, N2);

    for (int e = 0; e < 6; ++e) 
    {
        std::ofstream file4("./psi_squared" + std::to_string(e) + ".csv");
        for (int i = 0; i < N2; i++) 
        {
            for (int j = 0; j < N2; j++)
            {
                int k = i * N2 + j;

                double psi = 0.0;
                for (int l = 0; l < n * n; l++)
                 {
                    double alpha_x = 1.0 / (wx * m);
                    double alpha_y =  1.0 / (wy * m);
    
                    psi += wektory_wlasne(l, e) * gauss(xk4[k], yk4[k], xk3[l], yk3[l], alpha_x, alpha_y);
                }
                psi_squared(i, j) = psi * psi; 

            }
        }
        file4 << psi_squared;
        file4.close();
    }


    std::vector<std::vector<double>> energy;

    Eigen::MatrixXd H_matrix2(n * n, n * n);
    Eigen::MatrixXd S_matrix2(n * n, n * n);
    Eigen::MatrixXd K_matrix2(n * n, n * n);
    Eigen::MatrixXd V_matrix2(n * n, n * n);
    
    H_matrix2.setZero();
    S_matrix2.setZero();
    K_matrix2.setZero();
    V_matrix2.setZero();

    std::ofstream file4("./energy.csv");
        for (double wx2 = 0.1 / Eh; wx2 <= 0.51 / Eh; wx2 += 0.005 / Eh)  
        {
            for (int k = 0; k < n * n; k++) 
            {
                for (int l = 0; l < n * n; l++) 
                {
                    S_matrix2(k, l) = overlap(xk3[k], yk3[k], xk3[l], yk3[l], 1.0 / (wx2 * m), 1.0 / (wy * m));
                    K_matrix2(k, l) = K(xk3[k], yk3[k], xk3[l], yk3[l], 1.0 / (wx2 * m), 1.0 / (wy * m), m) * S_matrix2(k, l);
                    V_matrix2(k, l) = V(xk3[k], yk3[k], xk3[l], yk3[l], m, wx2, wy) * S_matrix2(k, l);
                    H_matrix2(k, l) = K_matrix2(k, l) + V_matrix2(k, l);
                }
            }
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_matrix2, S_matrix2);
            Eigen::VectorXd wartosci_wlasne2 = solver.eigenvalues();
            Eigen::MatrixXd wektory_wlasne2 = solver.eigenvectors();
            

            //energy.clear();
            file4 << wx2 * Eh; 
            for (int i = 0; i < 10; i++)
            {
                file4 << ", " << wartosci_wlasne2(i) * Eh;  
            }
            file4 << "\n";  
        }

    
    return 0;
}