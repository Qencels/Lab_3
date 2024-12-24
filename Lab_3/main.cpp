#include <iostream>
#include "functions.h"
using namespace std;
int main() {
    setlocale(LC_ALL, "rus");

    int n = 4, m = 4;
    double a = 0.0, b = 2.0, c = 0.0, d = 1.0;

    int NMax;
    double eps;

    std::cout << " ���������� ������ ������� ��� ������� ������ ������� ��� ��������� ��������" << std::endl;
    std::cout << " �������: u(x,y) = 1 � (x � 1)^2 � (y � 0.5)^2, x: [0, 2], y: [0, 1]\n";
    std::cout << " ������� ������������ ����� �������� = ";
    std::cin >> NMax;
    std::cout << " ������� ��������� �������� = ";
    std::cin >> eps;
    
    double** v = new double* [n + 1];
    for (int i = 0; i < n + 1; i++) {
        v[i] = new double[m + 1];
    }

    std::cout << " ��� ������� ���������� ����� � ������� ������ ������� (Nmax = " << NMax << ", eps = " << eps << ")" << std::endl;
    std::cout << " �� S = " << seidel(v, n, m, a, c, b, d, NMax, eps) << " �������� �������� ��������� �������: " << std::endl;
    for (int j = 0; j < m + 1; j++)
    {
        for (int i = 0; i < n + 1; i++)
        {
            std::cout << " v[" << i << "]" << "[" << j << "] = " << v[i][j] << std::endl;
        }
        std::cout << " --------------------------" << std::endl;
    }
    std::cout << " ��������� ����� ������� = " << residualNorm(v, n, m, a, c, b, d) << std::endl;
    std::cout << " ������������ ����� ����������� = " << errorStd(v, n, m, a, c, b, d) << std::endl;
}
