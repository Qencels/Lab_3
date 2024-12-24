#include <iostream>
#include "functions.h"
using namespace std;
int main() {
    setlocale(LC_ALL, "rus");

    int n = 3, m = 3;
    double a = 0.0, b = 1.0, c = 0.0, d = 0.2;

    int NMax;
    double eps;

    std::cout << " Реализация метода Зейделя для решения задачи Дирихле для уравнения Пуассона" << std::endl;
    std::cout << " Функция: u(x,y) = x^3 + y^3 + 2, x: [0, 1], y: [0, 0.2]\n";
    std::cout << " Укажите максимальное число итераций = ";
    std::cin >> NMax;
    std::cout << " Укажите требуемую точность = ";
    std::cin >> eps;
    
    double** v = new double* [n + 1];
    for (int i = 0; i < n + 1; i++) {
        v[i] = new double[m + 1];
    }

    std::cout << " При решении разностной схемы с помощью метода Зейделя (Nmax = " << NMax << ", eps = " << eps << ")" << std::endl;
    std::cout << " за S = " << seidel(v, n, m, a, c, b, d, NMax, eps) << " итераций получено численное решение: " << std::endl;
    for (int j = 0; j < m + 1; j++)
    {
        for (int i = 0; i < n + 1; i++)
        {
            std::cout << " v[" << i << "]" << "[" << j << "] = " << v[i][j] << std::endl;
        }
        std::cout << " --------------------------" << std::endl;
    }
    std::cout << " Евклидова норма невязки = " << residualNorm(v, n, m, a, c, b, d) << std::endl;
    std::cout << " Максимальная общая погрешность = " << errorStd(v, n, m, a, c, b, d) << std::endl;
}
