#include <iostream>
#include "functions.h"
using namespace std;
int main() {
    setlocale(LC_ALL, "rus");

    int n = 4, m = 4;
    double a = 0.0, b = 2.0, c = 0.0, d = 1.0;

    int NMax;
    double eps;

    std::cout << " –еализаци€ метода «ейдел€ дл€ решени€ задачи ƒирихле дл€ уравнени€ ѕуассона" << std::endl;
    std::cout << " ‘ункци€: u(x,y) = 1 Ц (x Ц 1)^2 Ц (y Ц 0.5)^2, x: [0, 2], y: [0, 1]\n";
    std::cout << " ”кажите максимальное число итераций = ";
    std::cin >> NMax;
    std::cout << " ”кажите требуемую точность = ";
    std::cin >> eps;
    
    double** v = new double* [n + 1];
    for (int i = 0; i < n + 1; i++) {
        v[i] = new double[m + 1];
    }

    std::cout << " ѕри решении разностной схемы с помощью метода «ейдел€ (Nmax = " << NMax << ", eps = " << eps << ")" << std::endl;
    std::cout << " за S = " << seidel(v, n, m, a, c, b, d, NMax, eps) << " итераций получено численное решение: " << std::endl;
    for (int j = 0; j < m + 1; j++)
    {
        for (int i = 0; i < n + 1; i++)
        {
            std::cout << " v[" << i << "]" << "[" << j << "] = " << v[i][j] << std::endl;
        }
        std::cout << " --------------------------" << std::endl;
    }
    std::cout << " ≈вклидова норма нев€зки = " << residualNorm(v, n, m, a, c, b, d) << std::endl;
    std::cout << " ћаксимальна€ обща€ погрешность = " << errorStd(v, n, m, a, c, b, d) << std::endl;
}
