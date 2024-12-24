#pragma once
#include <cmath>

double uxy(double x, double y)
{
    return (1 - pow((x - 1), 2) - pow((y - 0.5), 2));
}
double fxy(double x, double y)
{
    return -(6 * y + 2);
}
double uy0(double y)
{
    return (pow(y, 3) + 1);
}
double uy1(double y)
{
    return (pow(y, 3) + 2);
}
double ux0(double x)
{
    return (pow(x, 2) + 1);
}
double ux1(double x)
{
    return (pow(x, 2) + 2);
}

//вектор правой части размерности (n - 1) на (m - 1)
void vectorf(double** f, int n, int m, double xleft, double ybot, double xright, double ytop)
{
    double h = (xright - xleft) / n;
    double k = (ytop - ybot) / m;
    for (int i = 1; i < n; i++)
    {
        double xi = xleft + i * h;
        for (int j = 1; j < m; j++)
        {
            double yi = ybot + j * k;
            double res = 0;
            if (j == 1)
            {
                res += (1 / pow(k, 2)) * ux0(xi);
            }
            else if (j == m - 1)
            {
                res += (1 / pow(k, 2)) * ux1(xi);
            }
            if (i == 1)
            {
                res += (1 / pow(h, 2)) * uy0(yi);
            }
            else if (i == n - 1)
            {
                res += (1 / pow(h, 2)) * uy1(yi);
            }
            f[i][j] = -fxy(xi, yi) - res;
        }
    }
}

//известные значения сетчатой функции v в граничных узлах
void boundaryConditions(double** v, int n, int m, double xleft, double ybot, double xright, double ytop)
{
    double h = (xright - xleft) / n;
    double k = (ytop - ybot) / m;
    for (int i = 0; i < n + 1; i++)
    {
        double xi = xleft + i * h;
        v[i][0] = ux0(xi);
        v[i][m] = ux1(xi);
    }
    for (int j = 1; j < m; j++)
    {
        double yi = ybot + j * k;
        v[0][j] = uy0(yi);
        v[n][j] = uy1(yi);
    }
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            v[i][j] = 0;
        }
    }
}

//метод Зейделя
double seidel(double** v, int n, int m, double xleft, double ybot, double xright, double ytop, int Nmax, double eps)
{
    boundaryConditions(v, n, m, xleft, ybot, xright, ytop);
    int S = 0;
    double eps_max;
    double eps_curr;
    double h2 = -pow(n / (xright - xleft), 2);
    double k2 = -pow(m / (ytop - ybot), 2);
    double a2 = -2 * (h2 + k2);
    double v_old;
    double v_new;
    while (true)
    {
        eps_max = 0;
        if (S >= Nmax)
        {
            break;
        }
        for (int i = 1; i < n; i++)
        {
            double xi = xleft + i * (xright - xleft) / n;
            for (int j = 1; j < m; j++)
            {
                double yi = ybot + j * (ytop - ybot) / m;
                v_old = v[i][j];
                v_new = (-(h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1])) + fxy(xi, yi)) / a2;
                eps_curr = fabs(v_old - v_new);
                if (eps_curr > eps_max)
                {
                    eps_max = eps_curr;
                }
                v[i][j] = v_new;
            }
        }
        S++;
        if (eps_max < eps)
        {
            break;
        }
    }
    return S;
}

//расчет нормы невязки
double residualNorm(double** v, int n, int m, double xleft, double ybot, double xright, double ytop)
{
    double** vector_f = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        vector_f[i] = new double[m + 1];
    }
    vectorf(vector_f, n, m, xleft, ybot, xright, ytop);
    
    double h2 = pow(n / (xright - xleft), 2);
    double k2 = pow(m / (ytop - ybot), 2);
    double A = -2 * (h2 + k2);
    
    double resNorm = 0;
    
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            double resNormComponent;
            double leftSide;
            
            if (j != 1 && j != m - 1)
            {
                if (i != 1 && i != n - 1)
                {
                    leftSide = k2 * v[i][j - 1] + h2 * v[i - 1][j] + A * v[i][j] + h2 * v[i + 1][j] + k2 * v[i][j + 1];
                }
                else if (i == 1)
                {
                    leftSide = k2 * v[i][j - 1] + A * v[i][j] + h2 * v[i + 1][j] + k2 * v[i][j + 1];
                }
                else if (i == n - 1)
                {
                    leftSide = k2 * v[i][j - 1] + h2 * v[i - 1][j] + A * v[i][j] + k2 * v[i][j + 1];
                }
            }
            else if (j == 1)
            {
                if (i == 1)
                {
                    leftSide = A * v[i][j] + h2 * v[i + 1][j] + k2 * v[i][j + 1];
                }
                else if (i != n - 1)
                {
                    leftSide = h2 * v[i - 1][j] + A * v[i][j] + h2 * v[i + 1][j] + k2 * v[i][j + 1];
                }
                else if (i == n - 1)
                {
                    leftSide = h2 * v[i - 1][j] + A * v[i][j] + k2 * v[i][j + 1];
                }
            }
            else if (j == m - 1)
            {
                if (i == 1)
                {
                    leftSide = k2 * v[i][j - 1] + A * v[i][j] + h2 * v[i + 1][j];
                }
                else if (i != n - 1)
                {
                    leftSide = k2 * v[i][j - 1] + h2 * v[i - 1][j] + A * v[i][j] + h2 * v[i + 1][j];
                }
                else if (i == n - 1)
                {
                    leftSide = k2 * v[i][j - 1] + h2 * v[i - 1][j] + A * v[i][j];
                }
            }
            resNormComponent = abs(leftSide - vector_f[i][j]);
            resNorm += pow(resNormComponent, 2);
        }
    }
    resNorm = sqrt(resNorm);
    return resNorm;
}

//норма вектора погрешности ||z|| = max|zi|, норма Чебышева
double errorStd(double** v, int n, int m, double xleft, double ybot, double xright, double ytop)
{
    double max = 0;
    for (int i = 1; i < n; i++)
    {
        double xi = xleft + i * (xright - xleft) / n;
        for (int j = 1; j < m; j++)
        {
            double yi = ybot + j * (ytop - ybot) / m;
            double curr = abs(uxy(xi, yi) - v[i][j]);
            if (curr > max)
            {
                max = curr;
            }
        }
    }
    return max;
}

