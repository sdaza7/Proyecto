#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;
using fptr = double (double);

double oscilador (double w, double h, double y0, fptr ddx);
double ddx(double x);
const double M_PI = 3.14159265358979;
 const double w = M_PI;
 const double h = 0.001;
 double a = 3.21;

int main()
{
    cout << "t" << "         " << "X" << std::endl;
    oscilador(w, h, a, ddx);
    return 0;

}

double ddx(double x)
{
    double dderiv = -(w*w*x);
    return dderiv;
}

double oscilador (double w, double h, double y0, fptr ddx)
{
    double yb = 0; //yb es y_{-1}
    double y1 = 0;
    double yn = 0;
    yb = y0 + (h*h/2)*(ddx(y0));
    for (int ii=0; ii<=2000; ++ii)
    {
        cout << ii << "\t";
        if (ii == 0)
        {
            cout << y0 << "\n";
        }
        else
        {
        y1 = (h*h*ddx(y0)) + 2*y0 - yb;
        cout << y1 << "\n";
        yb = y0;
        y0 = y1;
        }
    }
}
