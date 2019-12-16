#include <iostream>
#include <cmath>
double R0 = 0.257453;
double z0 = 0.314687;
const double v = -0.5;
const long double u = 5*M_PI/4;
const double r0 = hypot(R0,z0);
const double r03 = pow(r0,3);
const double Q0 = 1 - (pow(((2*v/R0)+(R0/r03)),2));
const double dR0 =(sqrt(Q0))*cos(u);
const double dz0 = (sqrt(Q0))*sin(u);
double h = 0.001;
double range = 0.3;
double ddR(double R,double z);
double ddz (double R,double z);
double dphi (double R, double z, double s);
void NumMethod(double a,double b, double da, double db,double h);
double testNumMethod(double a, double da, double h, double s);
double testRk4(double a, double da, double h, double s);
int main(void)
{
    std::cout.precision(15); std::cout.setf(std::ios::scientific);
    NumMethod(R0, z0, dR0, dz0, h);
    return 0;
}
double ddR(double R,double z)
{
  double r=hypot(R,z);
  double dderiv=(((2*v)/R)+(R/(pow(r,3))))*((2*v/pow(R,2))+(3*R*R/(pow(r,5)))-(1/pow(r,3)));
  return dderiv;
}
double ddz (double R ,double z)
{
  double r = std::hypot(R,z);
  double dderiv=(((2*v)/R)+(R/(pow(r,3))))*(3*R*z/(pow(r,5)));
  return dderiv;
}
double dphi (double R, double z, double s)
{
    double r = std::hypot(R,z);
    double deriv = (((2*v)/R)+(R/(pow(r,3))))/R;
    return deriv;
}
void NumMethod(double a,double b, double da, double db,double h)
{
    double nt = (range - 0.0)/h;
    double Ri = 0.257453;
    double zi = 0.314687;
    double phi0 = 0.0;
    double phi= 0.0;
    double yR = a - h*da + (h*h*0.5)*ddR(a,b);
    double yz = b - h*db + (h*h*0.5)*ddz(a,b);
    std::cout << "            " << "S" << "                      " << "R" << "                         " << "Z" << "                     " << "Phi" << "                        " << "X" << "                        " << "Y" << "\n";
    for (int n = 0; n<=nt; ++n)
    {
        double s = 0.0 + n*h;
        if (s ==0)
            {
                std::cout << s << "  " << Ri << "    " << zi << "    " << phi0 << "    " << Ri*std::cos(phi0) << "    " << Ri*std::sin(phi0) <<"\n";
            }
        else
        {
        std::cout << s << " " << "\t";
        //Se aplica el metodo de Stormer
    double R = (h*h)*ddR(R0,z0) + (2*R0) - yR;
    double z = (h*h)*ddz(R0,z0) + (2*z0) - yz;
    //Se integra dphi dando uso al met+odo de Runge-Kutta
    double k1 = h*dphi(R0,z0, s);
    double k2 = h*dphi((R0+(k1*0.5)), (z0 + (k1*0.5)), (s + (h*0.5)));
    double k3 = h*dphi((R0+(k2*0.5)), (z0 + (k2*0.5)), (s + (h*0.5)));
    double k4 = h*dphi((R0+k3),(z0+k3), (s+h));
    phi += (((k1 + (2*k2) + (2*k3) + k4))/6.0);
    yR = R0;
    yz = z0;
    R0 = R;
    z0 = z;
    //En Coordenadas Cartesianas
    double x = R*std::cos(phi);
    double y = R*std::sin(phi);
    std::cout << R << "    " << z << "    " << phi << "    " << x << "    " << y << "\n";
        }
    }
}
