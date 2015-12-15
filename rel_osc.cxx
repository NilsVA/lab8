// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx, double* k1, double* k2, double* k3, double* k4);
void bs(double* yn, double* y0, double dx, double* k1, double* k2, double* k3, double* k4,double& teta);
//--------------------
using namespace std;
//--------------------

int main(void)
{
	ofstream out("solution");
  const int dim = 2;
	double dx = 0.25,x=0;
	const double L = 100;
	double teta;
	double yn[dim];
  double k1[dim], k2[dim], k3[dim], k4[dim];
	
  out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
  for(double p0=0;p0<5;p0+=0.2)
  {
    double y0[dim] = {p0, 0};
    x=0;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx,k1,k2,k3,k4);
		
		if(y0[1] > 0 && yn[1] < 0)
		{
		  bs(yn,y0,dx,k1,k2,k3,k4,teta);
		  
		  break;
		}
		//bs(yn,y0,dx,k1,k2,k3,k4);
    for(int i=0; i<dim; i++) y0[i] = yn[i];
		//out << x << "\t" << y0[0] << "\t" << y0[1] << endl;
	}
	out << p0 << "\t" << x+teta*dx-dx << "\t" << y0[0] << "\t" << y0[1] << endl;
  }
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx, double* k1, double* k2, double* k3, double* k4)
{
	const int dim = 2;
	
		
  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	
    
}
//-------------------
// Lorenz model
void f(double* const y0, const double x)
{
	
	double y[2] = { y0[0], y0[1]};

  
	  
	y0[0] = y[1];
	y0[1] = (-1.)*y[0]/sqrt(1.+pow(y[0],2));
}
 void bs(double* yn, double* y0, double dx, double* k1, double* k2, double* k3, double* k4,double& teta)
 {
        
	double b1,b2,b3,b4;
	double x1=y0[1];
	double x2=yn[1];
	double tetal=0;
	double tetar=1;
	teta = (tetal+tetar)/2.;
	double x3=1;
	while(abs(x3) > 1e-8 && teta>0)
	{
	b1=teta-(3./2.)*pow(teta,2)+(2./3.)*pow(teta,3);
	b2=pow(teta,2)-(2./3.)*pow(teta,3);
	b3=b2;
	b4=(-1./2.)*pow(teta,2)+(2./3.)*pow(teta,3);
	x3=x1+dx*(k1[1]*b1+k2[1]*b2+k3[1]*b3+k4[1]*b4);
	//x3=(x1+x2)/2.;
	
	  if(x3<0)
	  {
	    tetar=teta;
	  }
	  else
	  {
	    tetal=teta;
	  }
	 teta=(tetal+tetar)/2.; 
	}
	yn[1]=x3;
	
 }