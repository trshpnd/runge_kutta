#include <stdio.h>
#include <math.h>

const double a = 1.02; 		// constantes baseadas 
const double b = 0.159;		// no câncer ósseo (osteossarcoma)
const double vo = 0.01; 	// volume inicial em cm^3

const double to = 0; 			// tempo em dias
const double tf = 60;
const double h = 0.1;			// step

double f(double);
double df(double, double);
double euler_explicito(double, double);
double runge_kutta2(double, double);
double runge_kutta4(double, double);

int main(void) {

	FILE *fp = fopen("dados_saida.dat", "w");
	FILE *fe = fopen("dados_erro.dat","w");

	double euler, rk2, rk4, analitico, t;
	euler = rk2 = rk4 = vo;

	for(t = to; t < tf; t = t+h){
		euler = euler_explicito(t, euler);
		rk2 = runge_kutta2(t, rk2);
		rk4 = runge_kutta4(t, rk4);
		analitico = f(t);

		fprintf(fp, "%f %.4f %.4f %.4f %.4f\n", t, euler, rk2, rk4, analitico);
		fprintf(fe, "%f %.4f %.4f %.4f\n", t, fabs(euler-analitico)*100/analitico, fabs(rk2-analitico)*100/analitico, fabs(rk4-analitico)*100/analitico);

		//printf("%f %.2f %.2f %.2f %.2f\n", t, euler, rk2, rk4, analitico);
	}

	return 0;
}

double f(double t){
	return vo * (exp(a * (1 - exp(-b * t)) / b));
}

double df(double t, double v){ //dv/dt
	return a * v * exp(-b*t);
}

double euler_explicito(double t, double v){
	return v = v + df(t, v)*h;
}

double runge_kutta2(double t, double v){
	double k1 = df(t, v);
	double k2 = df(t + h, v + (k1 * h));
	return v + (h * (k1 + k2) / 2);
}

double runge_kutta4(double t, double v){
	double k1 = df(t, v);
	double k2 = df(t + (h/3), v + ((k1/3) * h));
	double k3 = df(t + (2 * h/3), v - (k1 * h/3) + (k2 * h));
	double k4 = df(t + h, v + (k1 * h) - (k2 * h) + (k3 * h));
	return v + (h * (k1 + (3 * k2) + (3 * k3) + k4) / 8);
}