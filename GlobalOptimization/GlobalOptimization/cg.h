#ifndef CG_H_INCLUDED
#define CG_H_INCLUDED

#include "report.h"
#include "testfun.h"

static const unsigned seed_p = std::chrono::system_clock::now().time_since_epoch().count();
static mt19937 gen_p(seed_p);

double inner_prod_dv(const int &n, point &dv1, point &dv2)
{
	double res = 0;
	for (int i = 0; i < n; ++i)
		res += dv1[i] * dv2[i];
	return res;
}

point mult_scalar_dv(const int &n, const double &alpha, point &dv)
{
	point res(n);
	for (int i = 0; i < n; ++i)
		res[i] = alpha * dv[i];
	return res;
}

point sum_dv(const int &n, point &dv1, point &dv2)
{
	point res(n);
	for (int i = 0; i < n; ++i)
		res[i] = dv1[i] + dv2[i];
	return res;
}

double GoldenSectionSearch(const int &fcode, const int &n, point &s, point &p)
{
	const double EPS = 1e-8;
	double a = 0;
	double b = 1e5;      
	double x0 = a + 0.5 * (3 - sqrt(5.0)) * (b - a);
	double x1 = b - x0 + a;

	while (abs(b - a) > EPS)
	{
		if (testfun::f(n, sum_dv(n, s, mult_scalar_dv(n, x0, p)), fcode) < testfun::f(n, sum_dv(n, s, mult_scalar_dv(n, x1, p)), fcode))
			b = x1;
		else
			a = x0;

		x1 = x0;
		x0 = b + a - x1;
	}
	return (a + b) / 2;
}

FullReport cg_FletcherRieves(const int &fcode, const int &n, const point &st, const int &max_call_f)
{
	FullReport fr;
	testfun::call_f = 0;
	uniform_real_distribution<double> distr_p(testfun::lb[fcode], testfun::rb[fcode]);
	const double eps = 1e-5;
	point x = st;
	double curVal = testfun::f(n, x, fcode);
	Report r(testfun::call_f, curVal);
	fr.insert_into_report(r);
	double prevVal = curVal;
	point antiGrad = mult_scalar_dv(n, -1.0, testfun::grad_f(n, x, fcode));	//antiGrad = (-1) * Grad 
	double gradSquare = inner_prod_dv(n, antiGrad, antiGrad);
	int numIter = 0;
	double alpha, beta, newGradSquare;
	point newGrad;
	while (testfun::call_f < max_call_f)
	{
		numIter++;
		alpha = GoldenSectionSearch(fcode, n, x, antiGrad);
		x = sum_dv(n, x, mult_scalar_dv(n, alpha, antiGrad));					//x = x + alpha * antiGrad
		for (int i = 0; i < n; ++i)
		{
			if (x[i] < testfun::lb[fcode])
				//x[i] = testfun::lb[fcode];
				x[i] = distr_p(gen_p);
			if (x[i] > testfun::rb[fcode])
				//x[i] = testfun::rb[fcode];
				x[i] = distr_p(gen_p);
		}

		newGrad = mult_scalar_dv(n, -1.0, testfun::grad_f(n, x, fcode));		//newGrad = (-1) * Grad
		newGradSquare = inner_prod_dv(n, newGrad, newGrad);
		if (numIter % (5 * n) == 0)												//обновление beta
			beta = 0;
		else
			beta = newGradSquare / gradSquare;
		antiGrad = sum_dv(n, newGrad, mult_scalar_dv(n, beta, antiGrad));		//antiGrad = newGrad + beta * antiGrad
		prevVal = curVal;	
		curVal = testfun::f(n, x, fcode);
		Report r(testfun::call_f, curVal);
		fr.insert_into_report(r);
		gradSquare = newGradSquare;
	} 
	// do ... while (gradSquare > eps);
	return fr;
}

#endif
