#include "report.h"
//#include "testfun.h"
#include "pt.h"
#include <iostream>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>

using namespace std;

typedef double db;

const db t = (3.0 - sqrt(5.0))*0.5, eps = 1e-10;
db b = 1000;

pt antigrad(const int &fcode, const int &n, pt &x)
{
	return x.grad_f(fcode, n)*(-1);
}

db find_border(const int &n, pt &x, pt &d)
{
	db k = 1e5;
	for (int i = 0; i < n; i++)
	{
		if (d.x[i]>0) k = min(k, (b - x.x[i]) / d.x[i]);
		if (d.x[i]<0) k = min(k, (-b - x.x[i]) / d.x[i]);
	}
	return max(0.0, k);
}

db R;
pt find_min(const int &fcode, const int &n, pt &x, pt &d)
{
	d = d.ort();
	db l = 0, k1 = l + (R - l)*t, k2 = R - (R - l)*t;
	db f1 = (x + d*k1).f(fcode, n), f2 = (x + d*k2).f(fcode, n);

	while (abs(l - R) > eps)
	{
		if (f1 < f2)
		{
			R = k2; k2 = k1; f2 = f1;
			k1 = l + (R - l)*t;
			f1 = (x + d*k1).f(fcode, n);
		}
		else
		{
			l = k1; k1 = k2; f1 = f2;
			k2 = R - (R - l)*t;
			f2 = (x + d*k2).f(fcode, n);
		}
	}
	return x + (d*l);
}

pt f_binary(const int &fcode, const int &n, pt &x, pt &d)
{
	pt newx = find_min(fcode, n, x, d);
	while (newx.f(fcode, n) > x.f(fcode, n))
	{
		R /= 2;
		newx = find_min(fcode, n, x, d);
	}
	return newx;
}

pt f_rand(const int &fcode, const int &n, pt &x, pt &d)
{
	srand(time(NULL));
	int cnt = 20;
	db newR = R, mn = x.f(fcode, n);
	pt newx = x;
	while (cnt--)
	{
		x = find_min(fcode, n, x, d);
		db value = x.f(fcode, n);
		if (mn > value) { mn = value; newx = x; }
		db num = rand() % 1000 + 1;
		num /= 1000;
		R = R*num;
	}
	return newx;
}

pt f_rand_cnt(const int &fcode, const int &n, pt &x, pt &d)
{
	srand(time(NULL));
	int cnt = 20;
	db newR = R, mn = x.f(fcode, n);
	pt newx = x;
	while (cnt--)
	{
		x = find_min(fcode, n, x, d);
		db value = x.f(fcode, n);
		if (mn > value) { mn = value; newx = x; }
		db num = rand() % 1000 + 1;
		num /= 1000;
		R = newR*num;
	}
	return newx;
}

FullReport cgm(const int &fcode, const int &n, const point &p, const int &max_call_f, const int &limit, const int &cgm_type)
{
	FullReport fr;
	pt x(n, p), r, d;
	b = limit;
	d = r = antigrad(fcode, n, x);
	int kol = 0;
	Report rep(kol, testfun::f(n, p, fcode));
	fr.insert_into_report(rep);
	do
	{
		R = find_border(n, x, d.ort());
		kol++;
		//x = find_min(fcode, n, x, d);
		//x = f_rand_cnt(fcode, n, x, d);
		//x = f_rand(fcode, n, x, d);
		x = f_binary(fcode, n, x, d);
		Report rep(kol, testfun::f(n, x.back(), fcode));
		fr.insert_into_report(rep);
		pt new_r = antigrad(fcode, n, x);
		db beta;
		switch (cgm_type)
		{
		case 1: {beta = max(0.0, (new_r*(new_r - r)) / (r*r)); break; }								    //Polak-Ribiere
		case 2: {beta = max(0.0, (new_r*(new_r - r)) / ((new_r - r)*d));	break; }				//Hestenes-Stiefel
		case 3: {beta = max(0.0, -(new_r*new_r) / (d*(new_r - r))); break; }						  //Dai-Yuan
		}
		d = new_r + (d*beta);
		r = new_r;
	} while (kol<max_call_f);
	return fr;
}
