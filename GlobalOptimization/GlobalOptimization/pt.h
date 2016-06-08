#ifndef PT_H_INCLUDED
#define PT_H_INCLUDED

#include <vector>
#include "testfun.h"

using namespace std;

class pt
{
public:
	int n;
	point x;
	pt() {}
	pt(const int & n_, const point & x_) : n(n_), x(x_) {}

	pt operator= (pt & p) {
		n = p.n; x = p.x;
		return p;
	}
	double operator* (pt & p) {
		double ans = 0;
		for (int i = 0; i < n; i++){
			ans += x[i] * p.x[i];
		}
		return ans;
	}
	pt operator* (double p) {
		pt ans = *this;
		for (int i = 0; i < n; i++)
			ans.x[i] *= p;
		return ans;
	}
	pt operator/ (double p) {
		pt ans = *this;
		for (int i = 0; i < n; i++)
			ans.x[i] /= p;
		return ans;
	}
	pt operator+ (pt & p){
		pt ans = *this;
		for (int i = 0; i < n; i++)
			ans.x[i] += p.x[i];
		return ans;
	}
	pt operator- (pt & p){
		pt ans = *this;
		for (int i = 0; i < n; i++)
			ans.x[i] -= p.x[i];
		return ans;
	}
	double mod() {
		return sqrt((*this)*(*this));
	}
	pt ort() {
		return (*this) / (*this).mod();
	}
	double f(const int & fcode, const int & n) {
		return testfun::f(n, x, fcode);
	}
	pt grad_f(const int & fcode, const int & n) {
		return pt(n, testfun::grad_f(n, x, fcode));
	}
	const point back() {
		return x;
	}
};
#endif
