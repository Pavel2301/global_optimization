#ifndef DIFEVO_H_INCLUDED
#define DIFEVO_H_INCLUDED

#include "person.h"
#include "population.h"

using namespace std;

static const unsigned difevo_seed = std::chrono::system_clock::now().time_since_epoch().count();
static mt19937 gen_difevo(difevo_seed);
static uniform_real_distribution<double> distr_alpha(0.5, 0.9);
static uniform_real_distribution<double> distr_alpha1(-1, 1);
static uniform_real_distribution<double> distr_alpha2(0, 1);
static uniform_real_distribution<double> distr_alpha3(0, 0.5);
static uniform_real_distribution<double> distr_alpha4(0.5, 0.9);

Person mutant_vector(const int &fcode, const Person &v1, const Person &v2, const Person &v3, const double &difevo_f)
{
	int n = v1.show_gene().size();
	point x(n, 0);
	for (int i = 0; i<n; ++i)
	{
		double alpha;
		alpha = (difevo_f == 0) ? distr_alpha(gen_difevo) : difevo_f;
		x[i] = v1.show_gene()[i] + alpha*(v2.show_gene()[i] - v3.show_gene()[i]);
	}
	Person p(fcode, x);
	return p;
}

Person mutant_vector(const int &fcode, const Person &v1, const Person &v2, const Person &v3, const double &difevo_f, const int &alpha_gen)
{
	uniform_real_distribution<double> distr_alpha(0.5, 0.9);
	int n = v1.show_gene().size();
	point x(n, 0);
	for (int i = 0; i<n; ++i)
	{
		double alpha;
		if (difevo_f == 0)
			switch (alpha_gen)
		{
			case 1: alpha = distr_alpha1(gen_difevo); break;
			case 2: alpha = distr_alpha2(gen_difevo); break;
			case 3: alpha = distr_alpha3(gen_difevo); break;
			case 4: alpha = distr_alpha4(gen_difevo); break;
		}
		else
			alpha = difevo_f;
		x[i] = v1.show_gene()[i] + alpha*(v2.show_gene()[i] - v3.show_gene()[i]);
	}
	Person p(fcode, x);
	return p;
}

FullReport difevo(const int &fcode, const int &n, const vector<point> &init_pop, const int &max_call_f, const double &mut_prob)
{
	FullReport fr;
	int psize = init_pop.size();
	vector<Person> w;
	for (int i = 0; i<psize; i++)
	{
		Person x(fcode, init_pop[i]);
		w.push_back(x);
	}
	Population p(w);
	uniform_int_distribution<int> difevo_distr(0, psize - 1);
	testfun::call_f = 0;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		vector<Person> next_pop;
		for (int k = 0; k < psize; k++)
		{
			int a = k, b = k, c = k;
			while ((a == k) || (b == k) || (c == k))
			{
				a = difevo_distr(gen_ga);
				b = difevo_distr(gen_ga);
				c = difevo_distr(gen_ga);
			}
			Person v = mutant_vector(fcode, p.show_persons()[a], p.show_persons()[b], p.show_persons()[c], 0);
			Person x = p.show_persons()[k];

			point q = pure_random(n);
			for (int i = 0; i<n; i++)
			{
				double coin = distr_coin(gen_ga);
				if (coin<mut_prob)
					q[i] = x.show_gene()[i];
				else
					q[i] = v.show_gene()[i];
			}
			Person res(fcode, q);

			if (person_cmp(res, x))
				next_pop.push_back(res);
			else
				next_pop.push_back(x);
		}
		Population h(next_pop);
		p = h;
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

FullReport difevo(const int &fcode, const int &n, const vector<point> &init_pop, const int &max_call_f, const double &mut_prob, const int &alpha_gen)
{
	FullReport fr;
	int psize = init_pop.size();
	vector<Person> w;
	for (int i = 0; i<psize; i++)
	{
		Person x(fcode, init_pop[i]);
		w.push_back(x);
	}
	Population p(w);
	uniform_int_distribution<int> difevo_distr(0, psize - 1);
	testfun::call_f = 0;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		vector<Person> next_pop;
		for (int k = 0; k < psize; k++)
		{
			int a = k, b = k, c = k;
			while ((a == k) || (b == k) || (c == k))
			{
				a = difevo_distr(gen_ga);
				b = difevo_distr(gen_ga);
				c = difevo_distr(gen_ga);
			}
			Person v = mutant_vector(fcode, p.show_persons()[a], p.show_persons()[b], p.show_persons()[c], 0, alpha_gen);
			Person x = p.show_persons()[k];

			point q = pure_random(n);
			for (int i = 0; i<n; i++)
			{
				double coin = distr_coin(gen_ga);
				if (coin<mut_prob)
					q[i] = x.show_gene()[i];
				else
					q[i] = v.show_gene()[i];
			}
			Person res(fcode, q);

			if (person_cmp(res, x))
				next_pop.push_back(res);
			else
				next_pop.push_back(x);
		}
		Population h(next_pop);
		p = h;
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

FullReport difevo(const int &fcode, const int &n, const int &psize, const int &max_call_f, const double &mut_prob)
{
	FullReport fr;
	uniform_int_distribution<int> difevo_distr(0, psize - 1);
	Population p(fcode, psize, n);
	testfun::call_f = 0;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		vector<Person> next_pop;
		for (int k = 0; k < psize; k++)
		{
			int a = k, b = k, c = k;
			while ((a == k) || (b == k) || (c == k))
			{
				a = difevo_distr(gen_ga);
				b = difevo_distr(gen_ga);
				c = difevo_distr(gen_ga);
			}
			Person v = mutant_vector(fcode, p.show_persons()[a], p.show_persons()[b], p.show_persons()[c], 0);
			Person x = p.show_persons()[k];

			point q = pure_random(n);
			for (int i = 0; i<n; i++)
			{
				double coin = distr_coin(gen_ga);
				if (coin<mut_prob)
					q[i] = x.show_gene()[i];
				else
					q[i] = v.show_gene()[i];
			}
			Person res(fcode, q);

			if (person_cmp(res, x))
				next_pop.push_back(res);
			else
				next_pop.push_back(x);
		}
		Population h(next_pop);
		p = h;
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

#endif