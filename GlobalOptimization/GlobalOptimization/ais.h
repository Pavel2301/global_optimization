#ifndef AIS_H_INCLUDED
#define AIS_H_INCLUDED

#include "population2.h"
#include "person.h"
#include "report.h"
#include "testfun.h"
using namespace std;

FullReport ais1(const int &fcode, const int &n, const int &psize, const int &max_call_f, const int &s_kol, const int &s,const int &d,
	const int &k, const int &beta, const int &nc, const int &r1, const double &mut_prob, const int mutcode)
/*параметры - номер функции, размерность, размер попул¤ции(количество имунных клеток в попул¤ции), 
макс. число вызовов ф-ции,
число клеток, выбираемых из попул¤ции в ходе селекции, 
число клеток с наихудшим состо¤нием приспособленности, максимальное число итераций.число вызовов функции,
параметр клонировани¤ в зависмости от выбора операции(beta|nc), параметр мутации, веро¤тность мутации
*/
{
	FullReport fr;
	Population2 p(fcode, psize, n); //генераци¤ начальной попул¤ции
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		Population2 first = p.clones(fcode, psize, s_kol, s, nc);//клонирование
		first.mutantgenesis(fcode, mutcode, psize, n, mut_prob, r1);//мутаци¤
		first.selector(s_kol, nc, psize);//селекци¤, формирование новой мутации
		p = first.reduction(psize);//уменьшение попул¤ции до исходного объема
		p.worst(fcode, d, psize,n);//обновление попул¤ции
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

FullReport ais1(const int &fcode, const int &n, const vector<point> &init_pop, const int &max_call_f, const int &s_kol, const int &s1, const int &d,
	const int &k, const int &beta, const int &nc, const int &r1, const double &mut_prob, const int mutcode)
	/*параметры - номер функции, размерность, размер попул¤ции(количество имунных клеток в попул¤ции),
	макс. число вызовов ф-ции,
	число клеток, выбираемых из попул¤ции в ходе селекции,
	число клеток с наихудшим состо¤нием приспособленности, максимальное число итераций.число вызовов функции,
	параметр клонировани¤ в зависмости от выбора операции(beta|nc), параметр мутации, веро¤тность мутации
	*/
{
	FullReport fr;
	int psize = init_pop.size();
	vector<Person> w;
	for (int i = 0; i<psize; i++)
	{
		Person x(fcode, init_pop[i]);
		w.push_back(x);
	}
	Population2 p(w);
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	double max_dist = 0;
	while (testfun::call_f < max_call_f)
	{
		Population2 first = p.clones(fcode, psize, s_kol, s1, nc);//клонирование
		first.mutantgenesis(fcode, mutcode, psize, n, mut_prob, r1);//мутаци¤
		first.selector(s_kol, nc, psize);//селекци¤, формирование новой мутации
		p = first.reduction(psize); //уменьшение попул¤ции до исходного объема
		p.worst(fcode, d, psize, n);//обновление попул¤ции
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

#endif
