#ifndef GA_H_INCLUDED
#define GA_H_INCLUDED

#include "population.h"
#include "person.h"

using namespace std;

FullReport ga(const int &fcode, const int &n, const int &psize,
	const int &max_call_f, const int &tournament_size, const int &mut_type, const double &mut_prob,
	const int &parents_size, const int &children_size, const int &num_k)
{
	/*параметры - номер функции, размерность, размер попул¤ции, макс. число вызовов ф-ции, размер турнира, тип мутации, веро¤тность мутации, число родителей, число детей, elitesize, число разрывов, стартовое число вызовов ф-ции*/
	FullReport fr;
	Population p(fcode, psize, n);              //генераци¤ начальной попул¤ции
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	while (testfun::call_f<max_call_f)
	{
		Population parents = p.tournament_selection(fcode, parents_size, tournament_size);              //селекци¤
		Population children = parents.popcrossover(fcode, children_size, num_k);                //скрещивание
		children.pmutation(fcode, mut_type, mprob);												//мутаци¤
		Population temp = children + p;
		p = temp.reduction(psize);                                                            //сокращение попул¤ции
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}

FullReport ga(const int &fcode, const int &n, const vector<point> &init_pop,
	const int &max_call_f, const int &tournament_size, const int &mut_type, const double &mut_prob,
	const int &parents_size, const int &children_size, const int &num_k)
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
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_avf());
	fr.insert_into_report(r);
	double max_dist = 0;
	while (testfun::call_f<max_call_f)
	{
		Population parents = p.tournament_selection(fcode, parents_size, tournament_size);
		Population children = parents.popcrossover(fcode, children_size, num_k);
		children.pmutation(fcode, mut_type, mprob);
		Population temp = children + p;
		p = temp.reduction(psize);
		Report r(testfun::call_f, p.calc_avf());
		fr.insert_into_report(r);
	}
	p.popsort();
	return fr;
}



#endif