#ifndef GAPSO_H_INCLUDED
#define GAPSO_H_INCLUDED

#include "popswarm.h"
#include "geneparticle.h"

using namespace std;

FullReport ga_pso(const int &fcode, const int &n, const int &psize,
	const int &max_call_f, const double &a1, const double &a2, const int &tournament_size, const int &mut_type, const double &mut_prob,
	const int &parents_size, const int &children_size, const int &num_k)
{
	/*параметры - номер функции, размерность, размер попул¤ции, макс. число вызовов ф-ции, размер турнира, тип мутации, веро¤тность мутации, число родителей, число детей, elitesize, число разрывов, стартовое число вызовов ф-ции*/
	FullReport fr;
	PopSwarm p(fcode, psize, n);              //генераци¤ начальной попул¤ции
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_av_fitness());
	fr.insert_into_report(r);
	while (testfun::call_f<max_call_f)
	{
		PopSwarm parents = p.popSelection(fcode, parents_size, tournament_size);              //селекци¤
		PopSwarm children = parents.popCrossover(fcode, children_size, num_k);                //скрещивание
		//mrad = children.pmutation(fcode, mut_type, mprob);    
		children.popMutation(fcode, mut_type, mprob);
		children.swarmMove(a1, a2, fcode);
		//children.pmutation2(fcode, mut_type, mprob);
		PopSwarm temp = children + p;
		p = temp.popReduction(psize);                                                            //сокращение попул¤ции
		//cout << testfun::call_f << ";" << point_to_comma(to_string(p.calc_avf())) << ";" << mprob << ";" << max_dist << ";" << mrad << endl;
		//cout << point_to_comma(to_string(p.calc_avf())) << ";";
		Report r(testfun::call_f, p.calc_av_fitness());
		fr.insert_into_report(r);
	}
	//cout << endl;
	p.popSort();
	return fr;
}

FullReport ga_pso(const int &fcode, const int &n, const vector<point> &init_pop,
	const int &max_call_f, const double &a1, const double &a2, const int &tournament_size, const int &mut_type, const double &mut_prob,
	const int &parents_size, const int &children_size, const int &num_k)
{
	/*параметры - номер функции, размерность, размер попул¤ции, макс. число вызовов ф-ции, размер турнира, тип мутации, веро¤тность мутации, число родителей, число детей, elitesize, число разрывов, стартовое число вызовов ф-ции*/
	FullReport fr;
	int psize = init_pop.size();
	vector<GeneParticle> w;
	for (int i = 0; i<psize; i++)
	{
		GeneParticle x(fcode, init_pop[i]);
		w.push_back(x);
	}
	PopSwarm p(w);              //генераци¤ начальной попул¤ции
	testfun::call_f = 0;
	double mprob = mut_prob;
	Report r(testfun::call_f, p.calc_av_fitness());
	fr.insert_into_report(r);
	while (testfun::call_f<max_call_f)
	{
		PopSwarm parents = p.popSelection(fcode, parents_size, tournament_size);              //селекци¤
		PopSwarm children = parents.popCrossover(fcode, children_size, num_k);                //скрещивание
		//mrad = children.pmutation(fcode, mut_type, mprob);    
		children.popMutation(fcode, mut_type, mprob);
		children.swarmMove(a1, a2, fcode);
		//children.pmutation2(fcode, mut_type, mprob);
		PopSwarm temp = children + p;
		p = temp.popReduction(psize);                                                            //сокращение попул¤ции
		//cout << testfun::call_f << ";" << point_to_comma(to_string(p.calc_avf())) << ";" << mprob << ";" << max_dist << ";" << mrad << endl;
		//cout << point_to_comma(to_string(p.calc_avf())) << ";";
		Report r(testfun::call_f, p.calc_av_fitness());
		fr.insert_into_report(r);
	}
	//cout << endl;
	p.popSort();
	return fr;
}

#endif