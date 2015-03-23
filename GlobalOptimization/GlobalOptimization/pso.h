#ifndef PSO_H_INCLUDED
#define PSO_H_INCLUDED

#include "particle.h"
#include "swarm.h"

using namespace std;

FullReport pso(const int &fcode, const int &n, const int &swarm_size, const int &max_call_f, const double &a1, const double &a2)
{
	FullReport fr;
	Swarm p_s(n, swarm_size, fcode);             //генерация роя частиц случайным образом
	testfun::call_f = 0;
	Report r(testfun::call_f, p_s.calc_averageValue());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		p_s.swarmMove(a1, a2, fcode);
		Report r(testfun::call_f, p_s.calc_averageValue());
		fr.insert_into_report(r);
	}
	return fr;
}

FullReport pso(const int &fcode, const int &n, vector<point> &init_swarm, const int &max_call_f, const double &a1, const double &a2)
{
	FullReport fr;
	int swarm_size = init_swarm.size();
	vector<Particle> w;
	for (int i = 0; i < swarm_size; i++)
	{
		Particle x(fcode, init_swarm[i]);
		w.push_back(x);
	}
	Swarm p_s(w);
	testfun::call_f = 0;
	Report r(testfun::call_f, p_s.calc_averageValue());
	fr.insert_into_report(r);
	while (testfun::call_f < max_call_f)
	{
		p_s.swarmMove(a1, a2, fcode);
		Report r(testfun::call_f, p_s.calc_averageValue());
		fr.insert_into_report(r);
	}
	return fr;
}


#endif