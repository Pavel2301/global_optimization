#ifndef POPULATION_H_INCLUDED
#define POPULATION_H_INCLUDED

#include "person.h"

static const unsigned seed_ga = std::chrono::system_clock::now().time_since_epoch().count();
static mt19937 gen_ga(seed_ga);
uniform_real_distribution<double> distr_coin(0.0, 1.0);

class Population
{
private:
	vector<Person> persons;
	int p_size;
public:
	Population(const vector<Person> &pop) :persons(pop)
	{
		p_size = persons.size();
	}
	Population(const Population &pop)
	{
		persons = pop.persons;
		p_size = pop.p_size;
	}
	Population(const int &fcode, const int &pop_size, const int &n) :p_size(pop_size)   //генераци¤ попул¤ции случайным образом
	{
		vector<Person> p;
		for (int i = 0; i<pop_size; ++i)
		{
			Person x(fcode, n);
			p.push_back(x);
		}
		persons = p;
	}

	/*Population(const Swarm &s, const int &fcode, const int &n)
	{
	p_size = s.show_swarm_size();
	for (int i = 0; i < p_size; i++)
	{
	Person p(fcode, s.show_ps()[i]);
	persons.push_back(p);
	}
	}*/

	int show_p_size()const
	{
		return p_size;
	}
	vector<Person> show_persons()const
	{
		return persons;
	}
	vector<point> show_points()const
	{
		vector<point> v;
		for (int i = 0; i<p_size; ++i)
			v.push_back(persons[i].show_gene());
		return v;
	}
	double calc_avf()const                         //средн¤¤ приспособленность попул¤ции
	{
		double avf = 0;
		for (int i = 0; i<p_size; ++i)
			avf = avf + persons[i].show_f();
		return avf / p_size;
	}

	double min_f()const                                 //минимальна¤ приспособленность в попул¤ции
	{
		double min = persons[0].show_f();
		for (int i = 1; i < p_size; i++)
		{
			if (persons[i].show_f() < min)
				min = persons[i].show_f();
		}
		return min;
	}

	Population tournament_selection(const int &fcode, const int &parents_size, const int &tournament_size)const   //селекци¤ (номер функции, число родителей, численность турнира) - турнирный отбор
	{
		vector<Person> p;
		uniform_int_distribution<int> parents_distr(0, p_size - 1);
		for (int i = 0; i<parents_size; ++i)
		{
			int s = parents_distr(gen_ga);
			double min_f = persons[s].show_f();         //минимальна¤ приспособленность = приспособленность случайной особи
			Person min_p = persons[s];                  //особь с минимальной приспособленностью
			for (int j = 1; j<tournament_size; ++j)     //нахождение особи с минимальной приспособленностью   
			{
				s = parents_distr(gen_ga);
				if (min_f>persons[s].show_f())
				{
					min_p = persons[s];
					min_f = persons[s].show_f();
				}
			}
			p.push_back(min_p);  //добавление в вектор родител¤
		}

		Population h(p);   //составление попул¤ции дл¤ скрещивани¤
		return h;
	}

	Population rank_selection(const int &fcode, const int &parents_size)            //рангова¤ селекци¤
	{
		vector<Person> p;
		popsort();
		for (int i = 0; i < parents_size; i++)
			p.push_back(persons[i]);
		Population h(p);
		return h;
	}

	Population popcrossover(const int &fcode, const int &children_size, const int &num_k)     //скрещивание (номер функции, количество детей)
	{
		vector<Person> p;         //попул¤ци¤
		uniform_int_distribution<int> children_distr(0, p_size - 1);
		for (int i = 0; i<children_size; ++i)
		{
			int a = children_distr(gen_ga);
			int b = children_distr(gen_ga);
			Person pa = persons[a];             // 1-¤ случайна¤ особь (из родительской попул¤ции)
			Person pb = persons[b];             // 2-¤ случайна¤ особь (из родительской попул¤ции)
			Person c = pa;
			crossover(fcode, num_k, pa, pb, c);     //скрещивание
			p.push_back(c);                     //добавление полученной особи в вектор p
		}
		Population q(p);
		return q;
	}

	void pmutation(const int &fcode, const int &mutcode, const double &mutprob)    //мутаци¤ (номер функции, тип мутации, веро¤тность мутации) 
	{
		for (int i = 0; i<p_size; ++i)
		{
			double coin = distr_coin(gen_ga);         //генерируем число от 0 до 1
			if (coin<mutprob)
			{
				double mutradius;
				mutation(fcode, mutcode, mutradius, persons[i]);
			}
		}
	}

	double pmutation_rand(const int &fcode, const int &mutcode)
	{
		for (int i = 0; i<p_size; ++i)
		{
			double mutprob = distr_coin(gen_ga);
			double coin = distr_coin(gen_ga);         //генерируем число от 0 до 1
			if (coin<mutprob)
			{
				double mutradius;
				mutation(fcode, mutcode, mutradius, persons[i]);
			}
		}
	}

	void pmutation2(const int &fcode, const int &mutcode, const double &mutprob)    //мутаци¤ (номер функции, тип мутации, веро¤тность мутации) 
	{
		for (int i = 0; i<p_size; ++i)
		{
			double coin = distr_coin(gen_ga);         //генерируем число от 0 до 1
			if (coin<mutprob)
				mutation2(fcode, mutcode, persons[i]);
			//mutation2b(fcode, mutcode, persons[i], 2);
		}
	}

	void popsort()               //сортировка попул¤ции
	{
		sort(persons.begin(), persons.end(), person_cmp);
	}


	Population reduction(const int &psize)          //сокращение попул¤ции (выбираем наиболее приспособленных особей)
	{
		if (psize <= p_size)
		{
			popsort();
			vector<Person>p;
			for (int i = 0; i<psize; ++i)
				p.push_back(persons[i]);
			Population h(p);
			return h;
		}
		else
		{
			return *this;
		}
	}

	friend const Population operator+(const Population &p1, const Population &p2)
	{
		vector<Person> merger;
		for (const auto& value : p1.persons)
			merger.push_back(value);
		for (const auto& value : p2.persons)
			merger.push_back(value);
		Population p(merger);
		return p;
	}

};

#endif

