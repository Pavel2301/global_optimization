#ifndef POPULATION2_H_INCLUDED
#define POPULATION2_H_INCLUDED

#include "person.h"

static const unsigned seed_ais = std::chrono::system_clock::now().time_since_epoch().count();
static mt19937 gen_ais(seed_ais);
uniform_real_distribution<double> distr_coin2(0.0, 1.0);

class Population2
{
private:
	vector<Person> persons2;
	int p_size;
public:
	Population2(const vector<Person> &pop) :persons2(pop)
	{
		p_size = persons2.size();
	}
	Population2(const Population2 &pop)
	{
		persons2 = pop.persons2;
		p_size = pop.p_size;
	}
	Population2(const int &fcode, const int &pop_size, const int &n) :p_size(pop_size)   //генераци¤ попул¤ции случайным образом
	{
		vector<Person> p;
		for (int i = 0; i<pop_size; ++i)
		{
			Person x(fcode, n);
			p.push_back(x);
		}
		persons2 = p;
	}

	int show_p_size()const
	{
		return p_size;
	}
	vector<Person> show_persons2()const
	{
		return persons2;
	}
	vector<point> show_points()const
	{
		vector<point> v;
		for (int i = 0; i<p_size; ++i)
			v.push_back(persons2[i].show_gene());
		return v;
	}
	double calc_avf()const                         //средн¤¤ приспособленность попул¤ции
	{
		double avf = 0;
		for (int i = 0; i<p_size; ++i)
			avf = avf + persons2[i].show_f();
		return avf / p_size;
	}

	double min_f()const                                 //минимальна¤ приспособленность в попул¤ции
	{
		double min = persons2[0].show_f();
		for (int i = 1; i < p_size; i++)
		{
			if (persons2[i].show_f() < min)
				min = persons2[i].show_f();
		}
		return min;
	}

	Population2 clones(const int &fcode, const int &psize, const int &s_kol, const int &s, const int &nc)//"клонирование особей, равномерный случай"
	{
		vector<Person> p;
		uniform_int_distribution<int> s_kol_distr(0, p_size - 1);
		
		for (int i = 0; i < psize ; ++i)
		{

			int k = 0;
			int s = s_kol_distr(gen_ais);
			double min_f = persons2[i].show_f();
			Person min_p = persons2[i];
			for (int j = i + 1; j < psize; ++j){
				s = s_kol_distr(gen_ais);
				if (persons2[i].show_f()>persons2[j].show_f()) {
					min_p = persons2[j];
					min_f = persons2[j].show_f(); swap(persons2[i], persons2[j]);
				}
			}
			p.push_back(min_p);
		}	
		for (int kk = 0; kk < s_kol; kk++){
			Person min_dd = p[kk];
			
		int mm = 0;
		while (mm != nc)
		{
			p.push_back(min_dd); mm++;
		}
		}
		//for (int i = 0; i < p.size(); i++) cout << p[i].show_f() << endl;
		Population2 h1(p);
		return h1;
	}
	void mutantgenesis(const int &fcode, const int &mutcode, const int &psize, const int &n, const int &mut_prob, const int &r1)
		//мутаци¤ дл¤ клонов AIS.
	{
		//cout << p_size << " ";
		//for(int i = 0; i < p_size; i++) cout << persons2[i].show_f() << " ";
		
		for (int i = psize; i < p_size; ++i)
		{
					double mutradius;
					mutation(fcode, mutcode, mutradius, persons2[i]);	
		}
		
	}
	void selector(const int &s_kol, const int &nc, const int &psize) //селекци¤
	{
		//vector<Person> p;
		//uniform_int_distribution<int> parents_distr(0, p_size - 1);
		for (int i = 0; i < s_kol; i++){
			for (int j = psize + i; j < psize + i + nc; j++){
				if (persons2[i].show_f() < persons2[j].show_f()) persons2[i] = persons2[j];
			}
		}
		
	}
	void worst(const int &fcode, const int &d, const int &psize, const int &n)//замена значений худших показателей
	{
		vector<Person> p;
		uniform_int_distribution<int> parents_distr(0, p_size - 1);
		
		popsort();
		for (int kk = psize - d; kk < psize; kk++)
		{
			Person x(fcode, n);
			persons2[kk] = x;
		}
	}

	void popsort()               //сортировка попул¤ции
	{
		sort(persons2.begin(), persons2.end(), person_cmp);
	}


	Population2 reduction(const int &psize)          //сокращение попул¤ции (выбираем наиболее приспособленных особей)
	{
		if (psize <= p_size)
		{
			popsort();
			vector<Person>p;
			for (int i = 0; i<psize; ++i)
				p.push_back(persons2[i]);
			Population2 h1(p);
			return h1;
		}
		else
		{
			return *this;
		}
	}

	friend const Population2 operator+(const Population2 &p1, const Population2 &p2)
	{
		vector<Person> merger;
		for (const auto& value : p1.persons2)
			merger.push_back(value);
		for (const auto& value : p2.persons2)
			merger.push_back(value);
		Population2 p(merger);
		return p;
	}

};

#endif

