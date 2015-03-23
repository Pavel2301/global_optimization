#ifndef PERSON_H_INCLUDED
#define PERSON_H_INCLUDED

#include <algorithm>

class Person
{
private:
	point gene;       //ген (особь)
	double fitness;   //приспособленность особи
public:
	Person(const int &fcode, const point &g) : gene(g)      //конструктор из объекта point
	{
		fitness = testfun::f(gene.size(), gene, fcode);
	}
	Person(const Person &p)   //конструктор копирования
	{
		gene = p.gene;
		fitness = p.fitness;
	}
	Person(const int &fcode, const int &n)  //конструктор случайным образом
	{
		point x = pure_random(n, testfun::lb[fcode], testfun::rb[fcode]);
		gene = x;
		fitness = testfun::f(gene.size(), gene, fcode);
	}

	/*Person(const int &fcode, const Particle &p)
	{
	gene = p._currentPosition;
	fitness = testfun::f(gene.size(), gene, fcode);
	}*/
	double recalc_f(const int &fcode)
	{
		fitness = testfun::f(show_n(), show_gene(), fcode);
		return fitness;
	}
	double show_f()const{ return fitness; }
	point show_gene()const{ return gene; }
	int show_n()const{ return gene.size(); }

};

bool person_cmp(const Person &lhs, const Person &rhs)
{
	return lhs.show_f()<rhs.show_f() ? true : false;
}


inline void crossover(const int &fcode, const int &k, const Person &p1, const Person &p2, Person &c1, Person &c2)
{
	/**k-число разрывов**/
	int n = p1.show_n();
	if (p2.show_n() == n)
	{
		vector<int> t(n - 1, 0);
		for (int i = 0; i < n - 1; ++i)
			t[i] = i;
		random_shuffle(t.begin(), t.end());
		vector<int> order(k, 0);
		for (int i = 0; i < k; ++i)
			order[i] = t[i];
		sort(order.begin(), order.end());
		bool flag = true;
		int s = 0;
		point a = p1.show_gene();
		point b = p2.show_gene();

		for (int i = 0; i < k; ++i)
		{
			for (int j = s; j < order[i]; ++j)
				if (flag)
				{
					a[j] = p2.show_gene()[j];
					b[j] = p1.show_gene()[j];
				}
			s = order[i];
			flag = !flag;
		}
		for (int j = s; j < n; ++j)
			if (flag)
			{
				a[j] = p2.show_gene()[j];
				b[j] = p1.show_gene()[j];
			}

		Person x1(fcode, a);
		Person x2(fcode, b);
		c1 = x1;
		c2 = x2;
	}
	else
	{
		c1 = p1;
		c2 = p2;
	}
}

inline void crossover(const int &fcode, const int &k, const Person &p1, const Person &p2, Person &c)    //скрещивание (номер функции, число разрывов, 1-я особь, 2-я особь, c) - многоточечный кроссинговер
{
	int n = p1.show_n();            //размерность 1-й особи
	if (p2.show_n() == n)           //если 2 особи имеют одинаковый размер
	{
		vector<int> t(n - 1, 0);
		for (int i = 0; i < n - 1; ++i)   //заполнение вектора порядковых номеров особей  
			t[i] = i;
		random_shuffle(t.begin(), t.end());   //случайная перетасовка вектора t 
		vector<int> order(k, 0);
		for (int i = 0; i<k; ++i)        //заносим в order элементы перетасованного вектора t
			order[i] = t[i];
		sort(order.begin(), order.end());    //сортируем order
		bool flag = true;
		int s = 0;
		point a = p1.show_gene();               //1-я особь     
		point b = p2.show_gene();               //2-я особь

		for (int i = 0; i<k; ++i)
		{
			for (int j = s; j<order[i]; ++j)
				if (flag)
				{
					a[j] = p2.show_gene()[j];
					b[j] = p1.show_gene()[j];
				}
			s = order[i];
			flag = !flag;
		}
		for (int j = s; j<n; ++j)
			if (flag)
			{
				a[j] = p2.show_gene()[j];
				b[j] = p1.show_gene()[j];
			}
		Person x1(fcode, a);
		Person x2(fcode, b);
		c = (person_cmp(x1, x2)) ? x1 : x2;
	}
	else
	{
		c = (person_cmp(p1, p2)) ? p1 : p2;
	}
}

void mutation(const int &fcode, const int &mutcode, double &mut_radius, Person &p)     //мутация (номер функции, тип мутации, радиус мутации, особь)
{
	if (mutcode == 1)
	{
		point q = p.show_gene();        //особь для мутации
		int k = q.size();               //размерность
		point q2 = stoh_inner_centered_ball_scaling(k, q, mut_radius);
		Person p2(fcode, q2);
		p = p2;
	}
	if (mutcode == 2)
	{
		point q = p.show_gene();        //особь для мутации
		int k = q.size();               //размерность
		point q2 = stoh_inner_centered_ball_scaling2(k, q, mut_radius);
		Person p2(fcode, q2);
		p = p2;
	}
	if (mutcode == 3)
	{
		static const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		static mt19937 gen(seed);
		point q = p.show_gene();
		int n = q.size();
		uniform_int_distribution<int> distr_c(0, n - 1);
		uniform_real_distribution<double> distr(testfun::lb[fcode], testfun::rb[fcode]);
		//for (int i = 0; i < 10; i++)
		//	{
		int c = distr_c(gen);
		q[c] = distr(gen);
		//}
		Person p2(fcode, q);
		p = p2;
	}

}

void mutation2(const int &fcode, const int &mutcode, Person &p)
{
	static const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	static mt19937 gen(seed);
	point q = p.show_gene();
	int n = q.size();
	uniform_int_distribution<int> distr_c(0, n - 1);
	uniform_real_distribution<double> distr(testfun::lb[fcode], testfun::rb[fcode]);
	int c = distr_c(gen);
	q[c] = distr(gen);
	Person p2(fcode, q);
	p = p2;
}

void mutation2b(const int &fcode, const int &mutcode, Person &p, const int percent)
{
	static const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	static mt19937 gen(seed);
	point q = p.show_gene();
	int n = q.size();
	double num_c = n * percent / 100;
	if (num_c > 0 && num_c < 0.5)
		num_c = 1;
	else
		num_c = floor(num_c + 0.5);
	uniform_real_distribution<double> distr(testfun::lb[fcode], testfun::rb[fcode]);
	for (int i = 0; i < num_c; i++)
	{
		uniform_int_distribution<int> distr_c(0, n - 1);
		int c = distr_c(gen);
		q[c] = distr(gen);
	}
	Person p2(fcode, q);
	p = p2;
}


inline double person_dist(const Person &lhs, const Person &rhs)
{
	point lp = lhs.show_gene();
	point rp = rhs.show_gene();
	int ln = lp.size();
	int rn = rp.size();
	return (ln == rn) ? dist_points(ln, lhs.show_gene(), rhs.show_gene()) : 0;
}

#endif

