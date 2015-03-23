#ifndef GENEPARTICLE_H_INCLUDED
#define GENEPARTICLE_H_INCLUDED

#include "point.h"

class GeneParticle
{
public:
	point _currentPoint;                     //текущая точка
	double _currentFitness;					//значение ЦФ в текущей точке
	point _localBestPoint;					//лучшая точка			     
	double _localBestFitness;				//значение ЦФ в лучшей точке
	vector<double> _speed;					//скорость

	GeneParticle() {}
	GeneParticle(const int &fcode, const point &p) : _currentPoint(p)                    //конструктор из объекта point
	{
		_currentFitness = testfun::f(_currentPoint.size(), _currentPoint, fcode);
		_localBestPoint = _currentPoint;
		_localBestFitness = _currentFitness;
		_speed.resize(_currentPoint.size());
		_speed.assign(_speed.size(), 0);
	}
	GeneParticle(const int &fcode, const point &p, const GeneParticle &old) : _currentPoint(p)                    //конструктор из объекта point
	{
		_currentFitness = testfun::f(_currentPoint.size(), _currentPoint, fcode);
		if (_currentFitness < old._currentFitness)
		{
			_localBestPoint = _currentPoint;
			_localBestFitness = _currentFitness;
		}
		else
		{
			_localBestPoint = old._currentPoint;
			_localBestFitness = old._currentFitness;
		}
		_speed = old._speed;	                                                               //?
	}
	GeneParticle(const int &n, const int &fcode)     //конструктор случайным образом
	{
		_currentPoint = pure_random(n, testfun::lb[fcode], testfun::rb[fcode]);
		_currentFitness = testfun::f(n, _currentPoint, fcode);
		_localBestPoint = _currentPoint;
		_localBestFitness = _currentFitness;
		_speed.resize(_currentPoint.size());
		_speed.assign(_speed.size(), 0);
	}
	void checkPoint(GeneParticle &x, const int &fcode)		//штраф отдельной частицы
	{
		for (int j = 0; j < x._currentPoint.size(); j++)
		{
			if (x._currentPoint[j] < testfun::lb[fcode])
				x._currentPoint[j] = testfun::lb[fcode];
			if (x._currentPoint[j] > testfun::rb[fcode])
				x._currentPoint[j] = testfun::rb[fcode];
		}
	}
};

bool point_cmp(const GeneParticle &p1, const GeneParticle &p2)
{
	return p1._currentFitness<p2._currentFitness ? true : false;
}

inline void crossoverPoint(const int &fcode, const int &k, const GeneParticle &p1, const GeneParticle &p2, GeneParticle &c1, GeneParticle &c2)
{
	/**k-число разрывов**/
	int n = p1._currentPoint.size();
	if (p2._currentPoint.size() == n)
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
		point a = p1._currentPoint;
		point b = p2._currentPoint;
		for (int i = 0; i < k; ++i)
		{
			for (int j = s; j < order[i]; ++j)
				if (flag)
				{
					a[j] = p2._currentPoint[j];
					b[j] = p1._currentPoint[j];
				}
			s = order[i];
			flag = !flag;
		}
		for (int j = s; j < n; ++j)
			if (flag)
			{
				a[j] = p2._currentPoint[j];
				b[j] = p1._currentPoint[j];
			}

		GeneParticle x1(fcode, a, p1);
		GeneParticle x2(fcode, b, p2);
		c1 = x1;
		c2 = x2;
	}
	else
	{
		c1 = p1;
		c2 = p2;
	}
}

inline void crossoverPoint(const int &fcode, const int &k, const GeneParticle &p1, const GeneParticle &p2, GeneParticle &c)    //скрещивание (номер функции, число разрывов, 1-я особь, 2-я особь, c) - многоточечный кроссинговер
{
	int n = p1._currentPoint.size();            //размерность 1-й особи
	if (p2._currentPoint.size() == n)           //если 2 особи имеют одинаковый размер
	{
		//std::srand(unsigned(std::time(0)));
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
		point a = p1._currentPoint;               //1-я особь     
		point b = p2._currentPoint;               //2-я особь

		for (int i = 0; i<k; ++i)
		{
			for (int j = s; j<order[i]; ++j)
				if (flag)
				{
					a[j] = p2._currentPoint[j];
					b[j] = p1._currentPoint[j];
				}
			s = order[i];
			flag = !flag;
		}
		for (int j = s; j<n; ++j)
			if (flag)
			{
				a[j] = p2._currentPoint[j];
				b[j] = p1._currentPoint[j];
			}
		GeneParticle x1(fcode, a, p1);
		GeneParticle x2(fcode, b, p2);
		c = (point_cmp(x1, x2)) ? x1 : x2;
	}
	else
	{
		c = (point_cmp(p1, p2)) ? p1 : p2;
	}
}

void mutation(const int &fcode, const int &mutcode, double &mut_radius, GeneParticle &p)     //мутация (номер функции, тип мутации, радиус мутации, особь)
{
	if (mutcode == 1)
	{
		GeneParticle q = p;        //особь для мутации
		int k = q._currentPoint.size();               //размерность
		point q2 = stoh_inner_centered_ball_scaling(k, q._currentPoint, mut_radius);
		GeneParticle p2(fcode, q2, q);
		p = p2;
	}
	if (mutcode == 2)
	{
		GeneParticle q = p;        //особь для мутации
		int k = q._currentPoint.size();                //размерность
		point q2 = stoh_inner_centered_ball_scaling2(k, q._currentPoint, mut_radius);
		GeneParticle p2(fcode, q2, q);
		p = p2;
	}
	if (mutcode == 3)
	{
		static const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		static mt19937 gen(seed);
		GeneParticle q = p;
		int n = q._currentPoint.size();
		uniform_int_distribution<int> distr_c(0, n - 1);
		uniform_real_distribution<double> distr(testfun::lb[fcode], testfun::rb[fcode]);
		//for (int i = 0; i < 10; i++)
		//	{
		int c = distr_c(gen);
		q._currentPoint[c] = distr(gen);
		//}
		GeneParticle p2(fcode, q._currentPoint, p);
		p = p2;
	}

}



#endif