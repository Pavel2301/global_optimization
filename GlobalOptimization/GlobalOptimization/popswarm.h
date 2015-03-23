#ifndef POPSWARM_H_INCLUDED
#define POPSWARM_H_INCLUDED

#include "geneparticle.h"

static const unsigned seed_gapso = std::chrono::system_clock::now().time_since_epoch().count();
static mt19937 gen_gapso(seed_gapso);
uniform_real_distribution<double> distr01(0.0, 1.0);

class PopSwarm
{
private:
	int _pop_size;                                      //размер попул¤ции
	vector<GeneParticle> _pswarm;                       //попул¤ци¤
	point _globalBestPoint;				                //лучша¤ точка из попул¤ции
	double _globalBestFitness;				            //значение ÷‘ в лучшей точке
public:
	PopSwarm(const vector<GeneParticle> &s) : _pswarm(s)
	{
		_pop_size = _pswarm.size();
		_globalBestPoint = _pswarm[0]._localBestPoint;
		_globalBestFitness = _pswarm[0]._localBestFitness;
		for (int i = 1; i < _pop_size; i++)
		{
			if (_pswarm[i]._localBestFitness < _globalBestFitness)
			{
				_globalBestPoint = _pswarm[i]._localBestPoint;
				_globalBestFitness = _pswarm[i]._localBestFitness;
			}
		}
	}
	PopSwarm(const vector<GeneParticle> &s, const PopSwarm &old) : _pswarm(s)
	{
		_pop_size = _pswarm.size();
		_globalBestPoint = old._globalBestPoint;
		_globalBestFitness = old._globalBestFitness;
	}
	PopSwarm(const int &fcode, const int &swarm_size, const int &n) : _pop_size(swarm_size)		//конструктор случайным образом
	{
		_pswarm.resize(swarm_size);
		_pswarm[0] = GeneParticle(n, fcode);
		_globalBestPoint = _pswarm[0]._localBestPoint;
		_globalBestFitness = _pswarm[0]._localBestFitness;
		for (int i = 1; i < swarm_size; i++)
		{
			_pswarm[i] = GeneParticle(n, fcode);
			if (_pswarm[i]._localBestFitness < _globalBestFitness)
			{
				_globalBestPoint = _pswarm[i]._localBestPoint;
				_globalBestFitness = _pswarm[i]._localBestFitness;
			}
		}
	}

	double calc_av_fitness()const                         //средн¤¤ приспособленность попул¤ции
	{
		double avf = 0;
		for (int i = 0; i<_pop_size; ++i)
			avf = avf + _pswarm[i]._currentFitness;
		return avf / _pop_size;
	}

	double calc_min_fitness()const                                 //минимальна¤ приспособленность в попул¤ции
	{
		double min = _pswarm[0]._currentFitness;
		for (int i = 1; i < _pop_size; i++)
		{
			if (_pswarm[i]._currentFitness < min)
				min = _pswarm[i]._currentFitness;
		}
		return min;
	}

	void popSort()               //сортировка попул¤ции
	{
		sort(_pswarm.begin(), _pswarm.end(), [](GeneParticle &p1, GeneParticle &p2) {return p1._currentFitness < p2._currentFitness; });
	}

	PopSwarm popSelection(const int &fcode, const int &parents_size, const int &tournament_size)const   //селекци¤ (номер функции, число родителей, численность турнира) - турнирный отбор
	{
		vector<GeneParticle> p;
		uniform_int_distribution<int> parents_distr(0, _pop_size - 1);
		for (int i = 0; i<parents_size; ++i)
		{
			int s = parents_distr(gen_gapso);
			double min_f = _pswarm[s]._currentFitness;         //минимальна¤ приспособленность = приспособленность случайной особи
			GeneParticle min_p = _pswarm[s];                  //особь с минимальной приспособленностью
			for (int j = 1; j<tournament_size; ++j)           //нахождение особи с минимальной приспособленностью   
			{
				s = parents_distr(gen_gapso);
				if (min_f>_pswarm[s]._currentFitness)
				{
					min_p = _pswarm[s];
					min_f = _pswarm[s]._currentFitness;
				}
			}
			p.push_back(min_p);  //добавление в вектор родител¤
		}
		PopSwarm h(p, *this);   //составление попул¤ции дл¤ скрещивани¤
		return h;
	}

	PopSwarm popCrossover(const int &fcode, const int &children_size, const int &num_k)     //скрещивание (номер функции, количество детей)
	{
		vector<GeneParticle> p;         //попул¤ци¤
		uniform_int_distribution<int> children_distr(0, _pop_size - 1);
		for (int i = 0; i<children_size; ++i)
		{
			int a = children_distr(gen_gapso);
			int b = children_distr(gen_gapso);
			GeneParticle pa = _pswarm[a];             // 1-¤ случайна¤ особь (из родительской попул¤ции)
			GeneParticle pb = _pswarm[b];             // 2-¤ случайна¤ особь (из родительской попул¤ции)
			GeneParticle c = pa;
			crossoverPoint(fcode, num_k, pa, pb, c);     //скрещивание
			p.push_back(c);                     //добавление полученной особи в вектор p
		}
		PopSwarm q(p, *this);
		return q;
	}

	void popMutation(const int &fcode, const int &mutcode, const double &mutprob)    //мутаци¤ (номер функции, тип мутации, веро¤тность мутации) 
	{
		double s = 0;
		for (int i = 0; i<_pop_size; ++i)
		{
			double coin = distr01(gen_gapso);         //генерируем число от 0 до 1
			if (coin<mutprob)
			{
				double mutradius;
				mutation(fcode, mutcode, mutradius, _pswarm[i]);
				if (_pswarm[i]._localBestFitness < _globalBestFitness)
				{
					_globalBestPoint = _pswarm[i]._localBestPoint;
					_globalBestFitness = _pswarm[i]._localBestFitness;
				}
			}
		}
	}

	PopSwarm popReduction(const int &psize)          //сокращение попул¤ции (выбираем наиболее приспособленных особей)
	{
		if (psize <= _pop_size)
		{
			popSort();
			vector<GeneParticle>p;
			for (int i = 0; i<psize; ++i)
				p.push_back(_pswarm[i]);
			PopSwarm h(p, *this);
			return h;
		}
		else
		{
			return *this;
		}
	}

	vector<double> changeLocalVelocity(GeneParticle &x, const double &a1, const double &a2) // изменение скорости одной частицы из ро¤
	{
		vector<double> v(x._currentPoint.size(), 0);
		double k = distr01(gen_gapso), a = a1 + a2;
		double compression_ratio = 2 * k / (fabs(2 - a - sqrt(a*a - 4 * a)));
		double xi1, xi2;
		for (int i = 0; i < v.size(); ++i)
		{
			xi1 = a1*distr01(gen_gapso)*(x._localBestPoint[i] - x._currentPoint[i]);
			xi2 = a2*distr01(gen_gapso)*(_globalBestPoint[i] - x._currentPoint[i]);
			//v[i] = compression_ratio*(v[i] + xi1 + xi2);
			v[i] = compression_ratio*(x._speed[i] + xi1 + xi2);
		}
		return v;
	}

	point particleMove(GeneParticle &x, const double &a1, const double &a2, const int &fcode)		//движение частицы + штраф
	{
		vector<double> v(x._currentPoint.size());
		v = changeLocalVelocity(x, a1, a2);
		x._speed = v;
		for (int i = 0; i<x._currentPoint.size(); ++i)
		{
			x._currentPoint[i] += v[i];
		}
		x.checkPoint(x, fcode);
		x._currentFitness = testfun::f(x._currentPoint.size(), x._currentPoint, fcode);
		if (x._currentFitness < x._localBestFitness)
		{
			x._localBestPoint = x._currentPoint;
			x._localBestFitness = x._currentFitness;
		}
		return x._currentPoint;
	}

	void swarmMove(const double &a1, const double &a2, const int &fcode)		//движение всего ро¤
	{
		if (a1 + a2 > 4)
		{
			for (int i = 0; i < _pop_size; i++)
			{
				particleMove(_pswarm[i], a1, a2, fcode);
				if (_pswarm[i]._localBestFitness < _globalBestFitness)
				{
					_globalBestPoint = _pswarm[i]._localBestPoint;
					_globalBestFitness = _pswarm[i]._localBestFitness;
				}
			}
		}
		else
			cout << "Error!";
	}

	void checkPopSwarm(const int &fcode)		// штраф всего ро¤ без движени¤
	{
		for (int i = 0; i < _pop_size; i++)
		{
			_pswarm[i].checkPoint(_pswarm[i], fcode);
		}
	}

	friend const PopSwarm operator+(const PopSwarm &p1, const PopSwarm &p2)
	{
		vector<GeneParticle> merger;
		for (const auto& value : p1._pswarm)
			merger.push_back(value);
		for (const auto& value : p2._pswarm)
			merger.push_back(value);
		PopSwarm p(merger);
		return p;
	}

};


#endif