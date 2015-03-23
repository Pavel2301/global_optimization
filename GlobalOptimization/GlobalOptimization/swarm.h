#ifndef SWARM_H_INCLUDED
#define SWARM_H_INCLUDED

#include "particle.h"

class Swarm
{
private:
	int _swarm_size;                            //размер роя частиц
	vector<Particle> _ps;                       //вектор частиц
	point _globalBestPosition;				    //лучшая позиция из роя
	double _globalBestValue;				    //значение ЦФ в лучшей позиции роя
public:
	Swarm(const vector<Particle> &s) : _ps(s)
	{
		_swarm_size = _ps.size();
		_globalBestPosition = _ps[0]._localBestPosition;
		_globalBestValue = _ps[0]._localBestValue;
		for (int i = 1; i < _swarm_size; i++)
		{
			if (_ps[i]._localBestValue < _globalBestValue)
			{
				_globalBestPosition = _ps[i]._localBestPosition;
				_globalBestValue = _ps[i]._localBestValue;
			}
		}
	}
	Swarm(const int &n, const int &swarm_size, const int &fcode) : _swarm_size(swarm_size)		//конструктор случайным образом
	{
		_ps.resize(swarm_size);
		_ps[0] = Particle(n, fcode);
		_globalBestPosition = _ps[0]._localBestPosition;
		_globalBestValue = _ps[0]._localBestValue;
		for (int i = 1; i < swarm_size; i++)
		{
			_ps[i] = Particle(n, fcode);
			if (_ps[i]._localBestValue < _globalBestValue)
			{
				_globalBestPosition = _ps[i]._localBestPosition;
				_globalBestValue = _ps[i]._localBestValue;
			}
		}
	}

	/*Swarm(const Population &p, const int &fcode, const int &n)
	{
	_swarm_size = p.show_p_size();
	for (int i = 0; i < _swarm_size; i++)
	{
	Particle prt(fcode, p.show_points()[i]);
	_ps.push_back(prt);
	}
	_globalBestPosition = _ps[0]._localBestPosition;
	_globalBestValue = _ps[0]._localBestValue;
	for (int i = 1; i < _swarm_size; i++)
	{
	if (_ps[i]._localBestValue < _globalBestValue)
	{
	_globalBestPosition = _ps[i]._localBestPosition;
	_globalBestValue = _ps[i]._localBestValue;
	}
	}
	}*/

	/*Swarm(const int &swarm_size) :
	_swarm_size(swarm_size)		//изменить размер роя
	{
	_ps.resize(swarm_size);
	}*/

	vector<Particle> show_ps()const                 //вернуть сам рой
	{
		return _ps;
	}

	int show_swarm_size()const		                //вернуть размер роя
	{
		return _swarm_size;
	}

	point show_globalBestPosition()const			//вернуть лучшую позицию роя
	{
		return _globalBestPosition;
	}

	double show_globalBestValue()const		           //вернуть значение ЦФ в лучшей позиции роя
	{
		return _globalBestValue;
	}

	vector<double> changeLocalVelocity(Particle &x, const double &a1, const double &a2) // изменение скорости одной частицы из роя
	{
		vector<double> v(x._currentPosition.size(), 0);
		double k = dc(gen_pso), a = a1 + a2;
		double compression_ratio = 2 * k / (fabs(2 - a - sqrt(a*a - 4 * a)));
		double xi1, xi2;
		for (int i = 0; i < v.size(); ++i)
		{
			xi1 = a1*dc(gen_pso)*(x._localBestPosition[i] - x._currentPosition[i]);
			xi2 = a2*dc(gen_pso)*(_globalBestPosition[i] - x._currentPosition[i]);
			//v[i] = compression_ratio*(v[i] + xi1 + xi2);
			v[i] = compression_ratio*(x._velocity[i] + xi1 + xi2);
		}
		return v;
	}

	point particleMove(Particle &x, const double &a1, const double &a2, const int &fcode)		//движение частицы + штраф
	{
		vector<double> v(x._currentPosition.size());
		v = changeLocalVelocity(x, a1, a2);
		x._velocity = v;
		for (int i = 0; i<x._currentPosition.size(); ++i)
		{
			x._currentPosition[i] += v[i];
		}
		x.checkParticle(x, fcode);
		x._currentValue = testfun::f(x._currentPosition.size(), x._currentPosition, fcode);
		if (x._currentValue < x._localBestValue)
		{
			x._localBestPosition = x._currentPosition;
			x._localBestValue = x._currentValue;
		}
		return x._currentPosition;
	}

	void swarmMove(const double &a1, const double &a2, const int &fcode)		//движение всего роя
	{
		if (a1 + a2 > 4)
		{
			for (int i = 0; i < _swarm_size; i++)
			{
				particleMove(_ps[i], a1, a2, fcode);
				if (_ps[i]._localBestValue < _globalBestValue)
				{
					_globalBestValue = _ps[i]._localBestValue;
					_globalBestPosition = _ps[i]._localBestPosition;
				}
			}
		}
		else
			cout << "Error!";
	}

	void checkSwarm(const int &fcode)		// штраф всего роя без движения
	{
		for (int i = 0; i < _swarm_size; i++)
		{
			_ps[i].checkParticle(_ps[i], fcode);
		}
	}

	double calc_averageValue()const		//среднее значение функции
	{
		double av = 0.0;
		for (int i = 0; i < _swarm_size; i++)
			av += _ps[i]._currentValue;
		return av / _swarm_size;
	}

	double calc_minValue()const		//минимальное значение функции
	{
		double mv = _ps[0]._currentValue;
		for (int i = 1; i < _swarm_size; i++)
			mv = min(mv, _ps[i]._currentValue);
		return mv;
	}

};



#endif