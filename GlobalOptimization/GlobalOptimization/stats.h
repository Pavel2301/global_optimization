#ifndef STAT_H_INCLUDED
#define STAT_H_INCLUDED

#include <iostream>
#include <vector>

using namespace std;

struct stats
{
	double st_min;
	double st_max;
	double st_mean;
	double st_sd;
};

inline stats get_stats(const vector<double> &v)
{
	stats s;
	double sum = 0;
	double sd_sum = 0;
	double s_min = v[0];
	double s_max = v[0];
	int n = v.size();
	for (int i = 0; i < n; ++i)
	{
		sum = sum + v[i];
		if (s_min > v[i])
			s_min = v[i];
		if (s_max < v[i])
			s_max = v[i];
	}
	for (int i = 0; i < n; ++i)
		sd_sum = sd_sum + (v[i] - sum / n)*(v[i] - sum / n);
	s.st_max = s_max;										/**максимальное**/
	s.st_min = s_min;										/**минимальное**/
	s.st_mean = sum / n;									/**среднее**/
	s.st_sd = sqrt(sd_sum / (n - 1));						/**стандартное отклонение**/
	return s;												//возвращаем статистику
}

inline double show_mean(const vector<double> &v)
{
	int n = v.size();
	double sum = 0;
	for (int i = 0; i<n; ++i)
		sum = sum + v[i];
	return sum / n;
}

inline void cout_stats(stats s)
{
	cout << "min" << ";" << point_to_comma(to_string(s.st_min)) << endl;
	cout << "max" << ";" << point_to_comma(to_string(s.st_max)) << endl;
	cout << "mean" << ";" << point_to_comma(to_string(s.st_mean)) << endl;
	cout << "sd" << ";" << point_to_comma(to_string(s.st_sd)) << endl;
}

#endif // STATS_H_INCLUDED