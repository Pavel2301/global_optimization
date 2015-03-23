#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <chrono>
#include <random>
#include <set>
#include <map>
#include <string>

#include "testfun.h"
#include "point.h"

using namespace std;
/*utilities*/
typedef vector<double> dv;

void gen_st_points(const int &fcode, const int &n, const int &psize)
{
	string s = "points\\" + testfun::namefun[fcode] + "-" + to_string(n);
	FILE *file = fopen(s.c_str(), "w");
	for (int i = 0; i < 100; i++)
	{
		for (int t = 0; t < psize; t++)
		{
			point x = pure_random(n, testfun::lb[fcode], testfun::rb[fcode]);
			for (int j = 0; j < n; j++)
				fprintf(file, "%f ", x[j]);
			fprintf(file, "\n");
		}
	}
	fclose(file);
}

void fprintf_vec(FILE *pFile, const vector<int> &v)
{
	for (const auto& value : v)
		fprintf(pFile, "%d;", value);
	fprintf(pFile, "\n");
}
void fprintf_vec(FILE *pFile, const vector<double> &v){
	for (const auto& value : v)
		fprintf(pFile, "%f;", value);
	fprintf(pFile, "\n");
}
template<class T>
string to_string(const T& t)
{
	ostringstream os;
	os << t;
	return os.str();
}
template<typename T>
inline ostream& operator<<(ostream &os, const vector<T> &v)
{
	for (unsigned int i = 0; i<v.size(); ++i)
		os << v[i] << ";";
	return os;
}
template<typename T>
inline ostream& operator<<(ostream &os, const vector<vector<T>> &w)
{
	for (unsigned int i = 0; i<w.size(); ++i)
		os << w[i] << endl;
	return os;
}
template<typename T, typename G>
inline ostream& operator<<(ostream &os, const map<T, G> &m)
{
	for (auto mit = m.begin(); mit != m.end(); ++mit)
		os << "(" << mit->first << ": " << mit->second << ")" << ";";
	return os;
}

inline string point_to_comma(string &s)
{
	bool flag = false;
	for (int i = 0; i<s.size() && !flag; i++)
	{
		if (s[i] == '.')
		{
			s[i] = ',';
			flag = true;
		}
	}
	return s;
}

inline string get_current_prefix(){
	using namespace std::chrono;
	duration<int, std::ratio<60 * 60 * 24> > one_day(1);
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	struct tm * timeinfo;
	char buffer[80];
	time(&tt);
	timeinfo = localtime(&tt);
	strftime(buffer, 80, "%Y%m%d_%H%M_", timeinfo);
	string start_time = to_string(buffer);
	return start_time;
}

template<class T>
inline void vv_linear_cout(const vector<vector<T>> &u){
	for (int i = 0; i<u.size(); ++i){
		for (int j = 0; j<u[i].size(); ++j){
			cout << u[i][i] << ";";
		}
	}
}

template<class T>
inline vector<double> aver_vec(const vector<T> &u, const int &x){
	int k = u.size();
	vector<double> v(k, 0);
	for (int i = 0; i<k; ++i) v[i] = u[i] / x;
	return v;
}

template<class T>
inline void transp_cout(const vector<vector<T>> &w){
	int k = w.size();
	int max_l = 0;

	vector<int>lens(k, 0);
	for (int i = 0; i<k; ++i){
		lens[i] = w[i].size();
		if (lens[i]>max_l) max_l = lens[i];
	}
	cout << lens << endl;
	for (int i = 0; i<max_l; ++i){
		for (int j = 0; j<k; ++j){
			if (i<lens[j]){ cout << w[j][i] << ";"; }
			else{ cout << ";"; }
		}
		cout << endl;
	}
}
#endif // UTILS_H_INCLUDED
