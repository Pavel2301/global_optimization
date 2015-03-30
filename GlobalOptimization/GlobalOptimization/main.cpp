#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <chrono>

#include "sa.h"
#include "ga.h"
#include "difevo.h"
#include "pso.h"
#include "gapso.h"
#include "point.h"
#include "report.h"
#include "testfun.h"
#include "stats.h"
#include "utils.h"

using namespace std;
using namespace testfun;

void check_points_file(const int &fcode, const int &n)
{
	ifstream p_file("points\\" + namefun[fcode] + "-" + to_string(n));
	if (!p_file)
		gen_st_points(fcode, n, 100);
	p_file.close();
}

void print_mean(vector<FullReport> &fr, const int &launches)
{
	long double temp = 0;
	for (int i = 0; i < fr[0].show_report_size(); i++)
	{
		for (int j = 0; j < launches; j++)
			temp += fr[j].show_report()[i].show_value_function();
		temp /= launches;
		cout << point_to_comma(to_string(temp)) << ";";
		temp = 0;
	}
}

void start_sa(const int &fcode, const int &n, const int &launches, const int &max_call_f)
{
	ifstream config_file("config\\sa.txt");
	string s;
	char c;
	double init_t, end_t, rate;
	int inner_iter;
	config_file >> s >> c >> init_t;
	config_file >> s >> c >> end_t;
	config_file >> s >> c >> inner_iter;
	config_file >> s >> c >> rate;
	config_file.close();
	check_points_file(fcode, n);
	ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
	string output_fname = "results\\SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(init_t) + "_" + to_string(end_t) + "_" + to_string(inner_iter) + "_" + to_string(rate) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	vector<FullReport> fr(launches);
	vector<double> final;
	point st(n);
	for (int j = 0; j < n; j++)
		points_file >> st[j];
	fr[0] = SimulatedAnnealing(fcode, n, st, max_call_f, init_t, end_t, inner_iter, rate);
	fr[0].print_full_report_cf();
	fr[0].print_full_report();
	final.push_back(fr[0].show_report().back().show_value_function());
	for (int i = 1; i < launches; i++)
	{
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		fr[i] = SimulatedAnnealing(fcode, n, st, max_call_f, init_t, end_t, inner_iter, rate);
		fr[i].print_full_report();
		final.push_back(fr[i].show_report().back().show_value_function());
	}
	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(init_t) + "_" + to_string(end_t) + "_" + to_string(inner_iter) + "_" + to_string(rate) + ".csv";
	FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
	cout_stats(stts);
	fclose(stats_file);
}

void start_ga(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\ga.txt");
	string s;
	char c;
	int tourn_size, mut_type, parents_size, children_size, num_k;
	double mut_prob;
	config_file >> s >> c >> tourn_size;
	config_file >> s >> c >> mut_type;
	config_file >> s >> c >> mut_prob;
	config_file >> s >> c >> parents_size;
	config_file >> s >> c >> children_size;
	config_file >> s >> c >> num_k;
	config_file.close();
	check_points_file(fcode, n);
	ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
	string output_fname = "results\\GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	vector<FullReport> fr(launches);
	vector<double> final;
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		w.push_back(st);
	}
	fr[0] = ga(fcode, n, w, max_call_f, tourn_size, mut_type, mut_prob, parents_size, children_size, num_k);
	fr[0].print_full_report_cf();
	fr[0].print_full_report();
	final.push_back(fr[0].show_report().back().show_value_function());
	for (int i = 1; i < launches; i++)
	{
		vector<point> w;
		for (int t = 0; t < psize; t++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			w.push_back(st);
		}
		fr[i] = ga(fcode, n, w, max_call_f, tourn_size, mut_type, mut_prob, parents_size, children_size, num_k);
		fr[i].print_full_report();
		final.push_back(fr[i].show_report().back().show_value_function());
	}
	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + ".csv";
	FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
	cout_stats(stts);
	fclose(stats_file);
}

void start_difevo(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\difevo.txt");
	string s;
	char c;
	int alpha_gen;
	double mut_prob;
	config_file >> s >> c >> mut_prob;
	config_file >> s >> c >> alpha_gen;
	config_file.close();
	check_points_file(fcode, n);
	ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
	string output_fname = "results\\DifEvo_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(mut_prob) + "_" + to_string(alpha_gen) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	vector<FullReport> fr(launches);
	vector<double> final;
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		w.push_back(st);
	}
	fr[0] = difevo(fcode, n, w, max_call_f, mut_prob, alpha_gen);
	fr[0].print_full_report_cf();
	fr[0].print_full_report();
	final.push_back(fr[0].show_report().back().show_value_function());
	for (int i = 1; i < launches; i++)
	{
		vector<point> w;
		for (int t = 0; t < psize; t++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			w.push_back(st);
		}
		fr[i] = difevo(fcode, n, w, max_call_f, mut_prob, alpha_gen);
		fr[i].print_full_report();
		final.push_back(fr[i].show_report().back().show_value_function());
	}
	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_DifEvo_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(mut_prob) + "_" + to_string(alpha_gen) + ".csv";
	FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
	cout_stats(stts);
	fclose(stats_file);
}

void start_pso(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\pso.txt");
	string s;
	char c;
	double a1, a2;
	config_file >> s >> c >> a1;
	config_file >> s >> c >> a2;
	config_file.close();
	check_points_file(fcode, n);
	ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
	string output_fname = "results\\PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	vector<FullReport> fr(launches);
	vector<double> final;
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		w.push_back(st);
	}
	fr[0] = pso(fcode, n, w, max_call_f, a1, a2);
	fr[0].print_full_report_cf();
	fr[0].print_full_report();
	final.push_back(fr[0].show_report().back().show_value_function());
	for (int i = 1; i < launches; i++)
	{
		vector<point> w;
		for (int t = 0; t < psize; t++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			w.push_back(st);
		}
		fr[i] = pso(fcode, n, w, max_call_f, a1, a2);
		fr[i].print_full_report();
		final.push_back(fr[i].show_report().back().show_value_function());
	}
	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
	cout_stats(stts);
	fclose(stats_file);
}

void start_gapso(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\ga_pso.txt");
	string s;
	char c;
	int tourn_size, mut_type, parents_size, children_size, num_k;
	double mut_prob, a1, a2;
	config_file >> s >> c >> tourn_size;
	config_file >> s >> c >> mut_type;
	config_file >> s >> c >> mut_prob;
	config_file >> s >> c >> parents_size;
	config_file >> s >> c >> children_size;
	config_file >> s >> c >> num_k;
	config_file >> s >> c >> a1;
	config_file >> s >> c >> a2;
	config_file.close();
	check_points_file(fcode, n);
	ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
	string output_fname = "results\\GA_PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	vector<FullReport> fr(launches);
	vector<double> final;
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		w.push_back(st);
	}
	fr[0] = ga_pso(fcode, n, w, max_call_f, a1, a2, tourn_size, mut_type, mut_prob, parents_size, children_size, num_k);
	fr[0].print_full_report_cf();
	fr[0].print_full_report();
	final.push_back(fr[0].show_report().back().show_value_function());
	for (int i = 1; i < launches; i++)
	{
		vector<point> w;
		for (int t = 0; t < psize; t++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			w.push_back(st);
		}
		fr[i] = ga_pso(fcode, n, w, max_call_f, a1, a2, tourn_size, mut_type, mut_prob, parents_size, children_size, num_k);
		fr[i].print_full_report();
		final.push_back(fr[i].show_report().back().show_value_function());
	}
	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_GAPSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
	cout_stats(stts);
	fclose(stats_file);
}

int main()
{
	string input_fname = "start.txt";
	FILE *input_file = freopen(input_fname.c_str(), "r", stdin);
	string s;
	char c;
	int mcode, fcode, n, launches, psize, max_call_f;
	cin >> s >> c >> mcode;
	cin >> s >> c >> fcode;
	cin >> s >> c >> n;
	cin >> s >> c >> launches;
	cin >> s >> c >> psize;
	cin >> s >> c >> max_call_f;
	fclose(input_file);
	switch (mcode)
	{
	case 0: {start_sa(fcode, n, launches, max_call_f); break; }					//метод имитации отжига
	case 1: {start_ga(fcode, n, launches, psize, max_call_f); break; }			//генетический алгоритм
	case 2: {start_difevo(fcode, n, launches, psize, max_call_f); break; }		//метод дифференциальной эволюции
	case 3: {start_pso(fcode, n, launches, psize, max_call_f); break; }			//метод роя частиц
	case 4: {start_gapso(fcode, n, launches, psize, max_call_f); break; }		
	}
	return 0;
}

