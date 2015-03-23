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
	string output_fname = "results\\SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(init_t) + "_" + to_string(end_t) + "_" + to_string(inner_iter) + "_" + to_string(rate) + ".csv";
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

	/*vector<double> temp;
	for (int i = 0; i < fr[0].show_report_size(); i++)
	{
	for (int j = 0; j < launches; j++)
	temp.push_back(fr[j].show_report()[i].show_value_function());
	cout << point_to_comma(to_string(show_mean(temp))) << ";";
	temp.clear();
	}*/

	print_mean(fr, launches);
	points_file.close();
	fclose(output_file);
	stats stts = get_stats(final);
	string stats_fname = "results\\Stats_SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(init_t) + "_" + to_string(end_t) + "_" + to_string(inner_iter) + "_" + to_string(rate) + ".csv";
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
	string output_fname = "results\\GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + ".csv";
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
	string stats_fname = "results\\Stats_GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + ".csv";
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
	string output_fname = "results\\DifEvo_" + namefun[fcode] + "_" + to_string(mut_prob) + "_" + to_string(alpha_gen) + ".csv";
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
	string stats_fname = "results\\Stats_DifEvo_" + namefun[fcode] + "_" + to_string(mut_prob) + "_" + to_string(alpha_gen) + ".csv";
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
	string output_fname = "results\\PSO_" + namefun[fcode] + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
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
	string stats_fname = "results\\Stats_PSO_" + namefun[fcode] + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
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
	string output_fname = "results\\GA_PSO_" + namefun[fcode] + "_" + to_string(tourn_size) + "_" + to_string(mut_type) + "_" + to_string(mut_prob) + "_" + to_string(parents_size) + "_" + to_string(children_size) + "_" + to_string(num_k) + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
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
	string stats_fname = "results\\Stats_GAPSO_" + namefun[fcode] + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
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
	case 0: {start_sa(fcode, n, launches, max_call_f); break; }
	case 1: {start_ga(fcode, n, launches, psize, max_call_f); break; }
	case 2: {start_difevo(fcode, n, launches, psize, max_call_f); break; }
	case 3: {start_pso(fcode, n, launches, psize, max_call_f); break; }
	case 4: {start_gapso(fcode, n, launches, psize, max_call_f); break; }
	}

	//const int n = 100;
	//const int k = 3;
	//const int launches = 100;
	//const int psize = 10;

	/*string input_fname = namefun[k] + "-" + to_string(n);
	FILE *input_file = freopen(input_fname.c_str(), "r", stdin);
	vector<FullReport> fr(launches);
	for (int i = 0; i < launches; i++)
	{
	ifstream f("input_sa.txt");
	point st(n);
	for (int j = 0; j < n; j++)
	cin >> st[j];
	while (!f.eof())
	{
	double t1, t2, rate;
	int inner_iter;
	f >> t1 >> t2 >> inner_iter >> rate;
	fr[i] = SimulatedAnnealing(k, n, st, 10000, t1, t2, inner_iter, rate);
	string output_fname = "SimulatedAnnealing_" + namefun[k] + "_" + to_string(t1) + "_" + to_string(t2) + "_" + to_string(inner_iter) + "_" + to_string(rate) + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	fr[i].print_full_report();
	}
	f.close();
	}*/

	/*string input_fname = namefun[k] + "-" + to_string(n) + "-pop20";
	FILE *input_file = freopen(input_fname.c_str(), "r", stdin);
	vector<vector<FullReport>> fr(10);

	for (int i = 0; i < launches; i++)
	{
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
	point st(n);
	for (int j = 0; j < n; j++)
	{
	cin >> st[j];
	}
	w.push_back(st);
	}
	ifstream f("input.txt");
	for (int q = 0; q < 10; q++)
	{
	double a1, a2;
	f >> a1 >> a2;
	//fr[q].push_back(pso(k, n, w, 5000, a1, a2));
	//string output_fname = "PSO_" + namefun[k] + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	fr[q].push_back(ga_pso(k, n, w, 5000, a1, a2, 3, 3, 0.9, 5, 5, 1, 3));
	string output_fname = "GA_PSO2" + namefun[k] + "_" + to_string(a1) + "_" + to_string(a2) + ".csv";
	//fr[q].push_back(ga(k, n, w, 5000, 3, 1, 0.9, 5, 5, 1, 3));
	//string output_fname = "GA_" + namefun[k] + ".csv";
	FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
	fr[q][i].print_full_report();
	}
	f.close();
	}

	/*const int swarm_size = 10;
	for (int i = 0; i < 10; i++)
	{
	pso(k, n, swarm_size, 1000, 2, 5);
	}*/

	/*ifstream f("input.txt");
	//const int psize = 20;
	int psize = 10;
	double y;
	//vector<double> sa_v(100);
	//string s = namefun[k] + "-" + to_string(n);
	string s;
	f >> s;
	if (s == "de")
	{
	int k;
	f >> k;
	//for (int z = 0; z < 2; z++)
	//{
	string s = namefun[k] + "-" + to_string(n) + "-pop20";
	FILE *file = freopen(s.c_str(), "r", stdin);
	for (int i = 0; i < 100; i++)
	{
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
	point st(n);
	for (int j = 0; j < n; j++)
	{
	cin >> st[j];
	}
	w.push_back(st);
	}
	difevo(k, n, w, 10000, 0.9, 1);
	}
	//psize += 10;
	//}
	}
	if (s == "ga")
	{
	int k;
	f >> k;
	//for (int z = 0; z < 2; z++)
	//	{
	string s = namefun[k] + "-" + to_string(n) + "-pop20";
	FILE *file = freopen(s.c_str(), "r", stdin);
	for (int i = 0; i < 100; i++)
	{
	vector<point> w;
	for (int t = 0; t < psize; t++)
	{
	point st(n);
	for (int j = 0; j < n; j++)
	{
	cin >> st[j];
	}
	w.push_back(st);
	}

	//ga(k, n, w, 1000, 3, 1, 0.9, 5, 5, 10, 2, 0);
	//ga(k, n, w, 1000, 3, 2, 0.9, 5, 5, 10, 2, 0);
	//ga(k, n, w, 1000, 3, 3, 0.9, 5, 5, 10, 2, 0);
	/*ga(k, n, w, 10000, 3, 1, 0.7, 5, 5, 10, 2, 0);
	ga(k, n, w, 10000, 3, 2, 0.7, 5, 5, 10, 2, 0);
	ga(k, n, w, 10000, 3, 3, 0.7, 5, 5, 10, 2, 0);
	ga(k, n, w, 10000, 3, 1, 0.5, 5, 5, 10, 2, 0);
	ga(k, n, w, 10000, 3, 2, 0.5, 5, 5, 10, 2, 0);
	ga(k, n, w, 10000, 3, 3, 0.5, 5, 5, 10, 2, 0);*/

	// ga(k, n, w, 1000, 3, 2, 0.9, 5, 5, 10, 3, 0);
	// ga(k, n, w, 1000, 3, 3, 0.9, 5, 5, 10, 3, 0);
	// ga(k, n, 10, 1000, 3, 3, 0.9, 5, 5, 10, 3, 0);
	//ga(k, n, w, 300000, 3, 1, 0.9, 5, 5, 10, 3, 0);


	//}
	//psize += 10;
	//	}

	//}
	//stats s_sa = get_stats(sa_v);
	//freopen("stat.csv", "a", stdout);
	//cout_stats(s_sa);
	//ga(0, 10, 20, 1000, 3, 1, 0.3, 10, 10, 10, 0);
	return 0;
}

