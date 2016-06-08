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
#include "ais.h"
#include "cg.h"
#include "cgm.h"
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
	int i, num_sets;
	config_file >> num_sets;
	vector<double> init_t(num_sets), end_t(num_sets), rate(num_sets);
	vector<int> inner_iter(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> init_t[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> end_t[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> inner_iter[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> rate[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(init_t[p]) + "_" + to_string(end_t[p]) + "_" + to_string(inner_iter[p]) + "_" + to_string(rate[p]) + ".csv";
		FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
		vector<FullReport> fr(launches);
		vector<double> final;
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		fr[0] = SimulatedAnnealing(fcode, n, st, max_call_f, init_t[p], end_t[p], inner_iter[p], rate[p]);
		fr[0].print_full_report_cf();
		fr[0].print_full_report();
		final.push_back(fr[0].show_report().back().show_value_function());
		for (i = 1; i < launches; i++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			fr[i] = SimulatedAnnealing(fcode, n, st, max_call_f, init_t[p], end_t[p], inner_iter[p], rate[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_SimulatedAnnealing_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(init_t[p]) + "_" + to_string(end_t[p]) + "_" + to_string(inner_iter[p]) + "_" + to_string(rate[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_ga(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\ga.txt");
	string s;
	char c;
	int i, num_sets;
	config_file >> num_sets;
	vector<int> tourn_size(num_sets), mut_type(num_sets), parents_size(num_sets), children_size(num_sets), num_k(num_sets);
	vector<double> mut_prob(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> tourn_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_type[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_prob[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> parents_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> children_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> num_k[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size[p]) + "_" + to_string(mut_type[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(parents_size[p]) + "_" + to_string(children_size[p]) + "_" + to_string(num_k[p]) + ".csv";
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
		fr[0] = ga(fcode, n, w, max_call_f, tourn_size[p], mut_type[p], mut_prob[p], parents_size[p], children_size[p], num_k[p]);
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
			fr[i] = ga(fcode, n, w, max_call_f, tourn_size[p], mut_type[p], mut_prob[p], parents_size[p], children_size[p], num_k[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		//print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_GeneticAlgorithm_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size[p]) + "_" + to_string(mut_type[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(parents_size[p]) + "_" + to_string(children_size[p]) + "_" + to_string(num_k[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_difevo(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\difevo.txt");
	string s;
	char c;
	int i, num_sets;
	config_file >> num_sets;
	vector<double> mut_prob(num_sets);
	vector<int> alpha_gen(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_prob[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> alpha_gen[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\DifEvo_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(mut_prob[p]) + "_" + to_string(alpha_gen[p]) + ".csv";
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
		fr[0] = difevo(fcode, n, w, max_call_f, mut_prob[p], alpha_gen[p]);
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
			fr[i] = difevo(fcode, n, w, max_call_f, mut_prob[p], alpha_gen[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_DifEvo_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(mut_prob[p]) + "_" + to_string(alpha_gen[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_pso(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\pso.txt");
	string s;
	char c;
	int i, num_sets;
	config_file >> num_sets;
	vector<double> a1(num_sets), a2(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> a1[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> a2[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(a1[p]) + "_" + to_string(a2[p]) + ".csv";
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
		fr[0] = pso(fcode, n, w, max_call_f, a1[p], a2[p]);
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
			fr[i] = pso(fcode, n, w, max_call_f, a1[p], a2[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(a1[p]) + "_" + to_string(a2[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_gapso(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\ga_pso.txt");
	string s;
	char c;
	int i, num_sets;
	config_file >> num_sets;
	vector<int> tourn_size(num_sets), mut_type(num_sets), parents_size(num_sets), children_size(num_sets), num_k(num_sets);
	vector<double> mut_prob(num_sets), a1(num_sets), a2(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> tourn_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_type[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_prob[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> parents_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> children_size[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> num_k[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> a1[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> a2[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\GA_PSO_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size[p]) + "_" + to_string(mut_type[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(parents_size[p]) + "_" + to_string(children_size[p]) + "_" + to_string(num_k[p]) + "_" + to_string(a1[p]) + "_" + to_string(a2[p]) + ".csv";
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
		fr[0] = ga_pso(fcode, n, w, max_call_f, a1[p], a2[p], tourn_size[p], mut_type[p], mut_prob[p], parents_size[p], children_size[p], num_k[p]);
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
			fr[i] = ga_pso(fcode, n, w, max_call_f, a1[p], a2[p], tourn_size[p], mut_type[p], mut_prob[p], parents_size[p], children_size[p], num_k[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		//print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_GAPSO3_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(tourn_size[p]) + "_" + to_string(mut_type[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(parents_size[p]) + "_" + to_string(children_size[p]) + "_" + to_string(num_k[p]) + "_" + to_string(a1[p]) + "_" + to_string(a2[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_ais(const int &fcode, const int &n, const int &launches, const int &psize, const int &max_call_f)
{
	ifstream config_file("config\\ais.txt");
	string s;
	char c;
	int i, num_sets;
	config_file >> num_sets;
	vector<int> s_kol(num_sets), s1(num_sets), d(num_sets), k(num_sets), beta(num_sets), nc(num_sets), r1(num_sets), mutcode(num_sets);
	vector<double> mut_prob(num_sets);
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> s_kol[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> s1[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> d[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> k[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> beta[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> nc[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> r1[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mut_prob[i];
	config_file >> s >> c;
	for (i = 0; i < num_sets; i++)
		config_file >> mutcode[i];
	config_file.close();
	check_points_file(fcode, n);
	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\ArtificialImmuneSystem_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(s_kol[p]) + +"_" + to_string(s1[p]) + "_" + to_string(d[p]) + "_" + to_string(k[p]) + "_" + to_string(beta[p]) + "_" + to_string(nc[p]) + "_" + to_string(r1[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(mutcode[p]) + ".csv";
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
		fr[0] = ais1(fcode, n, w, max_call_f, s_kol[p], s1[p], d[p], k[p], beta[p], nc[p], r1[p], mut_prob[p], mutcode[p]);
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
			fr[i] = ais1(fcode, n, w, max_call_f, s_kol[p], s1[p], d[p], k[p], beta[p], nc[p], r1[p], mut_prob[p], mutcode[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_ArtificialImmuneSystem_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(s_kol[p]) + +"_" + to_string(s1[p]) + "_" + to_string(d[p]) + "_" + to_string(k[p]) + "_" + to_string(beta[p]) + "_" + to_string(nc[p]) + "_" + to_string(r1[p]) + "_" + to_string(mut_prob[p]) + "_" + to_string(mutcode[p]) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

void start_cg(const int &fcode, const int &n, const int &launches, const int &max_call_f)
{
	check_points_file(fcode, n);
	//for (int p = 0; p < num_sets; p++)
	//{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\cg_FletcherRieves_" + namefun[fcode] + "_" + to_string(n)  + ".csv";
		FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
		std::vector<FullReport> fr(launches);
		std::vector<double> final;
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		fr[0] = cg_FletcherRieves(fcode, n, st, max_call_f);
		fr[0].print_full_report_cf();
		fr[0].print_full_report();
		final.push_back(fr[0].show_report().back().show_value_function());
		for (int i = 1; i < launches; i++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			fr[i] = cg_FletcherRieves(fcode, n, st, max_call_f);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}
		//print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_cg_FletcherRieves_" + namefun[fcode] + "_" + to_string(n) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	//}
}

void start_cgm(const int &fcode, const int &n, const int &launches, const int &max_call_f)
{
	ifstream config_file("config\\cgm.txt");
	string s;
	char c;
	int num_sets;
	config_file >> num_sets;
	std::vector<int> cgm_type(num_sets);
	config_file >> s >> c;
	for (int i = 0; i < num_sets; i++)
		config_file >> cgm_type[i];
	config_file.close();
	check_points_file(fcode, n);

	for (int p = 0; p < num_sets; p++)
	{
		ifstream points_file("points\\" + namefun[fcode] + "-" + to_string(n));
		string output_fname = "results\\cgm_" + to_string(cgm_type[p]) + "_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(launches) + ".csv";
		FILE *output_file = freopen(output_fname.c_str(), "a", stdout);
		vector<FullReport> fr(launches);
		vector<double> final;
		point st(n);
		for (int j = 0; j < n; j++)
			points_file >> st[j];
		fr[0] = cgm(fcode, n, st, max_call_f, rb[fcode], cgm_type[p]);
		fr[0].print_full_report_cf();
		fr[0].print_full_report();
		final.push_back(fr[0].show_report().back().show_value_function());

		for (int i = 1; i < launches; i++)
		{
			point st(n);
			for (int j = 0; j < n; j++)
				points_file >> st[j];
			fr[i] = cgm(fcode, n, st, max_call_f, rb[fcode], cgm_type[p]);
			fr[i].print_full_report();
			final.push_back(fr[i].show_report().back().show_value_function());
		}

		//print_mean(fr, launches);
		points_file.close();
		fclose(output_file);
		stats stts = get_stats(final);
		string stats_fname = "results\\Stats_cgm_" + to_string(cgm_type[p]) + "_" + namefun[fcode] + "_" + to_string(n) + "_" + to_string(launches) + ".csv";
		FILE *stats_file = freopen(stats_fname.c_str(), "a", stdout);
		cout_stats(stts);
		fclose(stats_file);
		final.clear();
	}
}

int main()
{
	string input_fname = "start.txt";
	FILE *input_file = freopen(input_fname.c_str(), "r", stdin);
	string s;
	char c;
	int nm, nf, n, launches, psize, max_call_f;
	cin >> nm >> s >> c;
	vector<int> mcodes(nm);
	for (int i = 0; i < nm; i++)
		cin >> mcodes[i];
	cin >> nf >> s >> c;
	vector<int> fcodes(nf);
	for (int i = 0; i < nf; i++)
		cin >> fcodes[i];
	cin >> s >> c >> n;
	cin >> s >> c >> launches;
	cin >> s >> c >> psize;
	cin >> s >> c >> max_call_f;
	fclose(input_file);
	for (int i = 0; i < nm; i++)
		for (int j = 0; j < nf; j++)
		{
			switch (mcodes[i])
			{
			case 0: {start_sa(fcodes[j], n, launches, max_call_f); break; }				//метод имитации отжига
			case 1: {start_ga(fcodes[j], n, launches, psize, max_call_f); break; }			//генетический алгоритм
			case 2: {start_difevo(fcodes[j], n, launches, psize, max_call_f); break; }		//метод дифференциальной эволюции
			case 3: {start_pso(fcodes[j], n, launches, psize, max_call_f); break; }			//метод роя частиц
			case 4: {start_gapso(fcodes[j], n, launches, psize, max_call_f); break; }
			case 5: {start_ais(fcodes[j], n, launches, psize, max_call_f); break;  }		//метод искусcтвенных имунных систем 
			case 6: {start_cg(fcodes[j], n, launches, max_call_f); break; }
			case 7: {start_cgm(fcodes[j], n, launches, max_call_f); break; }
			}
		}
	return 0;
}

