#ifndef TESTFUN_H_INCLUDED
#define TESTFUN_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <thread> 
#include <Windows.h>
#include <direct.h>

#include <random>
#include <chrono>
#include <vector>

using namespace std;
typedef vector<double> point;

#pragma once
namespace testfun
{
	static int call_f = 0; /*количество вызовов ÷‘*/
	const double pi = 3.1415926535897932384626433832795028841971;
	const double pi2 = 2 * pi;
	const double M_E = 2.71828182845904523536;		
	const vector<double> lb = { -5.12, -2.048, -1.28, -5.12, -500, -600, -30, -2, -5, -100, -10, -1, -1, -512, -10, -10, -10, -500, -pi, -10 ,-3}; /*левые границы области определения*/
	const vector<double> rb = {  5.12,  2.048,  1.28,  5.12,  500,  600,  30,  2,  5,  100,  10,  4,  1,  512,  10,  10,  10,  500,  pi,  10 , 3}; /*правые границы бласти определения*/ /*=> область определения*/
	const vector<string> namefun = { "Sphere", "Rosenbrock", "Gauss Quartic", "Rastrigin", "Schwefel", "Griewank", "Ackley", "RastriginNovgorod", "StepFunction", "Salomon", "Alpine1", "Brown", "Deb1", "Egg Holder", "Pinter", "Shubert3", "Streched V Sine Wave", "Trigonometric2", "Wavy", "Whitley", "Morse" };/*имена функций*/
	const vector<double> gauss_vec = { 0.0914399, -0.146639, 0.108227, 1, 0.598773, -1, 0.250188, -0.88703, -0.966537, -1, -0.526993, 0.160174, -1, -0.884465, 0.986888, -0.707977, 1, -1, 1, 0.623164, 1, -1, 0.831348, -1, 0.00521962, -0.066235, 1, 0.171589, 0.195258, -1, -0.00995592, -1, 0.395942, 1, 0.300302, -1, 1, 1, 0.996736, -1, -1, -0.91945, -0.777617, -0.284741, 0.0375635, -1, 0.720894, -0.473548, 0.420251, -1, 0.488805, -1, -1, -0.0976831, -0.206362, -0.396622, 0.509816, 0.0464454, -0.0725463, 0.752187, -0.060867, -0.00975769, -1, 0.505615, 0.278716, 0.559747, -0.893439, 1, 0.468256, 0.467339, -0.167235, 0.842858, -1, -1, 0.0635688, 1, 0.206896, -0.912863, 0.694523, -0.743222, 0.778778, 0.257005, 0.763326, 0.72237, 0.886368, -0.57182, -0.954208, 0.296418, 1, -0.127915, -0.241335, 0.27632, 0.2147, -0.595802, 0.789754, -0.3903, -0.632689, -1, 1, 0.728787, -1, 1, 0.359812, -0.794921, 1, -0.488669, 1, 0.718143, -0.979901, 0.667572, -0.881618, 0.573863, -0.508577, 0.701724, 0.449535, -1, -1, 0.0241156, 0.763896, -0.437631, -0.238105, 0.644866, 0.965775, 1, -0.676568, 0.812245, 0.00700361, -0.992248, 0.903787, 1, -0.970574, 0.843506, -0.34948, -0.211436, 1, -1, 0.879992, -1, 1, 0.00640636, 0.346083, 1, -0.466172, 1, 0.964828, 1, -0.620768, 0.243794, -0.366889, -0.845317, 0.558504, 0.952914, -0.737679, -0.88956, -0.607622, 1, 1, 1, -0.548275, 0.729677, 1, -0.7853, -0.0057842, -0.161104, 0.851031, 0.165215, 0.98558, 0.986174, -0.209502, -1, -0.0375259, -0.724831, -0.777656, -0.838853, 1, 0.220551, 0.383534, 0.477739, -0.0307834, -0.686389, -1, 1, -0.77841, -1, -0.668819, 0.486468, 0.816184, 0.609682, 0.53106, 1, -1, -1, 1, -0.815723, -1, -0.452409, -0.560971, -0.217311, 1, -0.449733, 0.447374, 0.850261, -0.385728, 0.525265, 1, 0.787643, -1, 0.661059, 0.422755, 1, -1, -0.20034, 0.398708, 0.912971, -0.440297, -1, 0.526395, 1, 0.827449, -0.523215, -0.405165, -0.811864, -0.354357, 0.272777, 0.615503, 0.0181488, -0.813896, -1, -0.491381, 0.836748, 0.215377, -0.10555, -0.17441, 0.743108, 0.97828, -0.764068, 0.421553, -0.51738, -0.924367, 0.551454, -1, 1, 1, 0.0871906, 0.469007, 0.93099, 1, -0.953593, -0.681067, 0.116124, 0.317598, -1, 1, -1, 0.551683, -1, -1, -0.428039, 0.0465474, -0.314273, -0.223779, 0.138863, 1, 1, -0.127786, -0.11236, -1, 0.522348, -0.284607, -0.418318, -0.141575, 1, 0.471972, -0.178323, 0.467209, -1, 0.0683227, 0.302949, 0.226416, 1, -0.365322, 1, 0.259011, -0.178799, -1, -0.950979, -1, 0.750239, 1, -0.755656, -0.01015, -0.0734345, 0.761021, -1, 0.946215, 0.000836779, 0.617069, 0.0891252, -0.823381, 1, 0.212234, -0.680532, -0.845184, 1, 0.0143709, 0.422468, 0.161331, -0.446941, 1, -0.713956, 0.55713, 0.268585, -1, 1, -1, 0.846447, -0.304294, -0.889344, -0.0321601, 0.0011341, -0.980085, -1, -1, -0.0303467, 0.148892, -0.464188, -1, -0.394032, 1, -1, -1, -0.138234, 0.570091, -0.352299, -0.512156, 1, 0.344292, -0.247441, 0.531572, 0.198788, 1, -1, -1, 1, -0.944122, -0.842723, -1, -1, 1, 0.231871, 0.289454, -0.00760637, 0.244116, 0.1894, 0.794509, -0.715875, 0.76217, -0.454719, -1, -0.69049, -0.243067, -1, -0.20672, 0.483317, -0.0185623, -0.689251, 0.745563, 0.427345, -0.951423, -0.928303, 1, -0.0985651, 1, 0.422297, 1, 0.298789, -0.990126, -1, -0.0246807, -0.555585, -1, -1, 0.558917, 0.351116, -0.748268, 0.546738, -1, -0.409396, -1, -1, -0.948496, 1, -1, -1, -0.707142, 1, 0.503063, -0.705092, 1, 1, -1, 0.533431, -0.949498, -1, 0.602081, -1, 0.180133, 1, 0.706843, -1, 0.892485, -0.293094, -0.696547, -0.0231491, -0.147953, -1, -0.8831, -0.0544744, 0.521512, -1, 0.840021, 0.832927, 0.543398, 0.551706, -1, -0.846502, -1, 0.977825, -1, -0.303345, 0.159131, 1, 1, -0.874298, 0.417851, -1, -1, -1, 0.426164, 0.259622, 0.25088, 0.0138074, 1, -1, 0.848442, -0.659153, -0.0280992, 0.20606, 0.784608, -0.263084, -0.0863187, -1, 0.521503, -1, -0.0631675, 1, 1, -0.291293, 1, 0.343724, 0.823164, 1, 0.220095, -0.406154, 0.301641, 0.00749308, 1, 0.555924, 1, -0.507266, -0.459688, 0.0184457, -0.473926, -1, 0.217024, 1, -1, 1, -1, 0.851793, 1, -1, 0.500908, -0.727952, -0.668661, -0.389491, -1, 0.275138, -0.752777, 0.3272, -1, 0.687884, -0.403383, -1, 1, -1, 0.586562, -0.350886, -0.516666, 0.82937, -0.647762, -1, 0.834452, -0.188635, -1, -1, 0.673823, 0.95759, -0.385446, 0.98135, -0.926381, -0.0453816, 0.671266, -1, -1, 0.371181, 1, -1, 1, 0.445529, -1, -1, 1, -1, 0.113333, 1, -1, -0.887383, -0.0149785, 0.360733, -0.461802, 0.661011, 0.678322, 1, -0.281698, -0.312592, -1, 0.214917, 0.475159, 1, 0.977707, -0.023947, -0.287094, 0.214615, 1, -1, -0.401195, -1, 0.280923, 1, -1, -0.503539, 0.611433, 1, -0.156979, 0.0270776, 0.461814, -0.332539, -0.0852935, 0.252258, -0.743706, -0.952287, -0.822598, 1, -0.75661, -0.917479, -1, -0.226822, 0.823081, 1, -0.384529, -0.722017, -1, 1, -0.364561, -0.431904, 1, -0.701875, -1, -0.236203, -0.477144, -0.108364, 1, -0.386026, -0.31065, 1, -1, 0.694604, 0.528554, 0.857443, 0.75115, 0.796866, -0.208725, 0.764144, 1, 0.978808, -1, 0.58187, -1, -1, 0.743453, 1, -0.511723, -0.168771, 0.562232, 0.588564, -0.641794, -0.291478, 0.956956, 1, 1, 0.052341, 0.617972, 0.305045, 1, 0.00183223, -0.798413, 0.48425, 1, -0.372767, -0.779387, 1, 1, -0.255661, -0.651594, 1, -0.537177, -0.237078, -0.419291, -0.420331, 0.171538, 0.529336, -0.432933, -0.479574, 1, 0.788725, -0.312968, 1, -0.320641, 0.748226, -0.141423, 1, -0.941076, -1, -1, -1, 1, 0.90258, 0.759454, -0.369247, 0.710074, 1, -0.539732, 0.181516, -1, -0.0744611, 1, 0.486817, -0.688828, 0.934153, -0.484456, 0.0766606, -0.266639, -0.0172779, 0.605863, 1, -1, -1, 1, 0.136194, -1, 1, 1, -1, -0.00105076, -0.489786, -0.994996, -0.529069, 1, -0.290865, -0.195814, 0.806083, 1, -1, 0.317196, -0.814989, 0.582772, -0.538554, 0.330113, 0.352756, 0.696136, -0.836488, 0.0121067, 1, 0.0848788, -0.807968, -0.314871, 1, -0.589831, -0.244633, 0.389061, 0.229403, -0.820473, -0.141237, 1, -1, 0.522161, 1, -0.434531, 0.656042, 0.651048, -1, 0.632303, -1, 0.238404, 0.261852, -0.608605, -0.090066, 0.331811, 0.216657, -0.121477, 0.355319, -0.530069, 1, -1, 0.479489, 0.563217, 0.428272, 0.27922, -0.264701, 0.609147, -0.00158745, -0.629053, 1, 0.983194, -0.064107, 1, 1, -1, -1, -1, -1, -1, -0.849376, -0.718464, -0.512918, 0.0380684, -1, -1, -0.213472, 0.596439, 1, -0.643606, -0.0638827, 0.6567, 0.0879948, 1, 0.432157, -0.541492, 1, -0.832868, 0.177751, -0.358543, 1, 1, -1, 0.837707, 0.241371, 0.778553, -1, 0.669127, 1, 0.29723, 0.741416, 0.139058, -0.912738, 0.0103664, -0.421024, -1, 1, 0.277515, -0.695842, -0.86309, -0.461908, 1, 0.10819, 0.645883, 1, -0.174938, -0.0809417, 0.694216, 0.472533, -0.350508, -0.808086, 1, 1, -0.475649, 0.881734, 1, -0.926209, 1, 0.249573, -0.796845, -1, -0.573845, 0.326703, 0.985515, -0.532394, -1, -0.0799123, -0.307396, 0.658809, 1, 0.10999, -1, -0.86283, -0.233218, 0.448796, 0.137781, -1, -0.155427, -0.464516, 1, -0.622774, -1, -0.305236, 0.629079, 0.54589, 0.281692, 0.0755439, 0.985078, -0.130483, 0.22512, 0.274668, 1, 0.682897, -0.984689, 0.538878, 0.27727, 1, 0.415508, 0.970822, 0.626247, 0.00144609, 1, 0.901946, -0.53738, 0.770024, 0.194323, -0.730869, 0.0731729, 0.629282, 1, -0.381206, -0.424663, -0.157238, 0.159857, 0.0444366, -0.784835, -1, -0.161918, 0.823979, 0.264592, 1, -0.0449488, -0.383353, -0.638476, 0.132545, -1, 0.0690439, -0.650745, -0.827241, 1, 0.415598, -0.379646, -0.985345, 0.309327, -1, 0.522755, 0.936742, -0.258587, -0.626271, 0.616002, -0.637355, -1, -0.692175, 0.315911, 0.385198, 0.26402, 0.667137, 0.252549, 0.536358, -0.672721, 1, -0.0167469, 0.464638, -1, -0.330613, 1, -1, 0.650636, -0.88624, 0.76483, 0.940916, -0.319173, 0.434664, 1, 1, -1, -0.370358, -1, 0.29734, 1, 1, -1, -1, -0.061959, 1, 0.835302, -1, 0.46406, 0.564174, 0.594174, 0.30045, 0.365947, -0.931851, -1, 0.67355, 0.0864398, 0.112111, -0.154632, 1, 0.0109375, 0.386785, -0.137266, -1, -0.429971, 0.786412, -1, -0.819866, 0.555753, 0.124084, -1, 0.422358, 0.656671, 0.341794, 0.0241202, 1, -1, -0.480025, -0.327384, 0.353583, -0.650522, -1, -0.70334, -0.855519, -0.672165, -0.0807826, 0.553136, 1, -1, -1, 0.911406, 0.986126, -1, 1, 1, -0.847747, -0.126584, 0.746313, 0.45466, -0.464976, -0.968192, -1, 0.935406, -1, -1, 0.365578, 0.0835101, 0.0473013, -0.0714998, 0.917451, 0.0348876, -0.401461, 0.144776, -0.433703, -0.869417, -0.924686, 0.378883, 1, -0.594084, 0.904527, 0.251976, -0.268591, 0.0202642, 1, -0.526496, 0.585918, 0.196808, 0.652261 };
	static const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	static mt19937 gen(seed);
	static normal_distribution < double > gauss_distr(0.0, 1.0);
	double morse_rho = 3;

	/* нормальное распределение - распределение Гаусса */
	double gauss()
	{
		double t = gauss_distr(gen);
		if (t>1.0) { return 1.0; }
		else if (t<-1.0) { return -1.0; }
		else { return t; }
	}
	double gauss(const int &i)
	{
		return gauss_vec[i];
	}

	double f_sphere(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += x[i] * x[i];
		return s;
	}

	double f_rosenbrock(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n - 1; ++i)
		{
			double k1 = x[i + 1] - x[i] * x[i];
			double k2 = x[i] - 1;
			s += (100 * k1*k1 + k2*k2);
		}
		return s;
	}

	double f_gauss(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)                        
		{
			double t = x[i] * x[i];
			s += (t*t + gauss(i));
		}
		return s;
	}

	double f_rastrigin(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)                        
		{
			s += (x[i] * x[i] - 10 * cos(pi2*x[i]));
		}
		return 10 * n + s;
	}

	double f_schwefel(const int &n, const point &x)
	{
		double s = 0;
		bool flag = true;
		for (int i = 0; i < n; ++i)
		{
			if ((x[i]>testfun::rb[4]) || (x[i]<testfun::lb[4]))
				flag = false;                 /*flag=false, если x[i] вышла за область определения*/
		}
		if (flag)
		{
			for (int i = 0; i < n; ++i)
				s += (-x[i] * sin(pow(abs(x[i]), 0.5)));
			return 418.9829*n + s;
		}
		else
		{
			return 418.9829*n + 500;
		}
	}

	double f_griewank(const int &n, const point &x)
	{
		double k = 1;
		for (int i = 0; i < n; ++i)
			k *= cos(x[i] / sqrt(i + 1));
		double s = 0;
		for (int i = 0; i < n; ++i)       //!
			s += (x[i] * x[i] / 4000);
		return 1 + s - k;
	}

	void exponents_for_ackley(const int &n, const point &x, double &a, double &b)
	{  /*вычисляет степени для функции Ackley*/
		double s1 = 0;
		double s2 = 0;
		for (int i = 0; i < n; ++i)
		{
			s1 += x[i] * x[i];
			s2 += cos(pi2*x[i]);
		}
		a = sqrt(s1 / n);
		b = s2;
	}

	double f_ackley(const int &n, const point &x)
	{
		double s1, s2;
		exponents_for_ackley(n, x, s1, s2);
		double a = -0.2*s1; /*степень*/
		double b = s2 / n; /*степень*/
		return 20 + exp(1) - 20 * exp(a) - exp(b);
	}

	double f_rastrigin_novgorod(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += (x[i] * x[i] - cos(18 * x[i] * x[i]));
		return n + s;
	}

	double f_step(const int &n, const point &x)
	{
		double s = 0;
		double s1 = 0;
		for (int i = 0; i < n; ++i)
			s1 += fabs(trunc(x[i]));
		if (s1 != 0)
		{
			for (int i = 0; i < n; ++i)
				s += trunc(x[i])*trunc(x[i]);
		}
		else
		{
			for (int i = 0; i < n; ++i)
				s += fabs(x[i]);
			s -= 1;
		}
		return s;
	}

	double f_salomon(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += x[i] * x[i];
		return 1 - cos(pi2 * sqrt(s)) + 0.1 * sqrt(s);
	}

	double f_alpine1(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += fabs(x[i] * sin(x[i]) + 0.1 * x[i]);
		return s;
	}

	double f_alpine2(const int &n, const point &x)
	{
		double s = 1;
		for (int i = 0; i < n; ++i)
			s *= sqrt(x[i])*sin(x[i]);
		return s;
	}

	double f_brown(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n - 1; ++i)
			s += pow(x[i] * x[i], x[i + 1] * x[i + 1] + 1) + pow(x[i + 1] * x[i + 1], x[i] * x[i] + 1);
		return s;
	}

	double f_csendes(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += pow(x[i], 6)*(2 + sin(1 / x[i]));
		return s;
	}

	double f_deb1(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += pow(sin(5 * pi * x[i]), 6);
		s *= -1;
		return s / n;
	}

	double f_deb3(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += pow(sin(5 * pi * (pow(x[i], 0.75) - 0.05)), 6);
		s *= -1;
		return s / n;
	}

	double f_egg_holder(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n - 1; ++i)
			s += -1 * (x[i + 1] + 47) * sin(sqrt(fabs(x[i + 1] + x[i] / 2 + 47))) - x[i] * sin(sqrt(fabs(x[i] - x[i + 1] - 47)));
		return s;
	}

	double f_pinter(const int &n, const point &x)
	{
		double a = x[n - 1] * sin(x[0]) + sin(x[1]);
		double b = x[n - 1] * x[n - 1] - 2 * x[0] + 3 * x[1] - cos(x[0]) + 1;
		double s1 = 0;
		for (int i = 0; i < n; ++i)
			s1 += (i + 1) * x[i] * x[i];
		double s2 = 0;
		s2 += 20 * sin(a) * sin(a);
		for (int i = 1; i < n - 1; ++i)
		{
			a = x[i - 1] * sin(x[i]) + sin(x[i + 1]);
			s2 += 20 * (i + 1) * sin(a) * sin(a);
		}
		a = x[n - 2] * sin(x[n - 1]) + sin(x[0]);
		s2 += 20 * n * sin(a) * sin(a);
		double s3 = 0;
		s3 += log10(1 + b * b);
		for (int i = 1; i < n - 1; ++i)
		{
			b = x[i - 1] * x[i - 1] - 2 * x[i] + 3 * x[i + 1] - cos(x[i]) + 1;
			s3 += (i + 1) * log10(1 + (i + 1) * b * b);
		}
		b = x[n - 2] * x[n - 2] - 2 * x[n - 1] + 3 * x[0] - cos(x[n - 1]) + 1;
		s3 += n * log10(1 + n * b * b);
		return s1 + s2 + s3;
	}

	double f_shubert3(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < 5; ++j)
				s += (j + 1)*sin((j + 2)*x[i] + j + 1);
		return s;
	}

	double f_shubert4(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < 5; ++j)
				s += (j + 1)*cos((j + 2)*x[i] + j + 1);
		return s;
	}

	double f_streched_v_sine_wave(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n - 1; ++i)
			s += pow(x[i + 1] * x[i + 1] + x[i] * x[i], 0.25)*(sin(50 * pow(x[i + 1] * x[i + 1] + x[i] * x[i], 0.1))* sin(50 * pow(x[i + 1] * x[i + 1] + x[i] * x[i], 0.1)) + 0.1);
		return s;
	}

	double f_trigonometric2(const int &n, const point &x)
	{
		double s = 0;
		double a = 6 * sin(14 * (x[0] - 0.9) * (x[0] - 0.9)) * sin(14 * (x[0] - 0.9) * (x[0] - 0.9));
		for (int i = 0; i < n; ++i)
			s += 8 * sin(7 * (x[i] - 0.9) * (x[i] - 0.9)) * sin(7 * (x[i] - 0.9) * (x[i] - 0.9)) + a + (x[i] - 0.9) * (x[i] - 0.9);
		return 1 + s;
	}

	double f_wavy(const int &n, const point &x)
	{
		double s = 0;
		double k = 10;        //?
		for (int i = 0; i < n; ++i)
			s += cos(k * x[i]) * exp(-x[i] * x[i] / 2);
		return 1 - s / n;
	}

	double f_whitley(const int &n, const point &x)
	{
		double s = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				s += pow(100 * (x[i] * x[i] - x[j]) * (x[i] * x[i] - x[j]) + (1 - x[j]) * (1 - x[j]), 2) / 4000 - cos(100 * (x[i] * x[i] - x[j]) * (x[i] * x[i] - x[j]) + (1 - x[j]) * (1 - x[j]) + 1);
		return s;
	}
	
	double rij(const int i, const int j, const point &x)	
	{
		return sqrt(
			(x[3 * i] - x[3 * j])*(x[3 * i] - x[3 * j]) +
			(x[3 * i + 1] - x[3 * j + 1])*(x[3 * i + 1] - x[3 * j + 1]) +
			(x[3 * i + 2] - x[3 * j + 2])*(x[3 * i + 2] - x[3 * j + 2])
			);
	}

	/* Morse metod */
	double grad_fMorse(const int i, const int j, const point &x, int dxId)
	{
		double F1;

		if (dxId == 3 * i)
			F1 = (2 * x[3 * i] - 2 * x[3 * j]);
		else 
			if (dxId == 3 * j)
				F1 = (2 * x[3 * j] - 2 * x[3 * i]);
		else 
			if (dxId == 3 * i + 1)
				F1 = (2 * x[3 * i + 1] - 2 * x[3 * j + 1]);
		else 
			if (dxId == 3 * j + 1)
				F1 = (2 * x[3 * j + 1] - 2 * x[3 * i + 1]);
		else 
			if (dxId == 3 * i + 2)
				F1 = (2 * x[3 * i + 2] - 2 * x[3 * j + 2]);
		else 
			if (dxId == 3 * j + 2)
				F1 = (2 * x[3 * j + 2] - 2 * x[3 * i + 2]);
		else 
			return 0;

		double Kor = rij(i, j, x);
		double epow = pow(M_E, morse_rho - morse_rho*Kor);

		if (Kor != 0)
		{
			return -1.5*F1*epow*(epow - 2.) / Kor - 1.5*epow*epow*F1 / Kor;
		}
		else
		{
			return 1e9;
		}		
	}

	double fMorse(const int i, const int j, const point &x)
	{
		double Kor = rij(i, j, x);
		double epow = pow(M_E, morse_rho - morse_rho*Kor);
		return epow*(epow - 2.);
	}

	double sumMorse(const int &n, const point &x, int dxId)
	{
		double sum = 0.;
		for (int i = 0; i < n/3; i++)
		{
			for (int j = i + 1; j < n/3; j++)
			{
				if (dxId == -1)
					sum += fMorse(i, j, x);
				else 
					sum += grad_fMorse(i, j, x, dxId);
			}
		}
		return sum;
	}

	double f_morse(const int &n, const point &x)
	{		
		return sumMorse(n, x, -1);
	}	

	double f(const int &n, const point &x, const int &fcode)
	{
		++call_f;
		switch (fcode)
		{
		case 0: return f_sphere(n, x);
		case 1: return f_rosenbrock(n, x);
		case 2: return f_gauss(n, x);
		case 3: return f_rastrigin(n, x);
		case 4: return f_schwefel(n, x);
		case 5: return f_griewank(n, x);
		case 6: return f_ackley(n, x);
		case 7: return f_rastrigin_novgorod(n, x);
		case 8: return f_step(n, x);
		case 9: return f_salomon(n, x);
		case 10: return f_alpine1(n, x);
		case 11: return f_brown(n, x);
		case 12: return f_deb1(n, x);
		case 13: return f_egg_holder(n, x);
		case 14: return f_pinter(n, x);
		case 15: return f_shubert3(n, x);
		case 16: return f_streched_v_sine_wave(n, x);
		case 17: return	f_trigonometric2(n, x);
		case 18: return f_wavy(n, x);
		case 19: return f_whitley(n, x);
		case 20: return f_morse(n, x);
		default: return 0;
		}
	}

	point grad_sphere(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		for (int i = 0; i < n; ++i)
			gr[i] = 2 * x[i];
		return gr;
	}

	point grad_rosenbrock(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0), sqx(n, 0);
		for (int i = 0; i < n; ++i)
			sqx[i] = x[i] * x[i];
		gr[0] = x[0] * (-400 * x[1] + 400 * sqx[0] + 2) - 2;
		for (int i = 1; i < n - 1; ++i)
			gr[i] = 100 * (x[i] * (1 - 4 * (x[i + 1] - sqx[i]) - sqx[i - 1])) - 2 * (1 - x[i]);
		gr[n - 1] = 200 * (x[n - 1] - sqx[n - 2]);
		return gr;
	}

	point grad_gauss(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		return gr;
	}

	point grad_rastrigin(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		for (int i = 0; i < n; ++i)
			gr[i] = 2 * x[i] + 10 * pi2 * sin(pi2 * x[i]);
		return gr;
	}

	point grad_schwefel(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0), sqrtx(n, 0);
		for (int i = 0; i < n; ++i)
			sqrtx[i] = sqrt(fabs(x[i]));
		for (int i = 0; i < n; ++i)
			gr[i] = -sin(sqrtx[i]) - sqrtx[i] * cos(sqrtx[i]) / 2;
		return gr;
	}

	point grad_griewank(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0), sqrti(n, 0), cxi(n, 0);
		for (int i = 0; i < n; ++i)
			sqrti[i] = sqrt(double(i + 1));
		for (int i = 0; i < n; ++i)
			cxi[i] = cos(x[i] / sqrti[i]);
		double s = 1;
		for (int i = 0; i < n; ++i)
			s *= cxi[i];
		for (int i = 0; i < n; ++i)
			gr[i] = x[i] / 2000 + (s / cxi[i])*(sin(x[i] / sqrti[i]) / sqrti[i]);
		return gr;
	}

	point grad_ackley(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		double s1, s2;
		exponents_for_ackley(n, x, s1, s2);
		double a = -0.2 * s1;
		double b = s2 / n;
		for (int i = 0; i<n; ++i) 
			gr[i] = (4 * exp(a)*x[i] / s1 + pi2*exp(b)*sin(pi2*x[i])) / n;
		return gr;
	}

	point grad_rastrigin_novgorod(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		for (int i = 0; i < n; ++i)
			gr[i] = 2 * x[i] + 36 * x[i] * sin(18 * x[i] * x[i]);
		return gr;
	}


	point grad_salomon(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		double s = 0;
		for (int i = 0; i < n; ++i)
			s += x[i] * x[i];
		s = sqrt(s);
		for (int i = 0; i < n; ++i)
			gr[i] = (0.1 * x[i]) / s + (pi2 * x[i] * sin(pi2 * s)) / s;
		return gr;
	}

	point grad_alpine1(const int &n, const point &x)
	{
		//++call_f;
		point gr(n, 0);
		for (int i = 0; i < n; ++i)
			gr[i] = (0.1 + x[i] * cos(x[i]) + sin(x[i])) * fabs(0.1 * x[i] + x[i] * sin(x[i]));
		return gr;
	}

	point grad_morse(const int &n, const point &x)
	{
		point gr(n,0);
		for (int i = 0; i < n; ++i){
			//++call_f;
			gr[i] = sumMorse(n, x, i);
		}
		return gr;
	}

	point grad_f(const int &n, const point &x, const int &fcode)
	{
		switch (fcode)
		{
		case 0: return grad_sphere(n, x);
		case 1: return grad_rosenbrock(n, x);
		case 2: return grad_gauss(n, x);
		case 3: return grad_rastrigin(n, x);
		case 4: return grad_schwefel(n, x);
		case 5: return grad_griewank(n, x);
		case 6: return grad_ackley(n, x);
		case 7:	return grad_rastrigin_novgorod(n, x);
		case 20: return grad_morse(n, x);
		default: return point(n, 0); 
		}
	}

}
#endif // TESTFUN_H_INCLUDED
