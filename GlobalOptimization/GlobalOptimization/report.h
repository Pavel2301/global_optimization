#ifndef REPORT_H_INCLUDED
#define REPORT_H_INCLUDED

#include "utils.h"

using namespace std;

class Report
{
private:
	int _call_function;                                       //число вызовов ЦФ
	double _value_function;                                   //значение ЦФ
public:
	Report(const int &call_function, const double &value_function) :_call_function(call_function), _value_function(value_function)
	{}
	void print_report_cf()const
	{
		cout << _call_function << ";";
	}
	void print_report()const
	{
		cout << point_to_comma(to_string(_value_function)) << ";";
		//cout << _call_function << ";" << point_to_comma(to_string(_value_function)) << endl;
	}
	int show_call_function()const
	{
		return _call_function;
	}
	double show_value_function()const
	{
		return  _value_function;
	}
};

class FullReport
{
private:
	vector<Report> _report;
public:
	FullReport() {}
	FullReport(const Report &report_item)
	{
		_report.push_back(report_item);
	}
	vector<Report> show_report()const
	{
		return _report;
	}
	int show_report_size()const
	{
		return _report.size();
	}
	int show_call_function(const int &k)const
	{
		if (k < _report.size())
			return _report[k].show_call_function();
		else
			return 0;
	}
	double show_value_function(const int &k)const
	{
		if (k < _report.size())
			return _report[k].show_value_function();
	}
	void print_report(const int &k)const
	{
		if (k < _report.size())
			_report[k].print_report();
	}
	void print_full_report()const
	{
		for (int i = 0; i < _report.size(); i++)
			_report[i].print_report();
		cout << endl;
	}
	void print_full_report_cf()const
	{
		for (int i = 0; i < _report.size(); i++)
			_report[i].print_report_cf();
		cout << endl;
	}
	void insert_into_report(const Report &report_item)
	{
		_report.push_back(report_item);
	}
	void merge_reports_after(const FullReport &after)
	{
		int size = after.show_report_size();
		for (int i = 0; i < size; i++)
			_report.push_back(after.show_report()[i]);
	}
};



#endif // REPORT_H_INCLUDED
