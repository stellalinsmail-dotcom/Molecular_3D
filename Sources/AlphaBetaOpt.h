#pragma once
#ifndef ALPHABETAAOPT_H
#define ALPHABETAAOPT_H

#include "types.h"
//#include "Atom.h"
//#include "Molecule.h"
//#include "Matrix.h"
#include "MMFF94Table.h"
#include "MMFF94Cal.h"

struct OptRecVal {
	double min_energy;
	double opt_time;
};

class OptRecLine {
private:
	string smiles;
	double min_energy;
	double opt_time;
public:
	OptRecLine(string s, double e, double t) :smiles(s), min_energy(e), opt_time(t) {}
	OptRecLine(string s, OptRecVal v) :smiles(s), min_energy(v.min_energy), opt_time(v.opt_time) {}
	OptRecLine(const vector<string>& info, const vector<int>& titlenum)
	{
		smiles = info[titlenum[0]];
		min_energy = atof(info[titlenum[1]].c_str());
		if (info[titlenum[2]] == "INF") opt_time = INFINITY;
		else opt_time = atof(info[titlenum[2]].c_str());
	}
	OptRecLine(const OptRecLine& r) :smiles(r.smiles), min_energy(r.min_energy), opt_time(r.opt_time) {}
	string GetIndex()const { return smiles; }
	OptRecVal GetVal()const { return { min_energy ,opt_time }; }
	double GetOptTime()const { return opt_time; }
	void Print(string sep = "\t")const
	{
		cout << smiles << sep << min_energy << sep << opt_time << endl;
	}

	friend ostream& operator<<(ostream& os, const OptRecLine& orc)
	{
		os << orc.smiles << "," << orc.min_energy << ",";
		if (orc.opt_time >= MAX_REC_TIME) os << "INF";
		else os << orc.opt_time;
		return os;
	}
};


//--- 优化相关函数 ---
template<typename T>
bool JudgeStop(const vector<T>& new_tb, const vector<T>& old_tb, const vector<T>& sep_tb)
{
	int ssize = sep_tb.size();
	for (int i = 0; i < ssize; i++)
	{
		double a = abs(old_tb[i] - new_tb[i]);
		if (a > sep_tb[i])
		{
			return false;
		}
	}
	return true;
}
template<typename T>
bool JudgeStop2(const vector<vector<T>>& new_tb, const vector<vector<T>>& old_tb, const vector<vector<T>>& sep_tb)
{
	int ssize = sep_tb.size();
	for (int i = 0; i < ssize; i++)
	{
		for (int j = 0; j < sep_tb[i].size(); j++)
		{
			double a = abs(old_tb[i][j] - new_tb[i][j]);
			if (a > sep_tb[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

vector<int> GenerateRandSeq(int size);
void CheckJsonExist(vector<OptRecLine>& optrec_tb);
V1_DTB AlphaOpt(bool has_circle, const V2_DTB& beta_tb, OptRecVal& alpha_result, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp, bool detail_print = false, bool savedata_yes = false, string now_smiles = "");
V2_DTB BetaOpt(bool has_circle, const V1_DTB& alpha_tb, OptRecVal& beta_result, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp, bool detail_print = false);

#endif // ALPHABETAAOPT_H
