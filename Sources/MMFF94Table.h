#pragma once

#ifndef MMFF94TABLE_H
#define MMFF94TABLE_H

#include "Matrix.h"

// --- MMFF94参数有效值类定义 ---
class BSVal {
public:
	double kb;
	double r0;
	BSVal(double kb_now = NAN, double r0_now = NAN) :kb(kb_now), r0(r0_now) {}
	friend ostream& operator<<(ostream& out, const BSVal& bsv)
	{
		out << bsv.kb << "," << bsv.r0;
		return out;
	}
	int GetPCount()const
	{
		return 2;
	}
};
class ABVal {
public:
	double ka_ijk;
	double theta_0;
	ABVal(double ka_now = NAN, double theta_now = NAN) :ka_ijk(ka_now), theta_0(theta_now) {}
	friend ostream& operator<<(ostream& out, const ABVal& abv)
	{
		out << abv.ka_ijk << "," << abv.theta_0;
		return out;
	}
	int GetPCount()const
	{
		return 2;
	}
};
class SBVal {
public:
	double kba_ijk;
	double kba_kji;

	SBVal(double kba_ijk_now = NAN, double kba_kji_now = NAN) :kba_ijk(kba_ijk_now), kba_kji(kba_kji_now) {}
	friend ostream& operator<<(ostream& out, const SBVal& sbv)
	{
		out << sbv.kba_ijk << "," << sbv.kba_kji;
		return out;
	}
	int GetPCount()const
	{
		return 2;
	}
};
class OPBVal {
public:
	double koop;
	OPBVal(double koop_now = NAN) :koop(NAN) {}
	friend ostream& operator<<(ostream& out, const OPBVal& opbv)
	{
		out << opbv.koop;
		return out;
	}
	int GetPCount()const
	{
		return 1;
	}
};
class TIVal {
public:
	double v1;
	double v2;
	double v3;
	TIVal(double v1_now = NAN, double v2_now = NAN, double v3_now = NAN) :v1(v1_now), v2(v2_now), v3(v3_now) {}
	friend ostream& operator<<(ostream& out, const TIVal& tiv)
	{
		out << tiv.v1 << "," << tiv.v2 << "," << tiv.v3;
		return out;
	}
	int GetPCount()const
	{
		return 3;
	}
};
class VDWVal {
public:
	double alpha_i;
	double n_i;
	double a_i;
	double g_i;
	VDWVal(double alpha_now = NAN, double n_now = NAN, double a_now = NAN, double g_now = NAN) :
		alpha_i(alpha_now), n_i(n_now), a_i(a_now), g_i(g_now) {
	}
	friend ostream& operator<<(ostream& out, const VDWVal& vdwv)
	{
		out << vdwv.alpha_i << "," << vdwv.n_i << "," << vdwv.a_i << vdwv.g_i;
		return out;
	}
	int GetPCount()const
	{
		return 4;
	}
};
class VDWProVal
{
public:
	double rv;
	double e;
	VDWProVal(double rv_ij = NAN, double epsilon_ij = NAN) :rv(rv_ij), e(epsilon_ij) {}
	friend ostream& operator<<(ostream& out, const VDWProVal& vdwv)
	{
		out << vdwv.rv << "," << vdwv.e;
		return out;
	}
	int GetPCount()const
	{
		return 2;
	}
};

// --- MMFF94参数表行类定义 ---
class BSLine
{
private:
	//string bt;
	int mtype_i;
	int mtype_j;
	double kb;
	double r0;
	//string source;
	int mtype_count = 2;
public:
	BSLine(const vector<string>& info, const vector<int>& titlenum)
	{
		//bt = info[titlenum[0]];
		mtype_i = atoi(info[titlenum[1]].c_str());
		mtype_j = atoi(info[titlenum[2]].c_str());
		kb = stod(info[titlenum[3]]);
		r0 = stod(info[titlenum[4]]);
	}

	MTYPE_INDEX GetMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_i, mtype_j };
	}
	MTYPE_INDEX GetRevMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1]) || (mtype_i == mtype[1] && mtype_j == mtype[0]);
	}
	BSVal GetVal()const { return BSVal(kb, r0); }

	int GetMTypeI()const { return mtype_i; }
	int GetMTypeJ()const { return mtype_j; }
	double GetKB()const { return kb; }
	double GetR0()const { return r0; }

	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(BS_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << mtype_j << sep << kb << sep << r0 << endl;
	}
};
class ABLine
{
private:
	//string angle_type;
	int mtype_i;
	int mtype_j;
	int mtype_k;
	double ka_ijk;
	double theta_0;
	//string comment_origin;
	int mtype_count = 3;
public:
	ABLine(const vector<string>& info, const vector<int>& titlenum)
	{
		//angle_type = info[titlenum[0]];
		mtype_i = atoi(info[titlenum[1]].c_str());
		mtype_j = atoi(info[titlenum[2]].c_str());
		mtype_k = atoi(info[titlenum[3]].c_str());
		ka_ijk = stod(info[titlenum[4]]);
		double theta_deg = stod(info[titlenum[5]]);
		theta_0 = DegToRad(theta_deg);
	}
	MTYPE_INDEX GetMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_i, mtype_j, mtype_k };
	}
	MTYPE_INDEX GetRevMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_k, mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0]);
	}
	ABVal GetVal()const { return ABVal(ka_ijk, theta_0); }

	int GetMTypeI()const { return mtype_i; }
	int GetMTypeJ()const { return mtype_j; }
	int GetMTypeK()const { return mtype_k; }
	double GetKA_IJK()const { return ka_ijk; }
	double GetTheta_0()const { return theta_0; }

	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(AB_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << mtype_j << sep << mtype_k << sep << ka_ijk << sep << theta_0 << endl;
	}
};
class SBLine
{
private:
	//string sbt;
	int mtype_i;
	int mtype_j;
	int mtype_k;
	double kba_ijk;
	double kba_kji;
	//string source;
	int mtype_count = 3;
public:
	SBLine(const vector<string>& info, const vector<int>& titlenum)
	{
		//sbt = info[titlenum[0]];
		mtype_i = atoi(info[titlenum[1]].c_str());
		mtype_j = atoi(info[titlenum[2]].c_str());
		mtype_k = atoi(info[titlenum[3]].c_str());
		kba_ijk = stod(info[titlenum[4]]);
		kba_kji = stod(info[titlenum[5]]);
	}
	MTYPE_INDEX GetMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_i, mtype_j, mtype_k };
	}
	MTYPE_INDEX GetRevMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_k, mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0]);
	}
	SBVal GetVal()const { return SBVal(kba_ijk, kba_kji); }

	int GetMTypeI()const { return mtype_i; }
	int GetMTypeJ()const { return mtype_j; }
	int GetMTypeK()const { return mtype_k; }
	double GetKBA_IJK()const { return kba_ijk; }
	double GetKBA_KJI()const { return kba_kji; }


	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(SB_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << mtype_j << sep << mtype_k << sep << kba_ijk << sep << kba_kji << endl;
	}
};
class OPBLine
{
private:
	int mtype_i;
	int mtype_j;
	int mtype_k;
	int mtype_l;
	double koop;
	//string source;
	int mtype_count = 4;
public:
	OPBLine(const vector<string>& info, const vector<int>& titlenum)
	{
		mtype_i = atoi(info[titlenum[0]].c_str());
		mtype_j = atoi(info[titlenum[1]].c_str());
		mtype_k = atoi(info[titlenum[2]].c_str());
		mtype_l = atoi(info[titlenum[3]].c_str());
		koop = stod(info[titlenum[4]]);
	}
	MTYPE_INDEX GetMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_i, mtype_j, mtype_k, mtype_l };
	}
	MTYPE_INDEX GetRevMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_k, mtype_j, mtype_i, mtype_l };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2] && mtype_l == mtype[3]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0] && mtype_l == mtype[3]);
	}
	OPBVal GetVal()const { return OPBVal(koop); }

	int GetMTypeI()const { return mtype_i; }
	int GetMTypeJ()const { return mtype_j; }
	int GetMTypeK()const { return mtype_k; }
	int GetMTypeL()const { return mtype_l; }
	double GetKoop()const { return koop; }


	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(OPB_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << mtype_j << sep << mtype_k << sep << mtype_l << sep << koop << endl;
	}
};
class TILine
{
private:
	//string tt;
	int mtype_i;
	int mtype_j;
	int mtype_k;
	int mtype_l;
	double v1;
	double v2;
	double v3;
	//string source;
	int mtype_count = 4;
public:
	TILine(const vector<string>& info, const vector<int>& titlenum)
	{
		//tt = info[titlenum[0]];
		mtype_i = atoi(info[titlenum[1]].c_str());
		mtype_j = atoi(info[titlenum[2]].c_str());
		mtype_k = atoi(info[titlenum[3]].c_str());
		mtype_l = atoi(info[titlenum[4]].c_str());
		v1 = stod(info[titlenum[5]]);
		v2 = stod(info[titlenum[6]]);
		v3 = stod(info[titlenum[7]]);
	}
	MTYPE_INDEX GetMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_i, mtype_j, mtype_k, mtype_l };
	}
	MTYPE_INDEX GetRevMTypeVec()const
	{
		return MTYPE_INDEX{ mtype_l, mtype_k, mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2] && mtype_l == mtype[3]) ||
			(mtype_i == mtype[3] && mtype_j == mtype[2] && mtype_k == mtype[1] && mtype_l == mtype[0]);
	}
	int GetMTypeI()const { return mtype_i; }
	int GetMTypeJ()const { return mtype_j; }
	int GetMTypeK()const { return mtype_k; }
	int GetMTypeL()const { return mtype_l; }
	double GetV1()const { return v1; }
	double GetV2()const { return v2; }
	double GetV3()const { return v3; }
	TIVal GetVal()const { return TIVal(v1, v2, v3); }

	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(TI_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << mtype_j << sep << mtype_k << sep << mtype_l << sep << v1 << sep << v2 << sep << v3 << endl;
	}
};
class VDWLine
{
private:
	int mtype_i;
	double alpha_i;
	double n_i;
	double a_i;
	double g_i;
	//string da;
	//string sym;
	//string origin;
	int mtype_count = 1;
public:
	VDWLine(const vector<string>& info, const vector<int>& titlenum)
	{
		mtype_i = atoi(info[titlenum[0]].c_str());
		alpha_i = stod(info[titlenum[1]]);
		n_i = stod(info[titlenum[2]]);
		a_i = stod(info[titlenum[3]]);
		g_i = stod(info[titlenum[4]]);
	}
	MTYPE_INDEX GetMTypeVec()const { return MTYPE_INDEX{ mtype_i }; }
	MTYPE_INDEX GetRevMTypeVec()const { return MTYPE_INDEX{ mtype_i }; }
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype[0] == mtype_i);
	}
	int GetMTypeI()const { return mtype_i; }
	double GetAlphaI()const { return alpha_i; }
	double GetNI()const { return n_i; }
	double GetAI()const { return a_i; }
	double GetGI()const { return g_i; }
	VDWVal GetVal()const
	{
		return VDWVal(alpha_i, n_i, a_i, g_i);
	}

	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(VDW_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << alpha_i << sep << n_i << sep << a_i << sep << g_i << endl;
	}

};

inline string GetSpheTbPath(int num)
{
	string filename = "M" + to_string(num) + ".csv";
	return SP_FOLDER + filename;
}

class MSymLine
{
private:
	string msym;
	int mtype;
public:
	MSymLine(const vector<string>& info, const vector<int>& titlenum)
	{
		msym = info[titlenum[0]];
		mtype = atoi(info[titlenum[1]].c_str());
	}
	string GetIndex() const { return msym; }
	int GetVal()const { return mtype; }
	void Print(string sep = "\t")const
	{
		cout << msym << sep << mtype << endl;
	}
};
class SP_SpLine
{
private:
	int seq;
	RadVec3 rvec;
public:
	SP_SpLine(const vector<string>& info, const vector<int>& titlenum)
	{
		seq = atoi(info[titlenum[0]].c_str());
		double r = stod(info[titlenum[1]]);
		double theta = DegToRad(stod(info[titlenum[2]]));
		double varphi = DegToRad(stod(info[titlenum[3]]));
		rvec = RadVec3(r, theta, varphi);
	}
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(SP_SPTB_TITLE, sep);
			cout << endl;
		}
		cout << seq << sep;
		rvec.Print(sep);
	}
	int GetSeq()const { return seq; }
	RadVec3 GetRVec() { return rvec; }
	void SpinTheta(double theta_rad)
	{
		rvec.theta += theta_rad;
		rvec.theta = fmod(rvec.theta, PI_DOUBLE);
	}
	void SpinVarphi(double varphi_rad)
	{
		rvec.varphi += varphi_rad;
		rvec.varphi = fmod(rvec.varphi, PI_DOUBLE);
	}
};

class SP_ReLine
{
private:
	int seq;
	Vec3 xyz;
public:
	SP_ReLine() :seq(-1), xyz() {}
	SP_ReLine(SP_SpLine line) :seq(line.GetSeq()), xyz(line.GetRVec()) {}
	SP_ReLine(const SP_ReLine& r) : seq(r.seq), xyz(r.xyz) {}
	SP_ReLine(int seq, Vec3 v) : seq(seq), xyz(v) {}
	int GetSeq() const { return seq; }
	Vec3 GetVec()const { return xyz; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(SP_RETB_TITLE, sep);
			cout << endl;
		}
		cout << seq << sep;
		xyz.Print(sep);
	}

};

class SP_SpTable {
private:
	vector<vector<SP_SpLine>> all_tb;
public:
	SP_SpTable()
	{
		vector<vector<SP_SpLine>> at;
		for (int i = 0; i <= MAX_SP_SIZE; i++)
		{
			if (i < 1) {
				at.push_back(vector<SP_SpLine>());
				continue;
			}
			vector<SP_SpLine> sp_sptb;
			ReadTableByTitle(GetSpheTbPath(i), SP_SPTB_TITLE, sp_sptb);
			at.push_back(sp_sptb);
		}
		all_tb = at;
	}
	SP_SpLine GetSPLine(int m, int row)const
	{
		return all_tb[m][row];
	}
	const vector<SP_SpLine>& operator[](int m)const
	{
		return all_tb[m];
	}
	void Print(string sep = "\t")
	{
		for (int i = 2; i <= MAX_SP_SIZE; i++)
		{
			string title = "M" + to_string(i) + ".csv";
			PrintCmdSepTitle(title);
			PrintTableTitle(SP_SPTB_TITLE);
			PrintSpecialVector(all_tb[i], sep);
		}
	}
};



struct EnergyFundTable
{
	vector<BSLine> bs_tb;
	vector<ABLine> ab_tb;
	vector<SBLine> sb_tb;
	vector<OPBLine> opb_tb;
	vector<TILine> ti_tb;
	vector<VDWLine> vdw_tb;
	MSYM_MAP msym_map;
};

EnergyFundTable ReadEnergySolidParam(bool print_yes = false, int max_row_count = -1);







#endif // DEBUG

