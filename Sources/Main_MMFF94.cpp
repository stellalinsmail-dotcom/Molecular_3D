#include "types.h"
#include "Atom.h"
#include "Molecule.h"
#include "Matrix.h"
//#include <map>

using namespace std;

// --- 文件路径处理 ---
string GetMMFFPath(string filename)
{
	return MMFF_CSV_FOLDER + filename;
}
string GetSpheTbPath(int num)
{
	string filename = "M" + to_string(num) + ".csv";
	return SP_FOLDER + filename;
}
string GetMNodeTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/MNodeTable.csv";
}
string GetAdjTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/BondTable.csv";
}
string GetXYZTbPath(string can_smiles, int n)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/XYZ_" + to_string(n) + ".csv";
}

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

// --- MMFF94参数表快表类定义 ---
template<typename TLine, typename TVal>
class FastMatrix1 {
private:
	vector<TVal> data;
public:
	FastMatrix1(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<TVal>(dsize, TVal());
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			if (i < 0 || i > max_mtype) continue;
			int ti = htb[i];
			if (ti < 0) continue;
			data[ti] = line.GetVal();
		}
	}
	TVal operator[](const FAST_TABLE_INDEX& m)const
	{
		return data[m[0]];
	}
	int GetSize()const
	{
		return data.size();
	}
};

template<typename TLine, typename TVal>
class FastMatrix2 {
private:
	vector<vector<TVal>> data;
public:
	FastMatrix2(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal()));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			if (min(i, j) < 0 || max(i, j) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			if (min(ti, tj) < 0) continue;
			data[ti][tj] = line.GetVal();
			data[tj][ti] = line.GetVal();
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]];
	}
};

template<typename TLine, typename TVal>
class FastMatrix3 {
private:
	vector<vector<vector<TVal>>> data;
public:
	FastMatrix3(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		//cout <<"Max MType: " << max_mtype << endl;
		data = vector<vector<vector<TVal>>>(dsize, vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal())));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			int k = line.GetMTypeK();
			if (Min3(i, j, k) < 0 || Max3(i, j, k) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			int tk = htb[k];
			if (Min3(ti, tj, tk) < 0) continue;
			data[ti][tj][tk] = line.GetVal();
			data[tk][tj][ti] = data[ti][tj][tk];
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]][m[2]];
	}
};

template<typename TLine, typename TVal>
class FastMatrix4 {
private:
	vector<vector<vector<vector<TVal>>>> data;
public:
	FastMatrix4(const vector<TLine> table, vector<int> htb, bool is_ij_solid = NO)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<vector<vector<vector<TVal>>>>(dsize, vector<vector<vector<TVal>>>(dsize, vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal()))));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			int k = line.GetMTypeK();
			int l = line.GetMTypeL();
			if (Min4(i, j, k, l) < 0 || Max4(i, j, k, l) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			int tk = htb[k];
			int tl = htb[l];
			if (Min4(ti, tj, tk, tl) < 0) continue;
			if (is_ij_solid)
			{
				data[ti][tj][tk][tl] = line.GetVal();
				data[tk][tj][ti][tl] = data[ti][tj][tk][tl];
			}
			else
			{
				data[ti][tj][tk][tl] = line.GetVal();
				data[tl][tk][tj][ti] = line.GetVal();
			}
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]][m[2]][m[3]];
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
			if (i < 2) {
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

class NeedCal {
public:
	bool eb;
	bool ea;
	bool eba;
	bool eoop;
	bool et;
	bool evdw;
	NeedCal(bool is_eb = true, bool is_ea = true, bool is_eba = true, bool is_eoop = true, bool is_et = true, bool is_evdw = true) :
		eb(is_eb), ea(is_ea), eba(is_eba), et(is_et), eoop(is_eoop), evdw(is_evdw) {
	}
};

struct EnergySolidParam
{
	F_BS bs_fast_tb;
	F_AB ab_fast_tb;
	F_SB sb_fast_tb;
	F_OPB opb_fast_tb;
	F_TI ti_fast_tb;
	F_VDW vdw_fast_tb;
	HASH_TB pro_htb;

	BFS2_TB bfs_tb;
	RE_TB re_tb;

	ADJ_LIST short_adj_list;
};
struct EnergyVaryParam
{
	R_TB r_tb;
	R_TB dr_tb;
	VAR_TB dvar_tb;
	PHI_TB phi_tb;
	CHI_TB chi_tb;
};

// --- MMFF94参数表字典类快表定义 ---
//
// --- 表格处理 ---
// template<typename T>class BSMap {
//	PARAM_MAP r0;
//	PARAM_MAP kb;
//public:
//	BSMap(const vector<BSLine>& table) :r0(CompressTableToMap(table, &BSLine::GetR0)),
//		kb(CompressTableToMap(table, &BSLine::GetKB)) {
//	}
//	double GetR0(const vector<int>& key) { return r0[key]; }
//	double GetKB(const vector<int>& key) { return kb[key]; }
//};
//class ABMap {
//private:
//	PARAM_MAP theta_0;
//	PARAM_MAP ka_ijk;
//public:
//	ABMap(const vector<ABLine>& table) :theta_0(CompressTableToMap(table, &ABLine::GetTheta_0)),
//		ka_ijk(CompressTableToMap(table, &ABLine::GetKA_IJK)) {
//	}
//	double GetTheta_0(const vector<int>& key) { return theta_0[key]; }
//	double GetKA_IJK(const vector<int>& key) { return ka_ijk[key]; }
//};
//class SBMap {
//private:
//	PARAM_MAP kba_ijk;
//	PARAM_MAP kba_kji;
//public:
//	SBMap(const vector<SBLine>& table) :kba_ijk(CompressTableToMap(table, &SBLine::GetKBA_IJK)),
//		kba_kji(CompressTableToMap(table, &SBLine::GetKBA_KJI)) {
//	}
//	double GetKBA_IJK(const vector<int>& key) { return kba_ijk[key]; }
//	double GetKBA_KJI(const vector<int>& key) { return kba_kji[key]; }
//};
//class OPBMap {
//private:
//	PARAM_MAP koop;
//public:
//	OPBMap(const vector<OPBLine>& table) :koop(CompressTableToMap(table, &OPBLine::GetKoop)) {}
//	double GetKoop(const vector<int>& key) { return koop[key]; }
//};
//class TIMap {
//private:
//	PARAM_MAP v1;
//	PARAM_MAP v2;
//	PARAM_MAP v3;
//public:
//	TIMap(const vector<TILine>& table) :v1(CompressTableToMap(table, &TILine::GetV1)),
//		v2(CompressTableToMap(table, &TILine::GetV2)),
//		v3(CompressTableToMap(table, &TILine::GetV3)) {
//	}
//	double GetV1(const vector<int>& key) { return v1[key]; }
//	double GetV2(const vector<int>& key) { return v2[key]; }
//	double GetV3(const vector<int>& key) { return v3[key]; }
//};
//class VDWMap {
//private:
//	PARAM_MAP alpha_i;
//	PARAM_MAP n_i;
//	PARAM_MAP a_i;
//	PARAM_MAP g_i;
//public:
//	VDWMap(const vector<VDWLine>& table) :alpha_i(CompressTableToMap(table, &VDWLine::GetAlphaI)),
//		n_i(CompressTableToMap(table, &VDWLine::GetNI)),
//		a_i(CompressTableToMap(table, &VDWLine::GetAI)),
//		g_i(CompressTableToMap(table, &VDWLine::GetGI)) {
//	}
//	double GetAlphaI(const vector<int>& key) { return alpha_i[key]; }
//	double GetNI(const vector<int>& key){ return n_i[key]; }
//	double GetAI(const vector<int>& key) { return a_i[key]; }
//	double GetGI(const vector<int>& key) { return g_i[key]; }
//};

// --- 辅助表格生成 ---

vector<int> GetHashTable(const MTYPE_SET& mtype_set)
{
	int max_mtype = *mtype_set.rbegin();
	vector<int> htb(max_mtype + 1, -1);
	int index = 0;
	for (auto mtype : mtype_set)
	{
		htb[mtype] = index;
		index++;
	}
	return htb;
}

vector<int> GetProHashTable(const vector<int>& mtype_tb, const vector<int>& htb)
{
	int tbsize = mtype_tb.size();
	vector<int> pro_htb(tbsize, -1);
	for (int i = 0; i < tbsize; i++)
	{
		pro_htb[i] = htb[mtype_tb[i]];
	}
	return pro_htb;
}

vector<double> GetAlphaSepTable(const ADJ_LIST short_adj_list)
{
	int nhc = short_adj_list.size();
	vector<double> a(nhc, NAN);
	for (int i = 0; i < nhc; i++)
	{
		a[i] = PI_HALF / (short_adj_list[i].GetBonds().size() - 1);
	}
	return a;

}

// --- 简单参数先导计算 ---

inline double xyz2r(const Vec3& v1, const Vec3& v2)
{
	return Vec3(v2 - v1).GetLen();
}
inline double xyz2vartheta_rad(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
	Vec3 ab = v2 - v1;
	Vec3 cb = v2 - v3;
	return acos(Cdot(ab, cb) / (ab.GetLen() * cb.GetLen()));
}
inline double xyz2phi_rad(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& v4)
{
	Vec3 a = v1 - v2;
	Vec3 b = v3 - v2;
	Vec3 c = v4 - v3;

	Vec3 n_ab = Cross(a, b);
	Vec3 n_bc = Cross(b, c);
	double phi_rad = acos(Cdot(n_ab, n_bc) / (n_ab.GetLen() * n_bc.GetLen()));

	//cout << endl;
	//n_ab.Print();
	//n_bc.Print();
	//cout << "val: " << Cdot(n_ab, n_bc) / (n_ab.GetLen() * n_bc.GetLen()) << endl;
	//cout << "phi_1: " << phi_rad << endl;

	double k = -Cdot(a, b) / Cdot(b, b);
	Vec3 d = a + k * b;
	bool judge = (Cdot(d, c) < 0);
	if (judge) phi_rad = PI - phi_rad;
	return phi_rad;
}
inline double xyz2chi_rad(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& v4)
{
	Vec3 c = v4 - v2;
	Vec3 n = Cross(v1 - v2, v3 - v2);
	double chi_rad = PI_HALF - acos(Cdot(c, n) / c.GetLen() / n.GetLen());
	return chi_rad;
}

// --- 参数表格先导计算 ---
VDWProVal PreCalRE(VDWVal vdw_i, VDWVal vdw_j)
{
	double alpha_i = vdw_i.alpha_i, alpha_j = vdw_j.alpha_i;
	double n_i = vdw_i.n_i, n_j = vdw_j.n_i;
	double a_i = vdw_i.a_i, a_j = vdw_j.a_i;
	double g_i = vdw_i.g_i, g_j = vdw_j.g_i;

	double B = 0.2;
	double beta = 12;
	double rv_ii = a_i * pow(alpha_i, 0.25);
	double rv_jj = a_j * pow(alpha_j, 0.25);
	double gamma_ij = (rv_ii - rv_jj) / (rv_ii + rv_jj);
	double rv_ij = 0.5 * (rv_ii + rv_jj) * (1 + B * (1 - exp(-beta * pow(gamma_ij, 2))));
	double epsilon_ij = 181.16 * g_i * g_j * alpha_i * alpha_j / (pow(alpha_i / n_i, 0.5) + pow(alpha_j / n_j, 0.5)) / pow(rv_ij, 6);
	return VDWProVal(rv_ij, epsilon_ij);
}
void PreCalRETb(RE_TB& re_tb, MTYPE_SET mtype_set, const F_VDW& vdw_fast_tb, const HASH_TB& pro_htb)
{
	int tsize = vdw_fast_tb.GetSize();
	int ac = pro_htb.size();
	re_tb.clear();
	re_tb = RE_TB(ac, vector<VDWProVal>(ac, VDWProVal()));
	vector<vector<VDWProVal>> re_m_tb(tsize, vector<VDWProVal>(tsize, VDWProVal()));
	for (int i = 0; i < tsize; i++)
	{
		for (int j = i; j < tsize; j++)
		{
			//cout << "***** fasttable_i,j: " << i << " " << j << endl;
			re_m_tb[i][j] = PreCalRE(vdw_fast_tb[{i}], vdw_fast_tb[{j}]);
			re_m_tb[j][i] = re_m_tb[i][j];
		}
	}
	for (int i = 0; i < ac; i++)
	{
		for (int j = i + 1; j < ac; j++)
		{
			re_tb[i][j] = re_m_tb[pro_htb[i]][pro_htb[j]];
			re_tb[j][i] = re_tb[i][j];
		}
	}
}
void PreCalDRTb(R_TB& r_tb, R_TB& dr_tb, const XYZ_TB& xyz_tb, ADJ_LIST short_adj_list, const F_BS& bs_fast_tb, const HASH_TB& pro_htb)
{
	int tsize = xyz_tb.size();
	int nhc = short_adj_list.size();
	r_tb.clear();
	dr_tb.clear();
	r_tb = R_TB(tsize, vector<double>(tsize, NAN));
	dr_tb = R_TB(tsize, vector<double>(tsize, NAN));
	for (int j = 0; j < tsize; j++)
	{
		for (int i = j + 1; i < tsize; i++)
		{
			r_tb[i][j] = xyz2r(xyz_tb[i], xyz_tb[j]);
			r_tb[j][i] = r_tb[i][j];
			//cout << i << "\t" << j<< "\tr_fcal: " << r_tb[i][j] << endl;
		}
		//PrintCommonVector(r_tb[j]);
		//cout << endl;
	}
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		for (int i = 0; i < nsize; i++)
		{
			int child = nb_node[i].GetDesSeq();
			if (child > parent)
			{
				FAST_TABLE_INDEX fti = { pro_htb[parent],pro_htb[child] };
				double r0 = bs_fast_tb[fti].r0;

				dr_tb[parent][child] = r_tb[parent][child] - r0;
				dr_tb[child][parent] = dr_tb[parent][child];
				//cout <<endl<< "r0: " <<r0<< endl;
				//cout << parent << "\t" << child << "\tr_cal: " << r_tb[parent][child] << endl;
				//cout <<parent<<"\t" << i << "\tdr: " << dr_tb[parent][i] << endl;
			}
		}
		//PrintCommonVector(dr_tb[j]);
		//cout << endl;
	}
}
void PreCalDVarTb(VAR_TB& dvar_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list, const F_AB& ab_fast_tb, const HASH_TB& pro_htb)
{
	int tsize = xyz_tb.size();
	int nhc = short_adj_list.size();

	dvar_tb.clear();
	dvar_tb = VAR_TB(nhc, vector<vector<double>>(MAX_ADJ_NODE_SIZE, vector<double>(MAX_ADJ_NODE_SIZE, NAN)));

	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		for (int i = 0; i < nsize; i++)
		{
			int child_i = nb_node[i].GetDesSeq();
			for (int k = i + 1; k < nsize; k++)
			{
				int child_k = nb_node[k].GetDesSeq();
				FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent],pro_htb[child_k] };
				double var = xyz2vartheta_rad(xyz_tb[child_i], xyz_tb[parent], xyz_tb[child_k]);
				double var0 = ab_fast_tb[fti].theta_0;
				dvar_tb[parent][i][k] = var - var0;
				dvar_tb[parent][k][i] = dvar_tb[parent][i][k];
				//cout<<parent <<" " << child_i << " " << child_k <<" " << var << "\t" << var0 << endl;
			}
			//PrintCommonVector(dvar_tb[parent][i]);
			//cout << endl;
		}
		//cout << endl;
	}
}
void PreCalPhiTb(PHI_TB& phi_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list)
{
	int nhc = short_adj_list.size();
	double sum_et = 0;
	phi_tb = PHI_TB(nhc, V3_DTB(nhc, V2_DTB(MAX_ADJ_NODE_SIZE, V1_DTB(MAX_ADJ_NODE_SIZE, NAN))));
	for (int j = 0; j < nhc; j++)
	{
		int parent_j = j;
		vector<PointTo> nb_node_j = short_adj_list[j].GetBonds();
		int jnsize = nb_node_j.size();
		for (int k = 0; k < jnsize; k++)
		{
			int parent_k = nb_node_j[k].GetDesSeq();
			if (parent_k > parent_j && parent_k < nhc)
			{
				vector<PointTo> nb_node_k = short_adj_list[parent_k].GetBonds();
				int knsize = nb_node_k.size();
				for (int i = 0; i < jnsize; i++)
				{
					int child_i = nb_node_j[i].GetDesSeq();
					if (child_i != parent_k)
					{
						for (int l = 0; l < knsize; l++)
						{
							int child_l = nb_node_k[l].GetDesSeq();
							if (child_l != parent_j)
							{
								//cout << endl;
								//cout << child_i << "\t" << parent_j << "\t" << parent_k << "\t" << child_l << endl;
								//FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent_j],pro_htb[parent_k],pro_htb[child_l] };
								phi_tb[parent_j][parent_k][i][l] = xyz2phi_rad(xyz_tb[child_i], xyz_tb[parent_j], xyz_tb[parent_k], xyz_tb[child_l]);

							}
						}
					}
				}
			}
		}
	}
}
void PreCalChiTb(CHI_TB& chi_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list)
{
	int nhc = short_adj_list.size();
	chi_tb.clear();
	chi_tb = PHI_TB(nhc, V3_DTB(MAX_ADJ_NODE_SIZE, V2_DTB(MAX_ADJ_NODE_SIZE, V1_DTB(MAX_ADJ_NODE_SIZE, NAN))));
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		if (nsize == 3)
		{
			for (int l = 0; l < 3; l++)
			{
				int child_l = nb_node[l].GetDesSeq();
				int i = (l + 1) % 3;
				int k = (i + 1) % 3;
				int child_i = nb_node[i].GetDesSeq();
				int child_k = nb_node[k].GetDesSeq();

				chi_tb[parent][l][i][k] = xyz2chi_rad(xyz_tb[child_i], xyz_tb[parent], xyz_tb[child_k], xyz_tb[child_l]);
			}
		}
	}
}
BFS2_TB Bfs2(int ac, const ADJ_LIST& short_adj_list)
{
	int nhc = short_adj_list.size();
	BFS2_TB all_rec(nhc, vector<int>());
	for (int j = 0; j < nhc; j++)
	{
		vector<bool> rec(ac, false);
		rec[j] = true;

		vector<PointTo> nb_node_a = short_adj_list[j].GetBonds();
		int nasize = nb_node_a.size();
		for (int na = 0; na < nasize; na++)
		{
			int child_a = nb_node_a[na].GetDesSeq();
			rec[child_a] = true;
			if (child_a < nhc)
			{
				vector<PointTo> nb_node_b = short_adj_list[child_a].GetBonds();
				int nbsize = nb_node_b.size();
				for (int nb = 0; nb < nbsize; nb++)
				{
					int child_b = nb_node_b[nb].GetDesSeq();
					rec[child_b] = true;
				}
			}
		}
		for (int i = j + 1; i < ac; i++)
		{
			if (!rec[i])
			{
				all_rec[j].push_back(i);
			}
		}
	}
	return all_rec;
}


// 注释：除R表外，各表大小为 (ac|nhc)*6*6*...

// --- MMFF94能量计算函数定义 ---
inline double GetEB(double dr_ij, double kb)
{
	return 143.9525 * kb / 2 * pow(dr_ij, 2) * (1 - 2 * dr_ij + 2.33333333333 * pow(dr_ij, 2));
}
inline double GetEA(double dvar_rad, double ka)
{
	return 0.043844 * ka / 2 * pow(dvar_rad, 2) * (1 + -0.000122 * dvar_rad);
}
inline double GetEBA(double dr_ij, double dr_kj, double dvar_rad, double kba_ijk, double kba_kji)
{
	return 2.51210 * (kba_ijk * dr_ij + kba_kji * dr_kj) * dvar_rad;
}
inline double GetET(double phi_rad, double v1, double v2, double v3)
{
	return 0.5 * (v1 * (1 + cos(phi_rad)) + v2 * (1 - cos(2 * phi_rad)) + v3 * (1 + cos(3 * phi_rad)));
}
inline double GetEOOP(double chi_ijkl_rad, double koop)
{
	return 0.034844 * koop / 2 * pow(chi_ijkl_rad, 2);
}
inline double GetEVDW(double r_ij, double rv_ij, double epsilon_ij)
{
	//cout << "result1: " << 1.07 * rv_ij / (r_ij + 0.07 * rv_ij) << endl;
	//cout << "result2: " << 1.12 * pow(rv_ij, 7) / (pow(r_ij, 7) + 0.12 * pow(rv_ij, 7)) << endl;
	return epsilon_ij * pow(1.07 * rv_ij / (r_ij + 0.07 * rv_ij), 7) * (1.12 * pow(rv_ij, 7) / (pow(r_ij, 7) + 0.12 * pow(rv_ij, 7)) - 2);
}

// --- MMFF94能量分步骤计算 ---
double CalEB(const ADJ_LIST& short_adj_list, const R_TB& dr_tb, const F_BS& bs_fast_tb, const HASH_TB& pro_htb)
{
	int nhc = short_adj_list.size();
	double sum_eb = 0;
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		for (int i = 0; i < nsize; i++)
		{
			int child = nb_node[i].GetDesSeq();
			if (child > parent)
			{
				FAST_TABLE_INDEX fti = { pro_htb[parent],pro_htb[child] };

				sum_eb += GetEB(dr_tb[parent][child], bs_fast_tb[fti].kb);
				//cout << "EB: " << GetEB(dr_tb[parent][i], bs_fast_tb[fti].kb) << endl;
			}
		}
	}
	return sum_eb;
}
double CalEA(const ADJ_LIST& short_adj_list, const VAR_TB& dvar_tb, const F_AB& ab_fast_tb, const HASH_TB& pro_htb)
{
	int nhc = short_adj_list.size();
	double sum_ea = 0;
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		for (int i = 0; i < nsize; i++)
		{
			int child_i = nb_node[i].GetDesSeq();
			for (int k = i + 1; k < nsize; k++)
			{
				int child_k = nb_node[k].GetDesSeq();
				FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent],pro_htb[child_k] };
				double ka_ijk = ab_fast_tb[fti].ka_ijk;

				//cout << dvar_tb[parent][i][k] << endl;
				sum_ea += GetEA(dvar_tb[parent][i][k], ka_ijk);
				//cout << "EA: " << GetEA(dvar_tb[parent][i][k], ka_ijk) << endl;
			}
		}
	}
	return sum_ea;
}
double CalEBA(const ADJ_LIST& short_adj_list, const R_TB& dr_tb, const VAR_TB& dvar_tb, const F_SB& sb_fast_tb, const HASH_TB& pro_htb)
{
	int nhc = short_adj_list.size();
	double sum_eba = 0;
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		for (int i = 0; i < nsize; i++)
		{
			int child_i = nb_node[i].GetDesSeq();
			for (int k = i + 1; k < nsize; k++)
			{
				int child_k = nb_node[k].GetDesSeq();
				FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent],pro_htb[child_k] };
				double kba_ijk = sb_fast_tb[fti].kba_ijk;
				double kba_kji = sb_fast_tb[fti].kba_kji;
				sum_eba += GetEBA(dr_tb[parent][child_i], dr_tb[parent][child_k], dvar_tb[parent][i][k], kba_ijk, kba_kji);

				//cout << endl << parent << "\t" << child_i << "\t" << child_k << endl;
				//cout << dr_tb[parent][child_i] << "\t" << dr_tb[parent][child_k] << dvar_tb[parent][i][k] << endl;
				//cout << "EBA: " << GetEBA(dr_tb[parent][i], dr_tb[parent][k], dvar_tb[parent][i][k], kba_ijk, kba_kji) << endl;
			}
			//cout << endl;
		}
		//cout << endl;
	}
	return sum_eba;
}
double CalEOOP(const ADJ_LIST& short_adj_list, const CHI_TB& chi_tb, const F_OPB& opb_fast_tb, const HASH_TB& pro_htb)
{
	int nhc = short_adj_list.size();
	double sum_eoop = 0;
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		if (nsize == 3)
		{
			for (int l = 0; l < 3; l++)
			{
				int child_l = nb_node[l].GetDesSeq();
				int i = (l + 1) % 3;
				int k = (i + 1) % 3;
				int child_i = nb_node[i].GetDesSeq();
				int child_k = nb_node[k].GetDesSeq();

				FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent],pro_htb[child_k],pro_htb[child_l] };
				double koop = opb_fast_tb[fti].koop;
				//double chi_ijkl = xyz2chi_rad(xyz_tb[child_i], xyz_tb[parent], xyz_tb[child_k], xyz_tb[child_l]);
				double chi_ijkl = chi_tb[parent][l][i][k];
				sum_eoop += GetEOOP(chi_ijkl, koop);
			}
		}
	}
	return sum_eoop;
}
double CalET(const ADJ_LIST& short_adj_list, const PHI_TB& phi_tb, const F_TI& ti_fast_tb, const HASH_TB& pro_htb)
{
	int nhc = short_adj_list.size();
	double sum_et = 0;
	for (int j = 0; j < nhc; j++)
	{
		int parent_j = j;
		vector<PointTo> nb_node_j = short_adj_list[j].GetBonds();
		int jnsize = nb_node_j.size();
		for (int k = 0; k < jnsize; k++)
		{
			int parent_k = nb_node_j[k].GetDesSeq();
			if (parent_k > parent_j && parent_k < nhc)
			{
				vector<PointTo> nb_node_k = short_adj_list[parent_k].GetBonds();
				int knsize = nb_node_k.size();
				for (int i = 0; i < jnsize; i++)
				{
					int child_i = nb_node_j[i].GetDesSeq();
					if (child_i != parent_k)
					{
						for (int l = 0; l < knsize; l++)
						{
							int child_l = nb_node_k[l].GetDesSeq();
							if (child_l != parent_j)
							{
								FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent_j],pro_htb[parent_k],pro_htb[child_l] };
								double v1 = ti_fast_tb[fti].v1;
								double v2 = ti_fast_tb[fti].v2;
								double v3 = ti_fast_tb[fti].v3;
								double phi_rad = phi_tb[parent_j][parent_k][i][l];
								sum_et += GetET(phi_rad, v1, v2, v3);

								//cout << endl;
								//cout << child_i << "\t" << parent_j << "\t" << parent_k << "\t" << child_l << endl;
								//cout << v1 << "\t" << v2 << "\t" << v3 << "\t" << phi_rad << endl;
								//cout << "ET: " << GetET(phi_rad, v1, v2, v3) << endl;
							}
						}
					}
				}
			}
		}
	}
	return sum_et;
}
double CalEVDW(const ADJ_LIST& short_adj_list, const BFS2_TB& bfs_tb, const R_TB& r_tb, const RE_TB& re_tb)
{
	double sum_evdw = 0;
	int nhc = bfs_tb.size();
	for (int j = 0; j < nhc; j++)
	{
		//cout << "BFS: " << j << endl;
		//PrintCommonVector(bfs_tb[j]);
		//cout << endl;

		int bnsize = bfs_tb[j].size();
		for (int i = 0; i < bnsize; i++)
		{
			int ci = bfs_tb[j][i];
			sum_evdw += GetEVDW(r_tb[j][ci], re_tb[j][ci].rv, re_tb[j][ci].e);

			//cout << endl;
			//cout << j << "\t" << ci << endl;
			//cout << r_tb[j][ci] << "\t" << re_tb[j][ci].rv << "\t" << re_tb[j][ci].e << endl;
			//cout << "EVDW: " << GetEVDW(r_tb[j][ci], re_tb[j][ci].rv, re_tb[j][ci].e) << endl;
		}

	}
	return sum_evdw;
}

double CalSumEnergy(bool print_yes, NeedCal need_cal,
	const EnergyVaryParam& evp, const EnergySolidParam& ft)
{
	double sum_eb = 0;
	double sum_ea = 0;
	double sum_eba = 0;
	double sum_eoop = 0;
	double sum_et = 0;
	double sum_evdw = 0;

	if (need_cal.eb) sum_eb = CalEB(ft.short_adj_list, evp.dr_tb, ft.bs_fast_tb, ft.pro_htb);
	if (need_cal.ea) sum_ea = CalEA(ft.short_adj_list, evp.dvar_tb, ft.ab_fast_tb, ft.pro_htb);
	if (need_cal.eb) sum_eba = CalEBA(ft.short_adj_list, evp.dr_tb, evp.dvar_tb, ft.sb_fast_tb, ft.pro_htb);
	if (need_cal.eoop) sum_eoop = CalEOOP(ft.short_adj_list, evp.chi_tb, ft.opb_fast_tb, ft.pro_htb);
	if (need_cal.et) sum_et = CalET(ft.short_adj_list, evp.phi_tb, ft.ti_fast_tb, ft.pro_htb);
	if (need_cal.evdw) sum_evdw = CalEVDW(ft.short_adj_list, ft.bfs_tb, evp.r_tb, ft.re_tb);

	double sum_E = sum_eb + sum_ea + sum_eba + sum_eoop + sum_et + sum_evdw;
	if (print_yes) PrintEnergy(sum_E, sum_eb, sum_ea, sum_eba, sum_eoop, sum_et, sum_evdw);

	return  sum_E;
}

//--- 初始三维坐标生成 ---
vector<Vec3> CalXYZ(const vector<double> alpha_tb, const  HASH_TB& pro_htb, const F_BS& bs_fast_tb, const ADJ_LIST& short_adj_list, SP_SpTable all_sp_tb)
{
	int ac = pro_htb.size();
	int nhc = alpha_tb.size();
	vector<Vec3> xyz_tb(ac, Vec3());
	vector<bool> rec_xyz(ac, false);
	for (int j = 0; j < nhc; j++)
	{
		int parent = j;
		//int m = mole_sp_tb[j];
		//if (m != 4) continue;
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();
		int m = nsize;

		//vector<SP_SpLine> new_m_arr(all_sp_tb[m]);
		vector<SP_ReLine> m_rect_arr(m, SP_ReLine());
		for (int sp_i = 0; sp_i < m; sp_i++)
		{
			SP_SpLine spline(all_sp_tb.GetSPLine(m, sp_i));
			spline.SpinTheta(alpha_tb[j]);
			SP_ReLine reline(spline);
			m_rect_arr[sp_i] = reline;
		}
		//cout << "MTB\n";
		//PrintTableTitle(SP_RETB_TITLE);
		//PrintSpecialVector(m_rect_arr);
		Matrix2<double> R;

		for (int i = 0; i < nsize; i++)
		{
			int child = nb_node[i].GetDesSeq();
			FAST_TABLE_INDEX fti = { pro_htb[parent],pro_htb[child] };

			double r0 = bs_fast_tb[fti].r0;
			//cout << "r0: " << r0 << endl;

			if (i == 0)
			{
				if (j == 0)
				{
					xyz_tb[child] = Vec3(0, 0, r0);
					rec_xyz[child] = true;
				}
				//Vec3 aaa(xyz_tb[child] - xyz_tb[parent]);
				//aaa.Print();
				CoordSys3 csys(xyz_tb[child] - xyz_tb[parent]);
				R = csys.GetMatrixR();
				//PrintCmdSepTitle("R");
				//R.Print();

			}
			else if (child > parent && rec_xyz[child] == false)
			{
				Vec3 adot = (m_rect_arr[i].GetVec() * R);
				xyz_tb[child] = r0 * adot + xyz_tb[parent];

				rec_xyz[child] = true;
			}

		}


	}
	return xyz_tb;

}

double CalSumEnergyByXYZ(bool print_yes, const NeedCal& need_cal, const vector<double> alpha_tb, SP_SpTable all_sp_tb,
	const EnergySolidParam& esp)
{
	vector<Vec3> xyz_tb = CalXYZ(alpha_tb, esp.pro_htb, esp.bs_fast_tb, esp.short_adj_list, all_sp_tb);


	R_TB dr_tb, r_tb;
	VAR_TB dvar_tb;
	PHI_TB phi_tb;
	CHI_TB chi_tb;

	if (need_cal.eb || need_cal.eba || need_cal.evdw)
	{
		PreCalDRTb(r_tb, dr_tb, xyz_tb, esp.short_adj_list, esp.bs_fast_tb, esp.pro_htb);
	}
	if (need_cal.ea || need_cal.eba)
	{
		PreCalDVarTb(dvar_tb, xyz_tb, esp.short_adj_list, esp.ab_fast_tb, esp.pro_htb);
	}
	if (need_cal.et)
	{
		PreCalPhiTb(phi_tb, xyz_tb, esp.short_adj_list);
	}
	if (need_cal.eoop)
	{
		PreCalChiTb(chi_tb, xyz_tb, esp.short_adj_list);
	}
	if (print_yes)
	{
		PrintCmdSepTitle("三维坐标");
		PrintSpecialVector(xyz_tb);
	}

	EnergyVaryParam evp{ r_tb,dr_tb,  dvar_tb,phi_tb,chi_tb };
	double sum_energy = CalSumEnergy(print_yes, need_cal, evp, esp);

	return sum_energy;
}



//--- 优化相关函数 ---
template<typename T>
bool JudgeStop(vector<T>new_tb, vector<T>old_tb, vector<T> sep_tb)
{
	int ssize = sep_tb.size();
	for (int i = 0; i < ssize; i++)
	{
		if (abs(old_tb[i] - new_tb[i]) > sep_tb[i])
		{
			return false;
		}
	}
	return true;
}

//类型识别



int main()
{
	PrintCmdSepTitle("基本参数表");
	vector<BSLine> bs_tb;
	vector<ABLine> ab_tb;
	vector<SBLine> sb_tb;
	vector<OPBLine> opb_tb;
	vector<TILine> ti_tb;
	vector<VDWLine> vdw_tb;

	string bs_filepath = GetMMFFPath(BS_FILENAME);
	string ab_filepath = GetMMFFPath(AB_FILENAME);
	string sb_filepath = GetMMFFPath(SB_FILENAME);
	string opb_filepath = GetMMFFPath(OPB_FILENAME);
	string ti_filepath = GetMMFFPath(TI_FILENAME);
	string vdw_filepath = GetMMFFPath(VDW_FILENAME);


	ReadTableByTitle(bs_filepath, BS_TB_TITLE, bs_tb);
	ReadTableByTitle(ab_filepath, AB_TB_TITLE, ab_tb);
	ReadTableByTitle(sb_filepath, SB_TB_TITLE, sb_tb);
	ReadTableByTitle(opb_filepath, OPB_TB_TITLE, opb_tb);
	ReadTableByTitle(ti_filepath, TI_TB_TITLE, ti_tb);
	ReadTableByTitle(vdw_filepath, VDW_TB_TITLE, vdw_tb);

	// 
	int max_row_count = 2;
	PrintTableTitle(BS_TB_SHORT_TITLE);
	PrintSpecialVector(bs_tb, "\t", max_row_count);

	PrintTableTitle(AB_TB_SHORT_TITLE);
	PrintSpecialVector(ab_tb, "\t", max_row_count);

	PrintTableTitle(SB_TB_SHORT_TITLE);
	PrintSpecialVector(sb_tb, "\t", max_row_count);

	PrintTableTitle(OPB_TB_SHORT_TITLE);
	PrintSpecialVector(opb_tb, "\t", max_row_count);

	PrintTableTitle(TI_TB_SHORT_TITLE);
	PrintSpecialVector(ti_tb, "\t", max_row_count);

	PrintTableTitle(VDW_TB_SHORT_TITLE);
	PrintSpecialVector(vdw_tb, "\t", max_row_count);


	PrintCmdSepTitle("杂化坐标表");

	vector<SP_SpLine> sp_sptb;
	ReadTableByTitle(GetSpheTbPath(3), SP_SPTB_TITLE, sp_sptb, YES);
	//PrintTableTitle(SP_SPTB_TITLE);
	//PrintSpecialVector(sp_sptb);

	SP_SpTable all_sp_tb;
	//all_sptb.Print();

	PrintCmdSepTitle("分子相关表");
	string now_smiles = "CCC";

	//map<vector<int>, BSVal> bs_kb_map = CompressTableToMap<BSLine, BSVal>(bs_tb);
	//PrintParamMap(bs_kb_map, "\t", max_row_count);

	vector<MNode> mnode_tb;
	string mnode_tb_path = GetMNodeTbPath(now_smiles);
	ReadTableByTitle(mnode_tb_path, MNODE_TB_TITLE, mnode_tb);
	PrintTableTitle(MNODE_TB_TITLE);
	PrintSpecialVector(mnode_tb);

	// --- 邻接表格式转换 ---
	vector<AdjLine> adj_tb;
	string adj_tb_path = GetAdjTbPath(now_smiles);
	ReadTableByTitle(adj_tb_path, ADJ_TB_TITLE, adj_tb);
	PrintTableTitle(ADJ_TB_TITLE);
	PrintSpecialVector(adj_tb);

	PrintCmdSepTitle("扩展邻接表");
	vector<NodeBonds> long_adj_list = AdjTbToBondTb(adj_tb, mnode_tb, NO);
	PrintSpecialVector(long_adj_list);

	PrintCmdSepTitle("普通邻接表");
	vector<NodeBonds> short_adj_list = AdjTbToBondTb(adj_tb, mnode_tb);
	PrintSpecialVector(short_adj_list);


	//---待修改的部分：根据分子节点类型，生成对应的MMFF94类型表---
	vector<int> mtype_tb = { 1,1,1,5,5,5,5,5,5,5,5 };




	MTYPE_SET mtype_set(mtype_tb.begin(), mtype_tb.end());

	vector<int> htb = GetHashTable(mtype_set);
	vector<int> pro_htb = GetProHashTable(mtype_tb, htb);


	//PrintSet(mtype_set);

	PrintCmdSepTitle("htb");
	PrintCommonVector(htb);

	PrintCmdSepTitle("pro_htb");
	PrintCommonVector(pro_htb);


	int ac = mtype_tb.size();
	int nhc = short_adj_list.size();

	F_BS bs_fast_tb(bs_tb, htb);
	F_AB ab_fast_tb(ab_tb, htb);
	F_SB sb_fast_tb(sb_tb, htb);
	F_OPB opb_fast_tb(opb_tb, htb);
	F_TI ti_fast_tb(ti_tb, htb);
	F_VDW vdw_fast_tb(vdw_tb, htb);

	BFS2_TB bfs_tb = Bfs2(ac, short_adj_list);

	RE_TB re_tb;
	PreCalRETb(re_tb, mtype_set, vdw_fast_tb, pro_htb);


	EnergySolidParam esp = { bs_fast_tb ,ab_fast_tb,sb_fast_tb,opb_fast_tb,ti_fast_tb,vdw_fast_tb,pro_htb,bfs_tb,re_tb,short_adj_list };

	PrintCmdSepTitle("三维坐标生成及能量计算-记时");


	double a3 = 10 / 180.0 * PI;

	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;
	QueryPerformanceCounter(&start);
	//------------------------------**计时开始**------------------------------

	vector<double> old_alpha_tb = { 0,0, a3 };
	vector<Vec3> xyz_tb = CalXYZ(old_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);

	PrintCmdSepTitle("三维坐标");
	PrintSpecialVector(xyz_tb);

	WriteTable(GetXYZTbPath(now_smiles, 1), XYZ_TB_TITLE, xyz_tb);

	R_TB dr_tb, r_tb;
	PreCalDRTb(r_tb, dr_tb, xyz_tb, short_adj_list, bs_fast_tb, pro_htb);

	VAR_TB dvar_tb;
	PreCalDVarTb(dvar_tb, xyz_tb, short_adj_list, ab_fast_tb, pro_htb);

	PHI_TB phi_tb;
	PreCalPhiTb(phi_tb, xyz_tb, short_adj_list);

	CHI_TB chi_tb;
	PreCalChiTb(chi_tb, xyz_tb, short_adj_list);

	EnergyVaryParam evp{ r_tb,dr_tb,dvar_tb,phi_tb,chi_tb };

	double sum_energy = CalSumEnergy(YES, NeedCal(),evp,esp);



	//------------------------------**根据能量优化参数**------------------------------

	// --- alpha 第一次优化 粗略 ---

	const double acc_energy = 1e-4;
	const double acc_angle = 1e-4;

	double asep = PI_HALF;
	vector<double> new_alpha_tb(nhc, PI);
	new_alpha_tb[0] = 0;
	vector<double> asep_tb = GetAlphaSepTable(short_adj_list);

	PrintCmdSepTitle("sp_tb");
	PrintCommonVector(asep_tb);
	NeedCal alpha_need_cal(false, false, false, false, true, true);


	while (!JudgeStop(new_alpha_tb, old_alpha_tb, asep_tb))
	{
		old_alpha_tb = new_alpha_tb;
		for (int i = 1; i < nhc; i++)
		{
			double now_asep = asep_tb[i];
			double min_energy = INFINITY;
			double min_angle = 0, now_angle = 0;

			do {
				new_alpha_tb[i] = now_angle;

				double now_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, new_alpha_tb, all_sp_tb,esp);
				if (now_energy < min_energy)
				{
					min_energy = now_energy;
					min_angle = now_angle;
				}

				now_angle += now_asep;

			} while (now_angle < PI_DOUBLE);
			new_alpha_tb[i] = min_angle;
		}
	}
	vector<Vec3> a_xyz_tb = CalXYZ(new_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);

	PrintCmdSepTitle("第一次优化后三维坐标");
	PrintSpecialVector(a_xyz_tb);
	WriteTable(GetXYZTbPath(now_smiles, 2), XYZ_TB_TITLE, a_xyz_tb);

	// --- alpha 第二次优化 精细 ---

	double old_energy = INFINITY;
	double new_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, new_alpha_tb, all_sp_tb,esp);

	cout << "\nE: " << new_energy << endl;

	old_alpha_tb = new_alpha_tb;

	while (!JudgeStop<double>({ new_energy }, { old_energy }, { acc_energy }))
	{
		old_energy = new_energy;
		for (int i = 1; i < nhc; i++)
		{
			vector<double>now_alpha_tb = new_alpha_tb;
			double left_angle = new_alpha_tb[i] - asep_tb[i];
			double right_angle = new_alpha_tb[i] + asep_tb[i];

			now_alpha_tb[i] = left_angle;

			double left_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

			now_alpha_tb[i] = right_angle;

			double right_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

			double mid_energy = INFINITY;

			while (!JudgeStop<double>({ left_energy }, { right_energy }, { acc_energy }) &&
				!JudgeStop<double>({ left_angle }, { right_angle }, { acc_angle }))
			{
				double sep_angle = (right_angle - left_angle) * OPT_RATIO;
				double mid_left_angle = left_angle + sep_angle;
				double mid_right_angle = right_angle - sep_angle;

				now_alpha_tb[i] = mid_left_angle;
				double mid_left_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

				now_alpha_tb[i] = mid_right_angle;
				double mid_right_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

				if (left_energy < right_energy)
				{
					right_angle = mid_left_angle;
					right_energy = mid_left_energy;
				}
				else
				{
					left_angle = mid_right_angle;
					left_energy = mid_right_energy;
				}
			}
			new_alpha_tb[i] = (left_angle + right_angle) / 2;
		}
		new_energy = CalSumEnergyByXYZ(NO, alpha_need_cal, new_alpha_tb, all_sp_tb, esp);
	}

	vector<Vec3> b_xyz_tb = CalXYZ(new_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);
	PrintCmdSepTitle("第二次优化三维坐标");
	PrintSpecialVector(a_xyz_tb);
	cout << "\nE: " << new_energy << endl;
	WriteTable(GetXYZTbPath(now_smiles, 3), XYZ_TB_TITLE, b_xyz_tb);

	//------------------------------**计时结束**------------------------------
	QueryPerformanceCounter(&stop);

	double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
	std::cout << "\n(^-^)Time taken by function: " << duration << " ms" << std::endl;



	//WriteTable(GetXYZTbPath(now_smiles), XYZ_TB_TITLE, xyz_tb);
	//PrintEnergy(sum_E, sum_eb, sum_ea, sum_eba, sum_eoop, sum_et, sum_evdw);

}


