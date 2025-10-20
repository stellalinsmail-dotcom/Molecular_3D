//#include <array>
#include <map>

#include "types.h"
#include "Atom.h"
#include "Molecule.h"

//#include "Atom.h"
//#include "


#define BS_TB_TITLE "BT,MTYPE_I,MTYPE_J,KB,R0,SOURCE"
#define AB_TB_TITLE "ANGLE_TYPE,MTYPE_I,MTYPE_J,MTYPE_K,KA_IJK,THETA_0,COMMENT_ORIGIN"
#define SB_TB_TITLE "SBT,MTYPE_I,MTYPE_J,MTYPE_K,KBA_IJK,KBA_KJI,SOURCE"
#define OPB_TB_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,KOOP,Source"
#define TI_TB_TITLE	"TT,MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,V1,V2,V3,Source"
#define VDW_TB_TITLE "MTYPE,ALPHA_I,N_I,A_I,G_I,DA,SYMBOL,ORIGIN"

#define BS_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,KB,R0"
#define AB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,KA_IJK,THETA_0"
#define SB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,KBA_IJK,KBA_KJI"
#define OPB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE,KOOP"
#define TI_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,V1,V2,V3"
#define VDW_TB_SHORT_TITLE "MTYPE,ALPHA_I,N_I,A_I,G_I"

#define MMFF_CSV_FOLDER "File/Values/MMFF94_CSV/"
#define MMFF_OUTPUT_FOLDER "File/Output/"

#define BS_FILENAME  "6_MMFFBOND.csv"
#define AB_FILENAME "8_MMFFANG.csv"
#define SB_FILENAME "9_MMFFSTBN.csv"
#define OPB_FILENAME "11_MMFFOOP.csv"
#define TI_FILENAME "12_MMFFTOR.csv"
#define VDW_FILENAME "13_MMFFVDW.csv"


#define PI 3.14159265358979323846
#define cosd(x) cos(x * PI / 180.0)
#define sind(x) sin(x * PI / 180.0)

using namespace std;


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

	vector<int> GetMTypeVec()const
	{
		return vector<int>{ mtype_i, mtype_j };
	}
	vector<int> GetRevMTypeVec()const
	{
		return vector<int>{ mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1]) || (mtype_i == mtype[1] && mtype_j == mtype[0]);
	}

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
		theta_0 = stod(info[titlenum[5]]);
	}
	vector<int> GetMTypeVec()const
	{
		return vector<int>{ mtype_i, mtype_j, mtype_k };
	}
	vector<int> GetRevMTypeVec()const
	{
		return vector<int>{ mtype_k, mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0]);
	}
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
	vector<int> GetMTypeVec()const
	{
		return vector<int>{ mtype_i, mtype_j, mtype_k };
	}
	vector<int> GetRevMTypeVec()const
	{
		return vector<int>{ mtype_k, mtype_j, mtype_i };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0]);
	}
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
	vector<int> GetMTypeVec()const
	{
		return vector<int>{ mtype_i, mtype_j, mtype_k, mtype_l };
	}
	vector<int> GetRevMTypeVec()const
	{
		return vector<int>{ mtype_k, mtype_j, mtype_i, mtype_l };
	}
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype_i == mtype[0] && mtype_j == mtype[1] && mtype_k == mtype[2] && mtype_l == mtype[3]) ||
			(mtype_i == mtype[2] && mtype_j == mtype[1] && mtype_k == mtype[0] && mtype_l == mtype[3]);
	}
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
	vector<int> GetMTypeVec()const
	{
		return vector<int>{ mtype_i, mtype_j, mtype_k, mtype_l };
	}
	vector<int> GetRevMTypeVec()const
	{
		return vector<int>{ mtype_l, mtype_k, mtype_j, mtype_i };
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
	vector<int> GetMTypeVec()const { return vector<int> {mtype_i}; }
	vector<int> GetRevMTypeVec()const { return vector<int> {mtype_i}; }
	bool IsMTypeMatch(vector<int> mtype)const
	{
		if (mtype.size() < mtype_count) return false;
		return (mtype[0] == mtype_i);
	}
	double GetAlphaI()const { return alpha_i; }
	double GetNI()const { return n_i; }
	double GetAI()const { return a_i; }
	double GetGI()const { return g_i; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(VDW_TB_SHORT_TITLE, sep);
			cout << endl;
		}
		cout << mtype_i << sep << alpha_i << sep << n_i << sep << a_i << sep << g_i << endl;
	}

};



// --- 表格处理 ---
template<typename T>
map<vector<int>, double> CompressTableToMap(const vector<T>& table, double (T::* func)()const)
{
	map<vector<int>, double> table_map;
	for (auto& line : table)
	{
		table_map[line.GetMTypeVec()] = (line.*func)();
		table_map[line.GetRevMTypeVec()] = (line.*func)();
	}
	return table_map;
}

void AdjTbToBondTb(vector<NodeBonds>& bond_tb, const vector<AdjLine>& adj_tb)
{
	bond_tb.clear();
	int max_node_seq = 0;
	for (auto& line : adj_tb)
	{
		if (line.GetSeqI() > max_node_seq) max_node_seq = line.GetSeqI();
		if (line.GetSeqJ() > max_node_seq) max_node_seq = line.GetSeqJ();
	}
	int max_size = max_node_seq + 1;
	bond_tb.resize(max_size, NodeBonds(0));
	for (int i = 0; i < max_size; i++)
	{
		bond_tb[i] = NodeBonds(i);
	}
	for (auto& line : adj_tb)
	{
		PointTo pto_j(line.GetBondSym(), line.GetSeqJ());
		bond_tb[line.GetSeqI()].AddBond(pto_j);
		//PointTo pto_i(line.GetBondSym(), line.GetSeqI());
		//bond_tb[line.GetSeqJ()].AddBond(pto_i);
	}
}

// ---三维向量处理 ---
class Vec3 {
private:
	double len; //长度缓存，-1表示未计算
		
public:
	double x;
	double y;
	double z;

	Vec3() :x(0), y(0), z(0),len(0) {}
	Vec3(double xval, double yval, double zval) :x(xval), y(yval), z(zval),len(-1){}
	double GetLen()
	{
		if (len==-1)
		{
			len = sqrt(x * x + y * y + z * z);
		}
		return len;
	}
	friend double Cdot(const Vec3& a, const Vec3& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	friend Vec3 Cross(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
	}
	friend Vec3 operator-(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	friend Vec3 operator+(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	friend double GetAngle_Rad(Vec3& a,Vec3& b)
	{
		double theta = Cdot(a, b) / (a.GetLen() * b.GetLen());
		if (theta > 1.0) theta = 1.0;
		if (theta < -1.0) theta = -1.0;
		return acos(theta);
	}
	bool Normalize()
	{
		if (len == 0) return false;
		x /= len;
		y /= len;
		z /= len;
		len = 1.0;
		return true;
	}
	void Print(string sep = "\t")const
	{
		cout << x << sep << y << sep << z << endl;
	}
	friend ostream& operator<<(ostream& out, const Vec3& v)
	{
		out << v.x << "," << v.y << "," << v.z;
		return out;
	}
};

class CoordSys3 {
public: 
	Vec3 ox;
	Vec3 oy;
	Vec3 oz;
	CoordSys3(Vec3 oz)
	{

	}
};

// --- MMFF94能量计算函数定义 ---
double GetEB(double dr_ij, double kb)
{
	return 143.9525 * kb / 2 * pow(dr_ij, 2) * (1 - 2 * dr_ij + 2.33333333333 * pow(dr_ij, 2));
}
double GetEA(double dvar, double ka)
{
	return 0.043844 * ka / 2 * pow(dvar, 2) * (1 + -0.000122 * dvar);
}
double GetEBA(double dr_ij, double dr_kj, double dvar, double kba_ijk, double kba_kji)
{
	return 2.51210 * (kba_ijk * dr_ij + kba_kji * dr_kj) * dvar;
}
double GetET(double phi, double v1, double v2, double v3)
{
	return 0.5 * (v1 * (1 + cosd(phi)) + v2 * (1 - cosd(2 * phi)) + v3 * (1 + cosd(3 * phi)));
}
double GetEOOP(double chi_ijkl, double koop)
{
	return 0.034844 * koop / 2 * pow(chi_ijkl, 2);
}
double GetEVDW(double r_ij, double rv_ij, double epsilon_ij)
{
	return epsilon_ij * pow(1.07 * rv_ij / (r_ij + 0.07 * rv_ij), 7) * (1.12 * pow(rv_ij, 7) / (pow(r_ij, 7) + 0.12 * pow(rv_ij, 7)) - 2);
}

string GetMMFFPath(string filename)
{
	return MMFF_CSV_FOLDER + filename;
}
string GetMNodeTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/MNodeTable.csv";
}
string GetAdjTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/BondTable.csv";
}

//--- 字典打印输出 ---
void PrintMap(const map<vector<int>, double>& m, const string& sep = "\t", int max_row_count = -1)
{
	PrintTableTitle("MType,Keys", sep);
	int rowcount = 0;
	for (auto& pair : m)
	{
		cout << "(";
		PrintCommonVector(pair.first, ",");
		cout << ")" << sep << pair.second << endl;
		if (max_row_count > 0)
		{
			rowcount++;
			if (rowcount >= max_row_count) break;
		}
	}
}

int main()
{
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;
	QueryPerformanceCounter(&start);
	//------------------------------**计时开始**------------------------------

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


	string now_smiles = "CCC";

	map<vector<int>, double> bs_kb_map = CompressTableToMap(bs_tb, &BSLine::GetKB);
	PrintMap(bs_kb_map, "\t", max_row_count);

	vector<MNode> mnode_tb;
	string mnode_tb_path = GetMNodeTbPath(now_smiles);
	ReadTableByTitle(mnode_tb_path, MNODE_TB_TITLE, mnode_tb);
	PrintTableTitle(MNODE_TB_TITLE);
	PrintSpecialVector(mnode_tb);

	vector<AdjLine> adj_tb;
	string adj_tb_path = GetAdjTbPath(now_smiles);
	ReadTableByTitle(adj_tb_path, ADJ_TB_TITLE, adj_tb);
	PrintTableTitle(ADJ_TB_TITLE);
	PrintSpecialVector(adj_tb);

	vector<NodeBonds> bond_tb;
	AdjTbToBondTb(bond_tb, adj_tb);
	//PrintTableTitle(ADJ_TB_TITLE);
	PrintSpecialVector(bond_tb);


	//---待修改的部分：根据分子节点类型，生成对应的MMFF94类型表---
	vector<int> mtype_tb{ 1,1,1,5,5,5,5,5,5,5,5 };



	//------------------------------**计时结束**------------------------------
	QueryPerformanceCounter(&stop);

	double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
	std::cout << "\n(^-^)Time taken by function: " << duration << " ms" << std::endl;
}



