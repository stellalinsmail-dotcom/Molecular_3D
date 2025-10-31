#pragma once
#ifndef MMFF94CAL_H
#define MMFF94CAL_H

#include "MMFF94Table.h"
#include "Molecule.h"

class NeedCal {
public:
	bool eb;
	bool ea;
	bool eba;
	bool eoop;
	bool et;
	bool evdw;
	NeedCal(bool is_eb = true, bool is_ea = true, bool is_eba = true, bool is_eoop = false, bool is_et = true, bool is_evdw = true) :
		eb(is_eb), ea(is_ea), eba(is_eba), et(is_et), eoop(is_eoop), evdw(is_evdw) {
	}
	NeedCal(const NeedCal& nc) :eb(nc.eb), ea(nc.ea), eba(nc.eba), et(nc.et), eoop(nc.eoop), evdw(nc.evdw) {
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
	MNODE_TB mnode_tb;
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

vector<int> GetHashTable(const MTYPE_SET& mtype_set);

vector<int> GetProHashTable(const vector<int>& mtype_tb, const vector<int>& htb);


vector<double> GetAlphaSepTable(const ADJ_LIST& short_adj_list);

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
	Vec3 n_bc = Cross(c, b);
	double phi_rad = acos(Cdot(n_ab, n_bc) / (n_ab.GetLen() * n_bc.GetLen()));

	//cout << endl;
	//n_ab.Print();
	//n_bc.Print();
	//cout << "val: " << Cdot(n_ab, n_bc) / (n_ab.GetLen() * n_bc.GetLen()) << endl;
	//cout << "phi_1: " << phi_rad << endl;

	//double k = -Cdot(a, b) / Cdot(b, b);
	//Vec3 d = a + k * b;
	//bool judge = (Cdot(d, c) < 0);
	//if (judge) phi_rad = PI - phi_rad;
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
VDWProVal PreCalRE(VDWVal vdw_i, VDWVal vdw_j);

void PreCalRETb(RE_TB& re_tb, MTYPE_SET mtype_set, const F_VDW& vdw_fast_tb, const HASH_TB& pro_htb);
void PreCalDRTb(R_TB& r_tb, R_TB& dr_tb, const XYZ_TB& xyz_tb, ADJ_LIST short_adj_list, const F_BS& bs_fast_tb, const HASH_TB& pro_htb);

void PreCalDVarTb(VAR_TB& dvar_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list, const F_AB& ab_fast_tb, const HASH_TB& pro_htb);
void PreCalPhiTb(PHI_TB& phi_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list);
void PreCalChiTb(CHI_TB& chi_tb, const XYZ_TB& xyz_tb, const ADJ_LIST short_adj_list);
BFS2_TB Bfs2(int ac, const ADJ_LIST& short_adj_list);


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
inline double CalEB(const ADJ_LIST& short_adj_list, const R_TB& dr_tb, const F_BS& bs_fast_tb, const HASH_TB& pro_htb)
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

				double now_eb = GetEB(dr_tb[parent][child], bs_fast_tb[fti].kb);
				sum_eb += isnan(now_eb) ? 0 : now_eb;
				//cout << "EB: " << GetEB(dr_tb[parent][i], bs_fast_tb[fti].kb) << endl;
			}
		}

	}
	return sum_eb;
}
inline double CalEA(const ADJ_LIST& short_adj_list, const VAR_TB& dvar_tb, const F_AB& ab_fast_tb, const HASH_TB& pro_htb)
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

				//cout << "EA: " << GetEA(dvar_tb[parent][i][k], ka_ijk) << endl;
				double now_ea = GetEA(dvar_tb[parent][i][k], ka_ijk);
				sum_ea += isnan(now_ea) ? 0 : now_ea;
			}
		}
	}
	return sum_ea;
}
inline double CalEBA(const ADJ_LIST& short_adj_list, const R_TB& dr_tb, const VAR_TB& dvar_tb, const F_SB& sb_fast_tb, const HASH_TB& pro_htb)
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
			int nsize = nb_node.size();
			for (int k = i + 1; k < nsize; k++)
			{
				int child_k = nb_node[k].GetDesSeq();
				FAST_TABLE_INDEX fti = { pro_htb[child_i],pro_htb[parent],pro_htb[child_k] };
				double kba_ijk = sb_fast_tb[fti].kba_ijk;
				double kba_kji = sb_fast_tb[fti].kba_kji;
				double now_eba = GetEBA(dr_tb[parent][child_i], dr_tb[parent][child_k], dvar_tb[parent][i][k], kba_ijk, kba_kji);;
				sum_eba += isnan(now_eba) ? 0 : now_eba;
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
inline double CalEOOP(const ADJ_LIST& short_adj_list, const CHI_TB& chi_tb, const F_OPB& opb_fast_tb, const HASH_TB& pro_htb)
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
				double now_eoop = GetEOOP(chi_ijkl, koop);
				sum_eoop += isnan(now_eoop) ? 0 : now_eoop;
			}

		}

	}
	return sum_eoop;
}
inline double CalET(const ADJ_LIST& short_adj_list, const PHI_TB& phi_tb, const F_TI& ti_fast_tb, const HASH_TB& pro_htb)
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
								double now_et = GetET(phi_rad, v1, v2, v3);
								sum_et += isnan(now_et) ? 0 : now_et;

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
inline double CalEVDW(const ADJ_LIST& short_adj_list, const BFS2_TB& bfs_tb, const R_TB& r_tb, const RE_TB& re_tb)
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
			double now_evdw = GetEVDW(r_tb[j][ci], re_tb[j][ci].rv, re_tb[j][ci].e);
			sum_evdw += isnan(now_evdw) ? 0 : now_evdw;
			//cout << endl;
			//cout << j << "\t" << ci << endl;
			//cout << r_tb[j][ci] << "\t" << re_tb[j][ci].rv << "\t" << re_tb[j][ci].e << endl;
			//cout << "EVDW: " << GetEVDW(r_tb[j][ci], re_tb[j][ci].rv, re_tb[j][ci].e) << endl;
		}
	}
	return  sum_evdw;
}


inline ADJ_LIST CutAdjList(const ADJ_LIST& short_adj_list, int limitpos)
{
	int nhc = short_adj_list.size();
	int newsize = min(limitpos, nhc);
	ADJ_LIST a;

	for (int i = 0; i < newsize; i++)
	{
		vector<PointTo> ptb(short_adj_list[i].GetBonds());
		NodeBonds nb(i);
		for (auto& pt : ptb)
		{
			int nowseq = pt.GetDesSeq();
			if (nowseq < limitpos || nowseq>nhc)
			{
				nb.AddBond(pt);
			}
		}
		a.push_back(nb);
	}
	return a;
}


inline double CalSumEnergy(bool print_yes, const  NeedCal& need_cal,
	const EnergyVaryParam& evp, const EnergySolidParam& esp, int limitpos = -1)
{
	vector<NodeBonds> now_adj_list;
	if (limitpos != -1)
	{
		now_adj_list = now_adj_list = CutAdjList(esp.short_adj_list, limitpos);
	}
	else
	{
		now_adj_list = esp.short_adj_list;
	}
	double sum_eb = 0;
	double sum_ea = 0;
	double sum_eba = 0;
	double sum_eoop = 0;
	double sum_et = 0;
	double sum_evdw = 0;

	if (need_cal.eb) sum_eb = CalEB(now_adj_list, evp.dr_tb, esp.bs_fast_tb, esp.pro_htb);
	if (need_cal.ea) sum_ea = CalEA(now_adj_list, evp.dvar_tb, esp.ab_fast_tb, esp.pro_htb);
	if (need_cal.eb) sum_eba = CalEBA(now_adj_list, evp.dr_tb, evp.dvar_tb, esp.sb_fast_tb, esp.pro_htb);
	if (need_cal.eoop) sum_eoop = CalEOOP(now_adj_list, evp.chi_tb, esp.opb_fast_tb, esp.pro_htb);
	if (need_cal.et) sum_et = CalET(now_adj_list, evp.phi_tb, esp.ti_fast_tb, esp.pro_htb);
	if (need_cal.evdw) sum_evdw = CalEVDW(now_adj_list, esp.bfs_tb, evp.r_tb, esp.re_tb);

	//sum_et = (sum_et == NAN) ? 0 : sum_et;
	//sum_et = 0;
	double sum_E = sum_eb + sum_ea + sum_eba + sum_eoop + sum_et + sum_evdw;
	if (print_yes) PrintEnergy(sum_E, sum_eb, sum_ea, sum_eba, sum_eoop, sum_et, sum_evdw);

	return  sum_E;
}

//--- 初始三维坐标生成 ---
vector<Vec3> CalXYZ(const vector<double>& alpha_tb, const  HASH_TB& pro_htb, const F_BS& bs_fast_tb, const MNODE_TB& mnode_tb, const ADJ_LIST& short_adj_list, const SP_SpTable& all_sp_tb,bool print_yes=false);

double CalSumEnergyByXYZ(bool print_yes, XYZ_TB& xyz_tb, const NeedCal& need_cal, const vector<double>& alpha_tb, const SP_SpTable& all_sp_tb,
	const EnergySolidParam& esp, int limitpos = -1);

vector<string> GetMSymTb(int ac, const vector<MNode>& nodetb, const  ADJ_LIST& short_adj_list);

vector<int> GetMTypeTb(int ac, const vector<string>& msym_tb, const map<string, int>& msym_map);

EnergySolidParam GenEnergySolidParam(int& ac, const EnergyFundTable& eft, const Mole& d, bool print_yes = false, int max_row_count = -1);


#endif
