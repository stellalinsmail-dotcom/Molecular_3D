#include"MMFF94Cal.h"


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
vector<double> GetAlphaSepTable(const ADJ_LIST& short_adj_list)
{
	int nhc = short_adj_list.size();
	vector<double> a(nhc, NAN);
	for (int i = 0; i < nhc; i++)
	{
		a[i] = PI_HALF / max(short_adj_list[i].GetBonds().size() - 1, 1);
		//a[i] = min(PI_HALF / (short_adj_list[i].GetBonds().size() - 1), PI_HALF / 3);
	}
	return a;

}
// --- 简单参数先导计算 ---

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

//--- 初始三维坐标生成 ---
vector<Vec3> CalXYZ(const vector<double>& alpha_tb, const  HASH_TB& pro_htb, const F_BS& bs_fast_tb, const MNODE_TB& mnode_tb, const ADJ_LIST& short_adj_list, const SP_SpTable& all_sp_tb)
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
		int m = nsize, startseq = 0;

		string now_sym = mnode_tb[j].GetSym();
		if (now_sym == "O" || now_sym == "N")
		{
			m = 4;
		}
		rec_xyz[parent] = true;
		//vector<SP_SpLine> new_m_arr(all_sp_tb[m]);
		vector<SP_ReLine> m_rect_arr(m, SP_ReLine());
		if (m >= 1)
		{
			for (int sp_i = 0; sp_i < m; sp_i++)
			{
				SP_SpLine spline(all_sp_tb.GetSPLine(m, sp_i));
				spline.SpinTheta(alpha_tb[j]);
				SP_ReLine reline(spline);
				m_rect_arr[sp_i] = reline;
			}
		}
		//cout << "MTB\n";
		//PrintTableTitle(SP_RETB_TITLE);
		//PrintSpecialVector(m_rect_arr);
		Matrix2<double> R;
		Vec3 oz;
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
				oz = Vec3(xyz_tb[child] - xyz_tb[parent]);
				CoordSys3 csys(oz);
				R = csys.GetMatrixR();
				//PrintCmdSepTitle("R");
				//R.Print();

			}
			else if (rec_xyz[child] == true && i == 1)
			{
				Vec3 ox_ref(xyz_tb[child] - xyz_tb[parent]);
				CoordSys3 csys(oz, ox_ref);
				R = csys.GetMatrixR();
				//cout <<"重新建系：" << parent << "\t" << child << endl;
			}
			else if (child > parent && rec_xyz[child] == false)
			{
				Vec3 adot = (m_rect_arr[i].GetVec() * R);
				xyz_tb[child] = r0 * adot + xyz_tb[parent];

				rec_xyz[child] = true;
				//cout << parent << "\t" << child << endl;
			}

		}

	}
	return xyz_tb;

}

double CalSumEnergyByXYZ(bool print_yes, XYZ_TB& xyz_tb, const NeedCal& need_cal, const vector<double>& alpha_tb, const SP_SpTable& all_sp_tb,
	const EnergySolidParam& esp, int limitpos)
{
	vector<NodeBonds> now_adj_list;
	if (limitpos != -1)
	{
		now_adj_list = CutAdjList(esp.short_adj_list, limitpos);
	}
	else
	{
		now_adj_list = esp.short_adj_list;
	}
	xyz_tb = CalXYZ(alpha_tb, esp.pro_htb, esp.bs_fast_tb, esp.mnode_tb, esp.short_adj_list, all_sp_tb);


	R_TB dr_tb, r_tb;
	VAR_TB dvar_tb;
	PHI_TB phi_tb;
	CHI_TB chi_tb;

	if (need_cal.eb || need_cal.eba || need_cal.evdw)
	{
		PreCalDRTb(r_tb, dr_tb, xyz_tb, now_adj_list, esp.bs_fast_tb, esp.pro_htb);
	}
	if (need_cal.ea || need_cal.eba)
	{
		PreCalDVarTb(dvar_tb, xyz_tb, now_adj_list, esp.ab_fast_tb, esp.pro_htb);
	}
	if (need_cal.et)
	{
		PreCalPhiTb(phi_tb, xyz_tb, now_adj_list);
	}
	if (need_cal.eoop)
	{
		PreCalChiTb(chi_tb, xyz_tb, now_adj_list);
	}
	//if (print_yes)
	//{
	//	PrintCmdSepTitle("三维坐标");
	//	PrintSpecialVector(xyz_tb);
	//}

	EnergyVaryParam evp{ r_tb,dr_tb,  dvar_tb,phi_tb,chi_tb };
	double sum_energy = CalSumEnergy(print_yes, need_cal, evp, esp, limitpos);

	return sum_energy;
}

vector<string> GetMSymTb(int ac, const vector<MNode>& nodetb, const  ADJ_LIST& short_adj_list)
{
	//先求各个原子的MSym
	int nhc = short_adj_list.size();

	vector<string> ori_sym_tb(ac, "");

	for (int j = 0; j < nhc; j++)
	{
		vector<PointTo> nb_node = short_adj_list[j].GetBonds();
		int nsize = nb_node.size();

		string parent_sym = nodetb[j].GetSym();
		string now_sym = parent_sym;
		string h_cmpl_sym = "";
		for (int i = 0; i < nsize; i++)
		{
			int child = nb_node[i].GetDesSeq();

			//if (child==nhc)
			//{
			//	h_cmpl_sym=
			//}

			//非氢
			if (child < nhc)
			{
				string bond_sym = nb_node[i].GetBondSymbol();
				string child_sym = nodetb[child].GetSym();
				if (bond_sym == DOUBLE_BOND) now_sym += bond_sym + child_sym;
				else if (bond_sym == TRIPLE_BOND) now_sym += "%";

				//else
				//{
				//	if (parent_sym != "C")
				//	{
				//		now_sym += child_sym;
				//	}
				//}
			}
			//氢
			else
			{
				if (parent_sym != "C")
				{
					ori_sym_tb[child] = "H" + now_sym;
				}
				else
				{
					ori_sym_tb[child] = "H" + parent_sym;
				}
			}
		}
		if (now_sym == parent_sym)
		{
			now_sym += "R";
		}
		ori_sym_tb[j] = now_sym;
	}
	//for (int i = 0; i < ac; i++)
	//{
	//	cout << i << "\t" << ori_sym_tb[i] << endl;
	//}

	return ori_sym_tb;
}

vector<int> GetMTypeTb(int ac, const vector<string>& msym_tb, const map<string, int>& msym_map)
{
	vector<int> mtype_tb(ac, NAN);
	for (int i = 0; i < ac; i++)
	{
		auto iter = msym_map.find(msym_tb[i]);
		if (iter != msym_map.end())
		{
			mtype_tb[i] = iter->second;
		}
		else
		{
			cout << "Warning: 未找到MSym " << msym_tb[i] << " 对应的MType！" << endl;
		}
	}
	return mtype_tb;
}

EnergySolidParam GenEnergySolidParam(int& ac, const EnergyFundTable& eft, const Mole& d, bool print_yes, int max_row_count)
{
	MNODE_TB mnode_tb = d.GetNodeTable();
	ADJ_LIST old_bond_tb = d.GetBondTable();


	ADJ_LIST short_adj_list;
	ac = ExpandBondTb(short_adj_list, old_bond_tb, mnode_tb);

	vector<string> ori_sym_tb = GetMSymTb(ac, mnode_tb, short_adj_list);
	vector<int> mtype_tb = GetMTypeTb(ac, ori_sym_tb, eft.msym_map);


	MTYPE_SET mtype_set(mtype_tb.begin(), mtype_tb.end());

	vector<int> htb = GetHashTable(mtype_set);
	vector<int> pro_htb = GetProHashTable(mtype_tb, htb);


	F_BS bs_fast_tb(eft.bs_tb, htb);
	F_AB ab_fast_tb(eft.ab_tb, htb);
	F_SB sb_fast_tb(eft.sb_tb, htb);
	F_OPB opb_fast_tb(eft.opb_tb, htb);
	F_TI ti_fast_tb(eft.ti_tb, htb);
	F_VDW vdw_fast_tb(eft.vdw_tb, htb);

	BFS2_TB bfs_tb = Bfs2(ac, short_adj_list);

	RE_TB re_tb;

	PreCalRETb(re_tb, mtype_set, vdw_fast_tb, pro_htb);


	EnergySolidParam esp = { bs_fast_tb ,ab_fast_tb,sb_fast_tb,opb_fast_tb,ti_fast_tb,vdw_fast_tb,pro_htb,bfs_tb,re_tb,short_adj_list,mnode_tb };

	if (print_yes)
	{
		PrintCmdSepTitle("一般扩展邻接表");
		PrintSpecialVector(short_adj_list, "\t", max_row_count);

		PrintCmdSepTitle("MType 表格查找");
		cout << "Index\tMSymCal\tMType\n";
		for (int i = 0; i < min(ac, max_row_count); i++)
		{
			cout << i << "\t" << ori_sym_tb[i] << "\t" << mtype_tb[i] << endl;
		}
		PrintCmdSepTitle("htb");
		PrintCommonVector(htb);

		PrintCmdSepTitle("pro_htb");
		PrintCommonVector(pro_htb);
	}
	return esp;
}
