#include"AlphaBetaOpt.h"


inline vector<int> GenerateRandSeq(int size)
{
	//使用cstdlib和ctime生成随机数
	vector<int> seq(size);
	for (int i = 0; i < size; i++)
	{
		seq[i] = i;
	}
	unsigned seed = static_cast<unsigned>(time(0));
	shuffle(seq.begin(), seq.end(), default_random_engine(seed));
	return seq;
}


//粗略优化函数
//V1_DTB AlphaOptRough(bool has_circle, double& alpha_energy, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp)
//{
//	vector<double> alpha_tb;
//
//	return alpha_tb;
//}
V1_DTB AlphaOpt(bool has_circle, const V2_DTB& beta_tb, OptRecVal& alpha_result, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp, bool detail_print, bool savedata_yes, string now_smiles)
{
	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;
	QueryPerformanceCounter(&start);

	//------------------------------**计时开始**------------------------------

	vector<OptimizationRecord> opt_records;

	LARGE_INTEGER last_record_time = start;
	const double record_interval_ms = 0.1;  // 记录间隔（毫秒）
	long long total_iteration_count = 0;          // 总迭代次数
	double current_min_energy = INFINITY;   // 当前最低


	int nhc = esp.short_adj_list.size();


	int turn1count = 0;
	int rough_turncount = 0;
	int turn2count = 0;
	double asep = PI_HALF;

	ABOpt old_ab_opt(nhc);
	ABOpt new_ab_opt(nhc);
	ABOpt now_ab_opt(nhc);
	old_ab_opt.beta_tb = beta_tb;
	new_ab_opt.beta_tb = beta_tb;
	now_ab_opt.beta_tb = beta_tb;
	//vector<double> old_alpha_tb(nhc, 0);
	//vector<double> new_alpha_tb(nhc, PI);


	old_ab_opt.alpha_tb = V1_DTB(nhc, 0);
	new_ab_opt.alpha_tb = V1_DTB(nhc, PI);


	V1_DTB& new_alpha_tb = new_ab_opt.alpha_tb;
	V1_DTB& old_alpha_tb = old_ab_opt.alpha_tb;



	new_alpha_tb[0] = 0;


	int cyc = 0;
	for (int i = 0; i < nhc; i++)
	{
		//bool continue_yes = false;
		vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
		for (int j = 0; j < nb_node.size(); j++)
		{
			int child = nb_node[j].GetDesSeq();
			string bond_sym = nb_node[j].GetBondSymbol();
			if (child < i && IsSpecialBond(bond_sym))
			{
				//continue_yes = true;
				cyc++;
				break;
			}
		}
		//if (continue_yes) continue;
	}
	if (cyc == nhc - 1) return old_alpha_tb;


	vector<double> rec_asep_tb = GetAlphaSepTable(esp.short_adj_list);
	vector<double> asep_tb = rec_asep_tb;

	//PrintCmdSepTitle("sp_tb");
	//PrintCommonVector(asep_tb);
	NeedCal alpha_need_cal;

	if (!has_circle)
	{
		alpha_need_cal = NeedCal(false, false, false, false, true, true);
	}

	double start_energy = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

	const double acc_energy = 1e-2;
	const double acc_angle = 1e-3;

	const int angle_half_turn = min(nhc * 5, 20);
	const int print_turn = 20;

	double old_energy_a = INFINITY;
	double new_energy_a = start_energy;

	const int wait_turn_max = min(nhc * 2, 30);
	int wait_turn_count = 0;
	bool break_yes_a = false;

	const int max_rough_turn = 10;
	const int max_rough_big_turn = nhc * 200;

	double rough_min_energy = INFINITY;
	vector<double> rough_min_alpha_tb = old_alpha_tb;
	bool rough_big_break_yes = false;
AlphaOpt1:
	{
		// --- alpha 第一次优化 粗略 ---
		while (!JudgeStop(new_alpha_tb, old_alpha_tb, asep_tb))
		{
			if (wait_turn_count >= wait_turn_max)
			{
				if (detail_print) cout << "等待超过最大轮次，结束粗略优化。" << endl;
				new_energy_a = old_energy_a;
				new_alpha_tb = old_alpha_tb;

				break_yes_a = true;
				break;
			}
			if (new_energy_a <= old_energy_a)
			{
				//for (int i = 0; i < nhc; i++)
				//{
				//	new_alpha_tb[i] = fmod((new_alpha_tb[i] + old_alpha_tb[i]) / 2,PI_DOUBLE);
				//}
				old_alpha_tb = new_alpha_tb;
				old_energy_a = new_energy_a;
				wait_turn_count = 0;
			}
			else
			{
				new_energy_a = old_energy_a;
				new_alpha_tb = old_alpha_tb;
				wait_turn_count++;
			}

			total_iteration_count++;
			if (savedata_yes) {
				if (new_energy_a < current_min_energy) {
					current_min_energy = new_energy_a;
				}

				LARGE_INTEGER current_time;
				QueryPerformanceCounter(&current_time);
				double elapsed_ms = (current_time.QuadPart - last_record_time.QuadPart) * 1000.0 / frequency.QuadPart;

				if (elapsed_ms >= record_interval_ms) {
					double total_time_ms = (current_time.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
					opt_records.push_back(OptimizationRecord(total_time_ms, total_iteration_count, current_min_energy));
					last_record_time = current_time;
				}
			}
			//生成一个随机数列，打乱优化顺序
			vector<int> rand_seq;
			vector<bool>rec(nhc, false);
			if (has_circle || rough_turncount != 1) rand_seq = GenerateRandSeq(nhc);
			int continue_yes_count = 0;
			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle || rough_turncount != 1)i = rand_seq[rand_i];
				else i = rand_i;
				rec[i] = true;

				//bool continue_yes = false;
				//vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
				//for (int j = 0; j < nb_node.size(); j++)
				//{
				//	int child = nb_node[j].GetDesSeq();
				//	string bond_sym = nb_node[j].GetBondSymbol();
				//	if (child < nhc && rec[child] && IsSpecialBond(bond_sym))
				//	{
				//		new_alpha_tb[i] = old_alpha_tb[i];
				//		continue_yes = true;
				//		continue_yes_count++;
				//		break;
				//	}
				//}
				//if (continue_yes) continue;

				double now_asep = asep_tb[i];
				double min_energy = INFINITY;
				double min_angle = 0, now_angle = new_alpha_tb[i];
				double angle_addup = 0;
				do {
					new_alpha_tb[i] = now_angle;
					int limitpos = (has_circle) ? -1 : i;
					double now_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, new_ab_opt, all_sp_tb, esp, limitpos);
					if (now_energy < min_energy)
					{
						min_energy = now_energy;
						min_angle = now_angle;
					}

					angle_addup += now_asep;
					now_angle = fmod(now_angle + angle_addup, PI_DOUBLE);

				} while (angle_addup < PI_DOUBLE);
				new_alpha_tb[i] = min_angle;
			}
			if (continue_yes_count == nhc - 1)
			{
				cout << "由于全为特殊键，分子无需优化。" << endl;
				return old_alpha_tb;
			}
			new_energy_a = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);
			if (new_energy_a < rough_min_energy)
			{
				rough_min_energy = new_energy_a;
				rough_min_alpha_tb = new_alpha_tb;
			}
		}
		for (int i = 0; i < nhc; i++)
		{
			asep_tb[i] /= 2;
			if (asep_tb[i] < acc_angle) asep_tb[i] = acc_angle;
		}
		if (turn1count < max_rough_turn && !break_yes_a)
		{
			turn1count++;
			//cout << "粗略优化进行中，第 " << turn1count << " 轮..." << endl;
			goto AlphaOpt1;
		}
		rough_turncount++;
	}
	if (has_circle || nhc > 10)
	{
		if (rough_turncount < max_rough_big_turn)
		{
			if (new_energy_a < rough_min_energy)
			{
				rough_min_energy = new_energy_a;
				rough_min_alpha_tb = new_alpha_tb;
			}
			rough_turncount++;
			asep_tb = rec_asep_tb;
			turn1count = 0;
			//cout << "粗略优化进行中，第 " << turn1count << " 轮..." << endl;
			goto AlphaOpt1;
		}
	}
	else
	{
		rough_min_energy = new_energy_a;
		rough_min_alpha_tb = new_alpha_tb;
	}
	// --- alpha 第二次优化 精细 ---


	new_alpha_tb = rough_min_alpha_tb;



	double old_energy = INFINITY;

	double new_energy = CalSumEnergyByXYZ(detail_print, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

	//cout << "\nE: " << new_energy << endl;

	vector<double> final_alpha_tb_a = rough_min_alpha_tb;
	double final_energy_a = rough_min_energy;


	old_alpha_tb = new_alpha_tb;

	asep_tb = rec_asep_tb;
	//for (int i = 0; i < nhc; i++)
	//{
	//	asep_tb[i] = min(asep_tb[i],PI_HALF/3);
	//	if (asep_tb[i] < acc_angle) asep_tb[i] = acc_angle;
	//}

	if (detail_print) cout << "粗略优化" << total_iteration_count << "轮完成，进入精细优化阶段..." << endl;

	vector<double> acc_asep_tb(nhc, acc_angle);

	wait_turn_count = 0;

	//return new_alpha_tb;

AlphaOpt2:
	{
		do
		{
			//cout << turn2count << endl;
			if (wait_turn_count >= wait_turn_max)
			{
				new_energy = old_energy;
				new_alpha_tb = old_alpha_tb;
				if (detail_print) cout << "等待超过最大轮次，结束精细优化。" << endl;
				break;
			}
			for (int i = 0; i < nhc; i++)
			{
				new_alpha_tb[i] = fmod(new_alpha_tb[i] + PI_DOUBLE, PI_DOUBLE);
			}
			if (new_energy < old_energy)
			{
				old_alpha_tb = new_alpha_tb;
				old_energy = new_energy;
				wait_turn_count = 0;
			}
			else
			{
				new_energy = old_energy;
				new_alpha_tb = old_alpha_tb;
				wait_turn_count++;
			}

			total_iteration_count++;

			if (savedata_yes) {
				if (new_energy < current_min_energy) {
					current_min_energy = new_energy;
				}

				LARGE_INTEGER current_time;
				QueryPerformanceCounter(&current_time);
				double elapsed_ms = (current_time.QuadPart - last_record_time.QuadPart) * 1000.0 / frequency.QuadPart;

				if (elapsed_ms >= record_interval_ms) {
					double total_time_ms = (current_time.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
					opt_records.push_back(OptimizationRecord(total_time_ms, total_iteration_count, current_min_energy));
					last_record_time = current_time;
				}
			}

			vector<int> rand_seq;
			vector<bool> rec(nhc, false);
			if (has_circle) rand_seq = GenerateRandSeq(nhc);

			//int continue_yes_count = 0;
			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle)i = rand_seq[rand_i];
				else i = rand_i;

				//bool continue_yes = false;
				//vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
				//for (int j = 0; j < nb_node.size(); j++)
				//{
				//	int child = nb_node[j].GetDesSeq();
				//	string bond_sym = nb_node[j].GetBondSymbol();
				//	if (child < nhc && rec[child] && IsSpecialBond(bond_sym))
				//	{
				//		continue_yes = true;
				//		new_alpha_tb[i] = old_alpha_tb[i];
				//		//continue_yes_count++;
				//		break;
				//	}
				//}
				//if (continue_yes) continue;

				now_ab_opt.alpha_tb = new_alpha_tb;
				V1_DTB& now_alpha_tb = now_ab_opt.alpha_tb;
				double left_angle = new_alpha_tb[i] - asep_tb[i];
				double right_angle = new_alpha_tb[i] + asep_tb[i];

				now_alpha_tb[i] = left_angle;


				double left_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_ab_opt, all_sp_tb, esp);

				now_alpha_tb[i] = right_angle;

				double right_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_ab_opt, all_sp_tb, esp);

				double mid_energy = INFINITY;

				int turn3count = 0;
				while (!JudgeStop<double>({ left_energy }, { right_energy }, { acc_energy }) &&
					!JudgeStop<double>({ left_angle }, { right_angle }, { acc_angle }))
				{
					turn3count++;
					//cout <<"T3: " << turn3count << endl;
					//cout << fixed << setprecision(4) << left_angle - right_angle << "\t" << acc_angle << endl;
					//cout << left_energy - right_energy << "\t" << acc_energy << endl;
					double sep_angle = (right_angle - left_angle) * OPT_RATIO;
					double mid_left_angle = left_angle + sep_angle;
					double mid_right_angle = right_angle - sep_angle;

					now_alpha_tb[i] = mid_left_angle;
					double mid_left_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_ab_opt, all_sp_tb, esp);

					now_alpha_tb[i] = mid_right_angle;
					double mid_right_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_ab_opt, all_sp_tb, esp);

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
			//if (continue_yes_count == nhc)
			//{
			//	if (detail_print) cout << "本轮所有原子均因特殊键连接跳过优化，结束精细优化。" << endl;
			//	break;
			//}
			turn2count++;
			bool print_yes = false;
			if (detail_print)print_yes = (turn2count % print_turn == 0);
			if (turn2count % angle_half_turn == 0)
			{
				for (int i = 0; i < nhc; i++)
				{
					asep_tb[i] /= 2;
					if (asep_tb[i] < acc_angle) asep_tb[i] = acc_angle;
				}

			}
			if (print_yes)
			{
				cout << "\n精细优化进行中，第 " << turn2count << " 轮...\n" << endl;
			}
			new_energy = CalSumEnergyByXYZ(print_yes, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);
			if (print_yes)
			{
				cout << "\nEnergySep: " << new_energy - old_energy << endl;
				//print alphasep
				cout << "Alpha: ";
				for (int i = 0; i < nhc; i++)
				{
					cout << new_alpha_tb[i] << "\t";
				}
				cout << endl;

				cout << "AlphaSepReal ";
				for (int i = 0; i < nhc; i++)
				{
					double nowprintsep = abs(new_alpha_tb[i] - old_alpha_tb[i]);
					cout << nowprintsep << "\t";
				}
				cout << endl;

				cout << "AlphaSep ";
				for (int i = 0; i < nhc; i++)
				{
					cout << asep_tb[i] << "\t";
				}
				cout << endl;


			}

		} while (!JudgeStop<double>({ new_energy }, { old_energy }, { acc_energy }) ||
			!JudgeStop(new_alpha_tb, old_alpha_tb, acc_asep_tb));
	}

	if (final_energy_a < new_energy)
	{
		new_alpha_tb = final_alpha_tb_a;
		new_energy = final_energy_a;
	}


	if (detail_print)cout << "精细优化完成，共进行 " << turn2count << " 轮。" << endl;

	//vector<Vec3> b_xyz_tb = CalXYZ(new_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);


	//------------------------------**计时结束**------------------------------
	QueryPerformanceCounter(&stop);
	double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / 1000.0 / frequency.QuadPart;
	cout << "\n(^-^)Alpha优化用时： " << duration << " s" << endl << endl;

	alpha_result = { new_energy,duration };

	if (savedata_yes && !opt_records.empty()) {
		// 添加最终记录点
		double final_time_ms = (stop.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
		opt_records.push_back(OptimizationRecord(final_time_ms, total_iteration_count, current_min_energy));

		// 生成文件名（使用当前时间戳）
		string record_filename = OPT_TIME_FOLDER + now_smiles + "_"
			+ GetCurrentTimeString() + ".csv";

		// 保存到CSV文件
		WriteTable(record_filename, "Time_ms,Iteration,MinEnergy", opt_records);
		cout << "优化记录已保存至: " << record_filename << endl;
	}

	return new_alpha_tb;



	return new_alpha_tb;
}

V2_DTB BetaOpt(bool has_circle, const V1_DTB& alpha_tb, OptRecVal& beta_result, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp, bool detail_print)
{

	LARGE_INTEGER frequency;
	QueryPerformanceFrequency(&frequency);
	LARGE_INTEGER start, stop;
	QueryPerformanceCounter(&start);
	//------------------------------**计时开始**------------------------------

	int nhc = esp.short_adj_list.size();

	int turn1count = 0;
	int rough_turncount = 0;
	int turn2count = 0;
	double bsep = PI_HALF / 24;
	double big_bsep = PI_HALF / 6;
	double minval = -PI / 3;
	double maxval = PI_HALF / 6;

	ABOpt old_ab_opt(nhc);
	ABOpt new_ab_opt(nhc);
	ABOpt now_ab_opt(nhc);
	old_ab_opt.alpha_tb = alpha_tb;
	new_ab_opt.alpha_tb = alpha_tb;
	now_ab_opt.alpha_tb = alpha_tb;


	old_ab_opt.beta_tb = V2_DTB(nhc, V1_DTB(MAX_ADJ_NODE_SIZE, maxval));
	new_ab_opt.beta_tb = V2_DTB(nhc, V1_DTB(MAX_ADJ_NODE_SIZE, 0));
	V2_DTB& old_beta_tb = old_ab_opt.beta_tb;
	V2_DTB& new_beta_tb = new_ab_opt.beta_tb;
	//new_beta_tb[0][0] = 0;
	//old_beta_tb[0][0] = 0;

	//Beta优化不做特殊键判断
	//int cyc = 0;
	//for (int i = 0; i < nhc; i++)
	//{
	//	//bool continue_yes = false;
	//	vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
	//	for (int j = 0; j < nb_node.size(); j++)
	//	{
	//		int child = nb_node[j].GetDesSeq();
	//		string bond_sym = nb_node[j].GetBondSymbol();
	//		if (child < i && IsSpecialBond(bond_sym))
	//		{
	//			//continue_yes = true;
	//			cyc++;
	//			break;
	//		}
	//	}
	//	//if (continue_yes) continue;
	//}
	//if (cyc == nhc - 1) return old_beta_tb;


	V2_DTB rec_bsep_tb = GetBetaSepTable(esp.short_adj_list, bsep);
	V2_DTB bsep_tb = rec_bsep_tb;

	//PrintCmdSepTitle("sp_tb");
	//PrintCommonVector(asep_tb);
	NeedCal beta_need_cal;

	if (!has_circle)
	{
		beta_need_cal = NeedCal(false, true, true, true, true, true);
	}

	double start_energy = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

	const double acc_energy = 1e-2;
	const double acc_angle = 1e-3;

	const int angle_half_turn = min(nhc * 5, 20);
	const int print_turn = 20;

	double old_energy_a = INFINITY;
	double new_energy_a = start_energy;

	const int wait_turn_max = min(nhc * 2, 20);
	int wait_turn_count = 0;
	bool break_yes_a = false;

	const int max_rough_turn = 10;
	const int max_rough_big_turn = nhc * 200;

	double rough_min_energy = INFINITY;
	V2_DTB rough_min_beta_tb = old_beta_tb;
	bool rough_big_break_yes = false;
BetaOpt1:
	{
		// --- alpha 第一次优化 粗略 ---
		while (!JudgeStop2(new_beta_tb, old_beta_tb, bsep_tb))
		{
			if (wait_turn_count >= wait_turn_max)
			{
				if (detail_print) cout << "等待超过最大轮次，结束粗略优化。" << endl;
				new_energy_a = old_energy_a;
				new_beta_tb = old_beta_tb;

				break_yes_a = true;
				break;
			}
			if (new_energy_a <= old_energy_a)
			{
				//for (int i = 0; i < nhc; i++)
				//{
				//	new_alpha_tb[i] = fmod((new_alpha_tb[i] + old_alpha_tb[i]) / 2,PI_DOUBLE);
				//}
				old_beta_tb = new_beta_tb;
				old_energy_a = new_energy_a;
				wait_turn_count = 0;
			}
			else
			{
				new_energy_a = old_energy_a;
				new_beta_tb = old_beta_tb;
				wait_turn_count++;
			}
			//生成一个随机数列，打乱优化顺序
			vector<int> rand_seq;
			vector<bool>rec(nhc, false);
			if (has_circle) rand_seq = GenerateRandSeq(nhc);
			int continue_yes_count = 0;
			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle)i = rand_seq[rand_i];
				else i = rand_i;
				rec[i] = true;

				//bool continue_yes = false;
				vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
				//for (int j = 0; j < nb_node.size(); j++)
				//{
				//	int child = nb_node[j].GetDesSeq();
				//	string bond_sym = nb_node[j].GetBondSymbol();
				//	if (child < nhc && rec[child] && IsSpecialBond(bond_sym))
				//	{
				//		new_alpha_tb[i] = old_alpha_tb[i];
				//		continue_yes = true;
				//		continue_yes_count++;
				//		break;
				//	}
				//}
				//if (continue_yes) continue;
				for (int j = 0; j < nb_node.size(); j++)
				{
					int child = nb_node[j].GetDesSeq();
					if (child >= nhc) break;
					if (rec[child]) continue;

					rec[child] = true;
					double now_asep = bsep_tb[i][j];
					double min_energy = INFINITY;
					double min_angle = minval, now_angle = minval;
					double angle_addup = 0;
					do {
						new_beta_tb[i][j] = now_angle;
						int limitpos = (has_circle) ? -1 : i;
						double now_energy = CalSumEnergyByXYZ(NO, xyz_tb, beta_need_cal, new_ab_opt, all_sp_tb, esp);
						if (now_energy < min_energy)
						{
							min_energy = now_energy;
							min_angle = now_angle;
						}

						angle_addup += now_asep;
						now_angle += angle_addup;
						if (now_angle > maxval) now_angle = maxval;

					} while (angle_addup < maxval);
					new_beta_tb[i][j] = min_angle;
				}
			}
			if (continue_yes_count == nhc - 1)
			{
				cout << "由于全为特殊键，分子无需优化。" << endl;
				return old_beta_tb;
			}
			new_energy_a = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);
			if (new_energy_a < rough_min_energy)
			{
				rough_min_energy = new_energy_a;
				rough_min_beta_tb = new_beta_tb;
			}
		}
		for (int i = 0; i < nhc; i++)
		{
			for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
			{
				if (bsep_tb[i][j] != 0)
				{
					bsep_tb[i][j] /= 2;
					if (bsep_tb[i][j] < acc_angle) bsep_tb[i][j] = acc_angle;
				}

			}

		}
		if (turn1count < max_rough_turn && !break_yes_a)
		{
			turn1count++;
			//cout << "粗略优化进行中，第 " << turn1count << " 轮..." << endl;
			goto BetaOpt1;
		}
		rough_turncount++;
	}
	if (has_circle)
	{
		if (rough_turncount < max_rough_big_turn)
		{
			if (new_energy_a < rough_min_energy)
			{
				rough_min_energy = new_energy_a;
				rough_min_beta_tb = new_beta_tb;
			}
			rough_turncount++;
			bsep_tb = rec_bsep_tb;
			turn1count = 0;
			//cout << "粗略优化进行中，第 " << turn1count << " 轮..." << endl;
			goto BetaOpt1;
		}
	}
	else
	{
		rough_min_energy = new_energy_a;
		rough_min_beta_tb = new_beta_tb;
	}
	// ---beta 第二次优化 精细 ---


	new_beta_tb = rough_min_beta_tb;



	double old_energy = INFINITY;

	double new_energy = CalSumEnergyByXYZ(detail_print, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);
	//return rough_min_beta_tb;


	//cout << "\nE: " << new_energy << endl;

	V2_DTB final_beta_tb_a = rough_min_beta_tb;
	double final_energy_a = rough_min_energy;


	old_beta_tb = new_beta_tb;

	bsep_tb = GetBetaSepTable(esp.short_adj_list, big_bsep);
	//for (int i = 0; i < nhc; i++)
	//{
	//	asep_tb[i] = min(asep_tb[i],PI_HALF/3);
	//	if (asep_tb[i] < acc_angle) asep_tb[i] = acc_angle;
	//}

	if (detail_print) cout << "粗略优化" << rough_turncount << "轮完成，进入精细优化阶段..." << endl;

	V2_DTB acc_bsep_tb(nhc, V1_DTB(MAX_ADJ_NODE_SIZE, acc_angle));

	wait_turn_count = 0;

	//return new_alpha_tb;


BetaOpt2:
	{
		do
		{
			//cout << turn2count << endl;
			if (wait_turn_count >= wait_turn_max)
			{
				new_energy = old_energy;
				new_beta_tb = old_beta_tb;
				if (detail_print) cout << "等待超过最大轮次，结束精细优化。" << endl;
				break;
			}
			for (int i = 0; i < nhc; i++)
			{
				for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
				{
					new_beta_tb[i][j] = fmod(new_beta_tb[i][j] + PI_DOUBLE, PI_DOUBLE);
				}
			}
			if (new_energy < old_energy)
			{
				old_beta_tb = new_beta_tb;
				old_energy = new_energy;
				wait_turn_count = 0;
			}
			else
			{
				new_energy = old_energy;
				new_beta_tb = old_beta_tb;
				wait_turn_count++;
			}
			vector<int> rand_seq;
			vector<bool> rec(nhc, false);
			if (has_circle) rand_seq = GenerateRandSeq(nhc);

			//int continue_yes_count = 0;
			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle)i = rand_seq[rand_i];
				else i = rand_i;
				rec[i] = true;

				//bool continue_yes = false;

				vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
				//for (int j = 0; j < nb_node.size(); j++)
				//{
				//	int child = nb_node[j].GetDesSeq();
				//	string bond_sym = nb_node[j].GetBondSymbol();
				//	if (child < nhc && rec[child] && IsSpecialBond(bond_sym))
				//	{
				//		continue_yes = true;
				//		new_alpha_tb[i] = old_alpha_tb[i];
				//		//continue_yes_count++;
				//		break;
				//	}
				//}
				//if (continue_yes) continue;
				for (int j = 0; j < nb_node.size(); j++)
				{
					int child = nb_node[j].GetDesSeq();
					if (child >= nhc) break;
					if (rec[child]) continue;
					rec[child] = true;
					now_ab_opt.beta_tb = new_beta_tb;
					V2_DTB& now_beta_tb = now_ab_opt.beta_tb;
					double left_angle = max(new_beta_tb[i][j] - bsep_tb[i][j], minval);
					double right_angle = min(new_beta_tb[i][j] + bsep_tb[i][j], maxval);
					now_beta_tb[i][j] = left_angle;

					//now_ab_opt.beta_tb = new_beta_tb;
					//V2_DTB& now_beta_tb = now_ab_opt.beta_tb;

					//double left_angle = new_beta_tb[i] - bsep_tb[i];
					//double right_angle = new_beta_tb[i] + bsep_tb[i];

					//now_beta_tb[i] = left_angle;


					double left_energy = CalSumEnergyByXYZ(NO, xyz_tb, beta_need_cal, now_ab_opt, all_sp_tb, esp);

					now_beta_tb[i][j] = right_angle;

					double right_energy = CalSumEnergyByXYZ(NO, xyz_tb, beta_need_cal, now_ab_opt, all_sp_tb, esp);

					double mid_energy = INFINITY;

					int turn3count = 0;
					while (!JudgeStop<double>({ left_energy }, { right_energy }, { acc_energy }) &&
						!JudgeStop<double>({ left_angle }, { right_angle }, { acc_angle }))
					{
						turn3count++;
						//cout <<"T3: " << turn3count << endl;
						//cout << fixed << setprecision(4) << left_angle - right_angle << "\t" << acc_angle << endl;
						//cout << left_energy - right_energy << "\t" << acc_energy << endl;
						double sep_angle = (right_angle - left_angle) * OPT_RATIO;
						double mid_left_angle = left_angle + sep_angle;
						double mid_right_angle = right_angle - sep_angle;

						now_beta_tb[i][j] = mid_left_angle;
						double mid_left_energy = CalSumEnergyByXYZ(NO, xyz_tb, beta_need_cal, now_ab_opt, all_sp_tb, esp);

						now_beta_tb[i][j] = mid_right_angle;
						double mid_right_energy = CalSumEnergyByXYZ(NO, xyz_tb, beta_need_cal, now_ab_opt, all_sp_tb, esp);

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
					new_beta_tb[i][j] = (left_angle + right_angle) / 2;
				}
			}
			//if (continue_yes_count == nhc)
			//{
			//	if (detail_print) cout << "本轮所有原子均因特殊键连接跳过优化，结束精细优化。" << endl;
			//	break;
			//}
			turn2count++;
			bool print_yes = false;
			if (detail_print)print_yes = (turn2count % print_turn == 0);
			if (turn2count % angle_half_turn == 0)
			{
				for (int i = 0; i < nhc; i++)
				{
					for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
					{
						if (bsep_tb[i][j] != 0)
						{
							bsep_tb[i][j] /= 2;
							if (bsep_tb[i][j] < acc_angle) bsep_tb[i][j] = acc_angle;
						}
					}
				}

			}
			if (print_yes)
			{
				cout << "\n精细优化进行中，第 " << turn2count << " 轮...\n" << endl;
			}
			new_energy = CalSumEnergyByXYZ(print_yes, xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);
			if (print_yes)
			{
				cout << "\nEnergySep: " << new_energy - old_energy << endl;
				//print alphasep
				cout << "Beta: ";
				for (int i = 0; i < nhc; i++)
				{
					for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
					{
						cout << new_beta_tb[i][j] << "\t";
					}
					cout << endl;
				}
				cout << endl;

				cout << "BetaSepReal ";
				for (int i = 0; i < nhc; i++)
				{
					for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
					{
						double nowprintsep = abs(new_beta_tb[i][j] - old_beta_tb[i][j]);
						cout << nowprintsep << "\t";
					}
					cout << endl;
				}
				cout << endl;

				cout << "BetaSep ";
				for (int i = 0; i < nhc; i++)
				{
					for (int j = 0; j < MAX_ADJ_NODE_SIZE; j++)
					{
						cout << bsep_tb[i][j] << "\t";
					}
					cout << endl;
				}
				cout << endl;


			}

		} while (!JudgeStop<double>({ new_energy }, { old_energy }, { acc_energy }) ||
			!JudgeStop2(new_beta_tb, old_beta_tb, acc_bsep_tb));
	}

	if (final_energy_a < new_energy) {
		new_beta_tb = final_beta_tb_a;
		new_energy = final_energy_a;
	}


	if (detail_print)cout << "精细优化完成，共进行 " << turn2count << " 轮。" << endl;


	//vector<Vec3> b_xyz_tb = CalXYZ(new_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);

	QueryPerformanceCounter(&stop);

	//------------------------------**计时结束**------------------------------

	double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / 1000.0 / frequency.QuadPart;
	cout << "\n(^-^)Beta优化用时： " << duration << " s" << endl << endl;

	beta_result = { new_energy,duration };
	//xyz_tb = new_beta_tb;

	return new_beta_tb;
}

//检查表格中记录的smiles格式对应的json文件是否存在于目录中，若不存在则删除该行记录
void CheckJsonExist(vector<OptRecLine>& optrec_tb)
{
	vector<OptRecLine> valid_records;
	string output_csv_folder = OPT_OUTPUT_FOLDER;
	string json_folder = JSON_OUTPUT_FOLDER;

	WriteTable(output_csv_folder + OPT_COPY_FILENAME, OPT_REC_TB_TITLE, optrec_tb);
	for (const auto& record : optrec_tb) {
		string json_filename = json_folder + "/" + record.GetIndex() + ".json";
		//cout << "检查文件: " << json_filename << endl;
		if (FileExists(json_filename)) {
			valid_records.push_back(record);
		}
		else {
			cout << "记录 " << record.GetIndex() << " 对应的JSON文件不存在，已删除该记录。" << endl;
		}
	}
	optrec_tb = valid_records;
	WriteTable(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb);

}
