#include "types.h"
#include "Atom.h"
#include "Molecule.h"
#include "Matrix.h"
#include "MMFF94Table.h"
#include "MMFF94Cal.h"
#include "MySocket.h"

using namespace std;


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

vector<int> GenerateRandSeq(int size)
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
vector<double> AlphaOptRough(bool has_circle, double& alpha_energy, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp)
{
	vector<double> alpha_tb;

	return alpha_tb;
}
vector<double> AlphaOpt(bool has_circle, double& alpha_energy, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp, bool detail_print = false)
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
	double asep = PI_HALF;
	vector<double> old_alpha_tb(nhc, 0);
	vector<double> new_alpha_tb(nhc, PI);
	new_alpha_tb[0] = 0;


	vector<double> rec_asep_tb = GetAlphaSepTable(esp.short_adj_list);
	vector<double> asep_tb = rec_asep_tb;

	//PrintCmdSepTitle("sp_tb");
	//PrintCommonVector(asep_tb);
	NeedCal alpha_need_cal;

	if (!has_circle)
	{
		alpha_need_cal = NeedCal(false, false, false, false, true, true);
	}

	double start_energy = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_alpha_tb, all_sp_tb, esp);

	const double acc_energy = 1e-2;
	const double acc_angle = 1e-3;

	const int angle_half_turn = nhc * 5;
	const int print_turn = 20;

	double old_energy_a = INFINITY;
	double new_energy_a = start_energy;

	const int wait_turn_max = min(nhc * 2, 30);
	int wait_turn_count = 0;
	bool break_yes_a = false;

	const int max_rough_turn = 10;
	const int max_rough_big_turn = nhc*100;

	double rough_min_energy = INFINITY;
	vector<double> rough_min_alpha_tb = old_alpha_tb;

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
			//生成一个随机数列，打乱优化顺序
			vector<int> rand_seq;
			if (has_circle) rand_seq = GenerateRandSeq(nhc);

			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle)i = rand_seq[rand_i];
				else i = rand_i;

				double now_asep = asep_tb[i];
				double min_energy = INFINITY;
				double min_angle = 0, now_angle = new_alpha_tb[i];
				double angle_addup = 0;
				do {
					new_alpha_tb[i] = now_angle;
					int limitpos = (has_circle) ? -1 : i;
					double now_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, new_alpha_tb, all_sp_tb, esp, limitpos);
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
			new_energy_a = CalSumEnergyByXYZ(NO, xyz_tb, NeedCal(), new_alpha_tb, all_sp_tb, esp);
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
	if (has_circle)
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

	double new_energy = CalSumEnergyByXYZ(detail_print, xyz_tb, NeedCal(), new_alpha_tb, all_sp_tb, esp);

	//cout << "\nE: " << new_energy << endl;

	vector<double> final_alpha_tb_a = new_alpha_tb;
	double final_energy_a = new_energy;


	old_alpha_tb = new_alpha_tb;

	asep_tb = rec_asep_tb;
	//for (int i = 0; i < nhc; i++)
	//{
	//	asep_tb[i] = min(asep_tb[i],PI_HALF/3);
	//	if (asep_tb[i] < acc_angle) asep_tb[i] = acc_angle;
	//}

	if (detail_print) cout << "粗略优化" << turn1count << "轮完成，进入精细优化阶段..." << endl;

	vector<double> acc_asep_tb(nhc, acc_angle);

	wait_turn_count = 0;

AlphaOpt2:
	{
		do
		{
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
			vector<int> rand_seq;
			if (has_circle) rand_seq = GenerateRandSeq(nhc);

			for (int rand_i = 0; rand_i < nhc; rand_i++)
			{
				int i;
				if (has_circle)i = rand_seq[rand_i];
				else i = rand_i;

				vector<double>now_alpha_tb = new_alpha_tb;
				double left_angle = new_alpha_tb[i] - asep_tb[i];
				double right_angle = new_alpha_tb[i] + asep_tb[i];

				now_alpha_tb[i] = left_angle;


				double left_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

				now_alpha_tb[i] = right_angle;

				double right_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

				double mid_energy = INFINITY;

				while (!JudgeStop<double>({ left_energy }, { right_energy }, { acc_energy }) ||
					!JudgeStop<double>({ left_angle }, { right_angle }, { acc_angle }))
				{
					double sep_angle = (right_angle - left_angle) * OPT_RATIO;
					double mid_left_angle = left_angle + sep_angle;
					double mid_right_angle = right_angle - sep_angle;

					now_alpha_tb[i] = mid_left_angle;
					double mid_left_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

					now_alpha_tb[i] = mid_right_angle;
					double mid_right_energy = CalSumEnergyByXYZ(NO, xyz_tb, alpha_need_cal, now_alpha_tb, all_sp_tb, esp);

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
			new_energy = CalSumEnergyByXYZ(print_yes, xyz_tb, NeedCal(), new_alpha_tb, all_sp_tb, esp);
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

	if (final_energy_a < new_energy) new_alpha_tb = final_alpha_tb_a;


	if (detail_print)cout << "精细优化完成，共进行 " << turn2count << " 轮。" << endl;
	alpha_energy = new_energy;
	//vector<Vec3> b_xyz_tb = CalXYZ(new_alpha_tb, pro_htb, bs_fast_tb, long_adj_list, all_sp_tb);

	QueryPerformanceCounter(&stop);

	//------------------------------**计时结束**------------------------------

	double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / 1000.0 / frequency.QuadPart;
	cout << "\n(^-^)Alpha优化用时： " << duration << " s" << endl << endl;

	return new_alpha_tb;
}


// 计划（伪代码）：
// 1. 使用粒子群优化（PSO）搜索 alpha 向量的全局最优解。
// 2. 初始化：根据原子数 nhc 设定粒子数 swarm_size；为每个粒子随机生成位置（角度范围 [0, PI)）和小幅速度；初始化个体最优 pbest 和全局最优 gbest。
// 3. 迭代：在每一代中，按惯性权重 w（线性递减），学习因子 c1, c2 更新速度：
//    v = w*v + c1*r1*(pbest - pos) + c2*r2*(gbest - pos)
//    控制速度幅度并更新位置后规一化到 [0, PI)。
// 4. 每更新一次位置就计算能量（调用 CalSumEnergyByXYZ，print_yes = NO），如果比 pbest 好则更新 pbest；并更新 gbest。
// 5. 终止条件：达到最大迭代次数或全局收敛（能量变化小于阈值若干代）。
// 6. 结束后用 gbest 生成最终三维坐标（CalSumEnergyByXYZ(print_yes = YES)）并返回 gbest，同时通过引用设置 alpha_energy 和 xyz_tb。
// 7. 参数选择：w 逐渐从 0.9 线性降到 0.4，c1=c2=1.5，速度限幅为 PI，最大迭代与耐心轮数根据 nhc 设置。

vector<double> AlphaPsoOpt(double& alpha_energy, XYZ_TB& xyz_tb, const SP_SpTable& all_sp_tb, EnergySolidParam& esp)
{

	LARGE_INTEGER b_frequency;
	QueryPerformanceFrequency(&b_frequency);
	LARGE_INTEGER b_start, b_stop;
	QueryPerformanceCounter(&b_start);
	//------------------------------**计时开始**------------------------------
	int nhc = esp.short_adj_list.size();
	if (nhc <= 0) {
		alpha_energy = INFINITY;
		xyz_tb.clear();
		return vector<double>();
	}

	// PSO 参数
	int swarm_size = max(12, min(60, nhc * 4)); // 粒子数
	int max_iter = 400;                          // 最大迭代次数
	double w_max = 0.9, w_min = 0.4;            // 惯性权重线性衰减
	double c1 = 1.5, c2 = 1.5;                  // 学习因子
	double vel_max = PI_DOUBLE;                 // 速度限幅
	double energy_tol = 1e-6;                   // 能量收敛阈值
	int patience_limit = 40;                    // 无改进耐心代数

	std::mt19937_64 rng(static_cast<unsigned long>(time(nullptr)));
	std::uniform_real_distribution<double> uni_angle(0.0, PI_DOUBLE);
	std::uniform_real_distribution<double> uni_vel(-PI_DOUBLE / 8.0, PI_DOUBLE / 8.0);
	std::uniform_real_distribution<double> uni_01(0.0, 1.0);

	// 粒子群数据结构
	vector<vector<double>> pos(swarm_size, vector<double>(nhc));
	vector<vector<double>> vel(swarm_size, vector<double>(nhc));
	vector<vector<double>> pbest_pos(swarm_size, vector<double>(nhc));
	vector<double> pbest_energy(swarm_size, INFINITY);

	vector<double> gbest_pos(nhc);
	double gbest_energy = INFINITY;

	NeedCal need_cal; // 默认计算所有能量项

	// 初始化粒子
	for (int s = 0; s < swarm_size; ++s) {
		for (int i = 0; i < nhc; ++i) {
			pos[s][i] = uni_angle(rng);
			vel[s][i] = uni_vel(rng);
		}
		// 评估初始能量（不打印）
		XYZ_TB tmp_xyz;
		double e = CalSumEnergyByXYZ(NO, tmp_xyz, need_cal, pos[s], all_sp_tb, esp);
		pbest_pos[s] = pos[s];
		pbest_energy[s] = e;
		if (e < gbest_energy) {
			gbest_energy = e;
			gbest_pos = pos[s];
		}
	}

	int iter = 0;
	int no_improve_count = 0;
	double last_best_energy = gbest_energy;

	while (iter < max_iter && no_improve_count < patience_limit) {
		double w = w_max - (w_max - w_min) * (double(iter) / double(max_iter)); // 线性下降

		for (int s = 0; s < swarm_size; ++s) {
			// 更新速度与位置
			for (int i = 0; i < nhc; ++i) {
				double r1 = uni_01(rng);
				double r2 = uni_01(rng);
				double cognitive = c1 * r1 * (pbest_pos[s][i] - pos[s][i]);
				double social = c2 * r2 * (gbest_pos[i] - pos[s][i]);
				vel[s][i] = w * vel[s][i] + cognitive + social;
				// 限幅
				if (vel[s][i] > vel_max) vel[s][i] = vel_max;
				if (vel[s][i] < -vel_max) vel[s][i] = -vel_max;
				// 更新位置并映射到 [0, PI)
				pos[s][i] += vel[s][i];
				// 规范化到 [0, PI)
				pos[s][i] = fmod(pos[s][i], PI_DOUBLE);
				if (pos[s][i] < 0) pos[s][i] += PI_DOUBLE;
			}

			// 评估能量（不打印）
			XYZ_TB tmp_xyz;
			double e = CalSumEnergyByXYZ(NO, tmp_xyz, need_cal, pos[s], all_sp_tb, esp);

			// 更新个体最优
			if (e < pbest_energy[s]) {
				pbest_energy[s] = e;
				pbest_pos[s] = pos[s];
			}
			// 更新全局最优
			if (e < gbest_energy) {
				gbest_energy = e;
				gbest_pos = pos[s];
			}
		}

		// 早停检测
		if (abs(last_best_energy - gbest_energy) < energy_tol) {
			no_improve_count++;
		}
		else {
			no_improve_count = 0;
			last_best_energy = gbest_energy;
		}
		iter++;
	}

	// 计算并返回最终坐标（打印最终结果）
	XYZ_TB final_xyz;
	double final_e = CalSumEnergyByXYZ(NO, final_xyz, need_cal, gbest_pos, all_sp_tb, esp);
	alpha_energy = final_e;
	xyz_tb = final_xyz;

	//------------------------------**计时结束**------------------------------
	QueryPerformanceCounter(&b_stop);
	double b_duration = (b_stop.QuadPart - b_start.QuadPart) * 1000.0 / 1000.0 / b_frequency.QuadPart;
	std::cout << "\n(^-^)PSO优化用时: " << b_duration << " s" << std::endl;


	return gbest_pos;
}


class OptRec {
private:
	string smiles;
	double min_energy;
public:
	OptRec(string s, double e) :smiles(s), min_energy(e) {}
	OptRec(const vector<string>& info, const vector<int>& titlenum)
	{
		smiles = info[titlenum[0]];
		min_energy = atof(info[titlenum[1]].c_str());
	}
	OptRec(const OptRec& r) :smiles(r.smiles), min_energy(r.min_energy) {}
	string GetIndex()const { return smiles; }
	double GetVal()const { return min_energy; }
	void Print(string sep = "\t")const
	{
		cout << smiles << sep << min_energy << endl;
	}

	friend ostream& operator<<(ostream& os, const OptRec& orc)
	{
		os << orc.smiles << "," << orc.min_energy;
		return os;
	}
};




int main()
{
	EnergyFundTable eft = ReadEnergySolidParam();
	SmilesFundTable sft = ReadSmilesSolidParam();
	SP_SpTable all_sp_tb;
	const string output_folder = JSON_OUTPUT_FOLDER;
	const string output_csv_folder = OPT_OUTPUT_FOLDER;
	const bool smiles_detail_print = false;
	const bool mmff_detail_print = true;

	vector<OptRec> optrec_tb;
	if (ReadTableByTitle(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb) == -1) {
		WriteTable(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb);
	}
	map<string, double> optrec_map = ChangeTableToMap<OptRec, double>(optrec_tb);


	cout << "\n数据初始化已完成(o゜▽゜)o☆\n";

	// 启动Socket服务器
	if (!InitializeWinsock()) {
		cout << "Winsock初始化失败！" << endl;
		return -1;
	}
	if (!StartSocketServer(DEFAULT_PORT)) {
		cout << "Socket服务器启动失败！" << endl;
		return -1;
	}

	int turncount = 0;
	while (serverRunning)
	{
		// 检查是否有客户端连接
		SOCKET clientSocket = accept(serverSocket, NULL, NULL);
		if (clientSocket == INVALID_SOCKET) {
			int error = WSAGetLastError();
			if (error != WSAEWOULDBLOCK) {
				cout << "接受连接失败，错误码: " << error << endl;
			}
			// 没有连接时休眠，降低CPU占用
			Sleep(100);
			continue;
		}

		//cout << "\n========================================" << endl;
		//cout << "客户端已连接！" << endl;
		//cout << "========================================\n" << endl;

		// 接收HTTP请求
		string httpRequest = ReceiveHTTPRequest(clientSocket);

		if (httpRequest.empty()) {
			cout << "未接收到有效请求" << endl;
			closesocket(clientSocket);
			continue;
		}

		// 解析HTTP请求
		string method, path;
		string requestBody = ParseHttpRequest(httpRequest, method, path);

		cout << "收到请求: " << method << " " << path << endl;

		// 处理OPTIONS请求（CORS预检）
		if (method == "OPTIONS") {
			SendHTTPResponse(clientSocket, 200, "text/plain", "OK");
			closesocket(clientSocket);
			cout << "处理OPTIONS预检请求完成\n" << endl;
			continue;
		}

		// 处理GET /status 请求（连接状态检查）
		if (method == "GET" && path == "/status") {
			string statusResponse = "{\"status\":\"running\",\"message\":\"服务器运行中\"}";
			SendHTTPResponse(clientSocket, 200, "application/json", statusResponse);
			closesocket(clientSocket);
			//cout << "状态检查请求处理完成\n" << endl;
			continue;
		}

		// 处理POST /convert 请求（分子转换）
		if (method == "POST" && path == "/convert") {
			turncount++;
			PrintCmdSepTitle("第" + to_string(turncount) + "轮结构式解析");

			cout << "接收到 " << requestBody.length() << " 字节JSON数据" << endl;
	
			// 解析JSON并构建分子
			vector<SimpleMNode> sn_tb;
			vector<AdjLine> adj_list;
			if (!AnalysisJsonFile(requestBody, sn_tb, adj_list))
			{
				cout << "\nJSON数据解析失败，请检查格式！" << endl;
				string errorResponse = "{\"status\":\"error\",\"message\":\"JSON解析失败\"}";
				SendHTTPResponse(clientSocket, 400, "application/json", errorResponse);
				closesocket(clientSocket);
				continue;
			}

			cout << requestBody << endl;
			Mole c;

			string now_smiles;
			if (adj_list.size() == 0)
			{
				c = Mole(sft.atomtable, sn_tb[0].sym);
				now_smiles = sn_tb[0].sym;
			}
			else {
				c = Mole(sft.atomtable, sn_tb, adj_list);
				now_smiles = c.GenerateCanSmiles(sft.primetable);
			}

			if (smiles_detail_print)
			{
				PrintCmdSepTitle("简化分子节点表");
				PrintSpecialVector(sn_tb);
				PrintCmdSepTitle("简化分子邻接表");
				PrintSpecialVector(adj_list);

				PrintCmdSepTitle("分子节点提取");
				c.PrintNodeTable();

				PrintCmdSepTitle("分子节点邻接表");
				c.PrintBondTable();


				PrintCmdSepTitle("秩排序结果");
				c.PrintOriRank();
				c.PrintRank();
				cout << "\n是否含环: " << (c.HasCircle() ? "是" : "否") << endl;
			}

			PrintCmdSepTitle("唯一SMILES生成");
			//cout << "原SMILES：\n";
			//cout << c.GetComSmiles() << endl << endl;
			cout << "唯一SMILES：\n";
			cout << now_smiles << endl;

			Mole d(sft.atomtable, now_smiles);

			/*
			PrintCmdSepTitle("分子节点提取");
			d.PrintNodeTable();
			PrintCmdSepTitle("分子节点邻接表");
			d.PrintBondTable();
			*/
			const string folderpath = output_folder;

			//continue;

			bool rec_find = optrec_map.find(now_smiles) != optrec_map.end();
			double rec_energy = rec_find ? optrec_map[now_smiles] : INFINITY;

			if (rec_find)
			{
				string sent_json_rec = ReadJsonFile(folderpath + "/" + now_smiles + ".json");
				string responseJson_rec = "{\n\"status\":\"precal\",\n\"data\":" + sent_json_rec + "\n}";
				cout << "正在发送预计算3D结构数据到网页..." << endl;
				if (SendHTTPResponse(clientSocket, 200, "application/json", responseJson_rec)) {
					cout << "3D结构数据发送成功！" << endl;
					cout << "数据大小: " << responseJson_rec.length() << " 字节" << endl;
				}
				else {
					cout << "3D结构数据发送失败！" << endl;
				}
			}
			

			//---------------------------**MMFF识别**----------------------------------
			int ac;
			EnergySolidParam esp = GenEnergySolidParam(ac, eft, d, smiles_detail_print);
			int nhc = esp.short_adj_list.size();
			bool has_circle = c.HasCircle();

			//---------------------------**Alpha优化**----------------------------------
			vector<double> old_alpha_tb(nhc, 0);
			vector<Vec3> old_xyz_tb;

			bool alpha_successful = false;

			PrintCmdSepTitle("Alpha优化前");
			double before_energy = CalSumEnergyByXYZ(YES, old_xyz_tb, NeedCal(), old_alpha_tb, all_sp_tb, esp);
			//WriteTable(GetXYZTbPath(now_smiles, 1), XYZ_TB_TITLE, xyz_tb);

			PrintCmdSepTitle("Alpha优化进度-记时");
			cout << "优化中……\n";

			XYZ_TB new_xyz_tb;
			double alpha_energy;
			vector<double> new_alpha_tb = AlphaOpt(has_circle, alpha_energy, new_xyz_tb, all_sp_tb, esp, mmff_detail_print);

			cout << "优化完成！\n";

			PrintCmdSepTitle("Alpha优化后");
			//WriteTable(GetXYZTbPath(now_smiles, 2), XYZ_TB_TITLE, alpha_xyz_tb);

			double final_energy = CalSumEnergyByXYZ(YES, new_xyz_tb, NeedCal(), new_alpha_tb, all_sp_tb, esp);

			if (before_energy < final_energy)
			{
				final_energy = before_energy;
				new_xyz_tb = old_xyz_tb;
			}



			//---------------------------**JSON格式传输**----------------------------------

			string sent_json;
			if (final_energy < rec_energy)
			{
				vector<SXYZ_3D> sxyz_tb = GetSXYZTb(new_xyz_tb, esp);
				vector<AdjABS_3D> adjline_tb = GetAdjABSTb(esp);
				cout << "\n优化成功！最终能量: " << final_energy << " kcal/mol\n";

				alpha_successful = true;
				PrintCmdSepTitle("JSON数据生成");

				string sxyz_json_part = VectorToJsonPart("SXYZ", sxyz_tb);
				string adjab_json_part = VectorToJsonPart("ADJ_AB", adjline_tb);

				string sent_json_path = folderpath + "/" + now_smiles + ".json";
				sent_json = "{\n" + sxyz_json_part + "\n,\n" + adjab_json_part + "\n}\n";

				WriteJsonFile(sent_json_path, sent_json);

				optrec_map[now_smiles] = final_energy;
			}
			else
			{
				cout << "\n优化未达预期，采用优化前结构。最终能量: " << rec_energy << " kcal/mol\n";
				sent_json = ReadJsonFile(folderpath + "/" + now_smiles + ".json");
			}

			//----------------------------**文件导出**-----------------------------------

			//cout << "文件 3d.json 已储存！\n";


			//PrintCmdSepTitle("文件夹" + now_smiles + "创建");
			//if (!CreateFolder(folderpath)) {
			//	cout << "Warning: 文件夹创建失败！可能已存在。\n";
			//}
			//else
			//{
			//	cout << "文件夹创建成功！\n";
			//}

			//PrintCmdSepTitle("CSV文件导出");

			//string sxyz_tb_path = folderpath + "/SXYZ.csv";
			//WriteTable(sxyz_tb_path, SXYZ_TB_TITLE, sxyz_tb);

			//string adjab_filename = folderpath + "/ADJ_AB.csv";
			//WriteTable(adjab_filename, ADJ_AB_TB_TITLE, adjline_tb);

			//cout << "文件 SXYZ.csv 和 ADJ_AB.csv 已储存！\n";

			//---------------------------**HTTP响应发送**----------------------------------
			if (!rec_find)
			{
				PrintCmdSepTitle("HTTP响应发送");

				// 构建响应JSON（包含状态和3D数据）
				string responseJson = "{\n\"status\":\"success\",\n\"data\":" + sent_json + "\n}";

				cout << "正在发送3D结构数据到网页..." << endl;
				if (SendHTTPResponse(clientSocket, 200, "application/json", responseJson)) {
					cout << "3D结构数据发送成功！" << endl;
					cout << "数据大小: " << responseJson.length() << " 字节" << endl;
				}
				else {
					cout << "3D结构数据发送失败！" << endl;
				}

				cout << "\nCanSmiles: " << now_smiles << endl;

				// 关闭客户端连接
				closesocket(clientSocket);
				cout << "\n客户端连接已关闭，等待下一个连接...\n" << endl;
			}

			//储存OptRec记录
			if (alpha_successful)
			{
				PrintCmdSepTitle("优化记录更新");
				vector<OptRec> optrec_tb = ChangeMapToTable<OptRec, double>(optrec_map);
				WriteTable(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb);
			}

		}
		else {
			// 未知请求
			string errorResponse = "{\"status\":\"error\",\"message\":\"未知请求路径\"}";
			SendHTTPResponse(clientSocket, 404, "application/json", errorResponse);
			closesocket(clientSocket);
			cout << "未知请求: " << method << " " << path << "\n" << endl;
		}
	}

	// 清理
	StopSocketServer();
	return 0;
}

