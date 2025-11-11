#include "types.h"
#include "Atom.h"
#include "Molecule.h"
#include "Matrix.h"
#include "MMFF94Table.h"
#include "MMFF94Cal.h"
#include "MySocket.h"
#include "AlphaBetaOpt.h"

using namespace std;


int main()
{
	EnergyFundTable eft = ReadEnergySolidParam();
	SmilesFundTable sft = ReadSmilesSolidParam();
	SP_SpTable all_sp_tb;
	const string output_folder = JSON_OUTPUT_FOLDER;
	const string output_csv_folder = OPT_OUTPUT_FOLDER;
	const bool html_print = false;
	const bool smiles_detail_print = false;
	const bool mmff_detail_print = false;
	const bool opt_print =false;

	vector<OptRecLine> optrec_tb;
	if (ReadTableByTitle(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb) == -1) {
		WriteTable(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb);
	}
	CheckJsonExist(optrec_tb);
	//WriteTable(output_csv_folder + OPT_REC_FILENAME, OPT_REC_TB_TITLE, optrec_tb);

	map<string, OptRecVal> optrec_map = ChangeTableToMap<OptRecLine, OptRecVal>(optrec_tb);


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
		// ====================== 处理静态文件请求 ======================

		// 处理GET / 或 /Molecular3D.html 请求
		if (method == "GET" && (path == "/" || path == "/Molecular3D.html")) {
			string htmlPath = "File\\Web\\Molecular3D.html";
			string htmlContent = ReadFileContent(htmlPath);

			if (!htmlContent.empty()) {
				SendHTTPResponse(clientSocket, 200, "text/html; charset=utf-8", htmlContent);
				closesocket(clientSocket);
				cout << "已发送 Molecular3D.html (" << htmlContent.length() << " 字节)\n" << endl;
				continue;
			}
			else {
				string errorResponse = "{\"status\":\"error\",\"message\":\"HTML文件未找到\"}";
				SendHTTPResponse(clientSocket, 404, "application/json", errorResponse);
				closesocket(clientSocket);
				cout << "错误: Molecular3D.html 文件未找到\n" << endl;
				continue;
			}
		}

		// 处理GET /documentation.md 请求
		if (method == "GET" && path == "/documentation.md") {
			string mdPath = "File\\Web\\documentation.md";
			string mdContent = ReadFileContent(mdPath);

			if (!mdContent.empty()) {
				SendHTTPResponse(clientSocket, 200, "text/markdown; charset=utf-8", mdContent);
				closesocket(clientSocket);
				cout << "已发送 documentation.md (" << mdContent.length() << " 字节)\n" << endl;
				continue;
			}
			else {
				string errorResponse = "{\"status\":\"error\",\"message\":\"文档文件未找到\"}";
				SendHTTPResponse(clientSocket, 404, "application/json", errorResponse);
				closesocket(clientSocket);
				cout << "错误: documentation.md 文件未找到\n" << endl;
				continue;
			}
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

			//cout << requestBody << endl;
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

			if (html_print)
			{

				PrintCmdSepTitle("网页信息");
				cout << requestBody << endl;
				PrintCmdSepTitle("简化分子节点表", H2_SEP_WIDTH);
				PrintSpecialVector(sn_tb);
				PrintCmdSepTitle("简化分子邻接表", H3_SEP_WIDTH);
				PrintSpecialVector(adj_list);

			}
			if (smiles_detail_print)
			{
				MoleInfoPrint(c);
			}


			//cout << "原SMILES：\n";
			//cout << c.GetComSmiles() << endl << endl;
			//cout << "唯一SMILES：\n";
			//cout << now_smiles << endl;

			if (!c.HasCircle() || !c.HasSpecialBond())
			{
				//Mole mid(sft.atomtable, now_smiles);
				now_smiles = c.GenerateCanSmilesNoCircle(sft.primetable);
				//cout << "进行转换！" << endl;
			}
			PrintCmdSepTitle("唯一SMILES生成");
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

			//string errorResponse_part = "{\"status\":\"error\",\"message\":\"未知请求路径\"}";
			//SendHTTPResponse(clientSocket, 404, "application/json", errorResponse_part);
			//closesocket(clientSocket);
			//cout << "未知请求: " << method << " " << path << "\n" << endl;
			//continue;

			bool rec_find = optrec_map.find(now_smiles) != optrec_map.end();
			double rec_energy = rec_find ? optrec_map[now_smiles].min_energy : INFINITY;
			double rec_opttime = rec_find ? optrec_map[now_smiles].opt_time : INFINITY;
			bool has_sent = false;
			if (rec_find && rec_opttime >= MAX_CAL_TIME)
			{
				PrintCmdSepTitle("预计算数据发送");
				string sent_json_rec = ReadJsonFile(folderpath + "/" + now_smiles + ".json");
				ResultSent("precal", clientSocket, sent_json_rec, now_smiles, rec_energy);
				has_sent = true;
			}


			//---------------------------**MMFF识别**----------------------------------
			int ac;
			EnergySolidParam esp = GenEnergySolidParam(ac, eft, d, mmff_detail_print);
			int nhc = esp.short_adj_list.size();
			bool has_circle = c.HasCircle();


			XYZ_TB final_xyz_tb;
			double final_energy;
			//---------------------------**Alpha优化**----------------------------------
			ABOpt old_ab_opt(nhc);
			ABOpt new_ab_opt(nhc);

			vector<double> old_alpha_tb(nhc, 0);
			XYZ_TB old_alpha_xyz_tb;

			bool update_optrec_yes = false;

			PrintCmdSepTitle("Alpha优化前");
			double alpha_before_energy = CalSumEnergyByXYZ(YES, old_alpha_xyz_tb, NeedCal(), old_ab_opt, all_sp_tb, esp);
			//WriteTable(GetXYZTbPath(now_smiles, 1), XYZ_TB_TITLE, xyz_tb);

			PrintCmdSepTitle("Alpha优化进度-记时");
			cout << "Alpha优化中……\n";

			XYZ_TB new_alpha_xyz_tb;
			OptRecVal alpha_opt_val;
			//vector<double> new_alpha_tb = AlphaOpt(has_circle, alpha_energy, new_xyz_tb, all_sp_tb, esp, opt_print, now_smiles);
			new_ab_opt.alpha_tb = AlphaOpt(has_circle, new_ab_opt.beta_tb, alpha_opt_val, new_alpha_xyz_tb, all_sp_tb, esp, opt_print);

			cout << "优化完成！\n";

			PrintCmdSepTitle("Alpha优化后");
			//WriteTable(GetXYZTbPath(now_smiles, 2), XYZ_TB_TITLE, alpha_xyz_tb);

			double alpha_energy = CalSumEnergyByXYZ(YES, new_alpha_xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

			if (alpha_before_energy < alpha_energy)
			{
				final_energy = alpha_before_energy;
				final_xyz_tb = old_alpha_xyz_tb;
			}
			else
			{
				final_energy = alpha_energy;
				final_xyz_tb = new_alpha_xyz_tb;
			}
			//string errorResponse = "{\"status\":\"error\",\"message\":\"未知请求路径\"}";
			//SendHTTPResponse(clientSocket, 404, "application/json", errorResponse);
			//closesocket(clientSocket);
			//cout << "未知请求: " << method << " " << path << "\n" << endl;

			//---------------------------**Beta优化**----------------------------------

			XYZ_TB new_beta_xyz_tb;
			OptRecVal beta_opt_val;

			PrintCmdSepTitle("Beta优化进度-记时");
			cout << "优化中……\n";
			new_ab_opt.beta_tb = BetaOpt(has_circle, new_ab_opt.alpha_tb, beta_opt_val, new_beta_xyz_tb, all_sp_tb, esp, opt_print);
			//
			////new_beta_xyz_tb = CalXYZ(new_ab_opt, esp.pro_htb, esp.bs_fast_tb, esp.mnode_tb, esp.short_adj_list, all_sp_tb);
			cout << "优化完成！\n";
			PrintCmdSepTitle("Beta优化后");
			////WriteTable(GetXYZTbPath(now_smiles, 2), XYZ_TB_TITLE, alpha_xyz_tb);

			double beta_energy = CalSumEnergyByXYZ(YES, new_beta_xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

			if (beta_energy < final_energy)
			{
				final_energy = beta_energy;
				final_xyz_tb = new_beta_xyz_tb;
			}
			OptRecVal final_opt_val = { final_energy, alpha_opt_val.opt_time + beta_opt_val.opt_time };

			//cout << "当前最小能量：" << final_energy << " kcal/mol\n";
			//---------------------------**Alpha优化**----------------------------------

			//PrintCmdSepTitle("Alpha优化进度-记时");
			//cout << "Alpha优化中……\n";

			////XYZ_TB new_alpha_xyz_tb;
			////OptRecVal alpha_opt_val;
			////vector<double> new_alpha_tb = AlphaOpt(has_circle, alpha_energy, new_xyz_tb, all_sp_tb, esp, opt_print, now_smiles);
			//new_ab_opt.alpha_tb = AlphaOpt(has_circle, new_ab_opt.beta_tb, alpha_opt_val, new_alpha_xyz_tb, all_sp_tb, esp, opt_print);

			//cout << "优化完成！\n";

			//PrintCmdSepTitle("Alpha优化后");
			////WriteTable(GetXYZTbPath(now_smiles, 2), XYZ_TB_TITLE, alpha_xyz_tb);

			// alpha_energy = CalSumEnergyByXYZ(YES, new_alpha_xyz_tb, NeedCal(), new_ab_opt, all_sp_tb, esp);

			//if (alpha_before_energy < alpha_energy)
			//{
			//	final_energy = alpha_before_energy;
			//	final_xyz_tb = old_alpha_xyz_tb;
			//}
			//else
			//{
			//	final_energy = alpha_energy;
			//	final_xyz_tb = new_alpha_xyz_tb;
			//}

			//---------------------------**JSON格式传输**----------------------------------

			string sent_json;
			//final_energy < rec_energy
			if (final_energy < rec_energy)
			{
				rec_energy = final_energy;
				if (rec_opttime < MAX_REC_TIME)
				{
					optrec_map[now_smiles] = { rec_energy,max(rec_opttime,final_opt_val.opt_time) };
				}
				else
				{
					optrec_map[now_smiles] = final_opt_val;
				}
				vector<SXYZ_3D> sxyz_tb = GetSXYZTb(final_xyz_tb, esp);
				vector<AdjABS_3D> adjline_tb = GetAdjABSTb(esp);
				cout << "\n优化成功！最终能量: " << final_energy << " kcal/mol\n";

				update_optrec_yes = true;
				//PrintCmdSepTitle("JSON数据生成");

				string sxyz_json_part = VectorToJsonPart("SXYZ", sxyz_tb);
				string adjab_json_part = VectorToJsonPart("ADJ_AB", adjline_tb);

				string sent_json_path = folderpath + "/" + now_smiles + ".json";
				sent_json = "{\n" + sxyz_json_part + "\n,\n" + adjab_json_part + "\n}\n";

				WriteJsonFile(sent_json_path, sent_json);


				//cout << sent_json << endl;
			}
			else
			{
				cout << "\n优化未达预期，采用优化前结构。最终能量: " << rec_energy << " kcal/mol\n";
				sent_json = ReadJsonFile(folderpath + "/" + now_smiles + ".json");
				if (rec_opttime >= MAX_REC_TIME)
				{
					rec_opttime = final_opt_val.opt_time;
					optrec_map[now_smiles] = { rec_energy,rec_opttime };
					update_optrec_yes = true;
				}
				else if (final_opt_val.opt_time > rec_opttime)
				{
					optrec_map[now_smiles] = { rec_energy,final_opt_val.opt_time };
					update_optrec_yes = true;
				}
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
			if (!has_sent)
			{
				PrintCmdSepTitle("正常计算结果发送");

				// 构建响应JSON（包含状态和3D数据）
				ResultSent("success", clientSocket, sent_json, now_smiles, rec_energy);

				cout << "\nCanSmiles: " << now_smiles << endl;
			}
			// 关闭客户端连接
			closesocket(clientSocket);
			cout << "\n客户端连接已关闭，等待下一个连接...\n" << endl;
			//储存OptRec记录
			if (update_optrec_yes)
			{
				PrintCmdSepTitle("优化记录更新");
				vector<OptRecLine> optrec_tb = ChangeMapToTable<OptRecLine, OptRecVal>(optrec_map);
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

