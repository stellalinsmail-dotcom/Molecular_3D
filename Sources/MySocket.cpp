#include "MySocket.h"


SOCKET serverSocket = INVALID_SOCKET;
bool serverRunning = false;

// 初始化Winsock
bool InitializeWinsock() {
	WSADATA wsaData;
	int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (iResult != 0) {
		cout << "WSAStartup失败，错误码: " << iResult << endl;
		return false;
	}
	return true;
}

// 启动Socket服务器 (使用HTTP协议)
bool StartSocketServer(const char* port) {
	struct addrinfo* result = NULL, hints;
	ZeroMemory(&hints, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;
	hints.ai_flags = AI_PASSIVE;

	int iResult = getaddrinfo(NULL, port, &hints, &result);
	if (iResult != 0) {
		cout << "getaddrinfo失败，错误码: " << iResult << endl;
		WSACleanup();
		return false;
	}

	serverSocket = socket(result->ai_family, result->ai_socktype, result->ai_protocol);
	if (serverSocket == INVALID_SOCKET) {
		cout << "创建socket失败，错误码: " << WSAGetLastError() << endl;
		freeaddrinfo(result);
		WSACleanup();
		return false;
	}

	// 设置SO_REUSEADDR选项，避免端口占用问题
	int opt = 1;
	setsockopt(serverSocket, SOL_SOCKET, SO_REUSEADDR, (const char*)&opt, sizeof(opt));

	iResult = ::bind(serverSocket, result->ai_addr, (int)result->ai_addrlen);
	if (iResult == SOCKET_ERROR) {
		cout << "绑定端口失败，错误码: " << WSAGetLastError() << endl;
		cout << "请检查端口" << port << "是否被其他程序占用" << endl;
		freeaddrinfo(result);
		closesocket(serverSocket);
		WSACleanup();
		return false;
	}

	freeaddrinfo(result);

	iResult = listen(serverSocket, SOMAXCONN);
	if (iResult == SOCKET_ERROR) {
		cout << "监听失败，错误码: " << WSAGetLastError() << endl;
		closesocket(serverSocket);
		WSACleanup();
		return false;
	}

	// 设置为非阻塞模式
	u_long mode = 1;
	ioctlsocket(serverSocket, FIONBIO, &mode);

	serverRunning = true;
	cout << "\n=========================================================" << endl;
	cout << "Socket服务器已启动（HTTP模式）" << endl;
	cout << "监听端口: " << port << endl;
	cout << "网页地址: http://localhost:" << port <<"/Molecular3D.html" << endl;
	cout << "按下Ctrl同时点击上方网页地址，即可打开网页^-^"<<endl;
	cout << "=========================================================\n" << endl;
	cout << "等待网页连接中...\n" << endl;
	return true;
}

// 解析HTTP请求头
string ParseHttpRequest(const string& request, string& method, string& path) {
	size_t methodEnd = request.find(' ');
	if (methodEnd == string::npos) return "";

	method = request.substr(0, methodEnd);

	size_t pathStart = methodEnd + 1;
	size_t pathEnd = request.find(' ', pathStart);
	if (pathEnd == string::npos) return "";

	path = request.substr(pathStart, pathEnd - pathStart);

	// 找到请求体（在两个\r\n\r\n之后）
	size_t bodyStart = request.find("\r\n\r\n");
	if (bodyStart == string::npos) return "";

	return request.substr(bodyStart + 4);
}

// 发送HTTP响应
bool SendHTTPResponse(SOCKET clientSocket, int statusCode, const string& contentType, const string& body) {
	string statusText;
	switch (statusCode) {
	case 200: statusText = "OK"; break;
	case 400: statusText = "Bad Request"; break;
	case 404: statusText = "Not Found"; break;
	case 500: statusText = "Internal Server Error"; break;
	default: statusText = "Unknown"; break;
	}

	stringstream response;
	response << "HTTP/1.1 " << statusCode << " " << statusText << "\r\n";
	response << "Content-Type: " << contentType << "\r\n";
	response << "Content-Length: " << body.length() << "\r\n";
	response << "Access-Control-Allow-Origin: *\r\n";
	response << "Access-Control-Allow-Methods: GET, POST, OPTIONS\r\n";
	response << "Access-Control-Allow-Headers: Content-Type\r\n";
	response << "Connection: close\r\n";
	response << "\r\n";
	response << body;

	string responseStr = response.str();
	int iSendResult = send(clientSocket, responseStr.c_str(), (int)responseStr.length(), 0);

	if (iSendResult == SOCKET_ERROR) {
		cout << "发送响应失败，错误码: " << WSAGetLastError() << endl;
		return false;
	}

	return true;
}

// 接收完整的HTTP请求
string ReceiveHTTPRequest(SOCKET clientSocket) {
	string fullData = "";
	char recvbuf[DEFAULT_BUFLEN];
	int recvbuflen = DEFAULT_BUFLEN;
	int iResult;

	// 设置为阻塞模式接收
	u_long mode = 0;
	ioctlsocket(clientSocket, FIONBIO, &mode);

	// 设置接收超时
	int timeout = 30000; // 30秒
	setsockopt(clientSocket, SOL_SOCKET, SO_RCVTIMEO, (const char*)&timeout, sizeof(timeout));

	do {
		iResult = recv(clientSocket, recvbuf, recvbuflen, 0);
		if (iResult > 0) {
			fullData.append(recvbuf, iResult);

			// 检查是否接收完整HTTP请求
			// HTTP请求以\r\n\r\n分隔头和体
			size_t headerEnd = fullData.find("\r\n\r\n");
			if (headerEnd != string::npos) {
				// 检查Content-Length
				size_t lengthPos = fullData.find("Content-Length:");
				if (lengthPos != string::npos) {
					size_t lengthStart = lengthPos + 15;
					size_t lengthEnd = fullData.find("\r\n", lengthStart);
					string lengthStr = fullData.substr(lengthStart, lengthEnd - lengthStart);
					int contentLength = atoi(lengthStr.c_str());

					size_t bodyStart = headerEnd + 4;
					if (fullData.length() >= bodyStart + contentLength) {
						break; // 接收完整
					}
				}
				else {
					// 没有Content-Length，可能是GET请求
					break;
				}
			}
		}
		else if (iResult == 0) {
			cout << "连接关闭" << endl;
			break;
		}
		else {
			int error = WSAGetLastError();
			if (error == WSAETIMEDOUT) {
				cout << "接收超时" << endl;
			}
			else {
				cout << "接收失败，错误码: " << error << endl;
			}
			break;
		}
	} while (iResult > 0);

	return fullData;
}

// 关闭Socket服务器
void StopSocketServer() {
	if (serverSocket != INVALID_SOCKET) {
		closesocket(serverSocket);
		serverSocket = INVALID_SOCKET;
	}
	WSACleanup();
	serverRunning = false;
	cout << "\nSocket服务器已关闭" << endl;
}


bool WriteJsonFile(const string& filepath, const string& json)
{
	ofstream ofs(filepath);
	if (!ofs.is_open())
	{
		cout << "Error: 无法打开文件 " << filepath << " 进行写入！" << endl;
		return false;
	}
	ofs << json;
	ofs.close();
	return true;
}

string ReadJsonFile(const string& filepath)
{
	ifstream ifs(filepath);
	if (!ifs.is_open())
	{
		cout << "Error: 无法打开文件 " << filepath << " 进行读取！" << endl;
		return "";
	}
	stringstream ss;
	ss << ifs.rdbuf();
	ifs.close();
	return ss.str();
}
/*
示例格式：
{"Atom":[[5,"C"],[6,"C"],[7,"C"]],"Adj":[[5,6,"single"],[5,7,"single"]]}

部分解析代码示例：
		string bondsym;
		string bondname = fields[2];
		if (bondname == "single") bondsym = SINGLE_BOND;
		else if (bondname == "double") bondsym = DOUBLE_BOND;
		else if (bondname == "triple") bondsym = TRIPLE_BOND;
		else bondsym = SINGLE_BOND;
*/
bool AnalysisJsonFile(const string& json, vector<SimpleMNode>& sntb, vector<AdjLine>& adj_list)
{
	size_t atom_pos = json.find("\"Atom\"");
	if (atom_pos == string::npos)
	{
		cout << "Error: JSON文件中未找到Atom部分！" << endl;
		return false;
	}
	size_t adj_pos = json.find("\"Adj\"");
	if (adj_pos == string::npos)
	{
		cout << "Error: JSON文件中未找到Adj部分！" << endl;
		return false;
	}
	size_t atom_start = json.find("[", atom_pos);
	size_t atom_end = json.find("]]", atom_start);
	string atom_array_str = json.substr(atom_start + 1, atom_end - atom_start - 1);
	size_t adj_start = json.find("[", adj_pos);
	size_t adj_end = json.find("]]", adj_start);
	string adj_array_str = json.substr(adj_start + 1, adj_end - adj_start - 1);
	// 解析Atom数组
	vector<string> atom_entries;

	vector<int> seq_tb;
	vector<SimpleMNode> mid_sntb;
	vector<AdjLine> mid_adj_list;
	int max_seq = -1;

	Split(atom_entries, atom_array_str, "],");
	for (const string& entry : atom_entries)
	{
		string clean_entry = entry;
		clean_entry.erase(remove(clean_entry.begin(), clean_entry.end(), '['), clean_entry.end());
		clean_entry.erase(remove(clean_entry.begin(), clean_entry.end(), ']'), clean_entry.end());
		vector<string> fields;
		Split(fields, clean_entry, ",");
		if (fields.size() >= 2)
		{
			int seq = atoi(fields[0].c_str());
			string sym = fields[1];
			//去除sym中的引号
			sym.erase(remove(sym.begin(), sym.end(), '\"'), sym.end());

			if (sym == "H") {
				continue;
			}
			if (seq > max_seq) max_seq = seq;
			seq_tb.push_back(seq);
			mid_sntb.push_back(SimpleMNode(seq, sym));
		}
	}
	// 解析Adj数组
	vector<string> adj_entries;
	Split(adj_entries, adj_array_str, "],");
	for (const string& entry : adj_entries)
	{
		string clean_entry = entry;
		clean_entry.erase(remove(clean_entry.begin(), clean_entry.end(), '['), clean_entry.end());
		clean_entry.erase(remove(clean_entry.begin(), clean_entry.end(), ']'), clean_entry.end());
		vector<string> fields;
		Split(fields, clean_entry, ",");
		if (fields.size() >= 3)
		{
			int seq1 = atoi(fields[0].c_str());
			int seq2 = atoi(fields[1].c_str());

			if (seq1 > max_seq || seq2 > max_seq) continue;
			string bondname = fields[2];
			string bondsym;
			if (bondname == "\"single\"") bondsym = SINGLE_BOND;
			else if (bondname == "\"double\"") bondsym = DOUBLE_BOND;
			else if (bondname == "\"triple\"") bondsym = TRIPLE_BOND;
			else bondsym = "";

			if (seq1 > seq2)
			{
				swap(seq1, seq2);
			}
			mid_adj_list.push_back(AdjLine(seq1, seq2, bondsym));
		}
	}

	sort(mid_adj_list.begin(), mid_adj_list.end(),
		[](const AdjLine& a, const AdjLine& b) {
			if (a.GetSeqI() != b.GetSeqI())
				return a.GetSeqI() < b.GetSeqI();
			else
				return a.GetSeqJ() < b.GetSeqJ();
		});
	//通过并查集思想去除其他无关原子

	vector<int> rootfindset(max_seq + 1, -1);
	vector<int> rootfindcount(max_seq + 1, 0);
	int max_root = -1;
	int max_count = -1;

	for (auto& adj : mid_adj_list)
	{
		int seq1 = adj.GetSeqI();
		int seq2 = adj.GetSeqJ();
		if (seq1 > seq2) swap(seq1, seq2);

		if (rootfindset[seq1] == -1 && rootfindset[seq2]==-1)
		{
			rootfindset[seq1] = seq1;
			rootfindset[seq2] = seq1;
			rootfindcount[seq1] = 1;
			//cout << "Connecting: " << seq1 << "< - " << seq2 << endl;
		}
		else if (rootfindset[seq1] != -1 && rootfindset[seq2] != -1)
		{
			//两个都有根节点，合并
			if (rootfindset[seq1] != rootfindset[seq2])
			{
				int old_root = rootfindset[seq2];
				int new_root = rootfindset[seq1];
				for (int i = 0; i <= max_seq; i++)
				{
					if (rootfindset[i] == old_root)
					{
						rootfindset[i] = new_root;
						rootfindcount[new_root]++;
					}
				}
			}
			continue;
		}
		else if (rootfindset[seq1] != -1 && rootfindset[seq2] == -1)
		{
			rootfindset[seq2] = rootfindset[seq1];
			rootfindcount[rootfindset[seq1]]++;
			//cout << "Connecting: " << seq1 << "< - " << seq2 << endl;
			continue;
		}
		else if (rootfindset[seq1] == -1 && rootfindset[seq2] != -1)
		{
			rootfindset[seq1] = rootfindset[seq2];
			rootfindcount[rootfindset[seq2]]++;
			//cout << "Connecting: " << seq2 << "< - " << seq1 << endl;
			continue;
		}
		//rootfindset[seq2] = rootfindset[seq1];
		//rootfindcount[rootfindset[seq1]]++;


	}
	////找出最大的连通分量的根节点
	//PrintCmdSepTitle("连通分量统计结果");
	//PrintCommonVector(rootfindset);
	//PrintCmdSepTitle("根查找结果");
	//PrintCommonVector(rootfindcount);
	//cout << endl;
	for (int i = 0; i <= max_seq; i++)
	{
		if (rootfindcount[i] > max_count)
		{
			max_count = rootfindcount[i];
			max_root = i;
		}
	}

	//for (int i = 0; i <= max_seq; i++)
	//{
	//	if (rootfindset[i] != -1 && rootfindset[i] == max_root)
	//	{
	//		new_seq_tb[i] = count;
	//		count++;
	//	}
	//}
	//cout << "MaxSeq: " << max_seq << endl;
	//cout << seq_tb.size() << endl;
	int count = 0;
	vector<int> new_seq_tb(max_seq + 1, -1);
	for (int i = 0; i < seq_tb.size(); i++)
	{
		int r = seq_tb[i];

		if (rootfindset[r] != -1 && rootfindset[r] == max_root)
		{
			new_seq_tb[r] = count;
			count++;
		}
	}

	if (mid_adj_list.size() == 0)
	{
		new_seq_tb[max_seq] = 0;
	}

	cout << "NewAtomCount: " << count << endl;
	for (auto& node : mid_sntb)
	{
		int old_seq = node.seq;
		node.seq = new_seq_tb[old_seq];

		if (node.seq == -1) continue;

		sntb.push_back(node);
	}
	if (mid_adj_list.size() == 0)
	{
		adj_list = vector<AdjLine>();
	}
	//vector<AdjLine> adj_list_temp;
	for (auto& adj : mid_adj_list)
	{
		int old_seq1 = adj.GetSeqI();
		int old_seq2 = adj.GetSeqJ();

		int new_seq1 = new_seq_tb[old_seq1];
		int new_seq2 = new_seq_tb[old_seq2];

		if (new_seq1 == -1 || new_seq2 == -1) continue;
		adj_list.push_back(AdjLine(new_seq1, new_seq2, adj.GetBondSym()));
	}
	return true;
}

//--- JSON格式解析 ---

vector<SXYZ_3D> GetSXYZTb(const XYZ_TB& xyz_tb, const  EnergySolidParam& esp)
{
	int ac = xyz_tb.size();
	int nhc = esp.short_adj_list.size();
	vector<SXYZ_3D> sxyz_tb(ac, SXYZ_3D());

	for (int i = 0; i < ac; i++)
	{
		if (i < nhc)sxyz_tb[i] = SXYZ_3D(esp.mnode_tb[i].GetSym(), xyz_tb[i]);
		else sxyz_tb[i] = SXYZ_3D("H", xyz_tb[i]);
	}
	return sxyz_tb;


}

vector<AdjABS_3D> GetAdjABSTb(const EnergySolidParam& esp)
{
	int nhc = esp.short_adj_list.size();
	vector<AdjABS_3D> adjline_tb;
	for (int i = 0; i < nhc; i++)
	{
		vector<PointTo> nb_node = esp.short_adj_list[i].GetBonds();
		int nsize = nb_node.size();
		for (int j = 0; j < nsize; j++)
		{
			int des = nb_node[j].GetDesSeq();
			if (des > i)
			{
				adjline_tb.push_back(AdjABS_3D(i + 1, des + 1,nb_node[j].GetBondSymbol()));
			}
		}
	}
	return adjline_tb;
}

void ResultSent(string status, SOCKET clientSocket, string sentjson, string now_smiles, double energy, double opt_time)
{
	string smiles_json_part = "\"SMILES\":\"" + now_smiles + "\"";
	string energy_json_part = "\"Energy\":" + to_string(energy);
	string opttime_json_part = "\"OptTime\":" + to_string(opt_time);
	string responseJson_rec = "{\n\"Status\":\"" + status + "\",\n"
		+ smiles_json_part + ",\n"
		+ energy_json_part + ",\n"
		+ opttime_json_part + ",\n"
		+ "\"Data\":" + sentjson + "\n}";

	cout << "正在发送3D结构数据到网页..." << endl;
	if (SendHTTPResponse(clientSocket, 200, "application/json", responseJson_rec)) {
		cout << "3D结构数据发送成功！" << endl;
		cout << "数据大小: " << responseJson_rec.length() << " 字节" << endl;
	}
	else {
		cout << "3D结构数据发送失败！" << endl;
	}
}