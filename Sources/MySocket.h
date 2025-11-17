#pragma once

#ifndef MYSOCKET_H
#define MYSOCKET_H

#include "types.h"
#include "Molecule.h"
#include "Matrix.h"
#include "MMFF94Cal.h"

// --- Socket Server Functions ---
extern SOCKET serverSocket;
extern bool serverRunning;

// 初始化Winsock
bool InitializeWinsock();

// 启动Socket服务器 (使用HTTP协议)
bool StartSocketServer(const char* port);

// 解析HTTP请求头
string ParseHttpRequest(const string& request, string& method, string& path);

// 发送HTTP响应
bool SendHTTPResponse(SOCKET clientSocket, int statusCode, const string& contentType, const string& body);

// 接收完整的HTTP请求
string ReceiveHTTPRequest(SOCKET clientSocket);

// 关闭Socket服务器
void StopSocketServer();

template<typename T>
string VectorToJsonArray(string stdtitle, const vector<T>& v)
{
	//title_tb存储各个字段的标题
	vector<string> title_tb;
	Split(title_tb, stdtitle, ',');
	string json = "{\n";
	int n = v.size();
	//通过循环输入v中的每个元素，格式为"title":value，value的顺序由重载的operator<<决定，与title_tb顺序一致
	for (int i = 0; i < n; i++)
	{
		json += "\t\"" + to_string(i) + "\": {";
		stringstream ss;
		ss << v[i];
		string line = ss.str();
		vector<string> value_tb;
		Split(value_tb, line, ',');

		int m = title_tb.size();
		for (int j = 0; j < m; j++)
		{
			json += "\"" + title_tb[j] + "\": ";
			json += value_tb[j];
			if (j < m - 1) json += ", ";
		}
		json += "}";

		if (i < n - 1) json += ",\n";
		else json += "\n";
	}


	json += "}\n";
	return json;

}

template<typename T>
string VectorToJsonPart(string nametitle, const vector<T>& v)
{
	//title_tb存储各个字段的标题
	//vector<string> title_tb;
	//Split(title_tb, stdtitle, ',');

	string json = "\"" + nametitle + "\": [\n";
	size_t n = v.size();
	//通过循环输入v中的每个元素，格式为"title":value，value的顺序由重载的operator<<决定，与title_tb顺序一致
	for (size_t i = 0; i < n; i++)
	{
		json += "[";
		stringstream ss;
		ss << v[i];
		string line = ss.str();
		vector<string> value_tb;
		Split(value_tb, line, ',');

		size_t m = value_tb.size();
		for (size_t j = 0; j < m; j++)
		{
			json += value_tb[j];
			if (j < m - 1) json += ", ";
		}
		json += "]";

		if (i < n - 1) json += ",\n";
		else json += "\n";
	}

	json += "]";
	return json;
}


bool WriteJsonFile(const string& filepath, const string& json);

string ReadJsonFile(const string& filepath);
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
bool AnalysisJsonFile(const string& json, vector<SimpleMNode>& sntb, vector<AdjLine>& adj_list);

// --- 面向3D网页的导出 ---

class SXYZ_3D {
private:
	string sym;
	Vec3 xyz;
public:
	SXYZ_3D() : sym("None"), xyz(Vec3()) {}
	SXYZ_3D(string s, Vec3 v) :sym(s), xyz(v) {}
	string GetSym() const { return sym; }
	Vec3 GetXYZ() const { return xyz; }
	void Print(string sep = "\t")const
	{
		cout << sym << sep;
		xyz.Print(sep);
	}
	friend ostream& operator<<(ostream& os, const SXYZ_3D& sxyz)
	{
		os << "\"" << sxyz.sym << "\"," << sxyz.xyz;
		return os;
	}
};

class AdjABS_3D
{
private:
	int seq1;
	int seq2;
	string bondsym;
public:
	AdjABS_3D(const vector<string>& info, const vector<int>& titlenum)
	{
		seq1 = atoi(info[titlenum[0]].c_str()) + 1;
		seq2 = atoi(info[titlenum[1]].c_str()) + 1;
		bondsym = info[titlenum[2]];
	}
	AdjABS_3D(int s1 = -0, int s2 = 0,string sym="") :seq1(s1), seq2(s2), bondsym(sym) {}
	int GetIndex1() const { return seq1; }
	int GetIndex2() const { return seq2; }
	void Print(string sep = "\t")const
	{
		cout << seq1 << sep << seq2 << endl;
	}
	friend ostream& operator<<(ostream& os, const AdjABS_3D& adjline)
	{
		os << adjline.seq1 << "," << adjline.seq2<<",\""<<adjline.bondsym<<"\"";
		return os;
	}

};


//--- JSON格式解析 ---

vector<SXYZ_3D> GetSXYZTb(const XYZ_TB& xyz_tb, const  EnergySolidParam& esp);

vector<AdjABS_3D> GetAdjABSTb(const EnergySolidParam& esp);

void ResultSent(string status, SOCKET clientSocket, string sentjson, string now_smiles, double energy, double opt_time);



#endif


