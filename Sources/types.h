#pragma once

#ifndef TYPES_H
#define TYPES_H

// --- Socket通信相关 (必须在windows.h之前) ---
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#endif

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <cmath>
#include <windows.h>
#include <iomanip>
#include <set>
#include <map>
#include <functional>
#include <random>

using namespace std;


#define DEFAULT_PORT "8765"
#define DEFAULT_BUFLEN 65536

#define MAX_CAL_TIME 500 //单位s,
#define MAX_REC_TIME 604800 //单位s (7天)

// --- 常量定义 ---
#define YES true
#define NO false
#define TITLE_ON 1
#define TITLE_OFF 0

#define VAL_SEP '/'
#define TITLE_SEP ','

#define SINGLE_BOND "-"
#define DOUBLE_BOND "="
#define TRIPLE_BOND "#"
#define AROMA_BOND "A"
#define HYDROGEN_BOND "H"

#define CC_UNIT 1
#define NH_UNIT 1
#define AN_UNIT 3
#define CS_UNIT 1
#define CV_UNIT 1
#define CH_UNIT 1

#define MAX_CIRCLE_COUNT 1000

#define COMMON_SORT 1
#define INIT_SORT 2

#define H1_SEP_WIDTH 60
#define H2_SEP_WIDTH 50
#define H3_SEP_WIDTH 40
#define H4_SEP_WIDTH 30


#define SEP_SYMBOL '-'

#define ATOM_TB_TITLE "Symbol,Number,Mass,Valence"
#define MNODE_TB_TITLE "Seq,Element,IsAroma,CCount,CHCount,NHBC,CharSym,CharVal,MType"
#define ADJ_TB_TITLE "Seq1,Seq2,BondSym"
#define OPT_REC_TB_TITLE "Smiles,MinEnergy,OptTime"
// --- MMFF94 ---

#define BS_TB_TITLE "BT,MTYPE_I,MTYPE_J,KB,R0,SOURCE"
#define AB_TB_TITLE "ANGLE_TYPE,MTYPE_I,MTYPE_J,MTYPE_K,KA_IJK,THETA_0,COMMENT_ORIGIN"
#define SB_TB_TITLE "SBT,MTYPE_I,MTYPE_J,MTYPE_K,KBA_IJK,KBA_KJI,SOURCE"
#define OPB_TB_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,KOOP,Source"
#define TI_TB_TITLE	"TT,MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,V1,V2,V3,Source"
#define VDW_TB_TITLE "MTYPE,ALPHA_I,N_I,A_I,G_I,DA,SYMBOL,ORIGIN"
#define MSYM_TB_TITLE "MSYMBOL,MTYPE,DEFINITION"


#define BS_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,KB,R0"
#define AB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,KA_IJK,THETA_0"
#define SB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,KBA_IJK,KBA_KJI"
#define OPB_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE,KOOP"
#define TI_TB_SHORT_TITLE "MTYPE_I,MTYPE_J,MTYPE_K,MTYPE_L,V1,V2,V3"
#define VDW_TB_SHORT_TITLE "MTYPE,ALPHA_I,N_I,A_I,G_I"
#define	MSYM_TB_SHORT_TITLE "MSYMBOL,MTYPE"


#define SP_SPTB_TITLE "Seq,R,Theta,Varphi"
#define SP_RETB_TITLE "Seq,X,Y,Z"
#define XYZ_TB_TITLE "X,Y,Z"

#define SXYZ_TB_TITLE "Symbol,X,Y,Z"
#define ADJ_AB_TB_TITLE "A,B"


#define MMFF_CSV_FOLDER "File/SolidParam/MMFF94_CSV/"
#define SP_FOLDER "File/SolidParam/C_SP_CSV/"
#define MMFF_OUTPUT_FOLDER "File/Output/CSV/"
#define JSON_OUTPUT_FOLDER "File/OutputOpt/Json"
#define OPT_OUTPUT_FOLDER "File/OutputOpt/"
#define OPT_TIME_FOLDER "File/OutputOpt/TimeRecord/"


#define ETB_FILEPATH "File/SolidParam/EL_PN/ElementsTable.csv";
#define PTB_FILEPATH "File/SolidParam/EL_PN/PrimeNumber1000.csv";


#define BS_FILENAME  "6_MMFFBOND.csv"
#define AB_FILENAME "8_MMFFANG.csv"
#define SB_FILENAME "9_MMFFSTBN.csv"
#define OPB_FILENAME "11_MMFFOOP.csv"
#define TI_FILENAME "12_MMFFTOR.csv"
#define VDW_FILENAME "13_MMFFVDW.csv"

#define MSYM_FILENAME "1_MMFFSYMB.csv"
#define OPT_REC_FILENAME "AlphaOptRecord.csv"
#define OPT_COPY_FILENAME "CopyRecord.csv"

#define PI 3.14159265358979323846
#define PI_HALF 1.57079632679489661923
#define PI_DOUBLE 6.2831853071795864769252
#define PI_270 4.71238898038468985769
#define cosd(x) cos(x * PI / 180.0)
#define sind(x) sin(x * PI / 180.0)
#define acosd(x) (180.0 / PI) * acos(x)


#define Min3(a,b,c) min(min(a,b),c)
#define Max3(a,b,c) max(max(a,b),c)
#define Max4(a,b,c,d) max(max(a,b),max(c,d))
#define Min4(a,b,c,d) min(min(a,b),min(c,d))

//#define CONST_PARAM_MAP const map<vector<int>, double>


#define MSYM_MAP map<string, int>
#define MTYPE_TB vector<int>
#define MTYPE_INDEX vector<int>
#define MTYPE_SET set<int>


#define MNODE_TB vector<MNode>

#define ANGLE_TB vector<double>
#define HASH_TB vector<int>
#define ADJ_LIST vector<NodeBonds>
#define XYZ_TB vector<Vec3>

#define FAST_TABLE_INDEX vector<int>
#define R_TB vector<vector<double>>
#define VAR_TB vector<vector<vector<double>>>
#define CHI_TB vector<vector<vector<vector<double>>>>
#define PHI_TB vector<vector<vector<vector<double>>>>
#define RE_TB vector<vector<VDWProVal>>
#define BFS2_TB vector<vector<int>>



#define V4_DTB vector<vector<vector<vector<double>>>>
#define V3_DTB vector<vector<vector<double>>>
#define V2_DTB vector<vector<double>>
#define V1_DTB vector<double>


#define F_BS FastMatrix2<BSLine, BSVal>
#define F_AB FastMatrix3<ABLine, ABVal>
#define F_SB FastMatrix3<SBLine, SBVal>
#define F_OPB FastMatrix4<OPBLine, OPBVal>
#define F_TI FastMatrix4<TILine, TIVal>
#define F_VDW FastMatrix1<VDWLine, VDWVal>



//#define ALL_SP_SPTB vector<vector<SP_SPLine>>    

#define MAX_SP_SIZE 6
#define MAX_ADJ_NODE_SIZE 6
#define OPT_RATIO 0.618

// --- 优化记录结构体 ---
struct OptimizationRecord {
	double time_ms;           // 时间（毫秒）
	long long iteration_count;      // 优化次数
	double min_energy;        // 最低能量

	OptimizationRecord(double t, long long iter, double e)
		: time_ms(t), iteration_count(iter), min_energy(e) {
	}
};

// 输出运算符重载（用于CSV导出）
inline ostream& operator<<(ostream& os, const OptimizationRecord& rec) {
	os << rec.time_ms << "," << rec.iteration_count << "," << rec.min_energy;
	return os;
}

// --- 文件路径处理 ---
inline string GetMMFFPath(string filename)
{
	return MMFF_CSV_FOLDER + filename;
}

inline string GetMNodeTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/MNodeTable.csv";
}
inline string GetAdjTbPath(string can_smiles)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/BondTable.csv";
}
inline string GetXYZTbPath(string can_smiles, int n)
{
	return MMFF_OUTPUT_FOLDER + can_smiles + "/XYZ_" + to_string(n) + ".csv";
}


// --- 通用辅助函数 ---

template<typename T>
void CopyDynArray(T*& new_array, const T*& old_array, int asize)
{
	new_array = new T[asize];
	for (int i = 0; i < asize; i++) new_array[i] = old_array[i];
}

char CharUpper(char a);
string StringLower(string a);
int LCSArray(const string& a, const string& b);

template<typename T>
bool SearchSimilarity(const vector<T>& table, const string& des, int& maxid, string(*tbfun)(T)) {
	int ratio = 0, ratioid = -1;
	for (int j = 0; j < table.size(); j++) {
		int a = LCSArray(tbfun(table[j]), des);
		if (ratio < a) {
			ratio = a;
			ratioid = j;
		}
	}
	maxid = ratioid;
	return ratio > 0;
}

bool SearchSimilarity(const vector<string>& table, const string& des, int& maxid);

template<typename T>
void Split(vector<T>& v, string s, const char ch)
{
	int start = 0, end = 0;
	v.clear();
	int i = 0;
	if (s.length() == 0) return;
	if (s[0] == ch)
	{
		start = 1;
	}
	for (i = start; i < s.length(); i++)
	{
		end = i;
		if (s[i] != ch && i != s.length() - 1) { continue; }
		else
		{
			int len = end - start;
			if (i == s.length() - 1) len++;
			if (len)
			{
				string s_cut = s.substr(start, len);
				v.push_back(s_cut);
			}
			start = end + 1;
		}
	}
}
// --- 字符串处理函数 ---
template<typename T>
void Split(vector<T>& v, string s, string ch)
{
	int start = 0, end = 0;
	v.clear();
	int chlen = ch.length();
	for (int i = 0; i < s.length(); i++)
	{
		end = i;
		if (s.substr(i, chlen) != ch && i != s.length() - chlen) { continue; }
		else
		{
			int len = end - start;
			if (i == s.length() - chlen) len += chlen;
			if (len)
			{
				string s_cut = s.substr(start, len);
				v.push_back(s_cut);
			}
			start = end + chlen;
			i += chlen - 1;
		}
	}
}
void SplitToInt(vector<int>& v, string s, const char ch);
void GetTitleSortSeq(vector<int>& seq, string stdtitle, const vector<string>& tabletitle);
bool IsLetter(const char& ch);
bool IsLowerLetter(const char ch);
bool IsNumber(const char& ch);
bool IsCircleNumber(const char& ch);
int MatchPunc(const char leftch, string str, int startpos);
int FindSingleAtom(string str, int startpos);
int FindLastCircleNumber(string str, int startpos);
string DeleteExtraBracket(string s);

// --- 输出辅助函数 ---
void PrintTableTitle(const string& text, const string outsep = "\t", const char insep = TITLE_SEP);


void PrintCmdSepTitle(const string title, int sepwidth = H1_SEP_WIDTH, const char fillsym = SEP_SYMBOL);

template<class T>
void PrintCommonVector(const vector<T>& v, const string& sep = "\t",int max_row_count = -1)
{
	int count = 0;
	for (auto it = v.begin(); it != v.end(); it++)
	{
		if (count != 0)cout << sep;
		cout << *it;
		if (max_row_count > 0)
		{
			if (count >= max_row_count) break;
		}
		count++;
	}
}

template<class T>
void PrintSpecialVector(const vector<T>& v, const string& sep = "\t", int max_row_count = -1)
{
	int rowcount = 0;
	for (auto& a : v)
	{
		a.Print(sep);
		cout << endl;
		if (max_row_count > 0)
		{
			rowcount++;
			if (rowcount >= max_row_count) break;
		}
	}
}

template<typename T>
void PrintSet(const set<T>& s, const string& sep = "\t", int max_row_count = -1)
{
	int rowcount = 0;
	for (auto& element : s)
	{
		cout << element << sep << endl;
		if (max_row_count > 0)
		{
			rowcount++;
			if (rowcount >= max_row_count) break;
		}
	}
}

void PrintEnergy(double sum_E, double sum_eb, double sum_ea, double sum_eba, double sum_eoop, double sum_et, double sum_evdw);

template<typename TIndex, typename TVal>
void PrintCommonMap(const map<TIndex, TVal>& m, const string& sep = "\t", int max_row_count = -1)
{
	//int n = TVal().GetPCount();
	//string title = "MType" + repeat_str(",Keys", n);
	//PrintTableTitle(title, sep);
	int rowcount = 0;
	for (auto& pair : m)
	{
		cout << "(" << pair.first << ")" << sep << "(" << pair.second << ")" << endl;
		if (max_row_count > 0)
		{
			rowcount++;
			if (rowcount >= max_row_count) break;
		}
	}
}


string GetCurrentTimeString();

// ---文件读写函数 ---
bool CreateFolder(const string folderpath);

void DeleteFolder(const string& folderpath);

template<typename T>
int ReadTable(string filename, string stdtitle, vector<T>& table, bool printyes = false)
{
	ifstream ifs;
	ifs.open(filename, ios::in);
	if (!ifs.is_open())
	{
		cout << "表格打开失败！\n";
		return 0;
	}
	vector<string> item;		//用于存放文件中的一行数据
	string temp;				//临时存储文件中的一行数据
	while (getline(ifs, temp))  //利用 getline（）读取文件中的每一行，并放入到 item 中
	{
		item.push_back(temp);
	}
	long int k = -1;
	int colsize = 0;
	vector<string>tabletitle;
	//vector<int> seq;

	for (auto it = item.begin(); it != item.end(); it++)
	{
		if (k < 0)
		{
			k++;
			string oldstr;
			string str;
			istringstream istr(*it);
			do {
				getline(istr, str, ',');
				if (str == oldstr) break;
				tabletitle.push_back(str);
				colsize++;
				oldstr = str;
			} while (str != "\n" && str != "");
			//cout << "all colsize:" << colsize << endl;
			//cout << "表格内容：\n";
		}
		string str;
		istringstream istr(*it);
		vector<string> stemp;
		//将字符串流 istr 中的字符读入到 str 字符串中，读取方式是以逗号为分隔符 
		for (int j = 0; j < colsize; j++)
		{
			getline(istr, str, ',');
			stemp.push_back(str);
			//if (colsize==1)cout << str << endl;
		}
		T info(stemp);
		table.push_back(info);
		//cout << stemp[0] << endl;
		k++;
	}
	ifs.close();
	item.clear();
	if (printyes) cout << "文件" << filename << "加载完毕！\n共有 " << k << " 行 " << colsize << " 列记录!" << endl;
	return k;
}
template<typename T>
int ReadTableByTitle(string filename, string stdtitle, vector<T>& table, bool printyes = false)
{
	ifstream ifs;
	ifs.open(filename, ios::in);
	if (!ifs.is_open())
	{
		cout << "表格 "+filename + " 打开失败！\n";
		return -1;
	}
	vector<string> item;		//用于存放文件中的一行数据
	string temp;				//临时存储文件中的一行数据
	while (getline(ifs, temp))  //利用 getline（）读取文件中的每一行，并放入到 item 中
	{
		item.push_back(temp);
	}
	long int k = -1;
	int colsize = 0;
	vector<string>tabletitle;
	vector<int> seq;

	for (auto it = item.begin(); it != item.end(); it++)
	{
		if (k < 0)
		{
			k++;
			string oldstr;
			string str;
			istringstream istr(*it);
			do {
				getline(istr, str, ',');
				if (str == oldstr) break;
				tabletitle.push_back(str);
				colsize++;
				oldstr = str;
			} while (str != "\n" && str != "");
			//cout << "all colsize:" << colsize << endl;
			//cout << "表格内容：\n";
			GetTitleSortSeq(seq, stdtitle, tabletitle);
			continue;
		}
		string str;
		istringstream istr(*it);
		vector<string> stemp;
		//将字符串流 istr 中的字符读入到 str 字符串中，读取方式是以逗号为分隔符 
		for (int j = 0; j < colsize; j++)
		{
			getline(istr, str, ',');
			stemp.push_back(str);
			if (colsize == 1)cout << str << endl;
		}
		T info(stemp, seq);
		table.push_back(info);
		//info.Print();
		//cout << endl;
		k++;
	}
	ifs.close();
	item.clear();
	if (printyes) cout << "文件" << filename << "加载完毕！\n共有 " << k << " 行 " << colsize << " 列记录!" << endl;
	return k;
}

template<typename T>
void WriteTable(string filename, string stdtitle, const  vector<T>& table, bool writetitle = YES)
{
	ofstream outFile(filename, ios::out);
	if (writetitle) outFile << stdtitle << endl;
	for (auto& line : table)
	{
		outFile << line << endl;
	}
	outFile.close();
}

bool FileExists(const string& filename);

string ReadFileContent(const string& filepath);

// --- 表格处理函数 ---
template<typename TLine,typename TVal>
map<string, TVal> ChangeTableToMap(const vector<TLine>& table)
{
	map<string, TVal> table_map;
	for (auto& line : table)
	{
		table_map[line.GetIndex()] = line.GetVal();
	}
	return table_map;
}

template<typename TLine,typename TVal>
vector<TLine> ChangeMapToTable(const map<string, TVal>& table_map)
{
	vector<TLine> table;
	for (auto& pair : table_map)
	{
		TLine line(pair.first, pair.second);
		table.push_back(line);
	}
	return table;
}


//--- 字典打印输出 ---
// 
// 
// template<typename T>
// using PARAM_MAP = std::map<std::vector<int>, T>;
//string repeat_str(const char* str, int n)
//{
//	string result = "";
//	for (int i = 0; i < n; i++)
//	{
//		result += str;
//	}
//	return result;
//}
//template <typename TVal>
//void PrintParamMap(const map<vector<int>, TVal>& m, const string& sep = "\t", int max_row_count = -1)
//{
//	int n = TVal().GetPCount();
//	string title = "MType" + repeat_str(",Keys", n);
//	PrintTableTitle(title, sep);
//	int rowcount = 0;
//	for (auto& pair : m)
//	{
//		cout << "(";
//		PrintCommonVector(pair.first, ",");
//		cout << ")" << sep << "(" << pair.second << ")" << endl;
//		if (max_row_count > 0)
//		{
//			rowcount++;
//			if (rowcount >= max_row_count) break;
//		}
//	}
//}
// 


#endif // TYPES_H
