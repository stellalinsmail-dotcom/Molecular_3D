#ifndef TYPES_H
#define TYPES_H

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

using namespace std;

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

#define SEP_WIDTH 50
#define SEP_SYMBOL '-'

#define ATOM_TB_TITLE "Symbol,Number,Mass,Valence"
#define MNODE_TB_TITLE "Seq,Element,IsAroma,CCount,CHCount,NHBC,CharSym,CharVal,MType"
#define ADJ_TB_TITLE "Seq1,Seq2,BondSym"

// --- 通用辅助函数 ---

template<typename T>
void CopyDynArray(T*& new_array,const T*& old_array, int asize)
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
	for (int i = 0; i < s.length(); i++)
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

void PrintCmdSepTitle(const string title, int sepwidth = SEP_WIDTH, const char fillsym = SEP_SYMBOL);

template<class T>
void PrintCommonVector(const vector<T>& v, const string& sep = "\t")
{
	int count = 0;
	for (auto it = v.begin(); it != v.end(); it++)
	{
		if (count != 0)cout << sep;
		cout << *it;
		count++;
	}
}

template<class T>
void PrintSpecialVector(const vector<T>& v, const string& sep = "\t",int max_row_count=-1)
{
	int rowcount = 0;
	for (auto& a : v)
	{
		a.Print(sep);
		if (max_row_count>0)
		{
			rowcount++;
			if (rowcount >= max_row_count) break;
		}
	}
}

string GetCurrentTimeString();

// ---表格读写函数 ---
bool CreateFolder(const string folderpath);

template<typename T>
int ReadTable(string filename, string stdtitle, vector<T>& table,bool printyes=false)
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
int ReadTableByTitle(string filename, string stdtitle, vector<T>& table,bool printyes=false)
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
			if (colsize==1)cout << str << endl;
		}
		T info(stemp, seq);
		table.push_back(info);
			//info.Print();
			//cout << endl;
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
		outFile << line;
	}
	outFile.close();
}

#endif // TYPES_H
