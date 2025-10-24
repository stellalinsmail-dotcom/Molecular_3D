#include "types.h"

char CharUpper(char a)
{
	if (a >= 'a' && a <= 'z') return a + 'A' - 'a';
	else return a;
}

string StringLower(string a)
{
	for (int i = 0; i < a.length(); i++)
	{
		if (a[i] >= 'A' && a[i] <= 'Z') a[i] = a[i] - 'A' + 'a';
	}
	return a;
}

int LCSArray(const string& a, const string& b)
{
	int alen = a.length(), blen = b.length();
	vector<vector<int>> m(alen + 1, vector<int>(blen + 1, 0));
	for (int i = 1; i <= alen; i++)
	{
		for (int j = 1; j <= blen; j++)
		{
			if (CharUpper(a[i - 1]) == CharUpper(b[j - 1]))  m[i][j] = m[i - 1][j - 1] + 1;
			else m[i][j] = max(m[i][j - 1], m[i - 1][j]);
		}
	}
	return m[alen][blen];
}

bool SearchSimilarity(const vector<string>& table, const string& des, int& maxid)
{
	return SearchSimilarity<string>(table, des, maxid, [](string a)->string { return a; });
}

void PrintTableTitle(const string& text,const string outsep, const char insep)
{
	for (int i = 0; i < text.length(); i++)
	{
		if (text[i] == insep) cout << outsep;
		else cout << text[i];
	}
	cout << endl;
}

void SplitToInt(vector<int>& v, string s, const char ch)
{
	vector<string> vs;
	Split(vs, s, ch);

	for (int i = 0; i < vs.size(); i++)
		v.push_back(atoi(vs[i].c_str()));
	vs.clear();
}

void GetTitleSortSeq(vector<int>& seq, string stdtitle, const vector<string>& tabletitle)
{
	vector<string> std;
	vector<string> tb(tabletitle);
	Split(std, stdtitle, ',');
	seq.clear();

	//PrintCommonVector(std);
	//cout << endl;
	//PrintCommonVector(tb);
	//cout << endl;

	for (int i = 0; i < std.size(); i++)
	{
		int maxid;
		if (SearchSimilarity(tb, std[i], maxid))
		{
			seq.push_back(maxid);
			tb[maxid] = "";
		}
		else cout << "查找序列失败！\n";
	}
	std.clear();
	tb.clear();
}

bool IsLetter(const char& ch)
{
	if (ch >= 'a' && ch <= 'z' || ch >= 'A' && ch <= 'Z') return true;
	return false;
}

bool IsLowerLetter(const char ch)
{
	if (ch >= 'a' && ch <= 'z') return true;
	return false;
}

bool IsNumber(const char& ch)
{
	if (ch >= '0' && ch <= '9') return true;
	return false;
}

bool IsCircleNumber(const char& ch)
{
	if (IsNumber(ch) || ch == '%') return true;
	return false;
}

int MatchPunc(const char leftch, string str, int startpos)
{
	int count = 1;
	char rightch;
	switch (leftch)
	{
	case '[': rightch = ']'; break;
	case '(': rightch = ')'; break;
	case '{': rightch = '}'; break;
	case '<': rightch = '>'; break;
	default: return -1;
	}
	for (int i = min(startpos + 1, str.length() - 1); i < str.length(); i++)
	{
		if (str[i] == rightch) count--;
		if (str[i] == leftch) count++;
		if (!count) return i;
	}
}
int FindSingleAtom(string str, int startpos)
{
	//for (int i = startpos + 1; i <min(str); i++)
	//{
	//	if (!IsLetter(str[i])) return i - 1;
	//	cout << i << endl;
	//}
	return startpos;
}
int FindLastCircleNumber(string str, int startpos)
{
	for (int i = min(startpos + 1, str.length() - 1); i < str.length(); i++)
	{
		if (!IsCircleNumber(str[i])) return i - 1;
	}
	return str.length() - 1;
}

string DeleteExtraBracket(string s)
{
	int leftend = 0, rightend = s.length(), slen = s.length();
	while (leftend < rightend)
	{
		rightend--;
		if (s[rightend] == ')' && (rightend == slen - 1 || s[rightend + 1] == ')'))
		{
			while (leftend < rightend && s[leftend] != '(')
			{
				leftend++;
			}
			s.erase(rightend, 1);
			rightend--;
			s.erase(leftend, 1);
			slen = s.length();
		}

	}
	return s;
}

string GetCurrentTimeString() {
	SYSTEMTIME st;
	GetLocalTime(&st); // 获取本地时间

	ostringstream oss;

	// 格式: YYYYMMDD_HHMMSS
	oss << setfill('0')
		<< setw(4) << st.wYear
		<< setw(2) << st.wMonth
		<< setw(2) << st.wDay
		<< "_"
		<< setw(2) << st.wHour
		<< setw(2) << st.wMinute
		<< setw(2) << st.wSecond;

	return oss.str();
}

bool CreateFolder(const string folderpath)
{
	size_t convertedChars = 0;

	wchar_t* folderPath = new wchar_t[folderpath.length() + 1];
	mbstowcs_s(&convertedChars, folderPath, folderpath.length() + 1, folderpath.c_str(), _TRUNCATE);
	bool issuccessful = CreateDirectory(folderPath, NULL);
	delete[] folderPath;
	return issuccessful;
}

void PrintCmdSepTitle(const string title, int sepwidth, const char fillsym)
{
	int scount = max((sepwidth - (int)title.length() - 6) / 2, 0);
	string sidesep(scount, fillsym);
	string sum = "\n" + sidesep + " * " + title + " * " + sidesep + "\n\n";
	cout << sum;
}


void PrintEnergy(double sum_E, double sum_eb, double sum_ea, double sum_eba, double sum_eoop, double sum_et, double sum_evdw)
{
	PrintCmdSepTitle("Sum E");
	cout << fixed << setprecision(8);
	cout << "sum_EB: " << sum_eb << endl;
	cout << "sum_EA: " << sum_ea << endl;
	cout << "sum_EBA: " << sum_eba << endl;
	cout << "sum_EOOP: " << sum_eoop << endl;
	cout << "sum_ET: " << sum_et << endl;
	cout << "sum_EVDW: " << sum_evdw << endl;

	cout << "*sum_E: " << sum_E << endl;
}