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
void DeleteFolder(const string& folderpath) {
	if (folderpath.empty()) {
		cout << "错误：文件夹路径为空" << endl;
		return;
	}

	// 将 UTF-8 std::string 转为 wide string（WinAPI 使用 UTF-16）
	auto Utf8ToWide = [](const string& s)->wstring {
		if (s.empty()) return wstring();
		int req = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, nullptr, 0);
		if (req == 0) return wstring();
		wstring w(req, L'\0');
		MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, &w[0], req);
		if (!w.empty() && w.back() == L'\0') w.pop_back();
		return w;
		};

	// 递归删除宽字符路径
	function<void(const wstring&)> DeleteFolderWide = [&](const wstring& wfolder) {
		// 检查是否存在
		DWORD attr = GetFileAttributesW(wfolder.c_str());
		if (attr == INVALID_FILE_ATTRIBUTES) return;

		// 如果是符号链接/重解析点，直接删除目录（避免递归穿越）
		if (attr & FILE_ATTRIBUTE_REPARSE_POINT) {
			SetFileAttributesW(wfolder.c_str(), FILE_ATTRIBUTE_NORMAL);
			RemoveDirectoryW(wfolder.c_str());
			return;
		}

		// 构造查找模式：folder\*
		std::wstring search = wfolder;
		if (!search.empty() && (search.back() != L'\\' && search.back() != L'/')) search += L'\\';
		search += L'*';

		WIN32_FIND_DATAW ffd;
		HANDLE hFind = FindFirstFileW(search.c_str(), &ffd);
		if (hFind != INVALID_HANDLE_VALUE) {
			do {
				const std::wstring name = ffd.cFileName;
				if (name == L"." || name == L"..") continue;

				std::wstring full = wfolder;
				if (!full.empty() && (full.back() != L'\\' && full.back() != L'/')) full += L'\\';
				full += name;

				if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
					// 递归删除子目录
					DeleteFolderWide(full);
				}
				else {
					// 先清除只读等属性，再删除文件
					SetFileAttributesW(full.c_str(), FILE_ATTRIBUTE_NORMAL);
					if (!DeleteFileW(full.c_str())) {
						// 可选：打印错误代码
						// cout << "无法删除文件: " << WideCharToMultiByte(CP_UTF8,0,full.c_str(),-1,...) << endl;
					}
				}
			} while (FindNextFileW(hFind, &ffd));
			FindClose(hFind);
		}

		// 删除目录自身（清除属性后删除）
		SetFileAttributesW(wfolder.c_str(), FILE_ATTRIBUTE_NORMAL);
		RemoveDirectoryW(wfolder.c_str());
		};

	std::wstring wpath = Utf8ToWide(folderpath);
	// 去掉末尾斜杠（统一处理）
	if (!wpath.empty() && (wpath.back() == L'\\' || wpath.back() == L'/')) wpath.pop_back();

	// 检查是否存在并且确实是目录
	DWORD finalAttr = GetFileAttributesW(wpath.c_str());
	if (finalAttr == INVALID_FILE_ATTRIBUTES) {
		cout << "删除失败：路径不存在: " << folderpath << endl;
		return;
	}
	if (!(finalAttr & FILE_ATTRIBUTE_DIRECTORY)) {
		cout << "删除失败：不是目录: " << folderpath << endl;
		return;
	}

	DeleteFolderWide(wpath);
	cout << "已尝试删除文件夹及其内容: " << folderpath << endl;
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