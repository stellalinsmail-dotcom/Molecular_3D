/*
一、一些暂时忽略的问题
  1、据资料说string可以容纳几百万字符，所以理论上应该长SMILES不会溢出
  2、限制了环生成时的断键方法，要求必须为单键。
  3、暂时忽略GENE算法中对环的多重键优先的规则；
  %%3、二符号元素的判断
*/

#include<iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <windows.h>

//仅用于输出：
#include <iomanip>


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

#define ATOM_TITLE "Symbol,Number,Mass,Valence"
#define MNODE_TITLE "Seq,Element,IsAroma,CCount,CHCount,NHBC,CharSym,CharVal,MType"
#define ADJ_TITLE "Seq1,Seq2,BondSym"

using namespace std;

template<typename T>
void CopyDynArray(T*& new_array, T* old_array, int asize)
{
	new_array = new T[asize];
	for (int i = 0; i < asize; i++) new_array[i] = old_array[i];
}

template<class T>
void PrintCommonVector(const vector<T>& v, const string& sep = "\t")
{
	for (auto it = v.begin(); it != v.end(); it++)
	{
		cout << *it << sep;
	}
}
template<class T>
void PrintSpecialVector(const vector<T>& v, const string& sep = "\t")
{
	for (auto& a : v)
	{
		a.Print();
	}
}

//int LCS(const string &a,const string &b) {
//	int alen = a.length(), blen = b.length();
//	if (!alen || !blen) return 0;
//	if (a[0] == b[0]) return LCS(a.substr(1), b.substr(1))+1;
//	else return max(LCS(a.substr(0), b.substr(1)), LCS(a.substr(1), b.substr(0)));
//}

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

bool SearchSimilarity(const vector<string>& table, const string& des, int& maxid)
{
	return SearchSimilarity<string>(table, des, maxid, [](string a)->string { return a; });
}

//输入格式为insep，输出格式为outsep
void PrintTableTitle(const string& text, string outsep = "\t", const char insep = TITLE_SEP)
{
	for (int i = 0; i < text.length(); i++)
	{
		if (text[i] == insep) cout << outsep;
		else cout << text[i];
	}
	cout << endl;
}
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

class Atom {
private:
	string sym; //原子符号
	int num; //原子序数
	int mass; //相对原子质量
	vector<int> val; //原子常见化合价
public:
	Atom(string symbol = "None") :sym(symbol), num(0), mass(0), val(0, 0) {};
	Atom(const Atom& e) :sym(e.sym), num(e.num), mass(e.mass), val(e.val) {}
	Atom(string symbol, int number, int relamass, vector<int> valence) :num(number), mass(relamass), sym(symbol) {
		val.assign(valence.begin(), valence.end());
	}
	Atom(const vector<string>& atom_info)
	{
		sym = atom_info[0];
		num = atoi(atom_info[1].c_str());
		mass = atoi(atom_info[2].c_str());
		SplitToInt(val, atom_info[3], VAL_SEP);
		sort(val.begin(), val.end());
	}
	Atom(const vector<string>& atom_info, const vector<int>& titlenum)
	{
		sym = atom_info[titlenum[0]];
		num = atoi(atom_info[titlenum[1]].c_str());
		mass = atoi(atom_info[titlenum[2]].c_str());
		SplitToInt(val, atom_info[titlenum[3]], VAL_SEP);
		sort(val.begin(), val.end());
	}
	~Atom() { val.clear(); }
	int GetNum()const { return num; }
	string GetSym()const { return sym; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const
	{
		if (title_state) {
			PrintTableTitle(ATOM_TITLE);
			cout << endl;
			//cout << "symbol" << sep << "number" << sep << "rmass" << sep << "valence" << endl;
		}
		cout << sym << sep << num << sep << mass << sep;
		PrintCommonVector(val, " ");
		cout << endl;
	}
	friend Atom SearchAtom(const vector<Atom>& atb, string symbol, bool isprecise = NO)
	{
		int index = 0;
		while (index < symbol.length())
		{
			if (!IsLetter(symbol[index]))
				symbol.erase(index, 1);
			else
				index++;
		}
		if (isprecise)
		{
			for (auto& v : atb)
			{
				if (v.sym == symbol) return v;
			}
		}
		int maxid;
		if (SearchSimilarity<Atom>(atb, symbol, maxid, [](Atom a)->string { return a.GetSym(); }))
			return Atom(atb[maxid]);

		return Atom("None");
	}
	int GetHigherVal(int n)
	{
		if (val.empty()) return 0;
		for (auto it = val.begin(); it != val.end(); it++)
		{
			if (*it >= n) return *it;
		}
		return val.back();
	}

};

//分子中原子节点
class MNode {
private:
	int seq;
	Atom e; //元素原子
	bool isa;//是否芳香
	int ccount; //连接数(等效非氢键数，双键+2，三键+3)
	int nhbcount; //非氢键数
	int chcount; //连接氢数
	int minhcount;//最少氢数（来自SMILES格式限定）
	bool csym; //电荷符号
	int cval; //电荷绝对值
	unsigned int rank;
	int mtype;

public:
	MNode(int sequence, Atom elem, bool isaroma, int connectcount = 0, int nonhygconcount = 0, int conhygcount = 0, int chargesym = 0, int chargeval = 0, int mmfftype = -1) :
		seq(sequence), e(elem), isa(isaroma), ccount(connectcount + isaroma), nhbcount(nonhygconcount), csym(chargesym), cval(chargeval), chcount(conhygcount), mtype(mmfftype), rank(0), minhcount(0) {
	}
	~MNode() {}
	void AddNonHydBond(string bondsymbol)
	{
		int add = 0;
		if (bondsymbol == SINGLE_BOND || bondsymbol == AROMA_BOND) add = 1;
		else if (bondsymbol == DOUBLE_BOND) add = 2;
		else if (bondsymbol == TRIPLE_BOND) add = 3;
		nhbcount++;
		ccount += add;
		CalConHyd();
	}
	void AddHydAtom(int n)
	{
		minhcount += n;
	}
	int CalConHyd() {
		chcount = max(e.GetHigherVal(ccount + minhcount) - ccount - (2 * csym - 1) * cval, minhcount);
		//cout << "chc\tcc\tnhbc\n";
		//cout << chcount << "\t" << ccount << "\t"<< nhbcount<<endl;
		return chcount;
	}
	bool IsAroma() { return isa; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)
	{
		if (title_state) {
			PrintTableTitle(MNODE_TITLE);
			cout << endl;
		}
		cout << seq << sep << e.GetSym() << sep << isa << sep << ccount << sep << chcount << sep << nhbcount << sep << csym << sep << cval << sep << mtype << endl;
	}
	unsigned long int GetOriRank()
	{
		unsigned int r = 0;
		r = ccount;
		r = r * pow(10, NH_UNIT) + nhbcount;
		r = r * pow(10, AN_UNIT) + e.GetNum();
		r = r * pow(10, CS_UNIT) + csym;
		r = r * pow(10, CV_UNIT) + cval;
		r = r * pow(10, CH_UNIT) + chcount;
		rank = r;
		return r;
	}
	string GetSym()const
	{
		return e.GetSym();
	}
	//unsigned int GetRank() { return rank; }
	//void SetRank(unsigned int r)
	//{
	//	rank = r;
	//}
	friend ostream& operator<<(ostream& out, const MNode& m) {
		out << m.seq << "," << m.e.GetSym() << "," << m.isa << "," << m.ccount << "," << m.chcount << "," << m.nhbcount << "," << m.csym << "," << m.cval << "," << m.mtype << "," << endl;
		return out;
	}
};

//邻接表单行内单节点
class PointTo {
private:
	string bondsym;
	int dseq;
public:
	PointTo(string bondsymbol, int sequence) :bondsym(bondsymbol), dseq(sequence) {}
	~PointTo() {}
	void Print()const { cout << "-" << bondsym << "->" << dseq; }
	int GetDesSeq()const { return dseq; }
	string GetBondSymbol()const { return bondsym; }
	friend ostream& operator<<(ostream& out, const PointTo& p) {
		out << p.dseq << "," << p.bondsym;
		return out;
	}
};

//邻接表单行
class NodeBonds
{
	int nseq;
	vector<PointTo> d;
public:
	NodeBonds(int nodeseq) :nseq(nodeseq) { d.clear(); }
	~NodeBonds() { d.clear(); }
	void AddBond(PointTo bond)
	{
		d.push_back(bond);
	}
	int CountBound()
	{
		return d.size();
	}
	void Print()const
	{
		cout << nseq << ": ";
		for (auto& a : d)
		{
			a.Print();
			cout << " ";
		}
	}
	vector<PointTo> GetBonds()const
	{
		vector<PointTo> p(d);
		return p;
	}
	friend ostream& operator<<(ostream& out, const NodeBonds& n) {
		for (auto& a : n.d)
		{
			out << n.nseq << "," << a << endl;
		}
		return out;
	}
};

struct PreNode {
	int frontseq;
	string str;
	int startpos = 0;
};

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
}
void ProcessAtom(vector<MNode>& ntb, vector<NodeBonds>& btb, const vector<Atom>& etb, const string& c)
{
	if (c.empty()) return;
	int startpos = 0, endpos = c.length() - 1;
	bool chargesym = 0, chargeyes = 0;
	int chargeval = 0, chcount = 0;
	if (IsNumber(c[endpos])) { chargeval = int(c[endpos] - '0'); endpos--; chargeyes = true; }
	if (c[endpos] == '+') { chargesym = 1; chargeyes = true; }
	if (c[endpos] == '-') { chargesym = 0; chargeyes = true; }
	if (chargeyes && !chargeval)
	{
		while (c[endpos] == '+' || c[endpos] == '-')
		{
			if (c[endpos] == '+')chargeval++;
			else chargeval--;
			endpos--;
		}
	}
	if (IsNumber(c[endpos])) { chcount += int(c[endpos] - '0'); endpos--; }
	int hcount = 0;
	bool hexist = false;
	if (c[endpos] == 'H')
	{
		endpos--;
		hexist = true;
	}
	if (IsNumber(c[endpos]))
	{
		endpos--;
		hcount = int(c[endpos] - '0');
	}
	string esym = c.substr(startpos, endpos - startpos + 1);
	Atom e(SearchAtom(etb, esym));
	//cout << "searCH: " << c.substr(startpos, endpos - startpos+1) << endl;
	//etb[0].Print();
	//cout << "Elem: ";
	//e.Print();
	//cout << endl;
	int seq = ntb.size();
	bool isaroma = false;
	if (IsLowerLetter(esym[0])) isaroma = true;
	ntb.push_back(MNode(seq, e, isaroma, 0, 0, chcount, chargesym, chargeval));
	if (hexist) ntb[seq].AddHydAtom(hcount);
	btb.push_back(NodeBonds(seq));
	//ntb[seq].CalConHyd();
}
//sequence从0开始，-1为空
int CircleMatch(int*& table, int nodeseq, int circleseq)
{
	if (table[circleseq] == -1) { table[circleseq] = nodeseq; return -1; }
	else { return table[circleseq]; }
}

class BinRank {
private:
	unsigned long int accordval;
	int seq;
public:
	BinRank(const BinRank& bk) :accordval(bk.accordval), seq(bk.seq) {}
	BinRank(unsigned long int accordvalue, int nowseq) :accordval(accordvalue), seq(nowseq) {}
	unsigned long int GetAccordValue()const { return accordval; }
	int GetSeq()const { return seq; }

};
class PrimeNumber {
private:
	unsigned long int pn;
public:
	PrimeNumber(int n = 0) :pn(n) {}
	PrimeNumber(vector<int> v) { if (!v.empty()) pn = v[0]; else pn = 0; }
	PrimeNumber(vector<string> vs) { if (!vs.empty()) pn = atoi(vs[0].c_str()); else pn = 0; }
	PrimeNumber(vector<string> vs, vector<int> seq) :PrimeNumber(vs) {}
	~PrimeNumber() {}
	void Print()const
	{
		cout << pn << endl;
	}
	int GetPN()const { return pn; }
};
bool SortRankCmp(BinRank a, BinRank b)
{
	return a.GetAccordValue() < b.GetAccordValue();
}

class TreeNode {
private:
	string bondsym;
	int accordval;
	int seq;
	bool isfc;
public:
	TreeNode(string bondsymbol, int accordvalue, int sequence, bool isfirstchild = false) :bondsym(bondsymbol), accordval(accordvalue), seq(sequence), isfc(isfirstchild) {}
	TreeNode(const TreeNode& tn) :bondsym(tn.bondsym), accordval(tn.accordval), seq(tn.seq), isfc(tn.isfc) {}
	string GetBondSymbol()const { return bondsym; }
	int GetAccordValue()const { return accordval; }
	int GetSeq()const { return seq; }
	bool IsFirstChild()const { return isfc; }
	void SetFirstChild() { isfc = true; }
};

bool SortTreeNodeCmp(TreeNode a, TreeNode b)
{
	return a.GetAccordValue() < b.GetAccordValue();
}

void ProcessCircleNumber(vector<string>& s, vector<vector<int>>& v, const vector<int>& findcounttb)
{
	vector<string> str(v.size(), "");
	int circount = 1;
	vector<int> newcirnum(v.size(), -1);
	for (int i = 0; i < findcounttb.size(); i++)
	{
		int seq = findcounttb[i];
		sort(v[seq].begin(), v[seq].end());
		for (auto& a : v[seq])
		{
			if (newcirnum[a] == -1)
			{
				newcirnum[a] = circount;
				circount++;
			}
			int b = newcirnum[a];
			if (b > 9) str[seq] += "%";
			str[seq] += to_string(b);
		}
	}
	s = str;
	//vector<string> str(v.size(), "");
	//for (int i = 0; i < v.size(); i++)
	//{
	//	sort(v[i].begin(), v[i].end());
	//	for (auto& a : v[i])
	//	{
	//		if (a > 9) str[i] += "%";
	//		str[i] += to_string(a);
	//	}
	//}
	//s = str;
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

class Mole {
private:
	vector<MNode> nodetb;
	vector<NodeBonds> bondtb;
	vector<int> ranktb;//r[i]表示nodetb中第i个元素的排序
	vector<bool> sortedtb;//表示节点是否有序
	string com_smiles;
	string can_smiles;
public:
	Mole(const vector<Atom>& etb, const string& smiles) :com_smiles(smiles)
	{
		vector<PreNode> p;
		int* rec = new int[1000];
		memset(rec, -1, 1000);
		p.push_back({ -1, smiles });
		string bondsym = "";
		int count = 0;
		while (!p.empty())
		{
			//count++;
			//if (count > 10) break;
			PreNode& ptop = p.back();

			//cout << "Round: " << ptop.startpos << "/" << ptop.str.length() << endl;
			if (ptop.startpos < ptop.str.length())
			{
				const char nowch = ptop.str[ptop.startpos];
				const string nowstr = ptop.str;
				if (nowch == '[' || IsLetter(nowch))
				{
					int endpos;
					if (nowch == '[') { endpos = MatchPunc(nowch, ptop.str, ptop.startpos); ptop.startpos++; endpos--; }
					if (IsLetter(nowch)) endpos = FindSingleAtom(ptop.str, ptop.startpos);
					//cout << ptop.startpos << " " << endpos << " " << ptop.str.substr(ptop.startpos, endpos - ptop.startpos + 1) << endl;
					ProcessAtom(nodetb, bondtb, etb, ptop.str.substr(ptop.startpos, endpos - ptop.startpos + 1));
					//cout << "size1: " << nodetb.size() << endl;
					int aindex = ptop.frontseq, bindex = nodetb.size() - 1;
					if (aindex != -1)
					{
						if (nodetb[aindex].IsAroma() && nodetb[bindex].IsAroma()) bondsym = AROMA_BOND[0];
						else if (bondsym == "") bondsym = SINGLE_BOND;

						bondtb[aindex].AddBond({ bondsym,bindex });
						bondtb[bindex].AddBond({ bondsym,aindex });

						nodetb[aindex].AddNonHydBond(bondsym);
						nodetb[bindex].AddNonHydBond(bondsym);

						bondsym = SINGLE_BOND;
					}
					if (nowch == '[') endpos++;
					ptop.startpos = endpos + 1;
					ptop.frontseq = nodetb.size() - 1;
				}
				else if (nowch == '(')
				{
					int endpos = MatchPunc(nowch, ptop.str, ptop.startpos); ptop.startpos++;
					int startpos = ptop.startpos; ptop.startpos = endpos + 1;
					p.push_back(PreNode({ ptop.frontseq,ptop.str.substr(startpos,endpos - startpos),0 }));
				}
				else if (IsCircleNumber(nowch))
				{
					int endpos = FindLastCircleNumber(ptop.str, ptop.startpos);
					vector<int> cnum;
					string cstr = ptop.str.substr(ptop.startpos, endpos - ptop.startpos + 1);
					for (int i = 0; i < cstr.length(); i++)
					{
						cnum.push_back(cstr[i] - '0');
					}
					//SplitToInt(cnum, ptop.str.substr(ptop.startpos, endpos - ptop.startpos + 1), '%');
					int aindex = ptop.frontseq;
					for (int i = 0; i < cnum.size(); i++)
					{
						int bindex = CircleMatch(rec, ptop.frontseq, cnum[i]);
						//cout << "nseq: " << ptop.frontseq << endl << "cnum: " << cnum[i] << endl << "table: " << rec[cnum[i]]<<endl << "bindex: " << bindex << endl;
						if (bindex != -1)
						{
							string cbsym = SINGLE_BOND;
							if (nodetb[aindex].IsAroma() && nodetb[bindex].IsAroma()) cbsym = AROMA_BOND;

							bondtb[aindex].AddBond({ cbsym,bindex });
							bondtb[bindex].AddBond({ cbsym,aindex });

							nodetb[aindex].AddNonHydBond(cbsym);
							nodetb[bindex].AddNonHydBond(cbsym);
							//cout << "Match!\n";
						}
					}
					ptop.startpos = endpos + 1;
				}
				else if (nowch == TRIPLE_BOND[0] || nowch == DOUBLE_BOND[0])
				{
					bondsym = nowch;
					ptop.startpos++;
				}
				else {
					ptop.startpos++;
				}
			}
			else
			{
				p.pop_back();
				//cout << "Pop!\n";
			}
		}
		delete[]rec;
		//cout << "size2: " << nodetb.size() << endl;
		vector<bool> stb(nodetb.size(), 0);
		vector<int> rtb(nodetb.size(), 0);
		sortedtb = stb;
		ranktb = rtb;
	}

	void PrintOriRank()
	{
		int i = 0;
		for (auto& a : nodetb)
		{
			//vector<PointTo> b(a.GetBonds());
			if (i != 0) cout << "-";
			cout << a.GetOriRank();
			i++;
		}
		cout << endl;
	}
	void PrintNodeTable()
	{
		PrintTableTitle(MNODE_TITLE);
		//cout << "size3: " << nodetb.size() << endl;
		for (auto& a : nodetb)
		{
			a.Print();
		}
	}
	void PrintBondTable()
	{
		cout << "Seq\tBond\n";
		int i = 0;
		for (auto& a : bondtb)
		{
			cout << i << "\t";
			a.Print();
			cout << endl;
			i++;
		}
	}
	void PrintRank()
	{
		//cout << "size: " << ranktb.size() << endl;
		int i = 0;
		for (auto& a : ranktb)
		{
			//vector<PointTo> b(a.GetBonds());
			if (i != 0) cout << "-";
			cout << a;
			i++;
		}
		cout << endl;
	}
	void PrintIsSorted()
	{
		//cout << "size: " << ranktb.size() << endl;
		int i = 0;
		for (auto& a : sortedtb)
		{
			//vector<PointTo> b(a.GetBonds());
			if (i != 0) cout << "-";
			cout << a;
			i++;
		}
		cout << endl;
	}
	void PrintAtomTable() {
		int i = 0;
		for (auto& a : nodetb)
		{
			if (i != 0) cout << "-";
			cout << a.GetSym();
			i++;
		}
		cout << endl;
	}

	void GetOriRankTable(vector<BinRank>& accordtb)
	{
		accordtb.clear();
		int i = 0;
		for (auto& a : nodetb)
		{
			accordtb.push_back(BinRank(a.GetOriRank(), i));
			i++;
		}
	}
	bool SortByBinRank(vector<BinRank>& accordtb, int startseq = 0)
	{
		//cout << "Before: \n";
		//PrintRank();
		//PrintIsSorted();

		//vector<int> rtb(nodetb.size(), 0);
		sort(accordtb.begin(), accordtb.end(), SortRankCmp);
		int repeatcount = 0;
		int oldrseq = -1;
		for (int j = 0; j < accordtb.size(); j++)
		{
			int nowseq = accordtb[j].GetSeq();
			if (!sortedtb[nowseq])
			{
				if (j > 0 && accordtb[j].GetAccordValue() == accordtb[j - 1].GetAccordValue())
				{
					ranktb[nowseq] = oldrseq; repeatcount++;
					if (j == 1) repeatcount++;
				}
				else {
					ranktb[nowseq] = j + startseq; oldrseq = j + startseq;
					if (j < accordtb.size() - 1 && accordtb[j].GetAccordValue() != accordtb[j + 1].GetAccordValue() || j == accordtb.size() - 1)
						sortedtb[nowseq] = true;
				}
			}
		}
		//cout << "rtb_size: " << rtb.size() << endl;
		//ranktb = rtb;
		cout << "After: ";

		cout << "(Repeat:" << repeatcount << "/" << accordtb.size() << ")" << endl;

		if (repeatcount == accordtb.size())
		{
			//cout << "破同分!\n";
			int nowseq = accordtb.back().GetSeq();
			ranktb[nowseq] = accordtb.size() + startseq - 1;
			sortedtb[nowseq] = true;
			PrintRank();
			PrintIsSorted();

			return false;
		}
		PrintRank();
		PrintIsSorted();

		if (repeatcount == 0) return true;
		else return false;
	}
	int ComplexSortByRank(vector<PrimeNumber>& ptb)
	{
		vector<BinRank> srtb;
		int recordval = -1;
		bool isrecord = false;
		for (int i = 0; i < sortedtb.size(); i++)
		{
			if (!sortedtb[i] && !isrecord)
			{
				recordval = ranktb[i];
				isrecord = true;
			}
			if (isrecord && ranktb[i] == recordval)
			{
				unsigned long int sum = 1;
				vector<PointTo> b(bondtb[i].GetBonds());
				//cout << "NearPN: ";
				for (auto& a : b)
				{
					int dseq = a.GetDesSeq();
					sum *= ptb[ranktb[dseq]].GetPN();
					//cout << ptb[ranktb[dseq]].GetPN() << " ";
				}
				//cout << "Sum: " << sum;
				//cout << endl;
				srtb.push_back({ sum,i });
			}
		}
		if (!isrecord) return true;
		SortByBinRank(srtb, recordval);
		return false;
	}
	void MoleSortWithPN(vector<PrimeNumber>& ptb)
	{
		vector<BinRank> oritb;
		GetOriRankTable(oritb);
		unsigned long int maxcount = nodetb.size() * log2(nodetb.size()), count = 1;
		cout << "第" << count << "次排序：\n";
		if (SortByBinRank(oritb)) return;
		cout << "\n第" << count + 1 << "次排序：\n";
		while (count < maxcount && !ComplexSortByRank(ptb)) {
			count++;
			cout << "\n第" << count + 1 << "次排序：\n";
		}
		cout << "\n总排序次数：" << count << endl;
	}
	string GenerateCanSmiles()
	{
		int rootseq;
		for (rootseq = 0; rootseq < ranktb.size(); rootseq++)
			if (ranktb[rootseq] == 0)
				break;

		string cs = "";

		int circlecount = 0;
		vector<vector<int>> circleint(nodetb.size());
		vector<string> circlestr;
		vector<int> findcounttb(ranktb.size(), -1);
		//vector<bool> endtb(ranktb.size(), false);
		//一共两次遍历，第一次通过DFS遍历记录环的标准断开位置，第二次遍历生成can_smiles
		for (int i = 0; i < 2; i++)
		{
			if (i == 1)ProcessCircleNumber(circlestr, circleint,findcounttb);
			vector<TreeNode> treetb;
			treetb.push_back(TreeNode("-", 0, rootseq, true));

			vector<bool> recordtb(ranktb.size(), false);//用于记录节点是否已被遍历
			vector<int> parenttb(ranktb.size(), -1);//用于防止环的误判
			//recordtb[rootseq] = true;//

			int leftbracketcount = 0;
			
			//int oldseq = -1;
			int findcount = 0;
			while (!treetb.empty())
			{
				TreeNode top(treetb.back());
				int topseq = top.GetSeq();
				treetb.pop_back();

				if (recordtb[topseq]) continue;
				recordtb[topseq] = true;
				findcounttb[findcount] = topseq;
				findcount++;

				vector<PointTo> ptb(bondtb[topseq].GetBonds());
				vector<TreeNode> drtb;
				//int count = 0;

				for (auto& a : ptb)
				{
					int dseq = a.GetDesSeq();
					if (!recordtb[dseq])
					{
						parenttb[dseq] = topseq;
						drtb.push_back(TreeNode(a.GetBondSymbol(), ranktb[dseq], dseq));
					}
					else if (i == 0 && parenttb[topseq] != dseq)
					{
						circlecount++;
						circleint[topseq].push_back(circlecount);
						circleint[dseq].push_back(circlecount);
					}
				}
				//oldseq = topseq;
				//if (count) cout << endl;
				if (i == 1)
				{
					if (!top.IsFirstChild()) { cs += "("; leftbracketcount++; };
					//cout << top.GetBondSymbol() << endl;
					if (top.GetBondSymbol() == DOUBLE_BOND) cs += DOUBLE_BOND;
					else if (top.GetBondSymbol() == TRIPLE_BOND) cs += TRIPLE_BOND;
					if (nodetb[topseq].IsAroma()) cs += StringLower(nodetb[topseq].GetSym());
					else cs += nodetb[topseq].GetSym();

					cs += circlestr[topseq];
				}
				if (drtb.empty())
				{
					//cout << "Empty! LeftBracketCount: "<< leftbracketcount<<endl;
					//if (i == 1) endtb[topseq] = true;
					if (leftbracketcount != 0 && i == 1) { cs += ")"; leftbracketcount--; }
					continue;
				}
				sort(drtb.begin(), drtb.end(), SortTreeNodeCmp);
				drtb.back().SetFirstChild();//秩最高为主链
				//cout << "FirstChild: " << drtb.back().IsFirstChild() << endl;
				for (int j = drtb.size() - 1; j >= 0; j--)
				{
					treetb.push_back(drtb[j]);//先遍历秩最低
				}
				//cout <<"左括号数：" << leftbracketcount << endl;
			}
			if (i == 1 && leftbracketcount != 0) cs += string(leftbracketcount, ')');
		}

		can_smiles = DeleteExtraBracket(cs);

		return can_smiles;
	}
	string GetComSmiles()const { return com_smiles; }
	string GetCanSmiles()
	{
		if (empty(can_smiles)) can_smiles = GenerateCanSmiles();
		return can_smiles;
	}
	vector <MNode> GetNodeTable()const {
		return nodetb;
	}
	vector <NodeBonds> GetBondTable()const {
		return bondtb;
	}
};


bool AtomTableCmp(Atom a, Atom b)
{
	return a.GetNum() < b.GetNum();
}

//template<typename T>
//T StringTo(string str){return T(str) }

template<typename T>
int ReadTable(string filename, string stdtitle, vector<T>& table, bool readtitle = YES)
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
			} while (str != "\n" && str != "" && colsize < 4);
			//cout << "all colsize:" << colsize << endl;
			//cout << "表格内容：\n";
			if (readtitle)
			{
				GetTitleSortSeq(seq, stdtitle, tabletitle);
				continue;
			}
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
		if (readtitle)
		{
			T info(stemp, seq);
			table.push_back(info);
			//info.Print();
			//cout << endl;
		}
		else
		{
			T info(stemp);
			table.push_back(info);
			//cout << stemp[0] << endl;
		}
		k++;
	}
	ifs.close();
	item.clear();
	cout << "文件" << filename << "加载完毕！\n共有 " << k << " 行 " << colsize << " 列记录!" << endl;
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
void PrintCmdSepTitle(const string title, int sepwidth = SEP_WIDTH, const char fillsym = SEP_SYMBOL)
{
	int scount = max((sepwidth - title.length() - 6) / 2, 0);
	string sidesep(scount, fillsym);
	string sum = "\n" + sidesep + " * " + title + " * " + sidesep + "\n\n";
	cout << sum;
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
int main()
{

	//string a = "2/4/6";
	//vector<int> v;
	//SplitToInt(v, a, '/');
	//PrintVec<int>(v, " ");

	//string a = "abcdaa";
	//string b = "bcef1c31s2d1`23";
	//cout << LCS(a, b);

	string etb_filename = "File/Tables/ElementsTable.csv";
	string ptb_filename = "File/Tables/PrimeNumber1000.csv";
	string output_folder = "File/Output";

	PrintCmdSepTitle("元素周期表读取");

	vector<Atom> atomtable;
	ReadTable(etb_filename, ATOM_TITLE, atomtable, true);
	sort(atomtable.begin(), atomtable.end(), AtomTableCmp);

	PrintTableTitle(ATOM_TITLE);
	PrintSpecialVector(atomtable);

	PrintCmdSepTitle("素数表读取");

	vector<PrimeNumber> primetable;
	ReadTable(ptb_filename, "PrimeNumber", primetable, false);


	//PrintTableTitle("PrimeNumber");
	//PrintVector(primetable);

	//cout << endl;
	//for (int i = 0; i < atomtable.size(); i++)
	//{
	//	atomtable[i].Print();
	//	cout << endl;
	//}
	//cout << "ANode: \n";
	//ANode a(atomtable[0],1,1,0,0,3);
	//Atom a(SearchAtom(atomtable, "S", NO));
	//cout << endl;
	//a.Print();
	int turncount = 0;
	while (1)
	{
		turncount++;
		PrintCmdSepTitle("第" + to_string(turncount) + "轮结构式解析");

		cout << "请输入普通SMILES结构式：" << endl;
		string smiles;
		cin >> smiles;

		LARGE_INTEGER frequency;
		QueryPerformanceFrequency(&frequency);
		LARGE_INTEGER start, stop;
		QueryPerformanceCounter(&start);
		//------------------------------**计时开始**------------------------------

		PrintCmdSepTitle("分子节点提取");

		Mole c(atomtable, smiles);
		c.PrintNodeTable();

		PrintCmdSepTitle("分子节点邻接表");
		c.PrintBondTable();

		PrintCmdSepTitle("初始秩生成");
		c.PrintAtomTable();
		c.PrintOriRank();

		PrintCmdSepTitle("秩排序过程");
		c.MoleSortWithPN(primetable);
		cout << "\n排序结果：\n";
		c.PrintRank();

		PrintCmdSepTitle("唯一SMILES生成");

		cout << "原SMILES：\n";
		cout << c.GetComSmiles() << endl << endl;
		cout << "唯一SMILES：\n";
		cout << c.GenerateCanSmiles() << endl;

		//------------------------------**计时结束**------------------------------
		QueryPerformanceCounter(&stop);

		double duration = (stop.QuadPart - start.QuadPart) * 1000.0 / frequency.QuadPart;
		std::cout << "\n(^-^)Time taken by function: " << duration << " ms" << std::endl;

		cout << "是否要储存节点信息表及邻接表？(Y/N)";
		char savechoice;
		cin >> savechoice;
		if (savechoice == 'Y' || savechoice == 'y')
		{
			string cs = c.GetCanSmiles();
			Mole d(atomtable, cs);
			const string folderpath = output_folder + "/" + cs;
			if (!CreateFolder(folderpath)) {
				cout << "文件夹创建失败！可能已存在。\n";
				char savechoice2;
				cin >> savechoice2;
				cout << "是否要继续覆盖节点信息表及邻接表？(Y/N)";
				if (!(savechoice2 == 'Y' || savechoice2 == 'y'))
				{
					continue;
				}
			}
			string now_time = GetCurrentTimeString();
			string atomtable_filename = folderpath + "/AtomTable.csv";
			string bondtable_filename = folderpath + "/BondTable.csv";
			PrintCmdSepTitle("节点表储存");
			WriteTable<MNode>(atomtable_filename, MNODE_TITLE, d.GetNodeTable(), true);
			cout << "节点表已储存至 " << atomtable_filename << endl;
			PrintCmdSepTitle("邻接表储存");
			WriteTable<NodeBonds>(bondtable_filename, ADJ_TITLE, d.GetBondTable(), true);
			cout << "邻接表已储存至 " << bondtable_filename << endl;
		}
	}
	return 0;
}