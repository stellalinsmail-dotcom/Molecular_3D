#ifndef MOLECULE_H
#define MOLECULE_H

#include "types.h"
#include "Atom.h"

// 前向声明
class PointTo;
class NodeBonds;
class MNode;

// --- 辅助结构和类 ---

struct PreNode {
	int frontseq;
	string str;
	int startpos = 0;
};

class BinRank {
private:
	unsigned long int accordval;
	int seq;
public:
	BinRank(const BinRank& bk);
	BinRank(unsigned long int accordvalue, int nowseq);
	unsigned long int GetAccordValue()const { return accordval; }
	int GetSeq()const { return seq; }

};
bool SortRankCmp(BinRank a, BinRank b);

class PrimeNumber {
private:
	unsigned long int pn;
public:
	PrimeNumber(int n = 0);
	PrimeNumber(vector<int> v);
	PrimeNumber(vector<string> vs);
	PrimeNumber(vector<string> vs, vector<int> seq);
	~PrimeNumber() {}
	void Print()const;
	int GetPN()const { return pn; }
};

class TreeNode {
private:
	string bondsym;
	int accordval;
	int seq;
	bool isfc;
public:
	TreeNode(string bondsymbol, int accordvalue, int sequence, bool isfirstchild = false);
	TreeNode(const TreeNode& tn);
	string GetBondSymbol()const { return bondsym; }
	int GetAccordValue()const { return accordval; }
	int GetSeq()const { return seq; }
	bool IsFirstChild()const { return isfc; }
	void SetFirstChild() { isfc = true; }
};
bool SortTreeNodeCmp(TreeNode a, TreeNode b);

// --- 分子结构类 ---

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
	MNode(int sequence, Atom elem, bool isaroma, int connectcount = 0, int nonhygconcount = 0, int conhygcount = 0, int chargesym = 0, int chargeval = 0, int mmfftype = -1);
	MNode(const vector<string>& info, const vector<int>& titlenum);
	
	~MNode() {}
	void AddNonHydBond(string bondsymbol);
	void AddHydAtom(int n);
	int CalConHyd();
	bool IsAroma() const{ return isa; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const;
	unsigned long int GetOriRank();
	string GetSym()const {return e.GetSym();}
	int GetCHCount()const { return  chcount; }

	friend ostream& operator<<(ostream& out, const MNode& m);
};

//邻接表单行内单节点
class PointTo {
private:
	string bondsym;
	int dseq;
public:
	PointTo(string bondsymbol, int sequence);
	~PointTo() {}
	void Print()const { cout << "-" << bondsym << "->" << dseq; }
	int GetDesSeq()const { return dseq; }
	string GetBondSymbol()const { return bondsym; }
	friend ostream& operator<<(ostream& out, const PointTo& p);
};

class AdjLine
{
private:
	int seq_i;
	int seq_j;
	string bondsym;
public:
	AdjLine(const vector<string>& info, const vector<int>& titlenum);
	int GetSeqI()const { return seq_i; }
	int GetSeqJ()const { return seq_j; }
	string GetBondSym()const { return bondsym; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const;
};


//邻接表单行
class NodeBonds
{
	int nseq;
	vector<PointTo> d;
public:
	NodeBonds(int nodeseq);
	~NodeBonds() { d.clear(); }
	void AddBond(PointTo bond);
	int CountBound();
	void Print(string sep = "\t")const;
	vector<PointTo> GetBonds()const;
	friend ostream& operator<<(ostream& out, const NodeBonds& n);
};

// --- 分子核心类 ---

class Mole {
private:
	vector<MNode> nodetb;
	vector<NodeBonds> bondtb;
	vector<int> ranktb;//r[i]表示nodetb中第i个元素的排序
	vector<bool> sortedtb;//表示节点是否有序
	string com_smiles;
	string can_smiles;

	void ProcessAtom(vector<MNode>& ntb, vector<NodeBonds>& btb, const vector<Atom>& etb, const string& c);
	int CircleMatch(int*& table, int nodeseq, int circleseq);
	void ProcessCircleNumber(vector<string>& s, vector<vector<int>>& v, const vector<int>& findcounttb);

public:
	Mole(const vector<Atom>& etb, const string& smiles);
	Mole(const vector<MNode>& now_nodetb, const  vector<NodeBonds>& bondtb);

	void PrintOriRank();
	void PrintNodeTable();
	void PrintBondTable();
	void PrintRank();
	void PrintIsSorted();
	void PrintAtomTable();

	void GetOriRankTable(vector<BinRank>& accordtb);
	bool SortByBinRank(vector<BinRank>& accordtb, int startseq = 0, bool printyes = false);
	int ComplexSortByRank(vector<PrimeNumber>& ptb);
	void MoleSortWithPN(vector<PrimeNumber>& ptb, bool printyes = false);
	string GenerateCanSmiles();

	string GetComSmiles()const { return com_smiles; }
	string GetCanSmiles();
	vector <MNode> GetNodeTable()const {
		return nodetb;
	}
	vector <NodeBonds> GetBondTable()const {
		return bondtb;
	}
};
// --- 辅助函数 ---
vector<NodeBonds> AdjTbToBondTb(const vector<AdjLine>& ,const vector<MNode>&, bool is_short = YES);

string DeleteHydrogen(string smiles);
#endif // MOLECULE_H