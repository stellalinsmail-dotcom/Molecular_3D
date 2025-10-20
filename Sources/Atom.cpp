#include "Atom.h"
#include "types.h" // 确保能访问到辅助函数

Atom::Atom(string symbol) :sym(symbol), num(0), mass(0), val(0, 0) {}
Atom::Atom(const Atom& e) :sym(e.sym), num(e.num), mass(e.mass), val(e.val) {}
Atom::Atom(string symbol, int number, int relamass, vector<int> valence) :num(number), mass(relamass), sym(symbol) {
	val.assign(valence.begin(), valence.end());
}
Atom::Atom(const vector<string>& atom_info)
{
	sym = atom_info[0];
	num = atoi(atom_info[1].c_str());
	mass = atoi(atom_info[2].c_str());
	SplitToInt(val, atom_info[3], VAL_SEP);
	sort(val.begin(), val.end());
}
Atom::Atom(const vector<string>& atom_info, const vector<int>& titlenum)
{
	sym = atom_info[titlenum[0]];
	num = atoi(atom_info[titlenum[1]].c_str());
	mass = atoi(atom_info[titlenum[2]].c_str());
	SplitToInt(val, atom_info[titlenum[3]], VAL_SEP);
	sort(val.begin(), val.end());
}
Atom::~Atom() { val.clear(); }

void Atom::Print(string sep, bool title_state)const
{
	if (title_state) {
		PrintTableTitle(ATOM_TB_TITLE);
		cout << endl;
		//cout << "symbol" << sep << "number" << sep << "rmass" << sep << "valence" << endl;
	}
	cout << sym << sep << num << sep << mass << sep;
	PrintCommonVector(val, " ");
	cout << endl;
}

int Atom::GetHigherVal(int n)
{
	if (val.empty()) return 0;
	for (auto it = val.begin(); it != val.end(); it++)
	{
		if (*it >= n) return *it;
	}
	return val.back();
}

Atom SearchAtom(const vector<Atom>& atb, string symbol, bool isprecise)
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

bool AtomTableCmp(Atom a, Atom b)
{
	return a.GetNum() < b.GetNum();
}