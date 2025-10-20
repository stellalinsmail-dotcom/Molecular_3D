#ifndef ATOM_H
#define ATOM_H

#include "types.h"

class Atom {
private:
	string sym; //原子符号
	int num; //原子序数
	int mass; //相对原子质量
	vector<int> val; //原子常见化合价
public:
	Atom(string symbol = "None");
	Atom(const Atom& e);
	Atom(string symbol, int number, int relamass, vector<int> valence);
	Atom(const vector<string>& atom_info);
	Atom(const vector<string>& atom_info, const vector<int>& titlenum);
	~Atom();

	int GetNum()const { return num; }
	string GetSym()const { return sym; }
	void Print(string sep = "\t", bool title_state = TITLE_OFF)const;
	int GetHigherVal(int n);

	friend Atom SearchAtom(const vector<Atom>& atb, string symbol, bool isprecise = NO);
};

bool AtomTableCmp(Atom a, Atom b);


#endif // ATOM_H