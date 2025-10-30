#include "Molecule.h"
#include "types.h" // 确保能访问到辅助函数和常量

// --- MNode 实现 ---
MNode::MNode(int sequence, Atom elem, bool isaroma, int connectcount, int nonhygconcount, int conhygcount, int chargesym, int chargeval, int mmfftype) :
	seq(sequence), e(elem), isa(isaroma), ccount(connectcount + isaroma), nhbcount(nonhygconcount), csym(chargesym), cval(chargeval), chcount(conhygcount), mtype(mmfftype), rank(0), minhcount(0) {
}
MNode::MNode(const vector<string>& info, const vector<int>& titlenum)
{
	seq = atoi(info[titlenum[0]].c_str());
	e = Atom(info[titlenum[1]]);
	isa = (info[titlenum[2]] == "1") ? true : false;
	ccount = atoi(info[titlenum[3]].c_str());
	chcount = atoi(info[titlenum[4]].c_str());
	nhbcount = atoi(info[titlenum[5]].c_str());
	csym = (info[titlenum[6]] == "1") ? true : false;
	cval = atoi(info[titlenum[7]].c_str());
	mtype = atoi(info[titlenum[8]].c_str());
	minhcount = 0;
	rank = 0;
}

void MNode::AddNonHydBond(string bondsymbol)
{
	int add = 0;
	if (bondsymbol == SINGLE_BOND || bondsymbol == AROMA_BOND) add = 1;
	else if (bondsymbol == DOUBLE_BOND) add = 2;
	else if (bondsymbol == TRIPLE_BOND) add = 3;
	nhbcount++;
	ccount += add;
	CalConHyd();
}
void MNode::AddHydAtom(int n)
{
	minhcount += n;
}
int MNode::CalConHyd() {
	chcount = max(e.GetHigherVal(ccount + minhcount) - ccount - (2 * csym - 1) * cval, minhcount);
	//cout << "chc\tcc\tnhbc\n";
	//cout << chcount << "\t" << ccount << "\t"<< nhbcount<<endl;
	return chcount;
}
void MNode::Print(string sep, bool title_state)const
{
	if (title_state) {
		PrintTableTitle(MNODE_TB_TITLE);
		cout << endl;
	}
	cout << seq << sep << e.GetSym() << sep << isa << sep << ccount << sep << chcount << sep << nhbcount << sep << csym << sep << cval << sep << mtype << endl;
}
unsigned long int MNode::GetOriRank()
{
	unsigned int r = 0;
	r = ccount;
	r = r * (unsigned long int)pow(10, NH_UNIT) + nhbcount;
	r = r * (unsigned long int)pow(10, AN_UNIT) + e.GetNum();
	r = r * (unsigned long int)pow(10, CS_UNIT) + csym;
	r = r * (unsigned long int)pow(10, CV_UNIT) + cval;
	r = r * (unsigned long int)pow(10, CH_UNIT) + chcount;
	rank = r;
	return r;
}

ostream& operator<<(ostream& out, const MNode& m) {
	out << m.seq << "," << m.e.GetSym() << "," << m.isa << "," << m.ccount << "," << m.chcount << "," << m.nhbcount << "," << m.csym << "," << m.cval << "," << m.mtype << endl;
	return out;
}

// --- PointTo 实现 ---
PointTo::PointTo(string bondsymbol, int sequence) :bondsym(bondsymbol), dseq(sequence) {}
ostream& operator<<(ostream& out, const PointTo& p) {
	out << p.dseq << "," << p.bondsym;
	return out;
}

// --- AdjLine 实现 ---
AdjLine::AdjLine(const vector<string>& info, const vector<int>& titlenum)
{
	seq_i = atoi(info[titlenum[0]].c_str());
	seq_j = atoi(info[titlenum[1]].c_str());
	bondsym = info[titlenum[2]];
}
void AdjLine::Print(string sep, bool title_state)const
{
	if (title_state) {
		PrintTableTitle(ADJ_TB_TITLE, sep);
		cout << endl;
	}
	cout << seq_i << sep << seq_j << sep << bondsym << endl;
}


// --- NodeBonds 实现 ---
NodeBonds::NodeBonds(int nodeseq) :nseq(nodeseq) { d.clear(); }
void NodeBonds::AddBond(PointTo bond)
{
	d.push_back(bond);
}
int NodeBonds::CountBound()
{
	return d.size();
}
void NodeBonds::Print(string sep)const
{
	cout << nseq << ": " << sep;
	for (auto& a : d)
	{
		a.Print();
		cout << sep;
	}
	cout << endl;
}
vector<PointTo> NodeBonds::GetBonds()const
{
	vector<PointTo> p(d);
	return p;
}
ostream& operator<<(ostream& out, const NodeBonds& n) {
	for (auto& a : n.d)
	{
		out << n.nseq << "," << a << endl;
	}
	return out;
}

// --- PreNode 相关的 SMILES 解析辅助函数（已移入types.cpp/h） ---
// 这些函数 FindSingleAtom, MatchPunc, FindLastCircleNumber 现在应该在 types.cpp/h 中实现。

// --- BinRank 实现 ---
BinRank::BinRank(const BinRank& bk) :accordval(bk.accordval), seq(bk.seq) {}
BinRank::BinRank(unsigned long int accordvalue, int nowseq) :accordval(accordvalue), seq(nowseq) {}
bool SortRankCmp(BinRank a, BinRank b)
{
	return a.GetAccordValue() < b.GetAccordValue();
}

// --- PrimeNumber 实现 ---
PrimeNumber::PrimeNumber(int n) :pn(n) {}
PrimeNumber::PrimeNumber(vector<int> v) { if (!v.empty()) pn = v[0]; else pn = 0; }
PrimeNumber::PrimeNumber(vector<string> vs) { if (!vs.empty()) pn = atoi(vs[0].c_str()); else pn = 0; }
PrimeNumber::PrimeNumber(vector<string> vs, vector<int> seq) :PrimeNumber(vs) {}
void PrimeNumber::Print(const string sep)const
{
	cout << pn << endl;
}

// --- TreeNode 实现 ---
TreeNode::TreeNode(string bondsymbol, int accordvalue, int sequence, bool isfirstchild) :bondsym(bondsymbol), accordval(accordvalue), seq(sequence), isfc(isfirstchild) {}
TreeNode::TreeNode(const TreeNode& tn) :bondsym(tn.bondsym), accordval(tn.accordval), seq(tn.seq), isfc(tn.isfc) {}
bool SortTreeNodeCmp(TreeNode a, TreeNode b)
{
	return a.GetAccordValue() < b.GetAccordValue();
}

// --- Mole 内部实现 ---

void Mole::ProcessAtom(vector<MNode>& ntb, vector<NodeBonds>& btb, const vector<Atom>& etb, const string& c)
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
			if (endpos < 0) return;
		}
	}
	if (IsNumber(c[endpos])) { chcount += int(c[endpos] - '0'); endpos--; }
	int hcount = 0;
	bool hexist = false;
	if (c[endpos] == 'H')
	{
		endpos--;
		hexist = true;
		if (endpos < 0) return;
	}
	if (IsNumber(c[endpos]))
	{
		endpos--;
		hcount = int(c[endpos] - '0');
		if (endpos < 0) return;
	}
	string esym = c.substr(startpos, endpos - startpos + 1);
	Atom e(SearchAtom(etb, esym));

	//e.Print();
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

int Mole::CircleMatch(int*& table, int nodeseq, int circleseq)
{
	if (table[circleseq] == -1) { table[circleseq] = nodeseq; return -1; }
	else { return table[circleseq]; }
}

void Mole::ProcessCircleNumber(vector<string>& s, vector<vector<int>>& v, const vector<int>& findcounttb)
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


// --- Mole 核心实现 ---
Mole::Mole(const vector<Atom>& etb, const string& smiles) :com_smiles(smiles),has_circle(false)
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
				has_circle = true;
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
	if (nodetb.size() == 1) nodetb[0].CalConHyd();
	delete[]rec;
	//cout << "size2: " << nodetb.size() << endl;
	vector<bool> stb(nodetb.size(), 0);
	vector<int> rtb(nodetb.size(), 0);
	sortedtb = stb;
	ranktb = rtb;
	
}
Mole::Mole(const vector<MNode>& now_nodetb, const  vector<NodeBonds>& now_bondtb):
	nodetb(now_nodetb), bondtb(now_bondtb), sortedtb(vector<bool>(nodetb.size(), 0)), ranktb(vector<int>(nodetb.size(), 0)), has_circle(false)
{
}
Mole::Mole(const Mole& m) :
	nodetb(m.nodetb), bondtb(m.bondtb), ranktb(m.ranktb), sortedtb(m.sortedtb), com_smiles(m.com_smiles), can_smiles(m.can_smiles), has_circle(m.has_circle)
{
}
Mole::Mole(const vector<Atom>& etb, const vector<SimpleMNode>& sn_tb, const vector<AdjLine>& adj_list)
{
	nodetb.clear();
	bondtb.clear();
	sortedtb.clear();
	ranktb.clear();
	com_smiles = "";
	can_smiles = "";
	has_circle = false;
	//处理节点表

	for (const auto& sn : sn_tb)
	{
		string atomsym = sn.sym;
		string bondsym = "";
		bool isaroma = false;
		ProcessAtom(nodetb, bondtb, etb, atomsym);
	}
	for (const auto& al : adj_list)
	{
		int aindex = al.GetSeqI();
		int bindex = al.GetSeqJ();
		string bondsym = al.GetBondSym();
		bondtb[aindex].AddBond({ bondsym,bindex });
		bondtb[bindex].AddBond({ bondsym,aindex });
		nodetb[aindex].AddNonHydBond(bondsym);
		nodetb[bindex].AddNonHydBond(bondsym);
	}
	//处理邻接表

	sortedtb = vector<bool>(nodetb.size(), 0);
	ranktb = vector<int>(nodetb.size(), 0);
}

void Mole::PrintOriRank()
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
void Mole::PrintNodeTable()
{
	PrintTableTitle(MNODE_TB_TITLE);
	//cout << "size3: " << nodetb.size() << endl;
	for (auto& a : nodetb)
	{
		a.Print();
	}
}
void Mole::PrintBondTable()
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
void Mole::PrintRank()
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
void Mole::PrintIsSorted()
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
void Mole::PrintAtomTable()
{
	int i = 0;
	for (auto& a : nodetb)
	{
		if (i != 0) cout << "-";
		cout << a.GetSym();
		i++;
	}
	cout << endl;
}

void Mole::GetOriRankTable(vector<BinRank>& accordtb)
{
	accordtb.clear();
	int i = 0;
	for (auto& a : nodetb)
	{
		accordtb.push_back(BinRank(a.GetOriRank(), i));
		i++;
	}
}
bool Mole::SortByBinRank(vector<BinRank>& accordtb, int startseq, bool printyes)
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
	if (printyes) cout << "After: ";

	if (printyes) cout << "(Repeat:" << repeatcount << "/" << accordtb.size() << ")" << endl;

	if (repeatcount == accordtb.size())
	{
		//cout << "破同分!\n";
		int nowseq = accordtb.back().GetSeq();
		ranktb[nowseq] = accordtb.size() + startseq - 1;
		sortedtb[nowseq] = true;
		if (printyes) PrintRank();
		if (printyes) PrintIsSorted();

		return false;
	}
	if (printyes) PrintRank();
	if (printyes) PrintIsSorted();

	if (repeatcount == 0) return true;
	else return false;
}
int Mole::ComplexSortByRank(const vector<PrimeNumber>& ptb)
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
void Mole::MoleSortWithPN(const vector<PrimeNumber>& ptb, bool printyes)
{
	vector<BinRank> oritb;
	GetOriRankTable(oritb);
	unsigned long int maxcount = nodetb.size() * log2(nodetb.size()), count = 1;
	if (printyes) cout << "第" << count << "次排序：\n";
	if (SortByBinRank(oritb, 0, printyes)) return;
	if (printyes)cout << "\n第" << count + 1 << "次排序：\n";
	while (count < maxcount && !ComplexSortByRank(ptb)) {
		count++;
		if (printyes) cout << "\n第" << count + 1 << "次排序：\n";
	}
	if (printyes) cout << "\n总排序次数：" << count << endl;
}
string Mole::GenerateCanSmiles(const vector<PrimeNumber>& primetable)
{
	if (nodetb.size()<=1) return com_smiles;
	MoleSortWithPN(primetable);
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
		if (i == 1)ProcessCircleNumber(circlestr, circleint, findcounttb);
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

			int max_dseq = ptb.back().GetDesSeq();
			for (auto& a : ptb)
			{
				int dseq = a.GetDesSeq();
				if (!recordtb[dseq])
				{
					parenttb[dseq] = topseq;
					int addnum = IsSpecialBond(a.GetBondSymbol());
					if (addnum) addnum += max_dseq;
					drtb.push_back(TreeNode(a.GetBondSymbol(), ranktb[dseq] - addnum,dseq));
				}
				else if (i == 0 && parenttb[topseq] != dseq)
				{
					has_circle = true;
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
string Mole::GetCanSmiles()
{
	//if (can_smiles.empty()) can_smiles = GenerateCanSmiles();
	return can_smiles;
}

// --- 辅助函数实现 ---
vector<AdjLine> BondTbToAdjTb(const vector<NodeBonds> & bond_tb)
{
	vector<AdjLine> adj_tb;
	for (int i=0;i<bond_tb.size();i++)
	{
		vector<PointTo> pto_tb = bond_tb[i].GetBonds();
		for (auto& pto : pto_tb)
		{
			int jseq = pto.GetDesSeq();
			if (i < jseq)
			{
				AdjLine aline(i, jseq, pto.GetBondSymbol());
				adj_tb.push_back(aline);
			}
		}
	}
	return adj_tb;
}

int ExpandBondTb(ADJ_LIST& new_bond_tb, const ADJ_LIST& old_bond_tb, const vector<MNode>& mnode_tb)
{
	new_bond_tb = old_bond_tb;
	int nhc = mnode_tb.size();
	
	int now_hseq = nhc;
	for (int i = 0; i < nhc; i++)
	{
		int chc = mnode_tb[i].GetCHCount();
		//ac += chc;
		for (int j = 0; j < chc; j++)
		{
			PointTo ch("-", now_hseq);
			new_bond_tb[i].AddBond(ch);
			now_hseq++;
		}
		
	}
	return now_hseq;

}
vector<NodeBonds> AdjTbToBondTb(const vector<AdjLine>& adj_tb, const vector<MNode>& mnode_tb, bool is_short)
{
	vector<NodeBonds> bond_tb;
	int nhc = mnode_tb.size();

	int ac = nhc;

	if (!is_short) {
		for (int i = 0; i < nhc; i++)
		{
			ac += mnode_tb[i].GetCHCount();
		}
	}
	bond_tb.resize(ac, NodeBonds(0));
	for (int i = 0; i < ac; i++)
	{
		bond_tb[i] = NodeBonds(i);
	}
	for (auto& line : adj_tb)
	{
		PointTo pto_j(line.GetBondSym(), line.GetSeqJ());
		bond_tb[line.GetSeqI()].AddBond(pto_j);
		//PointTo pto_i(line.GetBondSym(), line.GetSeqI());
		//bond_tb[line.GetSeqJ()].AddBond(pto_i);
	}
	int now_hseq = nhc;
	for (int i = 0; i < nhc; i++)
	{
		int chc = mnode_tb[i].GetCHCount();
		for (int j = 0; j < chc; j++)
		{
			PointTo ch("-", now_hseq);
			bond_tb[i].AddBond(ch);
			if (!is_short)
			{
				PointTo hc("-", i);
				bond_tb[now_hseq].AddBond(hc);
			}
			now_hseq++;
		}
	}
	return bond_tb;
}
string DeleteHydrogen(string smiles)
{
	int count = 0;
	vector<int> delpos;
	for (int i = 0; i < smiles.length(); i++)
	{
		if (smiles[i] == '[') count++;
		if (smiles[i] == ']') count--;
		if (count == 0 && (smiles[i] == 'H' || smiles[i] == 'h')) delpos.push_back(i);
	}
	for (int i = delpos.size() - 1; i >= 0; i--)
	{
		smiles.erase(delpos[i], 1);
	}
	return smiles;
}


int IsSpecialBond(const string& bondsym)
{
	if (bondsym == DOUBLE_BOND|| bondsym == AROMA_BOND) return 2;
	if (bondsym == TRIPLE_BOND) return 3;
	return 0;
}

SmilesFundTable ReadSmilesSolidParam(bool print_yes, int max_row_count)
{
	string etb_filename = "File/Tables/ElementsTable.csv";
	string ptb_filename = "File/Tables/PrimeNumber1000.csv";


	vector<Atom> atomtable;
	vector<PrimeNumber> primetable;

	ReadTableByTitle(etb_filename, ATOM_TB_TITLE, atomtable, print_yes);
	sort(atomtable.begin(), atomtable.end(), AtomTableCmp);

	ReadTable(ptb_filename, "PrimeNumber", primetable, print_yes);


	if (print_yes)
	{
		PrintCmdSepTitle("元素周期表读取");
		PrintTableTitle(ATOM_TB_TITLE);
		PrintSpecialVector(atomtable);

		PrintCmdSepTitle("素数表读取");
		PrintSpecialVector(primetable, "\t", max_row_count);
	}
	SmilesFundTable sft = { atomtable, primetable };
	return sft;
}


