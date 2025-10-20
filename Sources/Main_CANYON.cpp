#include "types.h"
#include "Atom.h"
#include "Molecule.h"


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
	ReadTableByTitle(etb_filename, ATOM_TB_TITLE, atomtable);
	sort(atomtable.begin(), atomtable.end(), AtomTableCmp);

	PrintTableTitle(ATOM_TB_TITLE);
	PrintSpecialVector(atomtable);

	PrintCmdSepTitle("素数表读取");

	vector<PrimeNumber> primetable;
	ReadTable(ptb_filename, "PrimeNumber", primetable);


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
		//c.PrintBondTable();

		PrintCmdSepTitle("初始秩生成");
		//c.PrintAtomTable();
		//c.PrintOriRank();

		PrintCmdSepTitle("秩排序过程");
		c.MoleSortWithPN(primetable);
		cout << "\n排序结果：\n";
		//c.PrintRank();

		PrintCmdSepTitle("唯一SMILES生成");

		cout << "原SMILES：\n";
		cout << c.GetComSmiles() << endl << endl;
		cout << "唯一SMILES：\n";
		cout << c.GetCanSmiles() << endl;

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
				cout << "\nWarning: 文件夹创建失败！可能已存在。\n";
				char savechoice2;
				cout << "\n是否要继续覆盖节点信息表及邻接表？(Y/N)";
				cin >> savechoice2;
				if (!(savechoice2 == 'Y' || savechoice2 == 'y'))
				{
					continue;
				}
			}
			string now_time = GetCurrentTimeString();
			string atomtable_filename = folderpath + "/MNodeTable.csv";
			string bondtable_filename = folderpath + "/BondTable.csv";
			PrintCmdSepTitle("节点表储存");
			WriteTable<MNode>(atomtable_filename, MNODE_TB_TITLE, d.GetNodeTable(), true);
			cout << "节点表已储存至 " << atomtable_filename << endl;
			PrintCmdSepTitle("邻接表储存");
			WriteTable<NodeBonds>(bondtable_filename, ADJ_TB_TITLE, d.GetBondTable(), true);
			cout << "邻接表已储存至 " << bondtable_filename << endl;
		}
	}
	return 0;
}