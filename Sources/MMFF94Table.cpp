#include"MMFF94Table.h"

EnergyFundTable ReadEnergySolidParam(bool print_yes, int max_row_count)
{
	if (print_yes)PrintCmdSepTitle("基本参数表");
	vector<BSLine> bs_tb;
	vector<ABLine> ab_tb;
	vector<SBLine> sb_tb;
	vector<OPBLine> opb_tb;
	vector<TILine> ti_tb;
	vector<VDWLine> vdw_tb;
	vector<MSymLine> msym_tb;

	string bs_filepath = GetMMFFPath(BS_FILENAME);
	string ab_filepath = GetMMFFPath(AB_FILENAME);
	string sb_filepath = GetMMFFPath(SB_FILENAME);
	string opb_filepath = GetMMFFPath(OPB_FILENAME);
	string ti_filepath = GetMMFFPath(TI_FILENAME);
	string vdw_filepath = GetMMFFPath(VDW_FILENAME);
	string msys_filepath = GetMMFFPath(MSYM_FILENAME);

	ReadTableByTitle(bs_filepath, BS_TB_TITLE, bs_tb);
	ReadTableByTitle(ab_filepath, AB_TB_TITLE, ab_tb);
	ReadTableByTitle(sb_filepath, SB_TB_TITLE, sb_tb);
	ReadTableByTitle(opb_filepath, OPB_TB_TITLE, opb_tb);
	ReadTableByTitle(ti_filepath, TI_TB_TITLE, ti_tb);
	ReadTableByTitle(vdw_filepath, VDW_TB_TITLE, vdw_tb);
	ReadTableByTitle(msys_filepath, MSYM_TB_TITLE, msym_tb);


	MSYM_MAP msym_map = ChangeTableToMap(msym_tb);

	//
	if (print_yes)
	{
		//int max_row_count = 2;
		PrintTableTitle(BS_TB_SHORT_TITLE);
		PrintSpecialVector(bs_tb, "\t", max_row_count);

		PrintTableTitle(AB_TB_SHORT_TITLE);
		PrintSpecialVector(ab_tb, "\t", max_row_count);

		PrintTableTitle(SB_TB_SHORT_TITLE);
		PrintSpecialVector(sb_tb, "\t", max_row_count);

		PrintTableTitle(OPB_TB_SHORT_TITLE);
		PrintSpecialVector(opb_tb, "\t", max_row_count);

		PrintTableTitle(TI_TB_SHORT_TITLE);
		PrintSpecialVector(ti_tb, "\t", max_row_count);

		PrintTableTitle(VDW_TB_SHORT_TITLE);
		PrintSpecialVector(vdw_tb, "\t", max_row_count);

		PrintCmdSepTitle("MMFFType");
		PrintTableTitle(MSYM_TB_SHORT_TITLE);
		PrintSpecialVector(msym_tb, "\t", 10);
	}
	//if (print_yes) {
	//	PrintCmdSepTitle("杂化坐标表");
	//	all_sp_tb.Print();
	//}
	EnergyFundTable eft = { bs_tb, ab_tb, sb_tb, opb_tb, ti_tb, vdw_tb, msym_map };
	return eft;

}