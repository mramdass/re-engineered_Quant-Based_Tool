//	
//	This tool is meant to be used to generate evidence files for LAST
//	reQBT Evidence File Generator (C) 2015 Munieshwar Ramdass
//	
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with this program. If not, see <http://www.gnu.org/licenses/>.
//

//
//	reQBT Evidence File Generator
//	Munieshwar Ramdass	(C++ Programmer)
//	July 31, 2015
//

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <windows.h>
#include <stdio.h>
#include <winuser.h>
#include <conio.h>
#include <algorithm>

typedef std::string String;								// The following typedef used for CSV as a means of convinience
typedef std::vector<String> csv_row;
typedef std::vector<csv_row> csv_database;

using namespace std;

// Parameters: csv_database& db, string db_name
// Returns double (0.0 or 1.0)
double read_csv(csv_database& db, string db_name) {
	ifstream ifs(db_name);
	String csvLine;
	while (getline(ifs, csvLine)) {
		istringstream csvStream(csvLine);
		csv_row csvRow;
		String csvCol;
		while (getline(csvStream, csvCol, ','))
			csvRow.push_back(csvCol);
		db.push_back(csvRow);
	}
	ifs.close();
	return true;
}

// Parameters: const string& name
// Returns nothing
void create_files(const string& name) {
	csv_database db;
	read_csv(db, name);
	if (db.size() == 0)
		return;
	csv_row header(db[0]);
	header.push_back("Alleles");
	int file_num(1);
	while (true) {
		ofstream ofs("Evidence_" + to_string(file_num) + ".csv");
		csv_row row(db[file_num]);
		for (int i(0); i < header.size(); ++i) {
			if (!(header[i] == "REP" && (row[i] == "INC" || row[i] == "inc")))
				ofs << header[i] << ",";
		}
		ofs << endl;
		for (int i(0); i < row.size(); ++i) {
			if (header[i] == "REP" && (row[i] == "NEG" || row[i] == "neg"))
				ofs << "-1" << ",";
			else if (!(header[i] == "REP" && (row[i] == "INC" || row[i] == "inc")))
				ofs << row[i] << ",";
		}

		vector<double> alleles;
		for (int i(0); i < header.size(); ++i) {
			if (header[i] == "REP") {
				istringstream box(row[i]);
				string sub_box;
				while (getline(box, sub_box, ';'))
					alleles.push_back(atof(sub_box.c_str()));
			}
		}

		sort(alleles.begin(), alleles.end());
		alleles.erase(unique(alleles.begin(), alleles.end()), alleles.end());
		double neg(0.0);
		if (alleles[0] == neg)
			alleles.erase(alleles.begin(), alleles.begin() + 1);
		if (alleles.size() != 0) {
			for (int i(0); i < alleles.size(); ++i) {
				ofs << alleles[i];
				if (i != alleles.size() - 1)
					ofs << ";";
			}
		}
		else
			ofs << "-1";

		ofs << endl;
		ofs.close();
		++file_num;
		if (file_num == db.size())
			break;
	}
}

int main() {
	HWND stealth;
	AllocConsole();
	stealth = FindWindowA("ConsoleWindowClass", NULL);
	ShowWindow(stealth, 0);

	string name("case.csv");
	create_files(name);
}
