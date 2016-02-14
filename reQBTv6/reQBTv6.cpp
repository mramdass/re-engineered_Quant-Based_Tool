//	
//	This tool can be used to achieve likelihood ratio for suspect(knowns_pn) and/or victim(knowns_pn) in a DNA mixture using probabilistic genotyping
//	re-engineered Quant-Based Tool (C) 2015 Munieshwar Ramdass, Nicholas Corpuz, Khagay Nagdimov
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
//	Programmer: Munieshwar Ramdass
//	August 2015
//


//
//	Application security has not been done - Will rely on an external program/person to check file format and exceptions
//	Simple procedures like checking for const methods or variables, or passing by references may not have been done
//	Code was written to work without regards to runtime or performance - Time costing transposing of matrix was done
//	Unsafe changing of global variables was done
//
//	Entire coding logic developed by Munieshwar Ramdass and can be explained in detail by him
//	Design of user interface done in R by Nick Corpuz and Khagay Nagdimov and can be explained by them
//	Any questions concerning the programs discussed above contact the team at: LAST.help.questions@gmail.com
//	

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <list>
#include <set>
#include <algorithm>
#include <iterator>
#include <functional>
#include <limits>
#include <cassert>
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <ppl.h>
#include <thread>
#include <future>
#include <direct.h>
#include <tuple>

typedef std::string String;								// The following typedef used for CSV as a means of convinience
typedef std::vector<String> csv_row;
typedef std::vector<csv_row> csv_database;

using namespace concurrency;
using namespace std;

int FILE_NUM = 1;
int RACE = 1;
double HOM_CONST = 0.0;		// inbreed - THETA

double PC0 = 0.0;			// Drop-in rates
double PC1 = 0.0;			// Not constant
double PC2 = 0.0;

//double PHET0 = 0.0;			// Drop-out rates
//double PHET1 = 0.0;			// Not constant
//double PHET2 = 0.0;
//double PHOM0 = 0.0;
//double PHOM1 = 0.0;

bool DEDUCIBLE = false;
double QUANT = 25.0;
int CONTRIBUTORS_PN = 1;
int CONTRIBUTORS_PD = 1;

double MINIMUM_WILD_FREQUENCY = 0.001;	// Not constant

int PRODUCT_GROUP = 15;				// Number of file LRs to multiple together in a group
									//char RUN_TYPE = 'E';				// Full, Moderate or Express run

string DATE_TIME = "Some Time";

class Genotypes;
class Person;
class Report;
class Allele;

class Thread_Constants {
public:
	Thread_Constants()
		: PHET0(0.0), PHET1(0.0), PHET2(0.0), PHOM0(0.0), PHOM1(0.0), RUN_TYPE('E'), RACE(1) {}

	double PHET0;			// Drop-out rates
	double PHET1;			// Not constant
	double PHET2;
	double PHOM0;
	double PHOM1;

	char RUN_TYPE;
	int RACE;

private:
	double PN_PHET0;
	double PN_PHET1;
	double PN_PHET2;
	double PN_PHOM0;
	double PN_PHOM1;

	double PD_PHET0;
	double PD_PHET1;
	double PD_PHET2;
	double PD_PHOM0;
	double PD_PHOM1;

	friend Report;
};

class Timer {
public:
	Timer() : start(clock()) {}
	double elapsed() { return (clock() - start) / CLOCKS_PER_SEC; }
	void reset() { start = clock(); }
private:
	double start;
};

class Allele {
public:
	Allele()
		: locus("W"), length(-1), freq(1.0) {}

	Allele(string locus, double length, double freq)
		: locus(locus), length(length), freq(freq) {}

	bool operator==(const Allele& right) const { return (locus == right.locus && length == right.length && freq == right.freq); }

	bool operator!=(const Allele& right) const { return !(*this == right); }

private:
	string locus;
	double length;
	double freq;

	friend Genotypes;
	friend Person;
	friend Report;

	friend double generate_wild_allele_freq(vector<Allele>& alleles);
	friend bool get_data(csv_database& evidence_db, csv_database& allele_db, vector<vector<Allele>>& replicates, vector<Allele>& distinct_alleles, vector<Person>& knowns_pn, vector<Person>& knowns_pd, int& unknowns_pn, int& unknowns_pd, double& quant, string& case_name, string& locus, int& contributors_pn, int& contributors_pd, int& contributors, Thread_Constants& tc);
};

class Genotypes {
public:
	Genotypes(vector<Allele>& distinct_alleles)
		: alleles(distinct_alleles) {

		//if (distinct_alleles.size() == 0)
		//alleles.push_back(Allele("W", -1, 1.0));
		//else
		//for (int i(0); i < distinct_alleles.size(); ++i)
		//alleles.push_back(distinct_alleles[i]);

		generate_allele_pairs();
		generate_freqs();
	}

	void generate_freqs() {
		for (int i(0); i < alleles.size(); ++i) {
			for (int j(i); j < alleles.size(); ++j) {
				if (i == j)
					freqs.push_back(calc_hom(alleles[i]));
				else
					freqs.push_back(calc_het(alleles[i], alleles[j]));
			}
		}
	}

	void generate_allele_pairs() {
		for (int i(0); i < alleles.size(); ++i) {
			for (int j(i); j < alleles.size(); ++j)
				allele_comb.push_back(pair<Allele, Allele>(alleles[i], alleles[j]));
		}
	}

	double calc_hom(Allele a) const { return a.freq * a.freq + a.freq * HOM_CONST * (1 - a.freq); }

	double calc_het(Allele a, Allele b) const { return 2 * a.freq * b.freq; }


private:
	vector<Allele> alleles;							// Distinct alleles must be in here
	vector<pair<Allele, Allele>> allele_comb;		// Possible genotypes from distinct alleles are here in order
	vector<double> freqs;							// Associated frequencies from respective 'allele_comb' are in here

	friend Allele;
	friend Person;
	friend Report;
};

class Person {
public:
	Person()
		: a(Allele()), b(Allele()) {}

	Person(Allele a, Allele b)
		: a(a), b(b) {

		generate_freq();

		if (a == b) {
			hom = true;
			het = false;
		}
		else {
			hom = false;
			het = true;
		}
	}

	void generate_freq() {
		if (a == b)
			freq = calc_hom(a);
		else
			freq = calc_het(a, b);
	}

	double calc_hom(Allele a) const { return a.freq * a.freq + a.freq * HOM_CONST * (1 - a.freq); }

	double calc_het(Allele a, Allele b) const { return 2 * a.freq * b.freq; }

	bool operator==(const Person& right) const { return (a == right.a && b == right.b && freq == right.freq && hom == right.hom && het == right.het); }

	bool operator!=(const Person& right) const { return !(*this == right); }

	bool operator<(const Person& right) const {
		if (a == right.a) return b.length < right.b.length;
		else return a.length < right.a.length;
	}

private:
	Allele a;
	Allele b;
	double freq;
	bool hom;
	bool het;

	friend Allele;
	friend Genotypes;
	friend Report;
};

class Report {
public:
	Report(csv_database& drop_out_db, vector<Person>& knowns_pn, int unknowns_pn, vector<Person>& knowns_pd, int unknowns_pd, Genotypes& genotypes, vector<vector<Allele>>& replicates, const string& case_name, const string& locus, Thread_Constants& thread_constants)
		: knowns_pn(knowns_pn), knowns_pd(knowns_pd), unknowns_pn(unknowns_pn), unknowns_pd(unknowns_pd), genotypes(genotypes), replicates(replicates), population_pn(pow(genotypes.allele_comb.size(), (knowns_pn.size() + unknowns_pn))), population_pd(pow(genotypes.allele_comb.size(), (knowns_pd.size() + unknowns_pd))), locus(locus), case_name(case_name), tc(thread_constants) {

		Timer t;
		t.reset();

		ofstream ofs, efs;
		ofs.open(DATE_TIME + "/output.csv", ios::app);

		for (int i(0); i < case_name.size(); ++i)
			if (this->case_name[i] == '/')
				this->case_name[i] = '-';

		if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') {
			efs.open(DATE_TIME + "/Evidence_" + to_string(FILE_NUM) + "_" + to_string(tc.RACE) + "_" + locus + "_" + this->case_name + ".csv", ios::app);

			efs << "Case Name:," << case_name << "\nLocus:," << locus << "\nQuant:," << QUANT << "\nPn size:," << (knowns_pn.size() + unknowns_pn) << "\nPd size:," << (knowns_pd.size() + unknowns_pd) << endl;
			if (DEDUCIBLE) efs << "Dedicible?,YES" << endl;
			else efs << "Dedicible?,NO" << endl;
			efs << "Race:,";
			if (tc.RACE == 1)
				efs << "BLACK" << endl;
			else if (tc.RACE == 2)
				efs << "CAUCASIAN" << endl;
			else if (tc.RACE == 3)
				efs << "HISPANIC" << endl;
			else if (tc.RACE == 4)
				efs << "ASIAN" << endl;
			efs << "Known(s) Pn:" << endl;
			if (knowns_pn.size() == 0) efs << ",No Knowns Pn" << endl;
			efs << ",";
			for (int i(0); i < knowns_pn.size(); ++i) {
				efs << "Person " << (i + 1) << ",(" << knowns_pn[i].a.length << "; " << knowns_pn[i].b.length << ")" << endl;
			}
			efs << "Known(s) Pn:" << endl;
			if (knowns_pd.size() == 0) efs << ",No Knowns Pd" << endl;
			efs << ",";
			for (int i(0); i < knowns_pd.size(); ++i) {
				efs << "Person " << (i + 1) << ",(" << knowns_pd[i].a.length << "; " << knowns_pd[i].b.length << ")" << endl;
			}
			efs << "\nReplicate(s):" << endl;
			for (int i(0); i < replicates.size(); ++i) {
				efs << ",Replicate " << (i + 1) << ",";
				for (int j(0); j < replicates[i].size(); ++j) {
					efs << replicates[i][j].length << ";";
				}
				efs << endl;
			}
			efs << endl;

		}

		vector<int> tmp;	// Handling absence of distinct alleles
		for (int i(0); i < genotypes.allele_comb.size(); ++i)
			tmp.push_back(i);
		if (unknowns_pn + knowns_pn.size() > 0 && unknowns_pn + knowns_pn.size() <= tmp.size())
			permute_vector_driver(tmp, unknowns_pn + knowns_pn.size(), 'n');
		else
			indicies_pn.push_back(tmp);
		if (unknowns_pd + knowns_pd.size() > 0 && unknowns_pd + knowns_pd.size() <= tmp.size())
			permute_vector_driver(tmp, unknowns_pd + knowns_pd.size(), 'd');
		else
			indicies_pd.push_back(tmp);

		allocate_population_space();

		generate_persons();

		generate_populations();

		vector<Person> wild_vector;
		if (persons[0].a.length == -1 && persons[0].b.length == -1) {
			population_pn.clear();
			for (int i(0); i < knowns_pn.size() + unknowns_pn; ++i)
				wild_vector.push_back(persons[0]);
			population_pn.push_back(wild_vector);
			wild_vector.clear();
		}
		if (persons[0].a.length == -1 && persons[0].b.length == -1) {
			population_pd.clear();
			for (int i(0); i < knowns_pd.size() + unknowns_pd; ++i)
				wild_vector.push_back(persons[0]);
			population_pd.push_back(wild_vector);
			wild_vector.clear();
		}

		check_person();

		if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') {
			efs << "\nGenotype Freq:" << endl;
			for (int i(0); i < persons.size(); ++i)
				efs << ",(" << persons[i].a.length << "; " << persons[i].b.length << "),(" << persons[i].a.freq << "; " << persons[i].b.freq << ")," << persons[i].freq << endl;
		}

		bool pn_done(false), pd_done(false);
		if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M') {
			pn_done = generate_pn(drop_out_db);
			pd_done = generate_pd(drop_out_db);
		}
		else if (tc.RUN_TYPE == 'C') {
			pn_done = generate_pn_express_with_data(drop_out_db);
			pd_done = generate_pd_express_with_data(drop_out_db);
		}
		else {
			pn_done = generate_pn_express(drop_out_db);
			pd_done = generate_pd_express(drop_out_db);
		}

		lr = pn / pd;

		if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'C') generate_hypothetical_scenario(drop_out_db);

		if (tc.RACE == 1) {
			cout << case_name << " - " << locus << " BLACK - LR: " << lr << " - Pn: " << pn << " - Pd: " << pd << endl;
			ofs << case_name << "," << locus << ",BLACK,," << lr << "," << pn << "," << pd << ",," << t.elapsed()
				<< ",," << tc.PN_PHET0 << "," << tc.PN_PHET1 << "," << tc.PN_PHET2 << "," << tc.PN_PHOM0 << "," << tc.PN_PHOM1
				<< ",," << tc.PD_PHET0 << "," << tc.PD_PHET1 << "," << tc.PD_PHET2 << "," << tc.PD_PHOM0 << "," << tc.PD_PHOM1 << endl;
			if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') efs << "\nBLACK\n";
		}
		else if (tc.RACE == 2) {
			cout << case_name << " - " << locus << " CAUCASIAN - LR: " << lr << " - Pn: " << pn << " - Pd: " << pd << endl;
			ofs << case_name << "," << locus << ",CAUCASIAN,," << lr << "," << pn << "," << pd << ",," << t.elapsed()
				<< ",," << tc.PN_PHET0 << "," << tc.PN_PHET1 << "," << tc.PN_PHET2 << "," << tc.PN_PHOM0 << "," << tc.PN_PHOM1
				<< ",," << tc.PD_PHET0 << "," << tc.PD_PHET1 << "," << tc.PD_PHET2 << "," << tc.PD_PHOM0 << "," << tc.PD_PHOM1 << endl;
			if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') efs << "\nCAUCASIAN\n";
		}
		else if (tc.RACE == 3) {
			cout << case_name << " - " << locus << " HISPANIC - LR: " << lr << " - Pn: " << pn << " - Pd: " << pd << endl;
			ofs << case_name << "," << locus << ",HISPANIC,," << lr << "," << pn << "," << pd << ",," << t.elapsed()
				<< ",," << tc.PN_PHET0 << "," << tc.PN_PHET1 << "," << tc.PN_PHET2 << "," << tc.PN_PHOM0 << "," << tc.PN_PHOM1
				<< ",," << tc.PD_PHET0 << "," << tc.PD_PHET1 << "," << tc.PD_PHET2 << "," << tc.PD_PHOM0 << "," << tc.PD_PHOM1 << endl;
			if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') efs << "\nHISPANIC\n";
		}
		else if (tc.RACE == 4) {
			cout << case_name << " - " << locus << " ASIAN - LR: " << lr << " - Pn: " << pn << " - Pd: " << pd << endl;
			ofs << case_name << "," << locus << ",ASIAN,," << lr << "," << pn << "," << pd << ",," << t.elapsed()
				<< ",," << tc.PN_PHET0 << "," << tc.PN_PHET1 << "," << tc.PN_PHET2 << "," << tc.PN_PHOM0 << "," << tc.PN_PHOM1
				<< ",," << tc.PD_PHET0 << "," << tc.PD_PHET1 << "," << tc.PD_PHET2 << "," << tc.PD_PHOM0 << "," << tc.PD_PHOM1 << endl;
			if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') efs << "\nASIAN\n";
		}

		if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'M' || tc.RUN_TYPE == 'C') {

			efs << "\n,,LR:," << lr << ",,Pn:," << pn << ",,Pd:," << pd << endl << endl;

			if (tc.RUN_TYPE == 'F' || tc.RUN_TYPE == 'C') {
				for (int i(0); i < denominators.size(); ++i) {
					for (int j(0); j < numerators.size(); ++j) {
						efs << ",LR:," << numerators[j].first / denominators[i].first << ",,Pn:," << numerators[j].first << ",,Pd:," << denominators[i].first << ",,";
						for (int k(0); k < numerators[j].second.size(); ++k)
							efs << "(" << numerators[j].second[k].a.length << "; " << numerators[j].second[k].b.length << "),,";
						efs << unknowns_pn << ",Uknown(s) Pn,,";
						for (int k(0); k < denominators[i].second.size(); ++k)
							efs << "(" << denominators[i].second[k].a.length << "; " << denominators[i].second[k].b.length << "),,";
						efs << unknowns_pd << ",Uknown(s) Pd" << endl;
					}
					efs << endl;
				}
				efs << endl;
			}
			efs << "\n,Time Elapsed:," << t.elapsed() << ",Second(s)" << endl;
		}

		efs.close();
		ofs.close();
	}

	bool is_subset_pn(int ind) {
		if (equal(population_pn[ind].begin(), population_pn[ind].begin() + knowns_pn.size(), knowns_pn.begin()))
			return true;
		return false;
	}

	bool is_subset_pd(int ind) {
		if (equal(population_pd[ind].begin(), population_pd[ind].begin() + knowns_pd.size(), knowns_pd.begin()))
			return true;
		return false;
	}

	bool generate_pd(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PD, 'd');

		double product(1.0), sum(0.0);

		ofstream ofs;
		ofs.open(DATE_TIME + "/Evidence_" + to_string(FILE_NUM) + "_" + to_string(tc.RACE) + "_" + locus + "_" + case_name + ".csv", ios::app);

		ofs << "\nDenominator:" << endl;
		ofs << "\n,HET0:," << tc.PHET0 << ",,pC0:," << PC0 << "\n,HET1:," << tc.PHET1 << ",,pC1:," << PC1
			<< "\n,HET2:," << tc.PHET2 << ",,pC2:," << PC2 << "\n,HOM0:," << tc.PHOM0 << "\n,HOM1:," << tc.PHOM1
			<< ",,Theta:," << HOM_CONST << endl << endl;


		for (int i(0); i < population_pd[0].size(); ++i)
			if (i >= (unknowns_pd + knowns_pd.size()) - ((unknowns_pd + knowns_pd.size()) - knowns_pd.size()))
				ofs << ",Person " << (i + 1) << ",Genotype Freq,";
		ofs << ",";
		for (int i(0); i < population_pd[0].size(); ++i) {
			for (int j(0); j < replicates.size(); ++j)
				ofs << "Drop-Out,";
			ofs << ",";
		}
		for (int i(0); i < replicates.size(); ++i)
			ofs << "Drop-In,";
		ofs << ",Row Product" << endl;

		for (int i(0); i < population_pd.size(); ++i) {
			if (is_subset_pd(i)) {
				for (int j(0); j < population_pd[i].size(); ++j) {	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pd + knowns_pd.size()) - ((unknowns_pd + knowns_pd.size()) - knowns_pd.size())) {
						product *= population_pd[i][j].freq;
						ofs << ",(" << population_pd[i][j].a.length << "; " << population_pd[i][j].b.length << ")," << population_pd[i][j].freq << ",";
					}

					//ofs << population_pd[i][j].freq << ",";
				}
				ofs << ",";
				for (int j(0); j < population_pd[i].size(); ++j) {	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k) {
						double prod = drop_out(population_pd[i][j], replicates[k]);
						product *= prod;

						ofs << prod << ",";
					}
					ofs << ",";
				}
				for (int j(0); j < replicates.size(); ++j) {		// MULTIPLY DROP-IN
					double prod = drop_in(population_pd[i], replicates[j]);
					product *= prod;

					ofs << prod << ",";
				}
				sum += product;
				ofs << "," << product << endl;
				product = 1.0;
			}
		}
		pd = sum;

		ofs.close();
		return true;
	}

	bool generate_pn(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PN, 'n');

		double product(1.0), sum(0.0), target_pn(0.0);
		double other_product(1.0), other_sum(0.0);

		ofstream ofs;
		ofs.open(DATE_TIME + "/Evidence_" + to_string(FILE_NUM) + "_" + to_string(tc.RACE) + "_" + locus + "_" + case_name + ".csv", ios::app);

		ofs << "\nNumerator:" << endl;
		ofs << "\n,HET0:," << tc.PHET0 << ",,pC0:," << PC0 << "\n,HET1:," << tc.PHET1 << ",,pC1:," << PC1
			<< "\n,HET2:," << tc.PHET2 << ",,pC2:," << PC2 << "\n,HOM0:," << tc.PHOM0 << "\n,HOM1:," << tc.PHOM1
			<< ",,Theta:," << HOM_CONST << endl << endl;

		for (int i(0); i < population_pn[0].size(); ++i)
			if (i >= (unknowns_pn + knowns_pn.size()) - ((unknowns_pn + knowns_pn.size()) - knowns_pn.size()))
				ofs << ",Person " << (i + 1) << ",Genotype Freq,";
		ofs << ",";
		for (int i(0); i < population_pn[0].size(); ++i) {
			for (int j(0); j < replicates.size(); ++j)
				ofs << "Drop-Out P" << (i + 1) << " R" << (j + 1) << ",";
			ofs << ",";
		}
		for (int i(0); i < replicates.size(); ++i)
			ofs << "Drop-In R" << (i + 1) << ",";
		ofs << ",Row Product" << endl;

		for (int i(0); i < population_pn.size(); ++i) {
			if (is_subset_pn(i)) {
				for (int j(0); j < population_pn[i].size(); ++j) {	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pn + knowns_pn.size()) - ((unknowns_pn + knowns_pn.size()) - knowns_pn.size())) {
						product *= population_pn[i][j].freq;
						ofs << ",(" << population_pn[i][j].a.length << "; " << population_pn[i][j].b.length << ")," << population_pn[i][j].freq << ",";
					}

					//ofs << population_pn[i][j].freq << ",";
				}
				ofs << ",";
				for (int j(0); j < population_pn[i].size(); ++j) {	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k) {
						double prod = drop_out(population_pn[i][j], replicates[k]);
						product *= prod;

						ofs << prod << ",";
					}
					ofs << ",";
				}
				for (int j(0); j < replicates.size(); ++j) {		// MULTIPLY DROP-IN
					double prod = drop_in(population_pn[i], replicates[j]);
					product *= prod;

					ofs << prod << ",";
				}
				sum += product;
				ofs << "," << product << endl;
				product = 1.0;
			}
			//numerators.push_back(sum);
		}
		pn = sum;
		ofs.close();
		return true;
	}

	bool generate_pd_express_with_data(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PD, 'd');

		double product(1.0), sum(0.0);

		ofstream ofs;
		ofs.open(DATE_TIME + "/Evidence_" + to_string(FILE_NUM) + "_" + to_string(tc.RACE) + "_" + locus + "_" + case_name + ".csv", ios::app);
		ofs << "\nCase:," << case_name << "," << locus << endl;
		ofs << "\nDenominator:" << endl;
		ofs << "\n,HET0:," << tc.PHET0 << ",,pC0:," << PC0 << "\n,HET1:," << tc.PHET1 << ",,pC1:," << PC1
			<< "\n,HET2:," << tc.PHET2 << ",,pC2:," << PC2 << "\n,HOM0:," << tc.PHOM0 << "\n,HOM1:," << tc.PHOM1
			<< ",,Theta:," << HOM_CONST << endl << endl;

		for (int i(0); i < population_pd.size(); ++i) {
			if (is_subset_pd(i)) {
				for (int j(0); j < population_pd[i].size(); ++j) {	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pd + knowns_pd.size()) - ((unknowns_pd + knowns_pd.size()) - knowns_pd.size())) {
						product *= population_pd[i][j].freq;
					}
				}
				for (int j(0); j < population_pd[i].size(); ++j) {	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k) {
						double prod = drop_out(population_pd[i][j], replicates[k]);
						product *= prod;
					}
				}
				for (int j(0); j < replicates.size(); ++j) {		// MULTIPLY DROP-IN
					double prod = drop_in(population_pd[i], replicates[j]);
					product *= prod;
				}
				sum += product;
				product = 1.0;
			}
		}
		pd = sum;

		ofs.close();
		return true;
	}

	bool generate_pn_express_with_data(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PN, 'n');

		double product(1.0), sum(0.0), target_pn(0.0);
		double other_product(1.0), other_sum(0.0);

		ofstream ofs;
		ofs.open(DATE_TIME + "/Evidence_" + to_string(FILE_NUM) + "_" + to_string(tc.RACE) + "_" + locus + "_" + case_name + ".csv", ios::app);
		ofs << "\nCase:," << case_name << "," << locus << endl;
		ofs << "\nNumerator:" << endl;
		ofs << "\n,HET0:," << tc.PHET0 << ",,pC0:," << PC0 << "\n,HET1:," << tc.PHET1 << ",,pC1:," << PC1
			<< "\n,HET2:," << tc.PHET2 << ",,pC2:," << PC2 << "\n,HOM0:," << tc.PHOM0 << "\n,HOM1:," << tc.PHOM1
			<< ",,Theta:," << HOM_CONST << endl << endl;

		for (int i(0); i < population_pn.size(); ++i) {
			if (is_subset_pn(i)) {
				for (int j(0); j < population_pn[i].size(); ++j) {	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pn + knowns_pn.size()) - ((unknowns_pn + knowns_pn.size()) - knowns_pn.size())) {
						product *= population_pn[i][j].freq;
					}
				}
				for (int j(0); j < population_pn[i].size(); ++j) {	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k) {
						double prod = drop_out(population_pn[i][j], replicates[k]);
						product *= prod;
					}
				}
				for (int j(0); j < replicates.size(); ++j) {		// MULTIPLY DROP-IN
					double prod = drop_in(population_pn[i], replicates[j]);
					product *= prod;
				}
				sum += product;
				product = 1.0;
			}
			//numerators.push_back(sum);
		}
		pn = sum;
		ofs.close();
		return true;
	}

	bool generate_pd_express(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PD, 'd');

		double product(1.0), sum(0.0);

		//Concurrency::parallel_for(0, (int)population_pd.size(), [&] (int i) {
		for (int i(0); i < population_pd.size(); ++i) {
			if (is_subset_pd(i)) {
				for (int j(0); j < population_pd[i].size(); ++j)	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pd + knowns_pd.size()) - ((unknowns_pd + knowns_pd.size()) - knowns_pd.size()))
						product *= population_pd[i][j].freq;

				for (int j(0); j < population_pd[i].size(); ++j)	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k)
						product *= drop_out(population_pd[i][j], replicates[k]);

				for (int j(0); j < replicates.size(); ++j)		// MULTIPLY DROP-IN
					product *= drop_in(population_pd[i], replicates[j]);

				sum += product;
				product = 1.0;
			}
		}
		pd = sum;
		return true;
	}

	bool generate_pn_express(csv_database& drop_out_db) {
		set_drop_out(drop_out_db, CONTRIBUTORS_PN, 'n');

		double product(1.0), sum(0.0);

		//Concurrency::parallel_for(0, (int)population_pn.size(), [&](int i) {
		for (int i(0); i < population_pn.size(); ++i) {
			if (is_subset_pn(i)) {
				for (int j(0); j < population_pn[i].size(); ++j)	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pn + knowns_pn.size()) - ((unknowns_pn + knowns_pn.size()) - knowns_pn.size()))
						product *= population_pn[i][j].freq;

				for (int j(0); j < population_pn[i].size(); ++j)	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k)
						product *= drop_out(population_pn[i][j], replicates[k]);

				for (int j(0); j < replicates.size(); ++j)		// MULTIPLY DROP-IN
					product *= drop_in(population_pn[i], replicates[j]);

				sum += product;
				product = 1.0;
			}
		}
		pn = sum;
		return true;
	}

	double drop_out(Person p, vector<Allele>& alleles) {		// Determining drop-out rate for a specific genotype (Person p)
		if (p.hom) {
			for (int i(0); i < alleles.size(); ++i) {
				if (p.a == alleles[i] && p.b == alleles[i])
					return tc.PHOM0;
			}
			return tc.PHOM1;
		}
		else if (p.het) {
			bool present_first = false;
			bool present_second = false;
			for (int i(0); i < alleles.size(); ++i) {
				if (p.a == alleles[i])
					present_first = true;
				if (p.b == alleles[i])
					present_second = true;
			}
			if (present_first && present_second)
				return tc.PHET0;
			else if (!present_first && !present_second)
				return tc.PHET2;
			else if ((present_first && !present_second) || (!present_first && present_second))
				return tc.PHET1;
		}
	}

	double drop_in(vector<Person>& p, vector<Allele>& alleles) {		// Determining drop-in rates for a specific replicate
		vector<Allele> tmp;
		for (int i(0); i < p.size(); ++i) {
			tmp.push_back(p[i].a);
			tmp.push_back(p[i].b);
		}
		tmp.erase(unique(begin(tmp), end(tmp)), end(tmp));
		int c(0);

		for (int i(0); i < alleles.size(); ++i)
			if (find(tmp.begin(), tmp.end(), alleles[i]) == tmp.end())
				++c;

		if (c == 0)
			return PC0;
		else if (c == 1)
			return PC1;
		else
			return PC2;
	}

	void check_person() {		// Checks if knowns' alleles in the numerator have been observed, and acts accordingly
		for (int ii(0); ii < knowns_pn.size(); ++ii) {
			bool found_two = false;
			bool found_a = false;
			bool found_b = false;
			for (int i(0); i < replicates.size(); ++i) {
				for (int j(0); j < replicates[i].size(); ++j) {
					if (knowns_pn[ii].a == replicates[i][j]) {
						found_a = true;
					}
					if (knowns_pn[ii].b == replicates[i][j]) {
						found_b = true;
					}
				}
			}
			if (found_a && !found_b) {		// Change second allele to wild if it is not observed
				knowns_pn[ii].b = persons[persons.size() - 1].b;
				knowns_pn[ii].generate_freq();
			}
			if (!found_a && found_b) {		// Change first allele to wild if it is not found
				knowns_pn[ii].a = knowns_pn[ii].b;
				knowns_pn[ii].b = persons[persons.size() - 1].b;
				knowns_pn[ii].generate_freq();
			}
			if (!found_a && !found_b) {		// Chane both alleles to wild if it is not found
				knowns_pn[ii].b = persons[persons.size() - 1].b;
				knowns_pn[ii].a = persons[persons.size() - 1].a;
				knowns_pn[ii].generate_freq();
			}
			if (knowns_pn[ii].a == knowns_pn[ii].b) {		// Book keeping
				knowns_pn[ii].hom = true;
				knowns_pn[ii].het = false;
			}
			else {
				knowns_pn[ii].het = true;
				knowns_pn[ii].hom = false;
			}
		}
		for (int ii(0); ii < knowns_pd.size(); ++ii) {	// For Knowns in the denominator
			bool found_two = false;
			bool found_a = false;
			bool found_b = false;
			for (int i(0); i < replicates.size(); ++i) {
				for (int j(0); j < replicates[i].size(); ++j) {
					if (knowns_pd[ii].a == replicates[i][j]) {
						found_a = true;
					}
					if (knowns_pd[ii].b == replicates[i][j]) {
						found_b = true;
					}
				}
			}
			if (found_a && !found_b) {		// Change second allele to wild if it is not observed
				knowns_pd[ii].b = persons[persons.size() - 1].b;
				knowns_pd[ii].generate_freq();
			}
			if (!found_a && found_b) {		// Change first allele to wild if it is not found
				knowns_pd[ii].a = knowns_pd[ii].b;
				knowns_pd[ii].b = persons[persons.size() - 1].b;
				knowns_pd[ii].generate_freq();
			}
			if (!found_a && !found_b) {		// Chane both alleles to wild if it is not found
				knowns_pd[ii].b = persons[persons.size() - 1].b;
				knowns_pd[ii].a = persons[persons.size() - 1].a;
				knowns_pd[ii].generate_freq();
			}
			if (knowns_pd[ii].a == knowns_pd[ii].b) {		// Book keeping
				knowns_pd[ii].hom = true;
				knowns_pd[ii].het = false;
			}
			else {
				knowns_pd[ii].het = true;
				knowns_pd[ii].hom = false;
			}
		}
	}

	void generate_populations() {		// Generates all combinaions of genotypes separately for numerator and denominator
		for (int i(0); i < indicies_pn.size(); ++i) {
			for (int j(0); j < indicies_pn[i].size(); ++j) {
				Person temp(genotypes.allele_comb[indicies_pn[i][j]].first, genotypes.allele_comb[indicies_pn[i][j]].second);
				population_pn[i][j] = temp;
			}
		}
		for (int i(0); i < indicies_pd.size(); ++i) {
			for (int j(0); j < indicies_pd[i].size(); ++j) {
				Person temp(genotypes.allele_comb[indicies_pd[i][j]].first, genotypes.allele_comb[indicies_pd[i][j]].second);
				population_pd[i][j] = temp;
			}
		}
	}

	void generate_persons() {		// Generates all possible genotypes of people that can be found based on observed alleles
		for (int i(0); i < genotypes.alleles.size(); ++i)
			for (int j(i); j < genotypes.alleles.size(); ++j)
				persons.push_back(Person(genotypes.alleles[i], genotypes.alleles[j]));
	}

	void permute_vector(const vector<int>& number_items, unsigned int length, vector<unsigned int>& position, unsigned int depth, unsigned int perimeter, char population_ID) {
		if (depth >= length) {
			vector<int> index_vector;
			for (unsigned int i = 0; i < position.size(); ++i)
				index_vector.push_back(number_items[position[i]]);
			if (population_ID == 'n')
				indicies_pn.push_back(index_vector);
			else if (population_ID == 'd')
				indicies_pd.push_back(index_vector);
			return;
		}

		for (unsigned int i = 0; i < number_items.size(); ++i) { // CHANGE i TO perimeter TO MAKE THIS A COMBINATION WITH REPETITION
			position[depth] = i;
			permute_vector(number_items, length, position, depth + 1, i, population_ID);
		}
		return;
	}

	void permute_vector_driver(const vector<int>& number_items, unsigned int length, char population_ID) {
		assert(length > 0 && length <= number_items.size());
		vector<unsigned int> positions(length, 0);
		permute_vector(number_items, length, positions, 0, 0, population_ID);
	}

	void allocate_population_space() {
		for (int i(0); i < population_pn.size(); ++i) {
			vector<Person> tmp_persons((knowns_pn.size() + unknowns_pn));
			population_pn[i] = tmp_persons;
		}
		for (int i(0); i < population_pd.size(); ++i) {
			vector<Person> tmp_persons((knowns_pd.size() + unknowns_pd));
			population_pd[i] = tmp_persons;
		}
	}

	void generate_hypothetical_scenario(csv_database& drop_out_db) {
		for (int i(0); i < population_pn.size(); ++i) {
			vector<Person> temp;
			if (i > 0 && !equal(population_pn[i - 1].begin(), population_pn[i - 1].end() - unknowns_pn, population_pn[i].begin(), population_pn[i].end() - unknowns_pn)) {
				for (int j(0); j < population_pn[i].size() - unknowns_pn; ++j)
					temp.push_back(population_pn[i][j]);
				pn_combinations.push_back(temp);
			}
			else if (i == 0) {
				for (int j(0); j < population_pn[i].size() - unknowns_pn; ++j)
					temp.push_back(population_pn[i][j]);
				pn_combinations.push_back(temp);
			}
		}

		for (int i(0); i < population_pd.size(); ++i) {
			vector<Person> temp;
			if (i > 0 && !equal(population_pd[i - 1].begin(), population_pd[i - 1].end() - unknowns_pd, population_pd[i].begin(), population_pd[i].end() - unknowns_pd)) {
				for (int j(0); j < population_pd[i].size() - unknowns_pd; ++j)
					temp.push_back(population_pd[i][j]);
				pd_combinations.push_back(temp);
			}
			else if (i == 0) {
				for (int j(0); j < population_pd[i].size() - unknowns_pd; ++j)
					temp.push_back(population_pd[i][j]);
				pd_combinations.push_back(temp);
			}
		}

		cout << "\n\tGenerating Scenarios with " << pn_combinations.size() << " different numerator(s) and\n\t"
			<< pd_combinations.size() << " different denominator(s).\n\tThis will take some time." << endl;
		numerators = vector<pair<double, vector<Person>>>(pn_combinations.size());
		denominators = vector<pair<double, vector<Person>>>(pd_combinations.size());
		/*
		parallel_for(size_t(0), pn_combinations.size(), [&](size_t i) {
		generate_px(drop_out_db, pn_combinations[i], 'n', i);
		});
		*/
		for (int i(0); i < pn_combinations.size(); ++i)
			generate_px(drop_out_db, pn_combinations[i], 'n', i);

		/*
		parallel_for(size_t(0), pd_combinations.size(), [&](size_t i) {
		generate_px(drop_out_db, pd_combinations[i], 'd', i);
		});
		*/
		for (int i(0); i < pd_combinations.size(); ++i)
			generate_px(drop_out_db, pd_combinations[i], 'd', i);

		cout << "\tEnding..." << endl;
	}

	bool generate_px(csv_database& drop_out_db, vector<Person>& knowns, char known_ID, int index) {		// The generic function that chould have been used if last few lines are modified
		if (known_ID == 'n') swap_drop_out_rates('n');
		else if (known_ID == 'd') swap_drop_out_rates('d');

		double product(1.0), sum(0.0);

		for (int i(0); i < population_pn.size(); ++i) {
			if (is_subset_px(i, knowns, known_ID)) {
				for (int j(0); j < population_pn[i].size(); ++j)	// MULTIPLY PERSONS FREQ
					if (j >= (unknowns_pn + knowns.size()) - ((unknowns_pn + knowns.size()) - knowns.size()))
						product *= population_pn[i][j].freq;

				for (int j(0); j < population_pn[i].size(); ++j)	// MULTIPLY DROP-OUT
					for (int k(0); k < replicates.size(); ++k)
						product *= drop_out(population_pn[i][j], replicates[k]);

				for (int j(0); j < replicates.size(); ++j)		// MULTIPLY DROP-IN
					product *= drop_in(population_pn[i], replicates[j]);

				sum += product;
				product = 1.0;
			}
		}
		//if (known_ID == 'n') numerators.push_back(pair<double, vector<Person>>(sum, knowns));
		//else if (known_ID == 'd') denominators.push_back(pair<double, vector<Person>>(sum, knowns));
		if (known_ID == 'n') numerators[index] = pair<double, vector<Person>>(sum, knowns);
		else if (known_ID == 'd') denominators[index] = pair<double, vector<Person>>(sum, knowns);

		return true;
	}

	bool is_subset_px(int ind, vector<Person>& knowns, char known_ID) {		// The generic function that should have been used
		if (equal(population_pn[ind].begin(), population_pn[ind].begin() + knowns.size(), knowns.begin()) && known_ID == 'n')
			return true;
		else if (equal(population_pd[ind].begin(), population_pd[ind].begin() + knowns.size(), knowns.begin()) && known_ID == 'd')
			return true;
		return false;
	}

	void set_drop_out(csv_database& drop_out_db, int contributors, char term) {
		int quant(QUANT);

		string dnd("ND");
		if (DEDUCIBLE) dnd = "D";
		else dnd = "ND";

		double full(1.0);
		for (int i(0); i < drop_out_db.size(); ++i) {					// Getting proper drop-out rates based on quant and number of contributors
			for (int j(0); j < drop_out_db[i].size(); ++j) {
				if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 25 && quant < 50) {
					double b(0.0), l(25.0), h(50.0);
					double slope((atof(drop_out_db[i][j + 3].c_str()) - atof(drop_out_db[i][j + 2].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 2].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 50 && quant < 100) {
					double b(0.0), l(50.0), h(100.0);
					double slope((atof(drop_out_db[i][j + 4].c_str()) - atof(drop_out_db[i][j + 3].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 3].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 100 && quant < 150) {
					double b(0.0), l(100.0), h(150.0);
					double slope((atof(drop_out_db[i][j + 5].c_str()) - atof(drop_out_db[i][j + 4].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 4].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 150 && quant < 250) {
					double b(0.0), l(150.0), h(250.0);
					double slope((atof(drop_out_db[i][j + 6].c_str()) - atof(drop_out_db[i][j + 5].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 5].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 250 && quant <= 500) {
					double b(0.0), l(250.0), h(500.0);
					double slope((atof(drop_out_db[i][j + 7].c_str()) - atof(drop_out_db[i][j + 6].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 6].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 6.25 && quant < 12.5) { // Single Source
					double b(0.0), l(6.25), h(12.5);
					double slope((atof(drop_out_db[i][j + 9].c_str()) - atof(drop_out_db[i][j + 8].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 8].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}
				else if (drop_out_db[i][j] == "HOM1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 12.5 && quant < 25) {
					double b(0.0), l(12.5), h(25);
					double slope((atof(drop_out_db[i][j + 2].c_str()) - atof(drop_out_db[i][j + 9].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 9].c_str()) - slope * l;
					tc.PHOM1 = slope * quant + b;
					tc.PHOM0 = full - tc.PHOM1;

					if (term == 'n') {
						tc.PN_PHOM1 = tc.PHOM1;
						tc.PN_PHOM0 = tc.PHOM0;
					}
					else if (term == 'd') {
						tc.PD_PHOM1 = tc.PHOM1;
						tc.PD_PHOM0 = tc.PHOM0;
					}
				}

				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 25 && quant < 50) {
					double b(0.0), l(25.0), h(50.0);
					double slope((atof(drop_out_db[i][j + 3].c_str()) - atof(drop_out_db[i][j + 2].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 2].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 50 && quant < 100) {
					double b(0.0), l(50.0), h(100.0);
					double slope((atof(drop_out_db[i][j + 4].c_str()) - atof(drop_out_db[i][j + 3].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 3].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 100 && quant < 150) {
					double b(0.0), l(100.0), h(150.0);
					double slope((atof(drop_out_db[i][j + 5].c_str()) - atof(drop_out_db[i][j + 4].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 4].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 150 && quant < 250) {
					double b(0.0), l(150.0), h(250.0);
					double slope((atof(drop_out_db[i][j + 6].c_str()) - atof(drop_out_db[i][j + 5].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 5].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 250 && quant <= 500) {
					double b(0.0), l(250.0), h(500.0);
					double slope((atof(drop_out_db[i][j + 7].c_str()) - atof(drop_out_db[i][j + 6].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 6].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 6.25 && quant < 12.5) { // Single Source
					double b(0.0), l(6.25), h(12.5);
					double slope((atof(drop_out_db[i][j + 9].c_str()) - atof(drop_out_db[i][j + 8].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 8].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}
				else if (drop_out_db[i][j] == "HET1-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 12.5 && quant < 25) {
					double b(0.0), l(12.5), h(25.0);
					double slope((atof(drop_out_db[i][j + 2].c_str()) - atof(drop_out_db[i][j + 9].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 9].c_str()) - slope * l;
					tc.PHET1 = slope * quant + b;

					if (term == 'n') tc.PN_PHET1 = tc.PHET1;
					else if (term == 'd') tc.PD_PHET1 = tc.PHET1;
				}

				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 25 && quant < 50) {
					double b(0.0), l(25.0), h(50.0);
					double slope((atof(drop_out_db[i][j + 3].c_str()) - atof(drop_out_db[i][j + 2].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 2].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 50 && quant < 100) {
					double b(0.0), l(50.0), h(100.0);
					double slope((atof(drop_out_db[i][j + 4].c_str()) - atof(drop_out_db[i][j + 3].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 3].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 100 && quant < 150) {
					double b(0.0), l(100.0), h(150.0);
					double slope((atof(drop_out_db[i][j + 5].c_str()) - atof(drop_out_db[i][j + 4].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 4].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 150 && quant < 250) {
					double b(0.0), l(150.0), h(250.0);
					double slope((atof(drop_out_db[i][j + 6].c_str()) - atof(drop_out_db[i][j + 5].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 5].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 250 && quant <= 500) {
					double b(0.0), l(250.0), h(500.0);
					double slope((atof(drop_out_db[i][j + 7].c_str()) - atof(drop_out_db[i][j + 6].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 6].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 6.25 && quant < 12.5) { // Single Source
					double b(0.0), l(6.25), h(12.5);
					double slope((atof(drop_out_db[i][j + 9].c_str()) - atof(drop_out_db[i][j + 8].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 8].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
				else if (drop_out_db[i][j] == "HET2-" + to_string(contributors) + "-" + dnd && drop_out_db[i][j + 1] == locus && quant >= 12.5 && quant < 25) {
					double b(0.0), l(12.5), h(25.0);
					double slope((atof(drop_out_db[i][j + 2].c_str()) - atof(drop_out_db[i][j + 9].c_str())) / (h - l));
					b = atof(drop_out_db[i][j + 9].c_str()) - slope * l;
					tc.PHET2 = slope * quant + b;
					tc.PHET0 = full - (tc.PHET1 + tc.PHET2);

					if (term == 'n') {
						tc.PN_PHET2 = tc.PHET2;
						tc.PN_PHET0 = tc.PHET0;
					}
					else if (term == 'd') {
						tc.PD_PHET2 = tc.PHET2;
						tc.PD_PHET0 = tc.PHET0;
					}
				}
			}
		}
	}

	void swap_drop_out_rates(char term) {
		if (term == 'n') {
			tc.PHET0 = tc.PN_PHET0;
			tc.PHET1 = tc.PN_PHET1;
			tc.PHET2 = tc.PN_PHET2;
			tc.PHOM0 = tc.PN_PHOM0;
			tc.PHOM1 = tc.PN_PHOM1;
		}
		else if (term == 'd') {
			tc.PHET0 = tc.PD_PHET0;
			tc.PHET1 = tc.PD_PHET1;
			tc.PHET2 = tc.PD_PHET2;
			tc.PHOM0 = tc.PD_PHOM0;
			tc.PHOM1 = tc.PD_PHOM1;
		}
	}

private:
	vector<Person> knowns_pn;
	vector<Person> knowns_pd;
	int unknowns_pn;
	int unknowns_pd;
	Genotypes genotypes;
	vector<vector<Allele>> replicates;
	string locus;
	string case_name;

	vector<vector<Person>> population_pn;
	vector<vector<Person>> population_pd;
	vector<vector<int>> indicies_pn;
	vector<vector<int>> indicies_pd;
	vector<Person> persons;
	double pn;
	double pd;
	double lr;

	vector<vector<Person>> pn_combinations;
	vector<vector<Person>> pd_combinations;
	vector<pair<double, vector<Person>>> numerators;
	vector<pair<double, vector<Person>>> denominators;

	Thread_Constants tc;

	friend Allele;
	friend Genotypes;
	friend Person;

	friend tuple<double, double, double> run_LAST(csv_database& allele_database, csv_database& drop_out_db, csv_database& drop_in_db);
};

// Returns double
double generate_wild_allele_freq(vector<Allele>& alleles) {
	double c(1.0);
	for (int i(0); i < alleles.size(); ++i) {
		c -= alleles[i].freq;
	}
	if (c < 0)
		c = MINIMUM_WILD_FREQUENCY;
	return c;
}

// Returns double
double get_input_freq(csv_database& db, string& locus, double length, int race) {
	for (int i(0); i < db.size(); ++i) {
		for (int j(0); j < db[i].size(); ++j) {
			if (db[i][j] == locus && atof(db[i][j + 1].c_str()) == length && race == 1)
				return atof(db[i][j + 2].c_str());
			else if (db[i][j] == locus && atof(db[i][j + 1].c_str()) == length && race == 2)
				return atof(db[i][j + 3].c_str());
			else if (db[i][j] == locus && atof(db[i][j + 1].c_str()) == length && race == 3)
				return atof(db[i][j + 4].c_str());
			else if (db[i][j] == locus && atof(db[i][j + 1].c_str()) == length && race == 4)
				return atof(db[i][j + 5].c_str());
		}
	}
	cout << "Frequency not found - LR given will most likely be wrong!" << endl;
	return 0.0;
}

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

// Returns double (0.0 or 1.0)
double read_csv_x_rows(csv_database& db, string db_name, int rows) {		// Reads x number of rows rows only
	ifstream ifs(db_name);
	String csvLine;
	int i(1);
	while (getline(ifs, csvLine)) {
		istringstream csvStream(csvLine);
		csv_row csvRow;
		String csvCol;
		while (getline(csvStream, csvCol, ','))
			csvRow.push_back(csvCol);
		db.push_back(csvRow);
		if (i == rows) break;
		++i;
	}
	ifs.close();
	return true;
}

// Returns bool
bool get_data(csv_database& evidence_db, csv_database& allele_db, vector<vector<Allele>>& replicates, vector<Allele>& distinct_alleles, vector<Person>& knowns_pn, vector<Person>& knowns_pd, int& unknowns_pn, int& unknowns_pd, double& quant, string& case_name, string& locus, int& contributors_pn, int& contributors_pd, int& contributors, Thread_Constants& tc) {
	vector<Allele> one_replicate;
	vector<Allele> pn_tmp;			// Temporary container for numerator known alleles
	vector<Allele> pd_tmp;			// Temporary container for denominator known alleles

	for (int i(0); i < evidence_db[0].size(); ++i) {
		if (evidence_db[0][i] == "Case Name" || evidence_db[0][i] == "Name" || evidence_db[0][i] == "ID")
			case_name = evidence_db[1][i];
		else if (evidence_db[0][i] == "Locus" || evidence_db[0][i] == "Loci")
			locus = evidence_db[1][i];
		else if (evidence_db[0][i] == "Unknowns Pn" || evidence_db[0][i] == "Unknown Pn")
			unknowns_pn = atoi(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Unknowns Pd" || evidence_db[0][i] == "Unknown Pd")
			unknowns_pd = atoi(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Quant")
			quant = atof(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Product Group")
			PRODUCT_GROUP = atof(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Contributors Pn" || evidence_db[0][i] == "Contributor Pn")
			contributors_pn = atof(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Contributors Pd" || evidence_db[0][i] == "Contributor Pd")
			contributors_pd = atof(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Contributors" || evidence_db[0][i] == "Contributor")
			contributors = atof(evidence_db[1][i].c_str());
		else if (evidence_db[0][i] == "Run Type" || evidence_db[0][i] == "Run Speed" || evidence_db[0][i] == "Speed") {
			if (evidence_db[1][i] == "Full" || evidence_db[1][i] == "F" || evidence_db[1][i] == "f" || evidence_db[1][i] == "FULL")
				tc.RUN_TYPE = 'F';
			else if (evidence_db[1][i] == "Mod" || evidence_db[1][i] == "M" || evidence_db[1][i] == "m" || evidence_db[1][i] == "MOD")
				tc.RUN_TYPE = 'M';
			else if (evidence_db[1][i] == "Comb" || evidence_db[1][i] == "C" || evidence_db[1][i] == "c" || evidence_db[1][i] == "COMB")
				tc.RUN_TYPE = 'C';
			else
				tc.RUN_TYPE = 'E';
		}
		else if (evidence_db[0][i] == "D/ND" || evidence_db[0][i] == "Deducible") {
			if (evidence_db[1][i] == "d" || evidence_db[1][i] == "D" || evidence_db[1][i] == "yes")
				DEDUCIBLE = true;
			else if (evidence_db[1][i] == "nd" || evidence_db[1][i] == "ND" || evidence_db[1][i] == "no")
				DEDUCIBLE = false;
			else {
				cout << "Error(s) in the D/ND column or another column!" << endl;
				return false;
			}
		}
		else if (evidence_db[0][i] == "Known Pn" || evidence_db[0][i] == "Knowns Pn") {
			if (locus == "W") {
				cout << "Move '" << evidence_db[0][i] << "' column to the right of the 'Locus' column in the case file.";
				return false;
			}
			istringstream csvCell(evidence_db[1][i]);
			string sub_cell;
			while (getline(csvCell, sub_cell, ';')) {
				Allele tmp_a(locus, atof(sub_cell.c_str()), get_input_freq(allele_db, locus, atof(sub_cell.c_str()), tc.RACE));
				if (tmp_a.length != -1) {
					pn_tmp.push_back(tmp_a);
				}
			}
		}
		else if (evidence_db[0][i] == "Known Pd" || evidence_db[0][i] == "Knowns Pd") {
			if (locus == "W") {
				cout << "Move '" << evidence_db[0][i] << "' column to the right of the 'Locus' column in the case file.";
				return false;
			}
			istringstream csvCell(evidence_db[1][i]);
			string sub_cell;
			while (getline(csvCell, sub_cell, ';')) {
				Allele tmp_a(locus, atof(sub_cell.c_str()), get_input_freq(allele_db, locus, atof(sub_cell.c_str()), tc.RACE));
				if (tmp_a.length != -1) {
					pd_tmp.push_back(tmp_a);
				}
			}
		}
		else if (evidence_db[0][i] == "REP" || evidence_db[0][i] == "rep" || evidence_db[0][i] == "Replicate") {
			if (locus == "W") {
				cout << "Move '" << evidence_db[0][i] << "' column to the right of the 'Locus' column in the case file.";
				return false;
			}
			istringstream csvCell(evidence_db[1][i]);
			string sub_cell;
			while (getline(csvCell, sub_cell, ';')) {
				Allele tmp_a(locus, atof(sub_cell.c_str()), get_input_freq(allele_db, locus, atof(sub_cell.c_str()), tc.RACE));
				if (tmp_a.length != -1) {
					one_replicate.push_back(tmp_a);
				}
			}
			if (one_replicate.size() != 0) {
				//sort(one_replicate.begin(), one_replicate.end());		// Not necessary
				replicates.push_back(one_replicate);
			}
			one_replicate.clear();
		}
		else if (evidence_db[0][i] == "Alleles" || evidence_db[0][i] == "Distinct Alleles") {
			if (locus == "W") {
				cout << "Move '" << evidence_db[0][i] << "' column to the right of the 'Locus' column in the case file.";
				return false;
			}
			istringstream csvCell(evidence_db[1][i]);
			string sub_cell;
			while (getline(csvCell, sub_cell, ';')) {
				Allele tmp_a(locus, atof(sub_cell.c_str()), get_input_freq(allele_db, locus, atof(sub_cell.c_str()), tc.RACE));
				if (tmp_a.length != -1) {
					distinct_alleles.push_back(tmp_a);
				}
			}
		}
	}
	double freq = generate_wild_allele_freq(distinct_alleles);
	Allele wild(locus + "-W", -1, freq);
	distinct_alleles.push_back(wild);

	for (int i(0); i < pn_tmp.size(); i += 2) {
		Person known(pn_tmp[i], pn_tmp[i + 1]);
		knowns_pn.push_back(known);
	}
	for (int i(0); i < pd_tmp.size(); i += 2) {
		Person known(pd_tmp[i], pd_tmp[i + 1]);
		knowns_pd.push_back(known);
	}
	return true;
}

// Returns double
tuple<double, double, double> run_LAST(csv_database& allele_database, csv_database& drop_out_db, csv_database& drop_in_db) {
	Thread_Constants tc;
	tc.RACE = RACE;
	++RACE;
	if (RACE == 5) RACE = 1;

	double low_copy_pC0(0.0), low_copy_pC1(0.0), low_copy_pC2(0.0), high_copy_pC0(0.0), high_copy_pC1(0.0), high_copy_pC2(0.0), theta(0.0), min_w(0.0);
	for (int i(0); i < drop_in_db.size(); ++i) {
		for (int j(0); j < drop_in_db[i].size(); ++j) {
			if (drop_in_db[i][j] == "LPC0")
				low_copy_pC0 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "LPC1")
				low_copy_pC1 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "LPC2")
				low_copy_pC2 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "HPC0")
				high_copy_pC0 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "HPC1")
				high_copy_pC1 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "HPC2")
				high_copy_pC2 = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "THETA")
				theta = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "B-MIN-WILD-FREQ" && tc.RACE == 1)
				min_w = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "C-MIN-WILD-FREQ" && tc.RACE == 2)
				min_w = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "H-MIN-WILD-FREQ" && tc.RACE == 3)
				min_w = atof(drop_in_db[i][j + 1].c_str());
			else if (drop_in_db[i][j] == "A-MIN-WILD-FREQ" && tc.RACE == 4)
				min_w = atof(drop_in_db[i][j + 1].c_str());
		}
	}
	HOM_CONST = theta;
	MINIMUM_WILD_FREQUENCY = min_w;		// Setting minimum wild frequency in case "w" becomes negative

	csv_database evidence_database;
	read_csv_x_rows(evidence_database, "Evidence_" + to_string(FILE_NUM) + ".csv", 2);

	vector<Person> knowns_pn, knowns_pd;
	int unknowns_pn(0), unknowns_pd(0), contributors_pn(-1), contributors_pd(-1), contributors(-1);
	vector<vector<Allele>> replicates;
	vector<Allele> alleles;
	double quant(100.0);
	string dnd("ND"), case_name("<Name>"), locus("W");

	bool data_done = get_data(evidence_database, allele_database, replicates, alleles, knowns_pn, knowns_pd, unknowns_pn, unknowns_pd, quant, case_name, locus, contributors_pn, contributors_pd, contributors, tc);

	if (contributors != -1 && contributors_pn == -1) contributors_pn = contributors;
	if (contributors != -1 && contributors_pd == -1) contributors_pd = contributors;

	if (contributors_pn == -1) contributors_pn = knowns_pn.size() + unknowns_pn;
	if (contributors_pd == -1) contributors_pd = knowns_pd.size() + unknowns_pd;
	CONTRIBUTORS_PN = contributors_pn;
	CONTRIBUTORS_PD = contributors_pd;

	if (unknowns_pn == 0) unknowns_pn = contributors_pn - knowns_pn.size();
	if (unknowns_pd == 0) unknowns_pd = contributors_pd - knowns_pd.size();

	if (!data_done) {
		cout << "Error in case file. LR given will be 0." << endl;
		return tuple<double, double, double>(-1.0, -1.0, -1.0);		// Error in data
	}

	if (quant <= 100) {					// QUANT of 100 will be run as low copy
		PC0 = low_copy_pC0;
		PC1 = low_copy_pC1;
		PC2 = low_copy_pC2;
	}
	else {
		PC0 = high_copy_pC0;
		PC1 = high_copy_pC1;
		PC2 = high_copy_pC2;
	}

	if (quant < 6.25 && contributors_pn == 1 && contributors_pd == 1) quant = 6.25;
	else if (quant < 25.0 && (contributors_pn > 1 || contributors_pd > 1)) quant = 25.0;
	else if (quant > 500) quant = 500;

	QUANT = quant;

	Genotypes genotypes(alleles);
	Report report(drop_out_db, knowns_pn, unknowns_pn, knowns_pd, unknowns_pd, genotypes, replicates, case_name, locus, tc);

	return tuple<double, double, double>(report.lr, report.pn, report.pd);
}

int main() {
	Timer t, c;

	cout << "///////////////////////////////////////////////////////" << endl;
	cout << "//" << endl;
	cout << "//\t\tWELCOME TO reQBT!" << endl;
	cout << "//" << endl;
	cout << "//" << endl;

	const char* months[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };
	SYSTEMTIME st;
	GetSystemTime(&st);
	ostringstream ss;
	ss << st.wHour << "-" << st.wMinute << "-" << st.wSecond << "-" << st.wMilliseconds <<
		"-" << months[st.wMonth - 1] <<
		"-" << st.wDay << "-" << st.wYear % 1000;
	DATE_TIME = ss.str();
	_mkdir(DATE_TIME.c_str());

	string s_allele_freq = "copy Allele_Frequencies.csv " + DATE_TIME;
	system(s_allele_freq.c_str());
	string s_drop_out_rates = "copy Drop_Out_Rates.csv " + DATE_TIME;
	system(s_drop_out_rates.c_str());
	string s_drop_in_rates = "copy Drop_In_Rates.csv " + DATE_TIME;
	system(s_drop_in_rates.c_str());
	string s_case = "copy case.csv " + DATE_TIME;
	system(s_case.c_str());

	ofstream ofs(DATE_TIME + "/output.csv");
	ofs.close();
	ofs.open(DATE_TIME + "/output.csv", ios::app);
	ofs << "Case Name,Locus,Race,,LR,Pn,Pd,,Time (seconds),,Pn HET0,Pn HET1,Pn HET2,Pn HOM0,Pn HOM1,,Pd HET0,Pd HET1,Pd HET2,Pd HOM0,Pd HOM1" << endl;

	//double LR_b(1.0), LR_c(1.0), LR_h(1.0), LR_a(1.0);
	tuple<double, double, double> B(1.0, 1.0, 1.0), C(1.0, 1.0, 1.0), H(1.0, 1.0, 1.0), A(1.0, 1.0, 1.0);

	csv_database allele_database(6);
	read_csv(allele_database, "Allele_Frequencies.csv");
	csv_database drop_out_rates_database(8);
	read_csv(drop_out_rates_database, "Drop_Out_Rates.csv");
	csv_database drop_in_rates_database(2);
	read_csv(drop_in_rates_database, "Drop_In_Rates.csv");
	vector<double> LRs;
	bool run(true);

	while (true && run) {
		ifstream ifs("Evidence_" + to_string(FILE_NUM) + ".csv");
		if (ifs) {

			cout << "Running Analysis for Black for Evidence " << FILE_NUM << endl;
			std::future<tuple<double, double, double>> ret_b = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//std::future<double> ret_b = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//LR_b *= run_LAST(allele_database, drop_out_rates_database, drop_in_rates_database);
			//++RACE;

			cout << "Running Analysis for Caucasian for Evidence " << FILE_NUM << endl;
			std::future<tuple<double, double, double>> ret_c = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//std::future<double> ret_c = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//LR_c *= run_LAST(allele_database, drop_out_rates_database, drop_in_rates_database);
			//++RACE;

			cout << "Running Analysis for Hispanic for Evidence " << FILE_NUM << endl;
			std::future<tuple<double, double, double>> ret_h = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//std::future<double> ret_h = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//LR_h *= run_LAST(allele_database, drop_out_rates_database, drop_in_rates_database);
			//++RACE;

			cout << "Running Analysis for Asian for Evidence " << FILE_NUM << endl;
			std::future<tuple<double, double, double>> ret_a = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//std::future<double> ret_a = std::async(&run_LAST, allele_database, drop_out_rates_database, drop_in_rates_database);
			//LR_a *= run_LAST(allele_database, drop_out_rates_database, drop_in_rates_database);
			//RACE = 1;

			//LR_b *= ret_b.get();
			//LR_c *= ret_c.get();
			//LR_h *= ret_h.get();
			//LR_a *= ret_a.get();
			tuple<double, double, double> rb = ret_b.get();
			tuple<double, double, double> rc = ret_c.get();
			tuple<double, double, double> rh = ret_h.get();
			tuple<double, double, double> ra = ret_a.get();

			std::get<0>(B) *= std::get<0>(rb);
			std::get<1>(B) *= std::get<1>(rb);
			std::get<2>(B) *= std::get<2>(rb);
			std::get<0>(C) *= std::get<0>(rc);
			std::get<1>(C) *= std::get<1>(rc);
			std::get<2>(C) *= std::get<2>(rc);
			std::get<0>(H) *= std::get<0>(rh);
			std::get<1>(H) *= std::get<1>(rh);
			std::get<2>(H) *= std::get<2>(rh);
			std::get<0>(A) *= std::get<0>(ra);
			std::get<1>(A) *= std::get<1>(ra);
			std::get<2>(A) *= std::get<2>(ra);
		}
		else
			run = false;

		if (FILE_NUM % PRODUCT_GROUP == 0) {		// Calculate likelihood ratio once PRODUCT_GROUP-number of files (or loci) has been evaluated - this is one full case
													//cout << "\nBlack: " << LR_b << " Caucasian: " << LR_c << " Hispanic: " << LR_h << " Asian: " << LR_a << endl;
													//ofs << "\nOverall LR,," << PRODUCT_GROUP << ",Files" << endl;
													//ofs << "\nBlack:," << LR_b << ",,Caucasian:," << LR_c << ",,Hispanic:," << LR_h << ",,Asian:," << LR_a << endl << endl;
			cout << endl;
			cout << "Overall: " << PRODUCT_GROUP << " Files - BLACK - LR: " << std::get<0>(B) << " - Pn: " << std::get<1>(B) << " - Pd: " << std::get<2>(B) << endl;
			cout << "Overall: " << PRODUCT_GROUP << " Files - CAUCASIAN - LR: " << std::get<0>(C) << " - Pn: " << std::get<1>(C) << " - Pd: " << std::get<2>(C) << endl;
			cout << "Overall: " << PRODUCT_GROUP << " Files - HISPANIC - LR: " << std::get<0>(H) << " - Pn: " << std::get<1>(H) << " - Pd: " << std::get<2>(H) << endl;
			cout << "Overall: " << PRODUCT_GROUP << " Files - ASIAN: - LR: " << std::get<0>(A) << " - Pn: " << std::get<1>(A) << " - Pd: " << std::get<2>(A) << endl;
			double case_time = c.elapsed();
			ofs << "Overall:," << PRODUCT_GROUP << " Files,BLACK,," << std::get<0>(B) << "," << std::get<1>(B) << "," << std::get<2>(B) << ",," << case_time << endl;
			ofs << "Overall:," << PRODUCT_GROUP << " Files,CAUCASIAN,," << std::get<0>(C) << "," << std::get<1>(C) << "," << std::get<2>(C) << ",," << case_time << endl;
			ofs << "Overall:," << PRODUCT_GROUP << " Files,HISPANIC,," << std::get<0>(H) << "," << std::get<1>(H) << "," << std::get<2>(H) << ",," << case_time << endl;
			ofs << "Overall:," << PRODUCT_GROUP << " Files,ASIAN,," << std::get<0>(A) << "," << std::get<1>(A) << "," << std::get<2>(A) << ",," << case_time << endl;
			//LR_b = 1.0;
			//LR_c = 1.0;
			//LR_h = 1.0;
			//LR_a = 1.0;
			std::get<0>(B) = 1.0;
			std::get<1>(B) = 1.0;
			std::get<2>(B) = 1.0;
			std::get<0>(C) = 1.0;
			std::get<1>(C) = 1.0;
			std::get<2>(C) = 1.0;
			std::get<0>(H) = 1.0;
			std::get<1>(H) = 1.0;
			std::get<2>(H) = 1.0;
			std::get<0>(A) = 1.0;
			std::get<1>(A) = 1.0;
			std::get<2>(A) = 1.0;
			c.reset();
		}
		++FILE_NUM;
	}

	ofs << "\n,Total Time Elapsed:," << t.elapsed() << ",Second(s)" << endl;
	ofs.close();
	cout << "\n\t\treQBT Analysis Finished!\n" << endl;
	system("pause");
}