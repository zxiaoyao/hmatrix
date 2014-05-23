#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>


using namespace std;


/**
 * Load hbond info between every two conformers from "hb.dat", which is a binary file.
 * The first 4 bytes is an integer indicating the number of conformers, followed by n*n matrix.
 * The matrix element (i, j) being 1 means there is hbond between conformer i and conformer j,
 * with conformer i being the donor and conformer j the accepter. 0 means no hbond between conformer i and j.
 *
 * @param hinter  the matrix stores the hbond info between every two conformers.
 * @param hb_dat_file the name of the binary file which saves the hbond info between every two conformers.
 */
void load_hb_data(vector<vector<char> > &hinter, string hb_dat_file="hb.dat") {
	ifstream iHbfile(hb_dat_file, ios::in | ios::binary);
	int n_conf;

	iHbfile.read((char *) &n_conf, sizeof(int));
	hinter.resize(n_conf);
	for (int i=0; i<n_conf; i++) hinter[i].resize(n_conf);

	for (int i=0; i<n_conf; i++) {
		for (int j=0; j<n_conf; j++) {
			iHbfile.read((char *) &hinter[i][j], 1);
		}
	}
	iHbfile.close();
}


/**
 * Load all the microstates from file ms.dat, which is a binary file.
 * The first 4 bytes is an integer giving the number of key residues whose microstates are saved.
 * Followd by the names of the key residues in a format like ASPA0085.
 * Then comes a set of microstate records, which have the conformer id number of the chosen conformer
 * of each residue in each microstate.
 *
 * @param resName the names of key residues whose microstates are saved and only which are shown in the hbond network.
 * @param resinter  the number of hbonding between two residues among all the microstates.
 * @param hinter  the matrix stores the hbond info between every two conformers, in which the element is either 0 or 1.
 * @param ms_dat_file the name of the binary file which saves all the microstates.
 */
int load_ms_data(vector<string> &resName, vector<vector<int> > &resinter,
				vector<vector<char> > &hinter, string ms_dat_file="ms.dat") {

	ifstream iMsfile(ms_dat_file, ios::in | ios::binary);

	// n_spe gives the number of special residues
	int n_spe;
	iMsfile.read((char *) &n_spe, 4);
	cout << "There are " << n_spe << " key residues." << endl;

	// resName are the names of the special residues
	resName.clear();
	char buffer[9];
	for (int i_spe=0; i_spe<n_spe; i_spe++) {
		iMsfile.read(buffer, 8);
		buffer[8] = '\0';
		resName.push_back(string(buffer));
	}

	resinter.clear();
	resinter.resize(n_spe);
	for (size_t i=0; i<n_spe; i++) resinter[i].resize(n_spe);

	int totalState = 0;
	int totalRecords = 0;

	const int PRINT_INTERVAL = 10000;

	unsigned short conf;
	double H_state, Hsq_state;
	int count;
	vector <unsigned short> iconf;

	// load the microstates one by one to save memory, instead of loading all of them at once.
	while (iMsfile.read((char *) &conf, 2)) {
		iconf.push_back(conf);
		for (int i_spe=1; i_spe<n_spe; i_spe++) {
			iMsfile.read((char *) &conf, 2);
			iconf.push_back(conf);
		}
		iMsfile.read((char *) &H_state, 8);
		iMsfile.read((char *) &Hsq_state, 8);
		iMsfile.read((char *) &count, 4);

		totalRecords++;
		totalState += count;
		if (totalRecords % PRINT_INTERVAL == 0) cout << totalRecords << " records have been loaded." << endl;

		for (int i_res=0; i_res<n_spe; i_res++) {
			for (int j_res=0; j_res<n_spe; j_res++) {
				if (j_res == i_res) continue;
				resinter[i_res][j_res] += hinter[iconf[i_res]][iconf[j_res]] * count;
			}
		}
		iconf.clear();
	}
	iMsfile.close();

	cout << totalRecords << " records and " << totalState << " states has been loaded." << endl;

	return totalState;
}


/**
 * Write the probability of hbond between two residues into file "hb.txt".
 * @param resName the names of key residues whose microstates are saved and only which are shown in the hbond network.
 * @param resinter the number of hbonding between two residues among all the microstates.
 * @param totalState the number of all the microstates.
 * @param threshold the hbond connection is written if the probability is larger or equal to this threshold.
 * @param hb_txt_file the name of the text file which outputs the probability of hbonding between two residues.
 */
void output_hb_txt(vector<string> &resName, vector<vector<int> > &resinter, int totalState, float threshold=0.001,
					string hb_txt_file="hb.txt") {

	ofstream ofile(hb_txt_file, ios::out);
	for (int i=0; i<resinter.size(); i++) {
		for (int j=0; j<resinter.size(); j++) {
			if ((float) resinter[i][j]/totalState >= threshold) {
				ofile << resName[i] << '\t' << resName[j] << '\t' << fixed << setprecision(3)
					<< ((float) resinter[i][j]/totalState) << endl;;
			}
		}
	}
	ofile.close();
}


int main()
{
	time_t time_start = time(NULL);

	vector<vector<char> > hinter;
	load_hb_data(hinter);

	vector<string> resName;
	vector<vector<int> > resinter;
	int totalState = load_ms_data(resName, resinter, hinter);

	cout << "Total time to load ms.dat: " << time(NULL) - time_start << "." << endl;

	output_hb_txt(resName, resinter, totalState);
	cout << "Total time: " << time(NULL) - time_start << "." << endl;

    return 0;
}
