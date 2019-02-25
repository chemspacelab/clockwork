#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>

#include "kabsch.h"


int main(int argc,char **argv)
{

	std::string filename_alpha = "xyz/wat_longer.xyz";
	std::string filename_beta = "xyz/wat_longer_moved.xyz";

	std::ifstream file_alpha (filename_alpha);
	std::ifstream file_beta (filename_beta);

	OpenBabel::OBMol mol_alpha;
	OpenBabel::OBMol mol_beta;

	OpenBabel::OBConversion conv;
	conv.SetInFormat("XYZ");

	conv.Read(&mol_alpha, &file_alpha);
	conv.Read(&mol_beta, &file_beta);

	auto coord_alpha = mol_alpha.GetCoordinates();
	auto coord_beta = mol_beta.GetCoordinates();

	unsigned int n_atoms = mol_alpha.NumAtoms();

	auto centroid_alpha = kabsch::centroid(coord_alpha, n_atoms);
	auto centroid_beta = kabsch::centroid(coord_beta, n_atoms);

	// Recenter
	for(unsigned int i = 0; i < n_atoms; i++) {
		for(unsigned int d = 0; d < 3; d++) {
			coord_alpha[3*i + d] -= centroid_alpha[d];
			coord_beta[3*i + d] -= centroid_beta[d];
		}
	}

	// rmsd
	std::cout << "rmsd: " << kabsch::rmsd(coord_alpha, coord_beta, n_atoms) << "\n";
	std::cout << "rotated rmsd: " << kabsch::kabsch_rmsd(coord_alpha, coord_beta, n_atoms) << "\n";

	delete[] centroid_alpha;
	delete[] centroid_beta;

	return 0;
}
