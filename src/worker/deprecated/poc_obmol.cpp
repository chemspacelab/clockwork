
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>


int main(int argc,char **argv)
{
	// Init
	OpenBabel::OBMol mol;

	mol.BeginModify();

	OpenBabel::OBAtom obatom;

	int atom_no = 6;

	double x = 0.05;
	double y = 1.55;
	double z = 6.00;

	obatom.SetAtomicNum(atom_no);
	obatom.SetVector(x, y, z);

	std::cout << " " << obatom.GetX();
	std::cout << " " << obatom.GetY();
	std::cout << " " << obatom.GetZ();

	std::cout << "\n";

	mol.AddAtom(obatom);

	mol.EndModify();

	std::cout << " " << mol.GetCoordinates()[0];
	std::cout << " " << mol.GetCoordinates()[1];
	std::cout << " " << mol.GetCoordinates()[2];

	std::cout << "\n";

	return 0;
}
