
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>

#include <cmath> // std:pow
#include <algorithm>    // std::min

OpenBabel::OBMol make_molecule(Archive archive, int idx)
{

	// Init
	OpenBabel::OBMol mol;
	mol.BeginModify();

	unsigned int n_atoms = archive.number_of_atoms[idx];
	unsigned int n_bonds = archive.number_of_bonds[idx];
	// unsigned int n_dihedrals = archive.number_of_dihedrals[idx];

	auto coord = archive.coordinates[idx];

	for (unsigned int i=0; i<n_atoms; i++)
	{
		OpenBabel::OBAtom obatom;

		auto atom_no = archive.element_numbers[idx][i];

		auto x = coord[3*i+0];
		auto y = coord[3*i+1];
		auto z = coord[3*i+2];

		obatom.SetAtomicNum(atom_no);
		obatom.SetVector(x, y, z);

		mol.AddAtom(obatom);
	}

	for (unsigned int i=0; i<n_bonds; i++)
	{
		auto bond_i = archive.bonds[idx][i*2] + 1;
		auto bond_j = archive.bonds[idx][i*2+1] + 1;
		auto bond_order = archive.bond_orders[idx][i];
		mol.AddBond(bond_i, bond_j, bond_order);
	}

	mol.EndModify();

	return mol;
}


double add_to_torsion(
	OpenBabel::OBMol *mol,
	int a,
	int b,
	int c,
	int d,
	double angle)
{

	double current = mol->GetTorsion(a+1, b+1, c+1, d+1);
	current += angle;

	double radian = current * 3.14159 / 180.0;

	mol->SetTorsion(
		mol->GetAtom(a+1),
		mol->GetAtom(b+1),
		mol->GetAtom(c+1),
		mol->GetAtom(d+1),
		radian);

	angle = mol->GetTorsion(a+1, b+1, c+1, d+1);

	return angle;
}


int scan_torsion(int i, int j, int n)
{
	double start = 360.0/std::pow(2, n+1);
	double step = std::min((double) 360.0/std::pow(2, n-1), (double) 180.0);
	unsigned int n_steps = std::pow(2, n-1);

	std::cout << "n_steps " << n_steps << "\n";
	
	for (unsigned int k; k < n_steps; k++)
	{
		std::cout << k << " " << start + step*k << "\n";
	}

	return 0;
}


void example_worker(Archive archive, int mol_index)
{

	// Needed such that openbabel does not try to parallelise
	omp_set_num_threads(1);

	// OpenBabel::OBForceField * ff = OpenBabel::OBForceField::FindForceField("MMFF94");
	OpenBabel::OBForceField * ff = OpenBabel::OBForceField::FindForceField("UFF");
	auto constraints = OpenBabel::OBFFConstraints();

	double angle;
	int a, b, c, d;

	auto coord = archive.coordinates[mol_index];
	auto n_atoms = archive.number_of_atoms[mol_index];

	auto mol = make_molecule(archive, mol_index);

	int n_steps_per_gradient = 5;
	int n_steps_per_round = 100;
	int n_steps_per_round_constraint = 1000;
	double threshold = 1e-4;

	// ff->SetLogFile(0);
	// ff->SetLogLevel(OBFF_LOGLVL_HIGH);

	angle = 90.0;

	std::cout << "starting" << "\n";

	for (unsigned int n=0; n<100; n++)
	{
		angle = 90.0;

		int i = 2;
		a = archive.dihedrals[mol_index][4*i + 0];
		b = archive.dihedrals[mol_index][4*i + 1];
		c = archive.dihedrals[mol_index][4*i + 2];
		d = archive.dihedrals[mol_index][4*i + 3];

		angle = add_to_torsion(&mol, a, b, c, d, angle);

		ff->Setup(mol);

		double energyA = 10000;
		double energyB = 1000;
		int ni = 0;

		ff->SetConstraints(constraints);
		constraints.AddTorsionConstraint(a+1, b+1, c+1, d+1, angle);
		ff->ConjugateGradientsInitialize(n_steps_per_gradient, threshold);

		while(energyA-energyB>threshold && ni < n_steps_per_round_constraint)
		{
			ni += 1;
			energyA = energyB;
			ff->ConjugateGradientsTakeNSteps(n_steps_per_gradient);
			energyB = ff->Energy();
		}

		constraints.DeleteConstraint(constraints.Size()-1);
		ff->ConjugateGradientsInitialize(n_steps_per_gradient, threshold);

		energyA = 10000;
		energyB = 1000;
		ni = 0;

		while(energyA-energyB>threshold && ni < n_steps_per_round)
		{
			ni += 1;
			energyA = energyB;
			ff->ConjugateGradientsTakeNSteps(n_steps_per_gradient);
			energyB = ff->Energy();
		}
	
		double compare = kabsch::kabsch_rmsd(coord, mol.GetCoordinates(), n_atoms);

		std::cout << "member " << angle << " " << compare << " " << constraints.Size() << "\n";

	}

}
