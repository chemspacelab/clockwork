
#include <iostream>

#include <GraphMol/GraphMol.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/MolOps.h>

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include <GraphMol/MolTransforms/MolTransforms.h>

#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/DistanceConstraint.h>
#include <ForceField/MMFF/AngleConstraint.h>
#include <ForceField/MMFF/TorsionConstraint.h>
#include <ForceField/MMFF/PositionConstraint.h>

#include <ForceField/ForceField.h>

#include <sstream>


std::string floattostring(double value)
{
	std::ostringstream strs;
	strs << value;
	std::string str = strs.str();
	return str;
}


void print_molecule(RDKit::ROMOL_SPTR mol, std::string name = "")
{
	// Print molecule in XYZ format

	unsigned int n_atoms = mol->getNumAtoms();
	RDKit::Conformer &conf = mol->getConformer();

	std::cout << n_atoms << "\n" << name << "\n";

	for (unsigned int i=0; i<n_atoms; i++)
	{
		RDKit::Atom *atom = mol->getAtomWithIdx(i);
		RDGeom::Point3D pos = conf.getAtomPos(i);
		std::cout << atom->getAtomicNum() << " ";
		std::cout << pos.x << " ";
		std::cout << pos.y << " ";
		std::cout << pos.z << " ";
		std::cout << "\n";
	}

}


double example_worker(int nexample)
{

	std::string sdf_filename = "sdf/pentane.sdf";

	// filename, sanitize, remove h
	RDKit::SDMolSupplier mol_supplier( sdf_filename, true, false);
	RDKit::ROMOL_SPTR mol( mol_supplier.next() );

	unsigned int n_atoms = mol->getNumAtoms();
	RDKit::Conformer &conf = mol->getConformer();

	double angle;
	double add_angle = 30.0;

	// Begin
	// angle = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);
	// print_molecule(mol, floattostring(angle));

	// Setup forcefield
	RDKit::MMFF::MMFFMolProperties *mmffMolProperties = new RDKit::MMFF::MMFFMolProperties(*mol);
    ForceFields::ForceField *ff = RDKit::MMFF::constructForceField(*mol, mmffMolProperties);
    ff->initialize();

	// std::cout << ff->contribs().size() << "\n";

	// angle = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);
	// print_molecule(mol, floattostring(angle));

	for (unsigned int n=0; n<nexample; n++)
	{

		angle = angle + add_angle;

		// Get current angle
		// double angle1 = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);

		// Set new angle
		MolTransforms::setDihedralDeg(conf, 0, 1, 5, 8, angle);

		// Set constraint
		ForceFields::MMFF::TorsionConstraintContrib *tc;
		tc = new ForceFields::MMFF::TorsionConstraintContrib(ff, 0, 1, 5, 8, angle, angle, 1.0e10);
		ff->contribs().push_back(ForceFields::ContribPtr(tc));

		// double angle2 = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);

		// Minimze
		// minimize(unsigned int maxIts = 200, double forceTol = 1e-4, double energyTol = 1e-6);
		ff->minimize(200, 1e-2, 1e-2);

		// double angle3 = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);

		// Remove last constraint
		ff->contribs().pop_back();

		// Minimize
		// ff->minimize(200, 1e-2, 1e-2);

		// double angle4 = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);

		// std::cout << "N " << n << " ";
		// std::cout << angle1 << " ";
		// std::cout << angle2 << " ";
		// std::cout << angle3 << " ";
		// std::cout << angle4 << " ";
		// std::cout << "\n";

	}

	// std::cout << ff->contribs().size() << "\n";

	// angle = MolTransforms::getDihedralDeg(conf, 0, 1, 5, 8);
	// print_molecule(mol, floattostring(angle));

	delete mmffMolProperties;
	// delete tc; // Delete via pop_back() ?

}


int main( int argc , char **argv ) {

	example_worker(1000);

	return 0;
}
