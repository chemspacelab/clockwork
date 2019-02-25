#include <iostream>
#include <vector>
#include <map>

//#define DEBUG

class Archive {
public:
	std::vector<unsigned int> molecule_ids;
	std::map<unsigned int, unsigned int> number_of_atoms;
	std::map<unsigned int, unsigned int*> element_numbers;
	std::map<unsigned int, double*> coordinates;
	std::map<unsigned int, unsigned int> number_of_bonds;
	std::map<unsigned int, unsigned int*> bonds;
	std::map<unsigned int, unsigned int*> bond_orders;
	std::map<unsigned int, unsigned int> number_of_dihedrals;
	std::map<unsigned int, unsigned int*> dihedrals;
	int read_archive(const char * filename);
};

int Archive::read_archive(const char * filename) {
	std::ifstream fh;
	fh.open(filename, std::ios::in | std::ios::binary);

	if (!fh.is_open()) {
		std::cout << "Unable to open file " << filename << std::endl;
		return 1;
	}

	char * buffer1 = new char[1];
	char * buffer4 = new char[4];
	char * buffer8 = new char[8];
	unsigned int atomcount, elementid, bondcount, dihedralcount;
	unsigned int molecule_id;
	
	while (fh.read(buffer4, 4) && !fh.eof()) {
		// read molecule id
		molecule_id = *(unsigned int*)buffer4;
		#ifdef DEBUG
			std::cout << "Reading molecule " << molecule_id << std::endl;
		#endif
		this->molecule_ids.push_back(molecule_id);

		// read number of atoms
		fh.read(buffer1, 1);
		atomcount = *(unsigned int*)buffer1;
		#ifdef DEBUG
			std::cout << " > has " << atomcount << " atoms" << std::endl;
		#endif
		number_of_atoms.insert(std::pair<unsigned int, unsigned int>(molecule_id, atomcount));

		// read element numbers
		unsigned int * numbers = new unsigned int[atomcount];
		#ifdef DEBUG
			std::cout << " > atoms have these atomic numbers:" << std::endl << " >> ";
		#endif
		for (unsigned int atom = 0; atom < atomcount; atom++) {
			fh.read(buffer1, 1);
			elementid = *(unsigned int*)buffer1;
			#ifdef DEBUG
				std::cout << elementid << ", ";
			#endif
			numbers[atom] = elementid;
		}
		#ifdef DEBUG
			std::cout << std::endl;
		#endif
		element_numbers.insert(std::pair<unsigned int, unsigned int*>(molecule_id, numbers));

		// read nuclear coordinates
		double * mcoordinates = new double[atomcount * 3];
		mcoordinates[0] = 0;
		mcoordinates[1] = 0;
		mcoordinates[2] = 0;
		#ifdef DEBUG
			std::cout << " > atoms have these coordinates: " << std::endl;
		#endif
		for (unsigned int atom = 1; atom < atomcount; atom++) {
			for (unsigned int dimension = 0; dimension < 3; dimension++) {
				fh.read(buffer8, 8);
				mcoordinates[atom * 3 + dimension] = *(double*)buffer8;
			}
			#ifdef DEBUG
				std::cout << " > atom " << atom << ": " << mcoordinates[atom * 3 + 0] << ", " << mcoordinates[atom * 3 + 1]<< ", " << mcoordinates[atom * 3 + 2] << std::endl;
			#endif
		}
		coordinates.insert(std::pair<unsigned int, double*>(molecule_id, mcoordinates));

		// read number of bonds
		fh.read(buffer1, 1);
		bondcount = *(unsigned int*)buffer1;
		#ifdef DEBUG
			std::cout << " > has " << bondcount << " bonds" << std::endl;
		#endif
		number_of_bonds.insert(std::pair<unsigned int, unsigned int>(molecule_id, bondcount));

		// read bond indices
		unsigned int * bondindices = new unsigned int[bondcount * 2];
		#ifdef DEBUG
			std::cout << " > bonds have these pairs: " << std::endl << " >> ";
		#endif
		for (unsigned int bond = 0; bond < bondcount; bond++) {
			fh.read(buffer1, 1);
			bondindices[bond * 2] = *(unsigned int*)buffer1;
			fh.read(buffer1, 1);
			bondindices[bond * 2 + 1] = *(unsigned int*)buffer1;
			#ifdef DEBUG
				std::cout << bondindices[bond * 2] << "-" << bondindices[bond * 2 + 1] << ", ";
			#endif
		}
		#ifdef DEBUG
			std::cout << std::endl;
		#endif
		bonds.insert(std::pair<unsigned int, unsigned int*>(molecule_id, bondindices));

		// read bond orders
		unsigned int * bondorders = new unsigned int[bondcount];
		#ifdef DEBUG
			std::cout << " > bonds have these orders: " << std::endl << " >> ";
		#endif
		for (unsigned int bond = 0; bond < bondcount; bond++) {
			fh.read(buffer1, 1);
			bondorders[bond] = *(unsigned int*)buffer1;
			#ifdef DEBUG
				std::cout << bondorders[bond] << ", ";
			#endif
		}
		#ifdef DEBUG
			std::cout << std::endl;
		#endif
		bond_orders.insert(std::pair<unsigned int, unsigned int*>(molecule_id, bondorders));

		// read number of dihedrals
		fh.read(buffer1, 1);
		dihedralcount = *(unsigned int*)buffer1;
		#ifdef DEBUG
			std::cout << " > has " << dihedralcount << " dihedrals" << std::endl;
		#endif
		number_of_dihedrals.insert(std::pair<unsigned int, unsigned int>(molecule_id, dihedralcount));

		// read dihedral indices
		unsigned int * dihedralindices = new unsigned int[dihedralcount * 4];
		#ifdef DEBUG
			std::cout << " > dihedrals have these tuples: " << std::endl << " >> ";
		#endif
		for (unsigned int dihedral = 0; dihedral < dihedralcount; dihedral++) {
			for (unsigned int pos = 0; pos < 4; pos++) {
				fh.read(buffer1, 1);
				dihedralindices[dihedral * 4 + pos] = *(unsigned int*)buffer1;	
			}
			
			#ifdef DEBUG
				std::cout << dihedralindices[dihedral * 4 + 0] << "-" << dihedralindices[dihedral * 4 + 1] << "-" << dihedralindices[dihedral * 4 + 2] << "-" << dihedralindices[dihedral * 4 + 3] << ", ";
			#endif
		}
		#ifdef DEBUG
			std::cout << std::endl;
		#endif
		dihedrals.insert(std::pair<unsigned int, unsigned int*>(molecule_id, dihedralindices));
	}

	delete[] buffer1;
	delete[] buffer4;
	delete[] buffer8;

	fh.close();
	return 0;
}