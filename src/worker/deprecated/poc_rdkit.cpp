#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

int main( int argc , char **argv ) {

	RDKit::ROMol *mol1 = RDKit::SmilesToMol( "Cc1ccccc1" );

	int N = mol1->getNumAtoms();

	std::cout << "num atoms " << N << "\n";

	return 0;
}
