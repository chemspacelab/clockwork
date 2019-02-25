
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>

#include "io.cpp"
#include "kabsch.h"
#include "molecule.cpp"

int main(int argc,char **argv)
{
	Archive archive;
	archive.read_archive("../../fixtures/sample.archive");

	// unsigned int idx = archive.molecule_ids[0];

	example_worker(archive, 0);

	return 0;
}
