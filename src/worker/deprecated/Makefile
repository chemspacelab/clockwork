CC=g++
OPT= -O3
CFLAGS=-Wall -g -I ../dep/dist/include/openbabel-2.0/ -I ../dep/dist/include/ -fopenmp -std=c++11
LDFLAGS=-L ../dep/dist/lib/ -lgomp -lopenbabel
NETWORK=-lhiredis
LDPOSTFLAGS=-lcblas -llapacke -llapack -lblas -lgfortran -llammps
MATHFLAGS=

RDBASE=../dep/rdkit
LDRDKIT= -L $(RDBASE)/build/lib -I ${RDBASE}/Code
FLAGSRDKIT= \
	-lRDKitDepictor -lRDKitSubstructMatch -lRDKitSubgraphs \
	-lRDKitGraphMol -lRDKitDataStructs -lRDKitRDGeometryLib \
	-lRDKitRDGeneral -lRDKitSmilesParse -lRDKitFileParsers \
	-lRDKitDistGeomHelpers -lRDKitForceFieldHelpers -lRDKitForceField \
	-lRDKitMolTransforms

all: worker

worker: main.cpp
	$(CC) $(CFLAGS) $(MATHFLAGS) $(LDFLAGS) main.cpp -o worker $(LDPOSTFLAGS) $(NETWORK)

test: poc_rmsd.o poc_obmol.o poc_conformations.o poc_torstion_scan.o

poc_rmsd.o: molecule.cpp
poc_obmol.o: molecule.cpp
poc_conformations.o: molecule.cpp
poc_torsion_scan.o: molecule.cpp

poc_lammps.o: poc_lammps.cpp
	mpicxx -o poc_lammps.o poc_lammps.cpp -L ../dep/dist/lib/ -I ../dep/dist/include/ -llammps -lmpi

poc_lammps_simple.o: poc_lammps_simple.cpp
	mpicxx -o poc_lammps_simple.o poc_lammps_simple.cpp -L ../dep/dist/lib/ -L ../dep/lammps/src -I ../dep/dist/include/ -I ../dep/lammps/src -llammps -lmpi



poc_rdkit.o: poc_rdkit.cpp
	$(CC) \
	-I ${RDBASE}/Code \
	-L $(RDBASE)/build/lib \
	-o poc_rdkit.o poc_rdkit.cpp \
	-lRDKitDepictor -lRDKitSubstructMatch -lRDKitSubgraphs \
	-lRDKitGraphMol -lRDKitDataStructs -lRDKitRDGeometryLib \
	-lRDKitRDGeneral -lRDKitSmilesParse 

molecule_rdkit.o: molecule_rdkit.cpp
	$(CC) $(LDRDKIT) -o $@ $< $(FLAGSRDKIT)

%.o: %.cpp
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(MATHFLAGS) $< $(LDPOSTFLAGS)

clean:
	rm -rf *.o worker
