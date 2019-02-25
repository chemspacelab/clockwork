#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <omp.h>
#include <ctime>
#include <cstring>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/obconversion.h>

#include "io.cpp"
#include "redis.cpp"
#include "kabsch.h"
#include "molecule.cpp"

int main(int argc,char **argv)
{
	omp_set_num_threads(1);

	Archive archive;
	archive.read_archive(argv[1]);
	std::cout << "Reading archive complete." << std::endl;
	
	redisContext *context = redis_connect();
	std::cout << "Connected to redis." << std::endl;

	unsigned char * payload = 0;
	size_t payloadsize;
	time_t packagestart;
	while (redis_fetch(context, &payload, &payloadsize) == 0) {
		packagestart = std::time(NULL);

		example_worker(archive, archive.molecule_ids[0]);

		if (redis_notify(context, payload, payloadsize, packagestart) != 0) {
			std::cout << "Problems in communicating to REDIS. Aborting." << std::endl;
		}

		delete[] payload;
	}

	redisFree(context);
	return 0;
}
