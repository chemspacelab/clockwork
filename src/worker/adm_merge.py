"""

merge jobs

"""

import gzip

import numpy as np
import scipy.special

import adm_fchl
import adm_molecule

def costfunc(clockres, n):
    """
    param:
        resolution list
        total number of torsions

    return
        cost

    """
    t = len(clockres)
    maxres = max(clockres)
    value = scipy.special.binom(n, t)
    no = 2**(maxres)
    value *= no
    return value


def info(jobstr):

    job = eval(jobstr)

    (midx, toridx, reslist), uniqueidx, energy, coord = job

    return midx, toridx, reslist, uniqueidx, energy, coord



def energyandfchl(jobs, atoms, decimals=1, threshold=0.98, reset_count=False, debug=False):
    """

    - sort energy
    - for all energy there are "unique", decimal=0.1
        - do qml

    """

    atoms = np.array(atoms)
    n_atoms = len(atoms)

    if type(jobs[0]) == str:
        jobs = [info(job) for job in jobs]

    energies = []
    costlist = []

    for job in jobs:

        energy = job[4]
        energy = float(energy)

        reslist = job[2]

        cost = costfunc(reslist, 30)

        energies.append(energy)
        costlist.append(cost)


    costlist = np.array(costlist)
    energies = np.round(energies, decimals=decimals)
    unique_energies = np.unique(energies)

    if debug:
        print("merge: found {:} energies".format(len(unique_energies)))

    keep = []
    keep_origin = []
    keep_merge = []

    for energy in unique_energies:

        jobidx, = np.where(energies == energy)

        if len(jobidx) == 1:

            keep.append(jobidx[0])
            keep_origin.append(jobidx[0])
            keep_merge.append(0)

            continue

        costs = costlist[jobidx]

        positions = []
        for i in jobidx:
            position = jobs[i][-1]
            position = np.array(position)
            position = position.reshape((n_atoms, 3))
            positions.append(position)


        if debug:
            print("merge: qml for {:} @ {:}".format(len(jobidx), energy))

        # Energy are the same, compare with FCHL
        representations = adm_fchl.get_representations_fchl(atoms, positions)
        representations = np.array(representations)

        similarity = adm_fchl.get_kernel_fchl(representations, representations)

        dealtwith = []

        for j, row in enumerate(similarity):

            if j in dealtwith: continue

            rowsim = np.argwhere(row >= threshold)
            rowsim = rowsim.flatten()
            rowcosts = costs[rowsim]

            # find the lowest cost, of all (incl dealtwith)
            costidx = np.argsort(rowcosts)
            keep_sim_cost = costidx[0]

            # want not dealt with index and lowest cost
            rowsim = [idx for idx in rowsim if idx not in dealtwith]
            costs_pos = costs[rowsim]
            costs_pos_low = np.argsort(costs_pos)[0]

            keep_sim_pos = rowsim[costs_pos_low]

            # convert similarity idx to jobidx
            keep_idx = jobidx[keep_sim_pos]
            keep_cost = jobidx[keep_sim_cost]

            keep.append(keep_idx)
            keep_origin.append(keep_cost)
            keep_merge.append(len(rowsim) - 1)

            dealtwith += list(rowsim)


    # correct jobs with origins
    rtnjobs = []

    for jidx, oidx, nmerge in zip(keep, keep_origin, keep_merge):

        job = jobs[jidx]
        job = list(job)

        if jidx != oidx:

            ojob = jobs[oidx]
            origin = ojob[2]
            job[2] = origin

        if reset_count:

            job[3] = nmerge

        else:

            job[3] += nmerge

        rtnjobs.append(job)

    return rtnjobs


def merge_energy(jobs, decimals=1):

    if type(jobs[0]) == str:
        jobs = [info(job) for job in jobs]

    energies = []
    costlist = []

    for job in jobs:

        energy = job[4]
        energy = float(energy)

        reslist = job[2]

        cost = costfunc(reslist, 30)

        energies.append(energy)
        costlist.append(cost)


    costlist = np.array(costlist)
    energies = np.round(energies, decimals=decimals)
    unique_energies = np.unique(energies)

    keep = []

    for energy in unique_energies:

        idx, = np.where(energies == energy)
        costs = costlist[idx]

        idxsort = np.argsort(costs)

        j = idxsort[0]
        keep.append(idx[j])


    jobs = [jobs[idx] for idx in keep]

    return jobs


def merge_fchl(jobs):
    """

    take a list of jobs and merge them

    """


    return jobs


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', type=str, help='', metavar='file')
    parser.add_argument('-s', '--sdf', type=str, help='', metavar='file')
    args = parser.parse_args()

    with open(args.filename) as f:
        lines = f.readlines()


    ext = args.sdf.split(".")[-1]

    if ext == "sdf":

        moldb = adm_molecule.load_sdf(args.sdf)

    elif ext == "gz":

        inf = gzip.open(args.sdf)
        moldb = adm_molecule.load_sdf_file(inf)

    else:

        print("unknown format")
        quit()


    # Select mol and find atoms
    job = info(lines[0])
    molidx = job[0]
    mol = moldb[molidx]

    atoms = [atom for atom in mol.GetAtoms()]
    atoms_str = [atom.GetSymbol() for atom in atoms]
    atoms_int = [atom.GetAtomicNum() for atom in atoms]


    jobs = energyandfchl(lines, atoms_int, reset_count=True)

    for job in jobs:
        print(job)

    return


if __name__ == '__main__':
    main()

