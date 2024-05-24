# Structure Suite YY
# this is a collection of tools for structure analysis of Molecular Dynamics Simulations results

# each atom's information is stored including atom_id, type, xyz position;
# the atom's index is based on the sequence that when they are read line by line from the lata file,
# NOT based on their atom_id (which can take more space because their
# atom_id not necessarily start from 1; also, when these atoms are throw to each sub-cell, we use their index, NOT atom id.
# when the connection graph is built, we also use their index, NOT atom ID.

from math import exp, log10
import sys, argparse, pdb
import random as rd
import numpy as np
from Atom import Atom
from SimulationBox import SimulationBox
from read import read_lata, read_velocity, readdx_upper, read_position_and_velocity, read_equi_position
from read import read_lata, read_velocity, readdx_upper, read_position_and_velocity, read_equi_position, readdx_binary_upper
from read import readdx_dynamical_matrix_gulp, read_binary_velocity, read_binary_position, read_binary_force
from read import read_binary_energy, read_binary_stress_tensor
from write import write_lata_atomic_density
# from cell import Cell
from Functions import throw_atom2cell
import Functions as fun
import time
import cluster_analysis as mycluster
import cProfile
from mpi4py import MPI

comm = MPI.COMM_WORLD
split = comm.Split(comm.Get_rank(), comm.Get_rank())

MAXDUMP = 50001


class allatoms:

    def __init__(self):
        self.atoms = {}  # store atoms information, use atom ID as key
        self.bonds = {}  # store bonds information, use atom ID as key
        self.atomnum = 0
        self.bondnum = 0
        self.maxID = -1

    # read single frame bonds information: type atom1-ID  atom2-ID
    # read bond info
    #  1  1  6660
    #  3  6655  8703
    def read_bonds(self, filename):
        with open(filename) as f:
            for index, line in enumerate(f):
                t = line.rstrip().split()
                if int(t[1]) not in self.bonds:
                    self.bonds[int(t[1])] = set([int(t[2])])
                else:
                    self.bonds[int(t[1])].add(int(t[2]))
                if int(t[2]) not in self.bonds:
                    self.bonds[int(t[2])] = set([int(t[1])])
                else:
                    self.bonds[int(t[2])].add(int(t[1]))
                self.bondnum += 1
        print("Import ", self.bondnum, " bonds!")
        # print(self.bonds)

    #  after delete atoms and add new atoms, the maximum ID of atoms might larger than self.atomnum
    #  because some ID has no atoms associated
    #  assuming there are only  Zn  H  C  N elements

    def outputxyz(self, xyzfile):
        Zn, H, C, N = {}, {}, {}, {}
        # nZn, nH, nC, nN = 1, 1, 1, 1
        j = 1
        # print(self.atoms[10].type_,"yes here")
        for i in self.atoms:
            if self.atoms[i].type_ == "Zn":
                Zn[j] = self.atoms[i].id_;
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "H":
                H[j] = self.atoms[i].id_;
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "C":
                C[j] = self.atoms[i].id_;
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "N":
                N[j] = self.atoms[i].id_;
                j += 1

        f = open(xyzfile, 'w')
        f.write('%d\n' % self.atomnum)
        f.write('add benzine to Zif4\n')
        for i in sorted(Zn.keys()):
            xyz = self.atoms[Zn[i]].xyz_
            # print(Zn[i], xyz)
            f.write("%d %s %f %f %f\n" % (i, self.atoms[Zn[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(H.keys()):
            xyz = self.atoms[H[i]].xyz_
            # print(xyz)
            f.write("%d %s %f %f %f\n" % (i, self.atoms[H[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(C.keys()):
            xyz = self.atoms[C[i]].xyz_
            f.write("%d %s %f %f %f\n" % (i, self.atoms[C[i]].type_, xyz[0], xyz[1], xyz[2]))
        for i in sorted(N.keys()):
            xyz = self.atoms[N[i]].xyz_
            f.write("%d %s %f %f %f\n" % (i, self.atoms[N[i]].type_, xyz[0], xyz[1], xyz[2]))
        f.close()

    def outputbond(self, bondfile):
        oldID2newID = {}
        newID2oldID = {}
        j = 1
        newbonds = {}
        bondtype = {"Zn-N": 1, "N-Zn": 1, "H-C": 2, "C-H": 2, "C-N": 3, "N-C": 3, "C-C": 4}
        for i in self.atoms:
            if self.atoms[i].type_ == "Zn":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "H":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "C":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        for i in self.atoms:
            if self.atoms[i].type_ == "N":
                oldID2newID[self.atoms[i].id_] = j
                newID2oldID[j] = self.atoms[i].id_
                j += 1
        # print("old2new ID:",oldID2newID)
        j = 1
        while j <= self.maxID:
            if j in self.bonds:
                newJ = oldID2newID[j]
                neighbour = self.bonds[j]
                for k in neighbour:
                    newK = oldID2newID[k]
                    assert newJ is not newK, "bonds error, same atoms!"
                    if newJ < newK:
                        if newJ in newbonds:
                            newbonds[newJ].add(newK)
                        else:
                            newbonds[newJ] = set([newK])
                    else:
                        if newK in newbonds:
                            newbonds[newK].add(newJ)
                        else:
                            newbonds[newK] = set([newJ])
            j += 1

        f = open(bondfile, 'w')
        # print(newbonds)
        for i in newbonds:
            for j in newbonds[i]:
                if i < j:
                    # print(newID2oldID[i], newID2oldID[j], self.atoms[newID2oldID[i]].type_, self.atoms[newID2oldID[j]].type_ )
                    typ = bondtype[self.atoms[newID2oldID[i]].type_ + "-" + self.atoms[newID2oldID[j]].type_]
                    f.write("%d %d %d\n" % (int(typ), i, j))
        f.close()

    # delete atoms and its associated bonds
    def delete_atom(self, atomID):
        assert atomID in self.atoms, "Atom %d does not exist!" % atomID
        neigh = self.bonds.pop(atomID)
        self.bondnum -= len(neigh)
        # print(self.bondnum)
        for i in neigh:
            if len(self.bonds[i]) == 1:
                self.bonds.pop(i)
            else:
                self.bonds[i].remove(atomID)
        self.atoms.pop(atomID)
        self.atomnum -= 1

    def add_atom(self, typ, coord):
        self.atomnum += 1
        newid = self.atomnum
        while newid in self.atoms:
            newid += 1
        self.atoms[newid] = Atom(newid, typ, coord[0], coord[1], coord[2])
        self.maxID = max(self.maxID, newid)
        return newid

    def add_bond(self, id1, id2):
        self.bondnum += 1
        assert id1 not in self.bonds or id2 not in self.bonds[id1], "Bond %d -> %d already exist!" % (id1, id2)

        if id1 not in self.bonds:
            self.bonds[id1] = set([id2])
        else:
            self.bonds[id1].add(id2)
        if id2 not in self.bonds:
            self.bonds[id2] = set([id1])
        else:
            self.bonds[id2].add(id1)

# input parameters
def parseOptions(comm):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("LATA", help="lata_file(typical format: ID type q x y z ...")
    parser.add_argument("out", help="output file prefix name")
    parser.add_argument("Cutoff", type=float, help="cutoff to determine if two points are close enough")
    parser.add_argument("Charge", type=int, help="if file contain charge [1/0]")
    parser.add_argument('-c', '--cellsize', type=float, help="cellsize to divide the simulation box")
    parser.add_argument('-ts', '--timestep', type=float, help="time_step(unit: femtosecond for fun-6 or picosecond for other) of each frame")
    parser.add_argument('-num', '--atom_num', type=int, help="number of atoms in the system for fun 4")
    parser.add_argument('-temp', '--temper', type=float, help="temperature of simulations for fun 6")
    parser.add_argument('-fnum', '--frame_num', type=int, help="number of frames in the system for fun 4 or 5 or 6")
    parser.add_argument('-cfnum', '--coor_frame_num', type=int, help="delta frames for correlation calculation")
    parser.add_argument('-sfnum', '--strong_frame_num', type=int, help="delta frames for strong correlation --correlation beyond this is ~zero")
    parser.add_argument('-mode_start', '--mode_start', type=int, help="for fun 6, first mode to calculate (first mode start from 1, not 0)")
    parser.add_argument('-mode_end', '--mode_end', type=int, help="for fun 6, last mode to calculate")
    parser.add_argument('-dym', '--dyn_matrix', type=str, help="dynamic matrix file generated by lammps () for fun 5")
    parser.add_argument('-Q0', '--Q0', type=str, help="Heat flux (symbol: J or Q) at each timestep for fun 6")
    parser.add_argument('-eng', '--energy', type=str, help="per-atom energy binary file generated by lammps () for fun 6")
    parser.add_argument('-S', '--stress_tensor', type=str, help="per-atom stress_tensor binary file generated by lammps () for fun 6")
    parser.add_argument('-equi_pos', '--equi_position', type=str, help="equilibrium position lata4olivia file by lammps for fun 5")
    parser.add_argument('-l', '--num_unit_cell_xyz', type=int, nargs=3, help="number of unit cell in X Y Z of the sample for fun 5")
    parser.add_argument('-fr', '--frame', type=int, help="use nth frame or start from nth frame")
    parser.add_argument('-f', '--fun', type=int, help="function: \n"
                                                      "1. number of nth neighbor (type Si) from a central atom (type Si)\n"
                                                      "2. density fluctuation using a spherical sample size(radius=R)\n"
                                                      "3. number of nth neighbor (type Si) from a central atom (type Si)"
                                                      "write down all coordination sequence, give histogram as well\n"
                                                      "4. get Spectral Energy Density from MD velocity\n"
                                                      "5. get Phonon lifetime from MD trajectory\n"
                                                      "6. GKMD analysis from normal mode and MD trajectory")
    '''
    parser.add_argument("NUM", type=int, help="Number of benzine to add")
    parser.add_argument("XYZ", help="XYZ_file(format: ID elementname X Y Z")
    parser.add_argument("BOND", help="Original bond list")
    parser.add_argument("OUTBOND", help="Output new bond list with benzine")
    parser.add_argument("OUTXYZ",  help="Output new xyz  file with benzine")
    parser.add_argument("Lx", type=float, help="Box length X")
    parser.add_argument("Ly", type=float, help="Box length Y")
    parser.add_argument("Lz", type=float, help="Box length Z")
    parser.add_argument("Nmin", type=int, help="minimum N index")
    parser.add_argument("Nmax", type=int, help="maximum N index")
    '''
    args = None
    try:
        if comm.Get_rank() == 0:
            args = parser.parse_args()
    finally:
        args = comm.bcast(args, root=0)
    if args is None:
        exit(0)
    return args

def main():
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = parseOptions(comm)

    if rank == 0:
        print("open files:", args.LATA)
    mysystem = [SimulationBox() for _ in range(MAXDUMP)]  # remember not to use [SimulationBox()] * 10 because this initiate 10 references to the same object
    print(mysystem[0].origin)

    print("Choose function: ", args.fun)
    print("Chose  frame:    ", args.frame)

    # function 1
    if args.fun == 1:
        if args.cellsize is None:
            print("Forget cellsize: -c XXX")
            sys.exit()
        (frame, nmols, starttime, endtime) = read_lata(mysystem, args.LATA, args.Charge, 500)
        print(frame, nmols, starttime, endtime)
        neighbor_list = None
        # only last frame is considered here
        frame = args.frame
        for i in range(frame, frame+1):
            fun.throw_atom2cell(mysystem, i, args.cellsize)
            neighbor_list = fun.get_connection_list(mysystem, i, 2, 1)
        # print(neighbor_list[0])
            graph = mycluster.Graph()
            graph.add_all_Edge(neighbor_list)
            neigh, dist, bond_to_current_neighbor_from_inner = graph.gen_nth_neighbor(mysystem, i)
            print("neighbor results:(index is nth order) ---->  ", end=" ")
            print(" ".join(str("{0:.3f}".format(x)) for x in neigh))
            # print(" ".join(str(x) for x in neigh))
            print("\nneighbor distance results:(index is nth order) ---->  ", end=" ")
            print(" ".join(str("{0:.3f}".format(x)) for x in dist))
            print("\nnum of bond to current level neighbor from inner level:(index is nth order) ---->  ", end=" ")
            print(" ".join(str("{0:.3f}".format(x)) for x in bond_to_current_neighbor_from_inner))

    # function 3
    # calculate CSQ sequence
    if args.fun == 3:
        if args.cellsize is None:
            print("Forget cellsize: -c XXX")
            sys.exit()
        (frame, nmols, starttime, endtime) = read_lata(mysystem, args.LATA, args.Charge, 500)
        print("total_frame, nmols, starttime, endtime")
        print(frame + 1, nmols, starttime, endtime)
        neighbor_list = None
        # only last frame is considered here
        frame = args.frame
        print("select [", frame, "] frame")
        # get Coordination sequence from all frames and also do average
        from collections import defaultdict
        CSQ = defaultdict(list)
        record_atom_density = defaultdict(list)
        for i in range(frame, frame+1):
            fun.throw_atom2cell(mysystem, i, args.cellsize)
            neighbor_list = fun.get_connection_list_Si_Si(mysystem, i, 2, 3.45)  # using 3.45Å as cutoff for Si-Si
            # print(neighbor_list[0])
            graph = mycluster.Graph()
            graph.add_all_Edge(neighbor_list)
            graph.gen_nth_neighbor_CSQ(mysystem, CSQ, i, record_atom_density)
            print(len(record_atom_density))
            write_lata_atomic_density(mysystem, "lata."+ args.out + ".N6_density." + str(frame), args.Charge, frame, record_atom_density)
        #if frame >= 0:
        #    for n in CSQ.items():
        #        CSQ[n][4] /= (frame + 1)
        f = open(args.out + ".Coordination_sequence_CSQ." + str(frame), 'w')
        f.write('#N1 is the first neighbor\n# CSQ sequence(N1,N2, ..., N6), N1[7], N2, N5, N6[10], Total_N1-N6(count_within_cluster)[11], '
                ' Total_N1-N2[12] Total_N3-N4[13] Frequency[14] cluster_density(N0+..+N6/2)/Vol[15]\n')
        for k, v in CSQ.items():
            # cluster density
            density =  v[4] / (4.0 / 3.0 * 3.1415926 * (v[8]/v[7]) ** 3)
            f.write('%1d %2d %2d %2d %2d %2d      %1d %2d   %2d %2d  %3.2f %3.2f %3.2f %3d  %1.5f\n' %
                    (k[0], k[1],k[2],k[3],k[4],k[5], v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], density))
        f.close()


    # function 2
    if args.fun == 2:
        if args.cellsize is None:
            print("Forget sphere radius: -c XXX")
            sys.exit()
        if args.frame is None:
            print("Forget frame number(starting from 1): -c XXX")
            sys.exit()
        (frame, nmols, starttime, endtime) = read_lata(mysystem, args.LATA, args.Charge, 500)
        # only last frame is considered here
        for i in range(frame, frame+1):
            fun.throw_atom2cell(mysystem, i, args.cellsize+0.000001)

    # function 4
    # get Spectral Energy Density (Unit: J s)
    if args.fun == 4:
        if args.atom_num is None:
            print("Forget atom number: -num XXX")
            sys.exit()
        (num_frame, nmols, starttime, endtime, velocity_array, mass_array) = read_velocity(args.LATA, args.Charge, args.atom_num, args.frame_num)
        res = fun.get_spectral_energy_density(velocity_array, mass_array, args.frame_num, args.timestep)
        f = open(args.out + ".spectral_energy_density", 'w')
        f.write('# omega/2Pi(THz) Phi(omega) k=0\n')
        for i in range(len(res)):
            f.write('%3.3f %3.5f\n' % (res[i][0], res[i][1]))
        f.close()

    # function 5
    # get Spectral Energy Density projection on each frequency (Unit: J s)
    if args.fun == 5:
        if args.atom_num is None:
            print("Forget atom number: -num XXX")
            sys.exit()

        # get normal mode
        # read dyn_matrix
        if args.dyn_matrix is None:
            print("Forget dynamic matrix file: -dyn FILENAME")
            sys.exit()


        import time
        from datetime import datetime
        now = datetime.now().time()
        print("now =", now)
        start_time = time.time()
        #dyn = readdx_upper(args.dyn_matrix, args.atom_num)
        #exit()
        #print(dyn)

        print(args.num_unit_cell_xyz)
        # read equilibrium position of each atom
        equi_position, unit_cell_size, cell_position = read_equi_position(args.equi_position, args.Charge, args.atom_num, args.num_unit_cell_xyz)

        # read velocity files
        (num_frame, nmols, velocity_array) = \
            read_binary_velocity(args.LATA, args.atom_num, args.frame_num, "SMALLBIG")
        print("velocity of atoms reading completed!!")

        # read dynamic matrix
        dyn = readdx_upper(args.dyn_matrix, args.atom_num)


        fun.get_phonon_life_time2(dyn, velocity_array, num_frame, nmols, args.timestep, args.out,
                                  args.num_unit_cell_xyz, unit_cell_size, cell_position)

        """
        # read position and velocity files
        print("reading the position and velocity of atoms in multiple trajectory ....")
        (num_frame, nmols, starttime, endtime, position_array, velocity_array) = \
            read_position_and_velocity(args.LATA, args.Charge, args.atom_num, args.frame_num, equi_position)
        print("Position and velocity of atoms reading completed!!")

        fun.get_phonon_life_time2(dyn, velocity_array, num_frame + 1, nmols, args.timestep, args.out,
                                 args.num_unit_cell_xyz, unit_cell_size, cell_position)
        """

        now = datetime.now().time()
        print("now =", now)
        print("Time used: ", time.time() - start_time, " seconds!")
        print("\nJob completed!!!\n")


    # function 6
    #  GKMA analysis
    if args.fun == 6:
        if args.atom_num is None:
            print("Forget atom number: -num XXX")
            sys.exit()

        if args.frame_num is None:
            print("Forget frame number: -fnum XXX")
            sys.exit()

        if args.coor_frame_num is None:
            print("Forget correlation frame number: -cfnum XXX")
            sys.exit()
        if args.strong_frame_num is None:
            print("Forget strong correlation frame number: -sfnum XXX")
            sys.exit()

        if args.temper is None:
            print("Forget temperature: -temp XXX[float]")
            sys.exit()

        if args.timestep is None:
            print("Forget timestep: -ts XXX[float](timestep per frame in femtosecond)")
            sys.exit()
        # get normal mode
        # read dyn_matrix
        if args.dyn_matrix is None:
            print("Forget dynamic matrix file: -dyn FILENAME")
            sys.exit()

        if args.energy is None:
            print("Forget per-atom energy file: -eng FILENAME")
            sys.exit()

        if args.stress_tensor is None:
            print("Forget per-atom stress_tensor file: -S FILENAME")
            sys.exit()

        if args.Q0 is None:
            print("Forget heat_flux_along_x/y/z at each timestep file: -Q0 FILENAME")
            sys.exit()

        import time
        from datetime import datetime
        now = datetime.now().time()
        print("now =", now)
        start_time = time.time()
        #dyn = readdx_upper(args.dyn_matrix, args.atom_num)
        #exit()
        #print(dyn)

        #print(args.num_unit_cell_xyz)
        # read equilibrium position of each atom
        # equi_position, unit_cell_size, cell_position = read_equi_position(args.equi_position, args.Charge, args.atom_num, args.num_unit_cell_xyz)

        #import pandas as pd
        # read Q_0 (heat flux at each step in x,y,z direction)
        Q_0 = np.genfromtxt(args.Q0, delimiter=None, dtype=float, usecols=(14,15,16), autostrip=True)
        #print(Q_0)
        #print(Q_0.shape)
        #print(Q_0[0,:])
        #exit()
        ##################################################
        # read dynamic matrix, regular or binary file
        dyn = readdx_upper(args.dyn_matrix, args.atom_num)
        #dyn = readdx_binary_upper(args.dyn_matrix, args.atom_num)

        # do diagonalization first
        # (1) numpy
        from numpy import linalg as LA
        print("Starting diagonalization of the dynamical matrix ....")
        start_time = time.time()
        w_value, e_vector = LA.eigh(dyn)
        print("w_value_shape, e_vector_shape", w_value.shape, e_vector.shape)
        print("Finishing diagonalization of the dynamical matrix :)")
        print("Time used: ", time.time() - start_time, " seconds!")
        # | # (2) scipy
        # | from scipy import linalg as spLA
        # | print("Starting diagonalization of the dynamical matrix ....")
        # | start_time = time.time()
        # | w_value, e_vector = spLA.eigh(dyn)
        # | print("w_value_shape, e_vector_shape", w_value.shape, e_vector.shape)
        # | print("Finishing diagonalization of the dynamical matrix :)")
        # | print("Time used: ", time.time() - start_time, " seconds!")
        # | exit()
        #############################################
        # read velocity*sqared_root_mass files
        (num_frame, nmols, velocity_ms_array, mass_sqr_array, volume) = \
            read_binary_velocity(args.LATA, args.atom_num, args.frame_num, "SMALLBIG")
        # read position files
        # (num_frame, nmols, position_array, mass_array) = \
        #     read_binary_position(args.LATA, args.atom_num, args.frame_num, "SMALLBIG")
        ## read force files
        #force_array= read_binary_force(args.LATA, args.atom_num, args.frame_num, "SMALLBIG")
        # read energy files
        energy_array = read_binary_energy(args.energy, args.atom_num, num_frame, "SMALLBIG")
        # read energy files
        stress_tensor_array = read_binary_stress_tensor(args.stress_tensor, args.atom_num, num_frame, "SMALLBIG")


        #print("velocity of atoms reading completed!!")
        small_flag = False
        if small_flag:
            dyn = np.array([[1,0,0],[0,1,0],[0,0,1]])
            #dyn = np.array([[1,2,3,0,1],[2,1,2,3,4,1],[3,2,2,0,2,1],[2,3,0,0,3,0],[0,4,2,3,1,0],[1,1,1,0,0,5]])
            mass_sqr_array = np.ones(3)
            velocity_ms_array = np.array([[1,-1,0],[1,0,1],[2,1,3]])
            energy_array = np.array([[1,0,2],[1,0,2],[1,0,2]])
            stress_tensor_array = np.ones(27).reshape(9,3)
            volume = 1.0
            args.temper = 1.0
            num_frame = 3
            nmols = 1
            args.timestep = 1
            args.out = "test"

        # group modes every 0.15 THz (X40 faster)
        fun.GKMA_fast(w_value, e_vector, mass_sqr_array, velocity_ms_array, energy_array, stress_tensor_array, Q_0, volume, args.temper,
                 num_frame, args.coor_frame_num, args.strong_frame_num, nmols, args.timestep, args.out)
        # calculate every single mode (slow)
        #fun.GKMA(w_value, e_vector, mass_sqr_array, velocity_ms_array, energy_array, stress_tensor_array, Q_0, volume, args.temper,
        #         num_frame, args.coor_frame_num, nmols, args.timestep, args.out, args.mode_start, args.mode_end)

        """
        # read position and velocity files
        print("reading the position and velocity of atoms in multiple trajectory ....")
        (num_frame, nmols, starttime, endtime, position_array, velocity_array) = \
            read_position_and_velocity(args.LATA, args.Charge, args.atom_num, args.frame_num, equi_position)
        print("Position and velocity of atoms reading completed!!")

        fun.get_phonon_life_time2(dyn, velocity_array, num_frame + 1, nmols, args.timestep, args.out,
                                 args.num_unit_cell_xyz, unit_cell_size, cell_position)
        """

        now = datetime.now().time()
        print("now =", now)
        print("Time used: ", time.time() - start_time, " seconds!")
    print("\nJob completed!!!!\n")

# yj remote
if __name__ == "__main__":
    main()
