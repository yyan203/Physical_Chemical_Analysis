
# read all frames from lata format: ID type q X Y Z ...
from Atom import Atom
from SimulationBox import SimulationBox

from ase.parallel import paropen
import struct
import time
import numpy as np


def read_lata(mysystem, filename, ifcharge, max_frame_to_read):

    # print(filename)
    f = open(filename)
    frame, nmols, starttime, endtime = -1, 0, 0, 0
    # print("ori", type(mysystem[0]), len(mysystem))

    while True:
        if frame + 1 >= max_frame_to_read:
            f.close()
            print("max frame:", frame, "reached")
            return frame, nmols, starttime, endtime
        line = f.readline()
        if not line:
            break
        t = line.split(' ')
        if t[0] == "ITEM:":
            frame += 1
            # print(t[0], t, frame)
            line = f.readline()
            t = line.split(' ')
            if frame == 0:
                starttime = int(t[0])
            endtime = int(t[0])
            line = f.readline()
            # number of atom
            line = f.readline()
            t = line.split(' ')
            nmols = int(t[0])
            # assign number of atoms
            mysystem[frame].nmols = nmols
            line = f.readline()

            # read box information
            line = f.readline()
            t = line.split(' ')
            mysystem[frame].origin[0] = float(t[0])
            # print("=", mysystem[0].origin[0], float(t[0]), frame)
            mysystem[frame].L[0] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            mysystem[frame].origin[1] = float(t[0])
            mysystem[frame].L[1] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            mysystem[frame].origin[2] = float(t[0])
            mysystem[frame].L[2] = float(t[1]) - float(t[0])

            line = f.readline()
            # read num atoms
            mysystem[frame].myatom = []
            for i in range(nmols):
                line = f.readline()
                t = line.split(' ')
                if ifcharge:
                    mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), float(t[3]), float(t[4]), float(t[5]), 0,0,0))
                else:
                    mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), float(t[2]), float(t[2]), float(t[4]), 0,0,0))
            # print("x=", mysystem[0].myatom[0].x)

    f.close()
    print("LATA file import:", nmols, " atoms!\n")
    return frame, nmols, starttime, endtime



# read velocity from lata format: ID type q X Y Z ...
#from Atom import Atom
from SimulationBox import SimulationBox
import numpy as np

# for fun 4 calculate spectral energy density
def read_velocity(filename, ifcharge, atom_num, frame_num):

    type_to_mass = {1: 15.999, 2: 28.0855}
    #type_to_mass = {1: 39.948}
    velocity_array = np.zeros([atom_num * 3, frame_num])
    mass_array = np.zeros(atom_num * 3)
    # print(filename)
    f = open(filename)
    frame, nmols, starttime, endtime = -1, 0, 0, 0
    # print("ori", type(mysystem[0]), len(mysystem))

    while True:
        if frame + 1 >= frame_num:
            f.close()
            print("max frame:", frame, "reached")
            # print(velocity_array[196 * 3 + 2])
            # print(mass_array[(196 * 3):(196 * 3 + 3)])
            #exit()
            return frame, nmols, starttime, endtime, velocity_array, np.transpose(mass_array)
            #return frame, nmols, starttime, endtime,
        line = f.readline()
        if not line:
            break
        t = line.split(' ')
        if t[0] == "ITEM:":
            frame += 1
            print("read ", frame, " frame")
            # print(t[0], t, frame)
            line = f.readline()
            t = line.split(' ')
            if frame == 0:
                starttime = int(t[0])
            endtime = int(t[0])
            line = f.readline()
            # number of atom
            line = f.readline()
            t = line.split(' ')
            nmols = int(t[0])
            # assign number of atoms
            #mysystem[frame].nmols = nmols
            line = f.readline()

            # read box information
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[0] = float(t[0])
            # print("=", mysystem[0].origin[0], float(t[0]), frame)
            #mysystem[frame].L[0] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[1] = float(t[0])
            #mysystem[frame].L[1] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[2] = float(t[0])
            #mysystem[frame].L[2] = float(t[1]) - float(t[0])

            line = f.readline()
            # read num atoms
            #mysystem[frame].myatom = []
            for i in range(nmols):
                line = f.readline()
                t = line.split(' ')
                if ifcharge:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[3]), float(t[4]), float(t[5])))
                    nth_atom = (int(t[0]) - 1) * 3
                    if frame == 0:
                        mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    velocity_array[nth_atom][frame] = float(t[3])
                    velocity_array[nth_atom + 1][frame] = float(t[4])
                    velocity_array[nth_atom + 2][frame] = float(t[5])
                else:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[2]), float(t[3]), float(t[4])))
                    nth_atom = (int(t[0]) - 1) * 3
                    if frame == 0:
                        mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    velocity_array[nth_atom][frame] = float(t[2])
                    velocity_array[nth_atom + 1][frame] = float(t[3])
                    velocity_array[nth_atom + 2][frame] = float(t[4])
            # print("x=", mysystem[0].myatom[0].x)
    f.close()
    print("LATA velocity file import:", nmols, " atoms!\n")
    return frame, nmols, starttime, endtime, velocity_array, np.transpose(mass_array)


#  fun 5, read dynamic matrix file
def readdx_upper(dynamic_matrix_filename, number_atom):

    import numpy as np
    desired_width = 640
    np.set_printoptions(linewidth=desired_width)

    #type_to_mass = {1: 15.999, 2: 28.0855}
    #mass_array = np.zeros(number_atom * 3)

    dyn = np.zeros([number_atom * 3, number_atom * 3], dtype=float)
    f = open(dynamic_matrix_filename, "r")
    if f is None:
        print("File %s does not exist!\n" % dynamic_matrix_filename);
        exit(0)
    #n2 = number_atom * number_atom * 9
    #n1 = number_atom * 3
    # for (int i=0;i< n2;i++){A[i]=0.0;}
    # double f1, f2, f3;
    # double avef1, avef2, avef3;
    # char templine[ONELINEMAX];
    # // ith atom  and jth atom
    # printf("\n\n>>>Please check if the dynamical matrix is symmetrical!\n"
    #                ">>>The following are [i(0,1,2)][i(0,1,2)] example terms "
    #                "along the diagonal:\n\n");
    print("Start reading the dynamic matrix")
    print(number_atom)
    #exit()
    #debug = True
    debug = False
    if debug is False:
        for i in range(number_atom):
            print("dyn_matrix", i, "/", number_atom)
            for k in range(3):
                for j in range(number_atom):
                    line = f.readline()
                    t = line.split(' ')
                    dyn[i*3+k][j*3+0] += float(t[0]) * 9648.5 # for unit conversion, so that w has unit of 1/ps or THz
                    dyn[i*3+k][j*3+1] += float(t[1]) * 9648.5
                    dyn[i*3+k][j*3+2] += float(t[2]) * 9648.5

        for i in range(3 * number_atom):
            print("dyn_matrix", i, "/", 3 * number_atom)
            for j in range(3 * number_atom):
                if i < j:
                    temp = 0.5 * (dyn[i][j] + dyn[j][i])
                    dyn[i][j] = dyn[j][i] = temp
    #  |  if debug:
    #  |      for i in range(3 * number_atom):
    #  |          print("dyn_matrix", i, "/", 3 * number_atom)
    #  |          for j in range(3 * number_atom):
    #  |              if i == j:
    #  |                  dyn[i][j] = 1.0
    #  |              else:
    #  |                  dyn[i][j] = dyn[j][i] = 0.0

    #print(dyn[0:12, 0:12])
    #from numpy import linalg as LA
    #w, e = LA.eigh(dyn)
    #print(w)
    f.close()
    print("Successfully imported the Dynamical matrix!\n")

    return dyn
# end of readdx

# read binary dynamic matrix file
def readdx_binary_upper(dynamic_matrix_filename, atom_num):
    dyn = np.empty((3,atom_num*atom_num*3))
    fileobj = paropen(dynamic_matrix_filename, "rb")
    #volume = 1.0
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)
######################
    start_time = time.time()
    while True:
        print("reading dynamic matrix:")
        try:
            for i in range(3):
                dyn[i,:] = read_variables("=" + str(atom_num*atom_num*3) + "d")
        except EOFError:
            break
    fileobj.close()
    #print("dynamic_matrix[0,0]:",dyn[0,0])
    #print(dyn.shape)
    dyn=dyn.reshape((atom_num*3, atom_num*3))
    dyn=0.5*(dyn + dyn.transpose())  # dyn symmetric
    if dyn[atom_num, atom_num-2] != dyn[atom_num-2, atom_num]:
        print(atom_num, dyn[atom_num, atom_num-2], dyn[atom_num-2, atom_num])
        print(dyn[0,0])
    assert dyn[atom_num, atom_num-2] == dyn[atom_num-2, atom_num]
    dyn = dyn * 9648.5 # for unit conversion, so that w has unit of 1/ps or THz
    #print(dyn.shape)
    #print(dyn[1000,1000])
    print("Finishing reading the dynamical matrix :)")
    print("Time used: ", time.time() - start_time, " seconds!")
    #exit()
    return dyn


# read velocity from lata format: ID type q X Y Z ...
#from Atom import Atom
from SimulationBox import SimulationBox
import numpy as np

from math import sqrt
# for fun 5 calculate phonon life time for each frequency
def read_position_and_velocity(filename, ifcharge, atom_num, frame_num, equi_position):

    type_to_mass = {1: sqrt(15.9994/atom_num), 2: sqrt(28.0855/atom_num)}
    ##type_to_mass = {1: 1, 2: 1}
    velocity_array = np.zeros([atom_num * 3, frame_num])
    position_array = np.zeros([atom_num * 3, frame_num])
    mass_array = np.zeros(atom_num * 3)
    # print(filename)
    f = open(filename)
    frame, nmols, starttime, endtime = -1, 0, 0, 0
    # print("ori", type(mysystem[0]), len(mysystem))

    while True:
        if frame + 1 >= frame_num:
            f.close()
            print("max frame:", frame, "reached")
            # print(velocity_array[196 * 3 + 2])
            # print(mass_array[(196 * 3):(196 * 3 + 3)])
            #exit()
            # print(position_array)
            # print(velocity_array)
            return frame, nmols, starttime, endtime, position_array, velocity_array
            #return frame, nmols, starttime, endtime,
        line = f.readline()
        if not line:
            break
        t = line.split(' ')
        if t[0] == "ITEM:":
            frame += 1
            print("read ", frame, " frame")
            # print(t[0], t, frame)
            line = f.readline()
            t = line.split(' ')
            if frame == 0:
                starttime = int(t[0])
            endtime = int(t[0])
            line = f.readline()
            # number of atom
            line = f.readline()
            t = line.split(' ')
            nmols = int(t[0])
            assert nmols == atom_num
            # assign number of atoms
            #mysystem[frame].nmols = nmols
            line = f.readline()

            # read box information
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[0] = float(t[0])
            # print("=", mysystem[0].origin[0], float(t[0]), frame)
            #mysystem[frame].L[0] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[1] = float(t[0])
            #mysystem[frame].L[1] = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            #mysystem[frame].origin[2] = float(t[0])
            #mysystem[frame].L[2] = float(t[1]) - float(t[0])

            line = f.readline()
            # read num atoms
            #mysystem[frame].myatom = []
            for i in range(nmols):
                line = f.readline()
                t = line.split(' ')
                if ifcharge:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[3]), float(t[4]), float(t[5])))
                    nth_atom = (int(t[0]) - 1) * 3
                    if frame == 0:
                        mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    position_array[nth_atom][frame]     = float(t[3]) - equi_position[nth_atom]
                    position_array[nth_atom + 1][frame] = float(t[4]) - equi_position[nth_atom + 1]
                    position_array[nth_atom + 2][frame] = float(t[5]) - equi_position[nth_atom + 2]
                    velocity_array[nth_atom][frame] = float(t[6])
                    velocity_array[nth_atom + 1][frame] = float(t[7])
                    velocity_array[nth_atom + 2][frame] = float(t[8])
                else:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[2]), float(t[3]), float(t[4])))
                    nth_atom = (int(t[0]) - 1) * 3
                    if frame == 0:
                        mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    position_array[nth_atom][frame]     = float(t[2]) - equi_position[nth_atom]
                    position_array[nth_atom + 1][frame] = float(t[3]) - equi_position[nth_atom + 1]
                    position_array[nth_atom + 2][frame] = float(t[4]) - equi_position[nth_atom + 2]
                    velocity_array[nth_atom][frame] = float(t[5])
                    velocity_array[nth_atom + 1][frame] = float(t[6])
                    velocity_array[nth_atom + 2][frame] = float(t[7])
            #  multiple by sqrt(m/N)
            position_array[:, frame] = np.multiply(position_array[:, frame], mass_array)
            velocity_array[:, frame] = np.multiply(velocity_array[:, frame], mass_array)
            # print("x=", mysystem[0].myatom[0].x)
    f.close()
    print("LATA position & velocity file import:", nmols, " atoms!\n")
    return frame, nmols, starttime, endtime, position_array, velocity_array



# for function 5, assuming no charge in lata4olivia file: id type vx vy vz
def read_binary_velocity(filename, atom_num, frame_num, intformat="SMALLBIG"):
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]
    colnames = ["id", "type", "vx", "vy", "vz"]
    fileobj = paropen(filename, "rb")
    volume = 1.0
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    #type_to_mass = {1: sqrt(15.999/atom_num), 2: sqrt(28.0855/atom_num)}
    #type_to_mass = {1: sqrt(39.948/atom_num)}
    #type_to_mass = {1: sqrt(15.9994/atom_num), 2: sqrt(28.0855/atom_num)}
    type_to_mass = {1: sqrt(15.9994), 2: sqrt(28.0855)}
    velocity_array = np.empty([atom_num * 3, frame_num], dtype=float)
    #velocity_array = np.empty([atom_num * 3, frame_num])
    mass_array = np.zeros(atom_num * 3)
    iframe = 0
    n_atoms = 0
    while True:
        print("reading vel frame: %d/%d" % (iframe, frame_num))
        try:
            # read header

            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")

            # lammps cells/boxes can have different boundary conditions on each
            # sides (makes mainly sense for different non-periodic conditions
            # (e.g. [f]ixed and [s]hrink for a irradiation simulation))
            # periodic case: b 0 = 'p'
            # non-peridic cases 1: 'f', 2 : 's', 3: 'm'
            #pbc = np.sum(np.array(boundary).reshape((3, 2)), axis=1) == 0
            #cell, celldisp = construct_cell(diagdisp, offdiag)

            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))
            ids = data[:, colnames.index("id")].astype(int)
            types = data[:, colnames.index("type")].astype(int)
            sort_order = np.argsort(ids)
            ids = ids[sort_order]
            types = types[sort_order]
            #print(ids[0:10])
            #exit()
            data = data[sort_order, :]
            velocity = data[:, [2, 3, 4]]
            #print("velocity shape:", velocity.shape)
            #print("velocity: ", velocity[:10, :])
            #print(np.array(velocity).reshape((1, -1)).shape)
            #print("velocity_array_iframe_shape:", velocity_array[:, iframe].shape)
            #print(np.array(velocity).reshape((1, -1))[0,:10])
            #exit()
            velocity_array[:, iframe] = np.array(velocity).reshape((1, -1))
            if iframe == 0:
                mass_array = np.array([type_to_mass[_] for _ in types])
                mass_array = np.repeat(mass_array, 3)
                (xlo,xhi,ylo,yhi,zlo,zhi) = diagdisp
                volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
                #print(volume,xlo,xhi)
            #print(velocity_array[-3:, iframe])

            velocity_array[:, iframe] = np.multiply(velocity_array[:, iframe], mass_array)

            iframe += 1
        except EOFError:
            break
    fileobj.close()
    print("binary velocity*squared_root_mass reading finished!")
    print("velocity shape:", velocity_array.shape)
    print(velocity_array[200, 2])
    print("mass shape:", mass_array.shape)
    print(mass_array[200])
    #exit()
    print(volume)
    return iframe, n_atoms, velocity_array, mass_array, volume

# for function 6, assuming no charge in lata4olivia file: id type x y z
def read_binary_position(filename, atom_num, frame_num, intformat="SMALLBIG"):
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]
    colnames = ["id", "type", "x", "y", "z"]
    fileobj = paropen(filename, "rb")
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    #type_to_mass = {1: sqrt(15.999/atom_num), 2: sqrt(28.0855/atom_num)}
    #type_to_mass = {1: sqrt(39.948/atom_num)}
    type_to_mass = {1: sqrt(15.9994), 2: sqrt(28.0855)}
    position_array = np.empty([atom_num * 3, frame_num], dtype=float)
    #velocity_array = np.empty([atom_num * 3, frame_num])
    mass_array = np.zeros(atom_num * 3)
    iframe = 0
    n_atoms = 0
    while True:
        print("reading frame: %d/%d" % (iframe, frame_num))
        try:
            # read header

            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")

            # lammps cells/boxes can have different boundary conditions on each
            # sides (makes mainly sense for different non-periodic conditions
            # (e.g. [f]ixed and [s]hrink for a irradiation simulation))
            # periodic case: b 0 = 'p'
            # non-peridic cases 1: 'f', 2 : 's', 3: 'm'
            #pbc = np.sum(np.array(boundary).reshape((3, 2)), axis=1) == 0
            #cell, celldisp = construct_cell(diagdisp, offdiag)

            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))
            ids = data[:, colnames.index("id")].astype(int)
            types = data[:, colnames.index("type")].astype(int)
            sort_order = np.argsort(ids)
            ids = ids[sort_order]
            types = types[sort_order]
            #print(ids[0:10])
            data = data[sort_order, :]
            position = data[:, [2, 3, 4]]
            #print(position.shape)
            #print(position[:10, :])
            position_array[:, iframe] = np.array(position).reshape((1, -1))
            if iframe == 0:
                mass_array = np.array([type_to_mass[_] for _ in types])
                mass_array = np.repeat(mass_array, 3)
            #print(position_array[-3:, iframe])
            #position_array[:, iframe] = np.multiply(position_array[:, iframe], mass_array)

            iframe += 1
            #print(ids[:10])
            #print(mass_array[:10])
            #exit()
        except EOFError:
            break
    fileobj.close()
    return iframe, n_atoms, position_array, mass_array

def read_binary_force(filename, atom_num, frame_num, intformat="SMALLBIG"):
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]
    colnames = ["id", "type", "fx", "fy", "fz"]
    fileobj = paropen(filename, "rb")
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    force_array = np.empty([atom_num * 3, frame_num], dtype=float)
    iframe = 0
    # n_atoms = 0
    while True:
        print("reading frame: %d/%d" % (iframe, frame_num))
        try:
            # read header

            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")
            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))
            ids = data[:, colnames.index("id")].astype(int)
            types = data[:, colnames.index("type")].astype(int)
            sort_order = np.argsort(ids)
            # ids = ids[sort_order]
            # types = types[sort_order]
            ## print(ids[0:10])
            data = data[sort_order, :]
            force = data[:, [2, 3, 4]]
            #print(force.shape)
            force_array[:, iframe] = np.array(force).reshape((1, -1))
            iframe += 1
            # exit()
        except EOFError:
            break
    fileobj.close()
    return force_array

def read_binary_energy(filename, atom_num, frame_num, intformat="SMALLBIG"):
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]
    colnames = ["id", "E"]
    fileobj = paropen(filename, "rb")
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    energy_array = np.empty([atom_num * 3, frame_num], dtype=float)
    iframe = 0
    n_atoms = 0
    while True:
        print("reading energy frame: %d/%d" % (iframe, frame_num))
        try:
            # read header

            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")
            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))
            ids = data[:, colnames.index("id")].astype(int)
            #types = data[:, colnames.index("type")].astype(int)
            sort_order = np.argsort(ids)
            data = data[sort_order, :]
            #print(data[:,[1]].shape)
            energy = np.repeat(data[:, [1]], 3)
            #print(energy.shape)
            #print(energy)
            energy_array[:, iframe] = np.array(energy).reshape((1, -1))
            iframe += 1
            #exit()
        except EOFError:
            break
    fileobj.close()
    print("binary energy reading finished!")
    print("energy shape:", energy_array.shape)
    print(energy_array[0,2])
    #exit()
    return energy_array

# read atomic stress tensor (3*3) defined by W_ab in https://docs.lammps.org/compute_stress_atom.html
# note W_ab is symmetric so that only six is needed
def read_binary_stress_tensor(filename, atom_num, frame_num, intformat="SMALLBIG"):
    tagformat, bigformat = dict(
        SMALLSMALL=("i", "i"), SMALLBIG=("i", "q"), BIGBIG=("q", "q")
    )[intformat]
    colnames = ["id", "xx", "yy", "zz", "xy", "xz", "yz"]
    fileobj = paropen(filename, "rb")
    def read_variables(string):
        obj_len = struct.calcsize(string)
        data_obj = fileobj.read(obj_len)
        if obj_len != len(data_obj):
            raise EOFError
        return struct.unpack(string, data_obj)

    stress_tensor_array = np.empty([atom_num * 9, frame_num], dtype=float)
    iframe = 0
    n_atoms = 0
    while True:
        print("reading stress frame: %d/%d" % (iframe, frame_num))
        try:
            # read header

            ntimestep, = read_variables("=" + bigformat)
            n_atoms, triclinic = read_variables("=" + bigformat + "i")
            boundary = read_variables("=6i")
            diagdisp = read_variables("=6d")
            if triclinic != 0:
                offdiag = read_variables("=3d")
            else:
                offdiag = (0.0,) * 3
            size_one, nchunk = read_variables("=2i")
            if len(colnames) != size_one:
                raise ValueError("Provided columns do not match binary file")
            data = []
            for _ in range(nchunk):
                # number-of-data-entries
                n_data, = read_variables("=i")
                # retrieve per atom data
                data += read_variables("=" + str(n_data) + "d")
            data = np.array(data).reshape((-1, size_one))
            ids = data[:, colnames.index("id")].astype(int)
            #types = data[:, colnames.index("type")].astype(int)
            sort_order = np.argsort(ids)
            data = data[sort_order, :]
            stress_tensor = data[:, [1, 4, 5, 4, 2, 6, 5, 6, 3]]
            #print(stress_tensor.shape)
            stress_tensor_array[:, iframe] = np.array(stress_tensor).reshape((1, -1))
            iframe += 1
            #exit()
        except EOFError:
            break
    fileobj.close()
    print("binary stress_tensor reading finished!")
    print("stress_tensor shape:", stress_tensor_array.shape)
    #print(stress_tensor_array[1000,2])
    #exit()
    return stress_tensor_array


# return equilibrium position of each atom
def read_equi_position(filename, ifcharge, atom_num, num_unit_cell):

    #type_to_mass = {1: sqrt(15.999/atom_num), 2: sqrt(28.0855/atom_num)}
    #type_to_mass = {1: 1, 2: 1}
    equi_position_array = np.zeros([atom_num * 3])
    cell_position_array = np.zeros([atom_num * 3], dtype=complex)
    #mass_array = np.zeros(atom_num * 3)
    frame = 0
    # print(filename)
    f = open(filename)
    while True:
        line = f.readline()
        if not line:
            break
        t = line.split(' ')
        if t[0] == "ITEM:":
            # print(t[0], t, frame)
            line = f.readline()
            t = line.split(' ')
            line = f.readline()
            # number of atom
            line = f.readline()
            t = line.split(' ')
            nmols = int(t[0])
            assert nmols == atom_num
            # assign number of atoms
            line = f.readline()
            # read box information
            line = f.readline()
            t = line.split(' ')
            xlo = float(t[0])
            lx = float(t[1]) - float(t[0])
            line = f.readline()
            t = line.split(' ')
            ylo = float(t[0])
            ly = float(t[1]) - float(t[0])
            line = f.readline()
            zlo = float(t[0])
            t = line.split(' ')
            lz = float(t[1]) - float(t[0])
            unit_cell_size = (lx/num_unit_cell[0], ly/num_unit_cell[1], lz/num_unit_cell[2])

            line = f.readline()
            # read num atoms
            #mysystem[frame].myatom = []
            for i in range(nmols):
                line = f.readline()
                t = line.split(' ')
                if ifcharge:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[3]), float(t[4]), float(t[5])))
                    nth_atom = (int(t[0]) - 1) * 3
                    #mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    equi_position_array[nth_atom:nth_atom + 3] = float(t[3]), float(t[4]), float(t[5])
                    cell_position_array[nth_atom:nth_atom + 3] = 1j * (xlo + int((float(t[3]) - xlo) / unit_cell_size[0]) * unit_cell_size[0]),\
                                                                 1j * (ylo + int((float(t[4]) - ylo) / unit_cell_size[1]) * unit_cell_size[1]),\
                                                                 1j * (zlo + int((float(t[5]) - zlo) / unit_cell_size[2]) * unit_cell_size[2])
                    #if nth_atom == 0:
                    #    print(xlo, float(t[3]), unit_cell_size[0], cell_position_array[0])
                    #    exit()
                else:
                    #mysystem[frame].myatom.append(Atom(int(t[0]), int(t[1]), 0.0,0.0,0.0, float(t[2]), float(t[3]), float(t[4])))
                    nth_atom = (int(t[0]) - 1) * 3
                    #mass_array[nth_atom:nth_atom + 3] = [type_to_mass[int(t[1])], type_to_mass[int(t[1])], type_to_mass[int(t[1])]]
                    equi_position_array[nth_atom:nth_atom + 3] = float(t[2]), float(t[3]), float(t[4])
                    cell_position_array[nth_atom:nth_atom + 3] = 1j * (xlo + int((float(t[2]) - xlo) / unit_cell_size[0]) * unit_cell_size[0]), \
                                                                 1j * (ylo + int((float(t[3]) - ylo) / unit_cell_size[1]) * unit_cell_size[1]), \
                                                                 1j * (zlo + int((float(t[4]) - zlo) / unit_cell_size[2]) * unit_cell_size[2])
            # print("x=", mysystem[0].myatom[0].x)
    f.close()
    print("LATA equilibrium position  import:", nmols, " atoms!\n")
    #print(equi_position_array[0:10])
    #print(len(equi_position_array))
    #exit()
    return equi_position_array, unit_cell_size, cell_position_array

#  read dynamical matrix from GULP output
def readdx_dynamical_matrix_gulp(dynamic_matrix_filename, number_atom):

    import numpy as np
    desired_width = 640
    np.set_printoptions(linewidth=desired_width)

    #type_to_mass = {1: 15.999, 2: 28.0855}
    #mass_array = np.zeros(number_atom * 3)
    #return dyn

    dyn = np.zeros([number_atom * 3, number_atom * 3], dtype=complex)
    f = open(dynamic_matrix_filename, "r")
    if f is None:
        print("File %s does not exist!\n" % dynamic_matrix_filename);
        exit(0)
    print("Start reading the dynamic matrix")
    print(number_atom)
    #exit()
    ii = 0
    line = f.readline()
    while line:
        t = line.split()
        for i in range(len(t)):
            row, col = ii // (number_atom * 3), ii % (number_atom * 3)
            dyn[row, col] = float(t[i]) * 9648.5
            ii += 1
        if ii % 100 == 0:
            print(ii)
        line = f.readline()
    print("compare here")
    #print(dyn[0:24, 0:24])
    #print(dyn[1000:1005, 1:6])
    #from numpy import linalg as LA
    #w, e = LA.eigh(dyn)
    #print(w)
    #exit()
    f.close()
    print("Successfully imported the Dynamical matrix!\n")
    return dyn
