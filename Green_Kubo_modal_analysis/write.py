
# read all frames from lata format: ID type q X Y Z ...
from Atom import Atom
from SimulationBox import SimulationBox

def write_lata_atomic_density(mysystem, filename, ifcharge, frame, record_atom_density):

    # print(filename)
    f = open(filename, 'w')

    nmols = mysystem[frame].nmols
    lx = mysystem[frame].origin[0]
    hx = lx + mysystem[frame].L[0]
    ly = mysystem[frame].origin[1]
    hy = ly + mysystem[frame].L[2]
    lz = mysystem[frame].origin[2]
    hz = lz + mysystem[frame].L[2]
    f.write('ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp\n%f %f\n%f %f\n%f %f\n'
            'ITEM: ATOMS id type q x y z CN_Si-Si N6_density N6_ave_distance  N3_N4 N1_N2\n' %
            (frame, nmols, lx, hx, ly, hx, lz, hz))
    for i in range(nmols):
        if i in record_atom_density:
            f.write('%d %d %d %f %f %f %f %f %f %f %d\n' % (i, mysystem[frame].myatom[i].type, 0, mysystem[frame].myatom[i].x[0],
                                               mysystem[frame].myatom[i].x[1], mysystem[frame].myatom[i].x[2],
                                               record_atom_density[i][0], record_atom_density[i][1],
                                               record_atom_density[i][2], record_atom_density[i][3],
                                               record_atom_density[i][4]))
        else:
            f.write('%d %d %d %f %f %f %f %f %f %f %d\n' % (i, mysystem[frame].myatom[i].type, 0, mysystem[frame].myatom[i].x[0],
                                               mysystem[frame].myatom[i].x[1], mysystem[frame].myatom[i].x[2],
                                               0, 0, 0, 0, 0))

    f.close()
    print("LATA file with N6 density written!\n")
    return
