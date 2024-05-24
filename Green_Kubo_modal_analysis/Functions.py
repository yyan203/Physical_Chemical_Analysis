from math import fmod
from basic_function import bondlen
from random import random
from math import sqrt
import time

# 2019-05-01  YJ
# throw all atoms into the cells of each frame


def throw_atom2cell(mysystem, iframe, cellsize):

    mysystem[iframe].SetupCell(cellsize)

    originx= mysystem[iframe].origin[0]
    originy= mysystem[iframe].origin[1]
    originz= mysystem[iframe].origin[2]
    print("originx=", originx)
    cellsizex=(float)(mysystem[iframe].L[0]/mysystem[iframe].ncell[0])
    print("cellsizex=", cellsizex)
    cellsizey=(float)(mysystem[iframe].L[1]/mysystem[iframe].ncell[1])
    cellsizez=(float)(mysystem[iframe].L[2]/mysystem[iframe].ncell[2])

    # loop of every atoms of one frame in dump file
    myatom = mysystem[iframe].myatom
    # print(type(myatom), len(mysystem[iframe].myatom), mysystem[iframe].nmols)
    for i in range(mysystem[iframe].nmols):
        x1 = myatom[i].x[0]
        x2 = myatom[i].x[1]
        x3 = myatom[i].x[2]
        # if i == 0: print("myatom:", myatom[i].x)
        # according to the atom's position, allocate them to the right cell. Since cellsizex is float number,
        # to make it safe, when the index from atom's position is calculated from cellsizex, bring it back to within the cell.
        p1=(int)((x1-originx)/cellsizex)
        if p1==mysystem[iframe].ncell[0]: p1 = p1 - 1
        p2=(int)((x2-originy)/cellsizey)
        if p2==mysystem[iframe].ncell[1]: p2 = p2 - 1
        p3=(int)((x3-originz)/cellsizez)
        if p3==mysystem[iframe].ncell[2]: p3 = p3 - 1
        memberid = mysystem[iframe].mycell[p1][p2][p3].nmember
        # print ("p1=%d p2=%d p3=%d\n",p1,p2,p3)
        # allocate the atom's id to the right cell
        # mysystem[iframe]->mycell[p1][p2][p3]->member[memberid]=myatom[i]->id
        mysystem[iframe].mycell[p1][p2][p3].member[memberid] = i
        # print("mycell%d%d%d->member[%d]:%d\n",p1,p2,p3,memberid,mysystem[iframe]->mycell[p1][p2][p3]->member[memberid])

        # if (p1==0 && p2==0 && p3==0){print("mycell-000 has atom(id): %d\n",mysystem[iframe]->mycell[p1][p2][p3]->member[mysystem[iframe]->mycell[p1][p2][p3]->nmember])}
        # if(i==400){print("nmember of mycell[%d][%d][%d] is %d\n",p1,p2,p3,mysystem[iframe]->mycell[p1][p2][p3]->nmember)}

        # increase the nmember of cell when a new atom is added into it.
        mysystem[iframe].mycell[p1][p2][p3].nmember += 1
    print("\nAtoms thrown to cells finished!!\n")

# get number density fluctuation by sampling N random points and count atoms belonging to sphere of radius=R locating at these random points

def get_density_fluctuation(mysystem, iframe, sphereR, sampleN):

    count = [ 0.0 for _ in range(sampleN)]
    originx= mysystem[iframe].origin[0]
    originy= mysystem[iframe].origin[1]
    originz= mysystem[iframe].origin[2]
    Lx= mysystem[iframe].origin[0]
    Ly= mysystem[iframe].origin[1]
    Lz= mysystem[iframe].origin[2]

    for i in range(sampleN):
        x1 = random()*Lx + originx
        x2 = random()*Ly + originy
        x3 = random()*Lz + originz
        p1=(int)((x1-originx)/Lx)
        if p1==mysystem[iframe].ncell[0]: p1 = p1 - 1
        p2=(int)((x2-originy)/Ly)
        if p2==mysystem[iframe].ncell[1]: p2 = p2 - 1
        p3=(int)((x3-originz)/Lz)
        if p3==mysystem[iframe].ncell[2]: p3 = p3 - 1
        memberid = mysystem[iframe].mycell[p1][p2][p3].nmember





# 2019-06-06  YJ
# nearest neighbor (type Si) list of all  Si atoms that are connected with each other via bridging Oxygen, i.e. -Si-O-Si-

# cutoff Si-Si is the distance for farest Si neighbor to consider, here we can set a value > 3.45 (~first peak of Si-Si in silica)
# cutoff Si-O is the distance for Si-O bond, typically 2.2 A

from collections import defaultdict

def get_connection_list(mysystem, iframe, Si_type, O_type, cutoff_Si_O = 2.2, cutoff_Si_Si = 5.0):
    mysystem_ = mysystem[iframe]
    myatom = mysystem[iframe].myatom
    res = defaultdict(list)

    lx, ly, lz = mysystem_.L[0], mysystem_.L[1], mysystem_.L[2]
    ncell = []
    # initialize array ncell[3]
    for iii in range(3):
        ncell.append(mysystem_.ncell[iii])

    neigh_O = set() # O atom
    neigh_S = set() # Si atom
    ################## do analysis
    # loop among different cell
    for cellx in range(ncell[0]):
        for celly in range(ncell[1]):
            for cellz in range(ncell[2]):
                # //printf("mysystem[%d]->mycell[%d][%d][%d]->nmember:%d\n",iframe,cellx,celly,cellz,mysystem[iframe]->mycell[cellx][celly][cellz]->nmember)
                nmember = mysystem_.mycell[cellx][celly][cellz].nmember
                # //loop among all the atoms contained in one cell
                for imember in range(nmember):
                    comparei = mysystem_.mycell[cellx][celly][cellz].member[imember]
                    itype = myatom[comparei].type # get type
                    if itype == Si_type:
                        neigh_O.clear()
                        neigh_S.clear()

                        # loop among all neighbour cells
                        for ii in range(27):
                            c = mysystem_.mycell[cellx][celly][cellz].nbr[ii]
                            ncellx = (int) (c / (ncell[1] * ncell[2]))
                            ncelly = (int) (fmod((int) (c / ncell[2]), ncell[1]))
                            ncellz = (int) (fmod(c, ncell[2]))
                            #printf("mysystem[%d]->mycell[%d][%d][%d]: ncellx:%d ncelly:%d ncellz:%d\n",iframe,cellx,celly,cellz,ncellx,ncelly,ncellz)
                            inmember = mysystem_.mycell[ncellx][ncelly][ncellz].nmember
                            # loop among all the neighbour's atoms
                            for j in range(inmember):
                                comparej = mysystem_.mycell[ncellx][ncelly][ncellz].member[j]
                                jtype = myatom[comparej].type  # get jtype
                                if comparei != comparej:
                                    dist = bondlen(myatom[comparei], myatom[comparej], lx, ly, lz)
                                    if jtype == Si_type and dist <= cutoff_Si_Si:
                                        neigh_S.add(comparej)
                                    if jtype == O_type and dist <= cutoff_Si_O:
                                        neigh_O.add(comparej)
                        res[comparei] = []
                        for it in neigh_S:
                            for jt in neigh_O:
                                if bondlen(myatom[it], myatom[jt], lx, ly, lz) <= cutoff_Si_O:
                                    assert myatom[it].type == Si_type
                                    res[comparei].append(it)
                                    break
                        # debug
                        # if comparei == 0:
                        #     print(neigh_O, neigh_S)
                        # if myatom[comparei].id == 1544:
                        #     for i in res[comparei]:
                        #         print(myatom[i].id)
    return res

# 2020-04-23  YJ
# nearest neighbor (type Si) list of all  Si atoms that are close with each other via a cutoff = 3.45
# purely based on Si-Si distance
# cutoff Si-Si is the distance for farest Si neighbor to consider, here we choose a value = 3.45 (~first peak of Si-Si in silica)

def get_connection_list_Si_Si(mysystem, iframe, Si_type, cutoff_Si_Si = 3.45):
    mysystem_ = mysystem[iframe]
    myatom = mysystem[iframe].myatom
    res = defaultdict(list)

    lx, ly, lz = mysystem_.L[0], mysystem_.L[1], mysystem_.L[2]
    ncell = []
    # initialize array ncell[3]
    for iii in range(3):
        ncell.append(mysystem_.ncell[iii])

    #neigh_O = set() # O atom
    neigh_S = set() # Si atom
    ################## do analysis
    # loop among different cell
    for cellx in range(ncell[0]):
        for celly in range(ncell[1]):
            for cellz in range(ncell[2]):
                # //printf("mysystem[%d]->mycell[%d][%d][%d]->nmember:%d\n",iframe,cellx,celly,cellz,mysystem[iframe]->mycell[cellx][celly][cellz]->nmember)
                nmember = mysystem_.mycell[cellx][celly][cellz].nmember
                # //loop among all the atoms contained in one cell
                for imember in range(nmember):
                    comparei = mysystem_.mycell[cellx][celly][cellz].member[imember]
                    itype = myatom[comparei].type # get type
                    if itype == Si_type:
                        # neigh_O.clear()
                        neigh_S.clear()

                        # loop among all neighbour cells
                        for ii in range(27):
                            c = mysystem_.mycell[cellx][celly][cellz].nbr[ii]
                            ncellx = (int) (c / (ncell[1] * ncell[2]))
                            ncelly = (int) (fmod((int) (c / ncell[2]), ncell[1]))
                            ncellz = (int) (fmod(c, ncell[2]))
                            #printf("mysystem[%d]->mycell[%d][%d][%d]: ncellx:%d ncelly:%d ncellz:%d\n",iframe,cellx,celly,cellz,ncellx,ncelly,ncellz)
                            inmember = mysystem_.mycell[ncellx][ncelly][ncellz].nmember
                            # loop among all the neighbour's atoms
                            for j in range(inmember):
                                comparej = mysystem_.mycell[ncellx][ncelly][ncellz].member[j]
                                jtype = myatom[comparej].type  # get jtype
                                if comparei != comparej:
                                    dist = bondlen(myatom[comparei], myatom[comparej], lx, ly, lz)
                                    if jtype == Si_type and dist <= cutoff_Si_Si:
                                        neigh_S.add(comparej)
                                    # if jtype == O_type and dist <= cutoff_Si_O:
                                    # neigh_O.add(comparej)
                        res[comparei] = []
                        if len(neigh_S) <= 1:
                            print(comparei, myatom[comparei].id)
                        for it in neigh_S:
                            if bondlen(myatom[it], myatom[comparei], lx, ly, lz) <= cutoff_Si_Si:
                                assert myatom[it].type == Si_type
                                res[comparei].append(it)

    return res



# 2020-09-28  YJ
# calculate spectral energy density from MD velocity of atoms
# refer to Eq. (4) in Thomas, McGaughey, PHYSICAL REVIEW B 81, 081411(R) 2010

from math import sin, cos
import numpy as np

def get_spectral_energy_density(velocity_array, mass_array, frame_num, time_step):
    time_array = np.zeros(frame_num)
    for i in range(frame_num):
        time_array[i] = i * time_step  # in picosecond
    time_array = np.transpose(time_array)
    #k = 0.0
    miu = 40.0  # Thz
    pi = 3.1415926
    N = 500  # number of miu
    N_T = 1.0  # number of unit cell # default is one
    tau_0 = time_step * frame_num
    prefactor = 1.0 / 4.0 / pi / N_T / tau_0
    print("dt = ", time_step)

    res = [[0.0, 0.0] for _ in range(N)]
    for j in range(N):
        print("Calculate progress：", j, "/", N)
        res[j][0] = miu / N * j  # THz
        omega = miu / N * j * 2 * pi    # THz * 2 * pi
        real_part = np.dot(velocity_array,  np.cos(time_array * omega))
        imag_part = np.dot(velocity_array,  np.sin(time_array * (-1) * omega))
        sqr_sum = (real_part ** 2 + imag_part ** 2) * (time_step ** 2)
        #  print("sqr_sum & time_array & mass_array")
        #  print(len(sqr_sum), len(time_array), len(mass_array))
        #  print(sqr_sum)
        #  print(time_array)
        #  print(mass_array)
        #  exit()
        matrix = np.dot(np.transpose(sqr_sum), mass_array)
        res[j][1] = np.sum(matrix) * prefactor

        #if j == 10:
        #    print("10th atom integral")
        #    print(res_real_x[j], res_imag_x[j],)
        #    print(res_real_y[j], res_imag_y[j],)
        #    print(res_real_z[j], res_imag_z[j],)
        #if j == 243:
        #    print("200th atom integral")
        #    print(res_real_x[200], res_imag_x[200],)
        #    print(res_real_y[200], res_imag_y[200],)
        #    print(res_real_z[200], res_imag_z[200],)
        print(res[j][0], res[j][1])
    return res

# for function 5 to get phonon life time from autocorrelation function of normal mode energy
def get_phonon_life_time(dyn, position_array, velocity_array, frame_num, atom_num, timestep, outputfilename,
                         num_unit_cell, unit_cell_size, cell_position):
    #print(position_array[0, 0:9])
    #print(velocity_array[0, 0:9])
    #exit()
    #|  for e in range(atom_num * 3):
    #|      if e == 0:
    #|          print("the mean position is: ", np.mean(position_array[e]))
    #|          print(position_array[0, 0:9])
    #|      position_array[e] = np.subtract(position_array[e], np.mean(position_array[e]))
    #|      if e == 0:
    #|          print(position_array[0, 0:9])
    from numpy import linalg as LA
    import time
    pi = 3.1415926
    # eigenvalue, eigenvector from diagonalization of dynamical matrix
    print("Starting diagonalization of the dynamical matrix ....")
    start_time = time.time()
    w_value, e_vector = LA.eigh(dyn)
    print("Finishing diagonalization of the dynamical matrix :)")
    print("Time used: ", time.time() - start_time, " seconds!")
    idx = np.argsort(-w_value)
    #print(w_value[0:10])
    w_value = np.array([w_value[idx]])
    print(w_value[0, 0:10])
    #print(e_vector[0:10, 0:10])
    e_vector = e_vector[:, idx]
    e_vector = np.conj(e_vector)
    #print(e_vector[0:10, 0:10])
    print(e_vector.shape)
    print(position_array.shape)
    print(velocity_array.shape)

    assert len(e_vector) == atom_num * 3
    #print(q[0:2,0:2])
    #print("something is here")
    #temp = np.transpose(w_value)
    #print(temp)
    #print(temp.shape)
    count = np.arange(frame_num, 0, -1)
    #print(count)
    num_unit_cell = [1, 1, 1]
    number_of_k = num_unit_cell[0] * num_unit_cell[1] * num_unit_cell[2]
    phonon_life_time = np.zeros((number_of_k, atom_num * 3), dtype=float)
    #phonon_life_time_subtract_mean = np.zeros(atom_num * 3, dtype=float)
    #flag = 1
    print(w_value[0, 0:10])
    print(w_value.shape)
    #print(w_value[0, 0:2])
    for ll in range(0, num_unit_cell[0]):
        for m in range(0, num_unit_cell[1]):
            for n in range(0, num_unit_cell[2]):
                wave_vector_index = ll + num_unit_cell[0] * m + num_unit_cell[0] * num_unit_cell[1] * n

                #l1, m1, n1 = ll / num_unit_cell[0], m / num_unit_cell[1], n / num_unit_cell[2]
                k = np.tile([[2 * pi / num_unit_cell[0] / unit_cell_size[0] * ll,
                              2 * pi / num_unit_cell[1] / unit_cell_size[1] * m,
                              2 * pi / num_unit_cell[2] / unit_cell_size[2] * n]], atom_num)
                print("wave_vector_index:", wave_vector_index)
                print(unit_cell_size, ll, m, n)
                print(cell_position[0:10])
                print(k[0, 0:10])
                exp_term = np.transpose(np.exp(np.multiply(k, cell_position)))
                print(np.transpose(exp_term)[0, 0:10])
                #time.sleep(10)
                position_array_temp = np.multiply(position_array, exp_term)
                velocity_array_temp = np.multiply(velocity_array, exp_term)
                q = np.dot(np.transpose(e_vector), position_array_temp)
                q_d = np.dot(np.transpose(e_vector), velocity_array_temp)

                # get potential energy for each mode
                print("get potential energy for each mode", "  wave vector index: ", ll, "-", m, "-", n)
                q = np.multiply(np.conj(q), q)
                q = np.multiply(q, np.transpose(w_value)).real
                    #print(q[0:2,0:2])
                q = np.multiply(q, 0.5)
                #print(q[0:2,0:2])
                #exit()

                # get kinetic energy for each mode
                print("get kinetic energy for each mode")
                q_d = np.multiply(np.conj(q_d), q_d).real
                q_d = np.multiply(q_d, 0.5)
                print("###################################\n################################\n\n\n\n")
                # total energy for each mode
                E = np.add(q, q_d)
                if ll == 0 and m == 0 and n == 0:
                    print(q[0, 0:30])  # potential
                    print(q_d[0, 0:30]) # kinetic
                    print(E[0, 0:30])  # total E
                    f = open("Energy_first_mode.dat.243", 'w')
                    f.write('#time(ps)   displacement_first_atom    total_energy_first_mode  Potential_E  Kinetic_E'
                            'frequency:%3.5f\n' % (sqrt(w_value[0, 243])/2.0/pi))
                    for i in range(len(E[243])):
                        f.write('%3.3f %3.5f %3.5f %3.5f %3.5f\n' %
                                (i*timestep, position_array[0, i]/sqrt(39.948/atom_num), E[243, i], q[243,i], q_d[243,i]))
                    f.close()
                    #print("output first mode energy done")
                    #exit()
                    #correlation_length = frame_num * 2 + 1


                for e in range(atom_num * 3):
                    if e == 243 and ll == m == n == 0:
                        autocorrelation = np.correlate(E[e], E[e], mode="full")
                        autocorrelation = autocorrelation[autocorrelation.size//2:]
                        autocorrelation = np.divide(autocorrelation, count)
                        normlize_factor = autocorrelation[0]
                        autocorrelation = np.divide(autocorrelation, normlize_factor)
                        f = open("Energy_autocorrelation.before_mean_subtract.dat.243", 'w')
                        f.write("#frame  time(ps)  Energy_autocorrelation\n")
                        for i in range(len(autocorrelation)):
                            f.write('%d %3.5f %3.5f\n' % (i, i*timestep, autocorrelation[i]))
                        f.close()

                    E[e] = np.subtract(E[e], np.mean(E[e]))
                    autocorrelation = np.correlate(E[e], E[e], mode="full")
                    #print(E[e])
                    #print("autocorrelation size:", autocorrelation.size)
                    autocorrelation = autocorrelation[autocorrelation.size//2:]
                    autocorrelation = np.divide(autocorrelation, count)
                    normlize_factor = autocorrelation[0]
                    autocorrelation = np.divide(autocorrelation, normlize_factor)
                    if e == 243 and ll == m == n == 0:
                        f = open("Energy_autocorrelation.dat.243", 'w')
                        f.write("#frame  time(ps)  Energy_autocorrelation\n")
                        for i in range(len(autocorrelation)):
                            f.write('%d %3.5f %3.5f\n' % (i, i*timestep, autocorrelation[i]))
                        f.close()
                    #print("timestep:", timestep)
                    #phonon_life_time[e] = np.sum(autocorrelation) * timestep
                    myindex = np.argmax(autocorrelation < 0.3679)
                    if myindex > 0:
                        phonon_life_time[wave_vector_index][e] = myindex * timestep
                    else:
                        phonon_life_time[wave_vector_index][e] = -1
                    #E[e] = np.subtract(E[e], np.mean(E[e]))
                    #autocorrelation = np.correlate(E[e], E[e], mode="full")
                    #autocorrelation = autocorrelation[autocorrelation.size//2:]
                    #autocorrelation = np.divide(autocorrelation, count)
                    #normlize_factor = autocorrelation[0]
                    #autocorrelation = np.divide(autocorrelation, normlize_factor)
                    #print(autocorrelation[0:50])
                    #phonon_life_time_subtract_mean[e] = np.sum(autocorrelation) * timestep

                    #if ll == m == n == 0:
                    print("[Eigenvector]:", e, "/", atom_num*3, " [time_step_per_frame]:", timestep,
                        " Phonon_life_time:", phonon_life_time[wave_vector_index][e])
                    #print("[Eigenvector]:", e, "/", atom_num*3, " [time_step_per_frame]:", timestep,
                    #      " Phonon_life_time:", phonon_life_time[e], " Phonon_life_time(subtract_mean):",
                    #      phonon_life_time_subtract_mean[e])
                    #exit()
                #print(phonon_life_time)
    print("phono_life_time_shape:", phonon_life_time.shape)
    #assert phonon_life_time.shape == (wave_vector_index, atom_num * 3)
    f = open(outputfilename + ".phonon_lifetime", 'w')
    #f.write('# omega/2Pi(THz) lifetime(ps) lifetime(ps)[Energy mean subtracted] k=0\n')
    f.write('# omega/2Pi(THz) lifetime(ps) num_unit_cell_X-Y-Z: %d %d %d \n' % (num_unit_cell[0], num_unit_cell[1], num_unit_cell[2]))
    for i in range(atom_num * 3):
        if w_value[0, i] < 0:
            f.write('%3.2f  ' % (-sqrt(-w_value[0, i])/2.0/pi))
        else:
            f.write('%3.2f  ' % (sqrt(w_value[0, i])/2.0/pi))
        for j in range(number_of_k):
            #if phonon_life_time[j][i] < 0.1:
            #    print(j, i, phonon_life_time[j][i])
            #    exit()
            #assert phonon_life_time[j][i] >= 0
            if phonon_life_time[j][i] < 0:
                continue
            f.write(' %3.3f' % (phonon_life_time[j][i]))
        f.write('\n')
    f.close()
    return

  # get eigenvalue
# def LAPACK_dsyev(int N, double *a, double *w):
#       MKL_INT  n = N, lda = N, info;
#           /* Executable statements */
#           printf( "\n#########################################\nDiagonalizing Dynamical Matrix using LAPACKE_dsyev, please wait for a few seconds ... ...\n\n\n" );
#           /* Solve eigenproblem */
#           info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, lda, w );
#           if( info > 0 ) {
#               printf( "The algorithm failed to compute eigenvalues.\n" );
#               exit( 1 );
#           }
#           /* Print eigenvalues */
#           //print_matrix( "Eigenvalues", 1, n, w, 1 );
#           /* Print eigenvectors */
#           //print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
#       printf("Print eigenvectors a[][1]\n");
#       for (int i=0; i < 10; i++){
#           printf("%lg ", a[n*i+1]);
#       }
#       printf("\n");
#   }

def lorentzian( x, x0, a, gam ):
    return a * gam**2 / ((x - x0)**2 + gam**2)

def lorentzian_log10( x, x0, a, gam ):
    return np.log10(a * gam**2 / ((x - x0)**2 + gam**2))

# function 5  get phonon life time from Lorentzian function
def get_phonon_life_time2(dyn, velocity_array, frame_num, atom_num, timestep, outputfilename,
                         num_unit_cell, unit_cell_size, cell_position):
    #|  for e in range(atom_num * 3):
    #|      if e == 0:
    #|          print("the mean position is: ", np.mean(position_array[e]))
    #|          print(position_array[0, 0:9])
    #|      position_array[e] = np.subtract(position_array[e], np.mean(position_array[e]))
    #|      if e == 0:
    #|          print(position_array[0, 0:9])
    from numpy import linalg as LA

    import time
    pi = 3.1415926
    # eigenvalue, eigenvector from diagonalization of dynamical matrix
    print("Starting diagonalization of the dynamical matrix ....")
    start_time = time.time()
    w_value, e_vector = LA.eigh(dyn)
    print("Finishing diagonalization of the dynamical matrix :)")
    print("Time used: ", time.time() - start_time, " seconds!")
    idx = np.argsort(-w_value)
    #print(w_value[0:10])
    w_value = np.array([w_value[idx]])
    print(w_value[0, 0:10])
    #print(e_vector[0:10, 0:10])
    e_vector = e_vector[:, idx]
    e_vector = np.conj(e_vector)
    print(velocity_array.shape)

    assert len(e_vector) == atom_num * 3
    count = np.arange(frame_num, 0, -1)
    #print(count)
    num_unit_cell = [1, 1, 1]
    number_of_k = num_unit_cell[0] * num_unit_cell[1] * num_unit_cell[2]
    phonon_life_time = np.zeros((atom_num * 3, number_of_k), dtype=float)
    #phonon_life_time_subtract_mean = np.zeros(atom_num * 3, dtype=float)
    print(w_value[0, 0:10])
    print(w_value.shape)
    #print(w_value[0, 0:2])
    for ll in range(0, num_unit_cell[0]):
        for m in range(0, num_unit_cell[1]):
            for n in range(0, num_unit_cell[2]):
                wave_vector_index = ll + num_unit_cell[0] * m + num_unit_cell[0] * num_unit_cell[1] * n

                #l1, m1, n1 = ll / num_unit_cell[0], m / num_unit_cell[1], n / num_unit_cell[2]
                k = np.tile([[2 * pi / num_unit_cell[0] / unit_cell_size[0] * ll,
                              2 * pi / num_unit_cell[1] / unit_cell_size[1] * m,
                              2 * pi / num_unit_cell[2] / unit_cell_size[2] * n]], atom_num)
                print("wave_vector_index:", wave_vector_index)
                print(unit_cell_size, ll, m, n)
                #print(cell_position[0:10])
                #print(k[0, 0:10])
                exp_term = np.transpose(np.exp(np.multiply(k, cell_position)))
                print(np.transpose(exp_term)[0, 0:10])
                #time.sleep(10)
                #position_array_temp = np.multiply(position_array, exp_term)
                velocity_array_temp = np.multiply(velocity_array, exp_term)
                #q = np.dot(np.transpose(e_vector), position_array_temp)
                q_d = np.dot(np.transpose(e_vector), velocity_array_temp)

                #| # get potential energy for each mode
                #| print("get potential energy for each mode", "  wave vector index: ", ll, "-", m, "-", n)
                #| q = np.multiply(np.conj(q), q)
                #| q = np.multiply(q, np.transpose(w_value)).real
                #| #print(q[0:2,0:2])
                #| q = np.multiply(q, 0.5)
                #| #print(q[0:2,0:2])
                #| #exit()
                from scipy.optimize import curve_fit
                import matplotlib
                #matplotlib.use('TkAgg')
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                n_col = q_d.shape[1]
                #q_d_fourier = np.zeros(n_col, dtype=float)
                #freq0 = np.arange(0, n_col, dtype=float)
                #for e in [0,10,20,30,200,500]:
                print("Do fourier transformation ... ")
                start_time = time.time()
                q_d_fourier = np.fft.fft(q_d, axis=-1)
                print("Fourier transformation completed!!!")
                print("Time used: ", time.time() - start_time, " seconds!")
                constants = 1 / 4.0 / pi / (n_col)
                freq = np.fft.fftfreq(n_col) * 2.0 * pi / timestep
                freq = freq[0:n_col//2 + 1]
                for e in range(atom_num * 3):
                    phi_all = np.multiply(np.conj(q_d_fourier[e, :]), q_d_fourier[e, :])
                    phi_all = np.multiply(phi_all, constants)
                    phi = phi_all[0:n_col//2 + 1].real
                    from lmfit import Model
                    model = Model(lorentzian)
                    params = model.make_params(x0=freq[np.argmax(phi)], a=np.amax(phi), gam=1.0)
                    params['gam'].min = 0.00000001
                    params['x0'].min = freq[np.argmax(phi)] - 1
                    params['x0'].max = freq[np.argmax(phi)] + 1
                    gam = 0.000001
                    try:
                        result = model.fit(phi, params, x=freq)
                        x0, a, gam = result.params['x0'].value, result.params['a'].value, result.params['gam'].value
                        phonon_life_time[e][wave_vector_index] = 1 / (2 * gam)
                    except:
                        #plt.plot(freq2, phi2, '-ok')
                        plt.plot(freq, phi, '-ok')
                        #plt.savefig("lorentzian_failed." + str(e) + ".png")
                        plt.savefig("lorentzian_long." + str(e) + "-" + str(x0) + "-" + str(a) + "-" + str(gam) + ".png")
                        #print(x0," ", a, " ", gam)
                        plt.clf()
                        #exit()
                    if e % atom_num == 0 or e == 512 or 1 / (2 * gam) > 300:
                    #if e > -1 or e % atom_num == 0 or e == 512 or 1 / (2 * gam) > 2000:
                        plt.yscale('log')
                        plt.plot(freq, phi, '-ok')
                        #plt.plot(freq2, phi2, '-ok')
                        x2 = np.linspace(x0 - 5, x0 + 5, 500)
                        plt.plot(x2, lorentzian(x2, x0, a, gam), c='r', lw=2, ls='dashed')
                        if 1 / (2 * gam) > 300:
                            plt.savefig("lorentzian_long." + str(e) + "-" + str(x0) + "-" + str(a) + "-" + str(gam) + ".png")
                            f = open("freq_phi.data." + str(e), 'w')
                            for iii in range(len(freq)):
                                f.write('%3.10f %3.10f\n' % (freq[iii], phi[iii]))
                            f.close()
                        else:
                            plt.savefig("lorentzian." + str(e) + "-" + str(x0) + "-" + str(a) + "-" + str(gam) + ".png")
                        plt.clf()
                        print(x0, " ", a, " ", gam)

                    print("[Eigenvector]:", e, "/", atom_num*3, " [time_step_per_frame]:", timestep,
                          " Phonon_life_time:", phonon_life_time[e][wave_vector_index])
                #exit()

    print("phono_life_time_shape:", phonon_life_time.shape)
    #assert phonon_life_time.shape == (wave_vector_index, atom_num * 3)
    f = open(outputfilename + ".phonon_lifetime", 'w')
    #f.write('# omega/2Pi(THz) lifetime(ps) lifetime(ps)[Energy mean subtracted] k=0\n')
    f.write('# omega/2Pi(THz) lifetime(ps) num_unit_cell_X-Y-Z: %d %d %d \n' % (num_unit_cell[0], num_unit_cell[1], num_unit_cell[2]))
    for i in range(atom_num * 3):
        if np.amin(phonon_life_time[i]) < 0:
            continue
        if w_value[0, i] < 0:
            f.write('%3.2f  ' % (-sqrt(-w_value[0, i])/2.0/pi))
        else:
            f.write('%3.2f  ' % (sqrt(w_value[0, i])/2.0/pi))
        for j in range(number_of_k):
            #if phonon_life_time[j][i] < 0.1:
            #    print(j, i, phonon_life_time[j][i])
            #    exit()
            #assert phonon_life_time[j][i] >= 0
            #if phonon_life_time[i][j] < 0:
            #    continue
            f.write(' %3.3f' % (phonon_life_time[i][j]))
        f.write('\n')
    f.close()
    return



def sum_rule(c, coor_frame_num):
    if len(c.shape) == 1:
        return c[:coor_frame_num].sum() - 0.5 * (c[0] + c[coor_frame_num-1])
    return c[:, :coor_frame_num].sum(axis=1) - 0.5 * (c[:, 0] + c[:, coor_frame_num-1])

def sum_rule_ave(c, coor_frame_num, middle):
    #assert middle < coor_frame_num
    if middle >= coor_frame_num:
        return sum_rule(c, coor_frame_num)
    max_repeat = coor_frame_num - middle
    mode_num = c.shape[0]
    count_xyz = np.arange(max_repeat, 0, -1)
    if len(c.shape) == 1:
        return c[:middle].sum() - 0.5 * c[0] + np.multiply(c[middle:coor_frame_num], count_xyz).sum() / (coor_frame_num - middle)
    return c[:,:middle].sum(axis=1) - 0.5 * c[:, 0] + np.multiply(c[:, middle:coor_frame_num], count_xyz).sum(axis=1) / (coor_frame_num - middle)


# for function 6 to do Green-Kubo Modal Analysis (GKMA), Wei Lv & Asegun Henry, Scientific Report, 2016, DOI: 10.1038/srep35720)
def GKMA_fast(dyn_w_value,dyn_e_vector, mass_sqr_array, velocity_ms_array, energy_array, stress_tensor_array, Q_0, volume, temperature,
         frame_num, coor_frame_num, strong_relation_frame, atom_num, time_step_per_frame, outputfilename):

    from numpy.fft import rfft, irfft
    # task:
    # 1. delete position_array
    print(Q_0.shape)
    assert Q_0.shape == (frame_num, 3)
    Q_0 = Q_0.transpose()
    from numpy import linalg as LA
    import time
    pi = 3.1415926
    # eigenvalue, eigenvector from diagonalization of dynamical matrix
    #print("Starting diagonalization of the dynamical matrix ....")
    #start_time = time.time()
    w_value, e_vector = dyn_w_value, dyn_e_vector
    #print("w_value, e_vector:", w_value, e_vector)
    #print("Finishing diagonalization of the dynamical matrix :)")
    #print("Time used: ", time.time() - start_time, " seconds!")
    idx = np.argsort(w_value)
    #idx = np.argsort(-w_value)
    #print(w_value[0:10])
    w_value = np.array([w_value[idx]])[0]
    e_vector = e_vector[:, idx]
    f_value = np.where(w_value < 0.000, -np.sqrt(-w_value)/2.0/pi, np.sqrt(w_value)/2.0/pi)
    #print(w_value[0, 0:10])
    #print(e_vector[0:10, 0:10])

    #freq_step = 0.15 # THz
    freq_step = 45 # THz
    start_freq = -1.0 # THz
    end_freq = 50  # Thz
    # combine mode within same frequency range
    new_f_vector = []
    new_e_vector = []

    #i = 0
    #while start_freq + i * freq_step < end_freq:
    #    fmin = start_freq + i * freq_step
    #    fmax = start_freq + (i+1) * freq_step
    #    id_list = np.where(np.logical_and(f_value >= fmin, f_value < fmax))[0]
    #    print("id_list:", len(id_list))
    #    if len(id_list) > 0:
    #        new_f_vector.append(np.sum(f_value[id_list[0]:(id_list[-1]+1)])/len(id_list))
    #        new_e_vector.append(np.sum(e_vector[:, id_list[0]:(id_list[-1]+1)], axis=1))
    #    i += 1.0
    #new_e_vector = np.array(new_e_vector).T
    #new_f_vector = np.array(new_f_vector)
    #mode_num = new_f_vector.shape[0]
    #print(new_f_vector.shape, new_e_vector.shape)
    #print(new_f_vector[:10])
    #print(new_e_vector.sum(axis=1)[:10])
    #print(e_vector.sum(axis=1)[:10])
    ##exit()
    #e_vector = new_e_vector

    new_f_value = []
    e_vector_conj = np.conj(e_vector)
    mode_num = e_vector.shape[1]

    print("Here")
    print(e_vector.shape, e_vector_conj.shape)
    #exit()
    #print(position_array.shape)
    print("velocity shape:", velocity_ms_array.shape)
    #mass_array = mass_array.reshape(atom_num*3, 1)
    stress_tensor_array = stress_tensor_array.transpose().reshape(frame_num, -1, 3, 3)
    print("stress tensor shape:", stress_tensor_array.shape)
    #exit()

    # Unit conversion
    # Warning!!!
    # In lammps, the unit of atomic stress from command(compute stress/atom) is pressure * volume.
    # So you need to check the unit of pressure and volume in the unit system you use
    # For example, if unit = metal, pressure * volume = bars * Å^3, this is not directly equal to the energy unit eV
    # Instead, pressure * volume = bars * Å^3 = 6.25e-7 eV ---> second_term_unit_conversion
    # this above gonna affect the 2nd term of the heat flux sum(S*v)
    second_term_unit_conversion = 1/1.6021765e6  # see above explantion
    kB=1.3806504e-23  # J/K
    eV2J = 1.602e-19
    A2m = 1.0e-10
    ps2s = 1.0e-12
    #timestep per frame in unit of fs
    dt2s = 1.0e-15
    T = temperature
    convert = eV2J*eV2J/ps2s/ps2s*dt2s/A2m
    debug = False
    #debug = True

    #Q_sum store heat total flux along x/y/z from all modes
    Q_sum = np.zeros((3,frame_num))
    Q_sum_convex = np.zeros((3,frame_num))
    # Start treat each mode separately
    print("w_value shape:", w_value.shape)
    #print(w_value[0,:].shape)
    #exit()
    #kappa = np.zeros(mode_num, dtype=float)
    kappa = []
    accum_kappa = 0.0
    #accum_kappa_xyz = [0.0, 0.0, 0.0]
    padding = np.zeros((3, frame_num))
    count_xyz = np.tile(np.arange(frame_num, 0, -1), (3, 1))
    Q_0_pad = np.concatenate((Q_0, padding), axis=1)
    bft = rfft(Q_0_pad, axis=1)
    f = open(outputfilename + ".GKMA.data", 'w')
    f.write('# mode_middle_frequency(THz) mode_omega/2Pi(THz) modal_kappa(W/m/k)  accumulated_kappa(W/m/k)  num_modes_within_range\n')
    # Q_nth_mode_along_time
    Q_n_t_all_x = []
    Q_n_t_all_y = []
    Q_n_t_all_z = []
    begin_i, end_i = 0, mode_num
    print("mode_num=", mode_num)
    #exit()
    #strong_relation_frame = 1000
    debug_time = False
    start_freq = f_value[0]
    freq_step = 0.20 # THz
    start_freq = -20.0 # THz
    end_freq = 50  # Thz
    num_step = int((end_freq - start_freq) // freq_step + 1)
    print("num_step", num_step)
    for j in range(num_step):
        cycle_start_t = time.time()
        fmin = start_freq + j * freq_step
        if fmin > end_freq:
            continue
        fmax = start_freq + (j+1) * freq_step
        id_list = np.where(np.logical_and(f_value >= fmin, f_value < fmax))[0]
        if len(id_list) == 0:
            continue
        ave_frequency = f_value[id_list[0]:(id_list[-1]+1)].sum()/len(id_list)
        middle_frequency = start_freq + (j+0.5) * freq_step
        #new_f_value.append(ave_frequency)
        new_f_value.append(middle_frequency)
        x_j_n_t = None
        for i in id_list:
            if debug_time:
                t1 = time.time()
                t11 = time.time()
                t12 = time.time()
                print("t12-t11:", t12-t11)
            p = np.divide(e_vector[:, i], mass_sqr_array)
            p_s = e_vector_conj[:, i]  # conjugate: p_s=p*
            #assert np.array_equal(e_vector[:, i], np.conj(p_s))
            if debug_time:
                t13 = time.time()
                print("t13-t12:", t13-t12)
            # p_s shapes become [1 x (number of atoms*3)]
            # get X_n_t_p
            #print("p_s shape:", p_s.shape)
            X_n_t_p = np.matmul(p_s, velocity_ms_array)  # matmul will be 1-D array
            if debug_time:
                t14 = time.time()
                print("t14-t13:", t14-t13)
            #print("X_n_t_p",X_n_t_p.shape)
            if debug:
                #print(p_s, velocity_ms_array)
                print("X_n_t_p:", X_n_t_p)
            assert X_n_t_p.shape == (frame_num,)  # Row=mode; Col=Time
            # x_j_n_t is x_dot_j(n, t)
            if x_j_n_t is None:
                x_j_n_t = np.multiply(p.reshape([atom_num * 3, 1]), X_n_t_p)
            else:
                x_j_n_t += np.multiply(p.reshape([atom_num * 3, 1]), X_n_t_p)
            #  print(x_j_n_t.shape)
            #  print(velocity_ms_array.shape)
            #  print(mass_sqr_array.shape)
            #  print(np.divide(velocity_ms_array[:, 0], mass_sqr_array)[:10])
            #  print(x_j_n_t[:10,0])
            #  print((p*mass_sqr_array)[:10])
            #  print(p_s[:10])
            #  print("mode_num=", mode_num)
            #  exit()
            if debug_time:
                t15 = time.time()
                print("t15-t14:", t15-t14)
                print("t15-t1(mode_sum):", t15-t1)
            #print("p", p)
            #exit()
            #print("x_j_n_t shape:", x_j_n_t.shape)
            if debug: print("x_j_n_t:", x_j_n_t)

        after_sum_x = time.time()
        print("sum_of_x time:", after_sum_x - cycle_start_t, " for #modes:", len(id_list))
        x_j_n_t_3 = x_j_n_t.transpose().reshape((frame_num, -1, 3)).repeat(3, axis=1).reshape(frame_num, -1, 3, 3)
        #print("x_j_n_t_3 shape:", x_j_n_t_3.shape)
        if debug: print("x_j_n_t_3 :", x_j_n_t_3)
        #exit()
        #print("mark")
        #print("energy shape:", energy_array.shape)
        # Q(n,t) = 1/Vol*[sum(E*v) - sum(S*v)] # I will not include 1/Vol but leave it to the integral
        #(1) sum(E*v)
        #first_term = np.multiply(energy_array, x_j_n_t).sum(axis=0)  # element-wise multiply + sum along axis=0
        first_term = np.multiply(energy_array, x_j_n_t).transpose().reshape((frame_num, atom_num, -1)).sum(axis=1)
        if debug_time:
            t16 = time.time()
            print("t16-t15:", t16-t15)
        first_term = first_term.T
        if debug:
            print(energy_array)
            print(x_j_n_t)
            print("first_term:", first_term)
            print(x_j_n_t.sum(axis=0))
            print("first_term:", first_term)
        #print("first_term shape:", first_term.shape)
        #print("stress tensor array", stress_tensor_array.shape)
        #(2) sum(S*v)
        second_term = np.multiply(stress_tensor_array, x_j_n_t_3).sum(axis=(1, 3)).transpose()
        #print("first_term & second_term\n", first_term, second_term)
        #exit()
        if debug_time:
            t17 = time.time()
            print("t17-t16:", t17-t16)

        if debug:
            print(stress_tensor_array)
            print(x_j_n_t_3)
            print("second_term:", second_term)
            exit()
        #print("second term shape:", second_term.shape)
        Q_n_t = (first_term - second_term * second_term_unit_conversion)
        if debug_time:
            t2 = time.time()
            print("t2-t1:", t2-cycle_start_t)
        #print(Q_n_t[0,:10])
        #exit()
        Q_n_t_all_x.append(Q_n_t[0, :])
        Q_n_t_all_y.append(Q_n_t[1, :])
        Q_n_t_all_z.append(Q_n_t[2, :])
        #Q_sum = Q_sum + Q_n_t
        #Q_sum_convex = Q_sum_convex + first_term
        #Q_n_t = (first_term - second_term)
        #Q_n_t = first_term
        #Q_n_t = second_term
        if debug:
            print("first_term & second_term", first_term, second_term)
        #print("Q_n_t shape:", Q_n_t.shape)
        if debug:
            print(Q_n_t)
            exit()
        ##In GK relation, the integral The integral is over the equilibrium flux autocovariance function;
        ## so we need to deduct the mean value of heat flux
        count = np.arange(len(Q_n_t[0]), 0, -1)
        autocorrelation = None

        # use fast fourier
        #Q_n_t_pad = np.concatenate((Q_n_t, padding), axis=1)
        aft = rfft(Q_n_t, n=(frame_num*2), axis=1)
        c = irfft(aft*np.conjugate(bft), axis=1)
        c = c[:, :frame_num]
        c = np.divide(c, count_xyz)
        #assert c.shape[0] == 3
        #autocorrelation = c.sum(axis=1) - 0.5 * c[:, 0]
        #autocorrelation = c.sum(axis=1)
        autocorrelation = c[:,:coor_frame_num].sum(axis=0)
        #print("autocor shape:", autocorrelation.shape)
        autocorrelation = autocorrelation / 3.0
        #print(autocorrelation)
        if debug_time:
            t3=time.time()
            print("t3-t2:", t3-t2)
        if debug:
            print("autocorr:", autocorrelation)
            print("autocorr.sum:", autocorrelation.sum())
        #autocorrelation = autocorrelation[:coor_frame_num]
        #kappa[i] = (autocorrelation.sum(axis=0) - 0.5 * autocorrelation[0]) * \
        kappa.append(sum_rule_ave(autocorrelation, coor_frame_num, strong_relation_frame) * \
                   time_step_per_frame/volume/kB/T/T*convert)  # in Unit of W/m/k
        #print("value:",time_step_per_frame/volume/kB/T/T*convert)
        #exit()

        accum_kappa += kappa[-1]
        if debug_time:
            t4=time.time()
            print("t4-t3:", t4-t3)
        #  | for k in range(3):
        #  |     accum_kappa_xyz[k] += (autocorrelation_xyz[k,:-1].sum() - 0.5 * autocorrelation_xyz[k,0]) * \
        #  |         time_step_per_frame/volume/kB/T/T*convert  # in Unit of W/m/k

        #assert w_value[0, i] != 0
        #print(w_value[0, i])
        print("Calculating_kappa_for_mode_nth:", j, "/", num_step, "  THz ", ave_frequency,  " kappa(W/m/k): ", kappa[-1], "mode_num:", len(id_list))
        print("cycle time:", time.time()-cycle_start_t, "heat_flux_time:", time.time()-after_sum_x)
        #exit()
        #||  if i==10 or i == 10 or i + 1 == e_vector_conj.shape[0]:
        #||      print(kappa[:10])
        #||      exit()
        #normlize_factor = autocorrelation[0]
        #autocorrelation = np.divide(autocorrelation, normlize_factor)
        #exit()
        #for i in range(atom_num * 3):
        f.write('%3.2f %3.2f  %3.3f  %3.3f %d\n' % (middle_frequency, ave_frequency, kappa[-1], accum_kappa, len(id_list)))
        #f.write(' %3.3f  %3.3f  %3.3f  %3.3f %3.3f' % (kappa[i], accum_kappa, accum_kappa_xyz[0], accum_kappa_xyz[1], accum_kappa_xyz[2]))
        if debug_time:
            t5=time.time()
            print("Full cycle t5-t1:", t5-t1)
        #exit()
    #print("Q_sum: ", Q_sum)
    f.close()
    with open('frequency_THz_increasing.npy', 'wb') as f:
        np.save(f, w_value)
    with open('grouped_ave_frequency_THz_increasing.npy', 'wb') as f:
        np.save(f, np.array(new_f_value))
    with open('Q_0_with_padding.npy', 'wb') as f:
        np.save(f, Q_0_pad)
        #with open('Q_n_t.'+str(mode_start)+'-'+str(mode_end)+'.npy', 'wb') as f:
        #np.save(f, Q_n_t)

    Q_n_t_all_x = np.array(Q_n_t_all_x)
    Q_n_t_all_y = np.array(Q_n_t_all_y)
    Q_n_t_all_z = np.array(Q_n_t_all_z)
    print(Q_n_t_all_x.sum(axis=0)[:10])
    print(Q_n_t_all_y.sum(axis=0)[:10])
    print(Q_n_t_all_z.sum(axis=0)[:10])
    import pickle as pk
    import gzip
    f_direction=pk.dump(Q_n_t_all_x, gzip.open('Q_n_t_all_x.gz', 'wb'))
    f_direction=pk.dump(Q_n_t_all_y, gzip.open('Q_n_t_all_y.gz', 'wb'))
    f_direction=pk.dump(Q_n_t_all_z, gzip.open('Q_n_t_all_z.gz', 'wb'))
    #with open('Q_n_t_all_y.npy', 'wb') as f:
    #    np.save(f, Q_n_t_all_y)
    #with open('Q_n_t_all_z.npy', 'wb') as f:
    #    np.save(f, Q_n_t_all_z)
    w_value_THz = np.where(w_value < 0, -np.sqrt(-w_value)/2.0/pi, np.sqrt(w_value)/2.0/pi)
    data=np.concatenate((np.arange(len(w_value)).reshape(1,len(w_value)), w_value.reshape(1,-1)), axis=0)
    np.savetxt('frequency_THz_increasing.dat', data.T, fmt='%.3f', delimiter=' ')

    print("job finished")
    return
    #print("Q_sum_convex: ", Q_sum_convex)
    autocor_xyz = np.zeros((3, frame_num))
    for j in range(3):
        autocor = np.correlate(Q_sum[j, :], Q_sum[j, :], mode="full")
        autocor = autocor[autocor.size//2:]
        #autocor = np.divide(autocor, count)
        autocor_xyz[j,:] = np.divide(autocor, count)
        print(autocor_xyz[j,:])
    #exit()

    # Note in lammps generated correlation file (e.g. profile.heatflux), the correlation is only up to (framenumber-1)
    # so to match kappa, I use autocor_xyz[:,:-1] or [:,:coor_frame_num]
    autocor_xyz = autocor_xyz[:,:coor_frame_num]
    temp = (autocor_xyz.sum(axis=1) - 0.5 * autocor_xyz[:,0]- 0.5 * autocor_xyz[:,-1]) * time_step_per_frame/volume/kB/T/T*convert
    #print(Q_sum)
    print("conversion:", time_step_per_frame/volume/kB/T/T*convert)  # in Unit of W/m/k
    print("convert:", convert)
    print("kappa from Q_sum (W/m/k) in x/y/z:", temp) # in Unit of W/m/k
    #exit()
    #################################
    ## calculate all modal kappa at one time
    ###################################
    #Q0 = np.tile(np.arange(frame_num, 0, -1), (Q_n_t_all_x.shape[0], 1))
    # for j in range(3):
    #     #Q_n_t[j, :] = Q_n_t[j, :] - Q_n_t[j,:].mean()
    #     #print(Q_n_t.shape, j, Q_n_t[0, :])
    #     autocor = np.correlate(Q_n_t[j, :], Q_0[j, :], mode="full")
    #     #print(autocor)
    #     #exit()
    #     autocor = autocor[autocor.size//2:]
    #     autocor = np.divide(autocor, count)
    #     if j == 0:
    #         autocorrelation = autocor
    #     else:
    #         autocorrelation += autocor
    #     #autocorrelation_xyz[j] = autocor
    # autocorrelation = autocorrelation / 3.0

    #exit()
    autocor_xyz = np.concatenate((np.arange(autocor_xyz.shape[1]).reshape(1,autocor_xyz.shape[1]), autocor_xyz), axis=0)
    np.savetxt('Q_sum.out', Q_sum.transpose(), fmt='%.2f', delimiter=' ')
    np.savetxt('Q_sum_convection.out', Q_sum_convex, fmt='%.2f', delimiter=' ')
    np.savetxt('Q_sum_autocor.out', autocor_xyz.transpose(), fmt='%3.2f', delimiter=' ')

    #print((energy_array[:, 0]*velocity_ms_array[:,0]/mass_sqr_array).sum())
    # calculate cross-correlation of heat flux
    print("Start calculate cross-correlation of heat flux\n")
    # padding with zeros
    padding = np.zeros((e_vector_conj.shape[0], frame_num))

    #Q_n_t_all_x = np.concatenate((Q_n_t_all_x, padding), axis=1)
    print(Q_n_t_all_x.shape)
    count = np.tile(np.arange(frame_num, 0, -1), (Q_n_t_all_x.shape[0], 1))
    # Q Q cross correlation
    Q_Q_coor = np.zeros((atom_num * 3, atom_num * 3))
    print_few_k = False
    if print_few_k:
        ka=0.0
        kb=0.0
    sum_all_QQ = 0.0
    # debug
    print(Q_n_t_all_x.sum(axis=0)[:10])
    print(Q_sum[0,:10])
    #exit()
    # use fast fourier transform to get the cross correlation along x/y/z directions
    aft = rfft(Q_n_t_all_x, n=frame_num*2, axis=1)
    bft = np.conjugate(aft)
    row_num = Q_n_t_all_x.shape[0]
    for i in range(row_num):
        print("calculate cross correlation ith:", i)
        c = irfft(aft[:row_num-i, ]*bft[i:, ], axis=1)
        c = c[:, :frame_num]
        #if i + 1 == Q_n_t_all_x.shape[0]:
        #    c = c[0:frame_num]
        #else:
        #    c = c[:, :frame_num]
        #print(c, count[i:,])
        c = np.divide(c, count[i:, ])

        # Note in lammps generated correlation file (e.g. profile.heatflux), the correlation is only up to (framenumber-1)
        # so to match kappa, I use [:,:coor_frame_num] or [:,:-1]
        cc_mat = c[:,:coor_frame_num].sum(axis=1) - 0.5 * c[:, 0]

        if i > 0:
            bottom = irfft(aft[i:, ]*bft[:row_num-i, ], axis=1)
            bottom = bottom[:, :frame_num]
            bottom = np.divide(bottom, count[i:, ])
            bot_mat = bottom[:,:coor_frame_num].sum(axis=1) - 0.5 * bottom[:, 0]

        if i == 0:
            sum_all_QQ += cc_mat.sum()
        else:
            sum_all_QQ += cc_mat.sum() + bot_mat.sum()
        if print_few_k:
            if i ==0:
                ka=cc_mat[0]
            if i ==1:
                kb=cc_mat[0]
        if i == 0:
            np.fill_diagonal(Q_Q_coor, cc_mat)
        else:
            np.fill_diagonal(Q_Q_coor[:-i, i:], cc_mat)
    print("Q_Q_coor shape:", Q_Q_coor.shape)
    with open('QQ_coorelation.npy', 'wb') as f:
        np.save(f, w_value[0, :])
        np.save(f, Q_Q_coor)
    np.savetxt('QQ_coorelation.dat', Q_Q_coor, fmt='%.2f', delimiter=' ')
    print("kappa from sum_all_QQ_along_x[This value can be different from calcualtion using"
          " lammps autocorrelation \nbecause only upper QQ correlation matrix is computed and diagonal"
          " is assumed but is not guaranteed] (W/m/k)\n:", sum_all_QQ * time_step_per_frame/volume/kB/T/T*convert,
          " (W/m/k) \n")  # in Unit of W/m/k
    print("convert:", convert)
    if print_few_k:
        print("Q[0,0]:",ka, " Q[0,1]:",kb)
    return



# for function 6 to do Green-Kubo Modal Analysis (GKMA), Wei Lv & Asegun Henry, Scientific Report, 2016, DOI: 10.1038/srep35720)
def GKMA(dyn_w_value,dyn_e_vector, mass_sqr_array, velocity_ms_array, energy_array, stress_tensor_array, Q_0, volume, temperature,
         frame_num, coor_frame_num, atom_num, time_step_per_frame, outputfilename, mode_start, mode_end):

    from numpy.fft import rfft, irfft
    # task:
    # 1. delete position_array
    print(Q_0.shape)
    assert Q_0.shape == (frame_num, 3)
    Q_0 = Q_0.transpose()
    from numpy import linalg as LA
    import time
    pi = 3.1415926
    # eigenvalue, eigenvector from diagonalization of dynamical matrix
    #print("Starting diagonalization of the dynamical matrix ....")
    #start_time = time.time()
    w_value, e_vector = dyn_w_value, dyn_e_vector
    #print("w_value, e_vector:", w_value, e_vector)
    #print("Finishing diagonalization of the dynamical matrix :)")
    #print("Time used: ", time.time() - start_time, " seconds!")
    idx = np.argsort(w_value)
    #idx = np.argsort(-w_value)
    #print(w_value[0:10])
    w_value = np.array([w_value[idx]])
    #print(w_value[0, 0:10])
    #print(e_vector[0:10, 0:10])
    e_vector = e_vector[:, idx]
    e_vector_conj = np.conj(e_vector)
    print("Here")
    print(e_vector.shape, e_vector_conj.shape)
    #print(position_array.shape)
    print("velocity shape:", velocity_ms_array.shape)
    #mass_array = mass_array.reshape(atom_num*3, 1)
    stress_tensor_array = stress_tensor_array.transpose().reshape(frame_num, -1, 3, 3)
    print("stress tensor shape:", stress_tensor_array.shape)
    #exit()

    # Unit conversion
    # Warning!!!
    # In lammps, the unit of atomic stress from command(compute stress/atom) is pressure * volume.
    # So you need to check the unit of pressure and volume in the unit system you use
    # For example, if unit = metal, pressure * volume = bars * Å^3, this is not directly equal to the energy unit eV
    # Instead, pressure * volume = bars * Å^3 = 6.25e-7 eV ---> second_term_unit_conversion
    # this above gonna affect the 2nd term of the heat flux sum(S*v)
    second_term_unit_conversion = 1/1.6021765e6  # see above explantion
    kB=1.3806504e-23  # J/K
    eV2J = 1.602e-19
    A2m = 1.0e-10
    ps2s = 1.0e-12
    #timestep per frame in unit of fs
    dt2s = 1.0e-15
    T = temperature
    convert = eV2J*eV2J/ps2s/ps2s*dt2s/A2m
    debug = False
    #debug = True

    #Q_sum store heat total flux along x/y/z from all modes
    Q_sum = np.zeros((3,frame_num))
    Q_sum_convex = np.zeros((3,frame_num))
    # Start treat each mode separately
    print("w_value shape:", w_value.shape)
    #print(w_value[0,:].shape)
    #exit()
    kappa = np.zeros(atom_num*3, dtype=float)
    accum_kappa = 0.0
    #accum_kappa_xyz = [0.0, 0.0, 0.0]
    padding = np.zeros((3, frame_num))
    count_xyz = np.tile(np.arange(frame_num, 0, -1), (3, 1))
    Q_0_pad = np.concatenate((Q_0, padding), axis=1)
    bft = rfft(Q_0_pad, axis=1)
    f = open(outputfilename + ".GKMA.data", 'w')
    f.write('# mode_omega/2Pi(THz) modal_kappa(W/m/k)  accumulated_kappa(W/m/k)\n')
    # Q_nth_mode_along_time
    Q_n_t_all_x = np.zeros((e_vector_conj.shape[0], frame_num))
    Q_n_t_all_y = np.zeros((e_vector_conj.shape[0], frame_num))
    Q_n_t_all_z = np.zeros((e_vector_conj.shape[0], frame_num))
    begin_i, end_i = mode_start, mode_end
    if begin_i is None:
        begin_i, end_i = 0, e_vector_conj.shape[0]
        print("calculate all modes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        begin_i -= 1
    debug_time = False
    for i in range(begin_i, end_i):
        cycle_start_t = time.time()
        if debug_time:
            t1 = time.time()
            t11 = time.time()
            t12 = time.time()
            print("t12-t11:", t12-t11)
        p = np.divide(e_vector[:, i], mass_sqr_array)
        p_s = e_vector_conj[:, i]  # conjugate: p_s=p*
        if debug_time:
            t13 = time.time()
            print("t13-t12:", t13-t12)
        if debug: print("p_s:", p_s)
        #exit()
        # p_s shapes become [1 x (number of atoms*3)]
        # get X_n_t_p
        #print("p_s shape:", p_s.shape)
        X_n_t_p = np.matmul(p_s, velocity_ms_array)  # matmul will be 1-D array
        if debug_time:
            t14 = time.time()
            print("t14-t13:", t14-t13)
        #print("X_n_t_p",X_n_t_p.shape)
        if debug:
            print(p_s, velocity_ms_array)
            print("X_n_t_p:", X_n_t_p)
        assert X_n_t_p.shape == (frame_num,)  # Row=mode; Col=Time
        # x_j_n_t is x_dot_j(n, t)
        x_j_n_t = np.multiply(p.reshape([atom_num * 3, 1]), X_n_t_p)
        if debug_time:
            t15 = time.time()
            print("t15-t14:", t15-t14)
        #print("p", p)
        #exit()
        #print("x_j_n_t shape:", x_j_n_t.shape)
        if debug: print("x_j_n_t:", x_j_n_t)
        x_j_n_t_3 = x_j_n_t.transpose().reshape((frame_num, -1, 3)).repeat(3, axis=1).reshape(frame_num, -1, 3, 3)
        #print("x_j_n_t_3 shape:", x_j_n_t_3.shape)
        if debug: print("x_j_n_t_3 :", x_j_n_t_3)
        #exit()
        #print("mark")
        #print("energy shape:", energy_array.shape)
        # Q(n,t) = 1/Vol*[sum(E*v) - sum(S*v)] # I will not include 1/Vol but leave it to the integral
        #(1) sum(E*v)
        #first_term = np.multiply(energy_array, x_j_n_t).sum(axis=0)  # element-wise multiply + sum along axis=0
        first_term = np.multiply(energy_array, x_j_n_t).transpose().reshape((frame_num, atom_num, -1)).sum(axis=1)
        if debug_time:
            t16 = time.time()
            print("t16-t15:", t16-t15)
        first_term = first_term.transpose()
        if debug:
            print(energy_array)
            print(x_j_n_t)
            print("first_term:", first_term)
            print(x_j_n_t.sum(axis=0))
            #exit()
            print("first_term:", first_term)
        #print("first_term shape:", first_term.shape)
        #print("stress tensor array", stress_tensor_array.shape)
        #(2) sum(S*v)
        second_term = np.multiply(stress_tensor_array, x_j_n_t_3).sum(axis=(1, 3)).transpose()
        #print("first_term & second_term\n", first_term, second_term)
        #exit()
        if debug_time:
            t17 = time.time()
            print("t17-t16:", t17-t16)

        if debug:
            print(stress_tensor_array)
            print(x_j_n_t_3)
            print("second_term:", second_term)
            exit()
        #print("second term shape:", second_term.shape)
        Q_n_t = (first_term - second_term * second_term_unit_conversion)
        if debug_time:
            t2 = time.time()
            print("t2-t1:", t2-t1)
        #Q_n_t_all_x[i, :] = Q_n_t.sum(axis=0)/3.0
        Q_n_t_all_x[i, :] = Q_n_t[0, :]
        Q_n_t_all_y[i, :] = Q_n_t[1, :]
        Q_n_t_all_z[i, :] = Q_n_t[2, :]
        Q_sum = Q_sum + Q_n_t
        Q_sum_convex = Q_sum_convex + first_term
        #Q_n_t = (first_term - second_term)
        #Q_n_t = first_term
        #Q_n_t = second_term
        if debug:
            print("first_term & second_term", first_term, second_term)
        #print("Q_n_t shape:", Q_n_t.shape)
        if debug:
            print(Q_n_t)
            exit()
        ##In GK relation, the integral The integral is over the equilibrium flux autocovariance function;
        ## so we need to deduct the mean value of heat flux
        count = np.arange(len(Q_n_t[0]), 0, -1)
        autocorrelation = None
        #autocorrelation_xyz = np.zeros((3, frame_num))

        # use fast fourier
        #Q_n_t_pad = np.concatenate((Q_n_t, padding), axis=1)
        aft = rfft(Q_n_t, n=(frame_num*2), axis=1)
        c = irfft(aft*np.conjugate(bft), axis=1)
        c = c[:, :frame_num]
        c = np.divide(c, count_xyz)
        #assert c.shape[0] == 3
        #autocorrelation = c.sum(axis=1) - 0.5 * c[:, 0]
        #autocorrelation = c.sum(axis=1)
        autocorrelation = c[:,:coor_frame_num].sum(axis=0)
        #print("autocor shape:", autocorrelation.shape)
        autocorrelation = autocorrelation / 3.0
        #print(autocorrelation)
        if debug_time:
            t3=time.time()
            print("t3-t2:", t3-t2)
        #for j in range(3):
        #    #Q_n_t[j, :] = Q_n_t[j, :] - Q_n_t[j,:].mean()
        #    #print(Q_n_t.shape, j, Q_n_t[0, :])
        #    autocor = np.correlate(Q_n_t[j, :], Q_0[j, :], mode="full")
        #    #print(autocor)
        #    #exit()
        #    autocor = autocor[autocor.size//2:]
        #    autocor = np.divide(autocor, count)
        #    if j == 0:
        #        autocorrelation = autocor
        #    else:
        #        autocorrelation += autocor
        #    #autocorrelation_xyz[j] = autocor
        #autocorrelation = autocorrelation / 3.0
        #print(autocorrelation)
        #exit()
        if debug:
            print("autocorr:", autocorrelation)
            print("autocorr.sum:", autocorrelation.sum())
            #sigma_square = autocorrelation[0]
            #print("covariance:", autocorrelation - sigma_square)
            #|  test = np.array([-2,3,0])
            #|  test = test - test.mean()
            #|  print("test:", test)
            #|  t2=np.correlate(test,test,mode="full")
            #|  t2 = t2[t2.size//2:]
            #|  t2 = np.divide(t2, count)
            #|  print("t2:", t2)
            #|  kappa[i] = autocorrelation.sum(axis=0) * time_step_per_frame/volume/1.0/T/T*1.0  # in Unit of W/m/k
        #autocovariance = autocorrelation - sigma_square
        #cal = (1.6e-19)**2/((1e-12)**2)*1e-15/1e-10/1.38e-23/90000
        #kappa[i] = autocorrelation.sum(axis=0) * time_step_per_frame  # in Unit of W/m/k
        #print("autocorre shape:", autocorrelation.shape)
        #print(autocorrelation[0], autocorrelation[-1])
        #print(autocorrelation[:5])
        #exit()
        # trapzoid rule to get sum, front and end only get half
        #print(autocorrelation)
        #print(autocorrelation.sum(axis=0) - 0.5 * (autocorrelation[0]+autocorrelation[-1]))
        #print(time_step_per_frame, volume, kB,T, convert)
        autocorrelation = autocorrelation[:coor_frame_num]
        #kappa[i] = (autocorrelation.sum(axis=0) - 0.5 * autocorrelation[0]) * \
        kappa[i] = (autocorrelation.sum(axis=0) - 0.5 * (autocorrelation[0]+autocorrelation[-1])) * \
                              time_step_per_frame/volume/kB/T/T*convert  # in Unit of W/m/k
        #print("value:",time_step_per_frame/volume/kB/T/T*convert)
        #exit()

        accum_kappa += kappa[i]
        if debug_time:
             t4=time.time()
             print("t4-t3:", t4-t3)
        #  | for k in range(3):
        #  |     accum_kappa_xyz[k] += (autocorrelation_xyz[k,:-1].sum() - 0.5 * autocorrelation_xyz[k,0]) * \
        #  |         time_step_per_frame/volume/kB/T/T*convert  # in Unit of W/m/k

        #assert w_value[0, i] != 0
        #print(w_value[0, i])
        if w_value[0, i] < 0:
            print("Calculating kappa for mode #:", i, " THz ", -sqrt(-w_value[0, i])/2.0/pi,  " kappa(W/m/k): ", kappa[i])
        else:
            print("Calculating kappa for mode #:", i, " THz ", sqrt(w_value[0, i])/2.0/pi,  " kappa(W/m/k): ", kappa[i])
        print("cycle time:", time.time()-cycle_start_t)
        #exit()
        #||  if i==10 or i == 10 or i + 1 == e_vector_conj.shape[0]:
        #||      print(kappa[:10])
        #||      exit()
        #normlize_factor = autocorrelation[0]
        #autocorrelation = np.divide(autocorrelation, normlize_factor)
        #exit()
        #for i in range(atom_num * 3):
        if w_value[0, i] < 0:
            f.write('%3.2f  ' % (-sqrt(-w_value[0, i]) / 2.0 / pi))
        else:
            f.write('%3.2f  ' % (sqrt(w_value[0, i])/2.0/pi))
        f.write(' %3.3f  %3.3f' % (kappa[i], accum_kappa))
        #f.write(' %3.3f  %3.3f  %3.3f  %3.3f %3.3f' % (kappa[i], accum_kappa, accum_kappa_xyz[0], accum_kappa_xyz[1], accum_kappa_xyz[2]))
        f.write('\n')
        if debug_time:
             t5=time.time()
             print("Full cycle t5-t1:", t5-t1)
        #exit()
    #print("Q_sum: ", Q_sum)
    f.close()
    if mode_start is not None:
        if mode_start == 1:
            with open('frequency_THz_increasing.npy', 'wb') as f:
                np.save(f, w_value[0, :])
            with open('Q_0_with_padding.npy', 'wb') as f:
                np.save(f, Q_0_pad)
        #with open('Q_n_t.'+str(mode_start)+'-'+str(mode_end)+'.npy', 'wb') as f:
            #np.save(f, Q_n_t)
        with open('Q_n_t_all_x.'+str(mode_start)+'-'+str(mode_end)+'.npy', 'wb') as f:
            np.save(f, Q_n_t_all_x[mode_start-1:mode_end, :])
        with open('Q_n_t_all_y.'+str(mode_start)+'-'+str(mode_end)+'.npy', 'wb') as f:
            np.save(f, Q_n_t_all_y[mode_start-1:mode_end, :])
        with open('Q_n_t_all_z.'+str(mode_start)+'-'+str(mode_end)+'.npy', 'wb') as f:
            np.save(f, Q_n_t_all_z[mode_start-1:mode_end, :])
        print("Only calculate mode range:", mode_start," to ", mode_end)
        print("job finished")
        return
    #print("Q_sum_convex: ", Q_sum_convex)
    autocor_xyz = np.zeros((3, frame_num))
    for j in range(3):
        autocor = np.correlate(Q_sum[j, :], Q_sum[j, :], mode="full")
        autocor = autocor[autocor.size//2:]
        #autocor = np.divide(autocor, count)
        autocor_xyz[j,:] = np.divide(autocor, count)
        print(autocor_xyz[j,:])
    #exit()

    # Note in lammps generated correlation file (e.g. profile.heatflux), the correlation is only up to (framenumber-1)
    # so to match kappa, I use autocor_xyz[:,:-1] or [:,:coor_frame_num]
    autocor_xyz = autocor_xyz[:,:coor_frame_num]
    temp = (autocor_xyz.sum(axis=1) - 0.5 * autocor_xyz[:,0]- 0.5 * autocor_xyz[:,-1]) * time_step_per_frame/volume/kB/T/T*convert
    #print(Q_sum)
    print("conversion:", time_step_per_frame/volume/kB/T/T*convert)  # in Unit of W/m/k
    print("convert:", convert)
    print("kappa from Q_sum (W/m/k) in x/y/z:", temp) # in Unit of W/m/k
    #exit()
    #################################
    ## calculate all modal kappa at one time
    ###################################
    #Q0 = np.tile(np.arange(frame_num, 0, -1), (Q_n_t_all_x.shape[0], 1))
    # for j in range(3):
    #     #Q_n_t[j, :] = Q_n_t[j, :] - Q_n_t[j,:].mean()
    #     #print(Q_n_t.shape, j, Q_n_t[0, :])
    #     autocor = np.correlate(Q_n_t[j, :], Q_0[j, :], mode="full")
    #     #print(autocor)
    #     #exit()
    #     autocor = autocor[autocor.size//2:]
    #     autocor = np.divide(autocor, count)
    #     if j == 0:
    #         autocorrelation = autocor
    #     else:
    #         autocorrelation += autocor
    #     #autocorrelation_xyz[j] = autocor
    # autocorrelation = autocorrelation / 3.0

    #exit()
    autocor_xyz = np.concatenate((np.arange(autocor_xyz.shape[1]).reshape(1,autocor_xyz.shape[1]), autocor_xyz), axis=0)
    np.savetxt('Q_sum.out', Q_sum.transpose(), fmt='%.2f', delimiter=' ')
    np.savetxt('Q_sum_convection.out', Q_sum_convex, fmt='%.2f', delimiter=' ')
    np.savetxt('Q_sum_autocor.out', autocor_xyz.transpose(), fmt='%3.2f', delimiter=' ')

    #print((energy_array[:, 0]*velocity_ms_array[:,0]/mass_sqr_array).sum())
    # calculate cross-correlation of heat flux
    print("Start calculate cross-correlation of heat flux\n")
    # padding with zeros
    padding = np.zeros((e_vector_conj.shape[0], frame_num))

    #Q_n_t_all_x = np.concatenate((Q_n_t_all_x, padding), axis=1)
    print(Q_n_t_all_x.shape)
    count = np.tile(np.arange(frame_num, 0, -1), (Q_n_t_all_x.shape[0], 1))
    # Q Q cross correlation
    Q_Q_coor = np.zeros((atom_num * 3, atom_num * 3))
    print_few_k = False
    if print_few_k:
        ka=0.0
        kb=0.0
    sum_all_QQ = 0.0
    # debug
    print(Q_n_t_all_x.sum(axis=0)[:10])
    print(Q_sum[0,:10])
    #exit()
    # use fast fourier transform to get the cross correlation along x/y/z directions
    aft = rfft(Q_n_t_all_x, n=frame_num*2, axis=1)
    bft = np.conjugate(aft)
    row_num = Q_n_t_all_x.shape[0]
    for i in range(row_num):
        print("calculate cross correlation ith:", i)
        c = irfft(aft[:row_num-i, ]*bft[i:, ], axis=1)
        c = c[:, :frame_num]
        #if i + 1 == Q_n_t_all_x.shape[0]:
        #    c = c[0:frame_num]
        #else:
        #    c = c[:, :frame_num]
        #print(c, count[i:,])
        c = np.divide(c, count[i:, ])

        # Note in lammps generated correlation file (e.g. profile.heatflux), the correlation is only up to (framenumber-1)
        # so to match kappa, I use [:,:coor_frame_num] or [:,:-1]
        cc_mat = c[:,:coor_frame_num].sum(axis=1) - 0.5 * c[:, 0]

        if i > 0:
            bottom = irfft(aft[i:, ]*bft[:row_num-i, ], axis=1)
            bottom = bottom[:, :frame_num]
            bottom = np.divide(bottom, count[i:, ])
            bot_mat = bottom[:,:coor_frame_num].sum(axis=1) - 0.5 * bottom[:, 0]

        if i == 0:
            sum_all_QQ += cc_mat.sum()
        else:
            sum_all_QQ += cc_mat.sum() + bot_mat.sum()
        if print_few_k:
            if i ==0:
                ka=cc_mat[0]
            if i ==1:
                kb=cc_mat[0]
        if i == 0:
            np.fill_diagonal(Q_Q_coor, cc_mat)
        else:
            np.fill_diagonal(Q_Q_coor[:-i, i:], cc_mat)
    print("Q_Q_coor shape:", Q_Q_coor.shape)
    with open('QQ_coorelation.npy', 'wb') as f:
        np.save(f, w_value[0, :])
        np.save(f, Q_Q_coor)
    np.savetxt('QQ_coorelation.dat', Q_Q_coor, fmt='%.2f', delimiter=' ')
    print("kappa from sum_all_QQ_along_x[This value can be different from calcualtion using"
          " lammps autocorrelation \nbecause only upper QQ correlation matrix is computed and diagonal"
          " is assumed but is not guaranteed] (W/m/k)\n:", sum_all_QQ * time_step_per_frame/volume/kB/T/T*convert,
          " (W/m/k) \n")  # in Unit of W/m/k
    print("convert:", convert)
    if print_few_k:
     print("Q[0,0]:",ka, " Q[0,1]:",kb)
    return
