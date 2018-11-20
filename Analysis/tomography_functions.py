# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:04:19 2018

@author: Jarnd
"""
from functools import reduce
import itertools as itt
import numpy as np
import operator as op
import math
#%% Basis generating functions
sq2 = np.sqrt(2)

I = np.mat([[1, 0], [0, 1]])
X = np.mat([[0, 1], [1, 0]])
Y = 1j*np.mat([[0, -1], [1, 0]])
Z = np.mat([[1, 0], [0, -1]])

paulis_1 = [I, X, Y, Z]

def get_pauli_basis(n, normalise=True):
    """
    Lists the Pauli matrices that form a basis for the space of
    operators on `n` qubits.
    """
    basis = [_ / sq2 for _ in paulis_1] if normalise else paulis_1
    
    P2 = []
    for Bde in itt.product(basis, repeat=n):
        P2.append(reduce(np.kron, Bde))
    return P2

def get_canon_basis(n):
    '''
    Get the canonical basis for a n-qubit system.
    Returns a list of length(2^n) with the ith element a zeros vector of length 2^n and the ith element a 1.
    '''
    basis = []
    for i in range(2**n):
        vect = np.mat(np.zeros((2**n, 1)))
        vect[i] = 1
        basis.append(vect)
    return basis


def get_max_ent_2n(n):
    '''
    Get the maximally entangled state on a system for 2n qubits.
    Returns a vector that is sum(i) 1/(2^n) |i> tensor |i>,
    with {|i>} the canonical basis for a n-qubit system,
    and the tensor product the kronecker tensor product.
    '''
    n_bas = get_canon_basis(n)
    d = len(n_bas)
    for i in range(d):
        if i == 0:
            psi_ome = np.kron(n_bas[i], n_bas[i])
        else:
            psi_ome += np.kron(n_bas[i], n_bas[i])
    return psi_ome/np.sqrt(d)

def get_choi_basis(n, chi_bas):
    """
    Note that you have to put in unitary operators in chi_bas in order
    to get normalised states out in the basis.
    """
    B_choi = []
    # I = np.array([[1, 0], [0, 1]])
    assert chi_bas[0].shape == (2**n, 2**n)
    psi_ome = get_max_ent_2n(n)
    for P in chi_bas:
        B_choi.append(np.kron(P, np.eye(2**n)) @ psi_ome)
    return B_choi

def get_pauli_list(n):
    p1names = ['I', 'X', 'Y', 'Z']
    '''
    Returns a list of tuples of the names of all n-qubit tensor products of Paulis,
    starting at I and ending at Z.
    The last pauli in the name is looped through first. e.g. for 2 qubits:
        The function returns [('I','I'),('I','X'),('I','Y'),('I','Z'),('X','I').....('Z','X'),('Z','Y'),('Z','Z')]
    '''
    return list(itt.product(p1names, repeat=n))


def get_pauli_names(n):
    '''
    Returns a list of the names of all n-qubit tensor products of Paulis,
    starting at I and ending at Z.
    The last pauli in the name is looped through first. e.g. for 2 qubits:
        The function returns [II,IX,IY,IZ,XI.....ZX,ZY,ZZ]
    '''
    p1names = ['I', 'X', 'Y', 'Z']
    pnames = []
    for p in itt.product(p1names, repeat=n):
        pnames.append(reduce(op.add, p))
    return pnames

#%% Mapping functions
def choi_to_chi(choi, B_choi, n):
    '''
    Maps the choi matrix to the chi matrix using the relation:
        chi(m,n) = <B_m|rho_choi|B_n>
    With the basis |B_m> = (B_m tensor I) * |Om>,
    and |Om> the maximally entangled state on 2n qubits, gotten from get_max_ent_2n()
    '''
    choi = np.mat(choi)
    chi = np.zeros(np.shape(choi), dtype='complex')
    for combi in itt.product(range(4**n), repeat=2):
        chi[combi] = complex(B_choi[combi[0]].H @ choi @ B_choi[combi[1]])
    chi /= np.trace(chi)                                                        # This allows for a non-normalized basis B_m
    return 2**n*chi                                                             # ensure chi matrix has trace d


def chi_to_choi(chi, B_choi, n):
    '''
    Maps the chi matrix to the choi matrix using the relation:
        choi = sum(m,n) chi(m,n)|B_m><B_n|
    With the basis |B_m> = (B_m tensor I) * |Om>,
    and |Om> the maximally entangled state on 2n qubits, gotten from get_max_ent_2n()
    '''
    choi = np.zeros(np.shape(chi), dtype='complex')
    for combi in itt.product(range(4**n), repeat=2):
        choi += chi[combi] * B_choi[combi[0]] @ B_choi[combi[1]].H
    return choi / np.trace(choi)                                                # ensure choi matrix has trace 1

#%% Helper functions for analysefunc.fit_chi_own()
def get_A_mat(B_prep, B_meas, B_chi):
    '''
    Calculate the A matrix that links the vectorized versions of the measurement vector lambda and the process matrix chi together.
    The matrix is calculated elementwise, with A(ij,kl):
        A(ij,mn) = tr(P_j P_m P_i P_n^(dagger)).
        
        P_j are the elements of B_meas
            These are the bases in which the measurements of lambda are performed
        P_i are the elements of B_prep
            These are the bases into which the input states of the measurements are reworked.
            See get_lambda_from_meas() for more details. That function only works for Pauli input states.
        P_m and P_n are the elements of B_chi
            This is the basis in which chi will be described. The Pauli basis is often chosen
        
        Pjm = P_j @ P_m is calculated at the start of every m loop,
        and then Pjmi = Pjm @ P_i is calculated at the start of every i loop.
        Finally, at every n loop tr(Pjmi @ P_n^(dagger)) is calculated.
        
        The indices of A are like:
        row    ij = j + (i*j_total)
        column mn = n + (m*n_total)
    '''
    jctot = len(B_meas)
    nctot = len(B_chi)
    A = np.zeros((len(B_prep)*len(B_meas), len(B_chi)**2), dtype='complex')
    for jc, j in enumerate(B_meas):
        for mc, m in enumerate(B_chi):
            jm = j @ m
            for ic, i in enumerate(B_prep):
                jmi = jm @ i
                for nc, n in enumerate(B_chi):
                    A[jc + (ic * jctot), nc + (mc * nctot)] = np.trace(jmi @ np.mat(n).H) # np.mat has a .H property that is the hermitian
    return A

def get_lambda_from_meas(tomo_set, meas_data, n):
    '''
    Calculates the lambda vector that relates to the chi matrix via the A matrix.
    The lambda vector is calculated from the measurement data and the specified tomography set.
    Lambda is a vector of length i_tot* j_tot, with i_tot and j_tot the number of different input and output bases respectively.
    Lambda is of the form:
        l(ij) = tr(P_j @ L(P_i)), for a input basis P_i, a measurement basis P_j and the system L() on which the tomography is performed.
    
    However, the measurement data is not of that form.
    The measurement data is first reworked using rework_data(tomo_set, meas_data) (see docstring of that function for details)
    The reworked measurement data is then a dictionary of all counts of the type:
        '_prep_X1(1)Y0(0)_meas_Y1(1)Z0(0)' : counts/#shots for
                                preperation of qubit(1) in the -1 eigenstate of X
                                preparation of qubit(0) in the +1 eigenstate of Y
                                measurement of qubit(1) in the Y basis and getting the -1 eigenvalue
                                measurement of qubit(0) in the Z basis and getting the +1 eigenvalue
    
    Then, for every element l(ij) = tr(P_j @ L(P_i)) is expanded into all these preparation states and measurement states.
    This is done by the function get_lamij_from_PiPj(). See that function for details.
    
    Furthermore, the function also calculates the standard deviation on lambda, lambdastddv.
    This is done for every P_i and P_j using the get_lamijstddv_from_PiPj() function.
    See that function for details.
    The function then returns:
        lambda, lampau, lambda_stddv
        with lampau a list with the ith element the pauli names corresponding to the ith element of lambda
    '''
    reworked_meas_data,  rew_data_stddv =  _rework_data_(tomo_set, meas_data)
    lam = []
    lampau = []
    lamstddv = []
    prep = get_pauli_list(n)
    meas = get_pauli_list(n)
    for Pi in prep:
        for Pj in meas:
            lam.append(_get_lamij_from_PiPj_(Pi, Pj, reworked_meas_data, n)) # Get the lamij corresponding to Pi&Pj 
            lampau.append([Pi, Pj])
            lamstddv.append(_get_lamijstddv_from_PiPj_(Pi, Pj, rew_data_stddv, n)) # Get the lamij_stddv corresponding to Pi&Pj
    return np.array(lam), lampau, np.array(lamstddv)

def get_lambda_from_meas_multiple(tomo_set, meas_data_all, n):
    '''
    Calculates the lambda vector that relates to the chi matrix via the A matrix for multiple sets of measurement data, assumed all to have the same number of shots.
    The lambda vector is calculated from the measurement data and the specified tomography set.
    Lambda is a vector of length i_tot* j_tot, with i_tot and j_tot the number of different input and output bases respectively.
    Lambda is of the form:
        l(ij) = tr(P_j @ L(P_i)), for a input basis P_i, a measurement basis P_j and the system L() on which the tomography is performed.
    
    However, the measurement data is not of that form.
    The measurement data is first reworked using rework_data(tomo_set, meas_data) (see docstring of that function for details)
    The reworked measurement data is then a dictionary of all counts of the type:
        '_prep_X1(1)Y0(0)_meas_Y1(1)Z0(0)' : counts/#shots for
                                preperation of qubit(1) in the -1 eigenstate of X
                                preparation of qubit(0) in the +1 eigenstate of Y
                                measurement of qubit(1) in the Y basis and getting the -1 eigenvalue
                                measurement of qubit(0) in the Z basis and getting the +1 eigenvalue
    
    Then, for every element l(ij) = tr(P_j @ L(P_i)) is expanded into all these preparation states and measurement states.
    This is done by the function get_lamij_from_PiPj(). See that function for details.
    
    Furthermore, the function also calculates the standard deviation on lambda, lambdastddv.
    This is done for every P_i and P_j using the get_lamijstddv_from_PiPj() function.
    See that function for details.
    The function then returns:
        lambda, lampau, lambda_stddv
        with lampau a list with the ith element the pauli names corresponding to the ith element of lambda
    '''
    nr_jobs = len(meas_data_all)
    reworked_meas_data_all_list = []
    reworked_meas_data_all = dict()
    reworked_meas_data_all_stddv = dict()
    for meas_data in meas_data_all:
        reworked_meas_data,  rew_data_stddv =  _rework_data_(tomo_set, meas_data)
        reworked_meas_data_all_list.append(reworked_meas_data)
    for jobnr, job in enumerate(reworked_meas_data_all_list):
        for key in job.keys():
            if jobnr == 0:
                reworked_meas_data_all[key] = reworked_meas_data[key]
            else:
                reworked_meas_data_all[key] = reworked_meas_data_all[key]+ reworked_meas_data[key]
    for key in reworked_meas_data_all.keys():
        reworked_meas_data_all[key] = reworked_meas_data_all[key]/nr_jobs
    for key in reworked_meas_data_all.keys():
        p = reworked_meas_data_all[key]
        reworked_meas_data_all_stddv[key] = np.sqrt(p*(1-p)/(nr_jobs*8192))
    lam = []
    lampau = []
    lamstddv = []
    prep = get_pauli_list(n)
    meas = get_pauli_list(n)
    for Pi in prep:
        for Pj in meas:
            lam.append(_get_lamij_from_PiPj_(Pi, Pj, reworked_meas_data_all, n)) # Get the lamij corresponding to Pi&Pj 
            lampau.append([Pi, Pj])
            lamstddv.append(_get_lamijstddv_from_PiPj_(Pi, Pj, reworked_meas_data_all_stddv, n)) # Get the lamij_stddv corresponding to Pi&Pj
    return np.array(lam), lampau, np.array(lamstddv)


def _rework_data_(tomo_set, meas_data):
    '''
    Reworks the data of meas_data into a reworked format that works with _get_lamij_from_PiPj()
    The measurement data provided by the IBM Q experience is of the form:
        '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)' : dict for 
                                preparation of qubit(1) in the -1 eigenstate of Y
                                preparation of qubit(0) in the +1 eigenstate of Z
                                measurement of qubit(1) in the Z basis
                                measurement of qubit(0) in the X basis
        Then dict is a dictionary containing elements 'counts', 'shots' and more.
        For 2 qubits, counts is like:
                                '00' : counts for measuring qubit(1) in +1 eigenstate and qubit(0) in +1 eigenstate
                                '01' : counts for measuring qubit(1) in +1 eigenstate and qubit(0) in -1 eigenstate
                                '10' : counts for measuring qubit(1) in -1 eigenstate and qubit(0) in +1 eigenstate
                                '11' : counts for measuring qubit(1) in -1 eigenstate and qubit(0) in -1 eigenstate
        
        This data is reworked into a dictionary reworked_meas_data of 2^n * the size of meas_data, with {keys : values} of the type:
        '_prep_Y1(1)Z0(0)_meas_Z0(1)X0(0)' : '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)'['counts']['00']/#shots
        '_prep_Y1(1)Z0(0)_meas_Z0(1)X1(0)' : '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)'['counts']['01']/#shots
        '_prep_Y1(1)Z0(0)_meas_Z1(1)X0(0)' : '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)'['counts']['10']/#shots
        '_prep_Y1(1)Z0(0)_meas_Z1(1)X1(0)' : '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)'['counts']['11']/#shots
        
    Furthermore, for every key of type '_prep_Y1(1)Z0(0)_meas_Z0(1)X0(0)' a dictionary with the standard deviation is calculated:
        Reworked_data_stddv : {'_prep_Y1(1)Z0(0)_meas_Z0(1)X0(0)' : sqrt(p*(1-p)/n),
        with p the success rate of the binomial distribution that are the counts. i.e.: p is the value in reworked_meas_data corresponding to the key
    
    The function then returns reworked_meas_data, rew_data_stddv
    '''
    reworked_meas_data = dict()
    rew_data_stddv = dict()
    labels = tomo_set['circuit_labels']                                         # All labels of the kind '_prep_Y1(1)Z0(0)_meas_Z(1)X(0)'
    for ind in range(len(labels)):
        shots = meas_data[ind]['shots']                                         # Get the number of shots associated with experiment 'ind'
        for j in [0, 1]:                                                        # Loop over +1 and -1 eigenstates of qubit 1
            for i in [0, 1]:                                                    # Loop over +1 and -1 eigenstates of qubit 0
                circuitlabel = labels[ind]                                      # Get the circuitlabel of the loop
                data = (meas_data[ind]['counts']                                # This is the total number of counts associated with the measurement of qubit 1 in j-eigenvalue...
                        [str(j)+str(i)])/shots                                  # and qubit 0 in i-eigenvalue, divided by total #shots
                updatedlabel = circuitlabel[:23] + \
                    str(j)+circuitlabel[23:27]+str(i)+circuitlabel[27:]         # The label is updated to the desired format
                reworked_meas_data[updatedlabel] = data;                        # Put the value for updatedlabel in the dictionary
                rew_data_stddv[updatedlabel] = np.sqrt(data*(1-data)/shots)     # The standard deviation on success rate p of a binomial variable is (p*(1-p)/n)^(1/2)
    return reworked_meas_data, rew_data_stddv

def _get_lamij_from_PiPj_(Pi, Pj, reworked_meas_data, n):
    '''
    This function calculates lamij from reworked_meas_data for a specific Pi and Pj,
    with Pi and Pj both a list with 2 elements 'I','X','Y','Z'.
    
    A lamij = tr[PjL(Pi)] is first decomposed into all combinations of all Pi's and Pj's involved.
    This works as follows:
        The function loops through the entire reworked_meas_data,
        and calculates the support in the eigendecomposition of PiPj on the elements.
        The support is calculated like this:
            For P_i = ['I','X'] and PJ = ['Y','Z'] the support is on all elements of reworked_meas_Data that:
                Have a +1 or -1 eigenvalue on qubit(1) prepared in 'I': there are only +1 eigenvalues, all preparations on qubit(1) are counted
                Have a +1 or -1 eigenvalue on qubit(0) prepared in 'X':
                    the +1 eigenvalues elements are of the kind: '_prep_??(1)X0(0)_meas_??(1)??(0)'
                    the -1 eigenvalues elements are of the kind: '_prep_??(1)X1(0)_meas_??(1)??(0)'
                Have a +1 or -1 eigenvalue on qubit(1) measured in 'Y':
                    the +1 eigenvalues elements are of the kind: '_prep_??(1)??(0)_meas_Y0(1)??(0)'
                    the -1 eigenvalues elements are of the kind: '_prep_??(1)??(0)_meas_Y1(1)??(0)'
                Have a +1 or -1 eigenvalue on qubit(0) measured in 'Z':
                    the +1 eigenvalues elements are of the kind: '_prep_??(1)??(0)_meas_??(1)Z0(0)'
                    the -1 eigenvalues elements are of the kind: '_prep_??(1)??(0)_meas_??(1)Z1(0)'
            Then, for all 4 paulis there is looped through 'namesP = ['X0', 'X1', 'Y0', 'Y1', 'Z0', 'Z1']'
            From this an circuitlabel is created of the kind '_prep_??(1)??(0)_meas_??(1)??(0)', for all 6^4 options:
                For every option , the value lam is calculated corresponding to the Pi's and the Pj's:
                    If the '_prep_??(1)??(0)_meas_??(1)??(0)' is a +1 eigenvalue of Pi and Pj then lam = 1
                    If the '_prep_??(1)??(0)_meas_??(1)??(0)' is a -1 eigenvalue of Pi and Pj then lam = -1
                    If the '_prep_??(1)??(0)_meas_??(1)??(0)' is a not an eigenstate of Pi and Pj then lam = 0
                Then, lam*meas is added to lamij, with meas the element of reworked_meas_data[_prep_??(1)??(0)_meas_??(1)??(0)']
                lam is calculated by using lookuptabels {'X': [1, -1, 0, 0, 0, 0], 'Y': [0, 0, 1, -1, 0, 0],
               'Z': [0, 0, 0, 0, 1, -1], 'I': [1/3, 1/3, 1/3, 1/3, 1/3, 1/3]}. Note that Identity is shared over all experiments.
    Then lamij is returned
    '''
    lookupP = {'X': [1, -1, 0, 0, 0, 0], 'Y': [0, 0, 1, -1, 0, 0],
               'Z': [0, 0, 0, 0, 1, -1], 'I': [1/3, 1/3, 1/3, 1/3, 1/3, 1/3]}
    namesP = ['X0', 'X1', 'Y0', 'Y1', 'Z0', 'Z1']
    Pi0list = lookupP[Pi[0]]
    Pi1list = lookupP[Pi[1]]
    Pj0list = lookupP[Pj[0]]
    Pj1list = lookupP[Pj[1]]
    lamij = 0                                                                   # Initialize lamij as 0

    for i0 in range(6):                                                         # Loop through every pauli eigenstate for preparation of qubit 0
        for i1 in range(6):                                                     # Loop through every pauli eigenstate for preparation of qubit 1
            for j0 in range(6):                                                 # Loop through every pauli eigenstate for measurement of qubit 0
                for j1 in range(6):                                             # Loop through every pauli eigenstate for measurement of qubit 1
                    lam = Pi0list[i0]*Pi1list[i1]*Pj0list[j0]*Pj1list[j1]       # Calculate lam from the lookup table and the elements of Pi and Pj
                    circuitlabel = '_prep_' + \
                        namesP[i1]+'(1)'+namesP[i0]+'(0)_meas_' + \
                        namesP[j1]+'(1)'+namesP[j0]+'(0)'                       # Creates a new circuitlable in the form of '_prep_??(1)??(0)_meas_??(1)??(0)'
                    if lam != 0:
                        meas = reworked_meas_data[circuitlabel]
                        lamij += lam*meas;                                      # Add lam*meas to the to be returned lamij
    return lamij

def _get_lamijstddv_from_PiPj_(Pi, Pj, rew_data_stddv, n):
    '''
    This function works much the same as _get_lamij_from_PiPj_():
        Here, the standard deviation on lamij is calculated as the square root of the sum of all (lam^2)*(meas^2).
    Check the docstring and source file of _get_lamij_from_PiPj() for more details.
    '''
    lookupP = {'X': [1, -1, 0, 0, 0, 0], 'Y': [0, 0, 1, -1, 0, 0],
               'Z': [0, 0, 0, 0, 1, -1], 'I': [1/3, 1/3, 1/3, 1/3, 1/3, 1/3]}
    namesP = ['X0', 'X1', 'Y0', 'Y1', 'Z0', 'Z1']
    Pi0list = lookupP[Pi[0]]
    Pi1list = lookupP[Pi[1]]
    Pj0list = lookupP[Pj[0]]
    Pj1list = lookupP[Pj[1]]
    lamijstddv2 = 0
    for i0 in range(6):
        for i1 in range(6):
            for j0 in range(6):
                for j1 in range(6):
                    lam = Pi0list[i0]*Pi1list[i1]*Pj0list[j0]*Pj1list[j1]
                    circuitlabel = '_prep_' + \
                        namesP[i1]+'(1)'+namesP[i0]+'(0)_meas_' + \
                        namesP[j1]+'(1)'+namesP[j0]+'(0)'
                    if lam != 0:
                        meas = rew_data_stddv[circuitlabel]
                        lamijstddv2 += (lam**2)*(meas**2)
    return np.sqrt(lamijstddv2)


###############################################################################
################################# Not used ####################################
###############################################################################
#
#def get_canonical_mat_basis(n):
#    E00 = np.mat([[1, 0], [0, 0]])
#    E01 = np.mat([[0, 1], [0, 0]])
#    E10 = np.mat([[0, 0], [1, 0]])
#    E11 = np.mat([[0, 0], [0, 1]])
#    E1 = [E00, E01, E10, E11]
#    E2 = []
#    for Bde in itt.product(E1, repeat=n):
#        B = 1
#        for i in Bde:
#            B = np.kron(B, i)
#        E2.append(B)
#    return E2
#
#def unitary_to_choi(U):
#    vect = np.mat(np.reshape(U, (-1, 1)))
#    return vect @ vect.H
#
#import Analysis.Paulifunctions as pf
#def get_A_mat_fast(nq):
#    indices = [0, 1, 2, 3]
#    d4 = (2**nq)**2
#    A = np.zeros((d4**2, d4**2), dtype='complex')
#    i = 0
#    for Bi in itt.product(indices, repeat=nq):
#        j = 0
#        for Bj in itt.product(indices, repeat=nq):
#            m = 0
#            for Bm in itt.product(indices, repeat=nq):
#                n = 0
#                for Bn in itt.product(indices, repeat=nq):
#                    A[j+i*d4, n+m*d4] = pf.calc_trace_P2prod([Bj, Bm, Bi, Bn])
#                    n += 1
#                m += 1
#            j += 1
#        i += 1
#    return A
#
#
#def get_A_mat_faster(nq):
#    indices = [0, 1, 2, 3]
#    d = 2**nq
#    row = np.empty((d**8, 1), dtype='complex')
#    i = 0
#    for Bs in itt.product(indices, repeat=nq*4):
#        row[i] = pf.calc_trace_P2prod([Bs[2:4], Bs[4:6], Bs[0:2], Bs[6:8]])
#        i += 1
#    A = np.reshape(row, ((d)**4, (d)**4))
#    return A