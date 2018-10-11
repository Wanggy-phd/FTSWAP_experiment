# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:04:19 2018

@author: Jarnd
"""
from functools import reduce
import itertools as itt
import numpy as np
import operator as op

import Analysis.Paulifunctions as pf

#%%

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

"""
def get_canonical_basis(n):
    E00 = np.mat([[1, 0], [0, 0]])
    E01 = np.mat([[0, 1], [0, 0]])
    E10 = np.mat([[0, 0], [1, 0]])
    E11 = np.mat([[0, 0], [0, 1]])
    E1 = [E00, E01, E10, E11]
    E2 = []
    for Bde in itt.product(E1, repeat=n):
        B = 1
        for i in Bde:
            B = np.kron(B, i)
        E2.append(B)
    return E2
"""

#%%
def get_canon_basis(n):
    basis = []
    for i in range(2**n):
        vect = np.mat(np.zeros((2**n, 1)))
        vect[i] = 1
        basis.append(vect)
    return basis


def get_max_ent_2n(n):
    n_bas = get_canon_basis(n)
    d = len(n_bas)
    for i in range(d):
        if i == 0:
            psi_ome = np.kron(n_bas[i], n_bas[i])
        else:
            psi_ome += np.kron(n_bas[i], n_bas[i])
    return psi_ome/np.sqrt(d)

#%%

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


#%% Mapping functions

def unitary_to_choi(U):
    vect = np.mat(np.reshape(U, (-1, 1)))
    return vect @ vect.H


def choi_to_chi(choi, B_choi, n):
    choi = np.mat(choi)
    chi = np.zeros(np.shape(choi), dtype='complex')
    for combi in itt.product(range(4**n), repeat=2):
        chi[combi] = complex(B_choi[combi[0]].H @ choi @ B_choi[combi[1]])
    chi /= np.trace(chi)
    return 2**n*chi # ensure chi matrix is trace-D


def chi_to_choi(chi, B_choi, n):
    choi = np.zeros(np.shape(chi), dtype='complex')
    for combi in itt.product(range(4**n), repeat=2):
        choi += chi[combi] * B_choi[combi[0]] @ B_choi[combi[1]].H
    return choi / np.trace(choi) # ensure choi matrix is trace-1

#%%


def get_pauli_list(n, P1names):
    return list(itt.product(P1names, repeat=n))


def get_pauli_names(n):
    p1names = ['I', 'X', 'Y', 'Z']
    pnames = []
    for p in itt.product(p1names, repeat=n):
        pnames.append(reduce(op.add, p))
    return pnames


def get_string_from_prepmeas(prep_paulis, eigP0, eigP1, meas_paulis):
    prep = prep_paulis[1] + str(eigP1) + '(1)' + prep_paulis[0] + str(eigP0) + '(0)'
    meas = meas_paulis[1] + '(1)' + meas_paulis[0] + '(0)'
    label = '_prep_' + prep + '_meas_' + meas
    return label


def get_A_mat_fast(nq):
    indices = [0, 1, 2, 3]
    d4 = (2**nq)**2
    A = np.zeros((d4**2, d4**2), dtype='complex')
    i = 0
    for Bi in itt.product(indices, repeat=nq):
        j = 0
        for Bj in itt.product(indices, repeat=nq):
            m = 0
            for Bm in itt.product(indices, repeat=nq):
                n = 0
                for Bn in itt.product(indices, repeat=nq):
                    A[j+i*d4, n+m*d4] = pf.calc_trace_P2prod([Bj, Bm, Bi, Bn])
                    n += 1
                m += 1
            j += 1
        i += 1
    return A


def get_A_mat_faster(nq):
    indices = [0, 1, 2, 3]
    d = 2**nq
    row = np.empty((d**8, 1), dtype='complex')
    i = 0
    for Bs in itt.product(indices, repeat=nq*4):
        row[i] = pf.calc_trace_P2prod([Bs[2:4], Bs[4:6], Bs[0:2], Bs[6:8]])
        i += 1
    A = np.reshape(row, ((d)**4, (d)**4))
    return A


def get_A_mat(B_prep, B_meas, B_chi):
    jctot = len(B_meas)
    nctot = len(B_chi)
    A = np.zeros((len(B_prep)*len(B_meas), len(B_chi)**2), dtype='complex')
    for jc, j in enumerate(B_meas):
        for mc, m in enumerate(B_chi):
            jm = j @ m
            for ic, i in enumerate(B_prep):
                jmi = jm @ i
                for nc, n in enumerate(B_chi):
                    A[jc + (ic * jctot), nc + (mc * nctot)] = np.trace(jmi @ np.mat(n).H)
    return A
#%%


def get_lambda_from_meas(tomo_set, meas_data, n):
    P1names = ['I', 'X', 'Y', 'Z']
    reworked_meas_data = rework_data(tomo_set, meas_data)
    lam = []
    lampau = []
    prep = get_pauli_list(n, P1names)
    meas = get_pauli_list(n, P1names)
    for Pi in prep:
        for Pj in meas:
            lam.append(get_lamij_from_PiPj(Pi, Pj, reworked_meas_data, n)[0])
            lampau.append([Pi, Pj])
    return lam, lampau


def rework_data(tomo_set, meas_data):
    reworked_meas_data = dict()
    labels = tomo_set['circuit_labels']
    for ind in range(len(labels)):
        for j in [0, 1]:
            for i in [0, 1]:
                circuitlabel = labels[ind]
                data = (meas_data[ind]['counts']
                        [str(j)+str(i)])/meas_data[ind]['shots']
                updatedlabel = circuitlabel[:23] + \
                    str(j)+circuitlabel[23:27]+str(i)+circuitlabel[27:]
                reworked_meas_data[updatedlabel] = data
    return reworked_meas_data


def get_lamij_from_PiPj(Pi, Pj, reworked_meas_data, n):
    lookupP = {'X': [1, -1, 0, 0, 0, 0], 'Y': [0, 0, 1, -1, 0, 0],
               'Z': [0, 0, 0, 0, 1, -1], 'I': [1/3, 1/3, 1/3, 1/3, 1/3, 1/3]}
    namesP = ['X0', 'X1', 'Y0', 'Y1', 'Z0', 'Z1']
    Pi0list = lookupP[Pi[0]]
    Pi1list = lookupP[Pi[1]]
    Pj0list = lookupP[Pj[0]]
    Pj1list = lookupP[Pj[1]]
    lamij = 0
    lams = []
    circuitlabels = []
    for i0 in range(6):
        for i1 in range(6):
            for j0 in range(6):
                for j1 in range(6):
                    lam = Pi0list[i0]*Pi1list[i1]*Pj0list[j0]*Pj1list[j1]
                    circuitlabel = '_prep_' + \
                        namesP[i1]+'(1)'+namesP[i0]+'(0)_meas_' + \
                        namesP[j1]+'(1)'+namesP[j0]+'(0)'
                    if lam != 0:
                        circuitlabels.append([circuitlabel, lam])
                        meas = reworked_meas_data[circuitlabel]
                        lams.append(lam)
                        lamij += lam*meas
    return lamij, lams, circuitlabels
