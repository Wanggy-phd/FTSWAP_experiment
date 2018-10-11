# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:22:48 2018

@author: Jarnd
"""
import qiskit.tools.qcvv.tomography as tomo
import Analysis.tomography_functions as tomoself
import numpy as np


def fit_tomodata(tomo_data, method=None):
    if method == 'Linear inversion':
        print('Error, no method yet!')
        choi_fit = []
    elif method == None:
        choi_fit = tomo.fit_tomography_data(tomo_data, options={'trace': 1})
    else:
        choi_fit = tomo.fit_tomography_data(
            tomo_data, method, options={'trace': 1})
    return choi_fit


def fit_chi_own(tomo_data, tomo_set, n):
    B_chi = tomoself.get_pauli_basis(n)
    B_prep = tomoself.get_pauli_basis(n, normalise=False)
    B_meas = tomoself.get_pauli_basis(n, normalise=False)

    lam, lampau = tomoself.get_lambda_from_meas(tomo_set, tomo_data['data'], n)

    A = tomoself.get_A_mat(B_prep, B_meas, B_chi)
    chivect = np.linalg.solve(A, lam)
    return np.reshape(chivect, ((2*n)**2, (2*n)**2))


def get_total_prob(tomo_data):
    meas_data = tomo_data['data']
    assert type(meas_data) == list
    counts = []
    for meas in meas_data:
        countsvalues = meas['counts'].values()
        counts.append(sum(list(countsvalues))/meas['shots'])
    return counts


def make_CP(chi, n):
    assert np.shape(chi) == ((2**n)**2, (2**n)**2)
    mineig = np.min(np.linalg.eigvals(chi))
    trace_chi = np.trace(chi)
    if mineig < 0:
        chiCP = np.add(chi, (-1)*mineig*np.eye((2*n)**2))
    return trace_chi*chiCP/np.trace(chiCP)


def get_TPsum(chi, B_chi):
    d2 = np.shape(chi)[0]
    iden = np.zeros_like(B_chi[0], dtype='complex')
    for m in range(d2):
        for n in range(d2):
            iden += chi[m, n]*np.mat(B_chi[n]).H@np.mat(B_chi[m])
    return iden


def check_TP(chi, B_chi, n):
    assert np.shape(chi) == ((2**n)**2, (2**n)**2)
    d = 2**n
    iden = np.eye(d, dtype='complex')
    TPsum = get_TPsum(chi, B_chi)
    diff = iden - TPsum
    if np.around(np.sum(diff), 0) == 0:
        return True
    else:
        return False

#%% Errors


def get_chi_error(chi, chi_bas, U, mode='p'):
    chi = np.mat(chi)
    U = np.mat(U)
    V = np.mat(np.zeros((len(chi_bas), len(chi_bas))), dtype='complex')
    mc = 0
    for i in range(len(chi_bas)):
        chi_bas[i] = np.mat(chi_bas[i])
    for m in chi_bas:
        nc = 0
        for n in chi_bas:
            if mode == 'p':
                V[mc, nc] = np.trace(m.H @ n @ U.H)
            if mode == 'n':
                V[mc, nc] = np.trace(m.H @ U.H @ n)
            nc += 1
        mc += 1
    return V @ chi @ V.H


def filter_chi_meas(chi, B_filter, nq):
    d2 = (2**nq)**2
    chivect = np.reshape(chi, (-1, 1))
    chifilter = np.linalg.solve(B_filter, chivect)
    return np.reshape(chifilter, (d2, d2))


def unit_to_choi(Unitary):
    vect = np.reshape(Unitary, (1, -1))
    return np.mat(vect).H @ np.mat(vect)


def get_measfiltermatind(chi_meas, nq):
    indices = [0, 1, 2, 3]
    d = 2**nq
    row = np.empty((d**8, 1), dtype='complex')
    i = 0
    chi_meas_row = np.reshape(chi_meas, (-1, 1))
    for Bs in tomoself.itt.product(indices, repeat=nq*4):
        ij = 0
        for Bij in tomoself.itt.product(indices, repeat=nq*2):
            tr1 = tomoself.pf.calc_trace_P2prod([Bs[0:2], Bij[0:2], Bs[4:6]])
            tr2 = tomoself.pf.calc_trace_P2prod([Bs[2:4], Bij[2:4], Bs[6:8]])
            row[i] += chi_meas_row[ij]*tr1*tr2
            ij += 1
        i += 1
        print('i is:', i)
    B = np.reshape(row, ((d)**4, (d)**4))
    return B


def get_measfiltermat(chi_meas, B_chi, nq):
    mc = 0
    nc = 0
    kc = 0
    lc = 0
    nctot = len(B_chi)
    lctot = len(B_chi)
    B = np.zeros((len(B_chi)*len(B_chi), len(B_chi)**2), dtype='complex')
    for m in B_chi:
        nc = 0
        for n in B_chi:
            kc = 0
            for k in B_chi:
                lc = 0
                for l in B_chi:
                    ic = 0
                    Bmnkl = 0
                    for i in B_chi:
                        jc = 0
                        for j in B_chi:
                            tr1 = np.trace(m@i@k)
                            tr2 = np.trace(n@l@j)
                            Bmnkl += chi_meas[ic, jc] * tr1 * tr2
                            jc += 1
                        ic += 1
                    B[nc+(mc*nctot), lc+(kc*lctot)] = Bmnkl
                    lc += 1
                kc += 1
            nc += 1
        mc += 1
    return B


def filter_meas(chi, chi_meas, n):
    k = n
    return k

def process_fidelity(chi_error, n):
    return chi_error[0,0]/(2**n);

def channel_fidelity(chi_error,B_choi, n):
    choi_error = tomoself.chi_to_choi(chi_error, B_choi, n)
    maxent = tomoself.get_max_ent_2n(n)
    return maxent.H @ choi_error @ maxent