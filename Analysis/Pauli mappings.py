# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:04:19 2018

@author: Jarnd
"""
import numpy as np
import itertools as itt

n = 2
def get_pauli_basis(n):
    I =    np.array([[1 , 0] , [0 , 1]]) * (1/2)
    X =    np.array([[0 , 1] , [1 , 0]]) * (1/2)
    Y = 1j*np.array([[0 ,-1] , [1 , 0]]) * (1/2)
    Z =    np.array([[1 , 0] , [0 ,-1]]) * (1/2)
    P1 = [I,X,Y,Z]
    P2 = []
    for Bde in itt.product(P1,repeat = n):
        B = 1
        for i in Bde:
            B = np.kron(B,i)
        P2.append(B)
    return P2

def get_canonical_basis(n):
    E00 = np.array([[1 , 0] , [0 , 0]]);
    E01 = np.array([[0 , 1] , [0 , 0]]);
    E10 = np.array([[0 , 0] , [1 , 0]]);
    E11 = np.array([[0 , 0] , [0 , 1]]);
    E1 = [E00,E01,E10,E11]
    E2 = []
    for Bde in itt.product(E1,repeat = n):
        B = 1
        for i in Bde:
            B = np.kron(B,i)
        E2.append(B)
    return E2




#%%
def get_canon_basis(n):
    basis = []
    for i in range(2**n):
        vect =  np.mat(np.zeros((2**n,1)))
        vect[i] = 1
        basis.append(vect)
    return basis

def get_max_ent_2n(n):
    n_bas = get_canon_basis(n)
    d = len(n_bas)
    for i in range(d):
        if i == 0:
            psi_ome = np.kron(n_bas[i],n_bas[i])
        else:
            psi_ome += np.kron(n_bas[i],n_bas[i])
    return psi_ome/np.sqrt(d)

#%%
def get_choi_basis(n, chi_bas):
    B_choi = []
    I =    np.array([[1 , 0] , [0 , 1]])
    psi_ome = get_max_ent_2n(n)
    for P in chi_bas:
        B_choi.append(np.kron(P,np.kron(I,I))@psi_ome)
    return B_choi

#%%
def unitary_to_choi(U):
    vect = np.mat(np.reshape(U,(-1,1)))
    return vect @ vect.T

def choi_to_chi(choi,B_choi,n):
    chi = np.zeros(np.shape(choi),dtype='complex')
    for combi in itt.product(range((2*n)**2), repeat = 2):
        chi[combi] = complex(B_choi[combi[0]].T @ choi @ B_choi[combi[1]])
    return chi

def chi_to_choi(chi,B_choi):
    choi = np.zeros(np.shape(chi),dtype='complex');
    for combi in itt.product(range((2*n)**2), repeat = 2):
        choi += chi[combi] * B_choi[combi[0]] @ B_choi[combi[1]].T
    return choi

#%%
def get_pauli_list(n,P1names):
    Pauli_list = []
    for i in range((2*n)):
        for j in range((2*n)):
            Pauli_list.append([P1names[i],P1names[j]])
    return Pauli_list

def get_string_from_prepmeas(prep_paulis,eigP0,eigP1,meas_paulis):
    prep = prep_paulis[1]+str(eigP1)+'(1)'+prep_paulis[0]+str(eigP0)+'(0)'
    meas = meas_paulis[1]+'(1)'+meas_paulis[0]+'(0)'
    label = '_prep_'+prep+'_meas_'+meas
    return label

    
def get_A_mat(B_prep,B_meas,B_chi):
    ic = 0;
    jc = 0;
    mc = 0;
    nc = 0;
    jctot = len(B_meas);
    nctot = len(B_chi);
    A = np.zeros((len(B_prep)*len(B_meas),len(B_chi)**2),dtype='complex')
    for i in B_prep:
        jc = 0;
        for j in B_meas:
            mc = 0;
            for m in B_chi:
                nc = 0;
                for n in B_chi:
                    A[jc+(ic*jctot),nc+(mc*nctot)] = np.trace(j@m@i@np.mat(n).H)
                    nc += 1
                mc += 1
            jc += 1
        ic += 1
    return A

def get_lambda_from_meas(meas_data, tomo_set,n):
    P1names = ['I','X','Y','Z']
    reworked_meas_data = rework_data(tomo_set,meas_data)
    lam = [];
    prep = get_pauli_list(n,P1names);
    meas = get_pauli_list(n,P1names);
    for Pi in prep:
        for Pj in meas:
            lam.append(get_lamij_from_PiPj(Pi,Pj,reworked_meas_data))
    return lam

def rework_data(tomo_set,meas_data):
    reworked_meas_data = dict()
    labels = tomo_set['circuit_labels']
    for ind in range(len(labels)):
        for j in [0,1]:
            for i in [0,1]:
                circuitlabel = labels[ind]
                data = meas_data[ind]['counts'][str(j)+str(i)]
                updatedlabel = circuitlabel[:23]+str(j)+circuitlabel[23:27]+str(i)+circuitlabel[27:]
                reworked_meas_data[updatedlabel] = data
    return reworked_meas_data

def get_lamij_from_PiPj(Pi,Pj,reworked_meas_data):
    lookupP = {{'X' : [1, -1, 0, 0, 0, 0] , 'Y' : [0, 0, 1, -1, 0, 0] , 'Z' : [0, 0, 0, 0, 1, -1] , 'I' : [1, 1, 1, 1, 1, 1]}}
    namesP = ['X0','X1','Y0','Y1','Z0','Z1'];
    Pi0list = lookupP[Pi[0]] ;
    Pi1list = lookupP[Pi[1]] ;
    Pj0list = lookupP[Pj[0]] ;
    Pj1list = lookupP[Pj[1]] ;
    lamij = 0;
    for i0 in range(5):
        for i1 in range(5):
            for j0 in range(5):
                for j1 in range(5):
                    lam = Pi0list[i0]*Pi1list[i1]*Pj0list[j0]*Pj1list[j1]
                    circuitlabel = '_prep_'+namesP[i1]+'(1)'+namesP[i0]+'(0)_meas_'+namesP[j1]+'(1)'+namesP[j0]+'(0)';
                    meas = reworked_meas_data[circuitlabel];
                    lamij += lam*meas;
    return lamij














