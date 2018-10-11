# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:35:42 2018

@author: Jarnd
"""
import sympy as smp
import numpy as np
import itertools as itt
from functools import reduce

ps, pm = smp.symbols('ps pm')

n=1
sq2 = np.sqrt(2)

I = np.mat([[1, 0], [0, 1]],dtype='complex')
X = np.mat([[0, 1], [1, 0]],dtype='complex')
Y = 1j*np.mat([[0, -1], [1, 0]],dtype='complex')
Z = np.mat([[1, 0], [0, -1]],dtype='complex')
H = smp.sqrt(smp.Rational(1,2))*(X+Z)
HI = np.kron(H,I)
IH = np.kron(I,H)
HH = HI * IH
paulis_1 = [I, X, Y, Z]
SWAP = np.mat([[1, 0, 0, 0],[0, 0, 1, 0],[0, 1, 0, 0],[0, 0, 0, 1]])



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

def get_decomp(U,bas,n):
    """
    Get the decomposition of unitary U into a basis specified by bas
    Returns a string of length 2^(2n) with the ith entry the weight of the ith matrix in bas
    """
    weights = np.zeros(len(bas),dtype = 'complex')
    for i,mat in enumerate(bas):
        weights[i] = np.trace(mat.conj().T*U)/np.trace(mat.conj().T*mat)
    return np.mat(weights)

def get_paulidecomp(U,n,normalisation = False):
    """
    Get the decomposition of unitary U into the pauli basis.
    Returns a string of length 2^(2n) with the ith entry the weight of the ith pauli
    """
    P = get_pauli_basis(n, normalise = normalisation)
    weights = np.zeros(len(P),dtype = 'complex')
    for i,pauli in enumerate(P):
        weights[i] = np.trace(pauli*U)/np.trace(pauli*pauli)
    return np.mat(weights)

def U_to_chi(U,n, normalisation):
    weights = get_paulidecomp(U,n, normalisation)
    return np.matmul(weights.conj().T, weights)

def get_chi_dep(p,n, normalized=False):
    d4 = 2**(2*n)
    supp = np.ones(d4,dtype='object')*p*smp.Rational(1,d4-1)
    supp[0] = 1-p
    return np.mat(np.diag(supp))*(2**n) if normalized else np.mat(np.diag(supp))

def left_ext_chi(chi,n,bas=None, normalise = False):
    if bas == None:
       bas = get_pauli_basis(n,normalise)
    I = np.eye(2**n)
    weights = get_decomp(I,bas,n)
    chi_I = np.matmul(weights.T, weights)
    return np.kron(chi_I,chi)

def right_ext_chi(chi,n,bas=None,normalise = False):
    if bas == None:
        bas = get_pauli_basis(n,normalise)
    I = np.eye(2**n)
    weights = get_decomp(I,bas,n)
    chi_I = np.matmul(weights.T, weights)
    return np.kron(chi,chi_I)

def get_supersuper_from_chi(chi,bas,n):
    SS = np.zeros((2**(2*n),2**(2*n)),dtype = 'object')
    for Bm, Bn in itt.product(enumerate(bas),repeat = 2):
        SS += np.kron(chi[Bm[0],Bn[0]],np.kron(Bn[1].conj(),Bm[1]))
    return SS

def get_chi_from_supersuper(SS,bas,n):
    chi = np.zeros((2**(2*n),2**(2*n)),dtype = 'object')
    for Bm, Bn in itt.product(enumerate(bas),repeat = 2):
        basvect = np.kron(Bn[1].T,Bm[1].conj().T)
        print(np.shape(chi))
        #print(basvect)
        print(Bm[0],Bn[0])
        chi[Bm[0],Bn[0]] = np.trace(basvect*SS)/np.trace(basvect*basvect.conj().T)
    return chi*(2**n)/(np.trace(chi)*np.trace(np.matmul(bas[0].conj().T,bas[0])))

def roundsym(expr):
    expr2 = expr
    for a in smp.preorder_traversal(expr):
        if isinstance(a, smp.Float):
            expr2 = expr2.subs(a, round(a, 1))
    return expr2

def get_chi_from_chis(chi_list, bas, n, normalisation):
    S = get_supersuper_from_chi(chi_list[0],bas,n)
    for chi in chi_list[1:]:
        S *= get_supersuper_from_chi(chi,bas,n)
    return get_chi_from_supersuper(S, bas, n)

normalisation = False

P1 = get_pauli_basis(1,normalise = normalisation)
P2 = get_pauli_basis(2,normalise = normalisation)
P3 = get_pauli_basis(3,normalise = normalisation)

c1_H = U_to_chi(H,1, normalisation)
c2_HI = U_to_chi(HI,2, normalisation)
c2_IH = U_to_chi(IH,2, normalisation)
c3_HII = right_ext_chi(c2_HI,2)
c3_IHI = right_ext_chi(c2_IH,2)
c3_IIH = left_ext_chi(c2_IH,2)


#%%
S2_HI = get_supersuper_from_chi(c2_HI, P2,2)
S2_IH = get_supersuper_from_chi(c2_IH, P2,2)
S2_HH = S2_HI * S2_IH

c2_HH_mul = get_chi_from_supersuper(S2_HH, P2, 2)
c2_HH = U_to_chi(HH,2,normalisation)
##%%
#S3_HII = get_supersuper_from_chi(c3_HII,P3, 3)
#S3_IHI = get_supersuper_from_chi(c3_IHI,P3, 3)
#S3_IIH = get_supersuper_from_chi(c3_IIH,P3, 3)
#
#S3_HHH = S3_HII * S3_IHI * S3_IIH
#
#c3_HHH = get_chi_from_supersuper(S3_HHH, P3, 3)
#
#HHH = np.kron(HH,H)
#c3_HHH_dir = U_to_chi(HHH,3,normalisation)
#
#c1_dep = get_chi_dep(ps,1, normalisation)
#c2_depI = right_ext_chi(c1_dep,1, normalise=normalisation)
#c2_Idep = left_ext_chi(c1_dep,1, normalise=normalisation)
#c3_depII = right_ext_chi(c2_depI,2, normalise=normalisation)
#c3_IdepI = left_ext_chi(c2_depI,2, normalise=normalisation)
#c3_IIdep = left_ext_chi(c2_Idep,2, normalise=normalisation)

