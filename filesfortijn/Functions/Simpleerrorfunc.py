# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 10:10:57 2018

@author: Jarnd
"""

"""
Created on Tue Oct  9 16:35:42 2018

@author: Jarnd
"""
import numpy as np
import itertools as itt
from functools import reduce
import operator as op

#%%
sq2 = np.sqrt(2)

I = np.mat([[1, 0], [0, 1]],dtype='complex')
X = np.mat([[0, 1], [1, 0]],dtype='complex')
Y = 1j*np.mat([[0, -1], [1, 0]],dtype='complex')
Z = np.mat([[1, 0], [0, -1]],dtype='complex')
II = np.kron(I,I)

paulis_1 = [I,X,Y,Z]


#%%
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

def get_pauli_names(n):
    p1names = ['I', 'X', 'Y', 'Z']
    pnames = []
    for p in itt.product(p1names, repeat=n):
        pnames.append(reduce(op.add, p))
    return pnames

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
    supp = np.ones(d4,dtype='complex')*p*(1/(d4-1))
    supp[0] = 1-p
    return np.mat(np.diag(supp))*(2**n) if normalized else np.mat(np.diag(supp))

def left_ext_chi(chi,n,bas=None, normalise = False):
    if bas == None:
       bas = get_pauli_basis(1,normalise)
    I = np.eye(2**1)
    weights = get_decomp(I,bas,1)
    chi_I = np.matmul(weights.T, weights)
    return np.kron(chi_I,chi)

def right_ext_chi(chi,n,bas=None,normalise = False):
    if bas == None:
        bas = get_pauli_basis(1,normalise)
    I = np.eye(2**1)
    weights = get_decomp(I,bas,1)
    chi_I = np.matmul(weights.T, weights)
    return np.kron(chi,chi_I)

def get_supersuper_from_chi(chi,bas,n):
    SS = np.zeros((2**(2*n),2**(2*n)),dtype = 'complex')
    for Bm, Bn in itt.product(enumerate(bas),repeat = 2):
        SS += np.kron(chi[Bm[0],Bn[0]],np.kron(Bn[1].conj(),Bm[1]))
    return np.mat(SS)

def get_chi_from_supersuper(SS,bas,n):
    chi = np.zeros((2**(2*n),2**(2*n)),dtype = 'complex')
    for Bm, Bn in itt.product(enumerate(bas),repeat = 2):
        basvect = np.kron(Bn[1].T,Bm[1].conj().T)
        chi[Bm[0],Bn[0]] = np.trace(basvect*SS)/np.trace(basvect*basvect.conj().T)
    return np.mat(chi*(2**n)/(np.trace(chi)*np.trace(np.matmul(bas[0].conj().T,bas[0]))))

def get_chi_from_chis(chi_list, bas, n, normalisation):
    S = get_supersuper_from_chi(chi_list[0],bas,n)
    for chi in chi_list[1:]:
        S = get_supersuper_from_chi(chi,bas,n)*S
    return get_chi_from_supersuper(S, bas, n)

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

def index(indices):
    indices.reverse()
    ind = 0
    for i in range(len(indices)):
        ind += indices[i]*(4**i)
    return ind

def swap_qubits_chi(chi,q1,q2,n):
    Perm = np.zeros_like(chi)
    for indices in itt.product(range(4),repeat = n):
        indices = list(indices)
        indicesnew = indices.copy()
        ind1 = index(indices)
        indicesnew[q2] = indices[q1]
        indicesnew[q1] = indices[q2]
        ind2 = index(indices)
        Perm[ind1,ind2] = 1
        Perm[ind2,ind1] = 1
    return Perm

def get_trace_S(n):
    Tr = np.zeros((2**(2*(n-1)),2**(2*n)),dtype = 'complex')
    Iten = 1
    for i in range(n-1):
        Iten = np.kron(Iten, I)
    for i in [0,1]:
        k = np.zeros((1,2),dtype = 'complex')
        k[0,i] = 1
        bas = np.kron(Iten,k)
        Tr += np.kron(bas,bas)
    return Tr

#def chi_trace_q(chi,n,q):
#    nontrace = list(range(n))
#    nontrace.remove(q)
#    dtraced = 2**(2*(n-1))
#    chi_traced = np.zeros((dtraced,dtraced))
#    for m in itt.product(range(4),repeat = (n-1)):
#        for n in itt.product(range(4), repeat = (n-1)):
#            
def chi3_trace_3(chi):
    chi_traced = np.zeros((16,16),dtype='complex')
    for m in itt.product(range(4),repeat = 2):
        mtraced = 4*m[0]+m[1]
        mold = 16*m[0]+4*m[1]
        for n in itt.product(range(4), repeat = 2):
           ntraced = 4*n[0]+n[1]
           nold = 16*n[0]+4*n[1]
           chi_traced[mtraced,ntraced] = chi[mold,nold]
    return chi_traced/np.trace(chi_traced)
           

#%%
normalisation = False

P1 = get_pauli_basis(1,normalise = normalisation)
P2 = get_pauli_basis(2,normalise = normalisation)
P3 = get_pauli_basis(3,normalise = normalisation)

P2norm = get_pauli_basis(2, normalise = True)
P3norm = get_pauli_basis(3, normalise = True)

H = (1/sq2)*(X+Z)

HI = np.kron(H,I)
IH = np.kron(I,H)
HH = HI * IH

CX = np.mat([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 0, 1],[0, 0, 1, 0]])
SWAP = np.mat([[1, 0, 0, 0],[0, 0, 1, 0],[0, 1, 0, 0],[0, 0, 0, 1]])
SWAPI = np.kron(SWAP,I)

ICX = np.kron(I,CX)
CXI = np.kron(CX,I)

ISWAP = np.kron(I,SWAP)
CIX = CXI * ISWAP


c1_H = U_to_chi(H,1, normalisation)
c2_HI = U_to_chi(HI,2, normalisation)
c2_IH = U_to_chi(IH,2, normalisation)


c2_SWAP = U_to_chi(SWAP, 2, normalisation)
c2_CX = U_to_chi(CX, 2, normalisation)

c3_HII = right_ext_chi(c2_HI,2)
c3_IHI = right_ext_chi(c2_IH,2)
c3_IIH = left_ext_chi(c2_IH,2)
c3_SWAPI = right_ext_chi(c2_SWAP,2)
c3_CXI = right_ext_chi(c2_CX,2)
c3_ICX = left_ext_chi(c2_CX,2)



S2_HI = get_supersuper_from_chi(c2_HI, P2,2)
S2_IH = get_supersuper_from_chi(c2_IH, P2,2)
S2_CX = get_supersuper_from_chi(c2_CX, P2,2)

S3_SWAPI = get_supersuper_from_chi(c3_SWAPI, P3, 3)
S3_HII = get_supersuper_from_chi(c3_HII, P3, 3)
S3_IHI = get_supersuper_from_chi(c3_IHI, P3, 3)
S3_IIH = get_supersuper_from_chi(c3_IIH, P3, 3)

S3_CXI = get_supersuper_from_chi(c3_CXI, P3, 3)
S3_ICX = get_supersuper_from_chi(c3_ICX, P3, 3)
S3_CIX = S3_SWAPI @ S3_ICX @ S3_SWAPI
