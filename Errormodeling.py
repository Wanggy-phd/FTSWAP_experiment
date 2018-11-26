# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 15:01:15 2018

@author: Jarnd
"""

#import Functions.Plotting as pt
#from Analysis.tomography_functions import get_pauli_names as paulilist
import numpy as np
import Functions.Simpleerrorfunc as sef
import Functions.Data_storage as store
import Analysis.tomography_functions as tomoself
import Analysis.Analysefunc as an
import Functions.Plotting as pt
import time
import datetime
#%%
chidict_FTSWAP = store.load_chi_last('FTSWAP')
chidict_NFTSWAP = store.load_chi_last('NFTSWAP')

chiFT = chidict_FTSWAP['chi_filtered']
chiNFT = chidict_NFTSWAP['chi_filtered']

B_chi = tomoself.get_pauli_basis(2)
B_choi = tomoself.get_choi_basis(2, B_chi)
choi_FT = tomoself.chi_to_choi(chiFT, B_choi, 2)
choi_NFT = tomoself.chi_to_choi(chiNFT, B_choi, 2)

chi_tot_meas = np.kron(chiNFT,chiFT)


Fidft_meas = an.get_chi_error(chiFT, B_chi, sef.SWAP)[0,0]/4
Fidnft_meas = an.get_chi_error(chiNFT, B_chi, sef.SWAP)[0,0]/4
Fidtot_meas = an.get_chi_error(chi_tot_meas, tomoself.get_pauli_basis(4,normalise = False), np.kron(sef.SWAP,sef.SWAP))[0,0]/16

gridnr = 20;

pslist = np.linspace(0, 0.05, num = gridnr)
pmlist = np.linspace(0, 0.10, num = gridnr)

fidnft = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')
fidft = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')

tdnft = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')
tdft = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')
tdtot = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')

fiddiff_tot = np.zeros((len(pslist),len(pmlist)),dtype = 'complex')

P2 = sef.P2
P3 = sef.P3
#%%
nr_loops = 0
average = 0
for sc, ps in enumerate(pslist):
    
    c1_dep = sef.get_chi_dep(ps,1, sef.normalisation)
    
    c2_depI = sef.right_ext_chi(c1_dep,1, normalise=sef.normalisation)
    c2_Idep = sef.left_ext_chi(c1_dep,1, normalise=sef.normalisation)
    
    c3_depII = sef.right_ext_chi(c2_depI,2,normalise = sef.normalisation)
    c3_IdepI = sef.left_ext_chi(c2_depI,2,normalise = sef.normalisation)
    c3_IIdep = sef.left_ext_chi(c2_Idep,2,normalise = sef.normalisation)
    
    S2_depI = sef.get_supersuper_from_chi(c2_depI, sef.P2, 2)
    S2_Idep = sef.get_supersuper_from_chi(c2_Idep, sef.P2, 2)
    
    S3_depII = sef.get_supersuper_from_chi(c3_depII, P3, 3)
    S3_IdepI = sef.get_supersuper_from_chi(c3_IdepI, P3, 3)
    S3_IIdep = sef.get_supersuper_from_chi(c3_IIdep, P3, 3)
    
    S2_fHH = S2_depI @ S2_Idep @ sef.S2_HI @ sef.S2_IH
    
    S3_fHHI = S3_depII @ S3_IdepI @ S3_IIdep @ sef.S3_HII @ sef.S3_IHI
    S3_fHIH = S3_depII @ S3_IdepI @ S3_IIdep @ sef.S3_HII @ sef.S3_IIH
    S3_fIHH = S3_depII @ S3_IdepI @ S3_IIdep @ sef.S3_IHI @ sef.S3_IIH
    
    for mc, pm in enumerate(pmlist):
        start = time.time()
        ######################################################################
        
        
        
        c2_dep2 = sef.get_chi_dep(pm,2,sef.normalisation)
        
        
        c3_dep2I = sef.right_ext_chi(c2_dep2,2,normalise = sef.normalisation)
        c3_Idep2 = sef.left_ext_chi(c2_dep2,2,normalise = sef.normalisation)


        #%%
        
        S2_dep2 = sef.get_supersuper_from_chi(c2_dep2, P2, 2)


        S3_dep2I = sef.get_supersuper_from_chi(c3_dep2I, P3, 3)
        S3_Idep2 = sef.get_supersuper_from_chi(c3_Idep2, P3, 3)
        S3_depI2 = sef.S3_SWAPI @ S3_Idep2 @ sef.S3_SWAPI

        
        #%%
        S2_fCX = S2_dep2 @ sef.S2_CX
        
        
        S3_fCXI = S3_IIdep @ S3_dep2I @ sef.S3_CXI
        S3_fCIX = S3_IdepI @ S3_depI2 @ sef.S3_CIX
        S3_fICX = S3_depII @ S3_Idep2 @ sef.S3_ICX
        
        
        S3_fSWAPI = S3_fCXI @ S3_fHHI @ S3_fCXI @ S3_fHHI @ S3_fCXI
        S3_fSWAIP = S3_fCIX @ S3_fHIH @ S3_fCIX @ S3_fHIH @ S3_fCIX
        S3_fISWAP = S3_fICX @ S3_fIHH @ S3_fICX @ S3_fIHH @ S3_fICX
        
        #%%
        S2nft_tot = S2_fCX @ S2_fHH @ S2_fCX @ S2_fHH @ S2_fCX
        S3ft_tot = S3_fISWAP @ S3_fSWAPI @ S3_fSWAIP
        
        
        c2nft_tot = sef.get_chi_from_supersuper(S2nft_tot, P2, 2)
        c3ft_tot = sef.get_chi_from_supersuper(S3ft_tot, P3, 3)
        c2ft_tot = sef.chi3_trace_3(c3ft_tot)
        
        choi_nft_comp = tomoself.chi_to_choi(c2nft_tot, B_choi, 2)
        choi_ft_comp = tomoself.chi_to_choi(c2ft_tot, B_choi, 2)
        
        c2nft_error = an.get_chi_error(c2nft_tot,sef.P2norm,sef.SWAP)
        c2ft_error = an.get_chi_error(c2ft_tot,sef.P2norm,sef.SWAP)
        
        c4tot = np.kron(c2nft_tot,c2ft_tot)
        c4tot_error = an.get_chi_error(c4tot,sef.get_pauli_basis(4,normalise=False),np.kron(sef.SWAP,sef.SWAP))
        
        
        fidnft[sc,mc] = c2nft_error[0,0]
        fidft[sc,mc] = c2ft_error[0,0]
        fiddiff_tot[sc,mc] = np.abs(c4tot_error[0,0]- Fidtot_meas)
        tdnft[sc,mc] = an.trace_dist(choi_nft_comp,choi_NFT)
        tdft[sc,mc] = an.trace_dist(choi_ft_comp,choi_FT)
        tdtot[sc,mc] = an.trace_dist(np.kron(choi_nft_comp,choi_ft_comp),np.kron(choi_NFT,choi_FT))
        

        #%%
        stop = time.time()
        diff = stop-start
        lapsesleft = gridnr**2 - (nr_loops+1)
        average = (nr_loops * average + diff)/(nr_loops+1)
        timeleft = average*lapsesleft
        nr_loops += 1
        percentage = nr_loops/(gridnr**2)

        print("Done: %d %%" %(percentage * 100))
        print('Time left:', time.strftime(format('H:%H M:%M S:%S'),time.gmtime(average*lapsesleft)))


#for name in dir():
#    if name.startswith('S') or name.startswith('c'):
#        del globals()[name]
#del name

#del S2_dep2, S2_depI, S2_fCX, S2_fHH, S2_Idep, S2nft_tot, S
#del ps, pm, mc, sc, P2, P3
#del average, nr_loops, timeleft, stop, diff, lapsesleft, percentage, start
#%%
data = {'NFTfid' : fidnft, 'FTfid' : fidft, 'totalfiddiff' : fiddiff_tot, 'NFTtd' : tdnft, 'FTtd' : tdft, 'Tottd' : tdtot, 'pslist' : pslist, 'pmlist' : pmlist}
store.save_chi(data, datetime.datetime.now(), 's','Simulations')   

pt.plot_contour_fid(pslist,pmlist,fidnft,'NFT Fidelities',levels=40)
pt.plot_contour_fid(pslist,pmlist,fidft,'FT Fidelities',levels=40)
pt.plot_contour_fid(pslist,pmlist,tdnft,'NFT Trace Distance',levels=40)
pt.plot_contour_fid(pslist,pmlist,tdft,'FT Trace Distance',levels=40)
#%%
pt.plot_contour_fid(pslist,pmlist,tdtot,'Distance total',levels=30)
pt.plot_contour_fid(pslist,pmlist,fiddiff_tot,'Summing fid distance total',levels=30)
#%%
pt.plot_surface(pslist,pmlist,tdtot)
pt.plot_surface(pslist,pmlist,fiddiff_tot)
#plt.contourf(pslist,pmlist, np.real(fid))
#%%
