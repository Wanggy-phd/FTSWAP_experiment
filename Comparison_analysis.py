# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 10:59:52 2018

@author: Jarnd
"""

import Functions.Data_storage as store
import Analysis.Analysefunc as an
import Analysis.tomography_functions as tomoself
import Functions.Plotting as pt
import numpy as np

FT_chi_dict = store.load_chi_last('FTSWAP')
NFT_chi_dict = store.load_chi_last('NFTSWAP')

FT_chi_perror = FT_chi_dict['chi_perror']
NFT_chi_perror = NFT_chi_dict['chi_perror']

FT_chi_perror_stddv = FT_chi_dict['chistddv_error']
NFT_chi_perror_stddv = NFT_chi_dict['chistddv_error']

#
#FT_single_weight, FT_single_weight_stddv = an.get_single_weight(FT_chi_perror,FT_chi_perror_stddv)
#FT_multi_weight, FT_multi_weight_stddv = an.get_single_weight(FT_chi_perror,FT_chi_perror_stddv)
#NFT_single_weight, NFT_single_weight_stddv = an.get_multi_weight(NFT_chi_perror,NFT_chi_perror_stddv)
#NFT_multi_weight, NFT_multi_weight_stddv = an.get_multi_weight(NFT_chi_perror,NFT_chi_perror_stddv)

FT_ratio, FT_ratio_stddv = an.get_multi_single_ratio(FT_chi_perror,FT_chi_perror_stddv)
NFT_ratio, NFT_ratio_stddv = an.get_multi_single_ratio(NFT_chi_perror,NFT_chi_perror_stddv)

print('FT ratio: %f±z%f' %(FT_ratio, FT_ratio_stddv))

print('NFT ratio: %f±z%f' %(NFT_ratio, NFT_ratio_stddv))

z = 1.96
FT_ratio_conf_upper_bound = FT_ratio + z*FT_ratio_stddv
NFT_ratio_conf_lower_bound = NFT_ratio - z*NFT_ratio_stddv
print('Upper bound on 95%% confidence interval FT ratio: ',FT_ratio_conf_upper_bound)
print('Lower bound on 95%% confidence interval NFT ratio: ',NFT_ratio_conf_lower_bound)
if FT_ratio_conf_upper_bound <= NFT_ratio_conf_lower_bound:
    print('Succes!')