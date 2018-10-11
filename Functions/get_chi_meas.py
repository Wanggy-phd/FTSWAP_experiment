# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:55:31 2018

@author: Jarnd
"""

#import sys
# sys.path.append('../')

import Functions.Data_storage as store
import Analysis.Analysis as an
import Analysis.tomography_functions as tomoself


fit_method = 'own'
n = 2

calcBfilter = True

#%% Gather the results from file

run_type = 'r'
circuit_name = 'Id'
timestamp = None
results_loaded = store.load_results(circuit_name, run_type, timestamp)
timestamp = results_loaded['Experiment time']


#%% Gather the tomo set and its outcomes from the results
tomo_set = results_loaded['Tomoset']
results = results_loaded['Results']
tomo_data = an.tomo.tomography_data(results, circuit_name, tomo_set)

#%% Tomography;
B_chi = tomoself.get_pauli_basis(n)
B_choi = tomoself.get_choi_basis(n, B_chi)

if fit_method == 'wizard':
    # Fitting choi with qiskit functions 'wizard' method and mapping choi to chi
    choi = an.fit_tomodata(tomo_data, method='wizard')
    chi = tomoself.choi_to_chi(choi, B_choi, n)
elif fit_method == 'leastsq':
    # Fitting choi with qiskit functions 'leastsq' method and mapping choi to chi
    choi = an.fit_tomodata(tomo_data, method='leastsq')
    choi = an.make_CP(choi, n)
    chi = tomoself.choi_to_chi(choi, B_choi, n)
elif fit_method == 'own':
    # Fitting choi with own functions and mapping chi
    chi = an.fit_chi_own(tomo_data, tomo_set, n)
    chi = an.make_CP(chi, n)
#%% Calculate B_filter
Bfilter = None
if calcBfilter:
    Bfilter = an.get_measfiltermat(chi, B_chi, n)
#%% Save to file
store.save_chi_meas(chi, timestamp, Bfilter)
store.save_chi_meas_last(chi, Bfilter)
