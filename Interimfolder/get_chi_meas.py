# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:55:31 2018

@author: Jarnd

This file gathers results of an Indentity channel tomography experiment,
and computes from that the process matrix chi of the identity channel,
as a representation of the SPAM errors of the chip.
The process matrix is then saved using the save_chi_meas and save_chi_meas_last method from Functions.Data_storage
Optionally, the filter matrix B_filter for the calculated process matrix can be computed.
    For this to happen, set calcBfilter to True
    
To be done: calculate the standard deviation on chi_meas, and then on B_filter.
"""


import Functions.Data_storage as store
import Analysis.Analysis as an
import Analysis.tomography_functions as tomoself



fit_method = 'own' # Specify the fithmethod: 'own', 'leastsq' and 'wizard'
n = 2

calcBfilter = True # True to calculate the B_filter

#%% Gather the results from file

run_type = 'r'
circuit_name = 'Id'
timestamp = None
results_loaded = store.load_results(circuit_name, run_type, timestamp)          # Load results from Id experiment after providing date and time
timestamp = results_loaded['Experiment time']                                   # Reload the timestamp from the results for saving the chi matrix dicitonary


#%% Gather the tomo set and its outcomes from the results
tomo_set = results_loaded['Tomoset']
results = results_loaded['Results']
tomo_data = an.tomo.tomography_data(results, circuit_name, tomo_set)

#%% Tomography;
B_chi = tomoself.get_pauli_basis(n)
B_choi = tomoself.get_choi_basis(n, B_chi)

if fit_method == 'wizard':
    # Fitting choi with qiskit functions 'wizard' method and mapping choi to chi
    choi = an.fit_tomodata(tomo_data, method='wizard')                          # The wizard method makes chi CP itself
    chi = tomoself.choi_to_chi(choi, B_choi, n)                                 # The chi matrix is needed, not the Choi
elif fit_method == 'leastsq':
    # Fitting choi with qiskit functions 'leastsq' method and mapping choi to chi
    choi = an.fit_tomodata(tomo_data, method='leastsq')
    choi = an.make_CP(choi, n)                                                  # The leastsq method does not provide CP, so it still has to be done
    chi = tomoself.choi_to_chi(choi, B_choi, n)                                 # The chi matrix is needed, not the Choi
elif fit_method == 'own':
    # Fitting choi with own functions and mapping chi
    chi = an.fit_chi_own(tomo_data, tomo_set, n)
    chi = an.make_CP(chi, n)                                                    # My own method does not provide CP, so it still has to be done
#%% Calculate B_filter
Bfilter = None
if calcBfilter:
    Bfilter = an.get_measfiltermat(chi, B_chi, n)                               # This will take quite some time, can be done seperately after reloading the chi_meas
#%% Save to file and save to last used chi file
store.save_chi_meas(chi, timestamp, Bfilter)                                    # Save to file
store.save_chi_meas_last(chi, Bfilter)                                          # Save to 'last' file
