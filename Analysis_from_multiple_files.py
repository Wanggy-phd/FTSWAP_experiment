# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 15:36:13 2018

@author: Jarnd
"""

"""
Created on Thu Aug 30 14:55:31 2018

@author: Jarnd
This file is used to analyse tomography data as gotten from the file Get_results after running a processtomo experiment and sending it to the IBM Q experience.
If direct = True, the last saved jobdata is used to gather the results. This only works if Get_results is ran on this specific job as the last too.
If direct = False, the specific experiment must be specified; the run type, the circuit name and the date and time is prompted from the user.

This file:
    -Calculates the chi matrix and choi matrix of the tomography experiment using the method specified by fit_method
    -Checks if chi and choi are CP and TP
    -Filters out the measurement with the use of the Bfilter matrix as calculated using an.get_measfiltermat() function.
        chi_meas is approximated by the chi matrix of an tomography experiment on the identity channel.
    -Calculates the error matrices using the an.get_chi_error() function
    -Calculates the process and channel fidelity and prints them
    -Saves the various process matrices to file
    -Plots city plots of the chi matrix and the error matrix
    -Plots the pauli twirl weights and calculates the multi- to single qubit ratio
"""

#%% Importing the needed functions and modules.
import Functions.Data_storage as store
import Analysis.Analysefunc as an
import Analysis.tomography_functions as tomoself
import Functions.Plotting as pt
import numpy as np

#%% Specify the fitting method and loading method
fit_method = 'own'


jobs = ['11_17-17_30_34']#,'11_17-17_31_06','11_17-17_31_30','11_17-17_32_20']#,'11_17-17_32_53']
circuit_name = 'FTSWAP'

#%% Load the results from file
meas_data_all = [];
for jobc,job in enumerate(jobs):
    run_type = 's'
    timestamp = job
    results_loaded = store.load_results(circuit_name, run_type, timestamp)


    #%% Gather the tomo set from the loaded results and its outcomes from the results
    tomo_set = results_loaded['Tomoset']
    results = results_loaded['Results']
    timestamp = results_loaded['Experiment time']
    Unitary = results_loaded['Unitary']
    tomo_data = an.tomo.tomography_data(results, circuit_name, tomo_set)
    meas_data_all.append(tomo_data['data'])
    n = int(np.log2(np.shape(Unitary)[0]))

#%% Tomography; obtaining chi and choi matrices, check CP and TP
B_chi = tomoself.get_pauli_basis(n, normalise = False)
B_choi = tomoself.get_choi_basis(n, B_chi)

((chi,chistddv),(choi,choistddv)) = an.fit_tomodata_multiple(meas_data_all, tomo_set, B_chi, B_choi, n, stddv=True)
print('CP:', an.check_CP(chi))
print('TP:', an.check_TP(chi, B_chi, n))




#%% Filter meas errors out
chi_filtered = chi
choi_filtered = choi
if run_type == 'r' and n == 2:
    Bfilter = store.load_chi_meas_last()['B_filter']
    chi_filtered = an.filter_chi_meas(chi, Bfilter, n)
    choi_filtered = tomoself.chi_to_choi(chi_filtered, B_choi, n)
    del Bfilter

#%% Calculate error matrices
chi_perror = an.get_chi_error(chi_filtered, B_chi, Unitary)
chi_unfiltered_perror = an.get_chi_error(chi, B_chi, Unitary)
chi_perror_stddv = an.get_chi_error(chistddv, B_chi, Unitary)

#%% Calculate traces and fidelities
process_fidelity = an.process_fidelity(chi_perror, n)
process_fidelity_unfiltered = an.process_fidelity(chi_unfiltered_perror, n)
channel_fidelity = an.channel_fidelity(chi_perror, B_choi, n)


print('Process fidelity from error matrix:', np.around(np.abs(process_fidelity),3))
print('Process fidelity from error matrix unr:', np.abs(process_fidelity))
print('Channel fidelity from Choi matrix:', np.around(np.abs(channel_fidelity),3))
print('Process fidelity from unfiltered error matrix:',
      np.around(np.abs(process_fidelity_unfiltered),3))

#%% Store the results
chi_dict = {'chi': chi , 'chi_filtered' : chi_filtered, 'chi_perror' : chi_perror
            , 'chistddv' : chistddv , 'chistddv_error' : chi_perror_stddv}

store.save_chi_last_mult(chi_dict, circuit_name)
store.save_chi_mult(chi_dict, timestamp, run_type, circuit_name)
#%% delete unwanted variables
del fit_method, results_loaded, run_type, timestamp,tomo_data, tomo_set
#%% Plotting
pt.plot_city(chi, tomoself.get_pauli_names(n), r'$\chi$', circuit_name+'chi', err = False)
pt.plot_city(chi_perror, tomoself.get_pauli_names(n), r'$\chi_{error}$', circuit_name+'chi_perror', err = True)
#%% Plotting with standard deviation
#pt.plot_city_with_var(chi,chistddv, tomoself.get_pauli_names(n),r'$\chi_{filtered}$')
#pt.plot_city_with_var(chi_perror,chi_perror_stddv, tomoself.get_pauli_names(n),r'$\chi_{error}$')
#%% Plotting of Pauli twirl and printing ratio
#pt.plot_pauli_weight(np.diag(chi_perror),2*np.diag(chi_perror_stddv),tomoself.get_pauli_names(n),n)
ratio = an.get_multi_single_ratio(chi_perror, chi_perror_stddv)
print('Ratio multi to single: %f +- %f' %(ratio[0],ratio[1]))