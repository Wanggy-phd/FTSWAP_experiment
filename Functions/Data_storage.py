# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 15:13:38 2018

@author: Jarnd
"""
import pickle as pickle
import os as os
from datetime import datetime


###############################################################################
def save_tomo_data(circuit_name, run_type, backend, tomo_data, nr_shots, notes=None):
    '''
    Save tomography data in a directory specified by run_type and circuit_name.
    Example: '/Experiment_data/Real_data/Hadamard' for a real experiment with name "Hadamard".
    The filename is then 
    '''
    timenow = datetime.now()
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'
    filepath = directory+circuit_name+'--' + \
        timenow.strftime("%m_%d-%H_%M_%S")+"--tomodata.txt"
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    data = {'Experiment time': timenow, 'Circuit name': circuit_name,
            'Type': run_type, 'Backend': backend, 'Data': tomo_data,
            'Shot number': nr_shots, 'Notes': notes}
    pickle.dump(data, open(filepath, "wb"))
    return data


def load_tomo_data(circuit_name, run_type, timestamp=None):
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'

    if timestamp == None:
        date = input('Date of experiment (mm_dd):')
        time = input('Time of experiment (hh_mm_ss):')
        filepath = directory+circuit_name+'--'+date+'-'+time+'--tomodata.txt'
        data_load = pickle.load(open(filepath, 'rb'))
    elif type(timestamp) == datetime:
        filepath = directory+circuit_name+'--' + \
            timestamp.strftime("%m_%d-%H_%M_%S")+"--tomodata.txt"
        data_load = pickle.load(open(filepath, "rb"))
    elif type(timestamp) == str:
        filepath = directory+circuit_name+'--'+timestamp+"--tomodata.txt"
        data_load = pickle.load(open(filepath, "rb"))
    return data_load


###############################################################################
def save_jobdata(circuit_name, jobs_data, tomo_set, backend, nr_shots, nr_batches, run_type, unitary, notes=None):
    '''
    Save data of the jobs created in a tomography experiment to a text file and a pickle format file.
    The directory and filenemes are specified by run_type and circuit_name and the time of saving.
    Example: '/Experiment_data/Real_data/Hadamard/Hadamard--MM_DD-HH_MM_SS-jobdata.pickle' (and .txt),
    for a real experiment with name "Hadamard".
    In the .pickle file is a dictionary of:
        -The circuit name
        -The experiement time (the time of saving)
        -The tomography set (as specified by the create_tomography functions)
        -The runtype ('r' or 's')
        -The used backend
        -The number of shots
        -The unitary operation of the circuit as a numpy matrix
        -The number of batches over which all tomography circuits are divided
        -A list of the actual jobs
        -Optional notes
    The function also saved a .txt file with the same directory and name that stores the jobids for quick reference
    '''
    timenow = datetime.now()
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'
    filepathtxt = directory+circuit_name+'/'+circuit_name + \
        '--'+timenow.strftime("%m_%d-%H_%M_%S")+"--jobdata.txt"
    filepathpickle = directory+circuit_name+'/'+circuit_name + \
        '--'+timenow.strftime("%m_%d-%H_%M_%S")+"--jobdata.pickle"
    os.makedirs(os.path.dirname(filepathtxt), exist_ok=True)
    fo = open(filepathtxt, 'w')
    fo.write('Job id\'s for experiment ran on '+jobs_data[0]['Date']+'\n')
    fo.write('Data saved on '+timenow.strftime("%m_%d-%H_%M_%S")+'\n')
    if run_type == 's':
        exptype = 'simulation'
    elif run_type == 'r':
        exptype = 'real'
    fo.write('Experiment type is '+exptype+'\n')
    fo.write('Number of batches is '+str(nr_batches)+'\n\n')
    for job in jobs_data:
        fo.write('Job nr %d/%d: \n' % (job['batchno']+1, nr_batches))
        fo.write('Date: '+job['Date']+'\n')
        fo.write('Job id:\n'+job['Jobid']+'\n\n')
    fo.close()
    datadict = {'Experiment time': timenow, 'Circuit name': circuit_name,
                'Type': run_type, 'Backend': backend, 'Data': jobs_data,
                'Tomoset': tomo_set, 'Shot number': nr_shots,
                'Unitary': unitary, 'Batchnumber': nr_batches, 'Notes': notes}
    fo = open(filepathpickle, 'wb')
    pickle.dump(datadict, fo)
    fo.close


def load_jobdata(circuit_name, run_type, timestamp=None):
    '''
    Loads the dictionary of jobdata stored in the .pickle file at location:
    'Experiment_data/Real_data/Circuitname/Circuitname--MM_DD-HH_MM_SS--jobdata.pickle'
    for a real experiment and conversely for a simulation experiment.
    The data should be stored using the save_jobdata() function
    '''
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'

    if timestamp == None:
        date = input('Date of experiment (mm_dd):')
        time = input('Time of experiment (hh_mm_ss):')
        filepath = directory+circuit_name+'/'+circuit_name + \
            '--'+date+'-'+time+'--jobdata.pickle'
        data_load = pickle.load(open(filepath, 'rb'))
    elif type(timestamp) == datetime:
        filepath = directory+circuit_name+'/'+circuit_name+'--' + \
            timestamp.strftime("%m_%d-%H_%M_%S")+"--jobdata.pickle"
        data_load = pickle.load(open(filepath, "rb"))
    elif type(timestamp) == str:
        filepath = directory+circuit_name+'/' + \
            circuit_name+'--'+timestamp+"--jobdata.pickle"
        data_load = pickle.load(open(filepath, "rb"))
    return data_load


###############################################################################
def save_last(circuit_name, jobdata, tomo_set, backend, nr_shots, nr_batches, run_type, unitary, notes=None):
    '''
    Acts the same as save_jobdata but stores the pickle file as 'Experiment_data/last--jobdata.pickle' for quick loading.
    Warning: this function ALWAYS overwrites a potential existing file.
    '''
    timenow = datetime.now()
    filepathpickle = "Experiment_data/last--jobdata.pickle"
    os.makedirs(os.path.dirname(filepathpickle), exist_ok=True)
    datadict = {'Experiment time': timenow, 'Circuit name': circuit_name,
                'Type': run_type, 'Backend': backend, 'Data': jobdata,
                'Tomoset': tomo_set, 'Shot number': nr_shots,
                'Unitary': unitary, 'Batchnumber': nr_batches, 'Notes': notes}
    fo = open(filepathpickle, 'wb')
    pickle.dump(datadict, fo)
    fo.close


def load_last():
    '''
    Load the jobdata dictionary stored in 'Experiment_data/last--jobdata.pickle'
    '''
    filepathpickle = "Experiment_data/last--jobdata.pickle"
    data_load = pickle.load(open(filepathpickle, 'rb'))
    return data_load

###############################################################################


def save_results(circuit_name, timestamp, run_type, backend, jobids, tomo_set, nr_batches, nr_shots, results, unitary, calibrationdata, notes=None):
    '''
    Save data of the results from a tomography experiment to a text file and a pickle format file.
    The directory and filenemes are specified by run_type and circuit_name and the time of saving.
    Example: '/Experiment_data/Real_data/Hadamard/Hadamard--MM_DD-HH_MM_SS-results.pickle' (and .txt),
    for a real experiment with name "Hadamard".
    In the .pickle file is a dictionary of:
        'Circuit name'-The circuit name
        'Experiment time'-The experiement time
        'Type'-The runtype ('r' or 's')
        'Backend'-The used backend
        'Jobids'-The original job ids
        'Tomoset'-The tomography set (as specified by the create_tomography functions)
        'Batchnumber'-The number of batches over which all tomography circuits are divided
        'Shot number'-The number of shots per circuit
        'Results'-The results of all jobs together as gotten from the backend
        'Unitary'-The unitary operation of the circuit as a numpy matrix
        'Calibration'-Data of the calibration of the quantum chip
        'Notes'-Optional notes
    '''
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'

    filepathpickle = directory+circuit_name+'/'+circuit_name + \
        '--'+timestamp.strftime("%m_%d-%H_%M_%S")+"--results.pickle"
    os.makedirs(os.path.dirname(filepathpickle), exist_ok=True)
    datadict = {'Experiment time': timestamp, 'Circuit name': circuit_name,
                'Type': run_type, 'Backend': backend, 'Results': results,
                'Tomoset': tomo_set, 'Shot number': nr_shots, 'Jobids': jobids,
                'Batchnumber': nr_batches, 'Calibration': calibrationdata,
                'Unitary': unitary, 'Notes': notes}
    fo = open(filepathpickle, 'wb')
    pickle.dump(datadict, fo)
    fo.close


def load_results(circuit_name, run_type, timestamp=None):
    '''
    Loads the dictionary of results data stored in the .pickle file at location:
    'Experiment_data/Real_data/Circuitname/Circuitname--MM_DD-HH_MM_SS--results.pickle'
    for a real experiment and conversely for a simulation experiment.
    The data should be stored using the save_results() function
    '''
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'
    if timestamp == None:
        date = input('Date of experiment (mm_dd):')
        time = input('Time of experiment (hh_mm_ss):')
        filepath = directory+circuit_name+'/'+circuit_name + \
            '--'+date+'-'+time+'--results.pickle'
        data_load = pickle.load(open(filepath, 'rb'))
    elif type(timestamp) == datetime:
        filepath = directory+circuit_name+'/'+circuit_name+'--' + \
            timestamp.strftime("%m_%d-%H_%M_%S")+'--results.pickle'
        data_load = pickle.load(open(filepath, "rb"))
    elif type(timestamp) == str:
        filepath = directory+circuit_name+'/' + \
            circuit_name+'--'+timestamp+'--results.pickle'
        data_load = pickle.load(open(filepath, "rb"))
    return data_load

###############################################################################


def save_chi_meas(chi_meas, timestamp=None, Bfilter=None):
    '''
    Save the process matrix from a tomography experiment on the identity channel.
    This process matrix can be used to calculate the filter matrix Bfilter (optionally included already)
    Using this filter matrix the SPAM errors can be filtered out.
    The process matrix and optional Bfilter matrix are stored at:
        'Experiment_data/Real_data/Id/Id--MM_DD-HH_MM_SS--chimeas.pickle'
    The dictionary contains:
        'chi_meas' - The process matrix chi
        'B_filter' - The optional Bfilter matrix
    '''
    directory = 'Experiment_data/Real_data/Id/Id'
    if timestamp == None:
        date = input('Date of experiment (mm_dd):')
        time = input('Time of experiment (hh_mm_ss):')
        filepath = directory+'--'+date+'-'+time+'--chimeas.pickle'
    elif type(timestamp) == datetime:
        filepath = directory+'--' + \
            timestamp.strftime("%m_%d-%H_%M_%S")+'--chimeas.pickle'
    elif type(timestamp) == str:
        filepath = directory+'--'+timestamp+'--chimeas.pickle'
    datadict = {'chi_meas': chi_meas, 'B_filter': Bfilter}
    fo = open(filepath, 'wb')
    pickle.dump(datadict, fo)
    fo.close


def load_chi_meas(timestamp=None):
    '''
    Loads the process matrix dictionary for the Identity channel stored in:
    'Experiment_data/Real_data/Id/Id--MM_DD-HH_MM_SS--chimeas.pickle'
    If a timestamp is provided the corresponding file is loaded, else the date and time is prompted
    '''
    directory = 'Experiment_data/Real_data/Id/Id'
    if timestamp == None:
        date = input('Date of experiment (mm_dd):')
        time = input('Time of experiment (hh_mm_ss):')
        filepath = directory+'--'+date+'-'+time+'--chimeas.pickle'
    elif type(timestamp) == datetime:
        filepath = directory+'--' + \
            timestamp.strftime("%m_%d-%H_%M_%S")+'--chimeas.pickle'
    elif type(timestamp) == str:
        filepath = directory+'--'+timestamp+'--results.pickle'
    return pickle.load(open(filepath, 'rb'))


def save_chi_meas_last(chi_meas, Bfilter=None):
    '''
    Save the process matrix from a tomography experiment on the identity channel.
    This process matrix can be used to calculate the filter matrix Bfilter (optionally included already)
    Using this filter matrix the SPAM errors can be filtered out.
    The process matrix and optional Bfilter matrix are stored at:
        'Experiment_data/Real_data/Id/chi_last.pickle'
        Warning: this function will ALWAYS overwrite a file stored there!
    The dictionary contains:
        'chi_meas' - The process matrix chi
        'B_filter' - The optional Bfilter matrix
    '''
    filepath = 'Experiment_data/Real_data/Id/chi_last.pickle'
    datadict = {'chi_meas': chi_meas, 'B_filter': Bfilter}
    fo = open(filepath, 'wb')
    pickle.dump(datadict, fo)
    fo.close


def load_chi_meas_last():
    '''
    Loads the process matrix as in save_chi_meas but from file 'Experiment_data/Real_data/Id/chi_last.pickle'
    '''
    return pickle.load(open('Experiment_data/Real_data/Id/chi_last.pickle', 'rb'))

###############################################################################

def save_chi(chi_dict, timestamp, run_type, circuit_name):
    '''
    Saves a dictionary of process matrices of any tomography experiment under 
    'Experiment_data/Real_data/circuit_name/circuit_name--MM_DD--HH_MM_SS-chimat.pickle'
    And /Simulation_data/ if the specified runtype is 's'
    '''
    if run_type == 's':
        directory = 'Experiment_data/Simulation_data/'
    elif run_type == 'r':
        directory = 'Experiment_data/Real_data/'
    filepath = directory+circuit_name+'/'+circuit_name + \
        '--'+timestamp.strftime("%m_%d-%H_%M_%S")+"--chimat.pickle"
    fo = open(filepath, 'wb')
    pickle.dump(chi_dict, fo)
    fo.close

def save_chi_last(chi_dict, circuit_name):
    '''
    Saves a dictionary of process matrices of any tomography experiment under 
    'Experiment_data/chi--circuit_name--.pickle'
    For quick reference use
    '''
    filepath = 'Experiment_data/'+'chi--'+circuit_name+'--.pickle'
    fo = open(filepath, 'wb')
    pickle.dump(chi_dict, fo)
    fo.close

def load_chi_last(circuit_name):
    '''
    Load the dictionary of process matrices from 'Experiment_data/chi--circuit_name--.pickle'
    as saved by save_chi_last()
    '''
    filepath = 'Experiment_data/'+'chi--'+circuit_name+'--.pickle'
    fo = open(filepath, 'rb')
    chi_dict = pickle.load(fo)
    fo.close
    return chi_dict