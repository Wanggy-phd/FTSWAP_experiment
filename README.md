# FTSWAP_experiment
Python files to run tomography experiments on the IBM Q Experience using qiskit. Specifically used for a Fault Tolerant Swap circuit.

The overall goal of the experiment is to compare the standard implementation of a SWAP circuit with a Fault Tolerant (FT) gadget that does the same thing.

The intended workflow:
## Setting up the experiment
- First specify the circuit and the unitary in a file saved in circuits using the circuit_template. The filename should be circuit_<circuitname>.py.
- Create a new instance of an experiment using the template_proccesstomo.py and save it to the root folder.
- Specify the run type and other variables in that file.
- Create an account at the IBM Q experience to get a API token and store that token in the file /IBM_Q_Experience/Qconfig.py. Note: process tomography is costly so doing these will spend many credits.

## Running the experiment
- Run the file FROM THE ROOT DIRECTORY. The experiments are send to the Q experience, and data on the jobs is automatically saved to /Experiment_Data/...

## Gathering the data
- Run the file get_results.py from the root directory. 
	Set direct to True if the last ran experiment is to be gathered, if direct = False the experinent type, name and time will be prompted. The file will print the status of all jobs within the experiment, and if they are all done will gather the results and save them to file.
	
## Getting data on SPAM errors for SPAM filtering
- The Id_processtomography.py experiment should be ran to obtain data for the SPAM error filtering methods. Instead of analyzing that data (after the results have been gathered) run the file get_chi_meas from the root folder. There is already data included for the SPAM errors on the ibmqx5 device, so this step can be skipped if that device is used.

## Analysing data
- Run the code Analysis_from_file.py from the root directory. 
	Here also there is a parameter 'direct', with the same behaviour as in get_results.py.
	Furthermore, there are other parameters to be set, with the most important one the 'fit_method'.
	See the file for details.
	This file now calculates the process matrix, the Choi matrix, the error matrix and the fidelities. Furthermore, it plots various representations and calculates the ratio of multi to single qubit errors.
	See the file for more details.

	
###############################################################################

There are some other files that have a more or less standalone function:
- The files in the directory /Errormodeling/ can be ran to obtain sumlation data and plots on the circuits in FTSWAP and NFTSWAP with a simple dephasing error model.

- The comparison analysis can be used to compare the FTSWAP and NFTSWAP experiments outcomes and calculates the multi-tosingle qubit error rates for both implementations.
