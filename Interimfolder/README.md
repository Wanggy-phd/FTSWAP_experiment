# FTSWAP_experiment
python qiskit files to run process tomography experiments on IBM Q experience quantum chips.

Run only files in the root folder for intended use.

This code requires a QConfig file in /IBM_Q_Experience to function properly. (You have to provide your own API!)
For more info see: https://qiskit.org/documentation/install.html
For immediate Qconfig file https://github.com/QISKit/qiskit-terra/blob/stable/Qconfig.py.default

############################################################################################################################
The 'circuitname'_processtomo.py files are to:
- create tomography circuits for the circuit in /circuits/'circuitname'
- send those to the server,
- save the jobids of all batches to /Experiment_data/'run_type'/circuitname/filenamewithtimestamp
- save (overwriting) the jobids of all batches to /Experiment_data/last--jobids

Parameters: 
run_type: 'r' for real experiment or 's' for simulation.

reg: set to rue to register at IBMQExperience.

nr_batches: The number of batches over which all tomography circuits are divided while sending to the IBM server.
      Different batches are different jobs for the server.
      
If you want to run a different circuit the file from which the circuit must be imported should be specified at line 50.

############################################################################################################################
The file 'analyse' (Dutch for analysis) reads the jobids (&other data) from the /Experiment_data/ directory,
and retrieves the jobs from the IBM server.
If all jobs are completed, simple tomography analysis is run on the data:
- fitting the choi matrix
- comparing to the perfect operation as defined in the circuit file (to do!)
- plotting of the choi matrix

Parameters:
direct: if True, the jobids are obtained from /Experiment_data/last--jobids.

        if False, the jobids are obtained from /Experiment_data/run_type/circuit_name/filenamewithtimestamp,
        where timestamp is searched for in the local variables, being either a string or a datetime.datetime object
        If no such variable exist, the date and time of the experiment is prompted in the console.
        
run_type: 'r' for real experiment or 's' for simulation.

fit_method: fitting method for tomo data fitting. Check qiskit.tools.qcvv.tomography.fit_tomography_data for more info.

circuit_name: the name of the circuit to load the jobids of the experiment from. Not needed when direct == True.
