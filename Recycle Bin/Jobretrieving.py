# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 14:15:48 2018

@author: Jarnd
"""

from IBMQuantumExperience import IBMQuantumExperience as IBMQ
import Qconfig
from qiskit.backends.ibmq.ibmqjob import IBMQJob
api = IBMQ(Qconfig.APItoken, Qconfig.config)


batchinfo = False
jobs_last = api.get_jobs(limit=9)
result_list = None

for job in jobs_last:
    print('Time is', job['creationDate'])
    print('Backend is ', job['backend'])
    print('Status is :', job['status'], '\n')
if batchinfo:
    batchnr = 1
    for batch in jobs_last:
        print('\nInfo of batch nr. %d \n' % batchnr)
        for qasm in batch['qasms']:
            print('Status:', qasm['status'])
            if qasm['status'] == 'DONE':
                print('Date:', qasm['result']['date'])
    batchnr += 1


result = IBMQJob.from_api(jobs_last[-1], api, False).result()
for i in range(7, -1, -1):
    result += IBMQJob.from_api(jobs_last[i], api, False).result()
