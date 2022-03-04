# g09log
Read the log file of g09, and gives some useful parameters.
Contact me: zhouych@foxmail.com


# compile it before use
gfortran g09log.f90 -o g09log

# How to use:
g09log logfile  jobs 

logfile: the name of gaussian output log filr 
jobs: 0-3
 0: chk converge energy
 1: obtain homo and lumo and its orbital number;
 2: write gjf from log, have to check the functional line
 3: get total energy
