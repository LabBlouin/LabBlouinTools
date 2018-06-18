#---------------------------Start program.job------------------------
#!/bin/bash

# The name of the job, can be whatever makes sense to you
#$ -N #NAME#

# The job should be placed into the queue 'all.q'.
#$ -q all.q

# Output, error stream files.
#$ -o stdout_#NAME#.dat
#$ -e stderr_#NAME#.dat

# The batchsystem should use the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batchsystem assumes to find the executable in this directory.
#$ -cwd

# Use the parallel environment "lam", which assigns two processes
# to one host. In this example, if there are not enough machines to run the
# mpi job on 120 processors the batchsystem can also use fewer than 120 but
# the job should not run on fewer than 30 processors.
#$ -pe mpich #NPMIN#-#NP#

# This is my email address for notifications. I want to have all notifications
# at the master node of this cluster.
#$ -M #EMAIL#

# Send me an email when the job is finished.
#$ -m e

# Set runtime length.
#$ -now y

# Pass external environment.
#$ -V

# This is the file to be executed.
#mpirun -np $NSLOTS 
#CMD#

#---------------------------End program.sge------------------------
