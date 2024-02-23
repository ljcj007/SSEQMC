#$ -S /bin/bash
#$ -pe orte 1
#$ -cwd
#

/opt/mpich2-intel/bin/mpirun -np ${NSLOTS} ./a.out
