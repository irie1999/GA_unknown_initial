#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=regular-o
#PJM -L node=1
#PJM --mpi proc=48
#PJM -L elapse=20:00:00
#PJM -g gu17
#PJM -j
#PJM -m e
#------- Program execution -------#
mpiexec ./main