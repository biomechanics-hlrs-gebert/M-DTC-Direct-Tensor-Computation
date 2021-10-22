#!/bin/bash
#--------------------------------------------------------------------------
# Set the environment for the Hazel Hen system.
# (A Cray XC40 HPC-System running the Cray Linux Environment CLE 6.0)
#--------------------------------------------------------------------------
#
module unload PrgEnv-cray
module load PrgEnv-gnu
module load cray-libsci
module load cray-petsc
#
