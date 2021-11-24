#!/bin/bash
###############################################################################
# Copyright 2021 Ralf Schneider
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################
#
# Unload MPT ------------------------------
module unload mpt/2.23
#
# Unload MPT ------------------------------
module load gcc/9.2.0
#
# MPI environment ------------------------
module load openmpi/4.0.4
#
mpi_prefix=/opt/hlrs/non-spack/mpi/openmpi/4.0.4-gcc-9.2.0/
export PATH=${mpi_prefix}/bin:$PATH
export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
#
# ----------------------------------------
# BLAS/LAPACK installation
# module load scalapack
export LAPACK_LIBPATH=$PWD/lib/lapack
# 
# export LAPACK_LIBPATH=/opt/hlrs/spack/rev-004_2020-06-17/scalapack/2.1.0-gcc-9.2.0-amna4d3j
#
# ----------------------------------------
# METIS installation
module load metis/5.1.0-int64
#
metis_prefix=/opt/hlrs/spack/rev-004_2020-06-17/metis/5.1.0-gcc-9.2.0-rdkkxlua
export METIS_INCPATH=${metis_prefix}/include
export METIS_LIBPATH=${metis_prefix}/lib
#
# ----------------------------------------
# PETSc installation
# module load petsc/3.12.2-int64
#
petsc_prefix=$PWD/petsc
export PETSC_INCPATH=${petsc_prefix}/include
export PETSC_LIBPATH=${petsc_prefix}/lib
export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
#
# ----------------------------------------
# Define std_out
export USE_STD_OUT=NO
#
# ----------------------------------------
# Root is a git repo?
export PROVIDES_GIT="NO"