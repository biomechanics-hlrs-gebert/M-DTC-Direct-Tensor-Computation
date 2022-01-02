#!/bin/bash
###############################################################################
# Copyright 2021 Johannes Gebert
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
#------------------------------------------------------------------------------
# Set the environment for the Julius I system.
#    Julius I is a Whiskey Lake Intel(R) Core(TM) i5-8365U CPU @ 1.60GHz
#    laptop running Arch Linux x86_64, Kernel > 5.14.12-arch1-1, zsh > 5.8
#------------------------------------------------------------------------------
#
# MPI environment ------------------------
mpi_prefix=/opt/mpi/openmpi-NO_F08-4.1.0
export PATH=${mpi_prefix}/bin:$PATH
export LD_LIBRARY_PATH=${mpi_prefix}/lib:$LD_LIBRARY_PATH
# ----------------------------------------
# BLAS/LAPACK installation
export LAPACK_LIBPATH=$PWD/lib/lapack
#
# ----------------------------------------
# METIS installation
metis_prefix=/opt/metis/metis-5.1.0
export METIS_INCPATH=${metis_prefix}/include
export METIS_LIBPATH=${metis_prefix}/lib
#
# ----------------------------------------
# PETSc installation
petsc_prefix=/opt/petsc/petsc-3.15
export PETSC_INCPATH=${petsc_prefix}/include
export PETSC_LIBPATH=${petsc_prefix}/lib
export LD_LIBRARY_PATH=${petsc_prefix}/lib:$LD_LIBRARY_PATH
#
# ----------------------------------------
# Define std_out
export USE_STD_OUT=YES
# ---------------------------------------------------------------------------
# Gnu Debugger - make tmpi available / check prerequisites
#
tmpi_prefix="/opt/tmpi"
#
# ----------------------------------------
# Root is a git repo?
export PROVIDES_GIT="YES"
#
export PATH=${tmpi_prefix}:$PATH
#
tools=( gdb tmpi tmux mpirun )
#
dbg_err=0
#
if [ "$NO_OUTPUT" != "YES" ]; then
    for tool in "${tools[@]}"; do
        echo -n "-- "
        if ! which ${tool} ; then # > /dev/null 2> /dev/null (to suppress cmd line output)
            echo "-- Please provide ${yellow}${tool}${nc} to use gdb with mpi."
            dbg_err=1
        fi
    done
    #
    if [[ $dbg_err == 0 ]]; then
        echo "--"
        echo "-- Usage of the GNU Debugger:"
        echo "-- »${yellow}tmpi $1 gdb --args mpirun n_cpus binary-input-file${nc}«"
        echo "-- After stopping gdb, [ctrl+b], [&], [y] and »exit« will get you "
        echo "-- back to the initial command line interface."
    fi
    #
    echo "-- "
fi
