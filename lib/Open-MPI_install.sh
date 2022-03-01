#!/bin/bash
# ----------------------------------------------------------------------------------------
# Open-MPI download, build and install script.
#
# Author:     Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created on: 14.09.2021
# Last edit:  27.02.2022
# ----------------------------------------------------------------------------------------
#
# Requirements:
# GCC
# GNU Make
# tar
# wget
# Ensure, your OS is allowed to run this script:
# sudo chmod +x Open-MPI_install.sh
#
# ----------------------------------------------------------------------------------------
#
# usage: ./<this script>
# The Script will ask for sudo. It's required to create the target directory at /opt
# ----------------------------------------------------------------------------------------
# Specify the Version of MPI
VERSION="4.1.2"
# ----------------------------------------------------------------------------------------
# Specify the general path of your mpi installation
base_prefix="/opt/mpi"
# ----------------------------------------------------------------------------------------
# Specify the integer kind of mpi (4/8, single/double precision) 
MPI_INT_KIND=4
# ----------------------------------------------------------------------------------------
# Specify whether mpi is build wih mpi_f08 or not (0=off, 1=on)
# ONLY USE MPI_F08 if MPI_INT_KIND=4 or c-integer-kind=8 as well!!
MPI_F08_SWITCH=0
# ----------------------------------------------------------------------------------------
#
#
#
# Modify only in case you know what you like to achieve :-)
# ----------------------------------------------------------------------------------------
# Check whether all programs and compilers are accesible
declare -a BIN=("gfortran" "gcc" "g++" "make" "tar" "wget")

echo "----------------------------------------------------------------------------------------"
for i in "${BIN[@]}"
do
   which "$i" > /dev/null 2> /dev/null
   if [ $? -eq 1 ]; then
	echo "Please make «"$i"» available."
	exit=1
   fi
done

filename=openmpi-$VERSION.tar.gz

if [ ! -f $PWD/$filename ]; then
	wget https://download.open-mpi.org/release/open-mpi/v4.1/$filename
	if [ $? -eq 1 ]; then
		echo ""
		echo "Download failed. Please have a look for the tarball of Open-MPI and place it in the directory of this script."
		echo "Take a look at: https://www.open-mpi.org/software/ompi/v4.1/"
		ABORT=1
	fi
fi

sudo mkdir -p $base_prefix >/dev/null 2>/dev/null &
sudo chown $USER:$USER -R $base_prefix

tar -xvzf $filename
cd openmpi-$VERSION

# Configure MPI.
# --prefix=<installation directory, independent of this one.>
# further options via ./configure --help
if [ $MPI_INT_KIND -eq 4 ]; then
	export F77=gfortran
	export FC=gfortran
	export CC=gcc
	export CXX=g++
	export FFLAGS="-m64"   # Real kind = 8
	export FCFLAGS="-m64"  # Real kind = 8
	export CFLAGS="-m64"   # Real kind = 8
	export CXXFLAGS="-m64" # Real kind = 8
	export mpifkkind_path="-I4"
elif [ $MPI_INT_KIND -eq 8 ]; then
	export F77=gfortran
	export FC=gfortran
	export CC=gcc
	export CXX=g++
	export FFLAGS="-m64 -fdefault-integer-8 --with-wrapper-fflags=-fdefault-integer-8 --with-wrapper-fcflags=-fdefault-integer-8"
	export FCFLAGS="-m64 -fdefault-integer-8"
	export CFLAGS="-m64"
	export CXXFLAGS="-m64"
	export mpifkkind_path="-I8"
fi
#
if [ $MPI_F08_SWITCH -eq 0 ]; then
	export mpif08_onoff="--enable-mpi-fortran=usempi"
	export mpif08_path=""
elif [ $MPI_F08_SWITCH -eq 1 ]; then
	export mpif08_onoff="--enable-mpi-fortran=usempif08"
	export mpif08_path="-MPIF08"
fi
installprefix=$base_prefix/openmpi$mpif08_path$mpifkkind_path-$VERSION

./configure --prefix=$installprefix 

make && make install

if [ $? -eq 0 ]; then
	rm -r ./openmpi-$VERSION* >/dev/null 2>/dev/null &
	echo ""
	echo "----------------------------------------------------------------------------------------"
	echo "Open-MPI successfully installed in:"
	echo ""
	echo $installprefix
	echo ""
	echo "Append "$installprefix"/bin to PATH and"
	echo "Append "$installprefix"/lib to LD_LIBRARY_PATH."
	echo "Often accomplished via an environment script, specifically written to a program."
else
	echo "The Installation of Open-MPI failed."
fi
echo "----------------------------------------------------------------------------------------"
