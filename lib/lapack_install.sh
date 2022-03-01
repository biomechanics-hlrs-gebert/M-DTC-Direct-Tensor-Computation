#!/bin/bash
# ----------------------------------------------------------------------------------------
# Lapack download, build and install script.
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
# »sudo chmod +x lapack_install.sh«
#
# ----------------------------------------------------------------------------------------
#
# usage: »./<this script>«
# The Script will ask for sudo. It's required to create the target directory at /opt
# ----------------------------------------------------------------------------------------
#
# Modify only in case you know what you like to achieve :-)
# ----------------------------------------------------------------------------------------
#
# Specify the Version of Lapack
VERSION="3.10.0"
# ----------------------------------------------------------------------------------------
#
# Check whether all programs and compilers are accesible
declare -a BIN=("gfortran" "gcc" "g++" "cmake" "make" "tar" "wget")
#
echo "----------------------------------------------------------------------------------------"
for i in "${BIN[@]}"
do
   which "$i" > /dev/null 2> /dev/null
   if [ $? -eq 1 ]; then
	echo "Please make »$i« available."
	exit
   fi
done

filename=v$VERSION.tar.gz

if [ ! -f "$PWD"/"$filename" ]; then
		wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/$filename
	if [ $? -eq 1 ]; then
		echo ""
		echo "Download failed. Please have a look for the tarball of lapack and place it in the directory of this script."
		echo "Take a look at: http://www.netlib.org/lapack/#_lapack_version_3_10_0_2"
		exit 1
	fi
fi

# Installation target i.e. /opt/lapack
PREPREF=/home/geb/09_trinity_site

sudo mkdir -p $PREPREF >/dev/null 2>/dev/null &
sudo chown "$USER":"$USER" -R $PREPREF

PREF=$PREPREF/lapack-$VERSION

tar -xf $PWD/v$VERSION.tar.gz

mkdir lapack-$VERSION-build

cd lapack-$VERSION-build || exit 2

cmake ../lapack-$VERSION -DCMAKE_INSTALL_PREFIX=$PREF

make install

cd .. || exit 2

if rm -r ./lapack-$VERSION* >/dev/null 2>/dev/null & then
	echo ""
	echo "----------------------------------------------------------------------------------------"
	echo "Lapack successfully installed in:"
	echo ""
	echo $PREF
	echo ""
	echo "Append $PREF/lib to LD_LIBRARY_PATH."
	echo "Often accomplished via an environment script, specifically written to a program."
else
	echo "The Installation of Lapack failed."
fi

echo "----------------------------------------------------------------------------------------"
