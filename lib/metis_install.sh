#!/bin/bash
# ----------------------------------------------------------------------------------------
# Metis download, build and install script.
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
# »sudo chmod +x Metis_install.sh«
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
# Specify the Version of Metis
VERSION="5.1.0"
# ----------------------------------------------------------------------------------------
#
# Check whether all programs and compilers are accesible
declare -a BIN=("gfortran" "gcc" "g++" "make" "tar" "wget")
#
echo "----------------------------------------------------------------------------------------"
for i in "${BIN[@]}"
do
    if ! which "$i" > /dev/null 2> /dev/null; then
	echo "Please make »$i« available."
	exit 1
   fi
done

filename=metis-$VERSION.tar.gz

PREF=/opt/metis/metis-$VERSION

sudo mkdir -p $PREF >/dev/null 2>/dev/null &

sudo chown $USER:$USER -R $PREF

cd $PREF || exit


if [ ! -f $PREF/$filename ]; then
	wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/$filename
	if [ $? -eq 1 ]; then
		echo ""
		echo "Download failed. Please have a look for the tarball of Metis and place it in the directory of this script."
		echo "Take a look at: http://glaros.dtc.umn.edu/gkhome/metis/metis/download"
		exit 1
	fi
fi

tar -xvzf $filename

cd $PREF/metis-$VERSION/build || exit 2

sed -i 's/REALTYPEWIDTH 32/REALTYPEWIDTH 64/g' $PREF/metis-$VERSION/include/metis.h
sed -i 's/IDXTYPEWIDTH 32/IDXTYPEWIDTH 64/g' $PREF/metis-$VERSION/include/metis.h

cmake -DGKLIB_PATH=../GKlib  -DCMAKE_INSTALL_PREFIX=$PREF ../

make install

cd .. || exit 2

# Grep after installation to provide a proper output to cli
echo "----------------------------------------------------------------------------------------"
echo

if ! grep "#define REALTYPEWIDTH 64" $PWD/include/metis.h; then
	echo "REALTYPEWIDTH 64 of metis.h may be errornous. Please check!"
	exit 1
fi


if ! grep "#define IDXTYPEWIDTH 64" $PWD/include/metis.h; then
	echo "IDXTYPEWIDTH 64 of metis.h may be errornous. Please check!"
	exit 1
fi

if rm -r "$PREF"/"$filename" >/dev/null 2>/dev/null & then
	echo "----------------------------------------------------------------------------------------"
	echo "Metis successfully installed in:"
	echo ""
	echo $PREF
	echo ""
	echo "Append $PREF/bin to PATH and"
	echo "Append $PREF/lib to LD_LIBRARY_PATH."
	echo "Often accomplished via an environment script, specifically written to a program."
else
	echo "The Installation of Metis failed."
fi
echo "----------------------------------------------------------------------------------------"
