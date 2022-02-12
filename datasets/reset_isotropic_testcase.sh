#!/bin/bash
#------------------------------------------------------------------------------
# Johannes Gebert - HLRS - NUM - gebert@hlrs.de - ddtc dataset reset
#------------------------------------------------------------------------------
# Reset a dataset for a subsequent computation
# Command arg 1 --> input dataset
# Command arg 2 --> output dataset.
#
# Only for use with scalar data @ .int4.st as ddtc does not process other.
#
# This simple procedure DOES NOT reset data contained in the stream file! 
# If you need 100% reliable data with validity up to the outermost voxels, 
# backup the initial pristine dataset.
#
# NO CHECKS. USE AT YOUR OWN RISK!
#------------------------------------------------------------------------------
#
bsnmin=$(basename $1 .meta)
bsnmout=$(basename $2 .meta)
#
if [ -f "$bsnmin.int4.st" ]; then
    mv $bsnmin.int4.st $bsnmin.raw
    echo " $bsnmin.int4.st resetted."
fi
#
fls=( ".real8.st" ".char.st" ".memlog" ".head" )
#
for fl in "${fls[@]}"; do
    if [ -f "./$bsnmin$fl" ]; then
        rm ./$bsnmin$fl
        echo " $bsnmin$fl resetted."
    fi
done
#
if [ -f ".$bsnmin.lock" ]; then
    rm .$bsnmin.lock
    echo ".$bsnmin.lock resetted."
fi
#
if [ -d "$bsnmout" ]; then
    rm -r $bsnmout*
fi
#
rm temporary.* > /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
    echo "temporary files resetted."
fi