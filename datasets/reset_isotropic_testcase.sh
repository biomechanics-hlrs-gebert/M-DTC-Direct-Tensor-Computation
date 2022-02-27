#!/bin/bash
#------------------------------------------------------------------------------
# Johannes Gebert - HLRS - NUM - gebert@hlrs.de - dtc dataset reset
#------------------------------------------------------------------------------
# Reset a dataset for a subsequent computation
# Command arg 1 --> input dataset
# Command arg 2 --> new app name
#
# Only for use with scalar data @ .int4.st as dtc does not process other.
#
# This simple procedure DOES NOT reset data contained in the stream file! 
# If you need 100% reliable data with validity up to the outermost voxels, 
# backup the initial pristine dataset.
#
# NO CHECKS. USE AT YOUR OWN RISK!
#------------------------------------------------------------------------------
#
# set -xv
#
#------------------------------------------------------------------------------
# Solution with two meta names:
#------------------------------------------------------------------------------
bsnmin=$(basename "$1" .meta)
bsnmout=$(basename "$2" .meta)
#------------------------------------------------------------------------------
# Solution wit app and input meta name
#------------------------------------------------------------------------------
bsnmin=$(basename "$2" .meta)
#
bsnmP1=$(echo "$bsnmin" | cut -d '_' -f 1 )
bsnmP2=$(echo "$bsnmin" | cut -d '_' -f 2 )
bsnmP5=$(grep 'NEW_BSNM_FEATURE' "$2" | tr -s ' ' |  cut -d ' '  -f 3 )
bsnmP3=$(grep 'NEW_BSNM_PURPOSE' "$2" | tr -s ' ' |  cut -d ' '  -f 3 )
bsnmP4="$1"
#
bsnmout="$bsnmP1""_""$bsnmP2""_""$bsnmP3""_""$bsnmP4""_""$bsnmP5"
# echo "$bsnmout"
#------------------------------------------------------------------------------
#
#
if [ -f "$bsnmin.int4.st" ]; then
    mv "$bsnmin".int4.st "$bsnmin".raw
    echo " $bsnmin.int4.st resetted."
fi
#
fls=( ".real8.st" ".char.st" ".memlog" ".head" )
#
for fl in "${fls[@]}"; do
    if [ -f "./$bsnmin$fl" ]; then
        rm ./"$bsnmin""$fl"
        echo " $bsnmin$fl resetted."
    fi
done
#
if [ -f ".$bsnmin.lock" ]; then
    rm ."$bsnmin".lock
    echo ".$bsnmin.lock resetted."
fi
#
if [ -d "$bsnmout" ]; then
    rm -r "$bsnmout"*
fi
#
if  rm temporary.* > /dev/null 2> /dev/null; then
    echo "temporary files resetted."
fi