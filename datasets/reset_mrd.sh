#!/bin/bash
#------------------------------------------------------------------------------
# Johannes Gebert - HLRS - NUM - gebert@hlrs.de - dtc dataset reset
#------------------------------------------------------------------------------
# Reset a crawled dataset.
# Command arg 1 --> input dataset
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
#------------------------------------------------------------------------------
#
fls=( ".covo" ".cr1" ".cr2" )
#
for fl in "${fls[@]}"; do
    if [ -f "./$bsnmin$fl" ]; then
        rm ./"$bsnmin""$fl"
        echo " $bsnmin$fl resetted."
    fi
    #
    if [ -f "./temporary$fl" ]; then
        rm ./"temporary""$fl"
        echo " temporary$fl resetted."
    fi
done
#
if [ -f "./.$bsnmin.lock" ]; then
    rm ./".$bsnmin.lock"
    echo ".$bsnmin.lock resetted."
fi