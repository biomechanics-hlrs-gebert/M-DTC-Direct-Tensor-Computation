#!/bin/bash
#------------------------------------------------------------------------------
# Johannes Gebert - HLRS - NUM - gebert@hlrs.de - dtc dataset reset
#------------------------------------------------------------------------------
# Reset a (crawled) dataset.
#------------------------------------------------------------------------------
#
declare -a basenames=(\
    "FH01-2_mu_Dev_CTIF_G31-Sig10" \
    "FH01-2_mu_dev_dtc_processors186" \
    "FH01-2_mu_dev_dtc_processors93" \
    "FH01-2_mu_dev_dtc_processors46" \
    )
#
declare -a files=(\
    "memlog" \
    "std_out" \
    "std_err" \
    "head" \
    "status" \
    "status*.loc*" \
    "mon" \
    "st" \
    )
#
for ((ii=0; ii<${#basenames[@]}; ii++));
do
    #
    #------------------------------------------------------------------------------
    # Delete files
    #------------------------------------------------------------------------------
    for suf in "${files[@]}"
    do
        file="$PWD/${basenames[ii]}.$suf"
        if [ -f "$file" ]; then
            echo "Resetting $file"
            rm "$file" >/dev/null 2>/dev/null &
        fi
    done
    #
    #------------------------------------------------------------------------------
    # Delete status file
    #------------------------------------------------------------------------------
    echo "Resetting $PWD/${basenames[ii]}.status"*".loc"*
    rm "$PWD/${basenames[ii]}.status"*".loc"* >/dev/null 2>/dev/null &
    #
    #------------------------------------------------------------------------------
    # Delete directory
    #------------------------------------------------------------------------------
    directory="$PWD/${basenames[ii]}"
    if [ -d "$directory" ]; then
        echo "Resetting $directory"
        rm -r "$directory" >/dev/null 2>/dev/null &
    fi
done