#!/bin/bash
# ----------------------------------------------------------------------------------------
# Extract max memory from log
# Author:      Johannes Gebert gebert@hlrs.de
# Created: on: 31.01.2022
# Last edit:   31.01.2022
# ----------------------------------------------------------------------------------------
if [[ -z $1 ]]; then
    echo "Please specify a file to read."
else
    if [[ ! -f $1 ]]; then
        echo "The file $1 does not exist. Please select a valid one."
    else
        GLBL_MAX=0
        num1=0
        #
        while IFS= read -r line
        do
            num1=$(echo "$line" | tr -s ' ' | cut -d " " -f 3 | tr -d Gi)
            # echo "$num1"
            #
            if [[ $num1 != *.* ]]  > /dev/null 2> /dev/null ; then
                num1=$num1".0"
            fi
            # echo "$num1"
            #
            if (( $(echo "$num1 > $GLBL_MAX" | bc -l) )); then
                GLBL_MAX=$num1
            fi
        done < "$1"

    fi
fi
#
echo "Maximum memory usage logged in file $1: $GLBL_MAX GiB."