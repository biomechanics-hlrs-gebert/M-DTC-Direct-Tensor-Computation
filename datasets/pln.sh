#!/bin/bash
# ----------------------------------------------------------------------------------------
# Print lines. I'm quite sure there are far more convenient ways to program that:-)
#
# Author:          Johannes Gebert gebert@hlrs.de
# Created: on :    13.04.2021
# Last edit:       12.02.2022
# ----------------------------------------------------------------------------------------
#
function echo_help {
    echo ""
    echo "Print a specific amount of lines."
    echo ""
    echo "Usage:"
    echo "pln      file  first-line  amount of lines"
    echo "pln  -a  file  first-line  last-line"
    echo "pln      file  start       last-line"
    echo "pln      file  first-line  end"
    echo ""
}
#
# Default:
calc_lines=1
#
# if [ -z $3 ]; then
#     echo_help
#     exit 0
# fi
#
for arg in "$@"
do
    if [ "$arg" == "help" ] || [ "$arg" == "--help" ] || [ "$arg" == "-h" ] || [ "$arg" == "h" ]; then
	    echo_help
        exit 0
    fi
    #
    # Use absolute line numbers:
    if [ "$arg" == "-a" ]; then
        calc_lines=0
    fi
    #
    # Print from line 1:
    if [ "$arg" == "-s" ] || [ "$arg" == "--start" ]; then
        file=$1
        sed -n '1,'"$3"'p' "$file"
        exit 0
    fi
    #
    # Print till last line:
    if [ "$arg" == "-e" ] || [ "$arg" == "--end" ]; then
        file=$1
        sed -n "$2"',$p' "$file"
        exit 0
    fi
done
#
if [ $calc_lines -eq 1 ]; then
    # $1 = -a
    file=$1
    start=$2
    end=$(("$2" + "$3"))
else
    if [ "$1" == "-a" ]; then
        file=$2
        #
        start=$3
        end=$4
    else
        echo_help
    fi
fi
#
sed -n "${start},${end}p" "$file"
