#!/bin/bash
# ----------------------------------------------------------------------------------------
# Manage the state of a subtree
#
# Author:    Johannes Gebert - HLRS - NUM - gebert@hlrs.de
# Created:   27.12.2021
# Last edit: 23.02.2022
# ----------------------------------------------------------------------------------------
SubtreeUrl=git@github.com:biomechanics-hlrs-gebert/A-CESO-Central_Sources.git
# ----------------------------------------------------------------------------------------
# Update the subtree
#
if ! which git > /dev/null 2> /dev/null ; then 
    echo "Git not found."
    echo "The program cannot get updates from this directory."
    echo "The program cannot compile if the »central_src« directory ist missing."
fi
#
if  ! ls -l "$PWD"/central_src > /dev/null 2> /dev/null ; then 
    operation="add"
else
    operation="pull"
fi
#
git subtree $operation --prefix central_src $SubtreeUrl main --squash