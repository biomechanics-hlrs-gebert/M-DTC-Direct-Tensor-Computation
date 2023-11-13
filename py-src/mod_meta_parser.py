#!/bin/python
#------------------------------------------------------------------------------
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#
# @description: 
#> Library of domain decomposition stuff
# -----------------------------------------------------------------------------
#
import os
import numpy as np
#
KCL     = 25   # Keyword character  length
UCL     = 8    # Unit    character  length
STD_SPC = 39   # Keyword standard space
#
def parse_meta(meta_file, program_requested):

    if not os.path.exists(meta_file):
        print("EE File »" + meta_file + "« not found.")
        exit(1)

    with open(meta_file) as file:
        lines = [line.rstrip() for line in file]

    keywords = {}

    # -----------------------------------------------------------------------------
    # Should be the same implementation like in Fortran (!)
    # -----------------------------------------------------------------------------
    programs_scope = False

    for line in lines:
        kywd_found = False
        values_found = False

        try:
            tag = line[0]
        except: 
            continue

        try:
            kywd_in = line[2:KCL-1]
            kywd_in.strip()
        except: 
            continue

        try:
            vals = line[KCL-1 : KCL+STD_SPC-1]
            values = vals.split()

            values_found = True
        except: 
            continue



        if tag in ['*', 'd', 'r', 'w']:
            keyword = kywd_in.strip()
            kywd_found = True

        if tag in ['p'] and kywd_in == program_requested:
            programs_scope = True
            continue

        elif tag in ['p'] and kywd_in != program_requested and programs_scope:
            programs_scope = False

            #############################################################
            ## Requested program can occur only once in the meta file! ##
            #############################################################
            break 

        if kywd_found and values_found:
            # A keyword gets overwritten by this method!
            keywords.update({keyword: values})

    return(keywords)
