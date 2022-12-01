#!/bin/python
#------------------------------------------------------------------------------
# create_directories.DTC.hawk.py
#------------------------------------------------------------------------------
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#
# @description: 
#> Create a file containing the task geometry for an HPE MPT batch job. 
#> Reads its information out of the meta file of the batch job.
#
# -----------------------------------------------------------------------------
# Import modules
# -----------------------------------------------------------------------------
import os, argparse
from random import randint
from time import sleep

# -----------------------------------------------------------------------------
# Set the parser for the command line arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Create directories for the struct process.')
parser.add_argument('-p_n_bsnm', action="store", dest="p_n_bsnm", help="Path and basename of the meta file.")
parser.add_argument('-worker_bsnm', action="store", dest="Worker_bsnm", help="Basename of the worker result tree.")
parser.add_argument('-ppd', action="store", dest="ppd", help="Parts per domain")
parser.add_argument('-size_mpi', action="store", dest="size_mpi", help="Size_mpi.")

# -----------------------------------------------------------------------------
# Parse the command line arguments
# -----------------------------------------------------------------------------
cmd_args = parser.parse_args()
 
# -----------------------------------------------------------------------------
# Create the directories for the struct process
# -----------------------------------------------------------------------------
for ii in range(1,int(cmd_args.size_mpi),int(cmd_args.ppd)):

    number = str(ii)

    Rank_directories = cmd_args.p_n_bsnm + "/" + cmd_args.Worker_bsnm + "_" + number.zfill(7)

    # sleep(randint(1,10))

    if not os.path.exists(Rank_directories):
        os.makedirs(Rank_directories)

        print("Create directory: " + Rank_directories)

