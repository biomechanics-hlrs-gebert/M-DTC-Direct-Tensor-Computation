#!/bin/python
#------------------------------------------------------------------------------
# write_task_geometry.DTC.hawk.py
#------------------------------------------------------------------------------
#> @author Johannes Gebert - HLRS - NUM - gebert@hlrs.de
#
# @description: 
#> Create a file containing the task geometry for an HPE MPT batch job. 
#
# -----------------------------------------------------------------------------
# Import modules
# -----------------------------------------------------------------------------
import os
import argparse
import math

#
# -----------------------------------------------------------------------------
# Hawk specific parameters (may be altered for another cluster/HPC system)
# -----------------------------------------------------------------------------
processors_per_node=128
master_rank=0                # Depends on the MPI implementation of DTC!

# -----------------------------------------------------------------------------
# Set the parser for the command line arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Create a file containing the task geometry for an HPE MPT batch job.')
parser.add_argument('-if', action="store", dest="in_file", help="Input batch file (PBS).")
parser.add_argument('-tf', action="store", dest="task_file", help="Input batch file (PBS).")

# -----------------------------------------------------------------------------
# Parse the command line arguments
# -----------------------------------------------------------------------------
cmd_args = parser.parse_args()
 
# -----------------------------------------------------------------------------
# Check if the cmd line args are valid/make sense
# -----------------------------------------------------------------------------
if not bool(cmd_args.in_file):
    print("Please specify an input file.")
    exit(1)

if not os.path.exists(cmd_args.in_file):
    print("Input file »" + cmd_args.in_file + "« does not exist.")
    exit(2)

if not bool(cmd_args.task_file):
    print("Please specify an input file.")
    exit(3)

if not os.path.exists(cmd_args.task_file):
    print("Input file »" + cmd_args.task_file + "« does not exist.")
    exit(4)


# -----------------------------------------------------------------------------
# Parse the batch file.
# -----------------------------------------------------------------------------
with open(cmd_args.in_file) as file:
    lines = [line.rstrip() for line in file]

    size_mpi=0
    ppd=0
    parse_keywords = False
    
    for line in lines:

        # -----------------------------------------------------------------------------
        # Parse the node ID.
        # -----------------------------------------------------------------------------
        keyword = int(line[2:23])

        # -----------------------------------------------------------------------------
        # Stop parsing if the meta file contains another program.
        # -----------------------------------------------------------------------------
        if int(line[0]) == "p" and parse_keywords:
            parse_keywords = False

        # -----------------------------------------------------------------------------
        # Start parsing if the keywords belong to DTC.
        # -----------------------------------------------------------------------------
        if keyword.strip() == "TENSOR-COMPUTATION":
            parse_keywords = True

        if keyword.strip() == "PROCESSORS" and parse_keywords:
            clipboard=line[24:65]

            size_mpi=int(clipboard.strip())

            print( "MM PROCESSORS of " + cmd_args.in_file + "=" + str(size_mpi) + "\n")


        if keyword.strip() == "MESH_PER_SUB_DMN" and parse_keywords:
            clipboard=line[24:65]

            ppd=int(clipboard.strip())

            print( "MM MESH_PER_SUB_DMN of " + cmd_args.in_file + "=" + str(ppd) + "\n")
       

    # -----------------------------------------------------------------------------
    # Check presence of the variables
    # -----------------------------------------------------------------------------
    if size_mpi == 0:
        print ("EE Keyword PROCESSORS not parsed properly.")
        exit(6)
    
    if ppd == 0:
        print ("EE Keyword MESH_PER_SUB_DMN not parsed properly.")
        exit(5)


# -----------------------------------------------------------------------------
# Calculate the task geometry:
# Example provided by comments.
# -----------------------------------------------------------------------------
if ppd > processors_per_node:

    # 5 = math.ceil(598/128)
    no_nodes_per_domain = math.ceil(ppd/processors_per_node)

    # 119 = math.floor(598/5)
    processes_per_node = math.floor(ppd/no_nodes_per_domain)

    # 595 = 119*5
    no_processes_per_domain = processes_per_node*no_nodes_per_domain

    # 3 = 598 - 595
    # 3 of 5 nodes must host an additional process
    overhead = ppd - no_processes_per_domain



# -----------------------------------------------------------------------------
# Concatenate the task geometry like {(0)(1,2,5)(4)}:
# -----------------------------------------------------------------------------
task_geom_string="{(" + str(master_rank) + ")"




task_geom_string = task_geom_string + "}"

# -----------------------------------------------------------------------------
# Write the *.task_geometry file
# -----------------------------------------------------------------------------
with open(cmd_args.task_file, 'w') as f:  # Just use 'w' mode in 3.x

