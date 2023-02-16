#!/bin/python
#------------------------------------------------------------------------------
# write_task_geometry.DTC.hawk.py
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
import os, argparse, math

#
# -----------------------------------------------------------------------------
# Hawk specific parameters (may be altered for another cluster/HPC system)
# -----------------------------------------------------------------------------
processors_per_node=128
master_rank=0          
#
FMT_STRING="-------------------------------------------------------------------------------"      # Depends on the MPI implementation of DTC!
#
# -----------------------------------------------------------------------------
# User output
# -----------------------------------------------------------------------------
print(
"""
-------------------------------------------------------------------------------
-- Creating hostfile | High-Performance Computing Center | Stuttgart (HLRS)
-------------------------------------------------------------------------------""")

# -----------------------------------------------------------------------------
# Set the parser for the command line arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Create a file containing the task geometry for an HPE MPT batch job.')
parser.add_argument('-if', action="store", dest="in_file", help="Input batch file (PBS).")
parser.add_argument('-nodelist', action="store", dest="nodelist", help="Nodefile of the job.")
parser.add_argument('-nodes', action="store", dest="nodes_available", help="Number of nodes of job submission.")
parser.add_argument('-ppn', action="store", dest="processors", help="Total number the number of processors used.")
parser.add_argument('-hf', action="store", dest="hostfile", help="Input batch file (PBS).")

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

if not bool(cmd_args.hostfile):
    print("Please specify an output file.")
    exit(1)

if os.path.exists(cmd_args.hostfile):
    print("The task_geometry file " + cmd_args.hostfile + " already exists.")
    # THIS IS NOT A BUG
    # exit(2)

if not bool(cmd_args.nodelist):
    print("Please specify a nodelist.")
    exit(1)

if not bool(cmd_args.nodes_available):
    print("Please specify the number of nodes available.")
    exit(3)

if not bool(cmd_args.processors):
    print("Please specify the number of processors used.")
    exit(3)
else:
    try:
        processors=int(cmd_args.processors.strip())

    except:
        pass
        
# -----------------------------------------------------------------------------
# Parse the batch file.
# -----------------------------------------------------------------------------
with open(cmd_args.in_file) as file:
    lines = [line.rstrip() for line in file]

    ppd=0
    parse_keywords = False
    
    for line in lines:

        # -----------------------------------------------------------------------------
        # Parse the node ID.
        # -----------------------------------------------------------------------------
        keyword = line[2:23]

        # -----------------------------------------------------------------------------
        # Stop parsing if the meta file contains another program.
        # -----------------------------------------------------------------------------
        try:
            if line[0] == "p" and parse_keywords:
                parse_keywords = False
        except:
            pass

        # -----------------------------------------------------------------------------
        # Start parsing if the keywords belong to DTC.
        # -----------------------------------------------------------------------------
        if keyword.strip() == "TENSOR_COMPUTATION":
            parse_keywords = True

        if keyword.strip() == "MESH_PER_SUB_DMN" and parse_keywords:
            try:
                clipboard=line[24:65]

                ppd=int(clipboard.strip())

            except:
                pass

    # -----------------------------------------------------------------------------
    # Check presence of the variables
    # -----------------------------------------------------------------------------
    if processors == 0:
        print ("Number of processors wrong.")
        exit(4)
    
    if ppd == 0:
        print ("EE Keyword MESH_PER_SUB_DMN not parsed properly.")
        print(FMT_STRING)
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
    
    processors_per_node=processes_per_node


    if overhead != 0:
        print ("WW The overhead is not 0. Some MPI processes may be placed unexpectedly.")
        print(FMT_STRING)


else:
    no_nodes_per_domain = 1
    overhead = 0
    processors_per_node=ppd




# -----------------------------------------------------------------------------
# Write the *.task_geometry file
# -----------------------------------------------------------------------------
out_file = open(cmd_args.hostfile, "w")

total_slots_available=0

with open(cmd_args.nodelist) as file:
    lines = [line.rstrip() for line in file]

    firstline= True
    for line in lines:

        if firstline:
            slots=1
            firstline=False
        else:
            slots=processors_per_node
        
        total_slots_available = total_slots_available + slots

        out_file.write(line + " slots=" + str(slots) + "\n")
 
out_file.close()
